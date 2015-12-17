classdef FirstOrderDynamics < dscomponents.ACoreFun
    
    properties(Dependent)
        UseFrequencyDetector;
    end
    
    properties
        MSLink;
        
        % The factors with which the primary and secondary affarent of the
        % spindle is multiplied before considered a "mean input current"
        % (then added to the external signal)
        %
        % @type rowvec<double> @default [1 1]*0.002
        SpindleAffarentWeights = sparse([1 1]*0.002);
    end
    
    properties(Access=private)
        mslinkfun;
        mslinkfunderiv;
        max_spindle_signals;
        
        fUseFD = false;
        nfibres;
        
        x0_motorunit;
    end
    
    properties(SetAccess=private)
        FrequencyDetector;
        freq_kexp;
        
        num_motoneuron_dof;
        num_sarco_dof;
        num_spindle_dof;
        
        input_motoneuron_link_idx;
        moto_input_noise_factors;
        sarco_output_idx;
        
        moto_sarco_link_moto_out;
        moto_sarco_link_sarco_in;
        
        % Motoneuron quantities
        moto_signal_input_pos;
        SomaInputFactors;
        Noise;
        IndepNoise_AP;
        normalized_cortex_signal;
        
        % The initial value.
        x0;
        
        % The offset for the sarcomere signals at t=0. Used to create
        % alpha(X,0) = 0
        %
        % See also: assembleX0
        sarco_mech_signal_offset;
    end
    
    methods
        function this = FirstOrderDynamics(sys)
            this = this@dscomponents.ACoreFun(sys);
            this.MSLink = models.motorunit.MSLink;
            
            s = load(fullfile(fileparts(which('fullmuscle.Model')),'FrequencyKexp.mat'));
            this.freq_kexp = s.kexp.toTranslateBase;
            
            % Get fibretype-independent noise from NoiseGenerator
            ng = models.motoneuron.NoiseGenerator;
            this.IndepNoise_AP = ng.AP;
        end
        
        function configUpdated(this)
            sys = this.System;
            mc = sys.Model.Config;
            if ~isempty(mc)
                ft = mc.FibreTypes;
                this.nfibres = length(ft);
                this.xDim = sys.NumTotalDofs;
                this.assembleX0;
                this.JSparsityPattern = this.computeSparsityPattern;
                
                this.FrequencyDetector = models.fullmuscle.FrequencyDetector(this.nfibres);
                
                %% Assemble noise signal for each fibre
                cs = mc.NormalizedCortexSignal;
                if isa(cs,'general.functions.AFunGen')
                    cs = cs.getFunction;
                elseif isempty(cs)
                    error('The NormalizedCortexSignal field of the models config is empty.');
                end
                this.normalized_cortex_signal = cs;
                
                % Load the extra spindle upper limit poly (larger than the
                % normal one by offset MaxExtraMeanSpindleSignal) and save
                % the maximum allowed spindle signals
                s = load(models.motoneuron.Model.FILE_UPPERLIMITPOLY);
                poly = s.upperlimit_poly + mc.MaxExtraMeanSpindleSignal;
                this.max_spindle_signals = polyval(poly, ft);
                
                ng = models.motoneuron.NoiseGenerator;
                ng.setFibreType(ft(1));
                thenoise = zeros(this.nfibres,length(ng.totalNoise));
                thenoise(1,:) = ng.totalNoise;
                for k=2:this.nfibres
                    ng.setFibreType(ft(k));
                    thenoise(k,:) = ng.totalNoise;
                end
                this.Noise = thenoise;
            
                % The noise signal added to the soma is multiplied with a
                % fibretype-dependent factor (constant 8).
                %this.SomaInputFactors = 1./(pi*ls.^2);
                mo = sys.Motoneuron;
                this.SomaInputFactors = 1./mo.getC8Value(ft);
                
                %% Component link indices
                % Get the positions where the input signal is mapped to the
                % motoneurons
                dm = sys.Motoneuron.Dims;
                dsa = sys.Sarcomere.Dims;
                this.input_motoneuron_link_idx = 2:dm:this.num_motoneuron_dof;
                this.moto_signal_input_pos = 2:dm:this.num_motoneuron_dof;
                this.moto_sarco_link_moto_out = sys.EndSecondOrderDofs + this.moto_signal_input_pos;
                this.moto_sarco_link_sarco_in = this.num_motoneuron_dof ...
                    + (1:dsa:this.num_sarco_dof);
            end
        end
       
        function prepareSimulation(this, mu)
            prepareSimulation@dscomponents.ACoreFun(this, mu);
            [this.mslinkfun, this.mslinkfunderiv] = this.MSLink.getFunction;
            
            % Register the ODE callback for the frequency integration
            sys = this.System;
            slv = sys.Model.ODESolver;
            if ~isa(slv,'solvers.MLWrapper')
                error('Only programmed to work with ML-builtin solvers so far!');
            end
            if this.fUseFD
                fd = this.FrequencyDetector;
                slv.odeopts = odeset(slv.odeopts,...
                    'OutputFcn',@(t,y,flag)fd.processSignal(t,y'),...
                    'OutputSel',this.moto_sarco_link_moto_out);
            else
                slv.odeopts.OutputFcn = [];
            end
        end

        function dy = evaluate(this, y, t)
            sys = this.System;
            nf = sys.nfibres;
            dy = zeros(this.fDim,1);
            off_mech = sys.EndSecondOrderDofs;
            
            %% Motoneurons
            mo = sys.Motoneuron;
            ymoto = y(off_mech + (1:this.num_motoneuron_dof));
            dy_m = mo.dydt(reshape(ymoto,mo.Dims,[]),t);
            dy(1:this.num_motoneuron_dof) = dy_m(:);
            
            %% Sacromeres
            sa = sys.Sarcomere;
            sarco_pos = this.num_motoneuron_dof + (1:this.num_sarco_dof);
            ys = reshape(y(off_mech + sarco_pos),sa.Dims,[]);
            dys = sa.dydt(ys, t);
            dy(sarco_pos) = dys(:);
            
            %% Link of motoneurons to sarcomeres
            moto_out = y(this.moto_sarco_link_moto_out);
            signal = this.mslinkfun(moto_out).*moto_out./sa.SarcoConst(1,:)';
            % Add signal to corresponding locations
            dy(this.moto_sarco_link_sarco_in) = ...
                dy(this.moto_sarco_link_sarco_in) + signal;
            
            %% External input signal
            % Get current external signal
            moto_signal = mo.checkMeanCurrent(this.mu(4)*this.normalized_cortex_signal(t));
            
            %% Spindles
            if sys.HasSpindle    
                sp = sys.Spindle;
                spindle_pos = this.num_motoneuron_dof ...
                    + this.num_sarco_dof + (1:this.num_spindle_dof);
                yspindle = reshape(y(off_mech + spindle_pos),sp.Dims,[]);

                %% Spindles -> Motoneurons
                % Get single spindle signals
                spindle_sig = this.SpindleAffarentWeights*sp.getAfferents(yspindle);
                % Compute the mean value over all signals
                spindle_sig = ones(1,nf)*mean(spindle_sig);
                % Compute the maximum possible effective spindle signal
                spindle_sig = min(spindle_sig,this.max_spindle_signals - moto_signal);
                % Add the signal to the current motoneuron signal
                moto_signal = moto_signal + spindle_sig;
                
                % Spindle self-dynamics
                % In the case of learned moto frequencies, we have an
                % instant feedback to the spindle. This is why we need the
                % spindle->moto signal first
                
                % Get motoneuron frequency
                if this.fUseFD
                    freq = this.FrequencyDetector.Frequency;
                else
                    freq = this.freq_kexp.evaluate([sys.Model.Config.FibreTypes; moto_signal]);
                end
                dys = sp.dydt(yspindle,t,freq,sys.f.lambda_dot,0);
                dy(spindle_pos) = dys(:);
            end
            
            %% Motoneuron exitation (external + spindle)
            indep_sig = moto_signal * this.IndepNoise_AP;
            fibre_dep_sig = 9*((moto_signal/9).^1.24635) .* this.Noise(:,round(t)+1)';
            dy(this.moto_signal_input_pos) = dy(this.moto_signal_input_pos) ...
                + ((indep_sig + fibre_dep_sig).*this.SomaInputFactors)';
        end
        
        function fx = evaluateCoreFun(this, x, t)
            error('Do not call directly; have custom evaluate method.');
        end
        
        function SP = computeSparsityPattern(this)
            sys = this.System;
            nf = this.nfibres;
            
            SP = sparse(this.fDim,this.xDim);
            [~, ~, SPLamDot] = sys.f.computeSparsityPattern;
            off_mech = sys.EndSecondOrderDofs;
            % Neuro
            dm = sys.Motoneuron.Dims;
            for k=1:nf
                pos = dm*(k-1)+(1:dm);
                SP(pos,off_mech + pos) = sys.Motoneuron.JSparsityPattern;%#ok
            end
            
            % Sarco
            J_sarco = sys.Sarcomere.JSparsityPattern;
            dsa = sys.Sarcomere.Dims;
            off = nf*dm;
            for k=1:nf
                pos = dsa*(k-1) + (1:dsa);
                SP(off+pos,off+off_mech + pos) = J_sarco;%#ok
            end
            
            % Spindle
            if sys.HasSpindle
                off_spindle = (dm+dsa)*nf;
                ds = sys.Spindle.Dims;
                JSp = sys.Spindle.JSparsityPattern;
                for k=1:nf
                    pos = ds*(k-1) + (1:ds);
                    SP(off_spindle+pos,off_spindle+off_mech + pos) = JSp;%#ok
                end
            end
            
            % Moto -> Sarco link
            % first entry of sarco gets 2nd output of motoneuron
            pos_sarco = nf*dm + (1:dsa:nf*dsa);
            pos_moto = off_mech + (1:dm:nf*dm) + 1;
            pos = sub2ind(size(SP),pos_sarco,pos_moto);
            SP(pos) = true; 
            
            if sys.HasSpindle
                sp = sys.Spindle;
                ds = sp.Dims;
                % Spindle -> Motoneuron link
                i = []; j = [];
                moto_pos = 2:dm:dm*nf;
                for k=1:nf
                    i = [i repmat(moto_pos',1,ds)];%#ok
                    j = [j repmat(ds*(k-1) + (1:ds),nf,1)];%#ok
                end
                SP(1:dm*nf,off_mech + (dm+dsa)*nf + (1:ds*nf)) = ...
                    sparse(i(:),j(:),ones(numel(i),1),dm*nf,ds*nf);
                
                % Mechanics -> spindle
                Jspin_Ldot = double(sp.JLdotSparsityPattern);
                Jspin_dmoto = double(sp.JMotoSparsityPattern);
                Jspin_Aff = double(sp.JAfferentSparsityPattern);
                for k=1:nf
                    spindle_pos = off_spindle + ds*(k-1) + (1:ds);
                    SP(spindle_pos,1:sys.NumStateDofs+sys.NumDerivativeDofs) = ...
                        logical(Jspin_Ldot*double(SPLamDot(k,:)));
                    % Create connecting link to self only when no frequency
                    % detector is used!
                    if ~this.fUseFD
                        spindle_pos_glob = off_mech + spindle_pos;
                        SP(spindle_pos,spindle_pos_glob) = ...
                            SP(spindle_pos,spindle_pos_glob) | logical(Jspin_dmoto*any(Jspin_Aff));
                    end
                end
            end
        end
        
        function J = getStateJacobian(this, y, t)
%             J = this.getStateJacobianFD(y,t);
%             return;

            sys = this.System;
            off_mech = sys.EndSecondOrderDofs;
            nf = sys.nfibres;
            
            %% Motoneuron
            mo = sys.Motoneuron;
            dm = mo.Dims;
            ymoto = reshape(y(off_mech + (1:this.num_motoneuron_dof)),dm,[]);
            Jdm = cell(1,nf);
            for k=1:nf
                Jdm{k} = mo.Jdydt(ymoto(:,k),t,k);
            end
            
            %% Sarcomeres
            sa = sys.Sarcomere;
            dsa = sa.Dims;
            sarco_pos = this.num_motoneuron_dof + (1:this.num_sarco_dof);
            ys = reshape(y(off_mech + sarco_pos),dsa,[]);
            if nf == 1
                Jsa{1} = sa.Jdydt(ys, t);
            else
                Jsa = cell(1,nf);
                [i,j] = find(sa.JSparsityPattern');
                % New, faster evaluation! Stick all into the Jdydt function
                % and create N sparse matrices from that!
                % Tricky pattern transpose and swap of i,j though
                Jsaall = sa.Jdydt(ys, t);
                for idx = 1:nf
                    Jsa{idx} = sparse(j,i,Jsaall(:,idx),dsa,dsa);
                end
            end
            
            %% Initialize Jacobian
            J = sparse(this.fDim,this.xDim);
            pos = 1:(dm+dsa)*nf;
            J(pos,off_mech+pos) = blkdiag(Jdm{:},Jsa{:});
            
            %% Motoneuron to Sarcomere coupling
            moto_out = y(this.moto_sarco_link_moto_out);
            dsignal_dmotoout = (this.mslinkfunderiv(moto_out).*moto_out + this.mslinkfun(moto_out))./sa.SarcoConst(1,:)';
            pos = sub2ind(size(J),this.moto_sarco_link_sarco_in,this.moto_sarco_link_moto_out);
            J(pos) = dsignal_dmotoout;
%             for k=1:nf
%                 J((k),(k)) = dsignal_dmotoout(k);
%             end
            
            %% External input signal
            % Get current external signal
            moto_signal = mo.checkMeanCurrent(this.mu(4)*this.normalized_cortex_signal(t));
            
            %% Spindle stuff
            if sys.HasSpindle
                sp = sys.Spindle;
                ds = sp.Dims;
                off_spindle = nf*(dm+dsa);
                spindle_pos = off_mech + off_spindle + (1:this.num_spindle_dof);
                yspindle = reshape(y(spindle_pos),ds,[]);
                Jsp = cell(1,nf);
                if this.fUseFD
                    freq = this.FrequencyDetector.Frequency;
                else   
                    % Get single spindle signals
                    spindle_sig = this.SpindleAffarentWeights*sp.getAfferents(yspindle);
                    % Compute the mean value over all signals
                    spindle_sig = ones(1,nf)*mean(spindle_sig);
                    % Use the upper bounded sum
                    spindle_sig = min(spindle_sig,this.max_spindle_signals - moto_signal);
                    moto_signal = moto_signal + spindle_sig;
                    freq_kexp_arg = [sys.Model.Config.FibreTypes; moto_signal];
                    freq = this.freq_kexp.evaluate(freq_kexp_arg);
                end
                K = this.System.f;
                i = []; j = []; s = [];
                moto_pos = 2:6:6*nf;
                for k=1:nf
                    spindle_pos = off_spindle + ds*(k-1) + (1:ds);
                    
                    %% Spindles by themselves
                    [Jsp{k}, Jspin_dLdot, Jspin_dmoto] = ...
                        sp.Jdydt(yspindle(:,k), t, freq(k), K.lambda_dot(k), 0);
                    
                    %% Mechanics to spindle coupling
                    J(spindle_pos,1:sys.NumStateDofs+sys.NumDerivativeDofs) = Jspin_dLdot'*K.JLamDot(k,:);
                    
                    %% Spindle to Motoneuron coupling
                    daffk_dy = this.SpindleAffarentWeights*sp.getAfferentsJacobian(yspindle(:,k));
                    %dnoise_daff = this.Noise(:,round(t)+1).*this.SomaInputFactors(:);
                    alpha = 1.24635;
                    dnoise_daff = (9^(1-alpha)*alpha*moto_signal.^(alpha-1))' .* this.Noise(:,round(t)+1);
                    i = [i repmat(moto_pos',1,9)];%#ok
                    j = [j repmat(9*(k-1) + (1:9),nf,1)];%#ok
                    s = [s dnoise_daff*daffk_dy/nf];%#ok
                    
                    %% Moto to Spindle coupling for learned frequencies
                    if ~this.fUseFD
                        kexp_Jac = this.freq_kexp.getStateJacobian(freq_kexp_arg(:,k));
                        
                        J(spindle_pos,spindle_pos) = J(spindle_pos,spindle_pos) ...
                            + Jspin_dmoto'*kexp_Jac(2)*daffk_dy;
                    end
                end
                J(1:this.num_motoneuron_dof,...
                  off_mech + off_spindle + (1:this.num_spindle_dof))...
                    = sparse(i(:),j(:),s(:),this.num_motoneuron_dof,this.num_spindle_dof);
            end
        end
        
        function value = get.UseFrequencyDetector(this)
            value = this.fUseFD;
        end
        
        function set.UseFrequencyDetector(this, value)
            if ~isequal(value,this.fUseFD) || isempty(this.fUseFD)
                this.fUseFD = value;
                % Update pattern (learned frequency causes different
                % spindle pattern)
                if ~isempty(this.System.Model.Config)
                    this.JSparsityPattern = this.computeSparsityPattern;
                end
            end
        end
        
        function assembleX0(this)
            sys = this.System;
            mc = sys.Model.Config;
            opt = mc.Options;
            
            % Load dynamic/affine x0 coefficients for moto/sarco system
            % from file
            s = load(fullfile(fileparts(which('models.motorunit.System')),...
                    sprintf('x0coeff%d.mat',opt.SarcoVersion)));
            this.x0_motorunit = dscomponents.AffineInitialValue;
            m = size(s.coeff,1);
            for k=1:m
                this.x0_motorunit.addMatrix(sprintf('polyval([%s],mu(1))',...
                    sprintf('%.14g ',s.coeff(k,:))),full(sparse(k,1,1,m,1)));
            end
            
            % Evaluate the dynamic ICs for the selected fibre types and use
            % them as x0
            ft = mc.FibreTypes;
            smoff = zeros(this.nfibres,1);
            x0_moto_sarco_spindle = zeros(sys.NumFirstOrderDofs,1);
            dm = sys.Motoneuron.Dims;
            dsa = sys.Sarcomere.Dims;
            for k=1:this.nfibres
                off = 0;
                x0ms = this.x0_motorunit.evaluate(ft(k));
                % add moto
                x0_moto_sarco_spindle(off + dm*(k-1) + (1:dm)) = x0ms(1:dm);
                % add sarco
                off = off + this.num_motoneuron_dof;
                x0_moto_sarco_spindle(off + dsa*(k-1) + (1:dsa)) = x0ms((dm+1):end);
                smoff(k) = x0ms(dm+53);
                if sys.HasSpindle
                    ds = sys.Spindle.Dims;
                    off = off + this.num_sarco_dof;
                    % add spindle
                    x0_moto_sarco_spindle(off + ds*(k-1) + (1:ds)) = sys.Spindle.y0;
                end
            end
            this.x0 = x0_moto_sarco_spindle;
            this.sarco_mech_signal_offset = smoff;
        end
        
        function updateDimensions(this, mc)
            nf = length(mc.FibreTypes);
            sys = this.System;
            
            dm = sys.Motoneuron.Dims;
            dsa = sys.Sarcomere.Dims;          
            this.num_motoneuron_dof = dm*nf;
            this.num_sarco_dof = dsa*nf;
            this.num_spindle_dof = 0;
            if sys.HasSpindle
                ds = sys.Spindle.Dims;
                this.num_spindle_dof = ds*nf;
            end
            this.fDim = (dm+dsa)*nf + this.num_spindle_dof;
            
            % This quantity is set here in updatedimensions as the
            % configUpdated method of the Dynamics class uses this value.
            % As the configUpdated for this class (i.e. the FO component)
            % is called after that (and in general there is no right order
            % to do so), it's put here.
            this.sarco_output_idx = this.num_motoneuron_dof ...
                + (53:dsa:this.num_sarco_dof);
        end
        
        function res = test_Jacobian(this, varargin)
            res = test_Jacobian@dscomponents.ACoreFun(this, varargin{:});
            res = res && models.fullmuscle.Spindle.test_Spindle_Jac;
        end
    end
end

