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
        
        % The upper limit for the mean input current fed to the motoneuron
        % soma. as in this current version this is the sum of spindle
        % feedback and external signal, the max value needs to be available
        % here to limit the sum.
        % The external signal is added at an higher level in the ode (B*u
        % component), but yet the sum of (spindle+ext_sig) <
        % max_moto_signal, which is why the external signal is accessed
        % here, too.
        max_moto_signals;
        
        freq_kexp;
        fUseFD = false;
        nfibres;
        
        x0_motorunit;
    end
    
    properties(SetAccess=private)
        num_motoneuron_dof;
        num_sarco_dof;
        num_spindle_dof;
%         off_moto;
%         off_sarco;
%         off_spindle;
        
        input_motoneuron_link_idx;
        moto_input_noise_factors;
        sarco_output_idx;
        
        moto_sarco_link_moto_out;
        moto_sarco_link_sarco_in;
        spindle_moto_link_moto_in;
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
        end
        
        function configUpdated(this)
            sys = this.System;
            mc = sys.Model.Config;
            if ~isempty(mc)
                this.nfibres = length(mc.FibreTypes);
                this.xDim = sys.NumTotalDofs;
                this.assembleX0;
                this.JSparsityPattern = this.computeSparsityPattern;
            end
        end
       
        
%         function projected = project(this, V, W)
%             projected = this.clone;
%             projected = project@dscomponents.ACoreFun(this, V, W, projected);
%         end
%         
%         function copy = clone(this)
%             % Create new instance
%             copy = models.muscle.Constraint(this.System);
%             copy = clone@dscomponents.ACoreFun(this, copy);
%             copy.ComputeUnassembled = this.ComputeUnassembled;
%             copy.fsys = this.fsys;
%             copy.idx_p_elems_unass = this.idx_p_elems_unass;
%             copy.Sigma = this.Sigma;
%             copy.fDim_unass = this.fDim_unass;
%         end

        function prepareSimulation(this, mu)
            prepareSimulation@dscomponents.ACoreFun(this, mu);
            [this.mslinkfun, this.mslinkfunderiv] = this.MSLink.getFunction;
            
            this.max_moto_signals = this.System.Motoneuron.getMaxMeanCurrents;
            
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
            ymoto = y(off_mech + (1:this.num_motoneuron_dof));
            mo = sys.Motoneuron;
            dy_m = mo.dydt(reshape(ymoto,mo.Dims,[]),t);
            dy(1:this.num_motoneuron_dof) = dy_m(:);
            
            %% Sacromeres
            sa = sys.Sarcomere;
            sarco_pos = this.num_motoneuron_dof + (1:this.num_sarco_dof);
            ys = reshape(y(off_mech + sarco_pos),sa.Dims,[]);
            dys = sa.dydt(ys, t);
            dy(sarco_pos) = dys(:);
            
            %% Link of motoneurons to sarcomeres
            moto_out = y(off_mech + this.moto_sarco_link_moto_out);
            signal = this.mslinkfun(moto_out).*moto_out./sa.SarcoConst(1,:)';
            % Add signal to corresponding locations
            dy(this.moto_sarco_link_sarco_in) = ...
                dy(this.moto_sarco_link_sarco_in) + signal;
            
            if sys.HasSpindle
                %% Spindles
                sp = sys.Spindle;
                spindle_pos = this.num_motoneuron_dof ...
                    + this.num_sarco_dof + (1:this.num_spindle_dof);
                yspindle = reshape(y(off_mech + spindle_pos),sp.Dims,[]);

                %% Spindles -> Motoneurons
                % Get single spindle signals
                spindle_sig = this.SpindleAffarentWeights*sp.getAfferents(yspindle);
                % Compute the mean value over all signals
                spindle_sig = ones(1,nf)*mean(spindle_sig);
                % Get current external signal
                ext_sig = sys.Inputs{2,1}(t);
                % Use the upper bounded sum
                eff_spindle_sig = min(spindle_sig,this.max_moto_signals - ext_sig);
                % Compute noisy signal
                noise_sig = mo.Noise(:,round(t)+1)'.*eff_spindle_sig.*mo.FibreTypeNoiseFactors;
    %             fprintf('Spindle->Neuron: adding %g at dy(%d)\n',noise_sig,this.spindle_moto_link_moto_in);
                dy(this.spindle_moto_link_moto_in) = ...
                    dy(this.spindle_moto_link_moto_in) + noise_sig';

                % Spindle actual
                % Get motoneuron frequency
                if this.fUseFD
                    freq = this.FrequencyDetector.Frequency;
                else
                    freq = this.freq_kexp.evaluate([sys.Model.Config.FibreTypes; eff_spindle_sig+ext_sig]);
                end
                dys = sp.dydt(yspindle,t,freq,sys.f.lambda_dot,0);
                dy(spindle_pos) = dys(:);
            end
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
                pos = dm*(k-1)+1:dm;
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
            SP(pos_sarco,pos_moto) = true;            
            
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
        
        function [J, JLamDot] = getStateJacobian(this, y, t)
%             J = this.getStateJacobianFD(y,t);
%             return;
            sys = this.System;
            
            % Get original constraint jacobian
            J = getStateJacobian@models.muscle.Constraint(this, y(1:sys.off_moto), t);
            mech_pressure_off = size(J,1);
            % append zeros behind for extra dofs
            J = [J sparse(mech_pressure_off,sys.num_extra_dof)];
            
            nf = sys.nfibres;
            
            %% Motoneuron
            mo = sys.Motoneuron;
            dm = mo.Dims;
            local_off_sarco = mech_pressure_off + nf*dm;
            for k=1:nf
                pos = dm*(k-1)+1:dm;
                J(mech_pressure_off+pos,sys.off_moto + pos) = ...
                    mo.Jdydt(y(moto_pos),t,k);
            end
            
            %% Sarcomeres
            sa = sys.Sarcomere;
            for k=1:nf
                sarco_pos = sys.off_sarco + 56*(k-1) + (1:56);
                J = blkdiag(J,sa.Jdydt(y(sarco_pos),t,k));
            end
            
            %% Motoneuron to Sarcomere coupling
            moto_out = y(this.moto_sarco_link_moto_out);
            fac = min(this.MSLink_MaxFactor,this.MSLinkFun(moto_out));
            dfac = this.MSLinkFunDeriv(moto_out);
            dsignal_dmotoout = (dfac .* moto_out + fac)./this.sarcoconst1';
            for k=1:nf
                J(this.moto_sarco_link_sarco_in(k),this.moto_sarco_link_moto_out(k)) = dsignal_dmotoout(k);
            end
            
            %% Sarcomere to mechanics coupling
            % The JS matrix is generated during the computation of the
            % mechanics jacobian, as the element/gauss loop is computed
            % there anyways. its not 100% clean OOP, but works for now.
            J(sys.NumStateDofs + (1:sys.NumDerivativeDofs), sys.off_sarco+(1:sys.num_sarco_dof)) = Jalpha;
            
            %% Spindle stuff
            if sys.HasSpindle
                sp = sys.Spindle;
                if this.fUseFD
                    freq = this.FrequencyDetector.Frequency;
                else
                    % For no detection, the current spindle signal is
                    % required
                    spindle_pos = sys.off_spindle + (1:sys.num_spindle_dof);
                    yspindle = reshape(y(spindle_pos),9,[]);
                    % Get single spindle signals
                    spindle_sig = this.SpindleAffarentWeights*sp.getAfferents(yspindle);
                    % Compute the mean value over all signals
                    spindle_sig = ones(1,nf)*mean(spindle_sig);
                    % Use the upper bounded sum
                    eff_spindle_sig = min(this.max_moto_signals, spindle_sig + sys.Inputs{2,1}(t));
                    freq_kexp_arg = [sys.Model.Config.FibreTypes; eff_spindle_sig];
                    freq = this.freq_kexp.evaluate(freq_kexp_arg);
                end
                
                i = []; j = []; s = [];
                moto_pos = 2:6:6*nf;
                for k=1:nf
                    spindle_pos = sys.off_spindle + 9*(k-1) + (1:9);
                    
                    %% Spindles by themselves
                    [Jspin, Jspin_dLdot, Jspin_dmoto] = sp.Jdydt(y(spindle_pos), t, freq(k), this.lambda_dot(k), 0);
                    J = blkdiag(J,Jspin);
                    
                    %% Mechanics to spindle coupling
                    J(spindle_pos,1:sys.NumStateDofs+sys.NumDerivativeDofs) = Jspin_dLdot'*JLamDot(k,:);
                    
                    %% Spindle to Motoneuron coupling
                    daffk_dy = this.SpindleAffarentWeights*sp.getAfferentsJacobian(y(spindle_pos));
                    dnoise_daff = mo.TypeNoise(:,round(t)+1).*mo.FibreTypeNoiseFactors(:);
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
                J(sys.off_moto + (1:sys.num_motoneuron_dof),...
                  sys.off_spindle + (1:sys.num_spindle_dof))...
                    = sparse(i(:),j(:),s(:),sys.num_motoneuron_dof,sys.num_spindle_dof);
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
            off = 0;
            for k=1:this.nfibres
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
            %this.nf = nf
            
            dm = sys.Motoneuron.Dims;
            dsa = sys.Sarcomere.Dims;
%              % Motoneurons are beginning after mechanics
%             % The current NumXXDofs are computed by the mechanics-only
%             % model, hence it's the correct offset.
%             this.off_moto = this.NumStateDofs+this.NumDerivativeDofs;
%             
%             % Sarcomeres are beginning after motoneurons
%             this.off_sarco = this.off_moto + fo.num_motoneuron_dof;
%             
%             % Spindles are beginning after sarcomeres
%             this.off_spindle = this.off_sarco + fo.num_sarco_dof;
            this.num_motoneuron_dof = dm*nf;
            this.num_sarco_dof = dsa*nf;
            this.num_spindle_dof = 0;
            if sys.HasSpindle
                ds = sys.Spindle.Dims;
                this.num_spindle_dof = ds*nf;
            end
            this.fDim = (dm+dsa)*nf + this.num_spindle_dof;
            
            %% Component link indices
            % Get the positions where the input signal is mapped to the
            % motoneurons
            this.input_motoneuron_link_idx = 2:dm:dm*nf;
            this.sarco_output_idx = this.num_motoneuron_dof + (53:dsa:dsa*nf);
            
            % same but improves readability
            this.moto_sarco_link_moto_out = 2:dm:dm*nf; 
            this.moto_sarco_link_sarco_in = this.num_motoneuron_dof + (1:dsa:dsa*nf);
            this.spindle_moto_link_moto_in = this.moto_sarco_link_moto_out;
        end
    end
    
    methods
        
    end
    
end

