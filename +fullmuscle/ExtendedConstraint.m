classdef ExtendedConstraint < models.muscle.Constraint
    
    properties(Access=private)
    
    end
    
    properties(SetAccess=private)
    
    end
    
    methods
        function this = ExtendedConstraint(sys)
            this = this@models.muscle.Constraint(sys);
        end
        
%         function configUpdated(this)
%             configUpdated@models.muscle.Constraint(this);
%         end
       
        
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

        function dy = evaluate(this, y, t)
            sys = this.System;
            dy = zeros(this.fDim,1);
            
            %% Mechanics
            uvp_pos = 1:sys.NumTotalDofs;
            % Use uvp as argument and also pass in s (=sarco forces)
            force = max(0,y(sys.sarco_output_idx)-sys.sarco_mech_signal_offset);
%             force = y(sys.sarco_output_idx);
%             force = zeros(length(sys.sarco_output_idx),1);
            uvps = [y(uvp_pos); force];
            dy(uvp_pos) = evaluate@models.muscle.Dynamics(this, uvps, t);
            
            %% Motoneurons
            mo = sys.Motoneuron;
            moto_pos = sys.off_moto+(1:sys.num_motoneuron_dof);
            dy_m = mo.dydt(reshape(y(moto_pos),6,[]),t);
            dy(moto_pos) = dy_m(:);
            
            %% Sacromeres
            sa = sys.Sarcomere;
            sarco_pos = sys.off_sarco + (1:sys.num_sarco_dof);
            ys = reshape(y(sarco_pos),56,[]);
            dys = sa.dydt(ys, t);
            dy(sarco_pos) = dys(:);
            
            %% Link of motoneurons to sarcomeres
            moto_out = y(this.moto_sarco_link_moto_out);
            fac = min(this.MSLink_MaxFactor,this.MSLinkFun(moto_out));
            signal = fac.*moto_out./this.sarcoconst1';
            % Add signal to corresponding locations
            dy(this.moto_sarco_link_sarco_in) = dy(this.moto_sarco_link_sarco_in) + signal;
            
            if sys.HasSpindle
                %% Spindles
                sp = sys.Spindle;
                spindle_pos = sys.off_spindle + (1:sys.num_spindle_dof);
                yspindle = reshape(y(spindle_pos),9,[]);

                % Link of spindle to motoneuron
                % Get single spindle signals
                spindle_sig = this.SpindleAffarentWeights*sp.getAfferents(yspindle);
                % Compute the mean value over all signals
                spindle_sig = ones(1,nf)*mean(spindle_sig);
                % Get current external signal
                ext_sig = sys.Inputs{2,1}(t);
                % Use the upper bounded sum
                eff_spindle_sig = min(spindle_sig,this.max_moto_signals - ext_sig);
                % Compute noisy signal
                noise_sig = mo.TypeNoise(:,round(t)+1)'.*eff_spindle_sig.*mo.FibreTypeNoiseFactors;
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
                dys = sp.dydt(yspindle,t,freq,this.lambda_dot,0);
                dy(spindle_pos) = dys(:);
            end
        end
        
        function SP = computeSparsityPattern(this)
            sys = this.System;
            
            % Get original constraint pattern
            SP = computeSparsityPattern@models.muscle.Constraint(this);
            mech_off = size(SP,1);
            % append zeros for extra dofs
            SP = [SP sparse(mech_off,sys.num_extra_dof)];
            
            [~, SPalpha, SPLamDot] = this.System.f.computeSparsityPattern;
            
            nf = sys.nfibres;
            
            % Neuro
            dm = sys.Motoneuron.Dims;
            for k=1:nf
                pos = dm*(k-1)+1:dm;
                SP(mech_off+pos,sys.off_moto + pos) = sys.Motoneuron.JSparsityPattern;%#ok
            end
            
            % Sarco
            J_sarco = sys.Sarcomere.JSparsityPattern;
            dsa = sys.Sarcomere.Dims;
            for k=1:nf
                pos = dm*nf + dsa*(k-1) + 1:dsa;
                SP(mech_off+pos,sys.off_sarco + pos) = J_sarco;%#ok
            end
            
            if sys.HasSpindle
                % Spindle
                JSp = sys.Spindle.JSparsityPattern;
                for k=1:nf
                    SP = blkdiag(SP, JSp);
                end
            end
            
            % Moto -> Sarco link
            for k=1:nf
                % first entry of sarco gets 2nd output of motoneuron
                off_sarco = sys.off_sarco + (k-1)*56 + 1;
                off_moto = sys.off_moto + (k-1)*6 + 2;
                SP(off_sarco,off_moto) = true;
            end
                    
            % Sarco -> Mechanics
            SP(sys.NumStateDofs + (1:sys.NumDerivativeDofs), sys.off_sarco+(1:sys.num_sarco_dof)) = SPalpha;
            
            if sys.HasSpindle
                % Spindle -> Motoneuron link
                i = []; j = [];
                moto_pos = 2:6:6*nf;
                for k=1:nf
                    i = [i repmat(moto_pos',1,9)];%#ok
                    j = [j repmat(9*(k-1) + (1:9),nf,1)];%#ok
                end
                SP(sys.off_moto + (1:sys.num_motoneuron_dof),...
                  sys.off_spindle + (1:sys.num_spindle_dof))...
                    = sparse(i(:),j(:),ones(numel(i),1),...
                    sys.num_motoneuron_dof,sys.num_spindle_dof);
                
                % Mechanics -> spindle
                Jspin_Ldot = double(sys.Spindle.JLdotSparsityPattern);
                Jspin_dmoto = double(sys.Spindle.JMotoSparsityPattern);
                Jspin_Aff = double(sys.Spindle.JAfferentSparsityPattern);
                for k=1:nf
                    spindle_pos = sys.off_spindle + 9*(k-1) + (1:9);
                    SP(spindle_pos,1:sys.NumStateDofs+sys.NumDerivativeDofs) = ...
                        logical(Jspin_Ldot*double(SPLamDot(k,:)));
                    % Create connecting link to self only when no frequency
                    % detector is used!
                    if ~this.fUseFD
                        SP(spindle_pos,spindle_pos) = ...
                            SP(spindle_pos,spindle_pos) | logical(Jspin_dmoto*any(Jspin_Aff));
                    end
                end
            end
        end
        
        function [J, JLamDot] = getStateJacobian(this, y, t)
%             J = this.getStateJacobianFD(y,t);
%             return;
            sys = this.System;
            
            %% Mechanics
            uvp_pos = 1:sys.NumTotalDofs;
            force = max(0,y(sys.sarco_output_idx)-sys.sarco_mech_signal_offset);
%             force = y(sys.sarco_output_idx);
%             force = zeros(length(sys.sarco_output_idx),1);
            uvps = [y(uvp_pos); force];
            [J, Jalpha, JLamDot] = getStateJacobian@models.muscle.Dynamics(this, uvps, t);
            
            %% Motoneuron
            mo = sys.Motoneuron;
            for k=1:nf
                moto_pos = sys.off_moto + 6*(k-1) + (1:6);
                J = blkdiag(J,mo.Jdydt(y(moto_pos),t,k));
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
    end
    
    methods(Access=private)
       
    end
    
end

