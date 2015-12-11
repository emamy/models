classdef Dynamics < models.motorunit.MotorunitBaseDynamics
    % Dynamics: Class for nonlinear dynamics of motorunit model
    %
    % @author Daniel Wirtz @date 2012-11-22
    %
    % @new{0,7,dw,2012-11-22} Added this class.
    %
    % This class is part of the framework
    % KerMor - Model Order Reduction using Kernels:
    % - \c Homepage http://www.agh.ians.uni-stuttgart.de/research/software/kermor.html
    % - \c Documentation http://www.agh.ians.uni-stuttgart.de/documentation/kermor/
    % - \c License @ref licensing
    
    properties
        % Stimulation frequency in kHz, duration in ms, ampiltude in µA/cm²
        % and start time for an external stimulus bypassing the
        % motoneuron. Stimulation pattern 
        %
        % @type double @default 0.05
        ExtFreq = 0.05;
        
        % @type double @default 0.5
        ExtDur = 0.5;
        
        % @type double @default 150
        ExtAmp = 150;
        
        % @type double @default 200
        ExtStart = 200;
                
        % @type double @default Inf
        ExtBurkeFreq = Inf;
        
        % @type double @default Inf
        ExtBurkeWidth = Inf;
        
        % @type double @default 0
        ExtAllenStart = 0;
        
        % Input mode flag.
        %
        % Turn off to use the "extInput" function instead of motoneuron
        % signal.
        %
        % @type logical @default true
        LinkSarcoMoto = true;
    end
    
    properties(Access=private)
        spm;
    end
    
    properties(Transient)
        hadpeak = false;
        signalweight;
    end
    
    methods
        function this = Dynamics(sys)
            this = this@models.motorunit.MotorunitBaseDynamics(sys);
            this.spm = sys.SinglePeakMode;
            this.fDim = sys.dsa + sys.dm;
            this.xDim = this.fDim;
            this.initJSparsityPattern;
        end
        
        function prepareSimulation(this, mu)
            prepareSimulation@models.motorunit.MotorunitBaseDynamics(this, mu);
            this.hadpeak = false;
        end
        
        function extIn = extInput(this, t)
            % For a bypass of the motoneuron an external stiumlus is
            % created here. Rectangular pulses of the frequency 'freq', the
            % duration 'dur' and the amplitude 'amp' start at the time
            % 'start'. Single twitches are realized by freqency = 0.
                                  
            freq = this.ExtFreq;
            dur = this.ExtDur;
            amp = this.ExtAmp;
            start = this.ExtStart; 
            BurkeFreq = this.ExtBurkeFreq;
            BurkeWidth = this.ExtBurkeWidth;
            AllenStart = this.ExtAllenStart;
            
            if freq == 0
                if t >=start && t<=start+dur
                    extIn = amp;
                else
                    extIn = 0;
                end
                
            else
                if t<= start
                    extIn  = 0;
                else
                    if t<= AllenStart
                        if (mod(t,1/freq) <= dur)
                            extIn  = amp;
                        else
                            extIn  = 0;
                        end
                    else
                        if (mod(t,1/freq) <= dur) && (mod(t-start,1/BurkeFreq) <= BurkeWidth)
                            extIn  = amp;
                        else
                            extIn  = 0;
                        end
                    end
                end
            end
        end
        
        function dy = evaluate(this, y, t)
            % Evaluates the nonlinear core function at given time and
            % state. Can evaluate vectorized arguments. Here, each column
            % represents one state.
            % The parts for motoneuron, spindle and sarcomeres are evaluated seperately in
            % respective rate functions and the rates are concatenated in this function. The sarcomere part is reshaped such
            % that a matrix with 58 rows and one column for each single sarcomere is passed to the function
            % for (single) sarcomere dynamics.
            %
            % Parameters:
            % y: @type matrix<double>, each column is one state vector
            % t: current time @type rowvec<double>
            % mu: fibre type parameter @type rowvec<double> 0 = slow twitch fibre, 1 = fast twitch fibre
            
            sys = this.System;
            dm = sys.dm;
            
            %% Inner dynamics (with partial nonlinear linking)
            dy = zeros(size(y));
            
            %% Neuro dynamics (last 2 arguments are 2 phases: phase_prim, phase_sec)
            dy(1:dm) = sys.moto.dydt(y(1:dm),t);
            
            %% Sarcomere dynamics
            sa = sys.sarco;
            dy(dm+1:end) = sa.dydt(y(dm+1:end),t);
            
            %% Link of motoneuron to sarcomer cell
            if sys.SinglePeakMode && sys.HadPeak
                signal = 0;
            elseif this.LinkSarcoMoto
                % "Link" at the only sarcomere
                signal = this.mslinkfun(y(2))*y(2)/sa.SarcoConst(1);
            else
                % Bypass of the motoneuron: external input for sarcomere
                signal = this.extInput(t);
            end
            linkpos = sys.MotoSarcoLinkIndex;
            dy(linkpos,:) = dy(linkpos,:) + signal;
        end
        
        function dy = evaluateCoreFun(this, y, t, mu)%#ok
            error('evaluate is overridden directly.');
        end
        
        function J = getStateJacobian(this, y, t)
            % Evaluates the full Jacobian of the nonlinear core function.
            %
            % Parameters:
            % y: @type matrix<double>, each column is one state vector
            % t: current time @type rowvec<double>
            % mu: fibre type parameter, 0 = slow twitch fibre, 1 = fast twitch fibre, @type rowvec<double>
            
%             J = this.getStateJacobianFD(y,t);
%             return;
            
            s = this.System;
            dm = s.dm;
            dsa = s.dsa;
            bs = dm + dsa;
            J = spalloc(bs,bs,45+303);  % 45: 40 entries in motoneuron and spindle part + reserve, 303: entries per single sarcomere
            
            %% Neuron jacobian
            J(1:dm,1:dm) = s.moto.Jdydt(y(1:dm),t,1);
            
            %% Sarcomere jacobian
            sa = s.sarco;
            J(dm+1:bs,dm+1:bs) = sa.Jdydt(y(dm+1:end), t);
            
            if this.LinkSarcoMoto
                % link position is dm+1
                J(s.MotoSarcoLinkIndex,2) = this.mslinkfun(y(2))/sa.SarcoConst(1);
            end
        end
        
        function copy = clone(this)
            % Create new instance
            copy = models.motorunit.Dynamics(this.System);
            
            % Call superclass clone (for deep copy)
            copy = clone@dscomponents.ACompEvalCoreFun(this, copy);
            
            % Copy local properties
            % No local properties are to be copied here, as so far
            % everything is done in the constructor.
        end
        
        function initJSparsityPattern(this)
            % Initializes the Sparsity pattern of the Jacobian
            %
            sys = this.System;
            dm = sys.dm;
            dsa = sys.dsa;
            size = dm + dsa;
            J = sparse(size,size);
            
            % Neuro
            J(1:dm,1:dm) = sys.moto.JSparsityPattern;
            
            % Sarco
            J(dm+1:size,dm+1:size) = sys.sarco.JSparsityPattern;
            
            % MS link
            if this.LinkSarcoMoto
                J(sys.MotoSarcoLinkIndex,2) = true;
            end
            this.JSparsityPattern = logical(J);
        end
    end
end
