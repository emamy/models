classdef MSLink < general.functions.AFunGen
    % Class to model the functional linking of the motoneuron model to the
    % sarcomere model
    
    properties
        % The `V_s` value of the motoneuron input at which the MSLink_MaxFactor should be attained
        %
        % @type double @default 40
        MSLink_MaxFactorSignal = 40;
        
        % The maximal factor with which the `V_s` value of the motoneuron should be amplified
        % when added to the sarcomere equations
        %
        % @type double @default 7
        MSLink_MaxFactor = 7;
        
        % The minimal factor at which the `V_s` value of the motoneuron should be amplified
        % when added to the sarcomere equations.
        %
        % Between both limits, an exponential Gaussian weight is applied with a certain radius,
        % so that the amplification factor has its minimum value at zero and maximum value at
        % MSLink_MaxFactorSignal.
        %
        % See also FibreDynamics.getLinkFactor
        %
        % @type double @default .3
        MSLink_MinFactor = .3;
    end
    
    methods
        
        function this = MSLink
            this.xLabel = 'Motoneuron signal [mV]';
            this.yLabel = 'Amplification factor [-]';
        end
        
        function [fhandle, dfhandle] = getFunction(this)
            maxf = this.MSLink_MaxFactor;
            minf = this.MSLink_MinFactor;
            diff = this.MSLink_MaxFactor-this.MSLink_MinFactor;
            ms = this.MSLink_MaxFactorSignal;
            fhandle = @(t)(t < ms).*(minf + diff*exp(-(t-ms).^2/150))...
                +(ms <= t)*maxf;
            dfhandle = @(t)(t < ms).*-exp(-(t-ms).^2/150).*(diff*(t-ms))/75;
        end
        
        function str = getConfigStr(this)
            str = sprintf('Min/max factors=%g/%g, max signal=%g',...
                this.MSLink_MinFactor,this.MSLink_MaxFactor,this.MSLink_MaxFactorSignal);
        end
        
        function pm = plot(this, varargin)
            varargin(end+1:end+2) = {'R' 0:.1:80};
            pm = plot@general.functions.AFunGen(this,varargin{:});
        end

    end
    
end