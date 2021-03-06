classdef BaseCoreFun < dscomponents.ACompEvalCoreFun
% BaseCoreFun: 
%
% @docupdate
%
% @author Daniel Wirtz @date 2013-03-14
%
% @new{0,7,dw,2013-03-14} Added this class.
%
% This class is part of the framework
% KerMor - Model Order Reduction using Kernels:
% - \c Homepage http://www.morepas.org/software/index.html
% - \c Documentation http://www.morepas.org/software/kermor/index.html
% - \c License @ref licensing
    
    properties(Constant)
        % The time in seconds needed before the activation function reaches its maximum value
        % and decays back to zero.
        %
        % Only applies for ActivationFunType == 2
        %
        % @default 30 @type double
        %
        % See also: MaxActivationTime ActivationFunType
        ActivationTransitionTime = 10; %[s]
        
        % The maximum time in seconds that spans the support of the piecewise activation
        % function.
        % It is composed of gaussian-shaped increase/decrease ends and a constant one in
        % between. For the default settings, we would have a maximum duration of
        % 500s-2*30s=440s for the activation rate of level one.
        %
        % The model parameter "atime" determines the level-one-activation time linearly between
        % [2*ActivationTransitionTime, MaxActivationTime]
        %
        % The activation function always starts at t=0.
        %
        % Only applies for ActivationFunType == 2
        %
        % @default 500 @type double
        %
        % See also: ActivationTransitionTime ActivationFunType
        MaxActivationTime = 400; %[s]
    end

    properties(Dependent)
        % Type of the activation function.
        %
        % Included for backward compatibility with older simulations.
        %
        % Admissible values:
        % 1: First version with gaussian activation over 300s
        % 2: Parameterized version
        %
        % @default 1 @type integer
        ActivationFunType;
    end
    
    properties(Access=protected)
        nodes;
        
        hlp;
        
        idxmat;
        
        Ji;
        Jj;
    end
    
    properties(Access=private)
        gaussian;
        fAFT = 1;
        fTau;
    end
    
    methods
        function this = BaseCoreFun(dynsys)
            this = this@dscomponents.ACompEvalCoreFun(dynsys);
            this.TimeDependent = true;
            
            this.fTau = dynsys.Model.tau;
            
            this.ActivationFunType = 2;
        end
        
        function copy = clone(this, copy)
            % Call superclass method
            copy = clone@dscomponents.ACompEvalCoreFun(this, copy);
            copy.gaussian = this.gaussian.clone;
            copy.fTau = this.fTau;
            copy.fAFT = this.fAFT;
        end
        
        function fx = evaluateCoreFun(this, x, t)
            error('This model overrides evaluate directly.');
        end
        
        function plotActivationFun(this, mu, pm)
            if nargin < 3
                pm = PlotManager;
                pm.LeaveOpen = true;
                if nargin < 2
                    mu = this.System.Model.getRandomParam;
                    
                end
            end
            
            h = pm.nextPlot('activation_fun','External input activation function','time','factor');
            t = linspace(0,min(this.System.Model.T,this.MaxActivationTime*1.2),2000);
            plot(h,t,this.activationFun(t/this.System.Model.tau,mu));
            
            if nargin < 3
                pm.done;
            end
        end
        
        function set.ActivationFunType(this, value)
            % Activation function setup
            if value == 1
                % \gamma=28 is such that we have K(0,27)~<.4 (27=150/tau)
                k = kernels.GaussKernel(28.206364723698);
            else
                k = kernels.GaussKernel;
                k.setGammaForDistance(this.ActivationTransitionTime/this.System.Model.tau,.001);
            end
            this.gaussian = k;
            this.fAFT = value;
        end
        
        function value = get.ActivationFunType(this)
            value = this.fAFT;
        end
    end
    
    methods(Access=protected)
        
        function idx = nodepos(this, nr)
            n = this.nodes;
            idx = zeros(1,n*length(nr));
            for k=1:length(nr)
                idx((k-1)*n+1:k*n) = (nr(k)-1)*this.nodes+1:nr(k)*this.nodes;
            end
        end
        
        function f = activationFun(this, t, mu)
            if this.fAFT == 1
                f = (this.gaussian.evaluateScalar(t-27)-.4).*(t<=54)/.6;
            else
                tau = this.fTau;
                ts = this.ActivationTransitionTime/tau;
                if isscalar(t)
                    te = ts+(mu(3)*(this.MaxActivationTime-2*this.ActivationTransitionTime))/tau;
                    if t > te+ts
                        f = 0;
                    elseif t > ts && t <= te
                        f = 1;
                    elseif t <= ts
                        f = (this.gaussian.evaluateScalar(t-ts)-.001)/.999;
                    elseif t>te && t<=te+ts
                        f = (this.gaussian.evaluateScalar(t-te)-.001)/.999;
                    else 
                        f = 0;
                    end
                else
                    te = ts+(mu(3,:)*(this.MaxActivationTime-2*this.ActivationTransitionTime))/tau;
                    f = (this.gaussian.evaluateScalar(t-ts)-.001).*(t<=ts)/.999 ...
                        + (t>ts).*(t<=te) ...
                        + (this.gaussian.evaluateScalar(t-te)-.001).*(t>te).*(t<=te+ts)/.999;
                end
            end
        end
    end
    
end