classdef RCLadder < models.BaseFullModel
% RCLadder: Model of a nonlinear resistor with independent current source
%
% This model has been used as benchmark in many papers. This implementation follows the details
% in \cite BS06
%
% This model is also used in several papers dealing with MOR of nonlinear systems, i.e.
% \cite RW03 or \cite CI04.
%
% @author Daniel Wirtz @date 2011-04-29
%
% @change{0,3,sa,2011-05-11} Implemented proeprty setter
%
% @new{0,3,dw,2011-04-29} Added this class.
%
% This class is part of the framework
% KerMor - Model Order Reduction using Kernels:
% - \c Homepage http://www.morepas.org/software/index.html
% - \c Documentation http://www.morepas.org/software/kermor/index.html
% - \c License @ref licensing
    
    properties
        % The target dimension of the RC ladder circuit
        %
        % @propclass{data}
        Dims;
    end
    
    methods
        
        function this = RCLadder(dims)
            
            this = this@models.BaseFullModel;
            
            if nargin == 0
                dims = 30;
            end
            this.Dims = dims;
            
            this.registerProps('Dims');
            
            this.Name = 'RC Ladder circuit model';
            
            this.T = 3;
            this.dt = .0025;
            this.tau = 1;
            this.SaveTag = sprintf('rcladder_d%d_T%g',dims,this.T);
            this.Data = data.ModelData(this);
            
            this.System = models.circ.RCLadderSys(this);
            
            % Only train with first input!
            this.DefaultInput = 1;
            this.DefaultMu = [];
            this.TrainingInputs = [2 3];
            
            this.Sampler = [];
            
            app = approx.KernelApprox(this.System);
            a = approx.algorithms.VKOGA;
            a.MaxRelErr = 1e-5;
            a.MaxAbsResidualErr = 1e-3;
            ec = kernels.config.ExpansionConfig;
            ec.StateConfig = kernels.config.GaussConfig('D',.3:.1:2);
            a.ExpConfig = ec;
            app.Algorithm = a;
                        
            t = data.selection.TimeSelector;
            t.Size = 12000;
            app.TrainDataSelector = t;
            this.Approx = app;
            
            s = spacereduction.PODReducer;
            s.Mode = 'abs';
            s.Value = 3;
            s.UseSVDS = false;
            this.SpaceReducer = s;
            
            s = solvers.SemiImplicitEuler(this);
            s.MaxStep = .05; % Stability constraint due to diffusion term
            this.ODESolver = s;
            
            this.ErrorEstimator = error.IterationCompLemmaEstimator;
        end
        
        function set.Dims(this, value)
            if ~isposintscalar(value)
                error('value must be a positive integer scalar');
            end
            this.Dims = value;
        end
    end
    
    methods(Static)
        function test_RCLadderModel
            m = models.circ.RCLadder(100);
            m.simulate([], 1);
            [t,y,~,x] = m.simulate([], 2);
            m.plot(t,y);
            m.plotState(t,x);
        end
    end
    
    methods(Static, Access=protected)
        function this = loadobj(this)
            if ~isa(this, 'models.circ.RCLadder')
                s = this;
                this = models.circ.RCLadder(s.Dims);
                % field "Dims" is set in constructor
                this = loadobj@models.BaseFullModel(this, s);
            else
                this = loadobj@models.BaseFullModel(this);
            end
        end
    end
end