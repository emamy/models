classdef Debug < models.fullmuscle.AMuscleConfig & matlab.unittest.TestCase
    % A set of simple configurations for debug purposes.
    % Select with "Version" option.
    % 
    % Versions:
    % 1-3: One element with increasing number of motor units / fibre types
    % 4: Setup without spindle
    % 5: One side is pulled
    %
    
    properties(Constant)
        RS = RandStream('mt19937ar','Seed',1);
    end
    
    methods
        function this = Debug(varargin)
            this = this@models.fullmuscle.AMuscleConfig(varargin{:});
            this.addOption('Spindle',true);
            %this.addOption('Pull',.1);
            this.addOption('ExtSig',true);
            this.addOption('NFibres',3);
            
            this.init;
        end
        
        function configureModel(this, m)
            configureModel@models.fullmuscle.AMuscleConfig(this, m);
            
            m.T = 50;
            m.dt = .1;

            m.ODESolver = solvers.MLode15i;
            m.ODESolver.RelTol = .001;
            m.ODESolver.AbsTol = .01;
            
%             o = this.Options;
        end
        
        function P = getBoundaryPressure(~, elemidx, faceidx)
            % Determines the neumann forces on the boundary.
            P = [];
            if elemidx == 1 && faceidx == 2
               P = 1;
            end
        end
    end
    
    methods(Access=protected)
        
        function geo = getGeometry(~)
            % Single cube with same config as reference element
            geo = fem.geometry.RegularHex27Grid([0 1],[0 1],[0 1]);
        end
        
        function [ft, ftw] = getFibreInfo(this)
            nf = this.Options.NFibres;
            % With this we are able to get a clean zero which wont happen
            % for "rand"
            ft = sort((this.RS.randperm(500,nf)-1)/499);
            ftw = this.getZeroFTWeights(length(ft));
            w = rand(1,nf);
            w = w./sum(w);
            for k = 1:nf
                ftw(:,k,:) = w(k);
            end
        end
        
        function sp = getSpindlePos(this)
            % Spindle position: first row element, second row gauss point
            % within element
            sp = [];
            if this.Options.Spindle
                nf = this.Options.NFibres;
                gp = this.RS.randperm(this.FEM.GaussPointsPerElem,nf);
                sp = [ones(1,nf); gp];
            end
        end
        
        function displ_dir = setPositionDirichletBC(this, displ_dir)
            geo = this.FEM.Geometry;
            displ_dir(:,geo.Elements(1,geo.MasterFaces(1,:))) = true;
        end
        
%         function [velo_dir, velo_dir_val] = setVelocityDirichletBC(this, velo_dir, velo_dir_val)
%             % Determines the dirichlet velocities.
%             %
%             % The unit for the applied quantities is [mm/ms] = [m/s]
%             geo = this.PosFE.Geometry;
%             switch this.Version
%             
%             case {4,5,6}
%                 % Move the whole front
%                 velo_dir(1,geo.Elements(1,geo.MasterFaces(2,:))) = true;
%                 velo_dir_val(velo_dir) = .05;
%             end
%         end
        
        function anull = seta0(~, anull)
            anull(1,:,:) = 1;
        end
    end
    
    properties(TestParameter)
        pull = {0 .1};
        cortex_sig = {0 3 6};
        spindle = {true false};
        nfibres = {1 4};
    end
    
    methods(Test)
        function TestFullMuscle(~, pull, cortex_sig, spindle, nfibres)
            mc = models.fullmuscle.examples.Debug(...
                'Pull',pull,'ExtSig',cortex_sig, 'Spindle', spindle,...
                'NFibres',nfibres);
            m = mc.createModel;
            [t,y] = m.simulate;
            df = m.getResidualForces(t,y);
            m.plot(t,y,'DF',df);
        end
    end
    
end

