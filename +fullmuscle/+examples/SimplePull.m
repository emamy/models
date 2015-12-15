classdef SimplePull < models.fullmuscle.AMuscleConfig
    % A full muscle example featuring a pulling on one side
    %
    % There are 7 inputs of various characteristics for testing and playing
    % around. Feel free to add you own in "getInputs" and invoke with
    % <myinstance>.simulate(mu,<myinputidx>)
    %
    % Versions (Option "Version"): 
    % 1:
    % 2:
    % 3:

    methods
        function this = SimplePull(varargin)
            % Creates a Debug simple muscle model configuration.
            %
            % Single cube with same config as reference element
            this = this@models.fullmuscle.AMuscleConfig(varargin{:});
            this.NeumannCoordinateSystem = 'global';
            % No external signals - whatsoever
            this.NormalizedCortexSignal = @(t)zeros(size(t));
            this.addOption('Version',1);
            this.addOption('X',[0 10]);
            this.addOption('Y',[0 10]);
            this.init;
        end
        
        function configureModel(this, m)
            configureModel@models.fullmuscle.AMuscleConfig(this, m);
            
            switch this.Options.Version
                case {1,2}
                    m.T = 100;
                    m.dt = .1;
                case 3
                    m.T = 1000;
                    m.dt = 1;
            end
            
            m.DefaultMu(1:4) = [1; 0; .15; 0];
            m.DefaultMu(13) = .250; %[MPa]
            m.EnableTrajectoryCaching = true;
        end
        
        function P = getBoundaryPressure(this, elemidx, faceidx)
            % Determines the neumann forces on the boundary.
            %
            % The unit for the applied quantities is megaPascal [MPa]
            %
            % In the default implementation there are no force boundary
            % conditions.
            P = [];
            elem = 1;
            if this.Options.Version == 3
                elem = 2;
            end
            if elemidx == elem && faceidx == 2
                P = 1;
            end
        end
        
        function u = getInputs(this)
            u{1} = this.getAlphaRamp(10,1);
            u{2} = this.getAlphaRamp(1,1);
            u{3} = this.getAlphaRamp(300,1);
            u{4} = this.getAlphaRamp(10,1,200);
            u{5} = this.getAlphaRamp(100,1,200);
            u{6} = this.getAlphaRamp(300,1,200);
            u{7} = this.getAlphaRamp(1000,1);
        end
        
    end
    
    methods(Access=protected)
        
        function geo = getGeometry(this)
            o = this.Options;
            geo = fem.geometry.RegularHex27Grid(o.X,o.Y,[0 10]);
        end
        
        function [ft, ftw] = getFibreInfo(this)
            switch this.Options.Version
                case 1
                   ft = 0;
                   ftw = this.getZeroFTWeights(1);
                   ftw(:,1,:) = 1;
                case 2
                   ft = [0 1];
                   ftw = this.getZeroFTWeights(2);
                   ftw(:,1,:) = .5;
                   ftw(:,2,:) = .5;
                case 3
                   ft = [0 .2 .4 .6 .8 1];
                   ftw = this.getZeroFTWeights(length(ft));
                   fac = exp((1:6)/2);
                   fac = fac / sum(fac);
                   ftw(:,1,:) = fac(1);
                   ftw(:,2,:) = fac(2);
                   ftw(:,3,:) = fac(3);
                   ftw(:,4,:) = fac(4);
                   ftw(:,5,:) = fac(5);
                   ftw(:,6,:) = fac(6);
            end
        end
        
        function sp = getSpindlePos(this)
            % Spindle position: first row element, second row gauss point
            % within element
            switch this.Options.Version
                case 1
                   sp = [1 9]';
                case 2
                   sp = [1 1
                         1 2];
                case 3
                   sp = [1 1 1 2 2 2
                         1 10 20 4 10 25]; 
            end
        end
        
        function displ_dir = setPositionDirichletBC(this, displ_dir)
            geo = this.FEM.Geometry;
            % Always fix back side
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
    
end

