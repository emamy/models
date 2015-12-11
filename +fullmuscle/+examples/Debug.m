classdef Debug < models.fullmuscle.AMuscleConfig
    % A simple configuration for Debug purposes.
    % 
    % Uses a single undeformed cube with triquadratic position shape
    % functions and trilinear shape functions for pressure.
    
    methods
        function this = Debug(varargin)
            this = this@models.fullmuscle.AMuscleConfig(varargin{:});
            this.addOption('Version',1);
            
            this.init;
        end
        
        function configureModel(this, m)
            configureModel@models.fullmuscle.AMuscleConfig(this, m);
            m.ODESolver = solvers.MLode15i;
            switch this.Options.Version
            case 1
                m.T = 100;
                m.dt = .1;
                
                m.ODESolver.RelTol = .001;
                m.ODESolver.AbsTol = .01;
            end
        end
        
        function P = getBoundaryPressure(this, elemidx, faceidx)
            % Determines the neumann forces on the boundary.
            %
            % The unit for the applied quantities is megaPascal [MPa]
            %
            % In the default implementation there are no force boundary
            % conditions.
            P = [];
            if this.Options.Version == 10
%                 if any(elemidx - (9:16) == 0) && faceidx == 2
                if elemidx == 1 && faceidx == 2
                   P = 2;
                end
            end
        end
    end
    
    methods(Access=protected)
        
        function geo = getGeometry(~)
            % Single cube with same config as reference element
            geo = fem.geometry.RegularHex27Grid([0 1],[0 1],[0 1]);
        end
        
        function [ft, ftw] = getFibreInfo(this)
%             ft = [0 .2 .4 .6 .8 1];
%             ft = [0 .4 1];
            ft = 0;
            ftw = this.getZeroFTWeights(length(ft));
            % Test: Use only slow-twitch muscles
            ftw(:,1,:) = .4;
%             ftw(:,2,:) = .05;
%             ftw(:,3,:) = .05;
%             ftw(:,4,:) = .1;
%             ftw(:,5,:) = .2;
%             ftw(:,6,:) = .2;
        end
        
        function sp = getSpindlePos(~)
            % Spindle position: first row element, second row gauss point
            % within element
            sp = [1; 1];
        end
        
        function displ_dir = setPositionDirichletBC(this, displ_dir)
            geo = this.FEM.Geometry;
            % Always fix back side
%             if this.Version == 10
%                 displ_dir(1,geo.Elements(1,geo.MasterFaces(1,:))) = true;
%             else
                displ_dir(:,geo.Elements(1,geo.MasterFaces(1,:))) = true;
%                 displ_dir(:,geo.Elements(1,geo.MasterFaces(2,:))) = true;
%                 displ_dir(:,geo.Elements(1,geo.MasterFaces(3,:))) = true;
%             end
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
        
        function anull = seta0(this, anull)
            switch this.Options.Version
            case {1}
                anull(1,:,:) = 1;
%             case {4,7}
%                 % No fibres
%             case {5,8}
%                 % Stretch along fibre direction
%                 anull(1,:,:) = 1;
%             case {6,9}
%                 % Stretch perpendicular to fibre direction
%                 anull(2,:,:) = 1;
%             case {10,11}
%                 anull(:,:,:) = 0;
            end
        end
    end
    
    methods(Static)
        function test_DebugConfig(version)
            if nargin < 1
                version = 1;
            end
            m = models.fullmuscle.Model(...
                models.fullmuscle.examples.Debug('Version',version));
            [t,y] = m.simulate;
            df = m.getResidualForces(t,y);
            
            m.plot('DF',df);
        end
    end
    
end

