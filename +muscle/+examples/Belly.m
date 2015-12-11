classdef Belly < models.muscle.AMuscleConfig
    
    properties(Access=private)
        NumParts;
    end
    
    methods
        function this = Belly(varargin)
            this = this@models.muscle.AMuscleConfig(varargin{:});
            this.init;
        end
        
        function configureModel(this, m)
            configureModel@models.muscle.AMuscleConfig(this, m);
            m.T = 150;
            m.dt = .1;
            m.DefaultMu(4) = 6;
            m.Plotter.DefaultArgs = {'Pool',true};
        end
    end
    
    methods(Access=protected)
        
        function geo = getGeometry(this)
            np = 4;
            geo = fem.geometry.Belly(linspace(0,10,np+1),'Radius',1,'InnerRadius',.5,'Gamma',2);
            this.NumParts = np;
        end
        
        function [types, ftw] = getFibreInfo(this)
            %% Muscle fibre weights
            types = [0 .2 .4 .6 .8 1];
            ftw = this.getZeroFTWeights(length(types));
            % Test: Use only slow-twitch muscles
            ftw(:,1,:) = .4;
            ftw(:,2,:) = .05;
            ftw(:,3,:) = .05;
            ftw(:,4,:) = .1;
            ftw(:,5,:) = .2;
            ftw(:,6,:) = .2;
        end
        
        function displ_dir = setPositionDirichletBC(this, displ_dir)
            %% Dirichlet conditions: Position (fix one side)
            geo = this.FEM.Geometry;
            for k = geo.NumElements-3:geo.NumElements
                displ_dir(:,geo.Elements(k,geo.MasterFaces(4,:))) = true;
            end
            for k = 1:4
                displ_dir(:,geo.Elements(k,geo.MasterFaces(3,:))) = true;
            end
        end
        
        function anull = seta0(~, anull)
            anull(2,:,:) = 1;
        end
    end
    
    methods(Static)
        function test_BellyModel
            m = models.muscle.Model(models.muscle.examples.Belly);
            m.dt = min(m.T/10,2);
            m.simulateAndPlot;
        end
    end
    
end

