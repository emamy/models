classdef Feedback < models.fullmuscle.AMuscleConfig & models.muscle.IExperimentModelConfig
    % Config for feedback experiments
    %
    
    properties(Constant)
        GeoLength = 30; % [mm]
        StretchTime = 100; % [ms]
    end
    
    properties
        StretchDist;
    end

    methods
        function this = Feedback(varargin)
            % Creates a Debug simple muscle model configuration.
            %
            % Single cube with same config as reference element
            this = this@models.fullmuscle.AMuscleConfig(varargin{:});
            this.NeumannCoordinateSystem = 'global';
            % No external signals - whatsoever
            this.NormalizedCortexSignal = @(t)zeros(size(t));
            this.addOption('Motorunits',4);
            this.addOption('SpindleGP',9);
            this.addOption('NumConf',4);
            
            this.init;
            this.initExperiment;
            o = this.Options;
            gl = this.GeoLength;
            % Set stretch distances so that lambdamax stretch is reached
            % for last config
            lambdamax = 1.4;
            this.StretchDist = linspace(0,lambdamax-1,o.NumConf)*gl; % [cm]
            
            this.VelocityBCTimeFun = ...
                general.functions.ConstantUntil(this.StretchTime);
            
            % IExperimentModelConfig values
            this.NumConfigurations = o.NumConf;
            this.NumOutputs = 1;
        end
        
        function configureModel(this, m)
            configureModel@models.fullmuscle.AMuscleConfig(this, m);
            
            m.T = 50;
            m.dt = .1;
                
            m.DefaultMu(1:4) = [1; 0; 0; 0];
            m.DefaultMu(13) = .250; %[MPa]
            m.DefaultInput = [];
            
            %this.ActivationRampMax = 1;
            
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
        
        function o = getOutputOfInterest(this, t, y)
            o = 1;
        end
        
        function x0 = getX0(this, x0)
            x0 = getX0@models.muscle.IExperimentModelConfig(this, x0);
        end
        
    end
    
    methods(Access=protected)
        
        function geo = getGeometry(~)
            %o = this.Options;
            geo = fem.geometry.RegularHex27Grid([0 3],[0 2],[0 2]);
        end
        
        function [ft, ftw] = getFibreInfo(this)
            nm = this.Options.Motorunits;
            ft = linspace(0,1,nm);       
            ftw = this.getZeroFTWeights(1);
            ftw_fc = general.functions.ExpDist(1,20);
            ftw_fh = ftw_fc.getFunction;
            w = ftw_fh(ft);
            w = w/sum(w);
            for k = 1:nm
                ftw(:,k,:) = w(k);
            end
        end
        
        function sp = getSpindlePos(this)
            % Put all spindles at same position, gauss point given by
            % option "SpindleGP"
            o = this.Options;
            sp = repmat([1; o.SpindleGP],1,o.Motorunits);
        end
        
        function displ_dir = setPositionDirichletBC(this, displ_dir)
            geo = this.FEM.Geometry;
            % Fix back side
            displ_dir(:,geo.Elements(1,geo.MasterFaces(1,:))) = true;
        end
        
        function [velo_dir, velo_dir_val] = setVelocityDirichletBC(this, velo_dir, velo_dir_val)
            % Determines the dirichlet velocities.
            %
            % The unit for the applied quantities is [mm/ms] = [m/s]
            geo = this.FEM.Geometry;
            
            % Move the whole front
            velo_dir(1,geo.Elements(1,geo.MasterFaces(2,:))) = true;
            v = this.StretchDist(this.CurrentConfigNr)/this.StretchTime;
            velo_dir_val(velo_dir) = v;
        end
        
        function anull = seta0(~, anull)
            anull(1,:,:) = 1;
        end
    end
    
end

