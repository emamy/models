classdef Model < models.BaseFullModel
% MotoModel: Motoneuron model
%
% The global time unit for this model is milliseconds [ms].
%
% Relevante Frequenzbereiche: 6-50Hz.
%
% @todo place citations
%
% @todo Motoneuron untersuchen und parameterbereich auf sinnvolle gr��en einschr�nken! dazu eigenes
% moto-model bauen und parameterbereich abtasten, dann parameterdom�nen-geometrie einbauen.
%
% @author Daniel Wirtz @date 2012-11-22
%
% @new{1,0,dw,2014-09-17} Added the upper limit polynomial for mean_current
% dependent on fibre type.
%
% @new{0,7,dw,2013-07-05} Added this class.
%

    properties(Constant)
        FILE_UPPERLIMITPOLY = fullfile(fileparts(mfilename('fullpath')),'data','upperlimitpoly.mat');
    end

    properties(Dependent)
        UseNoise;
    end
    
    properties(SetAccess=private)
        % Flag that determines if the mean current parameter is constrained
        % depending on the fibre type parameter.
        %
        % The actually applicable parameter domain of fibre_type and
        % mean_current is not a square if only "useful" frequencies in the
        % range of 10-60 Hz are required.
        %
        % An experimentally determined polynomial is used to limit the
        % mean_current.
        %
        % @type logical @default true
        %
        % See also: models.motoneuron.experiments.ParamDomainDetection
        FibreTypeDepMaxMeanCurrent = true;
    end
    
    methods
        function this = Model(limit_meancurrent)
            if nargin < 1
                limit_meancurrent = true;
            end
            this.FibreTypeDepMaxMeanCurrent = limit_meancurrent;
            % Creates a new motoneuron model
            %
            this.Name = 'Motoneuron model';
            this.System = models.motoneuron.System(this);
            this.TrainingInputs = 1;
            this.SaveTag = 'motoneuron';
            this.Data = data.ModelData(this);
            
            this.T = 150; % [ms]
            this.dt = .1; % [ms]
            
            this.ODESolver = solvers.MLWrapper(@ode15s);
            
            this.DefaultMu = [.1; 3];
            this.DefaultInput = 1;
        end
        
        function varargout = plot(this, varargin)
            % plots some interesting states of the model
            %
            % See also: musclefibres.MuscleFibreSystem
            [varargout{1:nargout}] = this.System.plot(varargin{:});
        end
        
        function [peaks, cov_isi, cov_dr, freq] = analyze(~, t, y)
            p = SignalAnalyzer;
            p.minV = 55;
            [peaks, cov_isi, cov_dr, freq] = p.analyze(t, y(2,:));
        end
    end
    
    methods
        function value = get.UseNoise(this)
            value = ~this.System.noiseGen.DisableNoise;
        end
        
        function set.UseNoise(this, value)
            this.System.noiseGen.DisableNoise = ~value;
        end
    end
    
    methods(Static)
        function res = test_Motoneuron
            m = models.motoneuron.Model;
            [t,y] = m.simulate;
            m.plot(t,y);
            res = true;
        end
    end
    
    methods(Static, Access=protected)
        function this = loadobj(this)
            if ~isa(this, 'models.motoneuron.Model')
                sobj = this;
                this = models.motoneuron.Model;
                this = loadobj@models.BaseFullModel(this, sobj);
            else
                this = loadobj@models.BaseFullModel(this);
            end
        end
    end
    
end