classdef Model < models.muscle.Model
% Model: 
%
% Features:
% - regionen mit summe der gewichtungen kleiner 1 ist region mit fett?!
% Scenarios:
% - breiter block an einer seite fest mit allen fasertypen über die länge drin & entsprechendem
% muskelverhalten
% - block mit links fest, rechts per neumann dran ziehen, dann
% spindel-feedback bei zu starker dehnung sorgt für konktraktion
% - dynamische aktivierung von aussen
% - dynamische aktivierung von aussen mit dynamischem zug
% - zyklischer zug
%
% Fragen:
% - Sollte der afferent-faktor zur signalübertragung für alle fasertypen gleich sein?
% - Lddot versuch? was wäre "realistisch"?
%
% @author Daniel Wirtz @date 2014-09-16
%
% @new{0,7,dw,2014-09-16} Added this class.
%
% This class is part of the framework
% KerMor - Model Order Reduction using Kernels:
% - \c Homepage http://www.agh.ians.uni-stuttgart.de/research/software/kermor.html
% - \c Documentation http://www.agh.ians.uni-stuttgart.de/documentation/kermor/
% - \c License @ref licensing

    properties(Constant)
        DataDir = fullfile(fileparts(mfilename('fullpath')),'data');
    end
    
    methods
        function this = Model(conf)
            if nargin < 1
                conf = models.fullmuscle.examples.Debug;
            end
            this = this@models.muscle.Model(conf);
        end
        
        function setConfig(this, value)
            if ~isa(value, 'models.fullmuscle.AMuscleConfig')
                error('Config must be a models.fullmuscle.AMuscleConfig instance');
            end
            setConfig@models.muscle.Model(this, value);
            this.Plotter = models.fullmuscle.MusclePlotter(this.System);
        end
        
        function varargout = plotGeometrySetup(this, varargin)
            varargin = [varargin {'GeoOnly',true}];
            [varargout{1:nargout}] = plotGeometrySetup@models.muscle.Model(this,varargin{:});
        end
        
    end
    
    methods(Access=protected)
        function subclassInit(this, config)
            
            this.System = models.fullmuscle.System(this, config);
            
            this.DefaultInput = 1;
            
            % Der tut auch wunderbar :-)
            slv = solvers.MLWrapper(@ode15s);
%             slv.odeopts = odeset('RelTol',1e-9,'AbsTol',1e-9);
            
%             slv = solvers.MLode15i(this);
%             slv.RelTol = .001;
%             slv.AbsTol = .01;
            this.ODESolver = slv;
        end
    end
    
end