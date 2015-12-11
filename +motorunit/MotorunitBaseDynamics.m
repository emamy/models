classdef MotorunitBaseDynamics < dscomponents.ACoreFun
% MotorunitBaseDynamics: 
%
% @docupdate
%
% @author Daniel Wirtz @date 2015-09-15
%
% @new{0,7,dw,2015-09-15} Added this class.
%
% This class is part of the framework
% KerMor - Model Order Reduction using Kernels:
% - \c Homepage http://www.morepas.org/software/
% - \c Documentation http://www.morepas.org/software/kermor/
% - \c License @ref licensing
    
    properties
        MSLink;
    end
    
    properties(SetAccess=private,GetAccess=protected)
        mslinkfun;
    end
    
    methods
        function this = MotorunitBaseDynamics(sys)
            this = this@dscomponents.ACoreFun(sys);
            this.MSLink = models.motorunit.MSLink;
        end
        
        function prepareSimulation(this, mu)
            prepareSimulation@dscomponents.ACoreFun(this, mu);
            this.mslinkfun = this.MSLink.getFunction;
        end
        
        function plotMotoSacroLinkFactorCurve(this)
            this.MSLink.plot;
        end
    end
    
end