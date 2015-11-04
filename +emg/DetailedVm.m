classdef DetailedVm < handle
    %VMGENERATOR A separate class used to compute detailed Vm signals for
    % defined parameters, to be used within the models.emg suite
    %
    
    properties(SetAccess=private)
        submesh_idx;
        % musclefibre model for full simulations
        musclefibremodel;
        mus_precomp;
        s;
    end
    
    methods
        
        function this = DetailedVm(len, dims, sv)
            dx = models.musclefibre.Model.dxDefault;
            N = round(len/dx);
            % Coarse index within fine resolution
            this.submesh_idx = round(linspace(1,N,dims));
            mf = models.musclefibre.Model(...
                'SarcoVersion',sv,...
                'N',N,'dx',dx,...
                'DynamicIC',true,'SPM',false,'OutputScaling',false,...
                'Spindle',false,'Noise',true,...
                'JunctionN',1);
            mf.EnableTrajectoryCaching = true;
            this.musclefibremodel = mf;
            
            % Hack for 
            mudir = fileparts(which('models.emg.Model'));
            s = load(fullfile(mudir,'data','mus.mat'));
            this.mus_precomp = s.mus;
            
            this.s = load('/data/local/musclefibre/precomp_Vms.mat');
        end
        
        function sig = computeSignal(this, t, mu)
            pos = Utils.findVecInMatrix(this.mus_precomp,mu);
            if ~isempty(this.s)
                if (t(2)-t(1) ~= this.s.t(2)-this.s.t(1))
                    error('Incompatible time-steps. Please check.');
                end
                sig = this.s.Vms{pos}(:,1:length(t));
            else
                if pos > 0
                    fprintf('Loading cached file ... ');
                    s = load(sprintf('/data/local/musclefibre/Vm_%d.mat',pos));
                    fprintf('done!\n');
                    fine_signal = s.Vm(:,1:2:2*length(t));
                else
                    mf = this.musclefibremodel;
                    mf.T = t(end); % infer from passed parameters
                    mf.dt = t(2)-t(1);
                    [~,~,~,x] = mf.simulate(mu,1);
                    % Find positions of Vm on each sarcomere
                    moto_off = mf.System.dm;
                    % Each first sarcomere model dimension contains the current
                    % Vm value
                    pos = moto_off + (1:mf.System.dsa:size(x,1)-moto_off);
                    fine_signal = x(pos,:);
                end
                % Only return the effectively needed signals from the
                % finer submodel simulation.
                sig = fine_signal(this.submesh_idx,:);
            end
        end
    end
end

