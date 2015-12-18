classdef MusclePlotter < models.muscle.MusclePlotter
% Extends the models.muscle.MusclePlotter functionality
%
% @docupdate
%
% @author Daniel Wirtz @date 2014-09-25
%
% @new{0,7,dw,2014-09-25} Added this class.
%
% This class is part of the framework
% KerMor - Model Order Reduction using Kernels:
% - \c Homepage http://www.agh.ians.uni-stuttgart.de/research/software/kermor.html
% - \c Documentation http://www.agh.ians.uni-stuttgart.de/documentation/kermor/
% - \c License @ref licensing
    
    methods
        
        function this = MusclePlotter(sys)
            this = this@models.muscle.MusclePlotter(sys);
        end
        
        function [pm, h_geo] = plot(this, t, y, varargin)
            opts = this.parsePlotArgs(varargin);
            
            mc = this.Config;
            sys = this.System;
            sel = opts.MU;
            fo = sys.FO;
            
            if isempty(opts.PM)
                if opts.GeoOnly
                    pm = PlotManager;
                else
                    pm = PlotManager(false,3,3);
                end
                pm.LeaveOpen = true;
            else
                pm = opts.PM;
            end
            
            [pd, t, y] = this.updatePlotData(this.plotdata, opts, t, y);
            this.plotdata = pd;
            
            if opts.Vid
                vw = VideoWriter(opts.Vid);
                vw.FrameRate = 25;
                vw.open;
            end
            
            if opts.Geo
                pmgeo = PlotManager;
                pmgeo.LeaveOpen = true;
%                 pmgeo = pm;
                h_geo = pmgeo.nextPlot('geo',sprintf('Deformation at t=%g',t(end)),'x [mm]','y [mm]');
                zlabel(h_geo,'z [mm]');
                axis(h_geo, this.getPlotBox(y));
                daspect([1 1 1]);
                view(h_geo, this.GeoView);
                hold(h_geo,'on');
            end
            
            if ~opts.GeoOnly
                ms = 10;
                
                if opts.Moto
                    h1 = pm.nextPlot('signal','Motoneuron signal','t [ms]','V_m');
                    plot(h1,t,pd.moto_vm(sel,:));
                    axis(h1,[0 t(end) min(pd.moto_vm(:)) max(pd.moto_vm(:))]);
                    hold(h1,'on');
                    ph1 = [];
                end

                if opts.Sarco
                    h2 = pm.nextPlot('force','Action potential','t [ms]','V_m');
                    plot(h2,t,pd.sarco_pot(sel,:));
                    axis(h2,[0 t(end) min(pd.sarco_pot(:)) max(pd.sarco_pot(:))]);
                    hold(h2,'on');
                    ph2 = [];

                    h3 = pm.nextPlot('force','Motorunit force signals','t [ms]','A_s');
                    force = pd.sarco_force(sel,:);
                    plot(h3,t,force,'r');
                    axis(h3,[0 t(end) min(force(:)) max(force(:))]);
                    hold(h3,'on');
                    ph3 = [];
                    
                    h3b = pm.nextPlot('force','Weighted Activation at Elem1,GP1','t [ms]','A_s');
                    force = mc.FibreTypeWeights(1,:,1)*pd.sarco_force;
                    plot(h3b,t,force,'b');
                    axis(h3b,[0 t(end) min(force(:)) max(force(:))]);
                    hold(h3b,'on');
                    ph3b = [];
                end

                if opts.Freq
                    hfreq = pm.nextPlot('frequency','Motoneuron frequency','t [ms]','Frequency [Hz]');
                    m = min(pd.freq(:));
                    M = max(pd.freq(:));
                    plot(hfreq,t,pd.freq(sel,:))
                    hold(hfreq,'on');
                    if opts.FreqDet && ~fo.UseFrequencyDetector
                        m = min(m, min(pd.freq_det(:)));
                        M = max(M, max(pd.freq_det(:)));
                        plot(hfreq,t,pd.freq_det(sel,:),'r--');
                    end
                    axis(hfreq,[0 t(end) m M]);
                    phfreq = [];
                end

                spos = mc.SpindlePositions;
                if opts.Spin
                    h_spin_l = pm.nextPlot('spindle_lambda',...
                        'Spindle lambda stretch','t [ms]','Lambda [L_0 = 1]');
                    plot(h_spin_l, t, pd.spindle_lambda(sel,:));
                    axis(h_spin_l,[0 t(end) min(pd.spindle_lambda(:)) max(pd.spindle_lambda(:))]);
                    hold(h_spin_l,'on');
                    phs = [];

                    if opts.Aff
                        h4 = pm.nextPlot('spindle','Afferents','t [ms]','aff');
                        affsel = [2*(sel-1)+1; 2*(sel-1)+2];        
                        plot(h4,t,pd.afferents(affsel(:),:)');
                        axis(h4,[0 t(end) min(pd.afferents(:)) max(pd.afferents(:))]);
                        hold(h4,'on');
                        ph4 = [];
                    end
                end
                
                h5 = pm.nextPlot('mean_current','Motoneuron input mean current','t [ms]','mean current');
                plot(h5,t,pd.ext_mean_current,'b--');
                hold(h5,'on');
                plot(h5,t,pd.eff_mean_current(sel,:),'r');
                if sys.HasSpindle
                    plot(h5,t,pd.spindle_mean_current,'g--');
                end
                axis(h5,[0 t(end) min(pd.eff_mean_current(:)) max(pd.eff_mean_current(:))+100*eps]);
                ph5 = [];

                if opts.Ext && ~isempty(sys.inputidx)
                    hext = pm.nextPlot('ext_neumann',...
                        'External pressure','t [ms]','Normal pressure [MPa]');
                    plot(hext, t, pd.uneum);
                    axis(hext,[0 t(end) min(pd.uneum(:)) max(pd.uneum(:))+eps]);
                    hold(hext,'on');
                    phext = [];
    %                 val = pd.uneum * pd.forcefactor;
    %                 hext_f = pm.nextPlot('ext_neumann',...
    %                     'External force','t [ms]','Normal force [mN]');
    %                 axis(hext_f,[0 t(end) min(val(:)) max(val(:))]);
    %                 hold(hext_f,'on');
                end
                
            end
                
            pm.done;
            
            fh = gcf;
            for ts = 1:length(t)
                % Quit if figure has been closed
                if ~ishandle(fh) || ~ishandle(h_geo)
                    break;
                end
                
                if opts.Geo
                    yf = pd.yfull(:,ts);
                    this.plotGeometry(h_geo, t(ts), yf, ts, opts);
                    if ~opts.GeoOnly && opts.Spin && ~isempty(spos)
                        % Plot spindle locations
                        for k = sel
                            u = yf(sys.idx_u_elems_local(:,:,spos(1,k)));
                            spindle_pt = u*pd.Ngp(:,spos(2,k));
                            plot3(h_geo,spindle_pt(1),spindle_pt(2),...
                                spindle_pt(3),'.',...
                                'MarkerSize',20,'Color',[.3 1 .3]);
                        end 
                    end
                end
                
                if ~opts.GeoOnly
                
                    if opts.Moto
                        delete(ph1);
                        ph1 = plot(h1,t(ts),pd.moto_vm(sel,ts),'r.',...
                            'MarkerSize',ms);
                    end

                    if opts.Sarco
                        delete(ph2);
                        ph2 = plot(h2,t(ts),pd.sarco_pot(sel,ts), ...
                            'r.','MarkerSize',ms);

                        delete(ph3);
                        ph3 = plot(h3,t(ts),pd.sarco_force(sel,ts),...
                            'r.','MarkerSize',ms);
                        
                        delete(ph3b);
                        walpha = mc.FibreTypeWeights(1,:,1) * pd.sarco_force(:,ts);
                        ph3b = plot(h3b,t(ts),walpha,'r.','MarkerSize',ms);
                    end

                    delete(ph5);
                    %t(ts),pd.ext_mean_current(ts),'b--',...
                    ph5 = plot(h5,t(ts),pd.eff_mean_current(sel,ts),'r.','MarkerSize',ms);
                    
                    if opts.Spin
                        delete(phs);
                        phs = plot(h_spin_l, t(ts), pd.spindle_lambda(sel,ts),...
                            'r.','MarkerSize',ms);
                        
%                         % Also add the spindle mean current
%                         plot(h5,time_part,pd.spindle_mean_current(1:ts),'g--');

                        if opts.Aff
                            delete(ph4);
                            ph4 = plot(h4,t(ts),pd.afferents(affsel(:),ts)',...
                                'r.','MarkerSize',ms);
                            
    %                         cla(h6);
    %                         plot(h6,time_part,pd.spindle_single_mean_current(sel,1:ts));
    %                         plot(h6,time_part,pd.eff_mean_current(sel,1:ts),'r--');
                        end
                    end
%                     axis(h5,'tight');

                    if opts.Freq
                        delete(phfreq);
                        phfreq = plot(hfreq,t(ts),pd.freq(sel,ts),...
                            'r.','MarkerSize',ms);
%                         if opts.FreqDet && ~fo.UseFrequencyDetector
%                             plot(hfreq,time_part,pd.freq_det(sel,1:ts),'r--');
%                         end
                    end

                    if opts.Ext && ~isempty(sys.inputidx)
                        delete(phext);
                        phext = plot(hext, t(ts), pd.uneum(ts), ...
                            'r.','MarkerSize',ms);
    %                     cla(hext_f);
    %                     plot(hext_f, time_part, pd.forcefactor*pd.uneum(1:ts));
                    end
                
                end
                
                if opts.Vid
                    vw.writeVideo(getframe(gcf));
                else
                    if opts.Pause == 0
                        drawnow;
                    else
                        pause(opts.Pause);
                    end
                end
            end
            
            if opts.Vid
                vw.close;
            end

            if isempty(opts.PM)
                pm.done;
            end
        end
    end
    
    methods(Access=protected)
        function [pd, t, y] = updatePlotData(this, pd, opts, t, y)
            mc = this.Config;
            sys = this.System;
            nf = length(mc.FibreTypes);
            fo = sys.FO;
            
            if (opts.Freq && fo.UseFrequencyDetector) || opts.FreqDet
                fd = fo.FrequencyDetector;
                nt = length(t);
                freq = zeros(nf,nt);
                fd.reset;
                fd.WindowSize = 4;
                for i=1:nt
                    fd.processSignal(t(i),y(fo.moto_sarco_link_moto_out,i)');
                    freq(:,i) = fd.Frequency;
                end
                if ~isempty(opts.F)
                    freq = freq(:,1:opts.F:end);
                end
                pd.freq = freq;
                if opts.FreqDet
                    pd.freq_det = freq;
                end
            end
            
            [pd, t, y] = updatePlotData@models.muscle.MusclePlotter(this, pd, opts, t, y);
            nt = length(t);
            
            pd.ext_mean_current = sys.mu(4)*fo.normalized_cortex_signal(t);
            
            if opts.Moto
                pos = sys.EndSecondOrderDofs + (2:6:6*nf);
                pd.moto_vm = y(pos,:);
            end
            
            if opts.Sarco
                off_sarco = sys.EndSecondOrderDofs + fo.num_motoneuron_dof;
                pos = off_sarco + (1:56:fo.num_sarco_dof);
                pd.sarco_pot = y(pos,:);
                
                pos = off_sarco + (53:56:fo.num_sarco_dof);
                force = bsxfun(@plus, -fo.sarco_mech_signal_offset, y(pos,:));
                force = bsxfun(@times,mc.forces_scaling,force);
                pd.sarco_force = force;
            end
            
            max_moto_signals = sys.Motoneuron.getMaxMeanCurrents(mc.FibreTypes);
            eff_mean_current = zeros(nf,nt);
            if opts.Spin
                off_spindle = sys.EndSecondOrderDofs ...
                    + fo.num_motoneuron_dof + fo.num_sarco_dof;
                pos = off_spindle + (9:9:fo.num_spindle_dof);
                pd.spindle_lambda = y(pos,:);
            
                % Freq also uses afferents if kernel expansions are used
                if opts.Aff || (opts.Freq && ~fo.UseFrequencyDetector)
                    afferents = zeros(2*nf,nt);
                    spindle_single_mean_current = zeros(nf,nt);
                    
                    for k=1:nf
                        spindle_pos = off_spindle + (k-1)*9 + (1:9);
                        af_pos = (k-1)*2 + (1:2);
                        yspindle = y(spindle_pos,:);
                        afferents(af_pos,:) = sys.Spindle.getAfferents(yspindle);
                        spindle_single_mean_current(k,:) = fo.SpindleAffarentWeights*afferents(af_pos,:);
                    end
                    spindle_mean_current = mean(spindle_single_mean_current,1);
                    for k=1:nf
                        eff_mean_current(k,:) = min(max_moto_signals(k),spindle_mean_current+pd.ext_mean_current);
                    end
                    pd.afferents = afferents;
                    pd.spindle_single_mean_current = spindle_single_mean_current;
                    pd.spindle_mean_current = spindle_mean_current;
                end
            else
                for k=1:nf
                    eff_mean_current(k,:) = min(max_moto_signals(k),pd.ext_mean_current);
                end
            end
            pd.eff_mean_current = eff_mean_current;
            
            if opts.Freq && ~fo.UseFrequencyDetector
                freq = zeros(nf,nt);
                for k=1:nf
                    x = [ones(size(pd.eff_mean_current(k,:)))*mc.FibreTypes(k);...
                        pd.eff_mean_current(k,:)];
                    freq(k,:) = fo.freq_kexp.evaluate(x);
%                     freq = reshape(freq,nf,[]);
                end
                pd.freq = freq;
            end
            
            % External neumann forces
            if opts.Ext && ~isempty(sys.inputidx)
                pd.uneum = sys.Inputs{1,sys.inputidx}(t)*sys.mu(3);
%                 pd.forcefactor = norm(pd.meanforces(:,1));
            end
        end
        
        function opts = parsePlotArgs(this, args)
            %       Freq  FreqDet Aff   Geo   Moto Sarco
%             defs = [false false false false true true false];
            defs = [true true true true true true true true false];
            opts = parsePlotArgs@models.muscle.MusclePlotter(this, args);
            i = inputParser;
            i.KeepUnmatched = true;
            i.addParamValue('Freq',defs(1),@(v)islogical(v));
            i.addParamValue('FreqDet',defs(2),@(v)islogical(v));
            i.addParamValue('Aff',defs(3),@(v)islogical(v));
            i.addParamValue('Geo',defs(4),@(v)islogical(v));
            i.addParamValue('Moto',defs(5),@(v)islogical(v));
            i.addParamValue('Sarco',defs(6),@(v)islogical(v));
            i.addParamValue('Spin',defs(7),@(v)islogical(v));
            i.addParamValue('Ext',defs(8),@(v)islogical(v));
            i.addParamValue('GeoOnly',defs(9),@(v)islogical(v));
            i.addParamValue('MU',1:length(this.Config.FibreTypes));
            i.parse(args{:});
            opts = Utils.copyStructFields(i.Results, opts);
            if ~this.System.HasSpindle
                opts.Spin = false;
            end
        end
    end
    
end