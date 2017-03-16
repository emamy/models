classdef Processor < handle
    %PROCESSOR Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        minV = -20;
    end
    
    methods
        
        function [T,V,invalid] = getVelocities(this, t, Vms, distances)
            nruns = length(Vms);
            if isscalar(distances)
                distances = ones(1,nruns)*distances;
            end
            V = cell(1,nruns);
            T = cell(1,nruns);
            
            pi = ProcessIndicator('Extracting velocities for %d runs',nruns,false,nruns);
            invalid = [];
            for idx = 1:nruns
                [peaks_junction, peaks_end] = this.getPeakLocations(t, Vms{idx});
                
                % Get the time instances the peak was recorded
                peaktimes = [];
                peaktimes(1,:) = t(peaks_junction);
                peaktimes(2,:) = t(peaks_end);
                % Compute velocities
                tdiff = diff(peaktimes,1);
                V{idx} = 10*distances(idx)./tdiff; % [cm/ms = 0.01m/0.001s = .1m/s]
                T{idx} = peaktimes(1,:);
                if length(T{idx}) < 2
                    invalid = [invalid idx];%#ok
                end
                pi.step;
            end
            pi.stop;
            % Remove invalid runs
            V(invalid) = [];
            T(invalid) = [];
        end
        
        function [T, P, invalid] = getPeakAmplitudes(this, t, Vms)
            nruns = length(Vms);
            P = cell(1,nruns);
            T = P;
            pi = ProcessIndicator('Extracting peak amplitudes from %d runs',nruns,false,nruns);
            invalid = [];
            for idx = 1:nruns
                Vm = Vms{idx}(1,:);
                peaks = this.getPeakIdx(t,Vm);
                if length(peaks) > 2
                    % Get the time instances the peak was recorded
                    T{idx} = t(peaks);
                    % get the Amplitudes
                    P{idx} = Vm(peaks);
                    % Interpolate data
%                     c = polyfit(peaktimes,peakvals,3);
%                     P(idx,:) = polyval(c,tgrid);
%                     plot(peaktimes,peakvals,'rx',tgrid,P(idx,:),'b');
                else
                    invalid = [invalid idx];%#ok
                end
                pi.step;
            end
            pi.stop;
            % Remove invalid runs
            P(invalid) = [];
            T(invalid) = [];
        end
        
        function [Vinterp, Vpoly] = getInterpolatedVelocities(~, T, V, tgrid, degree)
            if nargin < 5
                degree = 6;
            end
            %% Compute interpolated data
            nruns = length(T);
            tsteps = length(tgrid);
            Vinterp = zeros(nruns,tsteps);
            Vpoly = zeros(nruns,tsteps);
            for idx = 1:nruns
                Vi = interp1(T{idx},V{idx},tgrid,'pchip');
                Vinterp(idx,:) = Vi;
                c = polyfit(T{idx},V{idx},degree);
                Vpoly(idx,:) = polyval(c,tgrid);
            end
        end

        function [peaks_junction, peaks_end] = getPeakLocations(this, t, Vm)
            peaks_junction = this.getPeakIdx(t, Vm(1,:));
            peaks_end = this.getPeakIdx(t, Vm(2,:));

            % Trivial case
            ndiff = length(peaks_junction) - length(peaks_end);
            if ndiff > 0
                plot(t,Vm(1,:),'r',t(peaks_junction),Vm(1,peaks_junction),'bx',...
                    t,Vm(2,:),'b',t(peaks_end),Vm(2,peaks_end),'rx')
                peaks_junction(end-ndiff+1:end) = [];
                fprintf(2,'Check plot for correct removal of last junction peak!\n');
            end
        end
        
        function peaks = getPeakIdx(this, t, signal)
            % Find all locations with negative derivative
            negpos = find(diff(signal(1,:))<0);
            % Find the positions where the negative derivative starts
            firstpos = [1 find(diff(negpos) > 1)+1];
            % Check that the found locations are peaks and not in the lower noise
            ispeak = signal(1,negpos(firstpos)) > this.minV;
            % Remove unwanted
            firstpos(~ispeak) = [];
            peaks = negpos(firstpos);
            % Check for close consecutive peaks & remove
            tdiff = diff(t(peaks));
            closepeaks = tdiff < median(tdiff)/10;
            peaks(closepeaks) = [];
        end
    end
    
end

