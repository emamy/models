%%
clear classes;%#ok
%%
m = models.motoneuron.Model;
m.T = 3000;

%% Single simulation for test
[t,signal] = m.simulate;
peaks = p.getPeakIdx(t, signal(2,:));
plot(t,signal(2,:),'b',t(peaks),signal(2,peaks),'rx');

%% Sampling - entire domain
%s = sampling.GridSampler;
s = sampling.RandomSampler;
s.Samples = 2400;
s.Domain = models.motoneuron.ParamDomain;
mus = s.generateSamples(m);
%% Sampling - original setting
mus = [linspace(0,1,12); ones(1,12)];
%% Sampling - rough coverage of domain
mus = [0 0 .3 .3 .5 .5 .7 .7 .9 .9 .9 1
       1 3 1  3  1  3  1  3  1  3  5 7];

%% Compute
nmu = size(mus,2);
signals = zeros(nmu,length(m.Times));
PCPool.open;
parfor pidx = 1:nmu
    [t,y] = m.simulate(mus(:,pidx),1);%#ok
    signals(pidx,:) = y(2,:);
end
PCPool.close;

%% or load precomputed
% setup with fine grid of 314 mus over domain, t=400ms only
% load /home/dwirtz/papers/MWHR15' - EMG Signal generation'/matlab/data314mus_orig.mat
% setup with 12 fibre types for mc=3, T=3000 ms
% load /home/dwirtz/papers/MWHR15' - EMG Signal generation'/matlab/data12mus_mc3.mat

%%
nmu = size(mus,2);
p = SignalAnalyzer;
p.minV = 55;
cov = NaN(1,nmu);
np = cov;
t = m.Times;
for pidx = 1:nmu
    signal = signals(pidx,:);
    peaks = p.getPeakIdx(t, signal);
    if ~isempty(peaks)
        peaktimes = t(peaks);
        isi = diff(peaktimes)/1000;
        %dr = 1./(diff(peaktimes)/1000);
        cov(pidx) = std(isi)/mean(isi);
        np(pidx) = length(peaks);
    end
end
%%
%plot3(mus(1,:),mus(2,:),cov,'rx');
% col = np/max(np);
col = 1-cov/max(cov);
pm = PlotManager;
ax = pm.nextPlot('cov','Coefficients of variation for inter-spike intervals over parameter range','fibre type','mean current');
scatter3(ax,mus(1,:),mus(2,:),cov,10*ones(size(cov)),[1-col; col; zeros(size(col))]','filled');
hold(ax,'on');
[FT,MC] = meshgrid([0 1],[0 9]);
sh = surf(FT,MC,.15*ones(size(FT)),'Parent',ax,'FaceColor',[0 0 .5]);
alpha(sh,.3);
sh = surf(FT,MC,.20*ones(size(FT)),'Parent',ax,'FaceColor',[0 0 .5]);
alpha(sh,.3);
pm.done;