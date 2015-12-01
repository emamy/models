clear classes;
basedir = fileparts(mfilename('fullpath'));

file = fullfile(basedir, 'paramdomaindetection_withnoise.mat');
load(file);

%% Plot raw data
pm = PlotManager;
pm.FigureSize = [1024 768];
ax = pm.nextPlot('raw_data','Frequencies by detector','Fibre type \tau','Mean current \kappa');
tri = delaunay(ps(1,:),ps(2,:));
trisurf(tri,ps(1,:),ps(2,:),Hz,'Parent',ax,'FaceColor','interp','EdgeColor','interp');

%% assemble training data
m = min(ps,[],2);
M = max(ps,[],2);
n = 200;
ft = linspace(m(1),M(1),n);
mc = linspace(m(2),M(2),n);
[FT,MC] = meshgrid(ft,mc);
HZ = griddata(ps(1,:),ps(2,:),Hz,FT,MC,'nearest');

%% smooth out data
padsize = 6;
stencil = [0 0 1 0 0; 0 1 1 1 0; 1 1 1 1 1; 0 1 1 1 0; 0 0 1 0 0];
% stencil = [1 1 1; 1 1 1; 1 1 1];
% stencil = [0 1 0; 1 1 1; 0 1 0];
stencil = stencil/sum(stencil(:));
HZ = padarray(HZ,[padsize padsize],'replicate');
runs = 10;
for k = 1:runs
    HZ = conv2(HZ,stencil,'same');
end
HZ = HZ(padsize+1:end-padsize,padsize+1:end-padsize);
ax2 = pm.nextPlot('smoothed','Smoothed frequency surface','Fibre type \tau','Mean current \kappa');
surf(ax2,FT,MC,HZ,'FaceColor','interp','EdgeColor','interp');
axis(ax2,'tight');

%% Init
s = load(models.motoneuron.Model.FILE_UPPERLIMITPOLY);
maxv = polyval(s.upperlimit_poly,FT);
valid = MC <= maxv;
atd = data.ApproxTrainData([FT(valid) MC(valid)]',[],[]);
atd.fxi = HZ(valid)';

%% init algorithm
alg = approx.algorithms.VKOGA;
ec = kernels.config.ExpansionConfig;
comb = Utils.createCombinations([1 2 3 4],[1 2 3]);
sc = kernels.config.WendlandConfig('G',comb(1,:),'S',comb(2,:));
ec.StateConfig = sc;
alg.ExpConfig = ec;
alg.MaxExpansionSize = 1000;

%% compute
kexp = alg.computeApproximation(atd);

%% plot
afx = kexp.evaluate([FT(:) MC(:)]');
FX = reshape(afx,200,[]);
FX(~valid) = Inf;
ax3 = pm.nextPlot('approx','Learned smoothed frequency surface','Fibre type \tau','Mean current \kappa');
surf(ax3,FT,MC,FX,'FaceColor','interp','EdgeColor','interp');
axis(ax3,'tight');

%% Save stuff
pm.UseFileTypeFolders = false;
pm.ExportDPI = 150;
pm.FilePrefix = 'freq';
pm.NoTitlesOnSave = true;
zlim(ax2,zlim(ax));
zlim(ax3,zlim(ax));
pm.savePlots(basedir,'Format','jpg','Close',true);
file = fullfile(basedir, 'FreqencyLearning.mat');
save(file, 'alg', 'atd', 'HZ', 'FT', 'MC', 'valid');
file = fullfile(basedir, 'FrequencyKexp.mat');
save(file, 'kexp');