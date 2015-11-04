% Script to pre-compute several action potential shapes for faster use
% within the emg model. This is optional, as the correct shapes for the
% selected motor-units can be computed directly for the fibre types used
% in each emg model.
clear classes;%#ok
% Take 100 linearly spaced fibre types for sure
numsamples = 100;
fibretypes = linspace(0,1,numsamples);

% Also add the fibre types from "the" global parameter set (faster
% tests/experiments as precomp can be used when using the global param set)
s = load(fullfile(models.emg.Model.DataDir,'mus.mat'));
fibretypes = [s.mus(1,:) fibretypes];

numsamples = length(fibretypes);
PCPool.open;
for sv = 1:2
    m = models.motorunit.Shorten('SPM',true,'SarcoVersion',sv);
    % Wont need longer as the first peak is usually at about 6ms
    m.T = 60;
    % Fine time resolution to have sufficient detail for fine emg simulations.
    m.dt = .001;
    Shapes = cell(1,numsamples);
    Times = cell(1,numsamples);
    ctimes = zeros(1,numsamples);
    parfor i = 1:numsamples
        [Shapes{i}, Times{i}, ctimes(i)] = m.getActionPotentialShape(fibretypes(:,i));%#ok
    end
    save(sprintf('ShapeData_v%d.mat',sv),'Shapes','Times','fibretypes','ctimes');
end
PCPool.close;