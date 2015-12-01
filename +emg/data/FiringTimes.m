% Script to pre-compute several firing times for faster use
% within the emg model. This is optional, as the correct firing times for the
% selected motor-units can be computed directly for the fibre types used
% in each emg model.
%
% Note: For the experiments, this has not been used as the stand-alone
% firing times sometimes (due to discrete ODE solvers) do not match the
% firing times obtained via the direct computation of the full musclefibre
% model
clear classes;%#ok

m = models.motoneuron.Model;
% 10000 ms should be enough in terms of precomputed times
m.T = 10000;
% Fine time resolution to have sufficient detail for fine emg simulations.
m.dt = .1;
% Use samples with always full activation to have immediate peaks.
s = sampling.RandomSampler;
s.Domain = models.motoneuron.ParamDomain;
s.Samples = 200;
m.Sampler = s;

m.off1_createParamSamples;
mus = m.Data.ParamSamples;

% Also add the fibre types from "the" global parameter set (faster
% tests/experiments as precomp can be used when using the global param set)
s = load(fullfile(models.emg.Model.DataDir,'mus.mat'));
mus = [s.mus mus];
numsamples = size(mus,2);

Times = cell(1,numsamples);
ctimes = zeros(1,numsamples);
PCPool.open;
parfor i = 1:numsamples
    [t,y,ctimes(i)] = m.simulate(mus(:,i), 1);%#ok
    peaks = m.analyze(t,y);
    Times{i} = t(peaks);
end
PCPool.close;
save FiringTimes Times mus ctimes;