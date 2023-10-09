%% Simulate Interactive Parallel Model and plot results
% Top row are the results for the facilitatory parallel model
% Bottom row are the results for the inhibitory parallel model
%
% Columns should be: SIC_OR, SIC_AND, C_OR, C_AND, CCF

clear all
clc
close all
tic

colors = linspace(0, .6, 5)' * [1 1 1];

%%
lA = .02; % Low  Drift in channel A

pAB = [0 .25 .5 .75 1];
pBA = pAB;

cA = 10; % Criterion A
cB = 10; % Criterion B

t = 0:5000; % Time vector to evaluate

drift = [lA lA];

%% Get cdfs 
for i = 1:numel(pAB)
        inhibitory.cdf.and(:,i) = inh_and([drift, pAB(i), pBA(i), cA, cB]', t)';
        parallel.pdf.and(:,i) = diff([0; inhibitory.cdf.and(:,i)]);
end
    
%% Sample RTs from model at three levels of mixture
nsamples = 400;

parallel.pdf.and(parallel.pdf.and < 0) = 0;

Separation{1} = randsample(t, nsamples, true, parallel.pdf.and(:,2));
Separation{1} = Separation{1}'/1000;

Separation{2} = randsample(t, nsamples, true, parallel.pdf.and(:,3));
Separation{2} = Separation{2}'/1000;

Separation{3} = randsample(t, nsamples, true, parallel.pdf.and(:,4));
Separation{3} = Separation{3}'/1000;

save('inhibitoryParallel.mat', 'Separation')