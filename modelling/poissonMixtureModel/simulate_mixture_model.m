clear all
clc
close all

tic

%% SFT Functions
% Set up possible parameters
rate = .02; % Low  Drift in channel B
parallelRate = rate - .012;

pAB = 0;
pBA = pAB;

cA = 10; % Criterion A
cB = 10; % Criterion B

pcA = 10; % Criterion A
pcB = 10; % Criterion B

t = 0:3000; % Time vector to evaluate

% pset = 0:.25:1;
pset = [.1 .5 .9];

drift = [rate rate];

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get Parallel cdfs
parallel.drift.AND = [parallelRate, parallelRate];

parallel.cdf.AND = fac_and([parallel.drift.AND, 0, 0, pcA, pcB]', t)';
parallel.pdf.AND = diff([0; parallel.cdf.AND]);


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get serial pdfs
pX = .5;
    A = singChan([drift(1), cA]', t)'; % PDF channel A
    B = singChan([drift(2), cB]', t)'; % PDF channel B
    
    pdfA = diff([0; A]);
    pdfB = diff([0; B]);
    
    serial.pdf.AND = fftConv(pdfA, pdfB);
    serial.cdf.AND = cumsum(serial.pdf.AND);         % Exhaustive - convolve


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Mix Models
cnt = 1;
for p = pset
    %%
    mixedSerPar.pdf.AND(:,cnt) = p * serial.pdf.AND + (1 - p) * parallel.pdf.AND;
    mixedSerPar.cdf.AND(:,cnt) = cumsum(mixedSerPar.pdf.AND(:,cnt));
    
    cnt = cnt + 1;
end

toc;  

%% Sample RTs from model at three levels of mixture
nsamples = 400;

mixedSerPar.pdf.AND(mixedSerPar.pdf.AND < 0) = 0;

Separation{1} = randsample(t, nsamples, true, mixedSerPar.pdf.AND(:,1));
Separation{1} = Separation{1}'/1000;

Separation{2} = randsample(t, nsamples, true, mixedSerPar.pdf.AND(:,2));
Separation{2} = Separation{2}'/1000;

Separation{3} = randsample(t, nsamples, true, mixedSerPar.pdf.AND(:,3));
Separation{3} = Separation{3}'/1000;

save('mixedSerialParallel.mat', 'Separation')