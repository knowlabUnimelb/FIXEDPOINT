clear all
clc
close all

tic

%% SFT Functions
% Set up possible parameters
rate = .02; % Low  Drift in channel B
parallelRate = rate - .01;

pAB = 0;
pBA = pAB;

cA = 10; % Criterion A
cB = 10; % Criterion B

serial_cA = [9 10 11]; 
serial_cB = [10 11 12]; 

t = 0:3000; % Time vector to evaluate

% pset = [0 .25 .5 .75 1];
pset = [.1 .5 .9];

drift = [rate rate; rate+.005 rate+.005; rate+.01 rate+.01];

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get Parallel cdfs
parallel.drift.AND = [parallelRate, parallelRate];

parallel.cdf.AND = fac_and([parallel.drift.AND, 0, 0, cA, cB]', t)';
parallel.pdf.AND = diff([0; parallel.cdf.AND]);


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get serial pdfs
% serial 1
pX = .5;
A = singChan([drift(1,1), serial_cA(1)]', t)'; % PDF channel A
B = singChan([drift(1,2), serial_cB(1)]', t)'; % PDF channel B

pdfA = diff([0; A]);
pdfB = diff([0; B]);

serial.pdf.AND(:,1) = fftConv(pdfA, pdfB);
serial.cdf.AND(:,1) = cumsum(serial.pdf.AND(:,1));         % Exhaustive - convolve

% serial 2
pX = .5;
A = singChan([drift(2,1), serial_cA(2)]', t)'; % PDF channel A
B = singChan([drift(2,2), serial_cB(2)]', t)'; % PDF channel B

pdfA = diff([0; A]);
pdfB = diff([0; B]);

serial.pdf.AND(:,2) = fftConv(pdfA, pdfB);
serial.cdf.AND(:,2) = cumsum(serial.pdf.AND(:,2));         % Exhaustive - convolve

% serial 2
pX = .5;
A = singChan([drift(3,1), serial_cA(3)]', t)'; % PDF channel A
B = singChan([drift(3,2), serial_cB(3)]', t)'; % PDF channel B

pdfA = diff([0; A]);
pdfB = diff([0; B]);

serial.pdf.AND(:,3) = fftConv(pdfA, pdfB);
serial.cdf.AND(:,3) = cumsum(serial.pdf.AND(:,3));         % Exhaustive - convolve

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Mix Models
cnt = 1;
for p = pset
    %%
    mixedSerPar.pdf.AND(:,cnt) = p * serial.pdf.AND(:,cnt) + (1 - p) * parallel.pdf.AND;
    mixedSerPar.cdf.AND(:,cnt) = cumsum(mixedSerPar.pdf.AND(:,cnt));
    
    cnt = cnt + 1;
end
toc

%% Sample RTs from model at three levels of mixture
nsamples = 400;

mixedSerPar.pdf.AND(mixedSerPar.pdf.AND < 0) = 0;

Separation{1} = randsample(t, nsamples, true, mixedSerPar.pdf.AND(:,1));
Separation{1} = Separation{1}'/1000;

Separation{2} = randsample(t, nsamples, true, mixedSerPar.pdf.AND(:,2));
Separation{2} = Separation{2}'/1000;

Separation{3} = randsample(t, nsamples, true, mixedSerPar.pdf.AND(:,3));
Separation{3} = Separation{3}'/1000;

save('mixedSerialParallel_noSelectiveInfluence.mat', 'Separation')