clear all
clc
% close all

%% Load data file
subjectNumber = 6;
data = dlmread (fullfile('Data', sprintf('2018_MMT_FaceRules_subject%d.dat', subjectNumber))); %change subject number in data file

ploton = true; % Set to true to display plots
minrt = 200; % Minimum RT cutoff 200 ms 

%% Data file column names (some columns are necessary: 'sub', 'itm', 'rt'
cols = {'sub', 'con', 'sep', 'tri', 'itm', 'rsp', 'acc', 'rt'};

%% Remove overall long RTs
% cutoffPercentiles = repmat([99.9], 1, numel(subjectNumber)); % Throw out anything greater than this prctile (if 100 then don't throw out anything)
% This is useful for getting rid of extremely long RTs if the P, for instance, took a phonecall or fell asleep
% idxTooLong = find(data(:,end) > prctile(data(:,strcmp(cols, 'rt')), cutoffPercentiles(si)));
%
% 180404 DL: Changed to just throw out anything > 3 secs
idxTooLong = find(data(:,strcmp(cols, 'rt')) > 4);
data(idxTooLong, :) = []; % Delete long RTs
fprintf('Number of long RT trials removed = %d\n', numel(idxTooLong))

%% Remove timeouts
idx9 = find(data(:,strcmp('acc', cols)) == 9); % Remove timeouts
data(idx9,:) = []; % Delete timeouts
fprintf('Number of timeouts removed = %d\n', numel(idx9))
data(isnan(data(:,strcmp('rt', cols))), :) = []; % Delete nans

%% Remove error RTs
idx0 = find(data(:,strcmp('acc', cols)) == 0); % Remove errors
data(idx0,:) = []; % Delete errors
fprintf('Number of errors removed = %d\n', numel(idx0))
data(isnan(data(:,strcmp('rt', cols))), :) = []; % Delete error RTs

%% Convert RTs to msecs if not already in secs
% if max(data(:,end)) < 1000 
%     data(:,end) = data(:,end) * 1000;
% end

%% Extract RT data for all items in each separation condition
condition_number = data(:,2);
Separation{1} = data(condition_number == 2, strcmp(cols, 'rt'));
Separation{2} = data(condition_number == 3, strcmp(cols, 'rt'));
Separation{3} = data(condition_number == 4, strcmp(cols, 'rt'));

%% MCMC to find posterior gamma parameters
n.burnin = 1000; 
n.mcmc = 100000; 
n.chains = 4;
n.thin = 100; 

% Load the sample if they exist, if not, run the sampler then save the
% samples
if exist(sprintf('subject%d_samples_eg.mat', subjectNumber), 'file') == 2
    load(sprintf('subject%d_samples_eg.mat', subjectNumber))
else
    tic
    parpool('local')
    parfor i = 1:3
        samples{i} = sampleExGauss_fast(Separation{i}, n);
    end
    samples = cellfun(@(x)(exp(x)), samples, 'uni', 0);
    
    delete(gcp('nocreate'));
    toc
    save(sprintf('subject%d_samples_eg.mat', subjectNumber), 'samples')
end


%% Check convergence
figure('WindowStyle', 'docked');
titles = {'0', '100', '350'};
for i = 1:3
subplot(6,3,i); 
plot(samples{i}(n.burnin+1:n.thin:n.mcmc,:,1)); title(sprintf('Separation %s: %s parameter', titles{i}, '\mu')); xlabel('Iteration'); ylabel('Value');
subplot(6,3,i+3); 
plot(samples{i}(n.burnin+1:n.thin:n.mcmc,:,2)); title(sprintf('Separation %s: %s parameter', titles{i}, '\sigma')); xlabel('Iteration'); ylabel('Value');
subplot(6,3,i+6); 
plot(samples{i}(n.burnin+1:n.thin:n.mcmc,:,2)); title(sprintf('Separation %s: %s parameter', titles{i}, '\tau')); xlabel('Iteration'); ylabel('Value');
end

post{1} = samples{1}(n.burnin+1:n.thin:n.mcmc,:,:);
post{2} = samples{2}(n.burnin+1:n.thin:n.mcmc,:,:);
post{3} = samples{3}(n.burnin+1:n.thin:n.mcmc,:,:);

%% Sample from the posterior
nPostSamples = 1000;
rt = 0:.0005:4; 

f = @(x)(x(:));

for i = 1:3 % ndataset
for j = 1:3 % nparms
    parms{i}(:,j) = randsample(f(post{i}(:,:,j)), nPostSamples, false);
end
end

%%
dens1 = nan(nPostSamples, numel(rt)); dens2 = dens1; dens3 = dens1;
d12 = nan(nPostSamples, numel(rt)); d13 = d12; d23 = d12; 

tic
for i = 1:nPostSamples
   
    %% Now estimate the pdf based on the posterior
    for j = 1:3 % nData
        dens{j}(i,:) = exgausspdf(rt, parms{j}(i,1), parms{j}(i,2), parms{j}(i,3));
    end
    
    % Subtract gamma densities
    d12(i,:) = minus(dens{1}(i,:), dens{2}(i,:));
    d13(i,:) = minus(dens{1}(i,:), dens{3}(i,:));
    d23(i,:) = minus(dens{2}(i,:), dens{3}(i,:));
    
    %% Find the fixed points
    zc(i,:) = findZeroCrossing(rt, [d12(i,:); d13(i,:); d23(i,:)]);
end
toc 

% Find the mean and standard deviation of the bootstrapped estimates
for i = 1:3
    m{i} = mean(dens{i});
    s{i} = std(dens{i});
end

m12 = mean(d12);
s12 = std(d12);

m13 = mean(d13);
s13 = std(d13);

m23 = mean(d23);
s23 = std(d23);

%% Plot probability densities by condition number
colors = {[1 .5 .5], [.5 1 .5], [.5 .5 1]};
% fig = figure('WindowStyle', 'docked');
subplot(2,3,4)
for i = 1:3 % ndata
figh(i)=fill([rt, fliplr(rt)], [m{i} + 2 * s{i}, fliplr(m{i} - 2 * s{i})], colors{i}, 'FaceAlpha', .5, 'EdgeColor', colors{i});
hold on
plot(rt, m{i}, 'Color',  colors{i}, 'LineWidth', 2)
end
legend(figh, '0 pixel separation', '100 pixel separation', '350 pixel separation', 'Location', 'NorthEast')
xlabel('Time (sec)', 'FontSize', 14)
ylabel('Density', 'FontSize', 14)

%% Plot density difference functions
subplot(2,3,5)
f1=fill([rt, fliplr(rt)], [m12 + 2 * s12, fliplr(m12 - 2 * s12)], [1 .75 .5], 'FaceAlpha', .5, 'EdgeColor', [1 .75 .5]);
hold on
plot(rt, m12, 'Color',  [1 .75 .5], 'LineWidth', 2)
f2=fill([rt, fliplr(rt)], [m13 + 2 * s13, fliplr(m13 - 2 * s13)], [.75 1 .5], 'FaceAlpha', .5, 'EdgeColor', [.75 1 .5]);
plot(rt, m13, 'Color',  [.75 1 .5], 'LineWidth', 2)
f3=fill([rt, fliplr(rt)], [m23 + 2 * s23, fliplr(m23 - 2 * s23)], [.5 .75 1], 'FaceAlpha', .5, 'EdgeColor', [.5 .75 1]);
plot(rt, m23, 'Color',  [.5 .75 1], 'LineWidth', 2)
legend([f1, f2, f3], '0 pixel - 100 pixel', '0 pixel - 350 pixel', '100 pixel - 350 pixel', 'Location', 'NorthEast')
xlabel('Time (sec)', 'FontSize', 14)
ylabel('Density Difference', 'FontSize', 14)
line([0, 4], [0 0], 'Color', 'k', 'LineWidth', 2, 'LineStyle', '--')

%% Plot fixed point histogram
subplot(2,3,6)
bw = getbandwidth(zc(:));
for i = 1:3
    [density(:,i), bins(:,i)] = ksdensity(zc(:,i), 'kernel', 'epanechnikov', 'bandwidth', bw); % This uses KDE
%     h(i) = plot(bins(:,i), dens(:,i), 'Color', colors{i}, 'LineWidth', 2);
    h(i) = fill([bins(:,i); flipud(bins(:,i))], [density(:,i); zeros(numel(density(:,i)), 1)], colors{i}, 'FaceAlpha', .5, 'EdgeColor', colors{i});
    hold on
end
legend(h, '0 pixel - 100 pixel', '0 pixel - 350 pixel', '100 pixel - 350 pixel', 'Location', 'NorthEast')
xlabel('Fixed Point Estimate', 'FontSize', 14)
ylabel('Density', 'FontSize', 14)
set(gca,'XLim',[0.5,2.5]);

supertitle(sprintf('Subject %d', subjectNumber));

%% Compute the amount of overlap between 95% HDI's
cdfs = cumsum(density)./(ones(size(density,1), 1) * sum(density));
hdi = nan(2,3);
for i = 1:3
    hdi(1,i) = bins(find(cdfs(:,i)<=.025, 1, 'last'));
    hdi(2,i) = bins(find(cdfs(:,i)>=.985, 1, 'first'));
end

% XY coordinates of line 1: (hdi(1,1),0)-(hdi(2,1),0)
% XY coordinates of line 2: (hdi(1,2),0)-(hdi(2,2),0)
% XY coordinates of line 3: (hdi(1,3),0)-(hdi(2,3),0)
pairs = allcomb(1:3, 1:3); pairs(pairs(:,1) == pairs(:,2), :) = [];
pairs = unique(sort(pairs, 2), 'rows');
overlap = nan(1,3);
for i = 1:size(pairs, 1)
    % Get pairs
    x = hdi(:,pairs(i,:));
    
    % Sort so that the lines are ordered appropriately
    if x(1,1) ~= x(1,2)
        [a, b] = sort(x(1,:), 2);
        x = [a; x(2,b)];
    else
        [a, b] = sort(x(2,:), 2);
        x = [x(1,b); a];
    end
    
    % Check overlap
    if x(1,2) < x(2,1) && x(2,1) <= x(2,2) % Then the lines overlap
        overlap(i) = (x(2,1) - x(1,2))./(x(2,2) - x(1,1)); % Proportion overlap
    elseif x(1,2) < x(2,1) && x(2,1) > x(2,2) % Then they completely overlap
        overlap(i) = 1;
    else
        overlap(i) = 0; % If there's no overlap, then the proportion overlap is 0
    end
end
shared = [max(hdi(1,:)); min(hdi(2,:))];
sharedoverlap = zeros(1,3);
sharedoverlap(overlap > 0) = shared(2)-shared(1)./(hdi(2,overlap>0) - hdi(1,overlap>0)); 


%% 
conditions = [0 100 350];
figure
for i = 1:3
subplot(1,3,i)
t = linspace(0, 4, 50);
[c, e] = hist(Separation{i}, t);
b = bar(e, c./trapz(e,c), 'hist'); % * mean(tmp.dataacc == 1), 'hist');
set(b, 'FaceColor', [0 .95 .95])
hold on

for j = 1:50
    ppdens = exgausspdf(t, parms{i}(j,1), parms{i}(j,2), parms{i}(j,3));
    plot(t, ppdens, '-r'); 
end
set(gca,'XLim', [0, 4])
xlabel('time')
ylabel('density')
title(sprintf('Separation %d Pixels', conditions(i)))
end
% supertitle(sprintf('Posterior Predictive Fits of exGaussian to Subject %d', subjectNumber))

%% Summarize data
fprintf('Proportion Correct')
pcorrect = [cellfun(@(x)(numel(x)/400), Separation)]
% fprintf('Mean RT')
% meanRT = cellfun(@(x)(mean(x*1000)), Separation)
% fprintf('standard error')
% sterr = cellfun(@(x)(std(x*1000)./sqrt(numel(x))), Separation)

fprintf('Median RT')
meanRT = cellfun(@(x)(median(x*1000)), Separation)
fprintf('IQ Range')
sterr = cellfun(@(x)(prctile(x*1000, [25 75])), Separation, 'UniformOutput',false)