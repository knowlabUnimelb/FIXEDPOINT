function samples = sampleExGauss_fast(data, n)
% n.burnin = number of burnin samples
% n.mcmc = number of total mcmc

samples = nan(n.mcmc, n.chains, 3);
propvar = [.05, .05, .01]; % Variance for each of the parms

samples(1,:,:) = repmat(log([.8, .15, .18]), n.chains, 1);
for i = 2:n.mcmc
    samples(i,:,:) = samples(i-1,:,:);
    
    for j = 1:3
        % Propose a new value
        proposal = squeeze(samples(i,:,:));
        proposal(:,j) = proposal(:,j) + randn(n.chains,1) * propvar(j);
        
        % Get new likelihood
        newloglik = sum(log(exgausspdf(...
            data * ones(1,n.chains),...
            ones(numel(data),1) * exp(proposal(:,1))',...
            ones(numel(data),1) * exp(proposal(:,2))',...
            ones(numel(data),1) * exp(proposal(:,3))' )));
        
        % Get old likelihood
        oldloglik = sum(log(exgausspdf(...
            data * ones(1, n.chains),...
            ones(numel(data),1) * squeeze(exp(samples(i,:,1))),...
            ones(numel(data),1) * squeeze(exp(samples(i,:,2))),...
            ones(numel(data),1) * squeeze(exp(samples(i,:,3))) )));
        
        % Likelihood ratio
        ratio = exp(newloglik-oldloglik);
        
        % If ratio is high enough accept proposal
        p = rand(1,n.chains);
        samples(i, p < ratio, :) = proposal(p < ratio,:);
    end
end
