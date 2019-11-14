function [ccmat] = Train_NMF_seed(data,master_subset,numfactors)

ccmat = zeros(size(data,1));

% Pick out a subset of 80% samples
subset = master_subset;
subset(subset)=rand(sum(subset),1)>0.2;

% Run quick NMF
[genesigs,samplesigs] = nnmf_unmix_quick(data(:,subset),numfactors);
if size(genesigs,2)<numfactors
    disp('converged on fewer factors...')
    fprintf('found %d / %d \n',size(genesigs,2),numfactors);
end

% Normalize the gene factor loadings to a median of 1 
for i = 1:size(samplesigs,1)
    f = mean(genesigs(:,i));
    samplesigs(i,:)=samplesigs(i,:)*f;
    genesigs(:,i)=genesigs(:,i)/f;
end

% Pick out the exemplar genes from each factor 
sparsity = zeros(size(genesigs,1),1);
for j = 1:size(genesigs,2)
    sparsity(:,j)=genesigs(:,j)-max(genesigs(:,setdiff(1:size(genesigs,2),j)),[],2);
end
for j = 1:size(genesigs,2)
    sfilter=find(sparsity(:,j)>quantile(sparsity(:,j),1-50/size(genesigs,1)));
    ccmat(sfilter,sfilter) = ccmat(sfilter,sfilter)+1;
end
