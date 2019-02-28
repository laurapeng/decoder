function [samplesigsAll,genesigsAll,sparsityAll,markerGenesAll,topGenesAll]...
        = Save_cur_factor(samplesigsAll,genesigsAll,sparsityAll,markerGenesAll,topGenesAll,factorInd,samplesigs,genesigs,geneID)

for x = 1:numfactors
    samplesigsAll{factorInd,x} = samplesigs(x,:); 
    genesigsAll{factorInd,x} = genesigs(:,x);
    sparsity = genesigs(:,x)-max(genesigs(:,setdiff(1:size(genesigs,2),x)),[],2);
    sparsityAll{factorInd,x} = sparsity;
   
    % Find marker genes
    [sparsitySrt,order] = sort(sparsity, 'descend');
    geneIDSrt = geneID(order,:);
    sparsitySrt(:,2) = 1:1:length(sparsity);
    sparsitySrtPos = sparsitySrt(sparsitySrt(:,1) > 0, :);
    sparsitySrtPos(:,3) = sparsitySrtPos(:,1)/max(sparsitySrtPos(:,1));
    sparsitySrtPos(:,4) = sparsitySrtPos(:,2)/max(sparsitySrtPos(:,2));
    sparsitySrtPos(:,5) = abs(sparsitySrtPos(:,4) - sparsitySrtPos(:,3));
    [~,minInd] = min(sparsitySrtPos(:,5));
    markerGenesAll{factorInd,x} = geneIDSrt(1:minInd(1),1);

    % Fine top 250 genes
    [sparsitySrt,order] = sort(sparsity, 'descend');
    top = sparsitySrt>quantile(sparsitySrt,1-250/size(sparsitySrt,1));
    geneIDSrt = geneID(order,1);
    topGenesAll{factorInd,x} = geneIDSrt(top,1);
end
