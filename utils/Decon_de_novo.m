function Decon_de_novo(configFileName) 

%% ========== Preprocess ==================================================
tic
% Add path
binDECODER = fileparts(mfilename('fullpath'));
binDECODER = fileparts(binDECODER);
addpath(binDECODER)
addpath(fullfile(binDECODER,'data'))
addpath(fullfile(binDECODER,'utils'))

% Read and parse configure file 
configFile = fopen(configFileName);
configInfo = textscan(configFile,'%s\t%s');
fclose(configFile);

% Set parameters
tmpInd = find(strcmp(configInfo{1},'geneIDType'));
geneIDType = configInfo{2}{tmpInd};
tmpInd = find(strcmp(configInfo{1},'repTimes'));
repTimes = str2double(configInfo{2}{tmpInd});

tmpInd = find(strcmp(configInfo{1},'rangeK'));
rangeK = configInfo{2}{tmpInd};
switch rangeK
    case 'auto'
	minFactorNum = 2;
	maxFactorNum = 25;
    otherwise
	tmp = split(rangeK,':')
	minFactorNum = str2num(tmp{1});
	maxFactorNum = str2num(tmp{2});
end

tmpInd = find(strcmp(configInfo{1},'dataType'));
dataType = configInfo{2}{tmpInd};
switch dataType
    case 'RNAseq'
        numTrainGene=5000;
    case 'Microarray'
        numTrainGene=5000;
    case 'ATACseq'
        numTrainGene=8000;
    otherwise
        disp('Error: Data type not supported yet...');
        return
end

% Load data matrix
tmpInd = find(strcmp(configInfo{1},'dataMatrix'));
dataMatrix = configInfo{2}{tmpInd};
outDir = fileparts(dataMatrix);

tmpInd = find(strcmp(configInfo{1},'dataFormat'));
dataFormat = configInfo{2}{tmpInd};
[data,geneID,sampleID,master_subset] = Load_data(dataMatrix,dataFormat)

tmpInd = find(strcmp(configInfo{1},'logTransformed'));
logTransformed = configInfo{2}{tmpInd};
if strcmp(logTransformed,'no')
    data = log2(1+data);
end

% sample subset
original.data = data;
original.geneID = geneID;

% highly expressed genes
x = mean(data,2);
x = x<quantile(x,0.25); % upper 75%
data(x,:)=[];
geneID(x)=[];

% variable genes
x = std(2.^data(:,:),[],2);
x = x<quantile(x,1-numTrainGene/numel(x)); % upper numTrainGene genes
data(x,:)=[];
geneID(x)=[];

selected.data = data;
selected.geneID = geneID;

timePreproc = toc;
fprintf('Preprocessing completed, %3.0f minutes elapsed...\n',...
         timePreproc/60)
clear configFile configFileName 
clear dataFormat dataType dataMatrix logTransformed numTrainGene
clear rawTable x tmpInd tmp ans

%% ========== Run NMF through increasing number of K ======================
timeKLoop = timePreproc;
factorInfo = zeros((maxFactorNum-minFactorNum+1),3);
samplesigsAll = cell((maxFactorNum-minFactorNum+1),maxFactorNum);
genesigsAll = cell((maxFactorNum-minFactorNum+1),maxFactorNum);
sparsityAll = cell((maxFactorNum-minFactorNum+1),maxFactorNum);
topGenesAll = cell((maxFactorNum-minFactorNum+1),maxFactorNum);
markerGenesAll = cell((maxFactorNum-minFactorNum+1),maxFactorNum);
factorScore = zeros((maxFactorNum-minFactorNum+1),maxFactorNum);
factorMatch = zeros((maxFactorNum-minFactorNum+1),maxFactorNum);
factorInd = 0;
factorStartInd = 1;

for numfactors = minFactorNum:maxFactorNum
    factorInd = factorInd + 1;
    factorEndInd = factorInd;
    factorInfo(factorInd,1) = numfactors; 
    %%% ===== Train gene weight seed on 5K genes through repetitions ======
    timeIteration = timeKLoop;
    data = selected.data;
    geneID = selected.geneID ;
    rng(999);
    ccmat = zeros(size(data,1));
    for repetition = 1:repTimes
        tic;
        [ccmat] = Train_NMF_seed(data,master_subset,numfactors);
        timeIteration = timeIteration+toc;
        fprintf('K=%d: seed training - %3.0f iterations completed, %3.0f remaining, %3.0f minutes elapsed...\n',numfactors, repetition,repTimes-repetition,timeIteration/60)
    end
    save(fullfile(outDir, sprintf('K%d_ccmat_rep%d.mat',numfactors,repTimes)),'ccmat')
    
    %%% ===== Refine results and compute official 5k NMF ==================
    tic;
    % Pick genes which were factor-exemplars together and often
    genefilter = find(max(ccmat)>=max(1,quantile(max(ccmat),1-(50*numfactors)/size(data,1))));
    
    % Cluster these genes based on factor co-appearnace to get groups
    Consensus = ccmat(genefilter,genefilter)/repTimes;
    Consensus = (corr(Consensus)*-1)+1;
    z = linkage(squareform(Consensus),'average');
    % membership_y = cluster(z,'maxclust',numfactors);
    [~,membership_y,~]= dendrogram(z,numfactors,'orientation','left');
    [~,~,membership_y] = unique(membership_y);
    
    % Calculate cophenet correlation
    Y = pdist(z);
    Z = linkage(Y,'average');
    [coph,~] = cophenet(Z,Y);
    factorInfo(factorInd,2) = coph;
    
    % Rebalance genefilter list to have even coverage across groups
    old_genefilter = genefilter;
    for i = 1:numfactors
        j=old_genefilter(membership_y==i);
        k=find(mean(ccmat(j,j))<max(1,quantile(mean(ccmat(j,j)),1-(min(numel(j),25))/numel(j))));
        genefilter = setdiff(genefilter,j(k));
    end
    
    % Cluster these genes to find group thresholds
    Consensus = ccmat(genefilter,genefilter)/repTimes;
    Consensus = (corr(Consensus)*-1)+1;
    z = linkage(squareform(Consensus),'average');
    [tree,membership_y,~]= dendrogram(z,numfactors,'orientation','left');
    thresh = inf;
    for i = 1:length(tree)
        thresh = min(thresh,max(get(tree(i),'Xdata')));
    end
    
    % Cluster exemplar genes into their balanced categories
    [~,~,perm_y]= dendrogram(z,0,'orientation','right','colorthreshold',thresh);
    
    % Run last NMF on all samples using exemplar genes
    subset = master_subset;
    genesigs = ones(length(genefilter),numfactors);
    for i = 1:numfactors
        genesigs(membership_y(perm_y)==i,i) = 100;
    end
    
    % first on just the exemplar genes 100:1
    [genesigs,samplesigs] = nnmf(data(genefilter(perm_y),subset),numfactors,'w0',genesigs,'algorith','als');
    if ~isempty(find(all(genesigs==0)))
        tmpK = length(find(~all(genesigs==0)));
        data = original.data;
        geneID = original.geneID;
        if strcmp(autoStop,'yes') % stop if converged; cut tails (1|2) if < 0.5
            factorInd = factorInd - 1;
                factorEndInd = factorInd;
            if factorInfo(factorEndInd,3) < 0.5 && factorEndInd>=2
                    if factorInfo(factorEndInd-1,3) < 0.5
                factorEndInd = factorEndInd-2;
                    else
                        factorEndInd = factorEndInd-1;
                    end
            end
            fprintf('K=%d: Converged on %d factors, thus stopped...\n',numfactors, tmpK)
            timeKLoop = timeIteration + toc;
            break
        else
            fprintf('K=%d: Converged on %d factors, thus skipped...\n',numfactors, tmpK)
            timeKLoop = timeIteration + toc;
            continue
        end
    end
    for i = 1:size(samplesigs,1)
        f = mean(genesigs(:,i));
        samplesigs(i,:)=samplesigs(i,:)*f;
        genesigs(:,i)=genesigs(:,i)/f;
    end
    
    % then on 5k genes with better seeds
    temp = 0.01*ones(size(data,1),numfactors);
    temp(genefilter(perm_y),:) = genesigs;
    genesigs = temp;
    [genesigs,samplesigs] = nnmf(data(:,subset),numfactors,'w0',genesigs,'algorith','als');
    if ~isempty(find(all(genesigs==0)))
	tmpK = length(find(~all(genesigs==0)));
        data = original.data;
        geneID = original.geneID;
        if strcmp(autoStop,'yes')
            factorInd = factorInd - 1;
                factorEndInd = factorInd;
            if factorInfo(factorEndInd,3) < 0.5 && factorEndInd>=2
                    if factorInfo(factorEndInd-1,3) < 0.5
                        factorEndInd = factorEndInd-2;
                    else
                        factorEndInd = factorEndInd-1;
                    end
            end
            fprintf('K=%d: Converged on %d factors, thus stopped...\n',numfactors, tmpK)
            timeKLoop = timeIteration + toc;
            break
        else
            fprintf('K=%d: Converged on %d factors, thus skipped...\n',numfactors, tmpK)
            timeKLoop = timeIteration + toc;
            continue
        end
    end
    for i = 1:size(samplesigs,1)
        f = mean(genesigs(:,i));
        samplesigs(i,:)=samplesigs(i,:)*f;
        genesigs(:,i)=genesigs(:,i)/f;
    end
    clear ans Consensus repetition
    clear D Y Z f i j k
    clear ccmat genefilter old_genefilter perm_y membership_y
    clear temp thresh tree z 
    
    %%% ===== Project all genes and samples using NNLS ===================
    fprintf('K=%d: projecting using NNLS for all samples...\n',numfactors)
    sampleProjection = zeros(size(genesigs,2),size(data,2));
    for i = 1:size(sampleProjection,2) % for each sample
        [sampleProjection(1:size(genesigs,2),i),~] = lsqnonneg(genesigs , data(:,i) );
    end
    
    fprintf('K=%d: projecting using NNLS for all genes...\n',numfactors)
    data = original.data;
    geneID = original.geneID;
    geneProjection = zeros(numel(geneID),size(genesigs,2));
    for i = 1:size(geneProjection,1) % for each gene
        [geneProjection(i,:),~] = lsqnonneg(samplesigs',data(i,master_subset)');
    end
    
    genesigs = geneProjection; clear geneProjection
    samplesigs = sampleProjection; clear sampleProjection
    
    clear i
    save(fullfile(outDir,sprintf('K%d_res.mat',numfactors)),'genesigs','samplesigs','coph')
    fprintf('K=%d: factors identified...\n',numfactors)
    
    %%% ===== Process current factors =====================================
    fprintf('K=%d: saving current results...\n',numfactors)
    [samplesigsAll,genesigsAll,sparsityAll,markerGenesAll,topGenesAll]...
        = Id_cur_genes(samplesigsAll,genesigsAll,sparsityAll,markerGenesAll,topGenesAll,factorInd,samplesigs,genesigs,geneID);

    % calculate factor score and link factors
    fprintf('K=%d: factor linkages establishing with the previous run...\n',numfactors)
    [factorScore,factorMatch] = Link_previous_factor(factorInd,numfactors,factorInfo,topGenesAll);
    
    %%% ====== Evaluate current K =========================================
    fprintf('K=%d: current K evaluating...\n',numfactors)
    % calculate median of factor scores
    factorInfo(factorInd,3) = median(factorScore(factorInd,factorScore(factorInd,:)~=0));
    
    % Determine K to start and K to stop 
    if strcmp(autoStop,'yes')
        % find the first score that > 0.5
        if factorInd ~= 1
            if (factorStartInd==1) && (factorInfo(factorInd,3) >= 0.5)
                factorStartInd = factorInd;
            end
        end
        % stop if 4 continous < 0.5
        if factorStartInd~=1
            if factorInd-factorStartInd+1>=9
                if ((factorInfo(factorInd,3) < 0.5)...
                        && (factorInfo((factorInd-1),3) < 0.5)...
                        &&  (factorInfo((factorInd-2),3) < 0.5)...
                        &&  (factorInfo((factorInd-3),3) < 0.5))
                    factorEndInd = factorInd-4;
                    break
                    timeKLoop = timeIteration + toc;
                end
            end	    
        end
    end    
    timeKLoop = timeIteration + toc;
    clear coph genesigs samplesigs
end

tic
factorInfo = factorInfo(factorStartInd:factorEndInd,:);
samplesigsAll = samplesigsAll(factorStartInd:factorEndInd,:);
genesigsAll = genesigsAll(factorStartInd:factorEndInd,:);
sparsityAll = sparsityAll(factorStartInd:factorEndInd,:);
topGenesAll = topGenesAll(factorStartInd:factorEndInd,:);
markerGenesAll = markerGenesAll(factorStartInd:factorEndInd,:);
factorScore = factorScore(factorStartInd:factorEndInd,:);
factorMatch = factorMatch(factorStartInd:factorEndInd,:);

timeKLoop = timeKLoop + toc;
if strcmp(autoStop,'yes')
    if factorStartInd~=1
	fprintf('K=%d:%d(auto): NMF completed, %3.0f minutes elapsed...\n',factorInfo(1,1),factorInfo(end,1),timeKLoop/60)
    else
	fprintf('K=%d:%d(auto): warning: no robust K found...\n',factorInfo(1,1),factorInfo(end,1))
	fprintf('K=%d:%d(auto): NMF completed, %3.0f minutes elapsed...\n',factorInfo(1,1),factorInfo(end,1),timeKLoop/60)
    end
else
    fprintf('K=%d:%d(user defined): NMF completed, %3.0f minutes elapsed...\n',factorInfo(1,1),factorInfo(end,1),timeKLoop/60)
end
clear factorInd factorStartInd factorEndInd 
clear master_subset numfactors original repTimes selected subset autoStop

%% ========== Determine compartments from factors in multiple runs ========
tic
disp('Compartment selecting...')
% build linkages 
[factorLink,factorLinkScore] = Build_factor_link(factorInfo,factorMatch,factorScore);

% estimate threshold
[lowThre,highThre] = Est_score_threshold(factorLinkScore);

% pick compartment candidates
sltComp = {};
sltFactorLink = {};
sltFactorLinkScore = {};
sltFactorLinkFull= {};
sltFactorLinkFullScore = {};
n = 0;
for i = 2:size(factorLinkScore,1)
    [resIdx,factorType,nextFlag] = Pick_comp_candidates(factorLinkScore,i,lowThre,highThre);
    
    if nextFlag == 1
        continue
    else
        numfactors = factorInfoFlt(resIdx,1);
        factorInd = factorLink{i,resIdx};
        res = sprintf('%d.%d',numfactors,factorInd);
        if ~any(strcmp(sltComp,res))
            n = n+1;
            sltComp{n,1} = res;
            sltComp{n,2} = factorLinkScore{i,resIdx};
            sltComp{n,3} = factorType;
            sltFactorLink(n,1:resIdx) = factorLink(i,1:resIdx);
            sltFactorLinkScore(n,1:resIdx) = factorLinkScore(i,1:resIdx);
            sltFactorLinkFull(n,:) = factorLink(i,:);
            sltFactorLinkFullScore(n,:) = factorLinkScore(i,:);
        end   
    end
end
[sltComp] = Proc_comp_candidates(sltFactorLink,topGenesAll,sltComp);

% save 
for i = 1:size(sltComp,1)
    cur = sltComp{i,1};
    tmp = strsplit(cur,'.');
    samplesigs = samplesigsAll{factorInfo(:,1) == str2num(tmp{1,1}),str2num(tmp{1,2})};    
    sltComp{i,4} = samplesigs;
    genesigs = genesigsAll{factorInfo(:,1) == str2num(tmp{1,1}),str2num(tmp{1,2})};
    sltComp{i,5} = genesigs;
    sparsity = sparsityAll{factorInfo(:,1) == str2num(tmp{1,1}),str2num(tmp{1,2})};
    sltComp{i,6} = sparsity; 
    genes250 = topGenesAll{factorInfo(:,1) == str2num(tmp{1,1}),str2num(tmp{1,2})};
    sltComp{i,7} = genes250;
    genes = markerGenesAll{factorInfo(:,1) == str2num(tmp{1,1}),str2num(tmp{1,2})};
    sltComp{i,8} = genes;
    clear genes genes250 cur tmp spartsity genesigs samplesigs 
end

% unify genes
sltComp(:,end+1) = sltComp(:,end);
indPri = find(strcmp(sltComp(:,3), 'Primary'))';
for i = indPri
    genes1 = sltComp{i,end};
    for j = indPri(i+1:end)
        genes2 = sltComp{j,end};
        C1 = setdiff(genes1,genes2,'stable');
        sltComp{i,end} = C1;
        C2 = setdiff(genes2,genes1,'stable');
        sltComp{j,end} = C2;
    end
end

timeComp = timeKLoop + toc;
fprintf('%d compartments identified, %3.0f minutes elapsed...\n',size(sltComp,1),timeComp/60)
clear C1 C2 genes1 genes2 i j sparsity
clear ans ff flag i I j jj match scoreCur
clear numfactors res x y ind1 ind2 indPri numOvlp
clear str1 str2 l1 l2 m mm n nn
clear lowThre highThre factorType
clear resPoss1 resPoss2 resPoss3
clear points resLen resScore resInd
clear II III maxIdx tmpIdx minSt
clear endFlag factorInd genesigs
clear genes1 genes2 indTmp indTmp1 indTmp2
clear factorLink factorLinkScore factorScore
clear genesigsAll markerGenesAll topGenesAll samplesigsAll sparsityAll
clear sltFactorLink sltFactorLinkScore

%% ========== Annotate final compartments =================================
tic
disp('Compartment annotating by MSigDB...')
for i=1:size(sltComp,1)
    fprintf('Compartment annotation: compartment %s...\n',sltComp{i,1}); 
    switch dataType
	case 'ATACseq'
	    disp('Gene annotation for ATAC-seq data not supported yet...');
	otherwise    
	    [~,geneOrder] = sort(sltComp{i,6},'descend');
	    GSEAout = GSEA(geneID(geneOrder),geneIDType);
	    sltComp(i,10) = {GSEAout};
	    compLabel = GSEAout.GS_name{1};
	    compLabel((strfind(compLabel,'_')))='-';
	    sltComp(i,11) = {compLabel};  
	    fprintf('Compartment annotation: top MSigDB term determined as %s...\n',compLabel);
	    clear geneOrder GSEAout topGS compLabel  
    end
end
timeGSEA = timeComp + toc;
fprintf('Compartment annotation completed, %3.0f seconds elapsed...\n',timeGSEA/60)

%% ========== Plot factor tree ============================================
tic
disp('Plotting factor tree...')
% Prepare tree nodes
[treeNodes,treeInd] = Gen_tree_nodes(sltFactorLinkFull,factorInfo,factorMatch);
clear factorMatch sltFactorLinkFull sltFactorLinkFullScore

% plot
Plot_factor_tree(treeNodes,treeInd,sltComp)

timePlotTree = timeGSEA + toc;
fprintf('Factor tree plotting completed, %3.0f seconds elapsed...\n',timeGSEA/60)
clear treeNodes treeInd
clear factorMatch factorInfo sltFactorLinkFull

%% ========== Write results ===============================================
disp('Deconvoluted results saving...')
tic
% compartment info
disp('Writing final compartment information...')
fid = fopen(fullfile(outDir,'Final_compartment_info.txt'),'wt');
fprintf(fid,'ID\tScore\tType\tTopGSEA\n');
for j = 1:size(sltComp,1)
    fprintf(fid,'D%s\t',num2str(sltComp{j,1}));
    fprintf(fid,'%s\t',num2str(sltComp{j,2}));
    fprintf(fid,'%s\t',sltComp{j,3});
    fprintf(fid,'%s\n',sltComp{j,11});
end
fclose(fid);

header = strcat(sltComp(:,1),'_',sltComp(:,11))';
% write sample weights
disp('Writing final sample weights...')
fid = fopen(fullfile(outDir,'Final_sample_weights.txt'),'wt');
fprintf(fid,'%s\t','sampleID');
fprintf(fid,'%s\t',header{1,1:end-1});
fprintf(fid,'%s\n',header{1,end});
samplesigs = cell2mat(sltComp(:,4));
samplesigs = samplesigs';
for j = 1:size(samplesigs,1)
    fprintf(fid,'%s\t',sampleID{1,j});
    fprintf(fid,'%d\t',samplesigs(j,1:end-1));
    fprintf(fid,'%d\n',samplesigs(j,end));
end
fclose(fid); 

% write gene weights
disp('Writing final gene weights...')
fid = fopen(fullfile(outDir,'Final_gene_weights.txt'),'wt');
fprintf(fid,'%s\t','geneID');
fprintf(fid,'%s\t',header{1,1:end-1});
fprintf(fid,'%s\n',header{1,end});
genesigs = cell2mat(sltComp(:,5)');
for j = 1:size(genesigs,1)
    fprintf(fid,'%s\t',geneID{j,1});
    fprintf(fid,'%d\t',genesigs(j,1:end-1));
    fprintf(fid,'%d\n',genesigs(j,end));
end
fclose(fid); 

% write gene scores
disp('Writing final gene scores...')
fid = fopen(fullfile(outDir,'Final_gene_scores.txt'),'wt');
fprintf(fid,'%s\t','geneID');
fprintf(fid,'%s\t',header{1,1:end-1});
fprintf(fid,'%s\n',header{1,end});
genesigs = cell2mat(sltComp(:,6)');
for j = 1:size(genesigs,1)
    fprintf(fid,'%s\t',geneID{j,1});
    fprintf(fid,'%d\t',genesigs(j,1:end-1));
    fprintf(fid,'%d\n',genesigs(j,end));
end
fclose(fid); 

% write top250 genes
disp('Writing final top 250 genes...')
topGenes = cell(251, size(sltComp,1));
for i = 1:size(sltComp,1)
    topGenes{1,i} = header{1,i};
    tmpGenes = sltComp{i,7};
    topGenes(2:size(tmpGenes,1)+1,i) = tmpGenes(:,1);
end
fid = fopen(fullfile(outDir,'Final_top250_genes.txt'),'wt');
for j = 1:size(topGenes,1)
    fprintf(fid,'%s\t',topGenes{j,1:end-1});
    fprintf(fid,'%s\n',topGenes{j,end});
end
fclose(fid);

% write marker genes
disp('Writing final marker genes...')
markerGenes = {};
for i = 1:size(sltComp,1)
    markerGenes{1,i} = header{1,i};
    tmpGenes = sltComp{i,9};
    markerGenes(2:size(tmpGenes,1)+1,i) = tmpGenes(:,1);
end
fid = fopen(fullfile(outDir,'Final_marker_genes.txt'),'wt');
for j = 1:size(markerGenes,1)
    fprintf(fid,'%s\t',markerGenes{j,1:end-1});
    fprintf(fid,'%s\n',markerGenes{j,end});
end
fclose(fid);

% MSigDB
disp('Writing MSigDB gene sets annoations...')
for i=1:size(sltComp,1)
    GSEAout = sltComp{i,10};
    fname = fullfile(outDir,sprintf('Final_GSEA_D%s.txt',num2str(sltComp{i,1})));
    fid = fopen(fname,'wt');
    rows =size(GSEAout.pvals,1);
    fprintf(fid,'pValue\ttopAUC\tKSstat\tGSname\n');
    [~,GSorder] = sort(GSEAout.topAUC,'descend');
    for j=GSorder(:)'
        fprintf(fid,'%12.10f\t%6.4f\t%6.4f\t%s\n',...
           GSEAout.pvals(j),GSEAout.topAUC(j),GSEAout.KSstat(j),GSEAout.GS_name{j});
    end
    fclose(fid);
end

clear i j topGenes markerGenes fid fname rows 
clear GSorder header tmpGenes genesigs samplesigs
timeTotal = timePlotTree + toc;
fprintf('De novo deconvolution completed, %3.0f minutes elapsed in total.\n',timeTotal/60)
