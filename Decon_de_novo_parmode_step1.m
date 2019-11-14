function Decon_de_novo_parmode_step1(configFileName,numfactors) 

%% ========== Preprocess ==================================================
tic

% Add path
tmpInd = find(strcmp(configInfo{1},'binDECODER'));
binDECODER = configInfo{2}{tmpInd};
addpath(binDECODER)
addpath(fullfile(binDECODER,'data'))
addpath(fullfile(binDECODER,'utils'))

% Read and parse configure file 
configFile = fopen(configFileName);
if configFile == -1
    disp('Configure file not found...');
else
    configInfo = textscan(configFile,'%s\t%s');
    %configInfo = [configInfo{:}];
end
fclose(configFile);

% Set parameters
tmpInd = find(strcmp(configInfo{1},'geneIDType'));
geneIDType = configInfo{2}{tmpInd};
tmpInd = find(strcmp(configInfo{1},'repTimes'));
repTimes = str2double(configInfo{2}{tmpInd});

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
        disp('Error: Data type not supported...');
        return
end

% Load data matrix
tmpInd = find(strcmp(configInfo{1},'dataMatrix'));
dataMatrix = configInfo{2}{tmpInd};

tmpInd = find(strcmp(configInfo{1},'dataFormat'));
dataFormat = configInfo{2}{tmpInd};
[data,geneID,sampleID,master_subset] = Load_data(dataMatrix,dataFormat);

tmpInd = find(strcmp(configInfo{1},'logTransformed'));
logTransformed = configInfo{2}{tmpInd};
if strcmp(logTransformed,'no')
    data = log2(1+data);
end
outDir = fileparts(dataMatrix);

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
clear dataFormat dataMatrix logTransformed numTrainGene
clear rawTable x tmpInd tmp ans


%%% ===== Train gene weight seed on 5K genes through repetitions ======
timeIteration = timeKLoop;
data = selected.data;
geneID = selected.geneID ;
rng(999);
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

	fprintf('K=%d: Converged on %d factors, thus skipped...\n',numfactors, tmpK)
	timeKLoop = timeIteration + toc;
	break
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

	fprintf('K=%d: Converged on %d factors, thus skipped...\n',numfactors, tmpK)
	timeKLoop = timeIteration + toc;
	break
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

timeKLoop = timeIteration + toc;
fprintf('K=%d: deconvolution completed, %3.0f minutes elapsed...\n',numfactors,timeKLoop/60)

