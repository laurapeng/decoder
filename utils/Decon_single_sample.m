function [sampleWeights] = Decon_single_sample(configFile)
%function [sampleWeights] = Decon_single_sample(refSet,dataMatrix,dataFormat,geneIDType,logTransformed)


% Add path
binDECODER = fileparts(mfilename('fullpath'));
binDECODER = fileparts(binDECODER);
addpath(binDECODER)
addpath(fullfile(binDECODER,'data'))
addpath(fullfile(binDECODER,'utils'))

% Load reference
load(sprintf('%s.mat',refSet));
indPri = find(strcmp(sltComp(:,3),'Primary'));
sltComp = sltComp(indPri,:);
geneSigRef = cell2mat(sltComp(:,5)');

%%% if dataset too large
tmpFile = dir(fullfile(binDECODER,'data',sprintf('%s.*.mat',refSet)));
tmpFile = {tmpFile(:).name}';
if ~isempty(tmpFile)
    load(fullfile(binDECODER,'data',sprintf('%s.1.mat',refSet)))
    dataMerge = data;
    for i=2:size(tmpFile,1)
        load(fullfile(binDECODER,'data',sprintf('%s.%d.mat',refSet,i)))
        dataMerge = [dataMerge data];
    end
    data = dataMerge;
end

dataRef = log2(1+data);
switch geneIDType
    case 'TCGA'
        geneIDRef = geneID;
    case 'EntrezID'
        geneID = split(geneID,'|');
        geneIDRef = geneID(1:end,2);
    case 'geneSymbol'
        geneID = split(geneID,'|');
        geneIDRef = geneID(1:end,1);
end
fprintf('Reference %s loaded...\n',refSet);
clear data geneID geneSymbol EntrezID sampleID

% Load data for deconvolution
[data,geneID,sampleID] = Load_data(dataMatrix,dataFormat);
if strcmp(logTransformed,'no')
    data = log2(1+data);
end
[outDir,outPrefix,~] = fileparts(dataMatrix);

% Overlap genes
[~,ia,ib] = intersect(geneIDRef,geneID,'stable');
geneIDRef = geneIDRef(ia,:);
dataRef = dataRef(ia,:);
geneSigRef = geneSigRef(ia,:);
geneID = geneID(ib,:);
data = data(ib,:);

% Normalize data
maxRef = max(dataRef(:));
data = (maxRef/max(data(:)))*data;

% Calculate sample weights
sampleWeights = zeros(size(geneSigRef,2),size(data,2));
for i = 1:size(sampleWeights,2) % for each sample
    options = optimset('TolX',0.001);
     [sampleWeights(1:size(geneSigRef,2),i),~] = lsqnonneg(geneSigRef, data(:,i),options );
end

% Write results
outMat = sprintf('%s.sample_weights.mat',outPrefix);
save(fullfile(outDir,outMat),'sampleWeights','sampleID','sltComp')

outName = sprintf('%s.sample_weights.txt',outPrefix);
fid = fopen(fullfile(outDir,outName),'wt');
header = sltComp(:,11)';
fprintf(fid,'%s\t','sampleID');
fprintf(fid,'%s\t',header{1,1:end-1});
fprintf(fid,'%s\n',header{1,end});
samplesigs = sampleWeights';
for j = 1:size(samplesigs,1)
    fprintf(fid,'%s\t',sampleID{1,j});
    fprintf(fid,'%s',sprintf('%.3f\t',samplesigs(j,1:end-1)));
    fprintf(fid,'%s',sprintf('%.3f\n',samplesigs(j,end)));
end
fclose(fid); 

