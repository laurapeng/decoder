function [sampleWeights] = Decon_single_sample(configFileName)
%function [sampleWeights] = Decon_single_sample(refSet,dataMatrix,dataFormat,geneIDType,logTransformed)

outDir = pwd;

% Add path
binDECODER = fileparts(mfilename('fullpath'));
binDECODER = fileparts(binDECODER);
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

% Load reference
tmpInd = find(strcmp(configInfo{1},'refSet'));
refSet = configInfo{2}{tmpInd};
load(sprintf('%s.mat',refSet));
indPri = find(strcmp(sltComp(:,3),'Primary'));
sltComp = sltComp(indPri,:);
geneSigRef = cell2mat(sltComp(:,5)');

%%% if dataset too large
tmpFile = dir(fullfile(binDECODER,'data','*.mat'));
tmpFile = {tmpFile(:).name}';
for k = 1:size(tmpFile,1)
    if contains(tmpFile{k,1},refSet)
        tmpMatch = regexp(tmpFile{k,1},'(\d+)','match','once');
        
        if isempty(tmpMatch)
            continue
        elseif strcmp(tmpMatch,'1')
            load(tmpFile{k,1})
            %load(fullfile(binDECODER,'data',sprintf('%s.1.mat',refSet)))
            dataMerge = data;
        else
            load(tmpFile{k,1})
            %load(fullfile(binDECODER,'data',sprintf('%s.%d.mat',refSet,i)))
            dataMerge = [dataMerge data];
        end
    end
end

if exist('dataMerge')
    data = dataMerge;
end
dataRef = log2(1+data);

% parse geneID
tmpInd = find(strcmp(configInfo{1},'geneIDType'));
geneIDType = configInfo{2}{tmpInd};

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
tmpInd = find(strcmp(configInfo{1},'dataFormat'));
dataFormat = configInfo{2}{tmpInd};

tmpInd = find(strcmp(configInfo{1},'dataMatrix'));
dataMatrix = configInfo{2}{tmpInd};
[outDir,outPrefix,~] = fileparts(dataMatrix);

[data,geneID,sampleID] = Load_data(dataMatrix,dataFormat);

% log transformation
tmpInd = find(strcmp(configInfo{1},'logTransformed'));
logTransformed = configInfo{2}{tmpInd};

if strcmp(logTransformed,'no')
    data = log2(1+data);
end

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

