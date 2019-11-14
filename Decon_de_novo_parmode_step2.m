function Decon_de_novo_parmode_step2(configFileName) 

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
tmpInd = find(strcmp(configInfo{1},'dataType'));
dataType = configInfo{2}{tmpInd};

% Load data matrix
tmpInd = find(strcmp(configInfo{1},'dataMatrix'));
dataMatrix = configInfo{2}{tmpInd};

tmpInd = find(strcmp(configInfo{1},'dataFormat'));
dataFormat = configInfo{2}{tmpInd};
[data,geneID,sampleID,master_subset] = Load_data(dataMatrix,dataFormat);

outDir = fileparts(dataMatrix);
timePreproc = toc;
fprintf('Preprocessing completed, %3.0f minutes elapsed...\n',timePreproc/60)

%% ========== load saved data from increasing K ==================================================
tic

tmpFile = dir(fullfile(outDir,'*-factor_data.mat'));
facDat = {tmpFile(:).name}';
tmpFile = regexp(facDat(:,1), '-', 'split');
tmpFile = vertcat(tmpFile{:});
tmpFile = tmpFile(:,1);
tmpFile = cellfun(@str2double,tmpFile);
[facNum,I] = sort(tmpFile);
facDat = facDat(I,:);
factorInfo = facNum;
clear tmpFile I facNum

% find exemplar genes 
ff = size(factorInfo,1);
topGenes = cell(ff,factorInfo(end,1));
topGenesAll = cell(ff,factorInfo(end,1));
sparsityAll = cell(ff,factorInfo(end,1));
genesigsAll = cell(ff,factorInfo(end,1));
samplesigsAll = cell(ff,factorInfo(end,1));
sampleScoreAll = cell(ff,factorInfo(end,1));

for i = 1:ff
    % calculate sparsity for current numfactors
    load(facDat{i,1});
    factorInfo(i,2) = coph;
    
    for x = 1:size(genesigs,2)
        samplesigsAll{i,x} = samplesigs(x,:);
        %sampleScore = zscore(samplesigs(x,:),0,2);
	sampleScore = samplesigs(x,:);
        sampleScoreAll{i,x} = sampleScore;
        
        genesigsAll{i,x} = genesigs(:,x);
        sparsity = genesigs(:,x)-max(genesigs(:,setdiff(1:size(genesigs,2),x)),[],2);
        sparsityAll{i,x} = sparsity;
        % find dynamic top ones for current numfactors
        [sparsitySrt,order] = sort(sparsity, 'descend');
        geneIDSrt = geneID(order,:);
        sparsitySrt(:,2) = 1:1:length(sparsity);
        sparsitySrtPos = sparsitySrt(sparsitySrt(:,1) > 0, :);
        sparsitySrtPos(:,3) = sparsitySrtPos(:,1)/max(sparsitySrtPos(:,1));
        sparsitySrtPos(:,4) = sparsitySrtPos(:,2)/max(sparsitySrtPos(:,2));
        sparsitySrtPos(:,5) = abs(sparsitySrtPos(:,4) - sparsitySrtPos(:,3));
        [minValue,minInd] = min(sparsitySrtPos(:,5));
        topGenesAll{i,x} = geneIDSrt(1:minInd,1);
        clear sparsitySrt order sparsitySrtPos geneIDSrt
        % top 250
	[sparsitySrt,order] = sort(sparsity, 'descend');
        top = sparsitySrt>quantile(sparsitySrt,1-250/size(sparsitySrt,1));
	geneIDSrt = geneID(order,1);
        topGenes{i,x} = geneIDSrt(top,1);
    end
end
clear i x minInd minValue numfactors ylabels top sparsity coph
clear samplesigs genesigs sampleScore geneIDSrt sparsitySrt

% find corresponding factors in each run
gg = size(topGenes,1);
resMatch = zeros(ff,factorInfo(end,1));
resPct = zeros(ff,factorInfo(end,1));

for i = 1:gg  
    % find corresponding factor and calculate overlap
    numfactors = length(find(~cellfun(@isempty,topGenes(i,:))));
    for x = 1:numfactors
        if i == 1
            matchedFactor = x;
            maxPct = 1;
        else
            matchedFactor = 0;
            maxPct = 0;
            
            numfactorsPrev = length(find(~cellfun(@isempty,topGenes(i-1,:))));
            for y = 1:numfactorsPrev
                geneOvlp = size(intersect(topGenes{i,x},topGenes{i-1,y}),1);
                genePct = geneOvlp/size(topGenes{i-1,y},1);
                %geneCmb = size(unique([topGenes{i,x};topGenes{i-1,y}]),1);
                %genePct = geneOvlp/geneCmb;
                if(genePct > maxPct)
                    maxPct = genePct;
                    matchedFactor = y;
                end
            end
        end
        if maxPct < 0.1
            matchedFactor = 0;
        end
        resPct(i,x) = maxPct;
        resMatch(i,x) = matchedFactor;
        clear matchedFactor macOvlp
    end
end

clear ff gg geneOvlp i x y
clear numfactors numfactorsPrev 
clear genePct maxPct matchedFactor

% determine factors to filter
factorInfo(1,3) = 1;
for i = 2:size(resPct,1)
    tmp = resPct(i,resPct(i,:) ~= 0);
    tmp2 = median(tmp);
    factorInfo(i,3) = tmp2;
end

factorStartInd = 1;
numSkip = 0;
for factorInd = 1:size(factorInfo,1)
    if ((factorInfo(factorInd,3) >= 0.5) && (factorInd ~= 1))
        if factorStartInd==1
            factorStartInd = factorInd;
        end
    end
    
    % if break
    if factorInfo(factorInd,1)-factorInfo(factorInd-1,1)>3
        factorEndInd = factorInd-1;
        break
    end
    
    factorEndInd = factorInd;
    if (factorStartInd~=1)
        if factorInd-factorStartInd+1 >= 9
            if (factorInfo(factorInd,3) < 0.5)...
                    && (factorInfo((factorInd-1),3) < 0.5)...
                    &&  (factorInfo((factorInd-2),3) < 0.5)...
                    &&  (factorInfo((factorInd-3),3) < 0.5)
                factorEndInd = factorInd-4;
                break
            end
        end
    end
end

factorInfoFlt = factorInfo(factorStartInd:factorEndInd,:);
topGenesFlt = topGenes(factorStartInd:factorEndInd,:);
resMatchFlt = resMatch(factorStartInd:factorEndInd,:);
resPctFlt = resPct(factorStartInd:factorEndInd,:);

timeKLoop = timePreproc + toc;
fprintf('Load multiple runs completed, %3.0f minutes elapsed...\n',timeKLoop/60)

clear i ix tmp tmp2

%% ========== Determine compartments from factors in multiple runs ========
tic
disp('Compartment selecting...')
% build linkages 
[factorLink,factorLinkScore] = Build_factor_link(factorInfo,factorMatch,factorScore);

% estimate threshold
[lowThre,highThre] = Est_score_threshold(factorLinkScore,repTimes);

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
        numfactors = factorInfo(resIdx,1);
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
indPri = find(strcmp(sltComp(:,3), 'Major'))';
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
