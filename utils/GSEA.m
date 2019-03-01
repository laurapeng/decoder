function [GSEAout] = GSEA(rankedGenes,geneIDType)

binDECODER = fileparts(mfilename('fullpath'));
binDECODER = fileparts(binDECODER);
addpath(fullfile(binDECODER,'data'))

switch geneIDType
    case 'TCGA'
        load('msigdb.v3.1.entrez.mat')
        geneIDsplit = split(rankedGenes,'|');
        rankedGenes = geneIDsplit(1:end,2);
        rankedGenes = cellfun(@str2double,rankedGenes);
    case 'EntrezID'
        load('msigdb.v3.1.entrez.mat')
        if strcmp(class(rankedGenes),'cell')
            rankedGenes = cellfun(@str2double,rankedGenes);
        end
    case 'geneSymbol'
        load('msigdb.v3.1.symbols.mat')
end

% unify
GS_universe = unique(vertcat(GS_list{:}));
rankedGenes = rankedGenes(ismember(rankedGenes,GS_universe));

% initialize
numtests = size(GS_list,1);
GSEAout.pvals = ones(numtests,1);
GSEAout.KSstat = ones(numtests,1);
GSEAout.topAUC = zeros(numtests,1);
GSEAout.significant = false(numtests,1);

p_crit = (0.05/(numtests));

numgenes = size(rankedGenes,1);
topgenes = round(numgenes/10);
x0 = (1:numgenes)';
F0 = x0/numgenes;

hit = zeros(length(rankedGenes),1);

% check
fprintf('Compartment annotation: checking %5.0f gene sets (MSigDB v3.1)...\n',numtests);
for Gset = 1:numtests
	current_list = GS_list{Gset};
    	hit = ismember(rankedGenes,current_list);
	npos = sum(hit);
	if npos>=8 % Enough hits in the gene set to calculate stats
		[h,p,ksstat,cv]=kstest(x0(hit>0),[x0,F0],p_crit,'unequal');
		Fn = cumsum(hit/npos);
		GSEAout.topAUC(Gset) = sum(Fn(1:topgenes))/topgenes;
		GSEAout.pvals(Gset)  = p;
		GSEAout.KSstat(Gset) = ksstat;
    end
    %fprintf('Compartment annotation: %5.0f gene sets checked, %5.0f remaining...\n',Gset,numtests-Gset);
end



% Benjamin-Hochberg procedure to determine significance
[pvals,rankorder] = sort(GSEAout.pvals);
significant = false(numtests,1);
for k = 1:numtests
	p_crit = (k*0.05/(numtests));
	if pvals(k)<p_crit
		significant(k) = true;
	else
		break
	end
end
GSEAout.significant(rankorder) = significant;

% Rankd and get only top 50
[~,GSrank] = sort(GSEAout.topAUC.*GSEAout.significant,'descend');
GSEAout.topAUC = GSEAout.topAUC(GSrank(1:50));
GSEAout.pvals = GSEAout.pvals(GSrank(1:50));
GSEAout.KSstat = GSEAout.KSstat(GSrank(1:50));
GSEAout.significant = significant(GSrank(1:50));
GSEAout.GS_name = GS_name(GSrank(1:50));
