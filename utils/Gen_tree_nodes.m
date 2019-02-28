function [treeNodes,treeInd] = Gen_tree_nodes(sltFactorLinkFull,factorInfo,factorMatch)

% record involved factors
n = 0;
treeFactorList = cell(1,1);
for i = 1:size(sltFactorLinkFull,1)
    for j = 1:size(sltFactorLinkFull,2)
        if sltFactorLinkFull{i,j} == 0
            continue
        elseif isempty(sltFactorLinkFull{i,j})
            continue
        else
            numfactors = factorInfo(j,1);
            factorInd = sltFactorLinkFull{i,j};
            n = n+1;
            curFactor = sprintf('%d.%d',numfactors,factorInd);
            treeFactorList(n,1) = {curFactor};
        end
    end
end
[~,ia,~] = unique(treeFactorList(:,1));
treeFactorList = treeFactorList(ia,:);

% assign tree index to factors
ind = 0;
treeInd = cell(size(treeFactorList,1),1);
for i = 1:size(factorInfo,1)
    numfactors = factorInfo(i,1);
    for j = 1:numfactors
        curFactor = sprintf('%d.%d',numfactors,j);
        facIdx = find(strcmp(treeFactorList(:,1),curFactor));
        if ~isempty(facIdx)
            ind = ind+1;
            treeInd(ind,1) = {curFactor};
        end
    end
end

% generate node info
treeNodes = zeros(1,size(treeFactorList,1));
for i = 1:size(factorInfo,1)
    numfactors = factorInfo(i,1);    
    if i > 1
        matches = factorMatch(i,1:numfactors);  
    end
    for j = 1:numfactors
        queryName = sprintf('%d.%d',numfactors,j);
        nodeInd = find(cellfun(@(x)strcmp(x,queryName), treeInd(:,1)));
        if i == 1 % for the root
            treeNodes(nodeInd) = 0;
        else
            numfactorsPrev = factorInfo(i-1,1);
            if matches(j) == 0
                nodeName = 0;
            else
                queryName = sprintf('%d.%d',numfactorsPrev,matches(j));
                nodeName = find(cellfun(@(x)strcmp(x,queryName), treeInd(:,1)));
            end
            treeNodes(nodeInd) = nodeName;
        end 
    end
end