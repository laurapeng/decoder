function [sltComp] = Proc_comp_candidates(sltFactorLink,topGenesAll,sltComp)

% find link that's covered by another link
for m = 1:size(sltComp,1)
    l1 = cell2mat(sltFactorLink(m,:));
    for mm = 1:size(sltComp,1)
        if m == mm
            continue
        end
        l2 = cell2mat(sltFactorLink(mm,:));
        if length(l1) < length(l2)
            if (sum(l2(1:length(l1)) == l1) == length(l1))
                if strcmp(sltComp{m,3},'Major')
                    genes1 = topGenesAll{length(l1),l1(end)};
                    genes2 = topGenesAll{length(l2),l2(end)};
                    numOvlp = size(intersect(genes1,genes2),1);
                    if strcmp(sltComp{mm,3},'Major')...
                            && (numOvlp > 100)
                        sltComp{mm,3} = 'Dropped';
                    else 
                        sltComp{mm,3} = 'Minor';
                    end
                elseif strcmp(sltComp{m,3},'Unstable')
                    if strcmp(sltComp{mm,3},'Major')
                        sltComp{m,3} = 'Dropped';
                    elseif strcmp(sltComp{mm,3},'Unstable')
                        sltComp{mm,3} = 'Dropped';
                    end
                end
            end
        end
    end
end

% find 'Unstable' and 'Transient' that are the same with 'Primary'
indTmp1 = find(strcmp(sltComp(:,3),'Unstable'));
indTmp = sort(indTmp1);
indPri = find(strcmp(sltComp(:,3),'Major'));
for i = 1:length(indTmp)
    ind1 = indTmp(i);
    l1 = cell2mat(sltFactorLink(ind1,:));
    genes1 = topGenesAll{length(l1),l1(end)};
    for j = 1:length(indPri)
        ind2 = indPri(j);
        l2 = cell2mat(sltFactorLink(ind2,:));
        genes2 = topGenesAll{length(l2),l2(end)};
        numOvlp = size(intersect(genes1,genes2),1);
        if numOvlp > 100
            sltComp{ind1,3} = 'Dropped';
        end
    end
end

sltComp(strcmp(sltComp(:, 3), 'Dropped'), :) = [];
sltComp = sortrows(sltComp,2,'descend');
% sort by compt type
sltComp(find(strcmp(sltComp(:,3),'Major')),4) = {1};
sltComp(find(strcmp(sltComp(:,3),'Minor')),4) = {2};
sltComp(find(strcmp(sltComp(:,3),'Unstable')),4) = {3};
sltComp = sortrows(sltComp,4);
sltComp(:,4) = [];
