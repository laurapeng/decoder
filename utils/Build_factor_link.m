function [factorLink,factorLinkScore] = Build_factor_link(factorInfo,factorMatch,factorScore)

ff = size(factorInfo,1);
factorLink = {};
for i = 1:ff
    numfactors = factorInfo(i,1);
    factorLink{1,i} = sprintf('%d-factor',numfactors);
    if i == 1
        for j = 0:numfactors
            factorLink{j+2,i} = j;
        end
    else
        % the 0th factor
        indPrev = find(cellfun(@(x)isequal(x,0), factorLink(:,i-1)));
        factorLink{indPrev,i} = 0;
        % other factors
        match = factorMatch(i,1:numfactors);
        for j = 1:numfactors
            indPrev = find(cellfun(@(x)isequal(x,match(j)), factorLink(:,i-1)));
            if length(indPrev) > 1 % if already more than 2 dup 
                nrowCur = size(factorLink,1);
                % copy the rest to lower rows
                factorLink((indPrev(end)+2):nrowCur+1,1:i) ...
                    = factorLink((indPrev(end)+1):nrowCur,1:i);
                % copy the dup for one more line
                factorLink((indPrev(end)+1),1:i) ...
                    = factorLink(indPrev(end),1:i);
                factorLink{indPrev(end)+1,i} = j;
            elseif isempty(factorLink{indPrev,i}) % if first time
                factorLink{indPrev,i} = j;
            else % if already has one
                nrowCur = size(factorLink,1);
                factorLink((indPrev+2):nrowCur+1,1:i) ...
                    = factorLink((indPrev+1):nrowCur,1:i);
                factorLink(indPrev+1,1:i) ...
                    = factorLink(indPrev,1:i);
                factorLink{indPrev+1,i} = j;
            end
        end
    end
end

factorLinkScore = cell(size(factorLink));
for i = 1:size(factorLink,2)
    factorLinkScore(1,i) = factorLink(1,i);
    for j = 2:size(factorLink,1)
        if i == 1
            factorLinkScore{j,i} = 1;
        elseif factorLink{j,i} == 0
            factorLinkScore{j,i} = 1;
        else
            if ~isempty(factorLink{j,i})
                score = factorScore(i,factorLink{j,i});
                factorLinkScore{j,i} = round(score,4);
            end
        end
        
    end
end