function [factorScore,factorMatch] = Link_previous_factor(factorInd,numfactors,factorInfo,topGenesAll,factorScore,factorMatch)

for x = 1:numfactors
    switch factorInd
        case 1
            matchedFactor = x;
            maxPct = 1;
        otherwise
            matchedFactor = 0;
            maxPct = 0;
            numfactorsPrev = factorInfo((factorInd-1),1);
            for y = 1:numfactorsPrev
                geneOvlp = size(intersect(topGenesAll{factorInd,x},topGenesAll{(factorInd-1),y}),1);
                genePct = geneOvlp/250;
                if(genePct > maxPct)
                    maxPct = genePct;
                    matchedFactor = y;
                end
            end
            if maxPct < 0.1
                matchedFactor = 0;
            end
    end  
    factorScore(factorInd,x) = maxPct;
    factorMatch(factorInd,x) = matchedFactor;
    clear matchedFactor maxPct numfactorsPrev geneOvlp genePct
end
