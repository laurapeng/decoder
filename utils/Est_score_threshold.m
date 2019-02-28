function [lowThre,highThre] = Est_score_threshold(factorLinkScore)

n = 0;
scoreList = [];
for i = 2:size(factorLinkScore,1)
    scoreCur = factorLinkScore(i,:);
    scoreCur(cell2mat(cellfun(@(elem) elem == 1, scoreCur(:, :), 'UniformOutput', false))) = {0};
    scoreCur = cell2mat(scoreCur);
    for j = 1:size(factorLinkScore,2)
        if factorLinkScore{i,j} ~= 1
            n = n+1;
            scoreList(n) = factorLinkScore{i,j};
        end
    end
end
lowThre = quantile(scoreList,0.25);
if lowThre < 0.25
    lowThre = 0.25;
end
highThre = median(scoreList);
