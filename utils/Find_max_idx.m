function [resIdx,resLen,resScore] = Find_max_idx(startIdx,endIdx,highThre,scoreCur)

for j = endIdx:-1:startIdx
    if scoreCur(j) > highThre
        continue
    else
        break
    end
end

if endIdx~=j
    resLen = endIdx-startIdx+1;
    resScore = sum(scoreCur(j+1:endIdx))/resLen;
    tmpIdx = find(scoreCur(j+1:endIdx) == max(scoreCur(j+1:endIdx)));
    maxIdx = tmpIdx(1) + j;
    if maxIdx == startIdx && resLen <= 5  
        resIdx = endIdx;
    else
        resIdx = maxIdx;
    end
else
    resLen = 1;
    resScore = scoreCur(endIdx);
    resIdx = endIdx;
end
    