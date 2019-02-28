function [resIdx,nextFlag] = Pick_comp_candidates(factorLinkScore,i,lowThre,highThre)

nextFlag = 0;
factorType = 'Dropped';
scoreCur = factorLinkScore(i,:);
scoreCur(cell2mat(cellfun(@(elem) elem == 1, scoreCur(:, :), 'UniformOutput', false))) = {0};
scoreCur = cell2mat(scoreCur);
if sum(scoreCur) == 0
    resIdx = 0;
    nextFlag = 1;
    return
end
if size(find(scoreCur ~= 0),2) <=3
    resIdx = 0;
    nextFlag = 1;
    return
end    
[~,I] = sort(scoreCur,'descend');
if scoreCur(I(1)) <= highThre
    resIdx = 0;
    nextFlag = 1;
    return
end

% walk forward to build peak blocks
points = [];
resIdx = [];
resLen = [];
resScore = [];
m = 0;
endFlag = 0;

% detect the first low point p1
for j = 1:length(scoreCur)
    if scoreCur(j) == 0
        continue
    elseif scoreCur(j) > highThre
        continue
    else
        break
    end
end
points(1) = j;
resPoss1 = 0;

% if points(1) is the end
if points(1) == length(scoreCur)
    resPoss1 = 1;
% if points(1) is not the end
else
    % p1: if drop to low before above high
    % p1: or stay medium before above high
    for j = points(1):length(scoreCur)
        if scoreCur(j) <= lowThre
            resPoss1 = 1;
        elseif scoreCur(j) <= highThre
            continue
        else
            break
        end
    end
    if  j-points(1)>=2
        resPoss1 = 1;
    end

    % detect the next high point p2
    for j = points(1):length(scoreCur)
        if scoreCur(j) <= highThre
            continue
        else
            break
        end
    end
    points(2) = j;

    % if p2 is the end
    if points(2) == length(scoreCur)
        resPoss1 = 1;
    end
end

if resPoss1 == 1
    for j = points(1):-1:1
        if scoreCur(j) == 0
            break
        elseif scoreCur(j) <= highThre
            continue
        else
            break
        end
    end
    if scoreCur(j) ~= 0
        m=m+1;
        [resIdx(m),resLen(m),resScore(m)] = Find_max_idx(2,j,highThre,scoreCur);
    end    
end

% if p1 not the end
if points(1) ~= length(scoreCur)
    % if p2 not the end
    if  points(2) ~= length(scoreCur)
        % detect the next medium/low point p3 or end
        for j = points(2):length(scoreCur)
            if scoreCur(j) > highThre
                continue
            else
                break
            end
        end

        points(3) = j;
        resPoss2 = 0;
        % if points(3) is the end
        if points(3) == length(scoreCur)
            resPoss2 = 1;
        % if p3 not the end
        else
            % p3: if drop to low before above high
            % p3: or stay medium before above high
            for j = points(3):length(scoreCur)
                if scoreCur(j) <= lowThre
                    resPoss2 = 1;
                elseif scoreCur(j) <= highThre
                    continue
                else
                    break
                end
            end
            if  j-points(3)>=2
                resPoss2 = 1;
            end

            % detect the next high point p4
            for j = points(3):length(scoreCur)
                if scoreCur(j) <= highThre
                    continue
                else
                    break
                end
            end
            points(4) = j;

            % if p4 is the end
            if points(4) == length(scoreCur)
                resPoss2 = 1;
            end
        end

        if resPoss2 == 1
            for j = points(3):-1:points(2)
                if scoreCur(j) <= highThre
                    continue
                else
                    break
                end
            end
            m = m +1;
            if resPoss1 == 1
                [resIdx(m),resLen(m),resScore(m)] = Find_max_idx(points(2),j,highThre,scoreCur);
            else
                [resIdx(m),resLen(m),resScore(m)] = Find_max_idx(2,j,highThre,scoreCur);
            end
        end

        % if p3 not the end
        if  points(3) ~= length(scoreCur)
            % if p4 not the end
            if  points(4) ~= length(scoreCur)
                % detect the next medium/low point p5 or end
                for j = points(4):length(scoreCur)
                    if scoreCur(j) > highThre
                        continue
                    else
                        break
                    end
                end

                points(5) = j;
                resPoss3 = 0;
                % if points(5) is the end
                if points(5) == length(scoreCur)
                    resPoss3 = 1;
                % if p5 not the end
                else
                    % p5: if drop to low before above high
                    % p5: or stay medium before above high
                    for j = points(5):length(scoreCur)
                        if scoreCur(j) <= lowThre
                            resPoss3 = 1;
                        elseif scoreCur(j) <= highThre
                            continue
                        else
                            break
                        end
                    end
                    if  j-points(5)>=2
                        resPoss3 = 1;
                    end

                    % detect the next high point p6 or end
                    for j = points(5):length(scoreCur)
                        if scoreCur(j) <= highThre
                            continue
                        else
                            break
                        end
                    end
                    points(6) = j;

                    % if p6 is the end
                    if points(6) == length(scoreCur)
                        resPoss3 = 1;
                    end
                end

                if resPoss3 == 1
                    for j = points(5):-1:points(4)
                        if scoreCur(j) <= highThre
                            continue
                        else
                            break
                        end
                    end
                    m = m+1;
                    if resPoss2 == 1
                        [resIdx(m),resLen(m),resScore(m)] = Find_max_idx(points(4),j,highThre,scoreCur);
                    elseif resPoss1 == 1
                        [resIdx(m),resLen(m),resScore(m)] = Find_max_idx(points(2),j,highThre,scoreCur);
                    else
                        [resIdx(m),resLen(m),resScore(m)] = Find_max_idx(2,j,highThre,scoreCur);
                    end
                end

                % if p5 not the end     
                if points(5) ~= length(scoreCur)
                    % if p6 not the end     
                    if points(6) ~= length(scoreCur)
                        for j = points(6):length(scoreCur)
                            if scoreCur(j) > highThre
                                continue
                            else
                                break
                            end
                        end

                        points(7) = j;
                        resPoss4 = 0;
                        if points(7) == length(scoreCur)
                            resPoss4 = 1;
                        else
                            for j = points(7):length(scoreCur)
                                if scoreCur(j) <= lowThre
                                    resPoss4 = 1;
                                elseif scoreCur(j) <= highThre
                                    continue
                                else
                                    break
                                end
                            end
                            if  j-points(7)>=2
                                resPoss4 = 1;
                            end

                            % detect the next high point p8 or end
                            for j = points(7):length(scoreCur)
                                if scoreCur(j) <= highThre
                                    continue
                                else
                                    break
                                end
                            end
                            points(8) = j;
                            % if p8 is the end
                            if points(8) == length(scoreCur)
                                resPoss4 = 1;
                            end
                        end

                        if resPoss4 == 1
                            for j = points(7):-1:points(6)
                                if scoreCur(j) <= highThre
                                    continue
                                else
                                    break
                                end
                            end
                            m = m+1;
                            if resPoss3 == 1
                                [resIdx(m),resLen(m),resScore(m)] = Find_max_idx(points(6),j,highThre,scoreCur);
                            elseif resPoss2 == 1
                                [resIdx(m),resLen(m),resScore(m)] = Find_max_idx(points(4),j,highThre,scoreCur);
                            elseif resPoss1 == 1
                                [resIdx(m),resLen(m),resScore(m)] = Find_max_idx(points(2),j,highThre,scoreCur);
                            else
                                [resIdx(m),resLen(m),resScore(m)] = Find_max_idx(2,j,highThre,scoreCur);
                            end
                        end

                        if points(7) ~= length(scoreCur)
                            if points(8) ~= length(scoreCur)
                               for j = points(8):length(scoreCur)
                                    if scoreCur(j) > highThre
                                        continue
                                    else
                                        break
                                    end
                               end 

                                if j~= length(scoreCur)
                                    j = j-1;
                                else
                                    endFlag = 1;
                                end
                                m = m+1;   
                                if resPoss4 == 1
                                    [resIdx(m),resLen(m),resScore(m)] = Find_max_idx(points(8),j,highThre,scoreCur);
                                elseif resPoss3 == 1
                                    [resIdx(m),resLen(m),resScore(m)] = Find_max_idx(points(6),j,highThre,scoreCur);
                                 elseif resPoss2 == 1
                                    [resIdx(m),resLen(m),resScore(m)] = Find_max_idx(points(4),j,highThre,scoreCur);
                                elseif resPoss1 == 1
                                    [resIdx(m),resLen(m),resScore(m)] = Find_max_idx(points(2),j,highThre,scoreCur);
                                else
                                    [resIdx(m),resLen(m),resScore(m)] = Find_max_idx(points(2),points(8),highThre,scoreCur);
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end

if isempty(resIdx)
    resIdx = 0;
    nextFlag = 1;
    return
end

[~,III] = sort(resLen,'descend');
% lasting in less than 3 runs
if resLen(III(1)) < 3
    factorType = 'Unstable';
end

if strcmp(factorType,'Dropped')
    factorType = 'Primary';
end

if length(resIdx) > 1
    if (resLen(III(1)) == resLen(III(2)))...
            && (resIdx(III(1)) < resIdx(III(2)))
        resIdx = resIdx(III(2));
    elseif (resLen(III(1)) >= resLen(III(2))+1)...
            && (resIdx(III(1)) < resIdx(III(2)))...
            && (resLen(III(2)) >= 3)...
            && (endFlag == 1)
        resIdx = resIdx(III(2));
    else
        resIdx = resIdx(III(1));
    end
end