function Plot_factor_tree(treeNodes,treeInd,sltComp)

[x,y,~] = treelayout(treeNodes);
f = find(treeNodes~=0);
pp = treeNodes(f);
X = [x(f); x(pp); NaN(size(f))];
Y = [y(f); y(pp); NaN(size(f))];
X = X(:);
Y = Y(:);
n = length(treeNodes);

clf
plot(x, y, 'ko', X, Y, 'k-')
set(gca,'xtick',[])
set(gca,'ytick',[])
axis off
% re-color
hold on
tmpComp = find(strcmp('Major', sltComp(:,3)));
if length(tmpComp) ~= 0
    [~,ia,ib] = intersect(sltComp(tmpComp,1),treeInd(:,1));
    cmpIdxX = x(:,ib);
    cmpIdxy = y(:,ib);
    plot(cmpIdxX,cmpIdxy,'ro',...
        'MarkerSize', 9,...
        'LineWidth', 2,...
        'MarkerFaceColor', 'red',...
        'MarkerEdgeColor', 'red')
end

tmpComp = find(strcmp('Minor', sltComp(:,3)));
if length(tmpComp) ~= 0
    [~,ia,ib] = intersect(sltComp(tmpComp,1),treeInd(:,1));
    cmpIdxX = x(:,ib);
    cmpIdxy = y(:,ib);
    plot(cmpIdxX,cmpIdxy,'yo',...
        'MarkerSize', 9,...
        'LineWidth', 2,...
        'MarkerFaceColor', 'yellow',...
        'MarkerEdgeColor', 'yellow')
end

tmpComp = find(strcmp('Unstable', sltComp(:,3)));
if length(tmpComp) ~= 0
    [~,ia,ib] = intersect(sltComp(tmpComp,1),treeInd(:,1));
    cmpIdxX = x(:,ib);
    cmpIdxy = y(:,ib);
    plot(cmpIdxX,cmpIdxy,'bo',...
        'MarkerSize', 9,...
        'LineWidth', 2,...
        'MarkerFaceColor', 'blue',...
        'MarkerEdgeColor', 'blue')
end

% label
xx = x';
yy = y';
nodeName1 = string(treeInd(:,1));
text(xx(:,1), yy(:,1), nodeName1, 'VerticalAlignment','bottom','HorizontalAlignment','center','FontSize',10);
title({'DECODER factor tree'},'FontSize',15,'FontName','Times New Roman');
legend({'Filtered factor','Linkage','Major compartment','Unstable compartment','Minor compartment'},...
        'Location','northeastoutside','FontSize',10)
legend('boxoff')
print('-bestfit','Final_factor_tree','-dpdf')
close
