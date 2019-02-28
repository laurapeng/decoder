function [data,geneID,sampleID,master_subset] = Load_data(dataMatrix,dataFormat)

if exist(dataMatrix) ~= 2
    fprintf('Error: Data matrix %s not found !!!\n',dataMatrix);
    return
end

switch dataFormat
    case {'tsv','csv'}
        fprintf('Data matrix %s loading...\n',dataMatrix);
        rawTable = readtable(dataMatrix,'ReadVariableNames',false);
        geneID = table2cell(rawTable(2:end,1));
        sampleID = table2cell(rawTable(1,2:end));
        data = str2double(table2array(rawTable(2:end,2:end)));
        disp('Data matrix loaded...');
    case 'mat'
        fprintf('Data matrix %s loading...\n',dataMatrix);
        load(dataMatrix)
        disp('Data matrix loaded...');
    otherwise
        disp('Error: Data format not supported!!!')
        return
end

if ~exist(master_subset)
master_subset = logical(ones(size(data,2),1));
end

