function [data,geneID,sampleID,master_subset] = Load_data(dataMatrix,dataFormat)

if exist(dataMatrix) ~= 2
    fprintf('Error: Data matrix %s not found !!!\n',dataMatrix);
    return
end

switch dataFormat
    case {'tsv','csv'}
        fprintf('Data matrix %s loading...\n',dataMatrix);
        rawTable = importdata(dataMatrix);
        geneID = rawTable.textdata(2:end,1);
        sampleID = rawTable.textdata(1,2:end);
        data = rawTable.data;
        disp('Data matrix loaded...');
    case 'mat'
        fprintf('Data matrix %s loading...\n',dataMatrix);
        load(dataMatrix)
        disp('Data matrix loaded...');
    otherwise
        disp('Error: Data format not supported!!!')
        return
end

master_subset = logical(ones(size(data,2),1));

