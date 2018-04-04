function [GeneSymbol_SNV,Entrez_Gene_Id,...
    Variant_Classification,Sample_ID_SNV] = P01_load_data_SNV(filename, startRow, endRow)

delimiter = '\t';
if nargin<=1
    startRow = 1;
    endRow = inf;
end

formatSpec = '%s%s%*s%*s%*s%*s%*s%*s%s%*s%*s%*s%*s%*s%*s%s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%[^\n\r]';

fileID = fopen(filename,'r');

dataArray = textscan(fileID, formatSpec, endRow(1)-startRow(1)+1, 'Delimiter', delimiter, 'EmptyValue' ,NaN,'HeaderLines', startRow(1)-1, 'ReturnOnError', false);
for block=2:length(startRow)
    frewind(fileID);
    dataArrayBlock = textscan(fileID, formatSpec, endRow(block)-startRow(block)+1, 'Delimiter', delimiter, 'EmptyValue' ,NaN,'HeaderLines', startRow(block)-1, 'ReturnOnError', false);
    for col=1:length(dataArray)
        dataArray{col} = [dataArray{col};dataArrayBlock{col}];
    end
end

fclose(fileID);

GeneSymbol_SNV = dataArray{:, 1};
Entrez_Gene_Id = dataArray{:, 2};
Variant_Classification = dataArray{:, 3};
Sample_ID_SNV = dataArray{:, 4};

test_str = dataArray{1,1}(1);
if strcmp(test_str{1}(1),'#')
    GeneSymbol_SNV(1:2) = [];
    Entrez_Gene_Id(1:2) = [];
    Variant_Classification(1:2) = [];
    Sample_ID_SNV(1:2) = [];
elseif strcmp(test_str,'Hugo_Symbol')
    GeneSymbol_SNV(1) = [];
    Entrez_Gene_Id(1) = [];
    Variant_Classification(1) = [];
    Sample_ID_SNV(1) = [];
end

end

