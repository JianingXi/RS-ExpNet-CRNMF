function [GeneSymbol_EXP,Sample_ID_EXP,X_EXP_raw,EntrezID_EXP] = P02_load_data_EXP(input_txt_file_str)

fid = fopen(input_txt_file_str,'r');
tline_1 = fgetl(fid);
tline_2 = fgetl(fid);
fclose(fid);

str_header_1 = 'Hugo_Symbol	Entrez_Gene_Id	Cytoband	';
str_header_2 = 'Hugo_Symbol	Entrez_Gene_Id	';
str_1 = tline_1(1:length(str_header_1));
str_2 = tline_1(1:length(str_header_2));

if strcmp(str_header_1,str_1)
    Col_off = 3;
elseif strcmp(str_header_2,str_2)
    Col_off = 2;
else
    Col_off = 1;
end

str_row = 'Composite Element REF';
if strcmp(tline_2(1:length(str_row)),str_row)
    Row_off = 2;
else
    Row_off = 1;
end

M = dlmread(input_txt_file_str,'\t',Row_off,Col_off);
[~,N_col] = size(M);
X_EXP_raw = M';
clear M

fid = fopen(input_txt_file_str,'r');
C_text = textscan(fid,'%s',N_col+Col_off,'Delimiter','\t');
fclose(fid);


% Header_vec = {'Hugo_Symbol','Entrez_Gene_Id','Cytoband'};
% ind_head = strcmp(Header_vec,C_text{1,1}(1:3));
% if sum(ind_head) == 2

Sample_ID_EXP = C_text{1,1}((Col_off+1):end);

fid = fopen(input_txt_file_str,'r');
str_format = [repmat('\t%s',1,Col_off) repmat('\t%*s',1,N_col)];
str_format = str_format(3:end);
C_text = textscan(fid,str_format,'HeaderLines',Row_off);
GeneSymbol_EXP = C_text{1,1};
if Row_off >= 2
    EntrezID_EXP = C_text{1,2};
else
    EntrezID_EXP = 0;
end

fclose(fid);
end