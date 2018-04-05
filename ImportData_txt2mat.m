% The demo file of RS-ExpNet-CRNMF
% Author: Jianing Xi, USTC

bin_path = './bin';
addpath(bin_path);

load('./network/Adj_mat.mat', 'GeneSymbol_net')
LenGeneNet = length(GeneSymbol_net);

root_dir = './raw_input_data';
file_list = dir(root_dir);

save_dir = './input_cancer_data';
mkdir(save_dir);

% 3:length(file_list)
for i_file = 1:length(file_list)
    cur_name = file_list(i_file).name;
    if strcmp(cur_name,'.') || strcmp(cur_name,'..')
        continue;
    end
    disp(cur_name);
    cur_file_dir = [root_dir '/' cur_name ];
    
    input_txt_file_SNV = [cur_file_dir '/data_mutations_extended.txt'];        
    [GeneSymbol_SNV,EntrezID_SNV,Variant_Classification,...
        Sample_ID_SNV] = P01_load_data_SNV(input_txt_file_SNV);

    input_txt_file_EXP = [cur_file_dir '/data_expression_median.txt'];
    [GeneSymbol_EXP,Sample_ID_EXP,X_EXP_raw,EntrezID_EXP] = ...
        P02_load_data_EXP(input_txt_file_EXP);
    
    [SampleID_complete,X_SNV_Net,X_EXP_Net]...
        = P03_Intersect_Sample(Sample_ID_SNV, GeneSymbol_SNV, ...
        Sample_ID_EXP,GeneSymbol_EXP,X_EXP_raw,GeneSymbol_net,...
        LenGeneNet);

    save([save_dir '/' cur_name '.mat'],'X_EXP_Net','X_SNV_Net','SampleID_complete');

end

rmpath(bin_path);