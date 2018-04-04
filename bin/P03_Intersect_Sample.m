function [SampleID_complete,X_SNV_Net,X_EXP_Net]...
    = P03_Intersect_Sample(Sample_ID_SNV, GeneSymbol_SNV, ...
    Sample_ID_EXP,GeneSymbol_EXP,X_EXP_raw,GeneSymbol_net,LenGeneNet)

Sample_ID_SNV = cellfun(@(x) x(1:12),Sample_ID_SNV,'UniformOutput',0);
Sample_ID_EXP = cellfun(@(x) x(1:12),Sample_ID_EXP,'UniformOutput',0);

SampleID_complete = intersect(Sample_ID_EXP,Sample_ID_SNV);
LenSampleCom = length(SampleID_complete);


[Sample_ID_SNV_uni,~,Ind_Sample_uni] = unique(Sample_ID_SNV);
[GeneSymbol_SNV_uni,~,Ind_Gene_uni] = unique(GeneSymbol_SNV);

X_SNV_raw = sparse(Ind_Sample_uni,Ind_Gene_uni,1);
% Sample_ID_SNV_uni * GeneSymbol_SNV_uni

[~,Ind_comp,Ind_SNV] = intersect(SampleID_complete,Sample_ID_SNV_uni);
if norm(Ind_comp - (1:length(Ind_comp))',2) ~= 0
    error('SNV data sample ID mismatch.');
end
X_SNV_samp = X_SNV_raw(Ind_SNV,:);
% SampleID_complete * GeneSymbol_SNV_uni

X_SNV_Net = sparse(LenSampleCom,LenGeneNet);
[~,Ind_SNV,Ind_net] = intersect(GeneSymbol_SNV_uni,GeneSymbol_net);
X_SNV_Net(:,Ind_net) = X_SNV_samp(:,Ind_SNV);
% SampleID_complete * GeneSymbol_net



[~,Ind_comp,Ind_EXP] = intersect(SampleID_complete,Sample_ID_EXP);
if norm(Ind_comp - (1:length(Ind_comp))',2) ~= 0
    error('CNA data sample ID mismatch.');
end
X_EXP_samp = X_EXP_raw(Ind_EXP,:);
% SampleID_complete * GeneSymbol_EXP

X_EXP_Net = zeros(LenSampleCom,LenGeneNet);
[~,Ind_EXP,Ind_net] = intersect(GeneSymbol_EXP,GeneSymbol_net);
X_EXP_Net(:,Ind_net) = X_EXP_samp(:,Ind_EXP);
% SampleID_complete * GeneSymbol_net

end