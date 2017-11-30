CancerDataList = {'brca_tcga_pub';'coadread_tcga_pub';'gbm_tcga_pub'};


disp('Prior information: network');
network_dir = './network/Adj_mat.mat';
load(network_dir);
len_gene = length(GeneSymbol_net);

D_mat_half_inv = sparse(diag(sum(Adj_mat).^(-0.5)));
D_G = speye(len_gene);
W_G = D_mat_half_inv*Adj_mat*D_mat_half_inv;
L_G = D_G-W_G;

% Creating Directorys
output_save_dir = './output';
mkdir(output_save_dir);

for i_file = 1:length(CancerDataList)
    file_name_t = CancerDataList{i_file};
    input_mat_dir = ['./input_cancer_data/' file_name_t '.mat'];
    if ~exist(input_mat_dir,'file')
        continue;
    end
    
    disp([char(10) '-- -- File No.' num2str(i_file) ': ' file_name_t]);

    % read data
    disp('Prior information: mRNA expression pattern')
    temp_input_data = load(input_mat_dir);
    X_mut = (temp_input_data.X_SNV_Net~=0)+0;
    X_exp = temp_input_data.X_EXP_Net(:,sum(X_mut,1)>0);
    % sample_id = temp_input_data.SampleID_complete;
    N_sample = size(X_mut,1);

    clear temp_input_data
    
    PearsonCor = corr(X_exp');
    sigma_bandwidth = 1.0;
    W_S = exp(-((1-PearsonCor).^2)/(2*sigma_bandwidth^2));
    D_S = diag(sum(W_S,1));
    L_S = D_S - W_S;
    
    K_num = 4;
    lambda_S = 1.0;
    lambda_G = 1.0;
    
    % --- RS-CRNMF --- %
    disp('run RS-CRNMF ...');
    eps_t = 10^-5;
    
    % --- random initialization --- %
    S_init = max(rand(N_sample,K_num),eps_t);
    S_prev = S_init*diag(sum(S_init,1).^-1);
    G_init = max(rand(len_gene,K_num),eps_t);
    G_prev = G_init*diag(sum(S_init,1));
    delta_S = 1; delta_G = 1;
    cnt = 0;
    while delta_S > 10^-3 || delta_G > 10^-3 && cnt <= 500
        cnt = cnt + 1;
        G_numer = X_mut'*S_prev + lambda_G*W_G*G_prev;
        G_denum = G_prev*(S_prev'*S_prev) + lambda_G*D_G*G_prev + lambda_G*ones(len_gene,len_gene)*G_prev;
        G_new = G_prev.*G_numer./(G_denum + eps_t);

        S_numer = X_mut*G_new + lambda_S*W_S*S_prev;
        S_denum = S_prev*(G_new'*G_new) + lambda_S*D_S*S_prev + lambda_S*S_prev;
        S_new = S_prev.*S_numer./(S_denum + eps_t);

        NormFactor = sum(S_new,1);
        S_new = S_new*diag(NormFactor.^-1);
        G_new = G_new*diag(NormFactor);

        delta_S = norm(S_prev - S_new,'fro')^2/(norm(S_prev,'fro')^2);
        delta_G = norm(G_prev - G_new,'fro')^2/(norm(G_prev,'fro')^2);
        
        if ~mod(cnt,10)
            disp(['Iteration: ' num2str(cnt,'%d')]);
        end
        
        S_prev = S_new;
        G_prev = G_new;
    end
    disp(['Done!' char(10) 'Total iteration: ' num2str(cnt,'%d')]);
    save([output_save_dir '/result_' file_name_t '.mat'],'X_mut','S_new','G_new');
end