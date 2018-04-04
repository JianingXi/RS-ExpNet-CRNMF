CancerDataList = {'brca_tcga_pub';'coadread_tcga_pub';'gbm_tcga_pub'};


disp('Prior information: network');
network_dir = './network/Adj_mat.mat';
load(network_dir);
len_gene = length(GeneSymbol_net);

D_mat_half_inv = sparse(diag(sum(Adj_mat).^(-0.5)));
D_V = speye(len_gene);
W_V = D_mat_half_inv*Adj_mat*D_mat_half_inv;
L_V = D_V-W_V;

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
    W_U = exp(-((1-PearsonCor).^2)/(2*sigma_bandwidth^2));
    D_U = diag(sum(W_U,1));
    L_U = D_U - W_U;
    
    K_num = 4;
    lambda_LU = 1.0; lambda_RU = 1.0;
    lambda_LV = 1.0; lambda_RV = 1.0;
    
    % --- RS-CRNMF --- %
    disp('run RS-CRNMF ...');
    eps_t = 10^-5;
    
    % --- random initialization --- %
    U_init = max(rand(N_sample,K_num),eps_t);
    U_prev = U_init*diag(sum(U_init,1).^-1);
    V_init = max(rand(len_gene,K_num),eps_t);
    V_prev = V_init*diag(sum(U_init,1));
    delta_U = 1; delta_V = 1;
    cnt = 0;
    while delta_U > 10^-3 || delta_V > 10^-3 && cnt <= 500
        cnt = cnt + 1;
        V_numer = X_mut'*U_prev + lambda_LV*W_V*V_prev;
        V_denum = V_prev*(U_prev'*U_prev) + lambda_LV*D_V*V_prev + lambda_RV*ones(len_gene,len_gene)*V_prev;
        V_new = V_prev.*V_numer./(V_denum + eps_t);
        
        U_numer = X_mut*V_new + lambda_LU*W_U*U_prev;
        U_denum = U_prev*(V_new'*V_new) + lambda_LU*D_U*U_prev + lambda_RU*U_prev;
        U_new = U_prev.*U_numer./(U_denum + eps_t);
        
        NormFactor = sum(U_new,1);
        U_new = U_new*diag(NormFactor.^-1);
        V_new = V_new*diag(NormFactor);
        
        delta_U = norm(U_prev - U_new,'fro')^2/(norm(U_prev,'fro')^2);
        delta_V = norm(V_prev - V_new,'fro')^2/(norm(V_prev,'fro')^2);
        
        if ~mod(cnt,10)
            disp(['Iteration: ' num2str(cnt,'%d')]);
        end
        
        U_prev = U_new;
        V_prev = V_new;
    end
    disp(['Done!' char(10) 'Total iteration: ' num2str(cnt,'%d')]);
    
    [~, ind_gene] = sort(max(V_new,[],2),'descend');
    Candidates_list = GeneSymbol_net(ind_gene(1:200));

    save([output_save_dir '/result_' file_name_t '.mat'],'X_mut','U_new','V_new',...
        'Candidates_list');
end