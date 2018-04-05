# RS-ExpNet-CRNMF
Discovering mutated driver genes through a robust and sparse co-regularized matrix factorization framework with prior information from mRNA expression patterns and interaction network 
![image](https://github.com/JianingXi/RS-ExpNet-CRNMF/blob/master/bin/splash.jpg)

Developer: Jianing Xi <xjn@mail.ustc.edu.cn> from Health Informatics Lab, School of Information Science and Technology, University of Science and Technology of China

## Instructions to RS-ExpNet-CRNMF (version 1.0.0)

Requirement
------------------------
* 4GB memory
* MATLAB R2015a or later

Input data
------------------------
1. Somatic mutations of patients across multiple cancer types
The files `./input_cancer_data/gbm_tcga_pub.mat` `./input_cancer_data/coadread_tcga_pub.mat` and `./input_cancer_data/brca_tcga_pub.mat` contain the TCGA somatic mutation data matrix of three types of cancers, glioblastoma multiforme (GBM) [1], colon and rectal cancer (COADREAD) [2] and breast cancer (BRCA) [3], which is downloaded from [cBioPortal](http://www.cbioportal.org/data_sets.jsp) [4].

2. Prior information: mRNA expression data
The TCGA mRNA expression data of the samples of three types of cancers above are also provided in the files `./input_cancer_data/gbm_tcga_pub.mat` `./input_cancer_data/coadread_tcga_pub.mat` and `./input_cancer_data/brca_tcga_pub.mat`, downloaded from [cBioPortal](http://www.cbioportal.org/data_sets.jsp) [4].

3. Prior information: interaction network
The file `./network/Adj_mat.mat` contains the gene nodes and the edges of interaction network [iRefIndex 9](http://irefindex.org) [5].

run RS-ExpNet-CRNMF
------------------------
To apply RS-ExpNet-CRNMF, please run the Matlab script file `./demo.m` and the results will be automatically saved in file `./output/result_[cancer_file_name].mat` after the program is finished.

Output data
------------------------
In file `./output/result_[cancer_file_name].mat`, there are four output variables:

1. U_new
The variables U_new is the output sample representation matrix,
![image](https://github.com/JianingXi/RS-ExpNet-CRNMF/blob/master/bin/U_mat.PNG)
where K is the predefined dimension number of the latent representations. For k from 1 to K, the k-th vector u_{\*,k} indicates the assignment weights of the cancer cell sample to the k-th latent dimension. The i-th u_{i,\*} indicates the low-dimensional representations of the i-th cancer cell sample.

2. V_new
The matrix V_new is the gene representation matrix of the investigated genes,
![image](https://github.com/JianingXi/RS-ExpNet-CRNMF/blob/master/bin/V_mat.PNG)
where the k-th vector v_{\*,k} representing the weights of the tested genes in the k-th latent dimension. Each v_{j,\*} denotes the representations of the tested genes in the latent dimension. 

3. Candidates_list
This variable is a list containing driver gene candidates selected from the top 200 genes ranked by scores in gene representation matrix V_new.

Applying RS-ExpNet-CRNMF to other dataset
------------------------
To facilitate the applications to users' own dataset, we provide data a script file for data format transformation `./ImportData_txt2mat.m` to transform a txt-format dataset to mat-format data used in this package.
* Step 1: the users can create a new directory, and the directory name can be the name of the users' dataset. The txt-format files of somatic mutation and mRNA gene expression of the dataset should be included in this directory. An example of the application another dataset of kidney renal clear cell carcinoma (KIRC) is provided as two txt-format files `./raw_input_data/kirc_tcga_pub/data_mutations_extended.txt` (somatic mutation) and `./raw_input_data/kirc_tcga_pub/data_expression_median.txt` (mRNA gene expression) (the two files are downloaded from [cBioPortal](http://www.cbioportal.org/data_sets.jsp) [4]).
* Step 2: run the Matlab script file `./ImportData_txt2mat.m`, and the user's own dataset in directory `./raw_input_data/` will be automatically scanned and transformed to mat-format file.
* Step 3: the transformed mat-format file will be automatically saved in directory `./input_cancer_data`, which is the directory of mat-format input file for RS-ExpNet-CRNMF. Taking KIRC dataset as an example, if the script file `./ImportData_txt2mat.m` has been run, a mat-format file `./input_cancer_data/kirc_tcga_pub.mat` will be automatically created in directory `./input_cancer_data/`.
* Step 4: please run the Matlab script file `./demo.m` and the mat-format file of the users' dataset will be automatically scanned and analyzed by RS-ExpNet-CRNMF. The results as output file will be save after the program is finished. For KIRC example, when the script of RS-ExpNet-CRNMF is finished, the output file of KIRC result will be automatically saved in `./output/result_kirc_tcga_pub.mat`.


References
------------------------
[1] Cancer Genome Atlas Research Network and others: Comprehensive genomic characterization defines human glioblastoma genes and core pathways. Nature 455(7216), 1061 (2008)

[2] Cancer Genome Atlas Network and others: Comprehensive molecular characterization of human colon and rectal cancer. Nature 487(7407), 330-337 (2012)

[3] Cancer Genome Atlas Network and others: Comprehensive molecular portraits of human breast tumours. Nature 490(7418), 61-70 (2012)

[4] Gao, Jianjiong and Aksoy, Bu:lent Arman and Dogrusoz, Ugur and Dresdner, Gideon and Gross, Benjamin and Sumer, S Onur and Sun, Yichao and Jacobsen, Anders and Sinha, Rileen and Larsson, Erik and others: Integrative analysis of complex cancer genomics and clinical profiles using the cBioPortal. Science signaling 6(269), 1 (2013)

[5] Razick, S., Magklaras, G., Donaldson, I.M.: iRefIndex: a consolidated protein interaction database with provenance. BMC bioinformatics 9(1), 1 (2008)