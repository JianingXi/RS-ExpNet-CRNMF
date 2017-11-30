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

References
------------------------
[1] Cancer Genome Atlas Research Network and others: Comprehensive genomic characterization defines human glioblastoma genes and core pathways. Nature 455(7216), 1061 (2008)

[2] Cancer Genome Atlas Network and others: Comprehensive molecular characterization of human colon and rectal cancer. Nature 487(7407), 330-337 (2012)

[3] Cancer Genome Atlas Network and others: Comprehensive molecular portraits of human breast tumours. Nature 490(7418), 61-70 (2012)

[4] Gao, Jianjiong and Aksoy, Bu:lent Arman and Dogrusoz, Ugur and Dresdner, Gideon and Gross, Benjamin and Sumer, S Onur and Sun, Yichao and Jacobsen, Anders and Sinha, Rileen and Larsson, Erik and others: Integrative analysis of complex cancer genomics and clinical profiles using the cBioPortal. Science signaling 6(269), 1 (2013)

[5] Razick, S., Magklaras, G., Donaldson, I.M.: iRefIndex: a consolidated protein interaction database with provenance. BMC bioinformatics 9(1), 1 (2008)