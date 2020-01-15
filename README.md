## ICIBM_2019_SS_TSJ
# Pseudogene-gene functional networks are prognostic of patient survival in breast cancer.
Sasha Smerekanych<sup>1,2¶</sup>, Travis S Johnson<sup>2,3¶</sup>, Kun Huang<sup>3,4</sup>, Yan Zhang<sup>2,5*</sup>

<sup>1</sup> Kenyon College, Gambier, OH 43022, United States
<sup>2</sup> Department of Biomedical Informatics, College of Medicine, The Ohio State University, Columbus, OH 43210, United States
<sup>3</sup> Department of Medicine, School of Medicine, Indiana University, Indianapolis, IN 46202, USA 
<sup>4</sup> Regenstrief Institute, Indiana University, Indianapolis, IN 46262, USA
<sup>5</sup> The Ohio State University Comprehensive Cancer Center (OSUCCC – James), Columbus, OH 43210, United States
<sup>¶</sup> These authors contributed equally to this work.
<sup>*</sup> Correspondence: yan.zhang@osumc.edu

This github contains the relevant R scripts and data for the paper above.

Note that the GSE96058 (equivalent to GSE81538) data must be downloaded from GEO since the file size is too large to be hosted on github. Similarly the BRCA.rds file containing the pseudogene expression must be downloaded from the following link (https://www.dropbox.com/s/vi9mu46tfxr94df/BRCA.rds?dl=0) since it is too large to be hosted on github.

File descriptions:

ICIBM_Paper68_Supplement_code.rtfd.zip contains the R code for the analysis. The majority of the manuscript excluding the Swedish cohort (GSE81538) are included.

MANIFEST.txt contains the information of TCGA RNA expression data used for the genes in the analysis.

annot.Rdata is the gene and pseudogene annotation used during the analysis. This file is used int the R code above.

brca_gse81538.R is the analysis performed on the Swedish cohort of patients (GSE81538). This file was used to generate Supplementary_table_1.txt.

edgesblast_formatted.csv contains the pseudogene-gene functional networks as a list of edges. This file was used to generate the gene-pseudogene interaction terms in the R notebook file.
