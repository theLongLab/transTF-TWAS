## Users’ Manual of transTF-TWAS
## Overview
Transcriptome-wide association studies (TWAS) have been successful in identifying many disease susceptibility genes by integrating genetically predicted gene expression with genome-wide association study (GWAS) data. However, current TWAS models only consider cis-located variants. In this study, we introduce transTF-TWAS, which includes additional transcription factor (TF)-based trans-located variants for prediction model building. Our results show that transTF-TWAS outperforms other methods by significantly improving prediction models and identifying disease genes. Using data from the Genotype-Tissue Expression project, we predicted alternative splicing and gene expression and applied these models to large GWAS datasets for breast, prostate, and lung cancers. Our analysis revealed a total of 940 putative cancer susceptibility genes, including 492 that were previously unreported in GWAS loci and 143 in loci unreported by GWAS, at Bonferroni-corrected P < 0.05. These findings provide additional insight into the genetic susceptibility of human cancers and demonstrate that our approach can strengthen disease gene discovery.

![My Image](Fig1A_B.PNG)

## Methods
### 1. Prepare input data: 
**1) Description of the prioritized trans-variants used in this study:**\
We first prioritized functional variants by including TF-occupied variants that fell into functional regions. Functional regions include DNase I hypersensitive sites (DHSs), enhancer regions and promoter regions. DHSs was downloaded from https://www.meuleman.org/research/dhsindex/, which in total contain ~3.6M DHSs, and we only kept tissue specific DHS regions, i.e. we kept “Cancer/epithelial” related DHS regions for breast, prostate and lung cancers. Enhancers regions was downloaded from EpiMap http://compbio.mit.edu/epimap/ , which contains 2M non-tissue specific enhancers regions. For the promoter regions, we first utilized CAGE peak regions downloaded from Fantom https://fantom.gsc.riken.jp/5/data/, and also included transcription start site (TSS) +/- 2k for all genes as promoter regions.
For identifying trans-variants: 
Transcriptome factors (TFs), i.e. FOXA1, are transcripted and translated from genes that are named as TF-genes, (Fig. 1A). TF-genes are linked or regulated by variants called TFx. We first identified TFs that regulated target genes by including all TFs where their DNA binding variants sit within TSS+/-20K of the target gene. Then, for all the TFx that linked to the TF-genes are all considered as potential trans-variants for the target gene. To identify TFx for each TF-gene, we proposed four ways:
1)	By checking the significant variant-gene pairs downloaded from GTEx portal (https://www.gtexportal.org/home/datasets/), including both tissue specific eQTLs and eQTLs from whole blood. We also included all cis-eQTL downloaded from eQTLGen (https://www.eqtlgen.org/cis-eqtls.html). We would include all variant to TF-gene pairs with nominal p-value <0.05, and ensure that the variants sit within +/- 1M of the TF-gene gene body. 
2)	By looking for variants that link to the TF-gene from EpiMap enhancer-gene linking, where we would include variants within the enhancer regions that links to the target TF-gene as TFx. We used linked EpiMap enhancers based on expression-enhancer activity correlation across 833 cell-types that was downloaded from https://personal.broadinstitute.org/cboix/epimap/links/links_corr_only/ ; 
3)	Lastly, we would include variants that are within the promoter region of the target TF-gene (TSS +/- 2K). We also considered different choices for regions that harbor cis-regulatory variants, which includes the ones defined by linear windows surrounding gene start and end sites as well as by 3D genomics informed regions (i.e., Hi-C region). Hi-C measures the frequency at which two DNA fragments physically associate in 3D space, linking chromosomal structure directly to the genomic sequence.

For each gene, we prepared a csv file that contains all TF-based trans-variants (TFx) for all TF-genes that regulate this target gene with the format as below:

TF,CHR,LOC,GTEX-1117F,GTEX-1122O,GTEX-11EM3,GTEX-11EMC,GTEX-11GSP,GTEX-11I78,GTEX-11P81, … … \
GREB1_1,2,10494090,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, … … \
GREB1_2,2,10494743,2,1,1,1,1,2,0,2,0,0,2,0,1,1,1,1,1,1,0,1,0,0,1,1,2,0,1,1,1,1,1,2,0,0,0,0,1,1,0,0,1,0,1,1,1, … … \
GREB1_3,2,10494930,2,1,2,2,1,2,1,2,0,2,2,1,1,2,2,1,2,2,1,2,2,2,2,2,2,2,2,2,2,1,2,2,1,1,2,2,1,2,2,2,1,1,2,2,2, … … \
 … … \
FOXA1_1,14,37708623,0,0,0,1,1,1,0,0,0,0,0,0,0,1,0,1,0,0,1,0,0,0,0,0,1,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, … … \
FOXA1_2,14,37709692,0,1,0,1,0,1,0,1,1,2,2,0,2,1,2,1,0,0,2,0,0,0,1,2,2,1,0,2,1,1,0,2,0,1,1,1,0,1,0,1,1,1,0,2, … … \
FOXA1_3,14,37709729,2,1,2,1,1,1,1,1,1,0,1,2,0,0,1,1,2,1,0,2,2,1,1,0,1,1,2,0,1,2,1,0,2,1,1,1,1,1,1,1,0,1,2,2, … … \
 … … 

**2)	Gene expression file:** \
Take the breast cancer as an example. The fully processed, filtered and normalized gene expression matrices in bed format ("Breast_Mammary.v8.normalized_expression.bed") for prostate tissue was downloaded from GTEx portal (https://gtexportal.org/home/datasets). We included 151 samples in our analysis and removed sex chromosomes, by which we generated a new file named "Breast_Mammary.v8.normalized_expression.no_sex.bed". The covariates used in eQTL analysis, including top five genotyping principal components (PCs), were obtained from GTEx_Analysis_v8_eQTL_covariates.tar.gz, which was downloaded from GTEx portal (https://gtexportal.org/home/datasets). Then, we further performed a probabilistic estimation of expression residuals (PEER) analysis to adjust for top five genotyping PCs, age, and other potential confounding factors (PEERs)[1] for downstream prediction model building. There is a description of how to download and use the PEER tool here: https://github.com/PMBio/peer/wiki/Tutorial. The command that we used is shown as below: 

`Rscript ./code/Peer_Script.R`

According to the GTEx protocol, if the number of samples is between 150 and 250, 30 PEER factors should be used. For our study, the number of samples is 151, so we used 30 PEER factors. This command will generate a residual file named “GeneExpression_Breast_AfterRM_Residuals.csv”, and from this residual file, we generated the final gene expression data file named “Breast_Mammary.v8.normalized_expression.no_sex.rm_covariates.bed” as the input for our downstream predictive model. 

**3)	genotype file:**  
The whole genome sequencing file, GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_866Indiv.vcf, was downloaded from dbGaP (https://www.ncbi.nlm.nih.gov/projects/gap/cgi-bin/study.cgi?study_id=phs000424.v8.p2). The genotype dataset is quality controlled using the tool PLINK [2] (https://zzz.bwh.harvard.edu/plink/ ). Multiple QC steps were applied by excluding variants with missingness rate > 0.1, minor allele frequency < 0.01, high deviations from Hardy-Weinberg equilibrium at p<10-6, and removing samples with missingness rate > 0.1.

Take the breast cancer as an example. We prioritized putative regulatory variants based on their associations with breast cancer risk. For variants that bind to only one TF, we used the single TF beta value, and for other variants that bind to more than one TFs, we considered the largest beta values of the paired TFs. Once we obtained the beta values for all TF-occupied elements, we ranked those variants based on the beta values from largest to smallest, which illustrated with more important to less for breast cancer risk. As illustrated in our previous work [3], we only included top 50K TF-occupied regulatory variants. 

The tissue specific input genotype file ("genotype_file") with the format as below:

 CHR,LOC,GTEX-1117F,GTEX-1122O,GTEX-11EM3,GTEX-11EMC,GTEX-11GSP,GTEX-11I78,GTEX-11P81, … … \
 1,933303,0,0,0,0,0,0,0,0,1,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, … … \
 1,933411,1,2,2,2,2,2,2,2,2,2,2,2,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2, … … \
 1,933653,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,1,0,0,0, … … \
 … …

The input genotype file for calculating covariance ("covariance_genotype_file") should be non-tissue specific, with the same format as the tissue specific input genotype file.


**3)	SNP annotation file:** \
The input snp annotation file ("snp_annot_file"), contains only the top 50K regulatory variants with the format as below: 

SNP,varID,chr,pos,ref,effect \
rs1578391,chr1_629906_C_T_b38,1,629906,C,T \
rs6594029,chr1_630026_C_T_b38,1,630026,C,T \
rs114983708,chr1_778639_A_G_b38,1,778639,A,G \
rs71507461,chr1_827209_G_C_b38,1,827209,G,C \
rs71507462,chr1_827212_C_G_b38,1,827212,C,G \
… …


**4)	Gene annotation file:** \
The input gene annot file ("gene_annot_file") is downloaded from GENCODE: https://www.gencodegenes.org/human/release_26.html, in the GTF format and build in GRCh38.


**5)	GWAS file:** \
The input GWAS file ("gwas_file") contains "chr" and "position" columns, we just need to make sure that all the SNPs being trained in the Elastic Net model can be found in the GWAS dataset.


### 2. Training the predictive model:
We first trained a group lasso to select essential trans- groups. In total, there are G groups of TF-based trans-variants from each TF-gene. Once essential trans- groups are being selected, we built gene-expression prediction models for the sets of interest using processed genetic variants and normalized gene expression data generated in different tissue samples from GTEx v8. However, instead of using all flanking genetic variants (flanking ±1Mb region, MAF > 1%) to train the elastic-net model, like the regular TWAS does, we included the TF-occupied variants (50K) within flanking region of each gene plus all trans-variants from selected trans-groups. 

We processed one chromosome at a time by executing this code, take chromsome 1 as an example:\
`Rscript ./code/transTF_TWAS.R  1`

### References: 
1. Stegle, O., Parts, L., Piipari, M., Winn, J., and Durbin, R. (2012). Using probabilistic estimation of expression residuals (PEER) to obtain increased power and interpretability of gene expression analyses. Nat Protoc 7, 500-507.
2. Purcell, S., Neale, B., Todd-Brown, K., Thomas, L., Ferreira, M.A., Bender, D., Maller, J., Sklar, P., de Bakker, P.I., Daly, M.J., et al. (2007). PLINK: a tool set for whole-genome association and population-based linkage analyses. Am J Hum Genet 81, 559-575.
3. He, J., Wen, W., Beeghly, A., Chen, Z., Cao, C., Shu, X.O., Zheng, W., Long, Q., and Guo, X. (2022). Integrating transcription factor occupancy with transcriptome-wide association analysis identifies susceptibility genes in human cancers. Nat Commun 13, 7118. 10.1038/s41467-022-34888-0.

### Contacts
  Jingni He: jingni.he1@ucalgary.ca<br>
  Quan Long: quan.long@ucalgary.ca<br>
  Xingyi Guo: xingyi.guo@vumc.org<br>


