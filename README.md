## Users’ Manual of transTF-TWAS
## Overview
Transcriptome-wide association studies (TWAS) have been successful in identifying putative disease susceptibility genes by integrating cis-variants predicted gene expression with genome-wide association studies (GWAS) data, while trans-located variants for predicting gene expression remain largely unexplored. We introduce transTF-TWAS, by including transcription factor (TF)-linked trans-located variants to improve model building. Using data from the Genotype-Tissue Expression project and large disease GWAS datasets, we demonstrate that transTF-TWAS surpasses other existing TWAS approaches in both building gene prediction models and identifying disease-associated genes, from both simulations and real data analysis.  Our sTF-TWAS approach can significantly enhance the discovery of disease risk genes.

![My Image](Fig1A_B.png)

## Methods
### 1. Prepare input data: 
**1) Identifying cis-variants associate with TFs:**\
Step I - Identifying cis-variants associate with TFs: we firstly identified cis- variants located in STFCREs that potentially affect expression of a TF (e.g., FOXA1) by conducting cis-eQTL analysis and analyzing epigenetic data generated in breast-related cells. A set of the cis-variants regulating TF expression was determined based on the significant associations between TF gene expression and cis-variants (eQTL analysis), and/or variants with evidence of regulatory interactions with proximal promoters or distal enhancer-promoter regions (Fig. 1A). 

We downloaded approximately 3.6 million DNase I hypersensitive sites (DHSs) regions within human genome sequence (https://www.meuleman.org/research/dhsindex/). The enhancer regions were downloaded from EpiMap repository (http://compbio.mit.edu/epimap), which contains ~2M non-tissue specific enhancer regions. The CAGE peak regions were downloaded from FANTOM5 (https://fantom.gsc.riken.jp/5/data/), and we also included all regions within transcription start site (TSS) +/-2K for each gene as promoter regions. The eQTLs were downloaded from the GTEx portal (https://www.gtexportal.org/home/datasets/) and eQTLGen(https://www.eqtlgen.org/cis-eqtls.html). The Enhancer to gene link information across 833 cell-types was downloaded from EpiMap repository (https://personal.broadinstitute.org/cboix/epimap/links/links_corr_only/). We all used cell-type specific chromatin-chromatin interaction data from the 4D genomics(https://4dgenome.research.chop.edu/Tables/4DGenome_HomoSapiens_hg19.txt) and previous literature [1]. Finally, the TF-cis-regulatory-variants were identified based on the significant associations from eQTL results, and the regulatory evidence from the variants linked to the TF. Codes refers to `./code/WriteTFeQTL.py`.

Step II - TF-gene pair discovery:  we analyzed TF ChIP-seq data generated in breast cancer-related cells to characterize their genome-wide binding sites for susceptible TFs that have been identified in breast cancer from our prior work (Fig. 1B). We next characterized each gene potentially regulated by all possible susceptible TFs based on the evidence of the TF-DNA binding sites that are located in its flanking of transcription start sites (TSS, +/-20K; Fig. 1C). Cis-variants associated with TFs (subsequently termed TF-linked trans-variants) have the potential to regulate TF protein expression, which may consequently alter the gene expression of downstream targets. As a result, TF-linked trans-variants can influence the expression of genes located several megabases away or even on different chromosomes. Codes refers to `./code/WriteCSVFile.py` 

For each gene, we prepared a csv file that contains all TF-based trans-variants (TFx) for all TF-genes that regulate this target gene with the format as below:

TF,CHR,LOC,ID-1,ID-2,ID-3,ID-4,ID-5, … … \
GREB1_1,2,10494090,0,0,0,0,0,0,0,1, … … \
GREB1_2,2,10494743,2,1,1,1,1,2,0,2, … … \
 … … \
FOXA1_1,14,37708623,0,0,0,1,1,1,0, … … \
FOXA1_2,14,37709692,0,1,0,1,0,1,0, … … \
 … … 

**2)	Gene expression file:** \
Take the breast cancer as an example. The fully processed, filtered and normalized gene expression matrices in bed format ("Breast_Mammary.v8.normalized_expression.bed") for prostate tissue was downloaded from GTEx portal (https://gtexportal.org/home/datasets). We included 151 samples in our analysis and removed sex chromosomes, by which we generated a new file named "Breast_Mammary.v8.normalized_expression.no_sex.bed". The covariates used in eQTL analysis, including top five genotyping principal components (PCs), were obtained from GTEx_Analysis_v8_eQTL_covariates.tar.gz, which was downloaded from GTEx portal (https://gtexportal.org/home/datasets). Then, we further performed a probabilistic estimation of expression residuals (PEER) analysis to adjust for top five genotyping PCs, age, and other potential confounding factors (PEERs)[2] for downstream prediction model building. There is a description of how to download and use the PEER tool here: https://github.com/PMBio/peer/wiki/Tutorial. The command that we used is shown as below: 

`Rscript ./code/Peer_Script.R`

According to the GTEx protocol, if the number of samples is between 150 and 250, 30 PEER factors should be used. For our study, the number of samples is 151, so we used 30 PEER factors. This command will generate a residual file named “GeneExpression_Breast_AfterRM_Residuals.csv”, and from this residual file, we generated the final gene expression data file named “Breast_Mammary.v8.normalized_expression.no_sex.rm_covariates.bed” as the input for our downstream predictive model. 

**3)	genotype file:**  
The whole genome sequencing file, GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_866Indiv.vcf, was downloaded from dbGaP (https://www.ncbi.nlm.nih.gov/projects/gap/cgi-bin/study.cgi?study_id=phs000424.v8.p2). The genotype dataset is quality controlled using the tool PLINK [3] (https://zzz.bwh.harvard.edu/plink/ ). Multiple QC steps were applied by excluding variants with missingness rate > 0.1, minor allele frequency > 0.01, high deviations from Hardy-Weinberg equilibrium at p<10-6, and removing samples with missingness rate > 0.1.

Take the breast cancer as an example. We prioritized putative regulatory variants based on their associations with breast cancer risk. For variants that bind to only one TF, we used the single TF beta value, and for other variants that bind to more than one TFs, we considered the largest beta values of the paired TFs. Once we obtained the beta values for all TF-occupied elements, we ranked those variants based on the beta values from largest to smallest, which illustrated with more important to less for breast cancer risk. As illustrated in our previous work [4], we only included top 50K TF-occupied regulatory variants. 

The tissue specific input genotype file ("genotype_file") with the format as below:

 CHR,LOC,ID-1,ID-2,ID-3,ID-4,ID-5, … … \
 1,629906,0,0,0,0,0,0,0,0,1,0, … … \
 1,630026,1,0,0,2,1,0,0,0,2,1, … … \
 1,778639,0,0,0,0,0,0,0,0,0,0, … … \
 … …

The input genotype file for calculating covariance ("covariance_genotype_file") should be non-tissue specific, with the same format as the tissue specific input genotype file.


**3)	SNP annotation file:** \
The input snp annotation file ("snp_annot_file"), contains only the top 50K regulatory variants with the format as below: 

SNP,varID,chr,pos,ref,effect \
rs1578391,chr1_629906_C_T_b38,1,629906,C,T \
rs6594029,chr1_630026_C_T_b38,1,630026,C,T \
rs114983708,chr1_778639_A_G_b38,1,778639,A,G \
… …


**4)	Gene annotation file:** \
The input gene annot file ("gene_annot_file") is downloaded from GENCODE: https://www.gencodegenes.org/human/release_26.html, in the GTF format and build in GRCh38.


**5)	GWAS file:** \
The input GWAS file ("gwas_file") contains "chr" and "position" columns, we just need to make sure that all the SNPs being trained in the Elastic Net model can be found in the GWAS dataset.


### 2. Gene expression prediction model building based on trans-located variants:
We analyzed TF ChIP-seq data generated in target cancer-related cells to characterize their genome-wide binding sites for susceptible TFs using data from the Cistrome database (http://cistrome.org/).  We next characterized each gene potentially regulated by all possible susceptible TFs based on the evidence of their TF-DNA binding sites that are located in its flanking 20bk of TSS (i.e., number of G TFs; Fig. 1C).

Step III - Model training and disease-trait association analysis:  For each TF, we assessed the performance of a prediction model that utilized its TF linked trans-variants to predict expression of each target gene using Group Lasso method (i.e., number of G TFs; Fig. 1C). The Group Lasso’s property of encouraging between-group sparsity and within-group retainment aligns to our intention of selecting the actual functioning TFs and then retaining their cis-variants. The groups survive the regularization are corresponding to those of TF-linked trans-variants that may affect the expression of the gene. The final set of TF-linked trans-variants was identified for downstream gene expression model building by combining the groups from the significant models using Elastic Net (Fig. 1C). Under our sTF-TWAS framework, we next included the TF-linked trans-variants, together with the prioritized cis-variants to build gene expression prediction models. Here, we only focused on the set with 50K cis-variants (Fig. 1B), as the identified genes were highly overlapped among analyses with different number of variants (i.e., 50K vs 500K variants) in our prior work [4]. We conducted TWAS analyses by applying the gene expression prediction models, respectively, to GWAS summary statistics for breast, prostate, and lung cancers and other diseases to search for their susceptibility genes.

We processed one chromosome at a time by executing the below code, take chromsome 1 as an example:\
`Rscript ./code/transTF_TWAS.R  1`

### References: 
1. Rhie, S.K., Perez, A.A., Lay, F.D. et al. A high-resolution 3D epigenomic map reveals insights into the creation of the prostate cancer transcriptome. Nat Commun 10, 4154 (2019).
2. Stegle, O., Parts, L., Piipari, M., Winn, J., and Durbin, R. (2012). Using probabilistic estimation of expression residuals (PEER) to obtain increased power and interpretability of gene expression analyses. Nat Protoc 7, 500-507.
3. Purcell, S., Neale, B., Todd-Brown, K., Thomas, L., Ferreira, M.A., Bender, D., Maller, J., Sklar, P., de Bakker, P.I., Daly, M.J., et al. (2007). PLINK: a tool set for whole-genome association and population-based linkage analyses. Am J Hum Genet 81, 559-575.
4. He, J., Wen, W., Beeghly, A., Chen, Z., Cao, C., Shu, X.O., Zheng, W., Long, Q., and Guo, X. (2022). Integrating transcription factor occupancy with transcriptome-wide association analysis identifies susceptibility genes in human cancers. Nat Commun 13, 7118. 10.1038/s41467-022-34888-0.

### Contacts
  Jingni He: jingni.he1@ucalgary.ca<br>
  Quan Long: quan.long@ucalgary.ca<br>
  Xingyi Guo: xingyi.guo@vumc.org<br>


