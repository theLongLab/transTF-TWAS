## Users’ Manual of transTF-TWAS

## Overview

## Methods
### 1. Prepare input data: 
**1)	Functional regions include:
1. DNase I hypersensitive sites (DHSs) region is downloaded from https://www.meuleman.org/research/dhsindex/, and we only included “Cancer/epithelial” specific DHSs regions.
2. Enhancer region is downloaded from Epimap http://compbio.mit.edu/epimap/.
3. Cap Analysis Gene Expression (CAEG) peak is downloaded from fantom https://fantom.gsc.riken.jp/5/data/hg19.cage_peak_phase1and2combined_ann.txt.gz, indicating 4. promoter regions.
4. Promoter regions: we also include transcription start site (TSS) +/- 2k as promoter regions

**2) 
Hi-C group: If promoter region of the target gene (TSS+/-2K) overlapped with the region A of the interacted regions in 3D space, then all prioritized functional variants that located within the region B are all included as trans-variants for the target gene.

Enhancer-gene-link group: For each target gene, we included those prioritized functional variants that sit within the enhancer regions with the gene links

TF group: We first identified TFs that regulated target genes by including all TFs where their DNA binding variants sit within TSS+/-20K of the target gene. Then, for all the TFx that linked to the TF-genes are all considered as potential trans-variants for the target gene. In order to identify TFx for each TF-gene:
1) By checking the significant variant-gene pairs downloaded from GTEx portal (https://www.gtexportal.org/home/datasets/) , we would include all variant to TF-gene pairs with nominal p-value <0.05, and ensure that the variants sit within +/- 1M of the TF-gene gene body; 
2)We would include variants that have been utilized for generating imputed GREx for the TF-gene;
3)By looking for variants that link to the TF-gene from EpiMap enhancer-gene linking, where we would include variants within the enhancer regions that links to the target TF-gene as TFx;
4)Lastly, we would include variants that are within the promoter region of the target TF-gene (TSS +/- 20K).

**3）


