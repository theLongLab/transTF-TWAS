import sys
import os
import csv

TF_file=open("Twentytwo_TFs.txt",'r') ## TF      beta*   se      p, regarding to suscetibility TF towards cancer
lines=TF_file.readlines()
tf_gene_list=[]
for line in lines[1:]:
    tf_gene_list.append(line.split("\t")[0].strip())
TF_file.close()

gencode_file=open("/work/long_lab/jingni/project/trans_TWAS/Data_Dir/ori_data/TSS_Promoter/gencode.v26.GRCh38.genes.csv",'r')
gene_id_list=[]
geneid_geneid2_dict={}
geneid_name_dict={}
geneid_start_dict={}
geneid_end_dict={}
geneid_chr_dict={}
frank_region = 1000000
for line in gencode_file:
    gene_name=line.split(",")[8].split("\"")[1].strip()
    gene_id=line.split(",")[5].split("\"")[1].strip()
    if gene_name in tf_gene_list:
        curr_chr=line.split(",")[0].strip().split("chr")[1].strip()
        loc_start=int(line.split(",")[2].strip())
        loc_end=int(line.split(",")[3].strip())
        if loc_start < frank_region:
            loc_start = 1;
        else:
            loc_start = loc_start - frank_region
        loc_end = loc_end + frank_region
        geneid_start_dict[gene_id]=loc_start
        geneid_end_dict[gene_id]=loc_end
        gene_id_list.append(gene_id)
        geneid_geneid2_dict[gene_id.split(".")[0].strip()]=gene_id
        geneid_chr_dict[gene_id]=curr_chr
        geneid_name_dict[gene_id]=gene_name
gencode_file.close()
print(len(gene_id_list))

geneid_variant_dict={}
## Add trans by identified TF-gene associated eQTL
##a set of the TF-cis-regulatory-variants was determined based on the eQTL analysis in both target tissues and whole blood samples using data from GTEx portal and eQTLGen(at a nominal p-value < 0.05). 
with open("/PATH_TO/GTEx_Analysis_v8_eQTL/Breast_Mammary_Tissue.v8.signif_variant_gene_pairs.txt",'r') as eqtl_file:
    readeqtl = csv.reader(eqtl_file, delimiter='\t')
    header = next(readeqtl)
    for row in readeqtl:
        if (row[1].strip() in gene_id_list) and (float(row[6].strip())<0.05):
            if row[1].strip() not in geneid_variant_dict.keys():
                geneid_variant_dict[row[1].strip()]=[]
            if (int(row[0].strip().split("_")[1].strip()) >= geneid_start_dict[row[1].strip()]) and (int(row[0].strip().split("_")[1].strip()) <= geneid_end_dict[row[1].strip()]):
                geneid_variant_dict[row[1].strip()].append(row[0].strip())
eqtl_file.close()
print("eQTL_1")
print(len(geneid_variant_dict))

with open("/PATH_TO/GTEx_Analysis_v8_eQTL/Whole_Blood.v8.signif_variant_gene_pairs.txt",'r') as eqtl_file:
    readeqtl = csv.reader(eqtl_file, delimiter='\t')
    header = next(readeqtl)
    for row in readeqtl:
        if (row[1].strip() in gene_id_list) and (float(row[6].strip())<0.05):
            if row[1].strip() not in geneid_variant_dict.keys():
                geneid_variant_dict[row[1].strip()]=[]
            if (int(row[0].strip().split("_")[1].strip()) >= geneid_start_dict[row[1].strip()]) and (int(row[0].strip().split("_")[1].strip()) <= geneid_end_dict[row[1].strip()]):
                geneid_variant_dict[row[1].strip()].append(row[0].strip())
eqtl_file.close()
print("eQTL_2")
print(len(geneid_variant_dict))

with open("/PATH_TO/2019-12-11-cis-eQTLsFDR-ProbeLevel-CohortInfoRemoved-BonferroniAdded_b38.txt",'r') as eqtl_file:
    readeqtl = csv.reader(eqtl_file, delimiter='\t')
    header = next(readeqtl)
    for row in readeqtl:
        if row[7].strip() in geneid_geneid2_dict.keys():
            curr_geneid=geneid_geneid2_dict[row[7].strip()]
            if (curr_geneid in gene_id_list) and (float(row[0].strip())<0.05):
                if curr_geneid not in geneid_variant_dict.keys():
                    geneid_variant_dict[curr_geneid]=[]
                if (int(row[3].strip()) >= geneid_start_dict[curr_geneid]) and (int(row[3].strip()) <= geneid_end_dict[curr_geneid]):
                    tmp_chrpos="chr"+row[2].strip()+"_"+row[3].strip()+"_"+row[4].strip()+"_"+row[5].strip()+"_b38"
                    geneid_variant_dict[curr_geneid].append(tmp_chrpos)
eqtl_file.close()
print("eQTL_3")
print(len(geneid_variant_dict))
#print(geneid_variant_dict)

#Check overlapped eQTL with functional varants
#Functional variants are a group of TF-occupied variants that are located within DHSs, enhancer regions and promoter regions. 
sel_variant_dict={}
sel_variant_anno_dict={}
for curr_geneid in geneid_name_dict.keys():
    curr_chr=geneid_chr_dict[curr_geneid]
    if curr_chr !="X":
        func_file=open("/PATH_TO/TF_Combine_All_Function_chr"+curr_chr+".txt",'r')
        for line in func_file:
            curr_b38="chr"+line.split("\t")[0].strip()+"_b38"
            if curr_geneid in geneid_variant_dict.keys():
                if curr_b38 in geneid_variant_dict[curr_geneid]:
                    if curr_geneid not in sel_variant_dict.keys():
                        sel_variant_dict[curr_geneid]=[]
                        sel_variant_anno_dict[curr_geneid]=[]
                    sel_variant_dict[curr_geneid].append(curr_b38)
                    sel_variant_anno_dict[curr_geneid].append("eQTL:"+curr_b38)
        func_file.close()
print("Check_eQTL_overlap_with_functional_variants")
print(len(sel_variant_dict.keys()))

##we examined if these variants are located in the promoter region of a TF (TSS +/- 2K) 
geneid_tss_dict={}
tss_file=open("/PATH_TO/TSS_Promoter/TSS_Region.txt",'r')
for line in tss_file:
    curr_geneid=line.split("\t")[0].strip()
    if curr_geneid in geneid_name_dict.keys():
        curr_tss=int(line.split("\t")[4].strip())
        curr_start=curr_tss-2000
        curr_end=curr_tss+2000
        geneid_tss_dict[curr_geneid]=str(curr_start)+"_"+str(curr_end)
tss_file.close()

for curr_geneid in geneid_name_dict.keys():
    curr_chr=geneid_chr_dict[curr_geneid]
    if (curr_chr !="X") and (curr_geneid in geneid_tss_dict.keys()):
        tss_start=int(geneid_tss_dict[curr_geneid].split("_")[0].strip())
        tss_end=int(geneid_tss_dict[curr_geneid].split("_")[1].strip())
        func_file=open("/PATH_TO/TF_Combine_All_Function_chr"+curr_chr+".txt",'r')
        lines=func_file.readlines()
        for line in lines[1:]:
            curr_b38="chr"+line.split("\t")[0].strip()+"_b38"
            curr_pos=int(line.split("\t")[3].strip())
            if (curr_pos >= tss_start) and (curr_pos <= tss_end):
                if curr_geneid not in sel_variant_dict.keys():
                    sel_variant_dict[curr_geneid]=[]
                    sel_variant_anno_dict[curr_geneid]=[]
                sel_variant_dict[curr_geneid].append(curr_b38)
                sel_variant_anno_dict[curr_geneid].append("Promoter:"+curr_b38)
        func_file.close()

## Add Trans using enhancer_link
## Enhance region with an evidence of the enhancer linking to the TF based on expression-enhancer activity correlation across 833 cell-types from the EpiMap repository 
geneid_enhancer_dict={}
for curr_geneid in geneid_name_dict.keys():
    ## Enhancer_Gene_Links_* are coder in b38
    target_chr=geneid_chr_dict[curr_geneid]
    tss_file="/PATH_TO/Enhancer_Gene_Links_"+curr_geneid.split(".")[0]+".txt"
    if os.path.isfile(tss_file):
        en_tss_file=open(tss_file)
        en_lines=en_tss_file.readlines()
        loc_start=geneid_start_dict[curr_geneid]
        loc_end=geneid_end_dict[curr_geneid]
        enhancer_list=[]
        for line in en_lines[1:]:
            tmp_start=int(line.split("\t")[2].strip())
            tmp_end=int(line.split("\t")[3].strip())
            if (tmp_start >= loc_start) and (tmp_end <= loc_end):
                enhancer_list.append(str(tmp_start)+"_"+str(tmp_end))
        #print("enh")
        #print(len(enhancer_list))
        if (len(enhancer_list)!=0):
            geneid_enhancer_dict[curr_geneid]=enhancer_list
#print(len(geneid_enhancer_dict))

for curr_geneid in geneid_name_dict.keys():
    if curr_geneid in geneid_enhancer_dict.keys():
        curr_chr=geneid_chr_dict[curr_geneid]
        if curr_chr !="X":
            tmp_enhancer_list=geneid_enhancer_dict[curr_geneid]
            for curr_region in tmp_enhancer_list:
                tss_start=int(curr_region.split("_")[0].strip())
                tss_end=int(curr_region.split("_")[1].strip())
                func_file=open("/PATH_TO/TF_Combine_All_Function_chr"+curr_chr+".txt",'r')
                lines=func_file.readlines()
                for line in lines[1:]:
                    curr_b38="chr"+line.split("\t")[0].strip()+"_b38"
                    curr_pos=int(line.split("\t")[3].strip())
                    if (curr_pos >= tss_start) and (curr_pos <= tss_end):
                        if curr_geneid not in sel_variant_dict.keys():
                            sel_variant_dict[curr_geneid]=[]
                            sel_variant_anno_dict[curr_geneid]=[]
                        sel_variant_dict[curr_geneid].append(curr_b38)
                        sel_variant_anno_dict[curr_geneid].append("Enhancer:"+curr_b38)
        func_file.close()
print("Check_enhancer_overlap_with_functional_variants")
print(len(sel_variant_dict.keys()))

## Add Trans using 4D information
HiC_celltype_list=["Breast","HCC1954","HMEC","MCF7"]
geneid_4d_dict={}
for curr_geneid in geneid_name_dict.keys():
    target_chr=geneid_chr_dict[curr_geneid]
    chr_start_end_dict={}
    if (target_chr!="X") and (curr_geneid in geneid_tss_dict.keys()):
        loc_start=geneid_start_dict[curr_geneid]
        loc_end=geneid_end_dict[curr_geneid]
        tss_start=int(geneid_tss_dict[curr_geneid].split("_")[0].strip())
        tss_end=int(geneid_tss_dict[curr_geneid].split("_")[1].strip())
        with open("/PATH_TO/4D/4DGenome_HomoSapiens_hg38.final.txt") as region1_4d_file:
            read4d = csv.reader(region1_4d_file, delimiter='\t')
            header = next(read4d)
            for row in read4d:
                region2_chr=row[3].strip()
                if (row[0].strip() == "chr"+target_chr) and (region2_chr == "chr"+target_chr) and (row[9].strip() in HiC_celltype_list):
                    tmp_start=int(row[1].strip())
                    tmp_end=int(row[2].strip())
                    region2_start=row[4].strip()
                    region2_end=row[5].strip()
                    if (int(region2_start) >= loc_start) and (int(region2_end) <= loc_end): # Ensure the region2 are within cis region of the TF-GENE
                        if tmp_start<=tss_start and tmp_end>=tss_start and tmp_end <=tss_end:
                            #print("sel1")
                            #print(tmp_start,tmp_end)
                            if region2_chr not in chr_start_end_dict.keys():
                                chr_start_end_dict[region2_chr]=[]
                            sel_region2=region2_start+"_"+region2_end
                            if sel_region2 not in chr_start_end_dict[region2_chr]:
                                chr_start_end_dict[region2_chr].append(sel_region2)
                        elif tmp_start<=tss_start and tmp_end>=tss_end:
                            #print("sel2")
                            #print(tmp_start,tmp_end)
                            if region2_chr not in chr_start_end_dict.keys():
                                chr_start_end_dict[region2_chr]=[]
                            sel_region2=region2_start+"_"+region2_end
                            if sel_region2 not in chr_start_end_dict[region2_chr]:
                                chr_start_end_dict[region2_chr].append(sel_region2)
                        elif tmp_start>=tss_start and tmp_end>=tss_end and tmp_start <=tss_end:
                            #print("sel3")
                            #print(tmp_start,tmp_end)
                            if region2_chr not in chr_start_end_dict.keys():
                                chr_start_end_dict[region2_chr]=[]
                            sel_region2=region2_start+"_"+region2_end
                            if sel_region2 not in chr_start_end_dict[region2_chr]:
                                chr_start_end_dict[region2_chr].append(sel_region2)
                        elif tmp_start>=tss_start and tmp_end<=tss_end:
                            #print("sel4")
                            #print(tmp_start,tmp_end)
                            if region2_chr not in chr_start_end_dict.keys():
                                chr_start_end_dict[region2_chr]=[]
                            sel_region2=region2_start+"_"+region2_end
                            if sel_region2 not in chr_start_end_dict[region2_chr]:
                                chr_start_end_dict[region2_chr].append(sel_region2)
        region1_4d_file.close()
        with open("/PATH_TO/4D/4DGenome_HomoSapiens_hg38.final.txt") as region2_4d_file:
            read4d = csv.reader(region2_4d_file, delimiter='\t')
            header = next(read4d)
            for row in read4d:
                region1_chr=row[0].strip()
                if (row[3].strip() == "chr"+target_chr) and (region1_chr == "chr"+target_chr) and (row[9].strip() in HiC_celltype_list):
                    tmp_start=int(row[4].strip())
                    tmp_end=int(row[5].strip())
                    region1_start=row[1].strip()
                    region1_end=row[2].strip()
                    if (int(region1_start) >= loc_start) and (int(region1_end) <= loc_end):
                        if tmp_start<=tss_start and tmp_end>=tss_start and tmp_end <=tss_end:
                            #print("sel1")
                            #print(tmp_start,tmp_end)
                            if region1_chr not in chr_start_end_dict.keys():
                                chr_start_end_dict[region1_chr]=[]
                            sel_region1=region1_start+"_"+region1_end
                            if sel_region1 not in chr_start_end_dict[region1_chr]:
                                chr_start_end_dict[region1_chr].append(sel_region1)
                        elif tmp_start<=tss_start and tmp_end>=tss_end:
                            #print("sel2")
                            #print(tmp_start,tmp_end)
                            if region1_chr not in chr_start_end_dict.keys():
                                chr_start_end_dict[region1_chr]=[]
                            sel_region1=region1_start+"_"+region1_end
                            if sel_region1 not in chr_start_end_dict[region1_chr]:
                                chr_start_end_dict[region1_chr].append(sel_region1)
                        elif tmp_start>=tss_start and tmp_end>=tss_end and tmp_start <=tss_end:
                            #print("sel3")
                            #print(tmp_start,tmp_end)
                            if region1_chr not in chr_start_end_dict.keys():
                                chr_start_end_dict[region1_chr]=[]
                            sel_region1=region1_start+"_"+region1_end
                            if sel_region1 not in chr_start_end_dict[region1_chr]:
                                chr_start_end_dict[region1_chr].append(sel_region1)
                        elif tmp_start>=tss_start and tmp_end<=tss_end:
                            #print("sel4")
                            #print(tmp_start,tmp_end)
                            if region1_chr not in chr_start_end_dict.keys():
                                chr_start_end_dict[region1_chr]=[]
                            sel_region1=region1_start+"_"+region1_end
                            if sel_region1 not in chr_start_end_dict[region1_chr]:
                                chr_start_end_dict[region1_chr].append(sel_region1)
        region2_4d_file.close()
        print(curr_geneid)
        print(chr_start_end_dict)
        if "chr"+target_chr in chr_start_end_dict.keys():
            geneid_4d_dict[curr_geneid]=chr_start_end_dict["chr"+target_chr]

for curr_geneid in geneid_name_dict.keys():
    if curr_geneid in geneid_4d_dict.keys():
        curr_chr=geneid_chr_dict[curr_geneid]
        if (curr_chr !="X"):
            tmp_4d_list=geneid_4d_dict[curr_geneid]
            #print("4DDDD")
            #print(len(tmp_4d_list))
            for curr_region in tmp_4d_list:
                tss_start=int(curr_region.split("_")[0].strip())
                tss_end=int(curr_region.split("_")[1].strip())
                func_file=open("/PATH_TO/TF_Combine_All_Function_chr"+curr_chr+".txt",'r')
                lines=func_file.readlines()
                for line in lines[1:]:
                    curr_b38="chr"+line.split("\t")[0].strip()+"_b38"
                    curr_pos=int(line.split("\t")[3].strip())
                    if (curr_pos >= tss_start) and (curr_pos <= tss_end):
                        if curr_geneid not in sel_variant_dict.keys():
                            sel_variant_dict[curr_geneid]=[]
                            sel_variant_anno_dict[curr_geneid]=[]
                        sel_variant_dict[curr_geneid].append(curr_b38)
                        sel_variant_anno_dict[curr_geneid].append("4D:"+curr_b38)
        func_file.close()
print("Check_4D_overlap_with_functional_variants")
print(len(sel_variant_dict.keys()))

output=open("TFs_Trans_Variants_All_2k.txt",'w')
output.write("Gene_ID\tTF_name\tGene_Name\tVariants")
for gene_id in sel_variant_dict.keys():
    output.write("\n"+gene_id.strip()+"\t"+geneid_name_dict[gene_id]+"\t"+geneid_name_dict[gene_id])
    tmp_var_list=[]
    for variants in sel_variant_dict[gene_id]:
        if variants not in tmp_var_list:
            output.write("\t"+variants.strip())
            tmp_var_list.append(variants)
output.close()












