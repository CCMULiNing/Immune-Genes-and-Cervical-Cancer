#SMR analysis----
#SMR software (version 1.03)
#SMR software download website: https://yanglab.westlake.edu.cn/software/smr/#Download
smr --bfile reference_file --gwas-summary gwas.ma --beqtl-summary eqtl --out outcome


#Enrichment analysis: https://metascape.org/gp/index.html#/main/step1


#Multi-SNP based SMR analysis----
#SMR software (version 1.03)
smr --bfile reference_file --gwas-summary gwas.ma --beqtl-summary eqtl --out outcome --smr-multi


#MR analysis based on inverse variance weighted method----
#TwoSampleMR package (version 0.5.6) in R (version 4.1.2)
library(TwoSampleMR)
exp<-read.csv("exp_genetic_variants.csv")
out<-read.csv("out.csv")
exp1<-format_data(exp,snp_col = "SNP",
                  beta_col = "b", se_col = "SE",
                  effect_allele_col = "A1",other_allele_col = "A2",
                  eaf_col = "Freq",pval_col = "p")
out1<-format_data(out,type='outcome',snps=exp1$SNP,
                  snp_col = "SNP", effect_allele_col = "A1",
                  other_allele_col = "A2", eaf_col = "freq",
                  beta_col = "b", se_col = "se",
                  pval_col = "p")
har<-harmonise_data(exposure_dat = exp1,outcome_dat = out1)
res<-mr(har)


#Colocalization analysis----
#coloc package (version, 5.1.0.1) in R (version 4.1.2)

library(coloc)
#1. extract eQTL summary data for interesting genes using SMR
smr --beqtl-summary eqtl --query 1 --gene gene_id --cis-wind 1000 --out myquery
#GWAS summary data is in gcta format (SNP, A1, A2, frq, b, se, p, N)
#eQTL summary data is in SMR format (SNP, Chr, BP, A1, A2, Freq, Probe, Probe_Chr, Probe_bp, Gene, Orientation, b, SE, p)

#2. extract same SNPs in both datasets
gwas_data<-fread("gwasdata.ma",quote="")
geneeqtl_data<-read.table("myquery.txt", header=T, stringsAsFactor=F)
snps = intersect(gwas_data$SNP,geneeqtl_data$SNP)
smrfile<-smrfile[which(geneeqtl_data$SNP %in% snps),]
gwas_data1<-gwas_data[which(gwas_data$SNP %in% snps),]

#3. create list of input parameters for gwas data
MAF1=ifelse(gwas_data1$freq<0.5,gwas_data1$freq, 1-gwas_data1$freq)
gwas_data2=list("snp"=snps,"MAF"=MAF1,"beta"=gwas_data1$b,"varbeta"=gwas_data1$se^2,"pvalues"=gwas_data1$p,"N"=gwas_data1$n, "type"="cc", "s"=0.019126)
#note: "type"="cc"   Type of trait - cc for case/ctrl
#note: "s"= 0.019126   proportion of cases

#4. create list of input parameters for eqtl data
MAF2=ifelse(smrfile$Freq<0.5,smrfile$Freq, 1-smrfile$Freq)
smrfile1=list("snp"=smrfile$SNP,"MAF"=MAF2,"beta"=smrfile$b,"varbeta"=smrfile$SE^2,"pvalues"=smrfile$p,"N"=982, "type"="quant")
#note: "type"="quant"  quantitative trait

#5. formal colocalization analysis
res=coloc.abf(gwas_data2, smrfile1, p1 = 1e-04, p2 = 1e-04, p12 = 1e-05)


#SMR locus plot----
#1. use SMR software to generate a data file for plot
smr --bfile reference_file --gwas-summary gwas.ma --beqtl-summary eqtl --out myplot --plot --probe gene --probe-wind 1000 --gene-list glist-hg19

#2. use R (version 4.1.2) to draw the SMR locus plot
#download "ploy_SMR.r": https://yanglab.westlake.edu.cn/software/smr/#RscriptforSMRlocusplot
source("plot_SMR.r")
SMRData = ReadSMRData("myplot.gene.txt")
SMRLocusPlot(data=SMRData, smr_thresh=0.05, heidi_thresh=0.05, plotWindow=1000, max_anno_probe=5)