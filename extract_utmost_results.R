options(stringsAsFactors=F)
library("qqman")
library("data.table")
library("GWASTools")
args = commandArgs(trailingOnly=TRUE)
utmost_results = args[1] ## path to utmost joint test output
out_prefix = args[2] ## prefix to gene-level Manhattan plot/qq plot/significant gene table
gene_table = args[3] ## gene position info


## joint GBJ results summary
a = as.data.frame(fread(paste0(utmost_results), header=T, sep="\t"))
#gene    test_score      p_value
a = unique(a)
a = final[complete.cases(a), ]
thres = 0.05/nrow(a)
write.csv(a[a$pvalue<thres,], paste0(out_prefix,'/utmost_joint_sig_genes.csv'), row.names=F)


##### plot #####
gene_table = "/ysm-gpfs/pi/zhao/from_louise/bl537/GWAS/VA/metaxcan/GeneList.txt"
g = read.table(gene_table, header=T, sep="\t")
d = a[,c(1,3)]
colnames(d) = c("gene_name", "pvalue")  
dd = merge(d, g, all.x=T, by.x="gene_name", by.y="GENE")
  	
dd = dd[order(dd$CHR, dd$POS),]
a = dd
m = a[complete.cases(a),]
colnames(m) = c("Gene", "P", "CHR", "BP")
m$P = as.numeric(m$P)
mdir = paste0(out_prefix, "/utmost_joint_manhattan.png")
png(mdir, width=1000, height=600, res=100)
manhattan(m, chr = "CHR", bp = "BP", p = "P",
          col = c("#4682B4", "#B4464B"), chrlabs = as.character(c(1:22)), main=aname,
          suggestiveline = FALSE, genomewideline = -log10(thres),
          ylim = c(0, max(-log10(m$P), -log10(thres))+2),
          logp = TRUE,
          las = 2,
          cex = 0.5)
dev.off()

qdir = paste0(out_prefix, "/utmost_joint_qq.png")
png(qdir, width=1000, height=800, res=120)
qqPlot(m$P, truncate = FALSE, ylim = NULL, ci=TRUE, main=aname)
dev.off()



################################## Single-tissue results summary ##################################
#gf = c("csp572_bipol_eu", "csp572_cogni_eu", "csp572_schiz_eu")
#fd = paste0("/ysm-gpfs/pi/zhao/from_louise/yh367/VA/metaxcan_results1/", gf, '/')
#traits = c("bipol_eu", "cogni_eu", "schiz_eu")
#outfolder = paste0("/ysm-gpfs/pi/zhao/from_louise/yh367/VA/workflow/gene_level_analysis/MetaXcan/summary/",traits)
#sapply(outfolder, dir.create)
#thres = 2.43925e-07
#tt = "csp572"


## by tissue
#for(i in 1:length(gf)){
#	aname = traits[i]
#	folder = paste0(fd[i])
#	Files.csv = list.files(folder)
#	Files.csv = Files.csv[grepl(".csv", Files.csv)]
#	
#	pb = txtProgressBar(0, length(Files.csv), style=3)
#	for(j in 1:length(Files.csv)){
#	  d = as.data.frame(fread(paste0(folder, '/',Files.csv[j]), header=T, sep=","))
#	  #colnames(d) = c("gene_name", "pvalue")
#	  d = na.omit(d)
#	  d = unique(d)
#	  d = d[d$pvalue<thres, ]
#	  if(nrow(d)>0){
#	  	d['Tissue'] = Files.csv[j]
#	  }
#	  if (j==1) {a=d} else{
#	    a = rbind(a, d)
#	  }
#	  setTxtProgressBar(pb, j)
#	}
##	write.table(a, paste0(outfolder[i],'/',tt,'.',aname,'_by-tissue_sig_genes.txt'), quote=F, sep='\t', row.names=F, col.names=T)
#	write.table(a, paste0(outfolder[i],'/',tt,aname,'_by-tissue_sig_genes.txt'), quote=F, sep='\t', row.names=F, col.names=T)
#
#
#	##### plot #####
#	gene_table = "/ysm-gpfs/pi/zhao/from_louise/bl537/GWAS/VA/metaxcan/GeneList.txt"
#	g = read.table(gene_table, header=T, sep="\t")
#
#	pb = txtProgressBar(0, length(Files.csv), style=3)
#	for(j in 1:length(Files.csv)){
#	  d = as.data.frame(fread(paste0(folder, '/',Files.csv[j]), header=T, sep=","))
#	  d = na.omit(d)
#	  d = unique(d)
#	  d = d[, c("gene_name", "pvalue")]
#	  
#	  dd = merge(d, g, all.x=T, by.x="gene_name", by.y="GENE")
#	  
#	  if (j==1) {a=dd} else{
#	    a = rbind(a, dd)
#	  }
#	  setTxtProgressBar(pb, j)
#	}
#	
#	aa = a[order(a$CHR, a$POS),]
#	a = aa
#	m = a[complete.cases(a),]
#	colnames(m) = c("Gene", "P", "CHR", "BP")
#	  
#	mdir = paste0(outfolder[i], "/",tt,aname, "_by-tissue_manhattan.png")
##	mdir = paste0(outfolder[i], "/",tt,'.', aname, "_by-tissue_manhattan.png")
#	png(mdir, width=1000, height=600, res=100)
#	manhattan(m, chr = "CHR", bp = "BP", p = "P",
#	          col = c("#4682B4", "#B4464B"), chrlabs = as.character(c(1:22)), main=aname,
#	          suggestiveline = FALSE, genomewideline = -log10(thres),
#	          ylim = c(0, max(-log10(m$P), -log10(thres))+2),
#	          logp = TRUE,
#	          las = 2,
#	          cex = 0.5)
#	dev.off()
#
#}




#options(stringsAsFactors = F)
#tissue_helper <- function(tis){
#  tmp = strsplit(tis, '\\.')[[1]][1]
#  tmp = unlist(strsplit(tmp, '_'))
#  tmp = tmp[-1]
#  paste(tmp, collapse = '_')
#}
##xxx = c("alcamt", "pclc", "slprs")
##ttV = list()
##ttV[[1]] = c("alcamt.eu", "alcamt.g1", "alcamt.g2")
##ttV[[2]] = c("pclc.eu", "pclc.eu_cmbt", "pclc.eu_noncmbt")
##ttV[[3]] = c("slprs.13_eu", "slprs.64_eu")
##setwd("/ysm-gpfs/pi/zhao/from_louise/yh367/VA/workflow/gene_level_analysis/MetaXcan/summary/new_ones")
##tt = "new_ones.CasePclcEA"
##tt = "new_ones.slprEA"
##setwd("/ysm-gpfs/pi/zhao/from_louise/yh367/VA/workflow/gene_level_analysis/MetaXcan/summary/ptsd-exp-ea113017")
##tt = "ptsd-exp-ea113017"
#
##traits = c("e-c-ea", "une-b-ea", "une-c-ea")
##outfolder = paste0("/ysm-gpfs/pi/zhao/from_louise/yh367/VA/workflow/gene_level_analysis/MetaXcan/summary/",traits, '/')
##tt = "ptsd"
#
#traits = c("csp572_bipol_eu", "csp572_cogni_eu", "csp572_schiz_eu")
#outfolder = paste0("/ysm-gpfs/pi/zhao/from_louise/yh367/VA/workflow/gene_level_analysis/MetaXcan/summary/",traits, '/')
#
#for(i in 1:length(outfolder)){
#	setwd(outfolder[i])
#	aname = traits[i]
##	df = read.table(paste0(outfolder[i],'/',tt,'.',aname,'_by-tissue_sig_genes.txt'), header=T)
#	df = read.table(paste0(outfolder[i],'/',aname,'_by-tissue_sig_genes.txt'), header=T)
#	if(nrow(df)==0){
#	  print(paste0(tt, " empty!"))
#	}else{
#		df = df[,c(2,3,5,13)]
#		for(i in 1:nrow(df)){
#		  df[i,4] = tissue_helper(df[i,4])
#		}
#		df = df[order(df[,1]),]
#		write.table(df, paste0(aname,'.full.txt'), sep='\t', col.names=T, row.names=F, quote=F)
##		write.table(df, paste0(tt,'.',aname,'.full.txt'), sep='\t', col.names=T, row.names=F, quote=F)
#		gset = unique(df[,1])
#		most_sig = data.frame()
#		for(g in gset){
#		  chunk = df[df[,1]==g, ]
#		  insert = chunk[which.min(chunk[,3]),]
#		  most_sig = rbind(most_sig, insert)
#		}
#		write.table(most_sig, paste0(aname,'.most.sig.txt'), sep='\t', col.names=T, row.names=F, quote=F)
##		write.table(most_sig, paste0(tt,'.',aname,'.most.sig.txt'), sep='\t', col.names=T, row.names=F, quote=F)
#	}
#}
#
#