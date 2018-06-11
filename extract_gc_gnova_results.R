options(stringsAsFactors=F)

args = commandArgs(trailingOnly=TRUE)
ns_output = args[1] ## /gpfs/ysm/pi/zhao/from_louise/bl537/GWAS/VA/GNOVA/VA_UKB/0322/non-stratified/
out_prefix = args[2] ## prefix to summary table
trait = args[3] ## path to munged summary stats
ref_path = args[4] ##/gpfs/ysm/pi/zhao/from_louise/yh367/VA/workflow/GNOVA/sumstats/match_table_full.csv

u = read.csv(ref_path, header=T)
u = u[-76, -1] ## remove Apolipoprotein_A-I.Kettunen.2010 (problem with summary statistics)
list.disease = u$Code


#### No Need to Change ####
cutoff = 1e-3
fdr_cut = 0.01

res_table = as.data.frame(matrix(NA, length(list.disease), 6))
colnames(res_table) = c("Covariance", "Correlation", "P", "Sig_by_P", "FDR", "Sig_by_FDR")
Sig_by_P_list = c()
Sig_by_FDR_list = 1:length(list.disease)
pV = rep(-1, length(list.disease))

for (i in 1:length(list.disease)){
  gnova_output = paste0(ns_output, '/', trait, '_', list.disease[i], "_NS.txt")
  gnova_result = read.table(gnova_output, header=T)
  
  pval = gnova_result[1, "pvalue"]
  pV[i] = pval
  # covariance
  res_table[i, 1] = round(gnova_result[1, "rho"], digits = 3)
  
  # correlation 
  #if (trait[j]=='ipf')        {a[1, "h2_1"]  = 0.1281}
  if (list.disease[i]=='ALS') {gnova_result[1, "h2_2"]  = 0.085}
  res_table[i, 2] = round(gnova_result[1, "rho"]/sqrt(gnova_result[1, "h2_1"]*gnova_result[1, "h2_2"]), digits = 3)
  
  # p-value
  res_table[i, 3] = format(pval, scientific=T, digits=3)
  if (pval< cutoff){
    res_table[i, 4] = '*'
    Sig_by_P_list = c(Sig_by_P_list, i)
  }else
  {
    res_table[i, 4] = ''    
  }
  print(i)
}
fdrV = p.adjust(pV, method = "fdr")
Sig_by_FDR_list = Sig_by_FDR_list[fdrV<fdr_cut]
res_table[,5] = format(fdrV, scientific=T, digits=3)
res_table[fdrV<fdr_cut,6] = rep('*', sum(fdrV<fdr_cut))
res_table[fdrV>fdr_cut,6] = rep('', sum(fdrV>fdr_cut))
output_file = data.frame(list.disease, u$Categ, u$Trait, res_table)
colnames(output_file)[1:3] = c("Index", "Source", "Phenotype")
output_DIR = paste0(out_prefix, trait, ".csv")
#Sig_by_P_DIR = paste0(out_prefix, 'sig_by_p_', trait, ".txt")
#Sig_by_FDR_DIR = paste0(out_prefix, 'sig_by_fdr_', trait, ".txt")
write.csv(output_file, file=output_DIR, row.names = F)  
#write.table(output_file, file=output_DIR, quote = F, sep ="\t", row.names = F, col.names = T)  
#write.csv(u[Sig_by_P_list,], file=Sig_by_P_DIR, row.names = F)  
#write.csv(u[Sig_by_FDR_list,], file=Sig_by_FDR_DIR, row.names = F)  
print(paste0("GNOVA summary table written to ", output_DIR))

