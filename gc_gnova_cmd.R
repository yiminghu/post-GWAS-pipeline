rm(list = ls())
options(stringsAsFactors=F)


args = commandArgs(trailingOnly=TRUE)
aname = args[1] ## trait name
out_path = args[2] ## path to sumstats
sumstats = args[3] ## path to munged summary stats
gnova_path = args[4] ## gnova path: /ysm-gpfs/pi/zhao/from_louise/yh367/GNOVA/GNOVA-master
sqpath = args[5] ## output path for simple queue task list
ref_path = args[6] ##/gpfs/ysm/pi/zhao/from_louise/yh367/VA/workflow/GNOVA/sumstats/match_table_full.csv
bfile_path = args[7] ##/ysm-gpfs/pi/zhao/from_louise/yh367/GNOVA/genotype_1KG_eur_SNPmaf5/

mat_table_full = read.csv(ref_path, header = T)
dir.create(paste0(out_path, aname))
cmds = paste0('module load Python; python ', gnova_path, 'gnova.py ',sumstats,' ',mat_table_full$Curated_file,' --bfile ', bfile_path, 'eur_chr@_SNPmaf5 --out ',out_path, aname,'/',aname,'_',mat_table_full$Code,'_NS.txt;')    
writeLines(cmds, paste0(sqpath, aname, '_gnova_ns'))
