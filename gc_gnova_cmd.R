rm(list = ls())
options(stringsAsFactors=F)


args = commandArgs(trailingOnly=TRUE)
output_path = args[1] ## a folder to save all GNOVA output files (there will be ~2,500 files generated)
output_prefix = args[2] ## prefix of output files
sumstats = args[3] ## path to input munged summary stats
gnova_path = args[4] ## gnova path: /ysm-gpfs/pi/zhao/from_louise/yh367/GNOVA/GNOVA-master
sqpath = args[5] ## output path for simple queue task list
ref_path = args[6] ## reference table, header: Categ (source of published GWAS), 
## Trait (trait name), Curated_file (path to munged sumstats), Code (index for disease order)
## e.g. /gpfs/ysm/pi/zhao/from_louise/yh367/VA/workflow/GNOVA/sumstats/match_table_full.csv
## note this table requires updating if position of any munged sumstats changes. 
bfile_path = args[7] ## 1000 Genomes genotype data for estimating LD matrix
##e.g. /ysm-gpfs/pi/zhao/from_louise/yh367/GNOVA/genotype_1KG_eur_SNPmaf5/
print(args)

mat_table_full = read.csv(ref_path, header = T) ## read 
dir.create(output_path, showWarnings = FALSE) ## mkdir if not exist
existed_files <- dir(output_path)
system(paste0('rm -r ', sqpath, '*'))
system(paste0('rm -r ', 'core.*'))
cmds <- c()
for(i in 1:length(mat_table_full$Code)){
	if(paste0(output_prefix,'_',mat_table_full$Code[i], '_NS.txt') %in% existed_files){
		next
	}
	else{
		cmds <- c(cmds, paste0('source ~/.bashrc; python2 ', gnova_path, 'gnova.py ',sumstats,' ',mat_table_full$Curated_file[i],' --bfile ', bfile_path, 'eur_chr@_SNPmaf5 --out ',output_path,'/',output_prefix,'_',mat_table_full$Code[i],'_NS.txt;'))
	}
}
#print(cmds)
if(length(cmds) != 0){
	if(length(cmds) > 128*10){
		cuts_num <- ceiling(length(cmds)/11)
	}
	if(length(cmds) <= 128*10){
		cuts_num <- 128
	}
	newcmds <- split(cmds, ceiling((1:length(cmds))/cuts_num))
	for(i in 1:length(newcmds)){
		writeLines(newcmds[[i]], paste0(sqpath, '_part', i))
	}
	print(paste0('GNOVA task file written to ', sqpath))
}
if(length(cmds) == 0){
	print('GNOVA already finished!')
}
