'''
GNOVA 
SOFTWARE SETTINGS
'''
GNOVA_PATH=/ysm-gpfs/pi/zhao/from_louise/yh367/GNOVA/GNOVA-master #4
REF_TABLE_PATH=/gpfs/ysm/pi/zhao/from_louise/yh367/VA/workflow/GNOVA/sumstats/match_table_full.csv #6
## reference table, header: Categ (source of published GWAS), 
## Trait (trait name), Curated_file (path to munged sumstats), Code (index for disease order)
## e.g. /gpfs/ysm/pi/zhao/from_louise/yh367/VA/workflow/GNOVA/sumstats/match_table_full.csv
## note this table requires updating if position of any munged sumstats changes. 
REF_GENOME_PATH=/ysm-gpfs/pi/zhao/from_louise/yh367/GNOVA/genotype_1KG_eur_SNPmaf5/  #7
## 1000 Genomes genotype data for estimating LD matrix
##e.g. /ysm-gpfs/pi/zhao/from_louise/yh367/GNOVA/genotype_1KG_eur_SNPmaf5/

'''
INPUT & OUTPUT SETTINGS
'''

gnova_output_path="${TEMP}/gnova_ns/" #1 a folder to save all GNOVA output files (there will be ~2,500 files generated)
mkdir gnova_output_path
gnova_output_prefix=${TRAIT_NAME} #2 prefix of output files
sumstats=${MUNGED_PATH} #3 path to input munged summary stats
task_file="${TASK_FOLDER}/gnova_ns.task" #5 output path for simple queue task list

'''
TASK FILES SETTINGS (for parallell running purpose)
'''
Rscript --vanilla gc_gnova_cmd.R ${gnova_output_path} ${gnova_output_prefix} ${sumstats} ${GNOVA_PATH} ${task_file} ${REF_TABLE_PATH} ${REF_GENOME_PATH}
