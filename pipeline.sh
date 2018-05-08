'''
GLOBAL parameters
'''
PIPELINE_PATH=/input the path you cloned this repo here/
LOCUSZOOM_PATH=${PIPELINE_PATH}/post-GWAS-pipeline/locuszoom/
LDSC_PATH=${PIPELINE_PATH}/post-GWAS-pipeline/ldsc/
GNOVA_PATH=${PIPELINE_PATH}/post-GWAS-pipeline/gnova/
UTMOST_PATH=${PIPELINE_PATH}/post-GWAS-pipeline/UTMOST/
REF_TABLE_PATH=/gpfs/ysm/pi/zhao/from_louise/yh367/VA/workflow/GNOVA/sumstats/match_table_full.csv
'''
Input GWAS parameters
'''
TRAIT_NAME=example
sumstats=/gpfs/ysm/pi/zhao/from_louise/bl537/GWAS/VA/Pipeline/Grand/example/SumStats.txt
grandfolder=${PIPELINE_PATH}/post-GWAS-pipeline/grandfolder/ ##path for saving 

## Generate Standard GWAS / munge / LDSC task files ##
cd ${PIPELINE_PATH}/post-GWAS-pipeline/
python2.7 standardGWAS_munge_ldsc.py --study ${TRAIT_NAME} --sumstats ${sumstats} --grandfolder ${grandfolder} --ldsc_path ${LDSC_PATH} --locuszoom_path ${LOCUSZOOM_PATH}
## this step will generate the following files: 

## Generate GNOVA task files ##
REF_GENOME_PATH=${GNOVA_PATH}/genotype_1KG_eur_SNPmaf5/
GNOVA_results=${grandfolder}/GNOVA_results/
mkdir ${GNOVA_results}
gnova_output_path=${GNOVA_results}/${TRAIT_NAME}/ #1 a folder to save all GNOVA output files (there will be ~2,500 files generated)
mkdir ${gnova_output_path}
gnova_output_prefix=${TRAIT_NAME} #2 prefix of output files
sumstats=${MUNGED_PATH} #3 path to input munged summary stats
task_file="${TASK_FOLDER}/gnova_ns.task" #5 output path for simple queue task list

Rscript --vanilla gc_gnova_cmd.R ${grandfolder} ${gnova_output_prefix} ${sumstats} ${GNOVA_PATH} ${task_file} ${REF_TABLE_PATH} ${REF_GENOME_PATH}
