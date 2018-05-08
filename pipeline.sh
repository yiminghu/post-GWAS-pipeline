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
grandfolder=/gpfs/ysm/pi/zhao/from_louise/bl537/GWAS/VA/Pipeline/Grand ##path for saving 

## No need to change ##
cd /gpfs/ysm/pi/zhao/from_louise/bl537/GWAS/VA/Pipeline/Script
python2.7 Pipeline.py --study ${study} --sumstats ${sumstats} --grandfolder ${grandfolder}

