'''
GLOBAL parameters
'''


'''
Input GWAS parameters
'''
TRAIT_NAME=example
sumstats=/gpfs/ysm/pi/zhao/from_louise/bl537/GWAS/VA/Pipeline/Grand/example/SumStats.txt
grandfolder=/gpfs/ysm/pi/zhao/from_louise/bl537/GWAS/VA/Pipeline/Grand

## No need to change ##
cd /gpfs/ysm/pi/zhao/from_louise/bl537/GWAS/VA/Pipeline/Script
python2.7 Pipeline.py --study ${study} --sumstats ${sumstats} --grandfolder ${grandfolder}

