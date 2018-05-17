# post-GWAS-pipeline
This repository is a pipeline built for post-GWAS analysis. With GWAS summary statistics as input, the pipeline contains four major modules:
* *Standard GWAS analysis: Manhattan plot/QQ plot/[LocusZoom](https://genome.sph.umich.edu/wiki/LocusZoom_Standalone)*
* *[Heritability estimation](https://github.com/bulik/ldsc)/[annotation](http://genocanyon.med.yale.edu/GenoSkyline)-stratified enrichment analysis*
* *[GNOVA](https://github.com/xtonyjiang/GNOVA) genetic correlation estimation with 2419 UKB traits + 210 published GWAS*
* *[UTMOST](https://github.com/Joker-Jerome/UTMOST) cross-tissue gene-trait association analysis*

## Overview of the pipeline
<img src="./pipeline.png" width="900">

## Tutorial
### Clone the repo
```bash
PIPELINE_PATH=/input the path you want to install pipeline here/
cd ${PIPELINE_PATH}
git clone https://github.com/yiminghu/post-GWAS-pipeline.git
cd ${PIPELINE_PATH}/post-GWAS-pipeline
```
### 1. Install dependency
#### 1.1 R related (for Standard GWAS module and extracting results)
```bash
## with in R interface
install.packages('qqman')
install.packages('data.table')
install.packages('GWASTools')
```

#### 1.2 LocusZoom
Software [download](https://github.com/statgen/locuszoom-standalone) and [wiki](http://genome.sph.umich.edu/wiki/LocusZoom_Standalone).
```bash
cd ${PIPELINE_PATH}/post-GWAS-pipeline
wget https://statgen.sph.umich.edu/locuszoom/download/locuszoom_1.4.tgz
tar xvzf locuszoom_1.4.tgz
```
#### 1.3 LDSC
More info can be found on https://github.com/bulik/ldsc.
```bash
cd ${PIPELINE_PATH}/post-GWAS-pipeline
## install and setup python dependency
git clone git@github.com:bulik/ldsc.git
cd ./ldsc
conda env create --file environment.yml
source activate ldsc
```
Download ld scores
```bash
mkdir Input
mkdir Input/EUR
mkdir Input/EUR/genotype
mkdir Input/EUR/weights
cd Input/EUR
wget https://data.broadinstitute.org/alkesgroup/LDSCORE/eur_w_ld_chr.tar.bz2
tar -jxvf eur_w_ld_chr.tar.bz2
cd ../genotype
wget https://data.broadinstitute.org/alkesgroup/LDSCORE/1000G_Phase1_frq.tgz
tar xvzf 1000G_Phase1_frq.tgz
cd ../weights
wget https://data.broadinstitute.org/alkesgroup/LDSCORE/weights_hm3_no_hla.tgz
tar xvzf weights_hm3_no_hla.tgz
```

Download annotations
```bash
cd ${PIPELINE_PATH}/post-GWAS-pipeline/ldsc ## return to ldsc/
mkdir Annotations
mkdir Annotations/EUR
mkdir Annotations/EUR/Baseline
mkdir Annotations/EUR/GenoSkyline_Plus
## download baseline annotations
cd Annotations/EUR/Baseline
wget https://data.broadinstitute.org/alkesgroup/LDSCORE/1000G_Phase1_baseline_ldscores.tgz
## download GenoSkyline_Plus annotations
cd ../GenoSkyline_Plus
wget http://genocanyon.med.yale.edu/GenoSkylineFiles/GenoSkylinePlus/LD_score_files_1KGphase3.tar.gz
tar xvzf LD_score_files_1KGphase3.tar.gz
```

#### 1.4 GNOVA
More usage instructions can be found at https://github.com/xtonyjiang/GNOVA.
```bash
cd ${PIPELINE_PATH}/post-GWAS-pipeline/ ## return to post-GWAS-pipeline/
git clone git@github.com:xtonyjiang/GNOVA.git
## install python dependency if needed
pip install numpy --user
pip install scipy --user
pip install pandas --user
pip install sklearn --user
pip install bitarray --user
```
Download reference genome file
```bash
cd GNOVA
wget http://genocanyon.med.yale.edu/GNOVAFiles/genotype_1KG_eur_SNPmaf5.tar.gz
tar xvzf genotype_1KG_eur_SNPmaf5.tar.gz
```

#### 1.5 UTMOST
More usage instructions can be found at https://github.com/Joker-Jerome/UTMOST.
```bash
cd ${PIPELINE_PATH}/post-GWAS-pipeline/ ## return to post-GWAS-pipeline/
git clone https://github.com/Joker-Jerome/UTMOST.git
## install dependency
pip install numpy --user
pip install scipy --user
pip install pandas --user
pip install -Iv rpy2==2.8.6 --user
```
Download pre-calculated cross-tissue imputation models and covariance matrices
```bash
cd ./UTMOST
wget --load-cookies /tmp/cookies.txt "https://drive.google.com/uc?export=download&confirm=$(wget --quiet --save-cookies  /tmp/cookies.txt --keep-session-cookies --no-check-certificate 'https://drive.google.com/uc?export=download&id=1u8CRwb6rZ-gSPl89qm3tKpJArUT8XrEe' -O- | sed -rn 's/.*confirm=([0-9A-Za-z_]+).*/\1\n/p')&id=1u8CRwb6rZ-gSPl89qm3tKpJArUT8XrEe" -O sample_data.zip && rm -rf /tmp/cookies.txt
unzip sample_data.zip
cd sample_data
wget --load-cookies /tmp/cookies.txt "https://drive.google.com/uc?export=download&confirm=$(wget --quiet --save-cookies  /tmp/cookies.txt --keep-session-cookies --no-check-certificate 'https://drive.google.com/uc?export=download&id=1Kh3lHyTioKIXqCsREmsAyC-dS49KVO9G' -O- | sed -rn 's/.*confirm=([0-9A-Za-z_]+).*/\1\n/p')&id=1Kh3lHyTioKIXqCsREmsAyC-dS49KVO9G" -O covariance_tissue.tar.gz && rm -rf /tmp/cookies.txt
wget --load-cookies /tmp/cookies.txt "https://drive.google.com/uc?export=download&confirm=$(wget --quiet --save-cookies  /tmp/cookies.txt --keep-session-cookies --no-check-certificate 'https://drive.google.com/uc?export=download&id=1tqIW5Ms8p1StX7WXXWVa4TGKb5q58TPA' -O- | sed -rn 's/.*confirm=([0-9A-Za-z_]+).*/\1\n/p')&id=1tqIW5Ms8p1StX7WXXWVa4TGKb5q58TPA" -O covariance_joint.zip && rm -rf /tmp/cookies.txt
tar -zxvf covariance_tissue.tar.gz
unzip covariance_joint.zip
```

### 2. Setup global parameters (paths)
```bash
LOCUSZOOM_PATH=${PIPELINE_PATH}/post-GWAS-pipeline/locuszoom/
LDSC_PATH=${PIPELINE_PATH}/post-GWAS-pipeline/ldsc/
GNOVA_PATH=${PIPELINE_PATH}/post-GWAS-pipeline/gnova/
UTMOST_PATH=${PIPELINE_PATH}/post-GWAS-pipeline/UTMOST/
REF_TABLE_PATH=/gpfs/ysm/pi/zhao/from_louise/yh367/VA/workflow/GNOVA/sumstats/match_table_full.csv
REF_GENOME_PATH=${GNOVA_PATH}/genotype_1KG_eur_SNPmaf5/
grandfolder=${PIPELINE_PATH}/post-GWAS-pipeline/grandfolder/
mkdir ${grandfolder}
mkdir ${grandfolder}/joblist/
mkdir ${grandfolder}/LDSC/
mkdir ${grandfolder}/LDSC/sumstats/
mkdir ${grandfolder}/LDSC/Results/
mkdir ${grandfolder}/Standard/
mkdir ${grandfolder}/Standard/ManhattanPlot/
mkdir ${grandfolder}/Standard/qqPlot/
mkdir ${grandfolder}/Standard/Significant/
mkdir ${grandfolder}/Standard/LocusZoom/
mkdir ${grandfolder}/GNOVA/
mkdir ${grandfolder}/UTMOST/
```
**Note: REF_TABLE_PATH is a table contains paths to munged summary statistics of all reference GWAS (2,419 UKB + 210 publicly available GWAS). If any new reference GWAS are added in the future or any of the existing munged summary statistics are changed, please update this table accordingly!**

### 3. Reformatting Sumstats / QC
The pipeline assumes the following input format:

|CHROM  |POS    |A1     |A2     |BETA      |P        | SNP      |N      |
|:-----:|:-----:|:-----:|:-----:| :-----:| :-----:| :-----:| :-----:|
|1      |751343 |T      |A      |0.0477592 |0.0571189|rs28544273| 137044|
|1      |751756 |T      |C      |0.0476388 |0.0566586|rs28527770| 137044|
|1      |752894 |T      |C      |-0.0479547|0.0577561| rs3131971| 137044|
|1      |753405 |C      |A      |-0.0495965|0.0420519| rs3115860| 137044|
|1      |753425 |T      |C      |-0.0495346|0.0418289| rs3131970| 137044|

We also provide reformatting script for certain type of input (to be updated).

### 4. Generating task lists for all modules
```bash
jobfolder=${grandfolder}/joblist/${TRAIT_NAME}/
mkdir ${jobfolder}
## Generate Standard GWAS / munge / LDSC task files ##
cd ${PIPELINE_PATH}/post-GWAS-pipeline/
python2.7 standardGWAS_munge_ldsc.py --study ${TRAIT_NAME} --sumstats ${sumstats} --grandfolder grandfolder --ldsc_path ${LDSC_PATH} --locuszoom_path ${LOCUSZOOM_PATH}
### generates the following files in jobfolder:
### "standardGWAS_${TRAIT_NAME}"
### "munge_${TRAIT_NAME}"
### "annotation_enrichment_${TRAIT_NAME}_Tier1"
### "annotation_enrichment_${TRAIT_NAME}_Tier2"
### "annotation_enrichment_${TRAIT_NAME}_Tier3"
### "ldsc_${TRAIT_NAME}_Heritability"
## Generate GNOVA task files ##
GNOVA_results=${grandfolder}/GNOVA_results/
mkdir ${GNOVA_results}
gnova_output_path=${GNOVA_results}/${TRAIT_NAME}/
mkdir ${gnova_output_path}
gnova_output_prefix=${TRAIT_NAME} #2 prefix of output files
sumstats=${grandfolder}/LDSC/sumstats/${TRAIT_NAME}.sumstats.gz
task_file=${jobfolder}/gnova_ns.task
Rscript --vanilla gc_gnova_cmd.R ${grandfolder} ${gnova_output_prefix} ${sumstats} ${GNOVA_PATH} ${task_file} ${REF_TABLE_PATH} ${REF_GENOME_PATH}
```
### 5. Generating summary for GNOVA and UTMOST

