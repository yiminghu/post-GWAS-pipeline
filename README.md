# post-GWAS-pipeline
This repository is a pipeline built for post-GWAS analysis. With GWAS summary statistics as input, the pipeline contains four major modules:
* *Standard GWAS analysis: Manhattan plot/QQ plot/[LocusZoom](https://genome.sph.umich.edu/wiki/LocusZoom_Standalone)*
* *[Heritability estimation](https://github.com/bulik/ldsc)/[annotation](http://genocanyon.med.yale.edu/GenoSkyline)-stratified enrichment analysis*
* *[GNOVA](https://github.com/xtonyjiang/GNOVA) genetic correlation estimation with 2419 UKB traits + 210 published GWAS*
* *[UTMOST](https://github.com/Joker-Jerome/UTMOST) cross-tissue gene-trait association analysis*

## Overview of the pipeline
<img src="./pipeline.png" width="900">

## Quick-start
### Clone the repo
```bash
git clone https://github.com/yiminghu/post-GWAS-pipeline.git
cd post-GWAS-pipeline
```
### Install dependency
#### 1. R related (for Standard GWAS module and extracting results)
```bash
## with in R interface
install.packages('qqman')
install.packages('data.table')
install.packages('GWASTools')
```

#### 2. LocusZoom
Software [download](https://github.com/statgen/locuszoom-standalone) and [wiki](http://genome.sph.umich.edu/wiki/LocusZoom_Standalone).

#### 3. LDSC
More info can be found on https://github.com/bulik/ldsc.
```bash
## install and setup python dependency
git clone git@github.com:bulik/ldsc.git
cd ldsc
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
cd ../.. ## return to ldsc/
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

### GNOVA
See detailed instruction at https://github.com/xtonyjiang/GNOVA.

### UTMOST
See detailed instruction at https://github.com/Joker-Jerome/UTMOST.

## Setup paths

## Reformatting Sumstats/QC

## Generating task lists for all modules

## Generating summary for GNOVA and UTMOST
 
