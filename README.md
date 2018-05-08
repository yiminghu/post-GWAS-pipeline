# post-GWAS-pipeline
This repository is a pipeline built for post-GWAS analysis. With GWAS summary statistics as input, the pipeline will generate the following analysis:
* *Manhattan plot/QQ plot/[LocusZoom](https://genome.sph.umich.edu/wiki/LocusZoom_Standalone)*
* *[Heritability estimation](https://github.com/bulik/ldsc)/[annotation](http://genocanyon.med.yale.edu/GenoSkyline)-stratified enrichment analysis*
* *[GNOVA](https://github.com/xtonyjiang/GNOVA) genetic correlation estimation with 2419 UKB traits + 210 published GWAS*
* *[UTMOST](https://github.com/Joker-Jerome/UTMOST) cross-tissue gene-trait association analysis*

## Overview of the pipeline
<img src="./pipeline.png" width="900">

## Dependency
### R related (for Standard GWAS module and extracting results)
```bash
## with in R interface
install.packages('qqman')
install.packages('data.table')
install.packages('GWASTools')
```
### LocusZoom
Software [download](https://github.com/statgen/locuszoom-standalone)
Installation [wiki](http://genome.sph.umich.edu/wiki/LocusZoom_Standalone)

## Setup paths

## Reformatting Sumstats/QC

