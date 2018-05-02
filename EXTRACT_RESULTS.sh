'''
EXTRACT GNOVA & UTMOST RESULTS
Require: gnova_output_path, gnova_output_prefix, utmost_joint_output_path, utmost_joint_output_prefix
'''
GENE_TABLE=/ysm-gpfs/pi/zhao/from_louise/bl537/GWAS/VA/metaxcan/GeneList.txt
## generate GNOVA summary table
Rscript --vanilla extract_gc_gnova_results.R ${gnova_output_path} ${out_prefix} ${gnova_output_prefix} ${REF_TABLE_PATH}

## generate UTMOST gene-level Manhattan plot/qq plot/significant gene table
utmost_results="${utmost_joint_output_path}/${utmost_joint_output_prefix}_GTEx_1_17290.txt"
Rscript --vanilla extract_utmost_results.R ${utmost_results} ${TEMP} ${TRAIT_NAME} ${GENE_TABLE}
