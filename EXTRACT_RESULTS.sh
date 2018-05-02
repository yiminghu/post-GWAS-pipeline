'''
EXTRACT GNOVA & UTMOST RESULTS
Require: gnova_output_path, gnova_output_prefix, utmost_joint_output_path, utmost_joint_output_prefix
'''

## generate GNOVA summary table
Rscript --vanilla extract_gc_gnova_results.R ${gnova_output_path} ${out_prefix} ${gnova_output_prefix} ${REF_TABLE_PATH}

## generate UTMOST gene-level Manhattan plot/qq plot/significant gene table
Rscript --vanilla extract_utmost_results.R
