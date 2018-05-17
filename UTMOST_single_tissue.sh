'''
UTMOST
SOFTWARE SETTINGS
## Note: TISSUE_LIST need align with WEIGHT_DB_PATH/SINGLE_TISSUE_COV_PATH/JOINT_TISSUE_COV_PATH
'''
UTMOST_PATH=/ysm-gpfs/pi/zhao/from_louise/yh367/software/UTMOST
TISSUE_LIST=(Adipose_Subcutaneous Adipose_Visceral_Omentum Adrenal_Gland Artery_Aorta Artery_Coronary Artery_Tibial Brain_Anterior_cingulate_cortex_BA24 Brain_Caudate_basal_ganglia Brain_Cerebellar_Hemisphere Brain_Cerebellum Brain_Cortex Brain_Frontal_Cortex_BA9 Brain_Hippocampus Brain_Hypothalamus Brain_Nucleus_accumbens_basal_ganglia Brain_Putamen_basal_ganglia Breast_Mammary_Tissue Cells_EBV-transformed_lymphocytes Cells_Transformed_fibroblasts Colon_Sigmoid Colon_Transverse Esophagus_Gastroesophageal_Junction Esophagus_Mucosa Esophagus_Muscularis Heart_Atrial_Appendage Heart_Left_Ventricle Liver Lung Muscle_Skeletal Nerve_Tibial Ovary Pancreas Pituitary Prostate Skin_Not_Sun_Exposed_Suprapubic Skin_Sun_Exposed_Lower_leg Small_Intestine_Terminal_Ileum Spleen Stomach Testis Thyroid Uterus Vagina Whole_Blood) 
WEIGHT_DB_PATH=$UTMOST_PATH/sample_data/weight_db_GTEx/
SINGLE_TISSUE_COV_PATH=$UTMOST_PATH/sample_data/covariance_tissue/
JOINT_TISSUE_COV_PATH=$UTMOST_PATH/sample_data/covariance_joint/
GENE_INFO=$UTMOST_PATH/intermediate/gene_info.txt

'''
SET INPUT & OUTPUT PATH
'''
sumstats_path=/ysm-gpfs/pi/zhao/from_louise/yh367/VA/workflow/qc_passed_with_rs/
sumstats_file=set04162018A.txt
single_tissue_output_path=/ysm-gpfs/pi/zhao/from_louise/yh367/VA/workflow/gene_level_analysis/UTMOST/suicide/
joint_output_path=/ysm-gpfs/pi/zhao/from_louise/yh367/VA/workflow/gene_level_analysis/UTMOST/suicide_joint/
joint_output_prefix=suicide_GBJ

'''
SET TASK FILE FOR SINGLE-TISSUE & JOINT TEST (for parallell running purpose)
'''
single_tissue_task_file=/ysm-gpfs/pi/zhao/from_louise/yh367/VA/workflow/gene_level_analysis/UTMOST/suicide_single_tissue.task
joint_task_file=/ysm-gpfs/pi/zhao/from_louise/yh367/VA/workflow/gene_level_analysis/UTMOST/suicide_joint.task

for tissue in ${TISSUE_LIST[@]}
do
task_single_tissue="python2 ${UTMOST_PATH}/single_tissue_association_test.py --model_db_path ${WEIGHT_DB_PATH}/${tissue}.db --covariance ${SINGLE_TISSUE_COV_PATH}/${tissue}.txt.gz --gwas_folder ${sumstats_path} --gwas_file_pattern ${sumstats_file} --snp_column SNP --effect_allele_column ALT --non_effect_allele_column REF --beta_column ALT_EFFSIZE --pvalue_column PVALUE --output_file ${single_tissue_output_path}/${tissue}.csv"
echo $task_single_tissue >> ${single_tissue_task_file}
done
task_joint="python2 joint_GBJ_test.py --weight_db ${WEIGHT_DB_PATH}--output_dir ${joint_tissue_output_path} --cov_dir ${JOINT_TISSUE_COV_PATH} --input_folder ${single_tissue_output_path} --gene_info ${GENE_INFO} --output_name ${joint_output_prefix} --start_gene_index 1 --end_gene_index 17290"
echo $task_joint >> ${joint_task_file}

echo "SINGLE-TISSUE tasks written to ${task_single_tissue}!"
echo "JOINT TEST tasks written to ${joint_task_file}!"
