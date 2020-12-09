## coord step ##
python coord_genotypes.py --gf=/net/zhao/yh367/Kaiser/t2d/eur_t2d_sel2 --ssf=/net/zhao/yh367/PleioPred/t2d_cad_fg/t2d_merged.sumstats --ssf_format=BASIC --N=69033 --out=/net/zhao/yh367/Kaiser/t2d/t2d_coord_sel2 --vgf=/net/zhao/yh367/Kaiser/t2d/eur_t2d_sel2 --vbim=/net/zhao/yh367/Kaiser/t2d/eur_t2d_sel2.bim\
python coord_genotypes.py --gf=/net/zhao/yh367/Kaiser/t2d/eur_t2d_sel2 --ssf=/net/zhao/yh367/PleioPred/t2d_cad_fg/cad_merged.sumstats --ssf_format=BASIC --N=86995 --out=/net/zhao/yh367/Kaiser/t2d/cad_coord_sel2 --vgf=/net/zhao/yh367/Kaiser/t2d/eur_t2d_sel2 --vbim=/net/zhao/yh367/Kaiser/t2d/eur_t2d_sel2.bim\

python prior_generating.py --h5py_file=/net/zhao/yh367/Kaiser/t2d/t2d_coord_sel2 --LDSC_results_file=/net/zhao/ql68/Software/ldsc/Results/T2D/T2D_DIAGRAMv3_Curated_GC1_GS7_withBaseline.results --output_h2=/net/zhao/yh367/Kaiser/t2d/priors/h2_sel2_file.txt --output_pT=/net/zhao/yh367/Kaiser/t2d/priors/pT_sel2 --PS=1,0.3,0.1,0.03,0.01,0.003,0.001,0.0003,0.0001,0.00003,0.00001\
python PleioPrior_generating.py --h5py_file1=/net/zhao/yh367/Kaiser/t2d/t2d_coord_sel2 --LDSC_results_file1=/net/zhao/ql68/Software/ldsc/Results/T2D/T2D_DIAGRAMv3_Curated_GC1_GS7_withBaseline.results --h5py_file2=/net/zhao/yh367/Kaiser/t2d/cad_coord_sel2 --LDSC_results_file2=/net/zhao/ql68/Software/ldsc/Results/CAD/CAD_CARDIoGRAM_Curated_GC1_GS7_withBaseline.results --output_anno_h2=/net/zhao/yh367/Kaiser/t2d/priors/t2d_cad_anno_sel2_file.txt --output_ld_h2=/net/zhao/yh367/Kaiser/t2d/priors/t2d_cad_ld_sel2_file.txt\
python PleioPred_no_comp.py --coord_D1=/net/zhao/yh367/Kaiser/t2d/t2d_coord_sel2 --coord_D2=/net/zhao/yh367/Kaiser/t2d/cad_coord_sel2 --N1=69033 --N2=86995 --rho=0.1 --local_ld_prefix=/net/zhao/yh367/Kaiser/t2d/t2d_cad_sel2 --ld_radius=115 --hfile=/net/zhao/yh367/Kaiser/t2d/priors/t2d_cad_anno_sel2_file.txt --out=/net/zhao/yh367/Kaiser/t2d/output/pleio_inf/0.1/t2d_cad_sel2\

## PleioPred-non-infinitesimal ##
python PleioPred_bi_rho.py --coord_D1=/net/zhao/yh367/Kaiser/t2d/t2d_coord_sel1 --coord_D2=/net/zhao/yh367/Kaiser/t2d/cad_coord_sel1 --rho=0.5 --N1=69033 --N2=86995 --zero_jump_prob=0.05 --num_iter=500 --burn_in=200 --init_PV=0.25,0.25,0.25,0.25 --init_betas=/net/zhao/yh367/Kaiser/t2d/priors/t2d_cad_anno_initial_sel1.pickled.gz --alpha=1,1,1,1 --local_ld_prefix=/net/zhao/yh367/Kaiser/t2d/t2d_cad_sel1 --ld_radius=115 --hfile=/net/zhao/yh367/Kaiser/t2d/priors/t2d_cad_anno_sel1_file.txt --out=/net/zhao/yh367/Kaiser/t2d/output/pleio_non_inf/0.5/t2d_cad_sel1\
python PleioPred_bi_rho.py --coord_D1=/net/zhao/yh367/Kaiser/t2d/t2d_coord_sel1 --coord_D2=/net/zhao/yh367/Kaiser/t2d/cad_coord_sel1 --rho=0.5 --N1=69033 --N2=86995 --zero_jump_prob=0.05 --num_iter=500 --burn_in=200 --init_PV=0.25,0.25,0.25,0.25 --init_betas=/net/zhao/yh367/Kaiser/t2d/priors/t2d_cad_ld_initial_sel1.pickled.gz --alpha=1,1,1,1 --local_ld_prefix=/net/zhao/yh367/Kaiser/t2d/t2d_cad_sel1 --ld_radius=115 --hfile=/net/zhao/yh367/Kaiser/t2d/priors/t2d_cad_ld_sel1_file.txt --out=/net/zhao/yh367/Kaiser/t2d/output/pleio_non_inf_ld/0.5/t2d_cad_sel1\
python PleioPred_bi_mcmc.py --coord_D1=/net/zhao/yh367/Kaiser/t2d/t2d_coord_sel2 --coord_D2=/net/zhao/yh367/Kaiser/t2d/cad_coord_sel2 --N1=69033 --N2=86995 --zero_jump_prob=0.05 --num_iter=500 --burn_in=200 --init_PV=0.25,0.25,0.25,0.25 --init_betas=/net/zhao/yh367/Kaiser/t2d/priors/t2d_cad_anno_initial_sel2.pickled.gz --alpha=1,1,1,1 --local_ld_prefix=/net/zhao/yh367/Kaiser/t2d/t2d_cad_sel2 --ld_radius=115 --hfile=/net/zhao/yh367/Kaiser/t2d/priors/t2d_cad_anno_sel2_file.txt --out=/net/zhao/yh367/Kaiser/t2d/output/pleio_non_inf/t2d_cad_sel2\
python PleioPred_bi_mcmc.py --coord_D1=/net/zhao/yh367/Kaiser/t2d/t2d_coord_sel2 --coord_D2=/net/zhao/yh367/Kaiser/t2d/cad_coord_sel2 --N1=69033 --N2=86995 --zero_jump_prob=0.05 --num_iter=500 --burn_in=200 --init_PV=0.25,0.25,0.25,0.25 --init_betas=/net/zhao/yh367/Kaiser/t2d/priors/t2d_cad_ld_initial_sel2.pickled.gz --alpha=1,1,1,1 --local_ld_prefix=/net/zhao/yh367/Kaiser/t2d/t2d_cad_sel2 --ld_radius=115 --hfile=/net/zhao/yh367/Kaiser/t2d/priors/t2d_cad_ld_sel2_file.txt --out=/net/zhao/yh367/Kaiser/t2d/output/pleio_non_inf_ld/t2d_cad_sel2\


##### Kaiser T2D Single Prediction #####
pv = c('1.0','0.3','0.1','0.03','0.01','0.003','0.001','0.0003','0.0001','3e-05','1e-05')
c1 = paste0('source ~/.bashrc; cd /net/zhao/yh367/prediction_codes; python LDpred_no_output.py --coord=/net/zhao/yh367/Kaiser/t2d/t2d_coord_sel1 --ld_radius=115 --PS=',pv,' --local_ld_file_prefix=/net/zhao/yh367/Kaiser/t2d/t2d_sel1 --N=69033 --out=/net/zhao/yh367/Kaiser/t2d/output/ldpred_sel1/t2d_p_',pv)
c2 = paste0('source ~/.bashrc; cd /net/zhao/yh367/prediction_codes; python AnnoPred_h2.py --coord=/net/zhao/yh367/Kaiser/t2d/t2d_coord_sel1 --ld_radius=115 --local_ld_file_prefix=/net/zhao/yh367/Kaiser/t2d/t2d_sel1 --ps=',pv,' --PS=/net/zhao/yh367/Kaiser/t2d/priors/h2_sel1_file.txt --N=69033 --out=/net/zhao/yh367/Kaiser/t2d/output/h2_sel1/t2d_p_',pv)
c3 = paste0('source ~/.bashrc; cd /net/zhao/yh367/prediction_codes; python AnnoPred_pT.py --coord=/net/zhao/yh367/Kaiser/t2d/t2d_coord_sel1 --ld_radius=115 --local_ld_file_prefix=/net/zhao/yh367/Kaiser/t2d/t2d_sel1 --ps=',pv,' --PS=/net/zhao/yh367/Kaiser/t2d/priors/pT_sel1_',pv,'_file.txt --N=69033 --out=/net/zhao/yh367/Kaiser/t2d/output/pT_sel1/t2d_p_',pv)

writeLines(c(c1[-1],c2[-1],c3[-1]), '/net/zhao/yh367/Kaiser/SQ_Single_Trait_Prediction/t2d_sel1')
cd /usr/local/cluster/software/installation/SimpleQueue
./sqPBS.py zhao 3 t2d_single1 t2d_single1 /net/zhao/yh367/Kaiser/SQ_Single_Trait_Prediction/t2d_sel1 > /net/zhao/yh367/Kaiser/SQ_Single_Trait_Prediction/t2d_sel1.pbs
cd /net/zhao/yh367/Kaiser/SQ_Single_Trait_Prediction/
qsub t2d_sel1.pbs

dir.create('/net/zhao/yh367/Kaiser/t2d/output/ldpred_sel2'); dir.create('/net/zhao/yh367/Kaiser/t2d/output/h2_sel2'); dir.create('/net/zhao/yh367/Kaiser/t2d/output/pT_sel2')
pv = c('1.0','0.3','0.1','0.03','0.01','0.003','0.001','0.0003','0.0001','3e-05','1e-05')
c1 = paste0('source ~/.bashrc; cd /net/zhao/yh367/prediction_codes; python LDpred_no_output.py --coord=/net/zhao/yh367/Kaiser/t2d/t2d_coord_sel2 --ld_radius=115 --PS=',pv,' --local_ld_file_prefix=/net/zhao/yh367/Kaiser/t2d/t2d_sel2 --N=69033 --out=/net/zhao/yh367/Kaiser/t2d/output/ldpred_sel2/t2d_p_',pv)
c2 = paste0('source ~/.bashrc; cd /net/zhao/yh367/prediction_codes; python AnnoPred_h2.py --coord=/net/zhao/yh367/Kaiser/t2d/t2d_coord_sel2 --ld_radius=115 --local_ld_file_prefix=/net/zhao/yh367/Kaiser/t2d/t2d_sel2 --ps=',pv,' --PS=/net/zhao/yh367/Kaiser/t2d/priors/h2_sel2_file.txt --N=69033 --out=/net/zhao/yh367/Kaiser/t2d/output/h2_sel2/t2d_p_',pv)
c3 = paste0('source ~/.bashrc; cd /net/zhao/yh367/prediction_codes; python AnnoPred_pT.py --coord=/net/zhao/yh367/Kaiser/t2d/t2d_coord_sel2 --ld_radius=115 --local_ld_file_prefix=/net/zhao/yh367/Kaiser/t2d/t2d_sel2 --ps=',pv,' --PS=/net/zhao/yh367/Kaiser/t2d/priors/pT_sel2_',pv,'_file.txt --N=69033 --out=/net/zhao/yh367/Kaiser/t2d/output/pT_sel2/t2d_p_',pv)

writeLines(c(c1,c2,c3), '/net/zhao/yh367/Kaiser/SQ_Single_Trait_Prediction/t2d_sel2')
cd /usr/local/cluster/software/installation/SimpleQueue
./sqPBS.py zhao 3 t2d_single2 t2d_single2 /net/zhao/yh367/Kaiser/SQ_Single_Trait_Prediction/t2d_sel2 > /net/zhao/yh367/Kaiser/SQ_Single_Trait_Prediction/t2d_sel2.pbs
cd /net/zhao/yh367/Kaiser/SQ_Single_Trait_Prediction/
qsub t2d_sel2.pbs


bi_rho_cmds = function(coord1, coord2, rhoV, N1, N2, zj_p, iter, b_in, iPV, iBeta, alpha, localLD, LD_rad, hfile, outdir, dis_prefix){
	paste0('source ~/.bashrc; cd /net/zhao/yh367/prediction_codes; python PleioPred_bi_rho.py --coord_D1=',coord1,' --coord_D2=',coord2,' --rho=',rhoV,' --N1=',N1,' --N2=',N2,' --zero_jump_prob=',zj_p,' --num_iter=',iter,' --burn_in=',b_in,' --init_PV=',iPV,' --init_betas=',iBeta,' --alpha=',alpha,' --local_ld_prefix=',localLD,' --ld_radius=',LD_rad,' --hfile=',hfile,' --out=',outdir, rhoV,'/',dis_prefix)
}

##### Kaiser T2D Double Prediction #####
rhoV = seq(-0.9,0.9,0.1)
for(rv in rhoV){
	dir.create(paste0('/net/zhao/yh367/Kaiser/t2d/output/pleio_non_inf/',rv))
	dir.create(paste0('/net/zhao/yh367/Kaiser/t2d/output/pleio_non_inf_ld/',rv))
}
#anno1 = paste0('source ~/.bashrc; cd /net/zhao/yh367/prediction_codes; python PleioPred_bi_rho.py --coord_D1=/net/zhao/yh367/Kaiser/t2d/t2d_coord_sel1 --coord_D2=/net/zhao/yh367/Kaiser/t2d/cad_coord_sel1 --rho=',rhoV,' --N1=69033 --N2=86995 --zero_jump_prob=',zj_p,' --num_iter=',iter,' --burn_in=',b_in,' --init_PV=0.25,0.25,0.25,0.25 --init_betas=/net/zhao/yh367/Kaiser/t2d/priors/t2d_cad_anno_initial_sel1.pickled.gz --alpha=1,1,1,1 --local_ld_prefix=/net/zhao/yh367/Kaiser/t2d/t2d_cad_sel1 --ld_radius=115 --hfile=/net/zhao/yh367/Kaiser/t2d/priors/t2d_cad_anno_sel1_file.txt --out=/net/zhao/yh367/Kaiser/t2d/output/pleio_non_inf/',rhoV,'/t2d_cad_sel1')
#anno2 = paste0('source ~/.bashrc; cd /net/zhao/yh367/prediction_codes; python PleioPred_bi_rho.py --coord_D1=/net/zhao/yh367/Kaiser/t2d/t2d_coord_sel2 --coord_D2=/net/zhao/yh367/Kaiser/t2d/cad_coord_sel2 --rho=',rhoV,' --N1=69033 --N2=86995 --zero_jump_prob=',zj_p,' --num_iter=',iter,' --burn_in=',b_in,' --init_PV=0.25,0.25,0.25,0.25 --init_betas=/net/zhao/yh367/Kaiser/t2d/priors/t2d_cad_anno_initial_sel2.pickled.gz --alpha=1,1,1,1 --local_ld_prefix=/net/zhao/yh367/Kaiser/t2d/t2d_cad_sel2 --ld_radius=115 --hfile=/net/zhao/yh367/Kaiser/t2d/priors/t2d_cad_anno_sel2_file.txt --out=/net/zhao/yh367/Kaiser/t2d/output/pleio_non_inf/',rhoV,'/t2d_cad_sel2')
#ld1 = paste0('source ~/.bashrc; cd /net/zhao/yh367/prediction_codes; python PleioPred_bi_rho.py --coord_D1=/net/zhao/yh367/Kaiser/t2d/t2d_coord_sel1 --coord_D2=/net/zhao/yh367/Kaiser/t2d/cad_coord_sel1 --rho=',rhoV,' --N1=69033 --N2=86995 --zero_jump_prob=',zj_p,' --num_iter=',iter,' --burn_in=',b_in,' --init_PV=0.25,0.25,0.25,0.25 --init_betas=/net/zhao/yh367/Kaiser/t2d/priors/t2d_cad_ld_initial_sel1.pickled.gz --alpha=1,1,1,1 --local_ld_prefix=/net/zhao/yh367/Kaiser/t2d/t2d_cad_sel1 --ld_radius=115 --hfile=/net/zhao/yh367/Kaiser/t2d/priors/t2d_cad_ld_sel1_file.txt --out=/net/zhao/yh367/Kaiser/t2d/output/pleio_non_inf_ld/',rhoV,'/t2d_cad_sel1')
#ld2 = paste0('source ~/.bashrc; cd /net/zhao/yh367/prediction_codes; python PleioPred_bi_rho.py --coord_D1=/net/zhao/yh367/Kaiser/t2d/t2d_coord_sel2 --coord_D2=/net/zhao/yh367/Kaiser/t2d/cad_coord_sel2 --rho=',rhoV,' --N1=69033 --N2=86995 --zero_jump_prob=',zj_p,' --num_iter=',iter,' --burn_in=',b_in,' --init_PV=0.25,0.25,0.25,0.25 --init_betas=/net/zhao/yh367/Kaiser/t2d/priors/t2d_cad_ld_initial_sel2.pickled.gz --alpha=1,1,1,1 --local_ld_prefix=/net/zhao/yh367/Kaiser/t2d/t2d_cad_sel2 --ld_radius=115 --hfile=/net/zhao/yh367/Kaiser/t2d/priors/t2d_cad_ld_sel2_file.txt --out=/net/zhao/yh367/Kaiser/t2d/output/pleio_non_inf_ld/',rhoV,'/t2d_cad_sel2')
anno1 = bi_rho_cmds(coord1='/net/zhao/yh367/Kaiser/t2d/t2d_coord_sel1', coord2='/net/zhao/yh367/Kaiser/t2d/cad_coord_sel1', rhoV=seq(-0.9,0.9,0.1), N1=69033, N2=86995, zj_p=0.05, iter=250, b_in=100, iPV='0.25,0.25,0.25,0.25', iBeta='/net/zhao/yh367/Kaiser/t2d/priors/t2d_cad_anno_initial_sel1', alpha='1,1,1,1', localLD='/net/zhao/yh367/Kaiser/t2d/t2d_cad_sel1', LD_rad=115, hfile='/net/zhao/yh367/Kaiser/t2d/priors/t2d_cad_anno_sel1_file.txt', outdir='/net/zhao/yh367/Kaiser/t2d/output/pleio_non_inf/', dis_prefix='t2d_cad_sel1')
anno2 = bi_rho_cmds(coord1='/net/zhao/yh367/Kaiser/t2d/t2d_coord_sel2', coord2='/net/zhao/yh367/Kaiser/t2d/cad_coord_sel2', rhoV=seq(-0.9,0.9,0.1), N1=69033, N2=86995, zj_p=0.05, iter=250, b_in=100, iPV='0.25,0.25,0.25,0.25', iBeta='/net/zhao/yh367/Kaiser/t2d/priors/t2d_cad_anno_initial_sel2', alpha='1,1,1,1', localLD='/net/zhao/yh367/Kaiser/t2d/t2d_cad_sel2', LD_rad=115, hfile='/net/zhao/yh367/Kaiser/t2d/priors/t2d_cad_anno_sel2_file.txt', outdir='/net/zhao/yh367/Kaiser/t2d/output/pleio_non_inf/', dis_prefix='t2d_cad_sel2')
ld1 = bi_rho_cmds(coord1='/net/zhao/yh367/Kaiser/t2d/t2d_coord_sel1', coord2='/net/zhao/yh367/Kaiser/t2d/cad_coord_sel1', rhoV=seq(-0.9,0.9,0.1), N1=69033, N2=86995, zj_p=0.05, iter=250, b_in=100, iPV='0.25,0.25,0.25,0.25', iBeta='/net/zhao/yh367/Kaiser/t2d/priors/t2d_cad_ld_initial_sel1', alpha='1,1,1,1', localLD='/net/zhao/yh367/Kaiser/t2d/t2d_cad_sel1', LD_rad=115, hfile='/net/zhao/yh367/Kaiser/t2d/priors/t2d_cad_ld_sel1_file.txt', outdir='/net/zhao/yh367/Kaiser/t2d/output/pleio_non_inf_ld/', dis_prefix='t2d_cad_sel1')
ld2 = bi_rho_cmds(coord1='/net/zhao/yh367/Kaiser/t2d/t2d_coord_sel2', coord2='/net/zhao/yh367/Kaiser/t2d/cad_coord_sel2', rhoV=seq(-0.9,0.9,0.1), N1=69033, N2=86995, zj_p=0.05, iter=250, b_in=100, iPV='0.25,0.25,0.25,0.25', iBeta='/net/zhao/yh367/Kaiser/t2d/priors/t2d_cad_ld_initial_sel2', alpha='1,1,1,1', localLD='/net/zhao/yh367/Kaiser/t2d/t2d_cad_sel2', LD_rad=115, hfile='/net/zhao/yh367/Kaiser/t2d/priors/t2d_cad_ld_sel2_file.txt', outdir='/net/zhao/yh367/Kaiser/t2d/output/pleio_non_inf_ld/', dis_prefix='t2d_cad_sel2')
writeLines(c(anno1,anno2,ld1,ld2), paste0('/net/zhao/yh367/Kaiser/SQ_Double_Trait_Prediction/t2d_sparse_rho0'))
cd /usr/local/cluster/software/installation/SimpleQueue
./sqPBS.py zhao 1 t2d_sparse_rho0 t2d_sparse_rho0 /net/zhao/yh367/Kaiser/SQ_Double_Trait_Prediction/t2d_sparse_rho0 > /net/zhao/yh367/Kaiser/SQ_Double_Trait_Prediction/t2d_sparse_rho0.pbs
cd /net/zhao/yh367/Kaiser/SQ_Double_Trait_Prediction/
qsub t2d_sparse_rho0.pbs -lmem=300gb


##### NUgene T2D Single Prediction #####

##### NUgene T2D Double Prediction #####
rhoV = seq(-0.9,0.9,0.1)
for(rv in rhoV){
	dir.create(paste0('/net/zhao/yh367/Kaiser/t2d_nugene/output/pleio_non_inf/',rv))
	dir.create(paste0('/net/zhao/yh367/Kaiser/t2d_nugene/output/pleio_non_inf_ld/',rv))
}
anno1 = bi_rho_cmds(coord1='/net/zhao/yh367/PleioPred/t2d_cad_fg_qc/input/t2d_coord_sel1', coord2='/net/zhao/yh367/PleioPred/t2d_cad_fg_qc/input/cad_coord_sel1', rhoV=seq(-0.9,0.9,0.1), N1=69033, N2=86995, zj_p=0.05, iter=250, b_in=100, iPV='0.25,0.25,0.25,0.25', iBeta='/net/zhao/yh367/PleioPred/t2d_cad_fg_qc/priors/t2d_cad_anno_initial_sel1', alpha='1,1,1,1', localLD='/net/zhao/yh367/PleioPred/t2d_cad_fg_qc/tmp/t2d_cad_sel1', LD_rad=159, hfile='/net/zhao/yh367/PleioPred/t2d_cad/priors/h2_qc_sel1.txt', outdir='/net/zhao/yh367/Kaiser/t2d_nugene/output/pleio_non_inf/', dis_prefix='t2d_cad_sel1')
anno2 = bi_rho_cmds(coord1='/net/zhao/yh367/PleioPred/t2d_cad_fg_qc/input/t2d_coord_sel2', coord2='/net/zhao/yh367/PleioPred/t2d_cad_fg_qc/input/cad_coord_sel2', rhoV=seq(-0.9,0.9,0.1), N1=69033, N2=86995, zj_p=0.05, iter=250, b_in=100, iPV='0.25,0.25,0.25,0.25', iBeta='/net/zhao/yh367/PleioPred/t2d_cad_fg_qc/priors/t2d_cad_anno_initial_sel2', alpha='1,1,1,1', localLD='/net/zhao/yh367/PleioPred/t2d_cad_fg_qc/tmp/t2d_cad_sel2', LD_rad=159, hfile='/net/zhao/yh367/PleioPred/t2d_cad/priors/h2_qc_sel2.txt', outdir='/net/zhao/yh367/Kaiser/t2d_nugene/output/pleio_non_inf/', dis_prefix='t2d_cad_sel2')
ld1 = bi_rho_cmds(coord1='/net/zhao/yh367/PleioPred/t2d_cad_fg_qc/input/t2d_coord_sel1', coord2='/net/zhao/yh367/PleioPred/t2d_cad_fg_qc/input/cad_coord_sel1', rhoV=seq(-0.9,0.9,0.1), N1=69033, N2=86995, zj_p=0.05, iter=250, b_in=100, iPV='0.25,0.25,0.25,0.25', iBeta='/net/zhao/yh367/PleioPred/t2d_cad_fg_qc/priors/t2d_cad_ld_initial_sel1', alpha='1,1,1,1', localLD='/net/zhao/yh367/PleioPred/t2d_cad_fg_qc/tmp/t2d_cad_sel1', LD_rad=159, hfile='/net/zhao/yh367/PleioPred/t2d_cad/priors/ld_h2_qc_sel1.txt', outdir='/net/zhao/yh367/Kaiser/t2d_nugene/output/pleio_non_inf_ld/', dis_prefix='t2d_cad_sel1')
ld2 = bi_rho_cmds(coord1='/net/zhao/yh367/PleioPred/t2d_cad_fg_qc/input/t2d_coord_sel2', coord2='/net/zhao/yh367/PleioPred/t2d_cad_fg_qc/input/cad_coord_sel2', rhoV=seq(-0.9,0.9,0.1), N1=69033, N2=86995, zj_p=0.05, iter=250, b_in=100, iPV='0.25,0.25,0.25,0.25', iBeta='/net/zhao/yh367/PleioPred/t2d_cad_fg_qc/priors/t2d_cad_ld_initial_sel2', alpha='1,1,1,1', localLD='/net/zhao/yh367/PleioPred/t2d_cad_fg_qc/tmp/t2d_cad_sel2', LD_rad=159, hfile='/net/zhao/yh367/PleioPred/t2d_cad/priors/ld_h2_qc_sel2.txt', outdir='/net/zhao/yh367/Kaiser/t2d_nugene/output/pleio_non_inf_ld/', dis_prefix='t2d_cad_sel2')
writeLines(c(anno1[-15],anno2,ld1,ld2), paste0('/net/zhao/yh367/Kaiser/SQ_Double_Trait_Prediction/t2d_nugene_sparse_rho0'))
cd /usr/local/cluster/software/installation/SimpleQueue
./sqPBS.py general 1 t2d_nugene_sparse_rho0 t2d_nugene_sparse_rho0 /net/zhao/yh367/Kaiser/SQ_Double_Trait_Prediction/t2d_nugene_sparse_rho0 > /net/zhao/yh367/Kaiser/SQ_Double_Trait_Prediction/t2d_nugene_sparse_rho0.pbs
cd /net/zhao/yh367/Kaiser/SQ_Double_Trait_Prediction/
qsub t2d_nugene_sparse_rho0.pbs -lmem=200gb


##### Kaiser CD Single Prediction #####

##### Kaiser CD Double Prediction #####
rhoV = seq(-0.9,0.9,0.1)
for(rv in rhoV){
	dir.create(paste0('/net/zhao/yh367/Kaiser/cd/output/pleio_non_inf/',rv))
	dir.create(paste0('/net/zhao/yh367/Kaiser/cd/output/pleio_non_inf_ld/',rv))
}
anno1 = bi_rho_cmds(coord1='/net/zhao/yh367/PleioPred/cd/input/cd_coord_sel1', coord2='/net/zhao/yh367/PleioPred/cd/input/cel_coord_sel1', rhoV=seq(-0.9,0.9,0.1), N1=16730, N2=15283, zj_p=0.05, iter=250, b_in=100, iPV='0.25,0.25,0.25,0.25', iBeta='/net/zhao/yh367/PleioPred/cd/priors/cd_cel_anno_initial_sel1', alpha='1,1,1,1', localLD='/net/zhao/yh367/PleioPred/cd/tmp/cd_cel_sel1', LD_rad=61, hfile='/net/zhao/yh367/PleioPred/cd/priors/h2_cd_cel_sel1.txt', outdir='/net/zhao/yh367/Kaiser/cd/output/pleio_non_inf/', dis_prefix='cd_cel_sel1')
anno2 = bi_rho_cmds(coord1='/net/zhao/yh367/PleioPred/cd/input/cd_coord_sel2', coord2='/net/zhao/yh367/PleioPred/cd/input/cel_coord_sel2', rhoV=seq(-0.9,0.9,0.1), N1=16730, N2=15283, zj_p=0.05, iter=250, b_in=100, iPV='0.25,0.25,0.25,0.25', iBeta='/net/zhao/yh367/PleioPred/cd/priors/cd_cel_anno_initial_sel2', alpha='1,1,1,1', localLD='/net/zhao/yh367/PleioPred/cd/tmp/cd_cel_sel2', LD_rad=61, hfile='/net/zhao/yh367/PleioPred/cd/priors/h2_cd_cel_sel2.txt', outdir='/net/zhao/yh367/Kaiser/cd/output/pleio_non_inf/', dis_prefix='cd_cel_sel2')
ld1 = bi_rho_cmds(coord1='/net/zhao/yh367/PleioPred/cd/input/cd_coord_sel1', coord2='/net/zhao/yh367/PleioPred/cd/input/cel_coord_sel1', rhoV=seq(-0.9,0.9,0.1), N1=16730, N2=15283, zj_p=0.05, iter=250, b_in=100, iPV='0.25,0.25,0.25,0.25', iBeta='/net/zhao/yh367/PleioPred/cd/priors/cd_cel_ld_initial_sel1', alpha='1,1,1,1', localLD='/net/zhao/yh367/PleioPred/cd/tmp/cd_cel_sel1', LD_rad=61, hfile='/net/zhao/yh367/PleioPred/cd/priors/ld_h2_cd_cel_sel1.txt', outdir='/net/zhao/yh367/Kaiser/cd/output/pleio_non_inf_ld/', dis_prefix='cd_cel_sel1')
ld2 = bi_rho_cmds(coord1='/net/zhao/yh367/PleioPred/cd/input/cd_coord_sel2', coord2='/net/zhao/yh367/PleioPred/cd/input/cel_coord_sel2', rhoV=seq(-0.9,0.9,0.1), N1=16730, N2=15283, zj_p=0.05, iter=250, b_in=100, iPV='0.25,0.25,0.25,0.25', iBeta='/net/zhao/yh367/PleioPred/cd/priors/cd_cel_ld_initial_sel2', alpha='1,1,1,1', localLD='/net/zhao/yh367/PleioPred/cd/tmp/cd_cel_sel2', LD_rad=61, hfile='/net/zhao/yh367/PleioPred/cd/priors/ld_h2_cd_cel_sel2.txt', outdir='/net/zhao/yh367/Kaiser/cd/output/pleio_non_inf_ld/', dis_prefix='cd_cel_sel2')
writeLines(c(anno1,anno2,ld1,ld2), paste0('/net/zhao/yh367/Kaiser/SQ_Double_Trait_Prediction/cd_cel_sparse_rho0'))
cd /usr/local/cluster/software/installation/SimpleQueue
./sqPBS.py general 1 cd_cel_sparse_rho0 cd_cel_sparse_rho0 /net/zhao/yh367/Kaiser/SQ_Double_Trait_Prediction/cd_cel_sparse_rho0 > /net/zhao/yh367/Kaiser/SQ_Double_Trait_Prediction/cd_cel_sparse_rho0.pbs
cd /net/zhao/yh367/Kaiser/SQ_Double_Trait_Prediction/
qsub cd_cel_sparse_rho0.pbs -lmem=200gb

anno1 = bi_rho_cmds(coord1='/net/zhao/yh367/PleioPred/cd/input/cd_coord_sel1', coord2='/net/zhao/yh367/PleioPred/cd/input/uc_coord_sel1', rhoV=seq(-0.9,0.9,0.1), N1=16730, N2=26405, zj_p=0.05, iter=250, b_in=100, iPV='0.25,0.25,0.25,0.25', iBeta='/net/zhao/yh367/PleioPred/cd/priors/cd_uc_anno_initial_sel1', alpha='1,1,1,1', localLD='/net/zhao/yh367/PleioPred/cd/tmp/cd_uc_sel1', LD_rad=61, hfile='/net/zhao/yh367/PleioPred/cd/priors/h2_cd_uc_sel1.txt', outdir='/net/zhao/yh367/Kaiser/cd/output/pleio_non_inf/', dis_prefix='cd_uc_sel1')
anno2 = bi_rho_cmds(coord1='/net/zhao/yh367/PleioPred/cd/input/cd_coord_sel2', coord2='/net/zhao/yh367/PleioPred/cd/input/uc_coord_sel2', rhoV=seq(-0.9,0.9,0.1), N1=16730, N2=26405, zj_p=0.05, iter=250, b_in=100, iPV='0.25,0.25,0.25,0.25', iBeta='/net/zhao/yh367/PleioPred/cd/priors/cd_uc_anno_initial_sel2', alpha='1,1,1,1', localLD='/net/zhao/yh367/PleioPred/cd/tmp/cd_uc_sel2', LD_rad=61, hfile='/net/zhao/yh367/PleioPred/cd/priors/h2_cd_uc_sel2.txt', outdir='/net/zhao/yh367/Kaiser/cd/output/pleio_non_inf/', dis_prefix='cd_uc_sel2')
ld1 = bi_rho_cmds(coord1='/net/zhao/yh367/PleioPred/cd/input/cd_coord_sel1', coord2='/net/zhao/yh367/PleioPred/cd/input/uc_coord_sel1', rhoV=seq(-0.9,0.9,0.1), N1=16730, N2=26405, zj_p=0.05, iter=250, b_in=100, iPV='0.25,0.25,0.25,0.25', iBeta='/net/zhao/yh367/PleioPred/cd/priors/cd_uc_ld_initial_sel1', alpha='1,1,1,1', localLD='/net/zhao/yh367/PleioPred/cd/tmp/cd_uc_sel1', LD_rad=61, hfile='/net/zhao/yh367/PleioPred/cd/priors/ld_h2_cd_uc_sel1.txt', outdir='/net/zhao/yh367/Kaiser/cd/output/pleio_non_inf_ld/', dis_prefix='cd_uc_sel1')
ld2 = bi_rho_cmds(coord1='/net/zhao/yh367/PleioPred/cd/input/cd_coord_sel2', coord2='/net/zhao/yh367/PleioPred/cd/input/uc_coord_sel2', rhoV=seq(-0.9,0.9,0.1), N1=16730, N2=26405, zj_p=0.05, iter=250, b_in=100, iPV='0.25,0.25,0.25,0.25', iBeta='/net/zhao/yh367/PleioPred/cd/priors/cd_uc_ld_initial_sel2', alpha='1,1,1,1', localLD='/net/zhao/yh367/PleioPred/cd/tmp/cd_uc_sel2', LD_rad=61, hfile='/net/zhao/yh367/PleioPred/cd/priors/ld_h2_cd_uc_sel2.txt', outdir='/net/zhao/yh367/Kaiser/cd/output/pleio_non_inf_ld/', dis_prefix='cd_uc_sel2')
writeLines(c(anno1,anno2,ld1,ld2), paste0('/net/zhao/yh367/Kaiser/SQ_Double_Trait_Prediction/cd_uc_sparse_rho0'))
cd /usr/local/cluster/software/installation/SimpleQueue
./sqPBS.py general 1 cd_uc_sparse_rho0 cd_uc_sparse_rho0 /net/zhao/yh367/Kaiser/SQ_Double_Trait_Prediction/cd_uc_sparse_rho0 > /net/zhao/yh367/Kaiser/SQ_Double_Trait_Prediction/cd_uc_sparse_rho0.pbs
cd /net/zhao/yh367/Kaiser/SQ_Double_Trait_Prediction/
qsub cd_uc_sparse_rho0.pbs -lmem=200gb

cd /usr/local/cluster/software/installation/SimpleQueue
./sqPBS.py zhao 1 cd_uc_sparse_rho0 cd_uc_sparse_rho0 /net/zhao/yh367/Kaiser/SQ_Double_Trait_Prediction/cd_uc_sparse_rho0.REMAINING > /net/zhao/yh367/Kaiser/SQ_Double_Trait_Prediction/cd_uc_sparse_rho0.REMAINING.pbs
cd /net/zhao/yh367/Kaiser/SQ_Double_Trait_Prediction/
qsub cd_uc_sparse_rho0.REMAINING.pbs -lmem=200gb


##### Kaiser CEL Single Prediction #####

##### Kaiser CEL Double Prediction #####
rhoV = seq(-0.9,0.9,0.1)
for(rv in rhoV){
	dir.create(paste0('/net/zhao/yh367/Kaiser/cel/output/pleio_non_inf/',rv))
	dir.create(paste0('/net/zhao/yh367/Kaiser/cel/output/pleio_non_inf_ld/',rv))
}
anno1 = bi_rho_cmds(coord1='/net/zhao/yh367/PleioPred/cel/input/cel_coord_sel1', coord2='/net/zhao/yh367/PleioPred/cel/input/cd_coord_sel1', rhoV=seq(-0.9,0.9,0.1), N1=15283, N2=16370, zj_p=0.05, iter=250, b_in=100, iPV='0.25,0.25,0.25,0.25', iBeta='/net/zhao/yh367/PleioPred/cel/priors/cel_cd_anno_initial_sel1', alpha='1,1,1,1', localLD='/net/zhao/yh367/PleioPred/cel/tmp/cel_cd_sel1', LD_rad=73, hfile='/net/zhao/yh367/PleioPred/cel/priors/h2_cel_cd_sel1.txt', outdir='/net/zhao/yh367/Kaiser/cel/output/pleio_non_inf/', dis_prefix='cel_cd_sel1')
anno2 = bi_rho_cmds(coord1='/net/zhao/yh367/PleioPred/cel/input/cel_coord_sel2', coord2='/net/zhao/yh367/PleioPred/cel/input/cd_coord_sel2', rhoV=seq(-0.9,0.9,0.1), N1=15283, N2=16370, zj_p=0.05, iter=250, b_in=100, iPV='0.25,0.25,0.25,0.25', iBeta='/net/zhao/yh367/PleioPred/cel/priors/cel_cd_anno_initial_sel2', alpha='1,1,1,1', localLD='/net/zhao/yh367/PleioPred/cel/tmp/cel_cd_sel2', LD_rad=73, hfile='/net/zhao/yh367/PleioPred/cel/priors/h2_cel_cd_sel2.txt', outdir='/net/zhao/yh367/Kaiser/cel/output/pleio_non_inf/', dis_prefix='cel_cd_sel2')
ld1 = bi_rho_cmds(coord1='/net/zhao/yh367/PleioPred/cel/input/cel_coord_sel1', coord2='/net/zhao/yh367/PleioPred/cel/input/cd_coord_sel1', rhoV=seq(-0.9,0.9,0.1), N1=15283, N2=16370, zj_p=0.05, iter=250, b_in=100, iPV='0.25,0.25,0.25,0.25', iBeta='/net/zhao/yh367/PleioPred/cel/priors/cel_cd_ld_initial_sel1', alpha='1,1,1,1', localLD='/net/zhao/yh367/PleioPred/cel/tmp/cel_cd_sel1', LD_rad=73, hfile='/net/zhao/yh367/PleioPred/cel/priors/ld_h2_cel_cd_sel1.txt', outdir='/net/zhao/yh367/Kaiser/cel/output/pleio_non_inf_ld/', dis_prefix='cel_cd_sel1')
ld2 = bi_rho_cmds(coord1='/net/zhao/yh367/PleioPred/cel/input/cel_coord_sel2', coord2='/net/zhao/yh367/PleioPred/cel/input/cd_coord_sel2', rhoV=seq(-0.9,0.9,0.1), N1=15283, N2=16370, zj_p=0.05, iter=250, b_in=100, iPV='0.25,0.25,0.25,0.25', iBeta='/net/zhao/yh367/PleioPred/cel/priors/cel_cd_ld_initial_sel2', alpha='1,1,1,1', localLD='/net/zhao/yh367/PleioPred/cel/tmp/cel_cd_sel2', LD_rad=73, hfile='/net/zhao/yh367/PleioPred/cel/priors/ld_h2_cel_cd_sel2.txt', outdir='/net/zhao/yh367/Kaiser/cel/output/pleio_non_inf_ld/', dis_prefix='cel_cd_sel2')
writeLines(c(anno1,anno2,ld1,ld2), paste0('/net/zhao/yh367/Kaiser/SQ_Double_Trait_Prediction/cel_cd_sparse_rho0'))
cd /usr/local/cluster/software/installation/SimpleQueue
./sqPBS.py general 1 cel_cd_sparse_rho0 cel_cd_sparse_rho0 /net/zhao/yh367/Kaiser/SQ_Double_Trait_Prediction/cel_cd_sparse_rho0 > /net/zhao/yh367/Kaiser/SQ_Double_Trait_Prediction/cel_cd_sparse_rho0.pbs
cd /net/zhao/yh367/Kaiser/SQ_Double_Trait_Prediction/
qsub cel_cd_sparse_rho0.pbs -lmem=200gb

anno1 = bi_rho_cmds(coord1='/net/zhao/yh367/PleioPred/cel/input/cel_coord_sel1', coord2='/net/zhao/yh367/PleioPred/cel/input/uc_coord_sel1', rhoV=seq(-0.9,0.9,0.1), N1=15283, N2=26405, zj_p=0.05, iter=250, b_in=100, iPV='0.25,0.25,0.25,0.25', iBeta='/net/zhao/yh367/PleioPred/cel/priors/cel_uc_anno_initial_sel1', alpha='1,1,1,1', localLD='/net/zhao/yh367/PleioPred/cel/tmp/cel_uc_sel1', LD_rad=73, hfile='/net/zhao/yh367/PleioPred/cel/priors/h2_cel_uc_sel1.txt', outdir='/net/zhao/yh367/Kaiser/cel/output/pleio_non_inf/', dis_prefix='cel_uc_sel1')
anno2 = bi_rho_cmds(coord1='/net/zhao/yh367/PleioPred/cel/input/cel_coord_sel2', coord2='/net/zhao/yh367/PleioPred/cel/input/uc_coord_sel2', rhoV=seq(-0.9,0.9,0.1), N1=15283, N2=26405, zj_p=0.05, iter=250, b_in=100, iPV='0.25,0.25,0.25,0.25', iBeta='/net/zhao/yh367/PleioPred/cel/priors/cel_uc_anno_initial_sel2', alpha='1,1,1,1', localLD='/net/zhao/yh367/PleioPred/cel/tmp/cel_uc_sel2', LD_rad=73, hfile='/net/zhao/yh367/PleioPred/cel/priors/h2_cel_uc_sel2.txt', outdir='/net/zhao/yh367/Kaiser/cel/output/pleio_non_inf/', dis_prefix='cel_uc_sel2')
ld1 = bi_rho_cmds(coord1='/net/zhao/yh367/PleioPred/cel/input/cel_coord_sel1', coord2='/net/zhao/yh367/PleioPred/cel/input/uc_coord_sel1', rhoV=seq(-0.9,0.9,0.1), N1=15283, N2=26405, zj_p=0.05, iter=250, b_in=100, iPV='0.25,0.25,0.25,0.25', iBeta='/net/zhao/yh367/PleioPred/cel/priors/cel_uc_ld_initial_sel1', alpha='1,1,1,1', localLD='/net/zhao/yh367/PleioPred/cel/tmp/cel_uc_sel1', LD_rad=73, hfile='/net/zhao/yh367/PleioPred/cel/priors/ld_h2_cel_uc_sel1.txt', outdir='/net/zhao/yh367/Kaiser/cel/output/pleio_non_inf_ld/', dis_prefix='cel_uc_sel1')
ld2 = bi_rho_cmds(coord1='/net/zhao/yh367/PleioPred/cel/input/cel_coord_sel2', coord2='/net/zhao/yh367/PleioPred/cel/input/uc_coord_sel2', rhoV=seq(-0.9,0.9,0.1), N1=15283, N2=26405, zj_p=0.05, iter=250, b_in=100, iPV='0.25,0.25,0.25,0.25', iBeta='/net/zhao/yh367/PleioPred/cel/priors/cel_uc_ld_initial_sel2', alpha='1,1,1,1', localLD='/net/zhao/yh367/PleioPred/cel/tmp/cel_uc_sel2', LD_rad=73, hfile='/net/zhao/yh367/PleioPred/cel/priors/ld_h2_cel_uc_sel2.txt', outdir='/net/zhao/yh367/Kaiser/cel/output/pleio_non_inf_ld/', dis_prefix='cel_uc_sel2')
writeLines(c(anno1,anno2,ld1,ld2), paste0('/net/zhao/yh367/Kaiser/SQ_Double_Trait_Prediction/cel_uc_sparse_rho0'))
cd /usr/local/cluster/software/installation/SimpleQueue
./sqPBS.py general 1 cel_uc_sparse_rho0 cel_uc_sparse_rho0 /net/zhao/yh367/Kaiser/SQ_Double_Trait_Prediction/cel_uc_sparse_rho0 > /net/zhao/yh367/Kaiser/SQ_Double_Trait_Prediction/cel_uc_sparse_rho0.pbs
cd /net/zhao/yh367/Kaiser/SQ_Double_Trait_Prediction/
qsub cel_uc_sparse_rho0.pbs -lmem=200gb

cd /usr/local/cluster/software/installation/SimpleQueue
./sqPBS.py general 1 cel_uc_sparse_rho0 cel_uc_sparse_rho0 /net/zhao/yh367/Kaiser/SQ_Double_Trait_Prediction/cel_uc_sparse_rho0.REMAINING > /net/zhao/yh367/Kaiser/SQ_Double_Trait_Prediction/cel_uc_sparse_rho0.REMAINING.pbs
cd /net/zhao/yh367/Kaiser/SQ_Double_Trait_Prediction/
qsub cel_uc_sparse_rho0.REMAINING.pbs -lmem=200gb


get_res = function(out_dir, dis, rhoV){
	rep = length(rhoV)
	anno1 = c(paste0(out_dir,'/',rhoV,'/'))
	anno2 = c(paste0(out_dir,'/',rhoV,'/'))
	D1_anno1_fl_list = paste0(anno1,dis,'_sel1_auc__PleioPred_D1.txt')
	D1_anno2_fl_list = paste0(anno2,dis,'_sel2_auc__PleioPred_D1.txt')
	
	D1_anno1_y_list = paste0(anno1,dis,'_sel1_y__D1.txt') 
	D1_anno2_y_list = paste0(anno2,dis,'_sel2_y__D1.txt')
	D1_anno1_prs_list = paste0(anno1,dis,'_sel1_prs_PleioPred_D1.txt') 
	D1_anno2_prs_list = paste0(anno2,dis,'_sel2_prs_PleioPred_D1.txt')

	D1_anno1_auc = rep(0,length(anno1))
	D1_anno2_auc = rep(0,length(anno2))

	D1_anno1_cor = rep(0,length(anno1))
	D1_anno2_cor = rep(0,length(anno2))

	initial_PV = matrix(0,rep,4)
	avg_PV = matrix(0,rep,4)
	for(i in 1:rep){
	  if(file.exists(D1_anno1_fl_list[i])){fl = readLines(D1_anno1_fl_list[i]); aucL = unlist(strsplit(fl[4],': ')); D1_anno1_auc[i] = as.numeric(aucL[3])}else{D1_anno1_auc[i] = -1}
	  if(file.exists(D1_anno2_fl_list[i])){fl = readLines(D1_anno2_fl_list[i]); aucL = unlist(strsplit(fl[4],': ')); D1_anno2_auc[i] = as.numeric(aucL[3])}else{D1_anno2_auc[i] = -1}
	  initial_PV[i,] = as.numeric(readLines(paste0(anno1[i],dis,'_sel1_Initial_PV.txt')))	
	  avg_PV[i,] = as.numeric(readLines(paste0(anno1[i],dis,'_sel1_Avg_PV.txt')))	
	}
	for(i in 1:rep){
	  if(file.exists(D1_anno1_y_list[i])){y = as.numeric(readLines(D1_anno1_y_list[i])); prs = as.numeric(readLines(D1_anno1_prs_list[i])); D1_anno1_cor[i] = as.numeric(cor(y,prs))}else{D1_anno1_cor[i] = -1}
	  if(file.exists(D1_anno2_y_list[i])){y = as.numeric(readLines(D1_anno2_y_list[i])); prs = as.numeric(readLines(D1_anno2_prs_list[i])); D1_anno2_cor[i] = as.numeric(cor(y,prs))}else{D1_anno2_cor[i] = -1}
	  initial_PV[i,] = as.numeric(readLines(paste0(anno1[i],dis,'_sel2_Initial_PV.txt')))	
	  avg_PV[i,] = as.numeric(readLines(paste0(anno1[i],dis,'_sel2_Avg_PV.txt')))	
	}		
	cat('AUC1: ', D1_anno1_auc, '\n')
	cat('AUC2: ', D1_anno2_auc, '\n')
	cat('COR1: ', D1_anno1_cor, '\n')
	cat('COR2: ', D1_anno2_cor, '\n')
	cat('Initial PV: ', initial_PV, '\n')
	cat('Posterior PV mean: ', avg_PV, '\n')
#	print('Mean AUC')
#	print(cbind((D1_anno1_auc+D1_anno2_auc)/2,(D1_ld1_auc+D1_ld2_auc)/2))
#	print('Mean COR')
#	print(cbind((D1_anno1_cor+D1_anno2_cor)/2,(D1_ld1_cor+D1_ld2_cor)/2))
#	list(Anno1vsLD1= cbind(D1_anno1,D1_ld1), Anno2vsLD2=cbind(D1_anno2,D1_ld2), InitialPV = initial_PV, AvgPV = avg_PV, anno=(D1_anno1+D1_anno2)/2, ld=(D1_ld1+D1_ld2)/2)	
}

get_res('/net/zhao/yh367/Kaiser/cel/output/pleio_non_inf', 'cel_uc', seq(-0.9,0.9,0.1))
get_res('/net/zhao/yh367/Kaiser/cel/output/pleio_non_inf_ld', 'cel_uc', seq(-0.9,0.9,0.1))

get_res('/net/zhao/yh367/Kaiser/cd/output/pleio_non_inf', 'cd_uc', seq(-0.9,0.9,0.1))
get_res('/net/zhao/yh367/Kaiser/cd/output/pleio_non_inf_ld', 'cd_uc', seq(-0.9,0.9,0.1))


