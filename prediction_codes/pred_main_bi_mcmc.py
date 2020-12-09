try: 
    import scipy as sp
except Exception:
    print 'Using Numpy instead of Scipy.'
    import numpy as sp
    
from scipy import linalg 
import pdb
import plinkio
from plinkio import plinkfile
import random
import numpy as np
import time
import gzip
import itertools as it
from sklearn import metrics

import getopt
import sys
import traceback
import time
import os
import gzip
import itertools as it
import scipy as sp
import h5py
from scipy import stats
import cPickle
from progressbar import AnimatedMarker, Bar, BouncingBar, Counter, ETA, \
    FileTransferSpeed, FormatLabel, Percentage, \
    ProgressBar, ReverseBar, RotatingMarker, \
    SimpleProgress, Timer
import post_betas

chromosomes_list = ['chrom_%d'%(x) for x in range(1,23)]
chromosomes_list.append('chrom_X')


def pred_accuracy(y_true, y_pred):
    y_true = sp.copy(y_true)
    if len(sp.unique(y_true))==2:
        print 'dichotomous trait, calculating AUC'
        y_min = y_true.min()
        y_max = y_true.max()
        if y_min!= 0 or y_max!=1:
            y_true[y_true==y_min]=0
            y_true[y_true==y_max]=1
        fpr, tpr, thresholds = metrics.roc_curve(y_true, y_pred)
        auc = metrics.auc(fpr, tpr)
        return auc
    else:
        print 'continuous trait, calculating COR'
        cor = sp.corrcoef(y_true,y_pred)[0,1]
        return cor


def get_LDpred_ld_tables(snps, ld_radius=100, ld_window_size=0):
    """
    Calculates LD tables, and the LD score in one go...
    """
    
    ld_dict = {}
    m,n = snps.shape
    print m,n
    ld_scores = sp.ones(m)
    ret_dict = {}
    for snp_i, snp in enumerate(snps):
        # Calculate D
        start_i = max(0, snp_i - ld_radius)
        stop_i = min(m, snp_i + ld_radius + 1)
        X = snps[start_i: stop_i]
        D_i = sp.dot(snp, X.T) / n
        r2s = D_i ** 2
        ld_dict[snp_i] = D_i
        lds_i = sp.sum(r2s - (1-r2s) / (n-2),dtype='float32')
        #lds_i = sp.sum(r2s - (1-r2s)*empirical_null_r2)
        ld_scores[snp_i] =lds_i
    ret_dict['ld_dict']=ld_dict
    ret_dict['ld_scores']=ld_scores
    
    if ld_window_size>0:
        ref_ld_matrices = []
        for i, wi in enumerate(range(0, m, ld_window_size)):
            start_i = wi
            stop_i = min(m, wi + ld_window_size)
            curr_window_size = stop_i - start_i
            X = snps[start_i: stop_i]
            D = sp.dot(X, X.T) / n
            ref_ld_matrices.append(D)
        ret_dict['ref_ld_matrices']=ref_ld_matrices
    return ret_dict



def annopred_inf(beta_hats, pr_sigi, n=1000, reference_ld_mats=None, ld_window_size=100):
    """
    infinitesimal model with snp-specific heritability derived from annotation
    used as the initial values for MCMC of non-infinitesimal model
    """
    num_betas = len(beta_hats)
    updated_betas = sp.empty(num_betas)
    m = len(beta_hats)

    for i, wi in enumerate(range(0, num_betas, ld_window_size)):
        start_i = wi
        stop_i = min(num_betas, wi + ld_window_size)
        curr_window_size = stop_i - start_i
        Li = 1.0/pr_sigi[start_i: stop_i]
        D = reference_ld_mats[i]
        A = (n/(1))*D + sp.diag(Li)
        A_inv = linalg.pinv(A)
        updated_betas[start_i: stop_i] = sp.dot(A_inv / (1.0/n) , beta_hats[start_i: stop_i])  # Adjust the beta_hats

    return updated_betas

"""
res_dict = non_infinitesimal_mcmc(pval_derived_betas, Pi = prf_pi_chri_sorted, Sigi2=prf_sigi2_chri_sorted, sig_12=sig_12, h2=h2_chrom, n=n, ld_radius=ld_radius,
                        num_iter=num_iter, burn_in=burn_in, ld_dict=chrom_ld_dict[chrom_str],
                        start_betas=annopred_inf_chrom_dict[chrom_str], zero_jump_prob=zero_jump_prob)
pleiopred_inf(beta_hats1, beta_hats2, pr_sig1, pr_sig2, rho=0, n1=1000, n2=1000, ref_ld_mats1=None, ref_ld_mats2=None, ld_window_size=100)
pleiopred_genomewide(data_file_D1, data_file_D2, rho, ld_radius = None, ld_dict=None, out_file_prefix=None, n1=None, n2=None, PRF=None)            
"""


def pleiopred_genomewide(data_file_D1, data_file_D2, alpha, Pi, init_betas_prefix, ld_radius = None, ld_dict=None, out_file_prefix=None, n1=None, n2=None, PRF=None, num_iter=60, burn_in=10, zero_jump_prob=0.05, user_h1=None, user_h2=None):
    """
    Calculate LDpred for a genome
    """    
    prf_chr = PRF['chrom']
    prf_sids = PRF['sids']
    h2_D1 = PRF['h2_D1']
    h2_D2 = PRF['h2_D2']

    df1 = h5py.File(data_file_D1,'r')
    df2 = h5py.File(data_file_D2,'r')
    cord_data_g1 = df1['cord_data']
    cord_data_g2 = df2['cord_data']

    has_phenotypes1=False
    if 'y' in df1.keys():
        'Validation phenotypes of disease 1 found.'
        y1 = df1['y'][...]  # Phenotype
        num_individs1 = len(y1)
        prs_D1 = sp.zeros(num_individs1)
        has_phenotypes1=True

    has_phenotypes2=False
    if 'y' in df2.keys():
        'Validation phenotypes of disease 2 found.'
        y2 = df2['y'][...]  # Phenotype
        num_individs2 = len(y2)
        prs_D2 = sp.zeros(num_individs2)
        has_phenotypes2=True

    ld_scores_dict = ld_dict['ld_scores_dict']
    chrom_ld_dict = ld_dict['chrom_ld_dict']
    chrom_ref_ld_mats = ld_dict['chrom_ref_ld_mats']
    chrom_snps = ld_dict['chrom_snps']
    chrom_snpids = ld_dict['chrom_snpids']

    chrom_betas1 = ld_dict['chrom_betas1']
    chrom_betas2 = ld_dict['chrom_betas2']
        
    num_snps1 = 0
    sum_beta2s1 = 0
    num_snps2 = 0
    sum_beta2s2 = 0

    chr_list = list(set(cord_data_g1.keys()) & set(cord_data_g2.keys()))

    for chrom_str in chromosomes_list: 
        if chrom_str in chr_list:
            betas1 = chrom_betas1[chrom_str]
            n_snps1 = len(betas1)
            num_snps1 += n_snps1
            sum_beta2s1 += sp.sum(betas1 ** 2)
            betas2 = chrom_betas2[chrom_str]
            n_snps2 = len(betas2)
            num_snps2 += n_snps2
            sum_beta2s2 += sp.sum(betas2 ** 2)

    if user_h1 is None or user_h2 is None:        
        L1 = ld_scores_dict['avg_gw_ld_score']
        chi_square_lambda1 = sp.mean(n1 * sum_beta2s1 / float(num_snps1))
        print 'Genome-wide lambda inflation of D1:', chi_square_lambda1
        print 'Genome-wide mean LD score of D1:', L1
        gw_h2_ld_score_est1 = max(0.0001, (max(1, chi_square_lambda1) - 1) / (n1 * (L1 / num_snps1)))
        print 'Estimated genome-wide heritability of D1:', gw_h2_ld_score_est1
        
        assert chi_square_lambda1>1, 'Something is wrong with the GWAS summary statistics of D1.  Perhaps there were issues parsing of them, or the given GWAS sample size (N) was too small. Either way, lambda (the mean Chi-square statistic) is too small.  '
    
        L2 = ld_scores_dict['avg_gw_ld_score']
        chi_square_lambda2 = sp.mean(n2 * sum_beta2s2 / float(num_snps2))
        print 'Genome-wide lambda inflation of D2:', chi_square_lambda2
        print 'Genome-wide mean LD score of D2:', L2
        gw_h2_ld_score_est2 = max(0.0001, (max(1, chi_square_lambda2) - 1) / (n2 * (L2 / num_snps2)))
        print 'Estimated genome-wide heritability of D2:', gw_h2_ld_score_est2
        
        assert chi_square_lambda2>1, 'Something is wrong with the GWAS summary statistics of D2.  Perhaps there were issues parsing of them, or the given GWAS sample size (N) was too small. Either way, lambda (the mean Chi-square statistic) is too small.  '
    else:
        gw_h2_ld_score_est1 = user_h1
        gw_h2_ld_score_est2 = user_h2

    h2_new1 = sp.sum(h2_D1)
    sig_12_D1 = (1.0)/n1    
    pr_sig1 = {}

    h2_new2 = sp.sum(h2_D2)
    sig_12_D2 = (1.0)/n2 
    pr_sig2 = {}

    post_betas1 = {}
    post_betas2 = {}

    out1 = []
    out1.append('Estimated Genome-wide heritability: '+str(gw_h2_ld_score_est1)+'\n')
    out1.append('Posterior variance for each snp: '+str(sig_12_D1)+'\n')

    out2 = []
    out2.append('Estimated Genome-wide heritability: '+str(gw_h2_ld_score_est2)+'\n')
    out2.append('Posterior variance for each snp: '+str(sig_12_D2)+'\n')


## main calculation, chr by chr, posterior betas and prs ##
    
    beta1_current = chrom_betas1
    beta2_current = chrom_betas2

    for chrom_str in chromosomes_list:
        if chrom_str in chr_list:
            print 'Preparing annotation-based priors for Chromosome %s'%((chrom_str.split('_'))[1])           

            pval_derived_betas1 = chrom_betas1[chrom_str]
            pval_derived_betas2 = chrom_betas2[chrom_str]
            sids = chrom_snpids[chrom_str]

            n_snps_chrom = len(sids)

            chri = int(chrom_str.split('_')[1])
            prf_sids_chri = prf_sids[prf_chr==chri]
            h2_D1_chri = h2_D1[prf_chr==chri]
            h2_D2_chri = h2_D2[prf_chr==chri]
            if len(prf_sids_chri)==len(sids):
                if sum(prf_sids_chri==sids)==len(prf_sids_chri):
                    pr_sig1[chrom_str] = sp.copy(h2_D1_chri)
                    pr_sig2[chrom_str] = sp.copy(h2_D2_chri)
                else:
                    print 'sorting prior files'
                    pr_sig1[chrom_str] = sp.zeros(len(sids))
                    pr_sig2[chrom_str] = sp.zeros(len(sids))
                    for i, sid in enumerate(sids):
                        pr_sig1[chrom_str][i] = h2_D1_chri[prf_sids_chri==sid]
                        pr_sig2[chrom_str][i] = h2_D2_chri[prf_sids_chri==sid]
            else:
                print 'extracting prior files'
                pr_sig1[chrom_str] = sp.zeros(len(sids))
                pr_sig2[chrom_str] = sp.zeros(len(sids))
                for i, sid in enumerate(sids):
                    pr_sig1[chrom_str][i] = h2_D1_chri[prf_sids_chri==sid]
                    pr_sig2[chrom_str][i] = h2_D2_chri[prf_sids_chri==sid]

            pr_sig1[chrom_str] = gw_h2_ld_score_est1*pr_sig1[chrom_str]/h2_new1
            pr_sig2[chrom_str] = gw_h2_ld_score_est2*pr_sig2[chrom_str]/h2_new2

    ########################### using AnnoPred-baseline as initial values ###############################
    init_betas_path = '%s.pickled.gz'%init_betas_prefix
    if not os.path.isfile(init_betas_path): 
        print 'No initial values for mcmc found, generating ... '
        anno_post1 = {}
        anno_post2 = {}
        for chrom_str in chromosomes_list:
            if chrom_str in chr_list:
                pval_derived_betas1 = chrom_betas1[chrom_str]
                pval_derived_betas2 = chrom_betas2[chrom_str]
                annopred_betas1 = annopred_inf(
                    pval_derived_betas1, 
                    pr_sigi=pr_sig1[chrom_str], 
                    reference_ld_mats=chrom_ref_ld_mats[chrom_str], 
                    n=n1, 
                    ld_window_size=2*ld_radius
                    )
                annopred_betas2 = annopred_inf(
                    pval_derived_betas2, 
                    pr_sigi=pr_sig2[chrom_str], 
                    reference_ld_mats=chrom_ref_ld_mats[chrom_str], 
                    n=n2, 
                    ld_window_size=2*ld_radius
                    )
                anno_post1[chrom_str] = annopred_betas1
                anno_post2[chrom_str] = annopred_betas2
        init_betas = {'anno_post1':anno_post1, 'anno_post2':anno_post2}
        f = gzip.open(init_betas_path, 'wb')
        cPickle.dump(init_betas, f, protocol=2)
        f.close()
        print 'LD information is now pickled at %s'%init_betas_path
    else:    
        print 'Loading initial values for mcmc from file: %s'%init_betas_path
        f = gzip.open(init_betas_path, 'r')
        init_betas = cPickle.load(f)
        f.close()
    #### initial values ####
    print 'Preparing initial values for MCMC'
    beta1_current = init_betas['anno_post1']
    beta2_current = init_betas['anno_post2']
    avg_betas1 = {}
    avg_betas2 = {}
    avg_PV = sp.zeros(4)
    for chrom_str in chromosomes_list:
        if chrom_str in chr_list:
            avg_betas1[chrom_str] = sp.zeros(len(chrom_betas1[chrom_str]))
            avg_betas2[chrom_str] = sp.zeros(len(chrom_betas2[chrom_str]))
           
#    Pi = sp.random.dirichlet((alpha,alpha,alpha,alpha),1).flatten()
    print 'Initial PV: ('+str(Pi[0])+', '+str(Pi[1])+', '+str(Pi[2])+', '+str(Pi[3])+')'
    sp.savetxt('%s_Initial_PV'%(out_file_prefix)+'.txt',Pi)
    pb = 0
    pbar = ProgressBar(widgets=[Percentage(), ' ', Bar()," ", Timer()], maxval=num_iter*22).start()
    for k in range(num_iter):  #Big iteration
        A1 = 0
        A2 = 0
        A3 = 0
        A4 = 0    
        for chrom_str in chromosomes_list:
            if chrom_str in chr_list:
                posterior_betas = post_betas.bi_mcmc_all_chr(
                    chrom_betas1[chrom_str],
                    chrom_betas2[chrom_str], 
                    Pi=Pi,
                    pr_sig1=pr_sig1[chrom_str],
                    pr_sig2=pr_sig2[chrom_str], 
                    start_betas1=beta1_current[chrom_str],
                    start_betas2=beta2_current[chrom_str],
                    h2_D1=gw_h2_ld_score_est1 * (n_snps_chrom / float(num_snps1)),
                    n1=n1,
                    h2_D2=gw_h2_ld_score_est2 * (n_snps_chrom / float(num_snps2)),
                    n2=n2,
                    ld_radius=ld_radius, 
                    zj_p = zero_jump_prob,
                    ld_dict1=chrom_ld_dict[chrom_str],
                    ld_dict2=chrom_ld_dict[chrom_str]
                    )
                A1 += posterior_betas['A1']
                A2 += posterior_betas['A2']
                A3 += posterior_betas['A3']
                A4 += posterior_betas['A4']
                beta1_current[chrom_str] = posterior_betas['proposed_betas1']
                beta2_current[chrom_str] = posterior_betas['proposed_betas2']
                if k >= burn_in:
                    avg_betas1[chrom_str] += posterior_betas['curr_post_means1'] #Averaging over the posterior means instead of samples.
                    avg_betas2[chrom_str] += posterior_betas['curr_post_means2']
                pb = pb + 1
                pbar.update(pb)
        Pi = sp.random.dirichlet((alpha[0]+A1,alpha[1]+A2,alpha[2]+A3,alpha[3]+A4),1).flatten()
        if k >= burn_in:
            avg_PV += Pi
    pbar.finish()
    
## prs and auc ##
    avg_PV = avg_PV/float(num_iter-burn_in)
    print 'Posterior PV: ('+str(avg_PV[0])+', '+str(avg_PV[1])+', '+str(avg_PV[2])+', '+str(avg_PV[3])+')'
    sp.savetxt('%s_Avg_PV'%(out_file_prefix)+'.txt',avg_PV)

    for chrom_str in chromosomes_list:
        if chrom_str in chr_list:
            avg_betas1[chrom_str] = avg_betas1[chrom_str]/float(num_iter-burn_in)
            avg_betas2[chrom_str] = avg_betas2[chrom_str]/float(num_iter-burn_in)
            if has_phenotypes1:
                prs_chr_D1 = sp.dot(avg_betas1[chrom_str], chrom_snps[chrom_str])
                prs_D1 += prs_chr_D1
            if has_phenotypes2:
                prs_chr_D2 = sp.dot(avg_betas2[chrom_str], chrom_snps[chrom_str])
                prs_D2 += prs_chr_D2

############ PleioPred results #############
    corr_inf1 = sp.corrcoef(y1, prs_D1)[0, 1]
    r2_inf1 = corr_inf1 ** 2
    #results_dict[p_str]['r2_pd']=r2_inf
    print 'D1: the R2 prediction accuracy (observed scale) of PleioPred was: %0.4f (%0.6f)' % (r2_inf1, ((1-r2_inf1)**2)/num_individs1)
    out1.append('D1: the R2 prediction accuracy (observed scale) of PleioPred was: '+str(r2_inf1)+' ('+str(((1-r2_inf1)**2)/num_individs1)+')\n')
    
    if corr_inf1<0:
        prs_D1 = -1* prs_D1
    auc1 = pred_accuracy(y1,prs_D1)
    print 'D1: PleioPred AUC for the whole genome was: %0.4f'%auc1
    out1.append('D1: PleioPred AUC for the whole genome was: '+str(auc1)+'\n')
    out1.append('D1: PleioPred COR for the whole genome was: '+str(corr_inf1)+'\n')

    sp.savetxt('%s_y_'%(out_file_prefix)+'_D1.txt',y1)
    sp.savetxt('%s_prs'%(out_file_prefix)+'_PleioPred_D1.txt',prs_D1)

    #Now calibration                               
    ff_inf = open('%s_auc_'%(out_file_prefix)+'_PleioPred_D1.txt',"w")
    ff_inf.writelines(out1)
    ff_inf.close()

    corr_inf2 = sp.corrcoef(y2, prs_D2)[0, 1]
    r2_inf2 = corr_inf2 ** 2
    #results_dict[p_str]['r2_pd']=r2_inf
    print 'D2: the R2 prediction accuracy (observed scale) of PleioPred was: %0.4f (%0.6f)' % (r2_inf2, ((1-r2_inf2)**2)/num_individs2)
    out2.append('D2: the R2 prediction accuracy (observed scale) of PleioPred was: '+str(r2_inf2)+' ('+str(((1-r2_inf2)**2)/num_individs2)+')\n')
    
    if corr_inf2<0:
        prs_D2 = -1* prs_D2
    auc2 = pred_accuracy(y2,prs_D2)
    print 'D2: PleioPred AUC for the whole genome was: %0.4f'%auc2
    out2.append('D2: PleioPred AUC for the whole genome was: '+str(auc2)+'\n')
    out2.append('D2: PleioPred COR for the whole genome was: '+str(corr_inf2)+'\n')

    sp.savetxt('%s_y_'%(out_file_prefix)+'_D2.txt',y2)
    sp.savetxt('%s_prs'%(out_file_prefix)+'_PleioPred_D2.txt',prs_D2)

    #Now calibration                               
    ff_inf = open('%s_auc_'%(out_file_prefix)+'_PleioPred_D2.txt',"w")
    ff_inf.writelines(out2)
    ff_inf.close()

    f = gzip.open('%s_betas'%(out_file_prefix)+'_PleioPred_D1.pickled.gz', 'wb')
    cPickle.dump(avg_betas1, f, protocol=2)
    f.close()

    f = gzip.open('%s_betas'%(out_file_prefix)+'_PleioPred_D2.pickled.gz', 'wb')
    cPickle.dump(avg_betas2, f, protocol=2)
    f.close()


"""
p_dict = {'coord_D1':None, 'coord_D2':None, 'ld_radius':None, 'local_ld_prefix':None, 'hfile':None, 'PV':, 'gid':, 'out':None, 'N1':None, 'N2':None}
"""
def main(p_dict):
    local_ld_dict_file = '%s_ldradius%d.pickled.gz'%(p_dict['local_ld_prefix'], p_dict['ld_radius'])
    if not os.path.isfile(local_ld_dict_file):
        df1 = h5py.File(p_dict['coord_D1'])
        df2 = h5py.File(p_dict['coord_D2'])   

        chrom_ld_scores_dict1 = {}
        chrom_ld_dict1 = {}
        chrom_ref_ld_mats1 = {}
        ld_score_sum1 = 0
        num_snps1 = 0
        chrom_snps1 = {}
        chrom_betas1 = {}
        chrom_snpids = {}
        chrom_betas2 = {}

        print 'Calculating LD information w. radius %d'% p_dict['ld_radius']

        cord_data_g1 = df1['cord_data']
        cord_data_g2 = df2['cord_data']
        # find overlap of chrom list
        chr_list = list(set(cord_data_g1.keys()) & set(cord_data_g2.keys()))
        for chrom_str in chr_list:
            print 'Working on %s'%chrom_str
            print 'Sorting disease 1'
            g1 = cord_data_g1[chrom_str]
            if 'raw_snps_ref' in g1.keys():
                raw_snps1 = g1['raw_snps_ref'][...]
                snp_stds1 = g1['snp_stds_ref'][...]
                snp_means1 = g1['snp_means_ref'][...]
                betas1 = g1['betas'][...]
            #Filter monomorphic SNPs
            ok_snps_filter1 = snp_stds1>0
            ok_snps_filter1 = ok_snps_filter1.flatten()
            sids1 = g1['sids'][...]
            sids1 = sids1[ok_snps_filter1]

            print 'Sorting disease 2'
            g2 = cord_data_g2[chrom_str]
            if 'raw_snps_ref' in g2.keys():
                raw_snps2 = g2['raw_snps_ref'][...]
                snp_stds2 = g2['snp_stds_ref'][...]
                snp_means2 = g2['snp_means_ref'][...]
                betas2 = g2['betas'][...]
            #Filter monomorphic SNPs
            ok_snps_filter2 = snp_stds2>0
            ok_snps_filter2 = ok_snps_filter2.flatten()
            sids2 = g2['sids'][...]
            sids2 = sids2[ok_snps_filter2]

            print 'Extracting SNPs shared by both disease 1 and 2'
            ind1 = np.in1d(sids1,sids2)
            ind2 = np.in1d(sids2,sids1)
            sids_shared1 = sids1[ind1]
            sids_shared2 = sids2[ind2]
            raw_snps1 = raw_snps1[ok_snps_filter1][ind1]
            snp_means1 = snp_means1[ok_snps_filter1][ind1]
            snp_stds1 = snp_stds1[ok_snps_filter1][ind1]
            betas1 = betas1[ok_snps_filter1][ind1]
            betas2 = betas2[ok_snps_filter2][ind2]
            n_snps1 = len(raw_snps1)
            snp_means1.shape = (n_snps1,1)   
            snp_stds1.shape = (n_snps1,1)
            ### check order ###
            if sum(sids_shared1==sids_shared2)==len(sids_shared2):
                print 'Good!'
            else:
                print 'Shit happens, sorting sids1 and sids2'
                O1 = np.argsort(sids_shared1)
                O2 = np.argsort(sids_shared2)
                O3 = np.argsort(O2)
                sids_shared1 = sids_shared1[O1][O3]
                if sum(sids_shared1==sids_shared2)==len(sids_shared2):
                    raw_snps1 = raw_snps1[O1][O3]
                    snp_means1 = snp_means1[O1][O3]
                    snp_stds1 = snp_stds1[O1][O3]
                    betas1 = betas1[O1][O3]
                else:
                    print 'Stop! Problems with sorting!'

            # Normalize SNPs..
            chrom_snpids[chrom_str] = sids_shared1

            snps1 = sp.array((raw_snps1 - snp_means1)/snp_stds1,dtype='float32')
            assert snps1.shape==raw_snps1.shape, 'Array Shape mismatch'
            chrom_snps1[chrom_str] = snps1
            ret_dict1 = get_LDpred_ld_tables(snps1, ld_radius=p_dict['ld_radius'], ld_window_size=2*p_dict['ld_radius'])
            chrom_ld_dict1[chrom_str] = ret_dict1['ld_dict']
            chrom_ref_ld_mats1[chrom_str] = ret_dict1['ref_ld_matrices']
            ld_scores1 = ret_dict1['ld_scores']
            chrom_ld_scores_dict1[chrom_str] = {'ld_scores':ld_scores1, 'avg_ld_score':sp.mean(ld_scores1)}
            ld_score_sum1 += sp.sum(ld_scores1)
            num_snps1 += n_snps1

            chrom_betas1[chrom_str] = betas1
            chrom_betas2[chrom_str] = betas2


        avg_gw_ld_score1 = ld_score_sum1 / float(num_snps1)
        ld_scores_dict1 = {'avg_gw_ld_score': avg_gw_ld_score1, 'chrom_dict':chrom_ld_scores_dict1}    

        print 'Done calculating the LD table and LD score, writing to file:', local_ld_dict_file
        print 'Genome-wide average LD score was:', ld_scores_dict1['avg_gw_ld_score']
        ld_dict = {'ld_scores_dict':ld_scores_dict1, 'chrom_ld_dict':chrom_ld_dict1, 
        'chrom_ref_ld_mats':chrom_ref_ld_mats1, 'chrom_snps':chrom_snps1, 
        'chrom_betas1':chrom_betas1, 'chrom_betas2':chrom_betas2, 
        'chrom_snpids':chrom_snpids}

        f = gzip.open(local_ld_dict_file, 'wb')
        cPickle.dump(ld_dict, f, protocol=2)
        f.close()
        print 'LD information is now pickled at %s'%local_ld_dict_file
    else:
        print 'Loading LD information from file: %s'%local_ld_dict_file
        f = gzip.open(local_ld_dict_file, 'r')
        ld_dict = cPickle.load(f)
        f.close()
    
##################### using hfile as prior #######################
    print 'Starting calculation using h2 files as priors'
    print 'Loading prior information from file: %s'%p_dict['hfile']
    with open(p_dict['hfile']) as f:
        data = f.readlines()
    prf_chr = sp.empty(len(data),dtype='int8')
    prf_sids = []
    prf_h2_D1 = sp.zeros(len(data))
    prf_h2_D2 = sp.zeros(len(data))
    for i,line in enumerate(data):
        li = line.split()
        prf_chr[i] = int(li[0])
        prf_sids.append(li[1]) 
        prf_h2_D1[i] = float(li[2])
        prf_h2_D2[i] = float(li[3])  
    prf_sids = sp.array(prf_sids,dtype='str')
    prf = {}
    prf['chrom'] = prf_chr
    prf['sids'] = prf_sids
    prf['h2_D1'] = prf_h2_D1
    prf['h2_D2'] = prf_h2_D2
    H2_D1 = sp.sum(prf_h2_D1)
    H2_D2 = sp.sum(prf_h2_D2)

    pleiopred_genomewide(p_dict['coord_D1'], p_dict['coord_D2'], p_dict['alpha'], p_dict['init_PV'], init_betas_prefix=p_dict['init_betas'], out_file_prefix=p_dict['out'], ld_radius=p_dict['ld_radius'], ld_dict = ld_dict, n1=p_dict['N1'], n2=p_dict['N2'], PRF = prf, num_iter=p_dict['num_iter'], burn_in=p_dict['burn_in'], zero_jump_prob=p_dict['zero_jump_prob'], user_h1=p_dict['user_h1'], user_h2=p_dict['user_h2'])

           
 