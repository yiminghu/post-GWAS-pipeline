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
from numpy import random
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
from datetime import datetime
#from progressbar import AnimatedMarker, Bar, BouncingBar, Counter, ETA, \
#    FileTransferSpeed, FormatLabel, Percentage, \
#    ProgressBar, ReverseBar, RotatingMarker, \
#    SimpleProgress, Timer
#import post_betas

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



def initializing(beta_hats, h2V, n_indV, reference_ld_mats=None, ld_window_size=100):
    """
    infinitesimal model with snp-specific heritability derived from annotation
    used as the initial values for MCMC of non-infinitesimal model
    """
    num_betas, n_traits = beta_hats.shape
    updated_betas = np.zeros_like(beta_hats)
    m = len(beta_hats)
    for i, wi in enumerate(range(0, num_betas, ld_window_size)):
        start_i = wi
        stop_i = min(num_betas, wi + ld_window_size)
        curr_window_size = stop_i - start_i
        D = reference_ld_mats[i]
        for j in range(n_traits):
            A = n_indV[j]*D + (m / h2V[j]) * sp.eye(curr_window_size)
            A_inv = linalg.pinv(A)
            updated_betas[start_i: stop_i][:,j] = sp.dot(n_indV[j]*A_inv , beta_hats[start_i: stop_i][:,j])  # Adjust the beta_hats
    
    return updated_betas

def inner_iter(beta_hats1, n_indV, lambda1, lambda2, start_betas1=None, ld_radius=100, ld_dict1=None):
    """
    MCMC of non-infinitesimal model
    """
    #Pi = sp.random.dirichlet((alpha,alpha,alpha,alpha),1).flatten()
    m, n_traits = beta_hats1.shape
    
    curr_betas1 = sp.copy(start_betas1)
    # Iterating over effect estimates in sequential order
    iter_order = sp.arange(m)
    for i, snp_i in enumerate(iter_order):
        start_i = max(0, snp_i - ld_radius)
        focal_i = min(ld_radius, snp_i)
        stop_i = min(m, snp_i + ld_radius + 1)
        D1_i = ld_dict1[snp_i]
        local_betas1 = curr_betas1[start_i: stop_i]
        local_betas1[focal_i] = np.zeros(n_traits)
        num1 = beta_hats1[snp_i] - sp.dot(D1_i , local_betas1)
        bi_lasso1 = np.maximum(np.absolute(num1) - lambda1, 0)*np.sign(num1)
        bi_lasso_norm = np.sqrt(np.sum(bi_lasso1**2))
        curr_betas1[snp_i] = bi_lasso1*np.maximum(1 - lambda2/bi_lasso_norm, 0)
    return curr_betas1

def multi_glasso_train(beta1_current, n_indV, tune_idx, Y, X, lambda1, lambda2, num_iter=60):
    print "Starting training with lambda1 = %.5f and lambda2 = %.5f"%(lambda1, lambda2)
    n_tune = len(tune_idx)
    predicted0 = np.zeros(n_tune)
    for chrom_str in chromosomes_list:
        if chrom_str in chr_list:
            predicted0 += sp.dot(beta1_current[chrom_str][:,0], X[chrom_str][:,tune_idx])
    #print "%.2f kg = %.2f lb = %.2f gal = %.2f l" % (var1, var2, var3, var4)
    tune_cor_old = sp.corrcoef(Y[tune_idx], predicted0)[0, 1]
    print "Tuning COR of initial: %.3f"%tune_cor_old
    for k in range(num_iter):
        predicted1 = np.zeros(n_tune)
        for chrom_str in chromosomes_list:
            if chrom_str in chr_list:
                beta1_current[chrom_str] = inner_iter(
                    beta_hats1 = chrom_betas1[chrom_str],
                    n_indV = n_indV,
                    lambda1 = lambda1,
                    lambda2 = lambda2,
                    start_betas1=beta1_current[chrom_str],
                    ld_radius=ld_radius, 
                    ld_dict1=chrom_ld_dict[chrom_str],
                    )            
                predicted1 += sp.dot(beta1_current[chrom_str][:,0], X[chrom_str][:,tune_idx])
        #tune_err_new = np.mean((y1[tune_idx] - predicted)**2)
        tune_cor_new = sp.corrcoef(Y[tune_idx], predicted1)[0, 1]
        if np.isnan(tune_cor_new):
            break
        print "Tuning COR at %.1f step: %.3f"%(k,tune_cor_new)
        if tune_cor_new <= tune_cor_old:
            break
        else:
            tune_cor_old = tune_cor_new
    return beta1_current, tune_cor_new


def multi_glasso(data_file_D1, n_indV, init_betas_prefix, lambda1, lambda2, ld_radius = None, ld_dict=None, num_iter=60):
    df1 = h5py.File(data_file_D1,'r')
    cord_data_g1 = df1['cord_data']

    has_phenotypes1=False
    if 'y' in df1.keys():
        'Validation phenotypes of disease 1 found.'
        y1 = df1['y'][...]  # Phenotype
        num_individs1 = len(y1)
        prs_D1 = sp.zeros(num_individs1)
        has_phenotypes1=True

    ld_scores_dict = ld_dict['ld_scores_dict']
    chrom_ld_dict = ld_dict['chrom_ld_dict']
    chrom_ref_ld_mats = ld_dict['chrom_ref_ld_mats']
    chrom_snps = ld_dict['chrom_snps']
    chrom_snpids = ld_dict['chrom_snpids']
    chrom_betas1 = ld_dict['chrom_betas1']
        
    num_snps1 = 0
    chr_list = list(set(cord_data_g1.keys()))
    n_traits = chrom_betas1[chr_list[0]].shape[1]
    sum_beta2s1 = np.zeros(n_traits)
    for chrom_str in chromosomes_list: 
        if chrom_str in chr_list:
            betas1 = chrom_betas1[chrom_str]
            n_snps1 = len(betas1)
            num_snps1 += n_snps1
            sum_beta2s1 += np.sum(betas1 ** 2, axis = 0)


    L1 = ld_scores_dict['avg_gw_ld_score']
    chi_square_lambda1 = n_indV * sum_beta2s1 / float(num_snps1)
    print 'Genome-wide lambda inflation of D1:', chi_square_lambda1
    print 'Genome-wide mean LD score of D1:', L1
    gw_h2_ld_score_est1 = np.maximum(0.0001, (np.maximum(1, chi_square_lambda1) - 1) / (n_indV * (L1 / num_snps1)))
    print 'Estimated genome-wide heritability of multiple traits are', gw_h2_ld_score_est1
    
    assert chi_square_lambda1.any()>1, 'Something is wrong with the GWAS summary statistics of D1.  Perhaps there were issues parsing of them, or the given GWAS sample size (N) was too small. Either way, lambda (the mean Chi-square statistic) is too small.'

    ########################### using AnnoPred-baseline as initial values ###############################
    init_betas_path = '%s.pickled.gz'%init_betas_prefix
    if not os.path.isfile(init_betas_path): 
        print 'No initial values found, generating ... '
        init_betas = {}
        for chrom_str in chromosomes_list:
            if chrom_str in chr_list:
                pval_derived_betas1 = chrom_betas1[chrom_str]
                h1V = gw_h2_ld_score_est1*len(pval_derived_betas1)/num_snps1
                annopred_betas1 = initializing(
                    pval_derived_betas1, 
                    h2V=h1V, 
                    n_indV=n_indV,
                    reference_ld_mats=chrom_ref_ld_mats[chrom_str], 
                    ld_window_size=2*ld_radius
                    )
                init_betas[chrom_str] = annopred_betas1
        f = gzip.open(init_betas_path, 'wb')
        cPickle.dump(init_betas, f, protocol=2)
        f.close()
        print 'Initial values is now pickled at %s'%init_betas_path
    else:    
        print 'Loading initial values for mcmc from file: %s'%init_betas_path
        f = gzip.open(init_betas_path, 'r')
        init_betas = cPickle.load(f)
        f.close()

    print 'Starting cross validation!'
    ### divide into tuning and testing ###
    n_tune = num_individs1/2
    n_test = num_individs1 - n_tune
    tune_idx = random.choice(num_individs1, n_tune, replace=False)
    test_idx = np.delete(np.arange(num_individs1), tune_idx)
    cv_idx = np.array([tune_idx,test_idx])
    cv_ss = np.array([n_tune,n_test])

    final_results = np.zeros(2)
    cv_coef = {}
    cv_lambda = {}
    lambda1 = np.array([1e-5,3e-5,5e-5,8e-5,1e-4,3e-4,5e-4,8e-4,1e-3,3e-3,5e-3,8e-3,0.01,0.03,0.05,0.08,0.1,0.3,0.5,0.8,1,1.3,1.5,1.8,2,2.5,3,3.5,4,4.5,5,5.5,6,6.5,7,7.5,8,8.5,9,9.5,10,15,20,25,30,35,40,50,60,70,80,90,100])
    lambda2 = np.array([1e-5,3e-5,5e-5,8e-5,1e-4,3e-4,5e-4,8e-4,1e-3,3e-3,5e-3,8e-3,0.01,0.03,0.05,0.08,0.1,0.3,0.5,0.8,1,1.3,1.5,1.8,2,2.5,3,3.5,4,4.5,5,5.5,6,6.5,7,7.5,8,8.5,9,9.5,10,15,20,25,30,35,40,50,60,70,80,90,100])
    l1 = len(lambda1)
    l2 = len(lambda2)
    for fd in range(2):
        tune_cor_matrix = np.zeros((l1, l2))
        for tune1 in range(l1):
            for tune2 in range(l2):
                beta1_current = init_betas.copy()
                beta1, tune_cor = multi_glasso_train(beta1_current, n_indV, cv_idx[fd], y1, chrom_snps, lambda1[tune1], lambda2[tune2]) 
                num_non_zeros = 0.0
                for chrom_str in chromosomes_list:
                    if chrom_str in chr_list:
                        num_non_zeros += np.sum(beta1[chrom_str][:,0]!=0)
                if num_non_zeros > 0:
                    tune_cor_matrix[tune1,tune2] = tune_cor
                    predicted2 = sp.zeros(n_test)
                    for chrom_str in chromosomes_list:
                        if chrom_str in chr_list:
                            predicted2 += sp.dot(beta1[chrom_str][:,0], chrom_snps[chrom_str][:,cv_idx[(fd+1)%2]])
                    test_cor = sp.corrcoef(y1[cv_idx[(fd+1)%2]], predicted2)[0, 1]
                    print "Testing COR with lambda1 = %0.5f and lambda2 = %0.5f is %0.3f"%(lambda1[tune1], lambda2[tune2], test_cor)
                else:
                    break
        selected = np.unravel_index(tune_cor_matrix.argmax(), tune_cor_matrix.shape)
        beta1_current = init_betas.copy()
        beta1, tune_cor = multi_glasso_train(beta1_current, n_indV, cv_idx[fd], y1, chrom_snps, lambda1[selected[0]], lambda2[selected[1]])
        predicted2 = sp.zeros(n_test)
        for chrom_str in chromosomes_list:
            if chrom_str in chr_list:
                predicted2 += sp.dot(beta1[chrom_str][:,0], chrom_snps[chrom_str][:,cv_idx[(fd+1)%2]])
        test_cor = sp.corrcoef(y1[cv_idx[(fd+1)%2]], predicted2)[0, 1]
        print "Testing COR is %0.3f"%test_cor
        cv_lambda['cv'+str(fd)] = np.array([lambda1[selected[0]], lambda2[selected[1]]])
        cv_coef['cv'+str(fd)] = beta1
        final_results[fd] = test_cor

    avg_cor = np.mean(final_results)
    return avg_cor, cv_lambda, cv_coef

def main(p_dict):
    local_ld_dict_file = '%s_ldradius%d.pickled.gz'%(p_dict['local_ld_prefix'], p_dict['ld_radius'])
    if not os.path.isfile(local_ld_dict_file):
        df1 = h5py.File(p_dict['multi_coord'])
        chrom_ld_scores_dict1 = {}
        chrom_ld_dict1 = {}
        chrom_ref_ld_mats1 = {}
        ld_score_sum1 = 0
        num_snps1 = 0
        chrom_snps1 = {}
        chrom_betas1 = {}
        chrom_snpids = {}

        print 'Calculating LD information w. radius %d'% p_dict['ld_radius']

        cord_data_g1 = df1['cord_data']
        # find overlap of chrom list
        chr_list = list(set(cord_data_g1.keys()))
        for chrom_str in chr_list:
            print 'Working on %s'%chrom_str
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
            raw_snps1 = raw_snps1[ok_snps_filter1]
            snp_means1 = snp_means1[ok_snps_filter1]
            snp_stds1 = snp_stds1[ok_snps_filter1]
            betas1 = betas1[ok_snps_filter1]
            n_snps1 = len(raw_snps1)
            snp_means1.shape = (n_snps1,1)   
            snp_stds1.shape = (n_snps1,1)

            chrom_snpids[chrom_str] = sids1
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

        avg_gw_ld_score1 = ld_score_sum1 / float(num_snps1)
        ld_scores_dict1 = {'avg_gw_ld_score': avg_gw_ld_score1, 'chrom_dict':chrom_ld_scores_dict1}    

        print 'Done calculating the LD table and LD score, writing to file:', local_ld_dict_file
        print 'Genome-wide average LD score was:', ld_scores_dict1['avg_gw_ld_score']
        ld_dict = {'ld_scores_dict':ld_scores_dict1, 'chrom_ld_dict':chrom_ld_dict1, 
        'chrom_ref_ld_mats':chrom_ref_ld_mats1, 'chrom_snps':chrom_snps1, 
        'chrom_betas1':chrom_betas1, 'chrom_snpids':chrom_snpids}

        f = gzip.open(local_ld_dict_file, 'wb')
        cPickle.dump(ld_dict, f, protocol=2)
        f.close()
        print 'LD information is now pickled at %s'%local_ld_dict_file
    else:
        print 'Loading LD information from file: %s'%local_ld_dict_file
        f = gzip.open(local_ld_dict_file, 'r')
        ld_dict = cPickle.load(f)
        f.close()
    now = datetime.now()
    cv_res, cv_lambda, cv_coef = multi_glasso(p_dict["multi_coord"], p_dict['N'], p_dict['init_betas'], ld_radius = p_dict['ld_radius'], ld_dict=ld_dict, num_iter=60)
    print datetime.now() - now
    
    coef_file = "%s_coef_ldradius%d.pickled.gz"%(p_dict["out"], p_dict['ld_radius'])
    f = gzip.open(coef_file, 'wb')
    cPickle.dump(cv_coef, f, protocol=2)
    f.close()
    print "Coef. est. are now pickled at %s"%coef_file

    log_file = []
    log_file.append("cv_avg_cor = %0.3f\n"%cv_res)
    log_file.append("cv0_lambda = %0.6f, %0.6f\n"%(cv_lambda["cv0"][0], cv_lambda["cv0"][1]))
    log_file.append("cv1_lambda = %0.6f, %0.6f\n"%(cv_lambda["cv1"][0], cv_lambda["cv1"][1]))
    ff = open("%s_cv.log"%(p_dict["out"]), 'w')
    ff.writelines(log_file)
    ff.close()
    
               
     