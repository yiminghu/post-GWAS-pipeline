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
from scipy import linalg
from progressbar import AnimatedMarker, Bar, BouncingBar, Counter, ETA, \
    FileTransferSpeed, FormatLabel, Percentage, \
    ProgressBar, ReverseBar, RotatingMarker, \
    SimpleProgress, Timer


chromosomes_list = ['chrom_%d'%(x) for x in range(1,23)]
chromosomes_list.append('chrom_X')

p_dict1={'coord_D1':'/net/zhao/yh367/PleioPred/t2d_cad/input/T2D_coord_sel1', 'coord_D2':'/net/zhao/yh367/PleioPred/t2d_cad/input/CAD_coord_sel1', 'N1':69033, 'N2':86995,  
'local_ld_prefix':'/net/zhao/yh367/PleioPred/t2d_cad/tmp/t2d_cad_sel1','ld_radius':159, 'hfile':'/net/zhao/yh367/PleioPred/t2d_cad/priors/h2_sel1.txt', 'user_h1':None, 'user_h2':None, 'out':'/net/zhao/yh367/PleioPred/t2d_cad/priors/t2d_cad_anno_initial_sel1.pickled.gz'} 
bi_get_initial(p_dict1)

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


def bi_get_initial(p_dict):
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
        'chrom_betas3':chrom_betas3, 'chrom_snpids':chrom_snpids}

        f = gzip.open(local_ld_dict_file, 'wb')
        cPickle.dump(ld_dict, f, protocol=2)
        f.close()
        print 'LD information is now pickled.'
    else:
        print 'Loading LD information from file: %s'%local_ld_dict_file
        f = gzip.open(local_ld_dict_file, 'r')
        ld_dict = cPickle.load(f)
        f.close()

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
        #prf_pi[i] = p_dict['PS'][0]         
        prf_h2_D1[i] = float(li[2])
        prf_h2_D2[i] = float(li[3])  
    prf_sids = sp.array(prf_sids,dtype='str')
    prf = {}
    prf['chrom'] = prf_chr
    prf['sids'] = prf_sids
    prf['h2_D1'] = prf_h2_D1
    prf['h2_D2'] = prf_h2_D2
    
    data_file_D1=p_dict['coord_D1']
    data_file_D2=p_dict['coord_D2']
    out_file_prefix=p_dict['out']
    ld_radius=p_dict['ld_radius']
    ld_dict = ld_dict
    n1=p_dict['N1']
    n2=p_dict['N2']
    PRF = prf
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

    ld_scores_dict1 = ld_dict1['ld_scores_dict']
    chrom_ld_dict1 = ld_dict1['chrom_ld_dict']
    chrom_ref_ld_mats1 = ld_dict1['chrom_ref_ld_mats']
    chrom_snps1 = ld_dict1['chrom_snps']
    chrom_betas1 = ld_dict1['chrom_betas']
    chrom_snpids = ld_dict1['chrom_snpids']
    ld_scores_dict2 = ld_dict2['ld_scores_dict']
    chrom_ld_dict2 = ld_dict2['chrom_ld_dict']
    chrom_ref_ld_mats2 = ld_dict2['chrom_ref_ld_mats']
    chrom_snps2 = ld_dict2['chrom_snps']
    chrom_betas2 = ld_dict2['chrom_betas']
        
    #results_dict = {}
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
        
    L1 = ld_scores_dict['avg_gw_ld_score']
    chi_square_lambda1 = sp.mean(n1 * sum_beta2s1 / float(num_snps1))
    print 'Genome-wide lambda inflation of D1:', chi_square_lambda1
    print 'Genome-wide mean LD score of D1:', L1
    gw_h2_ld_score_est1 = max(0.0001, (max(1, chi_square_lambda1) - 1) / (n1 * (L1 / num_snps1)))
    print 'Estimated genome-wide heritability of D1:', gw_h2_ld_score_est1
    L2 = ld_scores_dict['avg_gw_ld_score']
    chi_square_lambda2 = sp.mean(n2 * sum_beta2s2 / float(num_snps2))
    print 'Genome-wide lambda inflation of D2:', chi_square_lambda2
    print 'Genome-wide mean LD score of D2:', L2
    gw_h2_ld_score_est2 = max(0.0001, (max(1, chi_square_lambda2) - 1) / (n2 * (L2 / num_snps2)))
    print 'Estimated genome-wide heritability of D2:', gw_h2_ld_score_est2
    h2_new1 = sp.sum(h2_D1)
    sig_12_D1 = (1.0)/n1    
    pr_sig1 = {}
    h2_new2 = sp.sum(h2_D2)
    sig_12_D2 = (1.0)/n2 
    pr_sig2 = {}
    anno_post1 = {}
    anno_post2 = {}
    post_betas1 = {}
    post_betas2 = {}
    ld_post1 = {}
    ld_post2 = {}
    ## main calculation, chr by chr, posterior betas and prs ##
    for chrom_str in chromosomes_list:
        if chrom_str in chr_list:
            print 'Calculating scores for Chromosome %s'%((chrom_str.split('_'))[1])           
            pval_derived_betas1 = chrom_betas1[chrom_str]
            pval_derived_betas2 = chrom_betas2[chrom_str]
            snps1 = chrom_snps1[chrom_str]
            snps2 = chrom_snps2[chrom_str]
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
            annopred_betas1 = annopred_inf(
                pval_derived_betas1, 
                pr_sigi=pr_sig1[chrom_str], 
                reference_ld_mats=chrom_ref_ld_mats1[chrom_str], 
                n=n1, 
                ld_window_size=2*ld_radius
                )
            annopred_betas2 = annopred_inf(
                pval_derived_betas2, 
                pr_sigi=pr_sig2[chrom_str], 
                reference_ld_mats=chrom_ref_ld_mats2[chrom_str], 
                n=n2, 
                ld_window_size=2*ld_radius
                )
            anno_post1[chrom_str] = annopred_betas1
            anno_post2[chrom_str] = annopred_betas2
    anno_post = {'anno_post1':anno_post1, 'anno_post2':anno_post2}

    f = gzip.open(out_file_prefix, 'wb')
    cPickle.dump(anno_post, f, protocol=2)
    f.close()

#f = gzip.open('/net/zhao/yh367/PleioPred/cel_cd_annopred_sel1_initial.pickled.gz', 'wb')
#cPickle.dump(anno_post, f, protocol=2)
#f.close()

#f = gzip.open('/net/zhao/yh367/PleioPred/cel_cd_ldpred_sel2_initial.pickled.gz', 'wb')
#cPickle.dump(anno_post, f, protocol=2)
#f.close()

#f = gzip.open('/net/zhao/yh367/PleioPred/cel_cd_annopred_sel2_initial.pickled.gz', 'wb')
#cPickle.dump(anno_post, f, protocol=2)
#f.close()#


#p_dict={'coord_D1':'/net/zhao/yh367/PleioPred/input2/cad_coord_sel1', 'coord_D2':'/net/zhao/yh367/PleioPred/input2/t2d_coord_sel1', 'N1':86995, 'N2':69033, 'PV':[0.1,0.1,0.2,0.6], 
#'local_ld_prefix':'/net/zhao/yh367/PleioPred/tmp/cad_t2d_sel1','ld_radius':167, 'hfile':'/net/zhao/yh367/PleioPred/priors/shit_sel1.txt', 'out':'/net/zhao/yh367/PleioPred/output/ldpred_sel1'} 
#
#p_dict={'coord_D1':'/net/zhao/yh367/PleioPred/input2/cad_coord_sel2', 'coord_D2':'/net/zhao/yh367/PleioPred/input2/t2d_coord_sel2', 'N1':86995, 'N2':69033, 'PV':[0.1,0.1,0.2,0.6], 
#'local_ld_prefix':'/net/zhao/yh367/PleioPred/tmp/cad_t2d_sel2','ld_radius':167, 'hfile':'/net/zhao/yh367/PleioPred/priors/shit_sel2.txt', 'out':'/net/zhao/yh367/PleioPred/output/ldpred_sel2'} 
#
#p_dict={'coord_D1':'/project/fas/zhao/yh367/CD_CEL/coord/cel_sel1', 'coord_D2':'/project/fas/zhao/yh367/CD_CEL/coord/cd_sel1', 'N1':15283, 'N2':16730, 'PV':[0.1,0.1,0.2,0.6], 
#'local_ld_prefix':'/net/zhao/yh367/PleioPred/cd_cel/tmp/cel_cd_sel1','ld_radius':137, 'hfile':'/net/zhao/yh367/PleioPred/cd_cel/priors/h2_sel1.txt', 'out':'/net/zhao/yh367/PleioPred/cd_cel/output/0.5/cel_cd_sel1'} 
#
#p_dict={'coord_D1':'/project/fas/zhao/yh367/CD_CEL/coord/cel_sel2', 'coord_D2':'/project/fas/zhao/yh367/CD_CEL/coord/cd_sel2', 'N1':15283, 'N2':16730, 'PV':[0.1,0.1,0.2,0.6], 
#'local_ld_prefix':'/net/zhao/yh367/PleioPred/cd_cel/tmp/cel_cd_sel2','ld_radius':137, 'hfile':'/net/zhao/yh367/PleioPred/cd_cel/priors/h2_sel2.txt', 'out':'/net/zhao/yh367/PleioPred/cd_cel/output/0.5/cel_cd_sel2'} 
#
#p_dict={'coord_D1':'/project/fas/zhao/yh367/CD_CEL/coord/cel_sel1', 'coord_D2':'/project/fas/zhao/yh367/CD_CEL/coord/cd_sel1', 'N1':15283, 'N2':16730, 'PV':[0.1,0.1,0.2,0.6], 
#'local_ld_prefix':'/net/zhao/yh367/PleioPred/cd_cel/tmp/cel_cd_sel1','ld_radius':137, 'hfile':'/net/zhao/yh367/PleioPred/cd_cel/priors/shit_sel1.txt', 'out':'/net/zhao/yh367/PleioPred/cd_cel/output/0.5/cel_cd_sel1'} 
#
#p_dict={'coord_D1':'/project/fas/zhao/yh367/CD_CEL/coord/cel_sel2', 'coord_D2':'/project/fas/zhao/yh367/CD_CEL/coord/cd_sel2', 'N1':15283, 'N2':16730, 'PV':[0.1,0.1,0.2,0.6], 
#'local_ld_prefix':'/net/zhao/yh367/PleioPred/cd_cel/tmp/cel_cd_sel2','ld_radius':137, 'hfile':'/net/zhao/yh367/PleioPred/cd_cel/priors/shit_sel2.txt', 'out':'/net/zhao/yh367/PleioPred/cd_cel/output/0.5/cel_cd_sel2'} 

## t2d_cad ##
# anno #
p_dict1={'coord_D1':'/net/zhao/yh367/PleioPred/t2d_cad/input/T2D_coord_sel1', 'coord_D2':'/net/zhao/yh367/PleioPred/t2d_cad/input/CAD_coord_sel1', 'N1':69033, 'N2':86995,  
'local_ld_prefix':'/net/zhao/yh367/PleioPred/t2d_cad/tmp/t2d_cad_sel1','ld_radius':159, 'hfile':'/net/zhao/yh367/PleioPred/t2d_cad/priors/h2_sel1.txt', 'out':'/net/zhao/yh367/PleioPred/t2d_cad/priors/t2d_cad_anno_initial_sel1.pickled.gz'} 
bi_get_initial(p_dict1)

p_dict2={'coord_D1':'/net/zhao/yh367/PleioPred/t2d_cad/input/T2D_coord_sel2', 'coord_D2':'/net/zhao/yh367/PleioPred/t2d_cad/input/CAD_coord_sel2', 'N1':69033, 'N2':86995,  
'local_ld_prefix':'/net/zhao/yh367/PleioPred/t2d_cad/tmp/t2d_cad_sel2','ld_radius':159, 'hfile':'/net/zhao/yh367/PleioPred/t2d_cad/priors/h2_sel2.txt', 'out':'/net/zhao/yh367/PleioPred/t2d_cad/priors/t2d_cad_anno_initial_sel2.pickled.gz'} 
bi_get_initial(p_dict2)

# ld #
options(stringsAsFactors=F)
pp1 = read.table('/net/zhao/yh367/PleioPred/t2d_cad/priors/h2_sel1.txt',header=F)
pp2 = read.table('/net/zhao/yh367/PleioPred/t2d_cad/priors/h2_sel2.txt',header=F)
pp1[,3] = rep(1,nrow(pp1))
pp1[,4] = rep(1,nrow(pp1))
pp2[,3] = rep(1,nrow(pp2))
pp2[,4] = rep(1,nrow(pp2))
write.table(pp1,'/net/zhao/yh367/PleioPred/t2d_cad/priors/ld_h2_sel1.txt',quote=F,row.names=F,col.names=F)
write.table(pp2,'/net/zhao/yh367/PleioPred/t2d_cad/priors/ld_h2_sel2.txt',quote=F,row.names=F,col.names=F)

p_dict1={'coord_D1':'/net/zhao/yh367/PleioPred/t2d_cad/input/T2D_coord_sel1', 'coord_D2':'/net/zhao/yh367/PleioPred/t2d_cad/input/CAD_coord_sel1', 'N1':69033, 'N2':86995,  
'local_ld_prefix':'/net/zhao/yh367/PleioPred/t2d_cad/tmp/t2d_cad_sel1','ld_radius':159, 'hfile':'/net/zhao/yh367/PleioPred/t2d_cad/priors/ld_h2_sel1.txt', 'out':'/net/zhao/yh367/PleioPred/t2d_cad/priors/t2d_cad_ld_initial_sel1.pickled.gz'} 
bi_get_initial(p_dict1)

p_dict2={'coord_D1':'/net/zhao/yh367/PleioPred/t2d_cad/input/T2D_coord_sel2', 'coord_D2':'/net/zhao/yh367/PleioPred/t2d_cad/input/CAD_coord_sel2', 'N1':69033, 'N2':86995,  
'local_ld_prefix':'/net/zhao/yh367/PleioPred/t2d_cad/tmp/t2d_cad_sel2','ld_radius':159, 'hfile':'/net/zhao/yh367/PleioPred/t2d_cad/priors/ld_h2_sel2.txt', 'out':'/net/zhao/yh367/PleioPred/t2d_cad/priors/t2d_cad_ld_initial_sel2.pickled.gz'} 
bi_get_initial(p_dict2)




## t2d_fg ##
# anno #
p_dict1={'coord_D1':'/net/zhao/yh367/PleioPred/t2d_cad_fg/input/t2d_coord_sel1', 'coord_D2':'/net/zhao/yh367/PleioPred/t2d_cad_fg/input/fg_coord_sel1', 'N1':69033, 'N2':46186,  
'local_ld_prefix':'/net/zhao/yh367/PleioPred/t2d_fg/tmp/t2d_fg_sel1','ld_radius':159, 'hfile':'/net/zhao/yh367/PleioPred/t2d_fg/priors/h2_sel1.txt', 'out':'/net/zhao/yh367/PleioPred/t2d_fg/priors/t2d_fg_anno_initial_sel1.pickled.gz'} 
bi_get_initial(p_dict1)

p_dict2={'coord_D1':'/net/zhao/yh367/PleioPred/t2d_cad_fg/input/t2d_coord_sel2', 'coord_D2':'/net/zhao/yh367/PleioPred/t2d_cad_fg/input/fg_coord_sel2', 'N1':69033, 'N2':46186,  
'local_ld_prefix':'/net/zhao/yh367/PleioPred/t2d_fg/tmp/t2d_fg_sel2','ld_radius':159, 'hfile':'/net/zhao/yh367/PleioPred/t2d_fg/priors/h2_sel2.txt', 'out':'/net/zhao/yh367/PleioPred/t2d_fg/priors/t2d_fg_anno_initial_sel2.pickled.gz'} 
bi_get_initial(p_dict2)

# ld #
options(stringsAsFactors=F)
pp1 = read.table('/net/zhao/yh367/PleioPred/t2d_fg/priors/h2_sel1.txt',header=F)
pp2 = read.table('/net/zhao/yh367/PleioPred/t2d_fg/priors/h2_sel2.txt',header=F)
pp1[,3] = rep(1,nrow(pp1))
pp1[,4] = rep(1,nrow(pp1))
pp2[,3] = rep(1,nrow(pp2))
pp2[,4] = rep(1,nrow(pp2))
write.table(pp1,'/net/zhao/yh367/PleioPred/t2d_fg/priors/ld_h2_sel1.txt',quote=F,row.names=F,col.names=F)
write.table(pp2,'/net/zhao/yh367/PleioPred/t2d_fg/priors/ld_h2_sel2.txt',quote=F,row.names=F,col.names=F)

p_dict1={'coord_D1':'/net/zhao/yh367/PleioPred/t2d_cad_fg/input/t2d_coord_sel1', 'coord_D2':'/net/zhao/yh367/PleioPred/t2d_cad_fg/input/fg_coord_sel1', 'N1':69033, 'N2':46186,  
'local_ld_prefix':'/net/zhao/yh367/PleioPred/t2d_fg/tmp/t2d_fg_sel1','ld_radius':159, 'hfile':'/net/zhao/yh367/PleioPred/t2d_fg/priors/ld_h2_sel1.txt', 'out':'/net/zhao/yh367/PleioPred/t2d_fg/priors/t2d_fg_ld_initial_sel1.pickled.gz'} 
bi_get_initial(p_dict1)

p_dict2={'coord_D1':'/net/zhao/yh367/PleioPred/t2d_cad_fg/input/t2d_coord_sel2', 'coord_D2':'/net/zhao/yh367/PleioPred/t2d_cad_fg/input/fg_coord_sel2', 'N1':69033, 'N2':46186,  
'local_ld_prefix':'/net/zhao/yh367/PleioPred/t2d_fg/tmp/t2d_fg_sel2','ld_radius':159, 'hfile':'/net/zhao/yh367/PleioPred/t2d_fg/priors/ld_h2_sel2.txt', 'out':'/net/zhao/yh367/PleioPred/t2d_fg/priors/t2d_fg_ld_initial_sel2.pickled.gz'} 
bi_get_initial(p_dict2)










