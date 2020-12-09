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

def bi_rho_all_chr(beta_hats1, beta_hats2, rho, Pi, pr_sig1, pr_sig2, zj_p, start_betas1=None, start_betas2=None, n1=1000, n2=1000, ld_radius=100, ld_dict1=None, ld_dict2=None, h2_D1=None, h2_D2=None):
    """
    MCMC of non-infinitesimal model
    """
    #Pi = sp.random.dirichlet((alpha,alpha,alpha,alpha),1).flatten()
    m = len(beta_hats1)
    
    curr_betas1 = sp.copy(start_betas1)
    curr_post_means1 = sp.zeros(m)
#    avg_betas1 = sp.zeros(m)

    curr_betas2 = sp.copy(start_betas2)
    curr_post_means2 = sp.zeros(m)

    # Iterating over effect estimates in sequential order
    h2_est1 = max(0.00001,sp.sum(curr_betas1 ** 2))
    h2_est2 = max(0.00001,sp.sum(curr_betas2 ** 2))
    shrink_factor = min(1-zj_p, h2_D1/h2_est1, h2_D2/h2_est2)

    rand_ps = sp.random.random(m)
    iter_order = sp.arange(m)
    for i, snp_i in enumerate(iter_order):
        if pr_sig1[snp_i]==0 or pr_sig2[snp_i]==0:
            if pr_sig1[snp_i]==0:
                curr_post_means1[snp_i] = 0
                curr_betas1[snp_i] = 0
            if pr_sig2[snp_i]==0:
                curr_post_means2[snp_i] = 0
                curr_betas2[snp_i] = 0
        else:
            start_i = max(0, snp_i - ld_radius)
            focal_i = min(ld_radius, snp_i)
            stop_i = min(m, snp_i + ld_radius + 1)
            D1_i = ld_dict1[snp_i]
            D2_i = ld_dict2[snp_i]
            local_betas1 = curr_betas1[start_i: stop_i]
            local_betas2 = curr_betas2[start_i: stop_i]
            local_betas1[focal_i] = 0
            local_betas2[focal_i] = 0
            num1 = beta_hats1[snp_i] - sp.dot(D1_i , local_betas1)
            num2 = beta_hats2[snp_i] - sp.dot(D2_i , local_betas2)

            v1 = pr_sig1[snp_i]/(Pi[0]+Pi[1])
            v2 = pr_sig2[snp_i]/(Pi[0]+Pi[2])
            v12 = rho*sp.sqrt(v1*v2)
            Sig = np.array([[v1,v12],[v12,v2]])
            detSig = v1*v2 - v12**2
            S = np.array([[v2/detSig+n1,-v12/detSig],[-v12/detSig,v1/detSig+n2]])
            detS = S[0,0]*S[1,1] - S[0,1]*S[1,0]
            S_inv = np.array([[S[1,1],-S[1,0]],[-S[0,1],S[0,0]]])/detS
            V1 = n1*num1
            V2 = n2*num2
            V = np.array([V1,V2])
            muV = np.dot(S_inv, V)
            VSV = np.dot(V.T, muV)/2.0
            inv_detSigS = 1.0/(detSig*detS)
            
            C1 = 1.0/(n1+1.0/v1) ##post var1
            C2 = 1.0/(n2+1.0/v2) ##post var2
            mu1 = V1*C1
            mu2 = V2*C2
            NC1 = C1*V1**2/2.0
            NC2 = C2*V2**2/2.0
            CR1 = 1.0/(n1*v1+1.0)
            CR2 = 1.0/(n2*v2+1.0)

            w11 = Pi[0]*sp.sqrt(inv_detSigS)
            w10 = Pi[1]*sp.sqrt(CR1)*sp.exp(NC1 - VSV)
            w01 = Pi[2]*sp.sqrt(CR2)*sp.exp(NC2 - VSV)
            w00 = Pi[3]*sp.exp(-VSV)
            wsum = w11 + w10 + w01 + w00
            postP = sp.array([w11, w11+w10, w11+w10+w01])*shrink_factor/wsum
            curr_post_means1[snp_i] = (w11*muV[0] + w10*mu1)/wsum
            curr_post_means2[snp_i] = (w11*muV[1] + w01*mu2)/wsum

            if rand_ps[i]<=postP[0]:
                [proposed_beta1,proposed_beta2] = np.random.multivariate_normal(muV, S_inv, size=1)[0]
            if rand_ps[i]>postP[0] and rand_ps[i]<=postP[1]:
                proposed_beta1 = stats.norm.rvs(mu1, C1, size=1)
                proposed_beta2 = 0
            if rand_ps[i]>postP[1] and rand_ps[i]<=postP[2]:
                proposed_beta1 = 0
                proposed_beta2 = stats.norm.rvs(mu2, C2, size=1)
            if rand_ps[i]>postP[2]:
                proposed_beta1 = 0
                proposed_beta2 = 0

            curr_betas1[snp_i] = proposed_beta1
            curr_betas2[snp_i] = proposed_beta2

    ########### update Pi ##########
    A1 = sp.sum((curr_betas1!=0) & (curr_betas2!=0))
    A2 = sp.sum((curr_betas1!=0) & (curr_betas2==0))
    A3 = sp.sum((curr_betas1==0) & (curr_betas2!=0))
    A4 = sp.sum((curr_betas1==0.0) & (curr_betas2==0.0))
    return {'proposed_betas1':curr_betas1, 'proposed_betas2':curr_betas2, 'curr_post_means1':curr_post_means1, 'curr_post_means2':curr_post_means2, 'A1':A1, 'A2':A2, 'A3':A3, 'A4':A4}


def bi_mcmc_all_chr(beta_hats1, beta_hats2, Pi, pr_sig1, pr_sig2, zj_p, start_betas1=None, start_betas2=None, n1=1000, n2=1000, ld_radius=100, ld_dict1=None, ld_dict2=None, h2_D1=None, h2_D2=None):
    """
    MCMC of non-infinitesimal model
    """
    #Pi = sp.random.dirichlet((alpha,alpha,alpha,alpha),1).flatten()
    m = len(beta_hats1)
    
    curr_betas1 = sp.copy(start_betas1)
    curr_post_means1 = sp.zeros(m)
#    avg_betas1 = sp.zeros(m)

    curr_betas2 = sp.copy(start_betas2)
    curr_post_means2 = sp.zeros(m)
#    avg_betas2 = sp.zeros(m)

#    Pi_traj = sp.zeros((4,num_iter+1))
#    Pi_traj[:,0] = Pi
#    s_traj1 = sp.zeros((m,num_iter+1))
#    s_traj2 = sp.zeros((m,num_iter+1))
#    s_traj1[:,0] = start_betas1
#    s_traj2[:,0] = start_betas2
#    m_traj1 = sp.zeros((m,num_iter))
#    m_traj2 = sp.zeros((m,num_iter))

    # Iterating over effect estimates in sequential order
    h2_est1 = max(0.00001,sp.sum(curr_betas1 ** 2))
    h2_est2 = max(0.00001,sp.sum(curr_betas2 ** 2))
    shrink_factor = min(1-zj_p, h2_D1/h2_est1, h2_D2/h2_est2)

    rand_ps = sp.random.random(m)
    iter_order = sp.arange(m)
    for i, snp_i in enumerate(iter_order):
        if pr_sig1[snp_i]==0 or pr_sig2[snp_i]==0:
            if pr_sig1[snp_i]==0:
                curr_post_means1[snp_i] = 0
                curr_betas1[snp_i] = 0
            if pr_sig2[snp_i]==0:
                curr_post_means2[snp_i] = 0
                curr_betas2[snp_i] = 0
        else:
            start_i = max(0, snp_i - ld_radius)
            focal_i = min(ld_radius, snp_i)
            stop_i = min(m, snp_i + ld_radius + 1)
            D1_i = ld_dict1[snp_i]
            D2_i = ld_dict2[snp_i]
            local_betas1 = curr_betas1[start_i: stop_i]
            local_betas2 = curr_betas2[start_i: stop_i]
            local_betas1[focal_i] = 0
            local_betas2[focal_i] = 0
            num1 = beta_hats1[snp_i] - sp.dot(D1_i , local_betas1)
            num2 = beta_hats2[snp_i] - sp.dot(D2_i , local_betas2)

            v1 = pr_sig1[snp_i]/(Pi[0]+Pi[1])
            v2 = pr_sig2[snp_i]/(Pi[0]+Pi[2])

            C1 = 1.0/(n1+1.0/v1) ##post var1
            C2 = 1.0/(n2+1.0/v2) ##post var2
            mu1 = n1*num1*C1
            mu2 = n2*num2*C2
            NC1 = C1*n1**2
            NC2 = C2*n2**2
            CR1 = 1.0/(n1*v1+1.0)
            CR2 = 1.0/(n2*v2+1.0)

            w11 = Pi[0]*sp.sqrt(CR1*CR2)
            w10 = Pi[1]*sp.sqrt(CR1)*sp.exp(-NC2*num2**2/2.0)
            w01 = Pi[2]*sp.sqrt(CR2)*sp.exp(-NC1*num1**2/2.0)
            w00 = Pi[3]*sp.exp(-NC2*num2**2/2.0-NC1*num1**2/2.0)
            wsum = w11 + w10 + w01 + w00
#            postP = [w11/wsum,(w11+w10)/wsum,(w11+w10+w01)/wsum]
#            outP = sp.array([w11,w10,w01,w00])/wsum
            postP = sp.array([w11, w11+w10, w11+w10+w01])*shrink_factor/wsum
            curr_post_means1[snp_i] = mu1*(w11+w10)/wsum
            curr_post_means2[snp_i] = mu2*(w11+w01)/wsum

            if rand_ps[i]<=postP[0]:
                proposed_beta1 = stats.norm.rvs(mu1, C1, size=1)
                proposed_beta2 = stats.norm.rvs(mu2, C2, size=1)
            if rand_ps[i]>postP[0] and rand_ps[i]<=postP[1]:
                proposed_beta1 = stats.norm.rvs(mu1, C1, size=1)
                proposed_beta2 = 0
            if rand_ps[i]>postP[1] and rand_ps[i]<=postP[2]:
                proposed_beta1 = 0
                proposed_beta2 = stats.norm.rvs(mu2, C2, size=1)
            if rand_ps[i]>postP[2]:
                proposed_beta1 = 0
                proposed_beta2 = 0

            curr_betas1[snp_i] = proposed_beta1
            curr_betas2[snp_i] = proposed_beta2

#            s_traj1[snp_i,k+1] = proposed_beta1
#            s_traj2[snp_i,k+1] = proposed_beta2
#
#            m_traj1[snp_i,k] = curr_post_means1[snp_i]
#            m_traj2[snp_i,k] = curr_post_means2[snp_i]

    ########### update Pi ##########
    A1 = sp.sum((curr_betas1!=0) & (curr_betas2!=0))
    A2 = sp.sum((curr_betas1!=0) & (curr_betas2==0))
    A3 = sp.sum((curr_betas1==0) & (curr_betas2!=0))
    A4 = sp.sum((curr_betas1==0.0) & (curr_betas2==0.0))
#    Pi = sp.random.dirichlet((alpha+A1,alpha+A2,alpha+A3,alpha+A4),1).flatten()
#    if k >= burn_in:
#        avg_betas1 += curr_post_means1 #Averaging over the posterior means instead of samples.
#        avg_betas2 += curr_post_means2
    return {'proposed_betas1':curr_betas1, 'proposed_betas2':curr_betas2, 'curr_post_means1':curr_post_means1, 'curr_post_means2':curr_post_means2, 'A1':A1, 'A2':A2, 'A3':A3, 'A4':A4}

def tri_mcmc_all_chr(beta_hats1, beta_hats2, beta_hats3, Pi, pr_sig1, pr_sig2, pr_sig3, zj_p, start_betas1=None, start_betas2=None, start_betas3=None, 
    n1=1000, n2=1000, n3=1000, ld_radius=100, ld_dict1=None, ld_dict2=None, ld_dict3=None, h2_D1=None, h2_D2=None, h2_D3=None):
    """
    MCMC of non-infinitesimal model
    """
    m = len(beta_hats1)
    
    curr_betas1 = sp.copy(start_betas1)
    curr_post_means1 = sp.zeros(m)

    curr_betas2 = sp.copy(start_betas2)
    curr_post_means2 = sp.zeros(m)

    curr_betas3 = sp.copy(start_betas3)
    curr_post_means3 = sp.zeros(m)

    h2_est1 = max(0.00001,sp.sum(curr_betas1 ** 2))
    h2_est2 = max(0.00001,sp.sum(curr_betas2 ** 2))
    h2_est3 = max(0.00001,sp.sum(curr_betas3 ** 2))
    shrink_factor = min(1-zj_p, h2_D1/h2_est1, h2_D2/h2_est2, h2_D3/h2_est3)

    # Iterating over effect estimates in sequential order
    rand_ps = sp.random.random(m)
    iter_order = sp.arange(m)
    for i, snp_i in enumerate(iter_order):
        start_i = max(0, snp_i - ld_radius)
        focal_i = min(ld_radius, snp_i)
        stop_i = min(m, snp_i + ld_radius + 1)
        D1_i = ld_dict1[snp_i]
        D2_i = ld_dict2[snp_i]
        D3_i = ld_dict3[snp_i]
        local_betas1 = curr_betas1[start_i: stop_i]
        local_betas2 = curr_betas2[start_i: stop_i]
        local_betas3 = curr_betas3[start_i: stop_i]
        local_betas1[focal_i] = 0
        local_betas2[focal_i] = 0
        local_betas3[focal_i] = 0
        num1 = beta_hats1[snp_i] - sp.dot(D1_i , local_betas1)
        num2 = beta_hats2[snp_i] - sp.dot(D2_i , local_betas2)
        num3 = beta_hats3[snp_i] - sp.dot(D3_i , local_betas3)

        v1 = pr_sig1[snp_i]/(Pi[0]+Pi[1]+Pi[2]+Pi[4])
        v2 = pr_sig2[snp_i]/(Pi[0]+Pi[1]+Pi[3]+Pi[5])
        v3 = pr_sig3[snp_i]/(Pi[0]+Pi[2]+Pi[3]+Pi[6])

        R1 = n1*v1
        R2 = n2*v2
        R3 = n3*v3

        db1 = 1 + 1.0/R1
        db2 = 1 + 1.0/R2
        db3 = 1 + 1.0/R3

        invC1 = sp.sqrt(1+R1)*sp.exp(-n1*num1**2/(2.0*db1))
        invC2 = sp.sqrt(1+R2)*sp.exp(-n2*num2**2/(2.0*db2))
        invC3 = sp.sqrt(1+R3)*sp.exp(-n3*num3**2/(2.0*db3))

        mu1 = num1/db1
        mu2 = num2/db2
        mu3 = num3/db3

        var1 = 1.0/(n1*db1)
        var2 = 1.0/(n2*db2)
        var3 = 1.0/(n3*db3)

        weights = np.zeros(8)
        weights[0] = Pi[0]
        weights[1] = Pi[1]*invC3
        weights[2] = Pi[2]*invC2
        weights[3] = Pi[3]*invC1
        weights[4] = Pi[4]*invC2*invC3
        weights[5] = Pi[5]*invC1*invC3
        weights[6] = Pi[6]*invC1*invC2
        weights[7] = Pi[7]*invC1*invC2*invC3
        wsum = np.sum(weights)
        postP = weights/wsum

        curr_post_means1[snp_i] = mu1*(postP[0]+postP[1]+postP[2]+postP[4])
        curr_post_means2[snp_i] = mu2*(postP[0]+postP[1]+postP[3]+postP[5])
        curr_post_means3[snp_i] = mu3*(postP[0]+postP[2]+postP[3]+postP[6])
        
        shrink_postP = np.cumsum(shrink_factor*postP[0:7])
        if rand_ps[i]<=shrink_postP[0]:
            proposed_beta1 = stats.norm.rvs(mu1, var1, size=1)
            proposed_beta2 = stats.norm.rvs(mu2, var2, size=1)
            proposed_beta3 = stats.norm.rvs(mu3, var3, size=1)
        if rand_ps[i]>shrink_postP[0] and rand_ps[i]<=shrink_postP[1]:
            proposed_beta1 = stats.norm.rvs(mu1, var1, size=1)
            proposed_beta2 = stats.norm.rvs(mu2, var2, size=1)
            proposed_beta3 = 0
        if rand_ps[i]>shrink_postP[1] and rand_ps[i]<=shrink_postP[2]:
            proposed_beta1 = stats.norm.rvs(mu1, var1, size=1)
            proposed_beta2 = 0
            proposed_beta3 = stats.norm.rvs(mu3, var3, size=1)
        if rand_ps[i]>shrink_postP[2] and rand_ps[i]<=shrink_postP[3]:
            proposed_beta1 = 0
            proposed_beta2 = stats.norm.rvs(mu2, var2, size=1)
            proposed_beta3 = stats.norm.rvs(mu3, var3, size=1)
        if rand_ps[i]>shrink_postP[3] and rand_ps[i]<=shrink_postP[4]:
            proposed_beta1 = stats.norm.rvs(mu1, var1, size=1)
            proposed_beta2 = 0
            proposed_beta3 = 0
        if rand_ps[i]>shrink_postP[4] and rand_ps[i]<=shrink_postP[5]:
            proposed_beta1 = 0
            proposed_beta2 = stats.norm.rvs(mu2, var2, size=1)
            proposed_beta3 = 0
        if rand_ps[i]>shrink_postP[5] and rand_ps[i]<=shrink_postP[6]:
            proposed_beta1 = 0
            proposed_beta2 = 0
            proposed_beta3 = stats.norm.rvs(mu3, var3, size=1)
        if rand_ps[i]>shrink_postP[6]:
            proposed_beta1 = 0
            proposed_beta2 = 0
            proposed_beta3 = 0

        curr_betas1[snp_i] = proposed_beta1
        curr_betas2[snp_i] = proposed_beta2
        curr_betas3[snp_i] = proposed_beta3
    ########### update Pi ##########
    A = np.zeros(8)
    A[0] = sp.sum((curr_betas1!=0) & (curr_betas2!=0) & (curr_betas3!=0))
    A[1] = sp.sum((curr_betas1!=0) & (curr_betas2!=0) & (curr_betas3==0))
    A[2] = sp.sum((curr_betas1!=0) & (curr_betas2==0) & (curr_betas3!=0))
    A[3] = sp.sum((curr_betas1==0) & (curr_betas2!=0) & (curr_betas3!=0))
    A[4] = sp.sum((curr_betas1!=0) & (curr_betas2==0) & (curr_betas3==0))
    A[5] = sp.sum((curr_betas1==0) & (curr_betas2!=0) & (curr_betas3==0))
    A[6] = sp.sum((curr_betas1==0) & (curr_betas2==0) & (curr_betas3!=0))
    A[7] = sp.sum((curr_betas1==0) & (curr_betas2==0) & (curr_betas3==0))
    return {'proposed_betas1':curr_betas1, 'proposed_betas2':curr_betas2, 'proposed_betas3':curr_betas3, 'curr_post_means1':curr_post_means1, 'curr_post_means2':curr_post_means2, 'curr_post_means3':curr_post_means3, 'A':A}
