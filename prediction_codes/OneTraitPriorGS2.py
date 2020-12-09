#!/usr/bin/env python

### generate prior file from h5py file directly ###
### generate_h2_pT generates two prior files from the results of LDSC and a fixed annotation file ###
### generate_h2_from_user generates one prior file from the user provided prior file ###
import h5py
import os
from collections import Counter
from collections import defaultdict
import datetime
import math
from argparse import ArgumentParser
from os.path import isfile, isdir, join
from sys import exit
import numpy as np

# Create the master argparser and returns the argparser object
def get_argparser():
  parser = ArgumentParser(prog="PriorGenerating", 
                          description="Generating two types of priors from Functional Annotations.")
  parser.add_argument('--h5py_file', required=True,
                      help="Path to coord file"
                           ", will be created if not present")
  parser.add_argument('--LDSC_results_file', required=True,
                      help="Path to corresponding LDSC results")
  parser.add_argument('--output_h2', required=True,
                      help="Path to generated h2 prior files")
  parser.add_argument('--output_pT', required=True,
                      help="Path and prefix to generated pT prior files")
  parser.add_argument('--PS', type=str)  

  return parser

def process_args(args):
  pdict = {}
  pdict['h5py_file'] = args.h5py_file
  pdict['LDSC_results_file'] = args.LDSC_results_file
  pdict['output_h2'] = args.output_h2
  pdict['output_pT'] = args.output_pT
  pdict['PS'] = [float(item) for item in args.PS.split(',')]
  return pdict


def generate_h2_pT(pdict):
    h5py_file = pdict['h5py_file']
    LDSC_results_file = pdict['LDSC_results_file']
    output_h2 = pdict['output_h2']
    PS = pdict['PS']
    output_pT = pdict['output_pT']
    # generate two types of prior files
    ### load the fixed input file ###
    if len(PS)==1:
        PS = [PS]
    h5f1 = h5py.File('/net/zhao/yh367/PleioPred/ref/GS2.h5','r')
    annot = h5f1['annot'][:]
    h5f1.close()
    h5f2 = h5py.File('/net/zhao/yh367/PleioPred/ref/1000G_SNP_info.h5','r')
    snp_chr = h5f2['snp_chr'][:]
    h5f2.close()
    ### get the snp list from h5py ###
    chromosomes_list = ['chrom_%d'%(x) for x in range(1,23)]
    chromosomes_list.append('chrom_X')
    
    df = h5py.File(h5py_file,'r')
    cord_data_g = df['cord_data']
    
    SNPids = []
    for chrom_str in chromosomes_list:
        if chrom_str in cord_data_g.keys():
            g = cord_data_g[chrom_str]
            #Filter monomorphic SNPs (SNPs with variance equal to 0)
            snp_stds = g['snp_stds_ref'][...]
            snp_stds = snp_stds.flatten()
            ok_snps_filter = snp_stds>0
            pval_derived_betas = g['betas'][...]
            pval_derived_betas = pval_derived_betas[ok_snps_filter]
            sids = g['sids'][...]
            SNPids = np.append(SNPids,sids[ok_snps_filter])
    num_snps = len(SNPids)
    ### overlap with SNP in annot files ###
    stt1 = np.in1d(snp_chr[:,2],SNPids)
    ant1 = annot[stt1]
    snp_chr1 = snp_chr[stt1]
    ### check order ###
    if sum(snp_chr1[:,2]==SNPids)==len(SNPids):
        print 'Good!'
    else:
        print 'Shit happens, sorting ant1 to have the same order as SNPids'
        O1 = np.argsort(snp_chr1[:,2])
        O2 = np.argsort(SNPids)
        O3 = np.argsort(O2)
        ant1 = ant1[O1][O3]

    ### load LDSC results ###
    LD_results = np.genfromtxt(LDSC_results_file,dtype=None,names=True)
    
    tau0 = LD_results['Coefficient']
  
    ### get heritability  ###
    sig2_0 = np.dot(ant1,tau0)
    
    ### adjust for minus terms ###
    sig2_0[sig2_0<0] = np.repeat(min(sig2_0[sig2_0>0]),np.sum(sig2_0<0))
    np.sum(sig2_0)
    
    ### save prior file (h2) ###
    h2_out = []
    for i in range(len(sig2_0)):
        h2_out.append(str(snp_chr1[:,0][i])+' '+str(snp_chr1[:,2][i])+' '+str(sig2_0[i])+'\n')
    #np.savetxt(output_h2,(snp_chr1[:,0],snp_chr1[:,1],sig2_0),fmt="%s")
    ff = open(output_h2,"w")
    ff.writelines(h2_out)
    ff.close()
    print 'h2 prior file saved at '+output_h2
    ### start calculating p_T ###
    M = np.empty(annot.shape[1])
    for i in range(len(M)):
        M[i] = np.sum(np.logical_and(annot[:,0],annot[:,i]))
    bgt = datetime.datetime.now()
    M_T = defaultdict(int)
    for i in range(annot.shape[0]):
        tup_i = tuple(annot[i])
        M_T[tup_i] += 1
    edt = datetime.datetime.now()
    print edt-bgt
    bgt = datetime.datetime.now()
    N_T = defaultdict(int)
    for i in range(ant1.shape[0]):
        tup_i = tuple(ant1[i])
        N_T[tup_i] += 1
    edt = datetime.datetime.now()
    print edt-bgt


    H0 = np.dot(M,tau0)
    N0 = float(len(SNPids))
    sig2V = np.dot(ant1,tau0)

    # N_T = {x:annotV1.count(x) for x in annotV1}
    
    M_TV = np.empty(ant1.shape[0])
    N_TV = np.empty(ant1.shape[0])
    for i in range(ant1.shape[0]):
        tup_i = tuple(ant1[i])
        M_TV[i] = M_T[tup_i]
        N_TV[i] = N_T[tup_i]

    for ps in PS:
        pr_p = (ps*N0/H0)*M_TV*sig2V/N_TV
        sig2 = M_TV*sig2V/N_TV
        m1 = min(pr_p[pr_p>0])
        m2 = min(sig2[sig2>0])
        pr_p[pr_p<0] = np.repeat(m1,np.sum(pr_p<0))
        sig2[sig2<0] = np.repeat(m2,np.sum(sig2<0))
        pr_p[pr_p>1] = np.repeat(1,np.sum(pr_p>1))
        pT_out = []
        for i in range(len(sig2)):
            pT_out.append(str(snp_chr1[:,0][i])+' '+str(snp_chr1[:,2][i])+' '+str(pr_p[i])+' '+str(sig2[i])+'\n')
        ff = open(output_pT+'_'+str(ps)+'_file.txt',"w")
#        ff = open(output_pT,"w")
        ff.writelines(pT_out)
        ff.close()
        print 'pT prior files saved at ' + output_pT+'_'+str(ps)+'_file.txt'

    print 'Suggested LD radius: ' + str(math.ceil(num_snps/3000.0)) 
    return math.ceil(num_snps/3000.0)

def main(pdict):
  print(pdict)
  generate_h2_pT(pdict) 

if __name__ == '__main__':
  args = get_argparser().parse_args()
  main(process_args(args))





