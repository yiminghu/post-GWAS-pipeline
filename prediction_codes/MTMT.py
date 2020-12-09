#!/usr/bin/env python

from argparse import ArgumentParser
from os.path import isfile, isdir, join
from sys import exit
import numpy as np
#from pleiopred import prior_generating, coord_trimmed, pre_sumstats
import multi_glasso_ver0

# Create the master argparser and returns the argparser object
def get_argparser():
  parser = ArgumentParser(prog="MTMT", 
                          description="A Multi-Task Learning approach to jointly predicting Multiple Traits")
  ## Input Files
  
  ## Parameters
  parser.add_argument('--N', required=True,
                      help="Sample size of training data of multiple traits, with the first one being the one of interest.")
  parser.add_argument('--init_betas', required=True,
                      help="path to initial values (AnnoPred-inf scores)")
  parser.add_argument('--num_iter', type=int, default=60, 
                      help="Number of iterations for MCMC, default to 60.")
  parser.add_argument('--local_ld_prefix', required=True,
                      help="A local LD file name prefix will be created if not present")
  parser.add_argument('--ld_radius', type=int,
                      help="If not provided, will use the number of SNPs in common divided by 3000")
  parser.add_argument('--multi_coord', required=True, 
                      help="Output H5 File for testing/reference genotypes and summary stats")
  parser.add_argument('--out', default="MTMT_out",
                      help="Output filename prefix for MTMT")
  return parser

def process_args(args):
  pdict = {}
  
  pdict['multi_coord'] = args.multi_coord
  pdict['N'] = [int(item) for item in args.N.split(',')]
  pdict['ld_radius'] = args.ld_radius
  pdict['local_ld_prefix'] = args.local_ld_prefix
  pdict['out'] = args.out
  pdict['num_iter'] = args.num_iter
  pdict['init_betas'] = args.init_betas
  return pdict

def main(pdict):
  print(pdict)
  multi_glasso_ver0.main(pdict) 

if __name__ == '__main__':
  args = get_argparser().parse_args()
  main(process_args(args))
