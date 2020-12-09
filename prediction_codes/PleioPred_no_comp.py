#!/usr/bin/env python

from argparse import ArgumentParser
from os.path import isfile, isdir, join
from sys import exit
import numpy as np
#from pleiopred import prior_generating, coord_trimmed, pre_sumstats
import pred_main_no_comparison

# Create the master argparser and returns the argparser object
def get_argparser():
  parser = ArgumentParser(prog="PleioPred", 
                          description="Genetic Risk Prediction by joint modeling of multiple diseases and functional annotation.")
  ## Input Files
  
  ## Parameters
  parser.add_argument('--N1', required=True, type=int,
                      help="Sample size of the first disease GWAS")
  parser.add_argument('--N2', required=True, type=int,
                      help="Sample size of the second disease GWAS")
  parser.add_argument('--rho', required=True, type=float,
                      help="Tuning parameter in (-1,1)"
                           ", the genetic correlation between diseases")
  parser.add_argument('--local_ld_prefix', required=True,
                      help="A local LD file name prefix"
                           ", will be created if not present")
  parser.add_argument('--hfile', required=True,
                      help="per-SNP heritability estimation")
  parser.add_argument('--ld_radius', type=int,
                      help="If not provided, will use the number of SNPs in" 
                           " common divided by 3000")
  parser.add_argument('--coord_D1', required=True, 
                      help="Output H5 File for coord_genotypes of D1")
  parser.add_argument('--coord_D2', required=True, 
                      help="Output H5 File for coord_genotypes of D2")
  parser.add_argument('--out', default="PleioPred_out",
                      help="Output filename prefix for AnnoPred")
  parser.add_argument('--user_h1', type=float,
                      help="User-provided heritability estimation for D1")
  parser.add_argument('--user_h2', type=float,
                      help="User-provided heritability estimation for D2")

  return parser

def process_args(args):
  pdict = {}
  
  pdict['coord_D1'] = args.coord_D1
  pdict['coord_D2'] = args.coord_D2
  pdict['N1'] = args.N1
  pdict['N2'] = args.N2

  if (args.rho>-1 and args.rho<1):
    pdict['rho'] = args.rho
  else:
    exit("Tuning parameter needs to be in (-1,1)!")

  pdict['ld_radius'] = args.ld_radius

  pdict['local_ld_prefix'] = args.local_ld_prefix
  pdict['hfile'] = args.hfile
  pdict['out'] = args.out
  pdict['user_h1'] = args.user_h1
  pdict['user_h2'] = args.user_h2
  return pdict

def main(pdict):
  print(pdict)
  pred_main_no_comparison.main(pdict) 

if __name__ == '__main__':
  args = get_argparser().parse_args()
  main(process_args(args))
