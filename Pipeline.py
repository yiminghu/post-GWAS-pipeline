#### Hope this rocks ####


import os
import os.path
import sys  
import argparse
    
parser = argparse.ArgumentParser()
parser.add_argument('--study', default=None, type=str,
                    help="Study name.")
parser.add_argument('--population', default='EUR', type=str,
                    help="Population.")
parser.add_argument('--sumstats', default=None, type=str,
                    help="Summary statistics.")
parser.add_argument('--grandfolder', default=None, type=str,
                    help="Grand folder.")

args = parser.parse_args()

if not os.path.exists(args.grandfolder+'/joblist'):
    os.makedirs(args.grandfolder+'/joblist')

if not os.path.exists(args.grandfolder+'/LDSC'):
    os.makedirs(args.grandfolder+'/LDSC')
if not os.path.exists(args.grandfolder+'/LDSC/sumstats'):
    os.makedirs(args.grandfolder+'/LDSC/sumstats')
if not os.path.exists(args.grandfolder+'/LDSC/Results'):
    os.makedirs(args.grandfolder+'/LDSC/Results')

if not os.path.exists(args.grandfolder+'/Standard'):
    os.makedirs(args.grandfolder+'/Standard')
if not os.path.exists(args.grandfolder+'/Standard/ManhattanPlot'):
    os.makedirs(args.grandfolder+'/Standard/ManhattanPlot')
if not os.path.exists(args.grandfolder+'/Standard/qqPlot'):
    os.makedirs(args.grandfolder+'/Standard/qqPlot')
if not os.path.exists(args.grandfolder+'/Standard/Significant'):
    os.makedirs(args.grandfolder+'/Standard/Significant')
if not os.path.exists(args.grandfolder+'/Standard/LocusZoom'):
    os.makedirs(args.grandfolder+'/Standard/LocusZoom')

Study                        = args.study
SumStats_file                = args.sumstats
jobfolder                    = args.grandfolder+'/joblist/'+Study
LDSC_SumStats_folder         = args.grandfolder+'/LDSC/sumstats/'+Study
LDSC_Results_folder          = args.grandfolder+'/LDSC/Results/'+Study
ManhattanPlot_folder         = args.grandfolder+'/Standard/ManhattanPlot/'+Study
qqPlot_folder                = args.grandfolder+'/Standard/qqPlot/'+Study
sig_folder                   = args.grandfolder+'/Standard/Significant/'+Study
LocusZoom_folder             = args.grandfolder+'/Standard/LocusZoom/'+Study
Population                   = args.population

    
if not os.path.exists(jobfolder):
    os.makedirs(jobfolder)
if not os.path.exists(LDSC_SumStats_folder):
    os.makedirs(LDSC_SumStats_folder)
if not os.path.exists(LDSC_Results_folder):
    os.makedirs(LDSC_Results_folder)
    os.makedirs(LDSC_Results_folder+'/Tier1')
    os.makedirs(LDSC_Results_folder+'/Tier2')
    os.makedirs(LDSC_Results_folder+'/Tier3')
    os.makedirs(LDSC_Results_folder+'/Heritability')
if not os.path.exists(ManhattanPlot_folder):
    os.makedirs(ManhattanPlot_folder)
if not os.path.exists(qqPlot_folder):
    os.makedirs(qqPlot_folder)
if not os.path.exists(sig_folder):
    os.makedirs(sig_folder)
if not os.path.exists(LocusZoom_folder):
    os.makedirs(LocusZoom_folder)

    
#### Step 1. Manhattan and QQ Plots ####
outdir = jobfolder+"/step1a_"+Study
outfile=open(outdir,'w')
a = "source ~/.bashrc; module load R; aname="+Study+"; SumStats_file="+SumStats_file+"; mfolder="+ManhattanPlot_folder+"; qfolder="+qqPlot_folder+"; sfolder="+sig_folder+"; lfolder="+LocusZoom_folder+"; jfolder="+jobfolder+"; Population="+Population+"; Rscript --vanilla /gpfs/ysm/pi/zhao/from_louise/bl537/GWAS/VA/Pipeline/Script/Manhattan_QQ_Significant.R ${aname} ${SumStats_file} ${mfolder} ${qfolder} ${sfolder} ${lfolder} ${jfolder} ${Population};"
outfile.write(a+'\n')
outfile.close()

#### Step2a. LDSC sumstats ####
outdir = jobfolder+"/step2a_"+Study
outfile = open(outdir,'w')
a = "source ~/.bashrc; SumStats_file="+SumStats_file+"; outfolder="+LDSC_SumStats_folder+"; aname="+Study+"; cd /ysm-gpfs/pi/zhao/from_louise/bl537/ldsc; python2.7 munge_sumstats.py --out ${outfolder}/${aname}  --sumstats ${SumStats_file};"        
outfile.write(a+'\n')
outfile.close()
    
#### Step2b. LDSC regression 3 tiers ####
Tier = "Tier1"
outdir = jobfolder+"/step2b_"+Study+"_"+Tier
outfile = open(outdir,'w')
BaselineDIR = "/ysm-gpfs/pi/zhao/from_louise/ql68/Software/ldsc/Annotations/EUR/Baseline/baseline"
AnnotDIR    = "/ysm-gpfs/pi/zhao/from_louise/ql68/Software/ldsc/Annotations/EUR/GenoSkyline_Plus"
InputDIR    = "/ysm-gpfs/pi/zhao/from_louise/ql68/Software/ldsc/Input"
a = "source ~/.bashrc; Tier="+Tier+"; outfolder="+LDSC_Results_folder+"/${Tier}; aname="+Study+"; sumstatsDIR="+LDSC_SumStats_folder+"; BaselineDIR="+BaselineDIR+"; AnnotDIR="+AnnotDIR+"; InputDIR="+InputDIR+"; cd /ysm-gpfs/pi/zhao/from_louise/bl537/ldsc; python2.7 ldsc.py --h2 ${sumstatsDIR}/${aname}.sumstats.gz --ref-ld-chr ${BaselineDIR}.,${AnnotDIR}/GenoSkyline_Plus_${Tier}. --w-ld-chr ${InputDIR}/EUR/weights/weights. --frqfile-chr ${InputDIR}/EUR/genotype/1000G.mac5eur. --overlap-annot --print-coefficients --out ${outfolder}/${aname};"
outfile.write(a+'\n')
outfile.close()
    
Tier = "Tier2"
outdir = jobfolder+"/step2b_"+Study+"_"+Tier
outfile = open(outdir,'w')
BaselineDIR = "/ysm-gpfs/pi/zhao/from_louise/ql68/Software/ldsc/Annotations/EUR/Baseline/baseline"
AnnotDIR    = "/ysm-gpfs/pi/zhao/from_louise/ql68/Software/ldsc/Annotations/EUR/GenoSkyline_Plus"
InputDIR    = "/ysm-gpfs/pi/zhao/from_louise/ql68/Software/ldsc/Input"
a = "source ~/.bashrc; Tier="+Tier+"; outfolder="+LDSC_Results_folder+"/${Tier}; aname="+Study+"; sumstatsDIR="+LDSC_SumStats_folder+"; BaselineDIR="+BaselineDIR+"; AnnotDIR="+AnnotDIR+"; InputDIR="+InputDIR+"; cd /ysm-gpfs/pi/zhao/from_louise/bl537/ldsc; python2.7 ldsc.py --h2 ${sumstatsDIR}/${aname}.sumstats.gz --ref-ld-chr ${BaselineDIR}.,${AnnotDIR}/GenoSkyline_Plus_${Tier}. --w-ld-chr ${InputDIR}/EUR/weights/weights. --frqfile-chr ${InputDIR}/EUR/genotype/1000G.mac5eur. --overlap-annot --print-coefficients --out ${outfolder}/${aname};"
outfile.write(a+'\n')
outfile.close()
    
Tier = "Tier3"
outdir = jobfolder+"/step2b_"+Study+"_"+Tier
outfile = open(outdir,'w')
BaselineDIR = "/ysm-gpfs/pi/zhao/from_louise/ql68/Software/ldsc/Annotations/EUR/Baseline/baseline"
AnnotDIR    = "/ysm-gpfs/pi/zhao/from_louise/ql68/Software/ldsc/Annotations/EUR/GenoSkyline_Plus"
InputDIR    = "/ysm-gpfs/pi/zhao/from_louise/ql68/Software/ldsc/Input"
a = "source ~/.bashrc; Tier="+Tier+"; outfolder="+LDSC_Results_folder+"/${Tier}; aname="+Study+"; sumstatsDIR="+LDSC_SumStats_folder+"; BaselineDIR="+BaselineDIR+"; AnnotDIR="+AnnotDIR+"; InputDIR="+InputDIR+"; cd /ysm-gpfs/pi/zhao/from_louise/bl537/ldsc; python2.7 ldsc.py --h2 ${sumstatsDIR}/${aname}.sumstats.gz --ref-ld-chr ${BaselineDIR}.,${AnnotDIR}/GenoSkyline_Plus_${Tier}. --w-ld-chr ${InputDIR}/EUR/weights/weights. --frqfile-chr ${InputDIR}/EUR/genotype/1000G.mac5eur. --overlap-annot --print-coefficients --out ${outfolder}/${aname};"
outfile.write(a+'\n')
outfile.close()
    
#### Step2c. LDSC regression heritability ####
sub = "Heritability"
outdir = jobfolder+"/step2c_"+Study+"_"+sub
outfile = open(outdir,'w')
InputDIR    = "/ysm-gpfs/pi/zhao/from_louise/ql68/Software/ldsc/Input"
a = "source ~/.bashrc; sub="+sub+"; outfolder="+LDSC_Results_folder+"/${sub}; aname="+Study+"; sumstatsDIR="+LDSC_SumStats_folder+"; InputDIR="+InputDIR+"; cd /ysm-gpfs/pi/zhao/from_louise/bl537/ldsc; python2.7 ldsc.py --h2 ${sumstatsDIR}/${aname}.sumstats.gz --ref-ld-chr ${InputDIR}/EUR/eur_w_ld_chr/ --w-ld-chr ${InputDIR}/EUR/eur_w_ld_chr/ --out ${outfolder}/${aname};"
outfile.write(a+'\n')
outfile.close()   
