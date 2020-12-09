## read hg file, partition file, bim file, summary stats file ##
get_inputs = function(dir_simu1, numr, output_dir_full, val_dir, h2_D1, h2_D2, h_partition1, h_partition2, anno_assign1, anno_assign2){
	## dir_simu1: path to simulated files
	## numr: index of replication, e.g. the first replication when numr = 1
	## output_dir_full: the path for saving generated files (e.g. reformmated summary statistics, prior files)
	## val_dir: path to validation files
	## h2_D1: heritability of D1
	## h2_D2: heritability of D2
	## h_partition1: the heritability in A, B, complement of union of A and B in D1
	## h_partition2: the heritability in A, B, complement of union of A and B in D2
	## anno_assign1: annotation assignment in D1
	## anno_assign2: annotation assignment in D2
	out_full = paste0(output_dir_full,numr) # the folder for the numr th replication
	dir.create(out_full) # create directory
	out_val = paste0(val_dir,numr) # the folder for the numr th validation data
	dir.create(out_val) # create directory
	setwd(paste0(dir_simu1,numr)) # change dir to simudated data of numr th replication
	## heritability partition of D1
	hgf1 = read.table(h_partition1,header=F)[,1]
	partf1 = read.table(anno_assign1,header=T)
	## heritability partition of D2
	hgf2 = read.table(h_partition2,header=F)[,1]
	partf2 = read.table(anno_assign2,header=F)
	## by writing the following codes I already assume the validataion data (undivided) are in paste0(dir_simu1,numr) and named test.
	bim_test = read.table('test.bim',header=F)
	bim_train_full = read.table('train.bim',header=F) # the snp information of the training data (get a snp list of summary stats)
	fam_test = read.table('test.fam',header=F)
	summary_stats_full = read.table('train.ss.assoc.linear',header=T) # summary statistics

	# reformatting summary stats
	summary_stats_full = cbind(summary_stats_full[,c(1,2)],bim_train_full[,c(5,6)],summary_stats_full[,c(3,7,9)])
	non_na_full = !is.na(summary_stats_full[,7])
	summary_stats_full = summary_stats_full[non_na_full,]
	summary_stats_full[,6] = exp(summary_stats_full[,6])
	summary_stats_full[,1] = paste0('chr',summary_stats_full[,1])
	# write the reformatted summary stats to paste0(out_full,'/','Input_sumstats.txt')
	write.table(summary_stats_full,paste0(out_full,'/','Input_sumstats.txt'),quote=F,row.names=F,col.names=c('hg19chrc','snpid','a1','a2','bp','or','p'))

	# per-SNP heritability for D1 (h2 prior file, the first prior)
	partf_mat1 = partf1[,c(2,3)]
	partf_mat1 = cbind(!(partf_mat1[,1]|partf_mat1[,2]),partf_mat1)
	partf_mat1 = as.matrix(partf_mat1)
	MM <- diag(colSums(partf_mat1))
	for(i in 1:(nrow(MM)-1)){
	  for(j in (i+1):ncol(MM)){
	    MM[i,j] <- sum(partf_mat1[,i]&partf_mat1[,j])
	    MM[j,i] <- MM[i,j]
	  }
	}
	tau1 = solve(MM,hgf1) # solving the linear equation I wrote today
	h_snp1 = partf_mat1%*%tau1
	h_snp1[h_snp1<0] = rep(min(h_snp1[h_snp1>0]),sum(h_snp1<0)) # if negative elements exist, substitute with smallest positive element
	h_snp_file1 = data.frame(bim_train_full[,1],bim_train_full[,2],h_snp1)
	write.table(h_snp_file1,paste0(out_full,'/','h_snp_D1.txt'),quote=F,row.names=F,col.names=F)
	## generating the pT files (the second priors)
	H0 <- h2
	N0 <- nrow(partf_mat1)
	sig2V <- partf_mat1%*%tau
	annot1V <- apply(partf_mat1,1,paste,collapse='')
	N_T <- table(annot1V)
	N_TV <- as.numeric(N_T[annot1V])
	H_TV <- N_TV*sig2V
	dir.create(paste0(output_dir_full,numr,'/p_priors'))
	dir.create(paste0(output_dir_half,numr,'/p_priors'))
	pv <- c(1,0.3,0.1,0.03,0.01,0.003,0.001,0.0003,0.0001,0.00003,0.00001)
	for(p0 in pv){
	  pr_p <- p0*N0*H_TV/(N_TV*H0)#(p0*N0/H0)*M_TV*sig2V/N_TV
	  sig2 <- sig2V
	  m1 <- min(pr_p[pr_p>0])
	  m2 <- min(sig2[sig2>0])
	  pr_p[pr_p>1] = rep(1,sum(pr_p>1))
	  pr_p[pr_p<0] = rep(m1,sum(pr_p<0))
	  sig2[sig2<0] = rep(m2,sum(sig2<0))
	  output <- data.frame(bim_train_full[,1],bim_train_full[,2],pr_p,sig2)
	  ## write the pT files to folder
	  write.table(output,paste0(output_dir_full,numr,'/p_priors/p0_',p0,'.txt'),quote=F,row.names=F,col.names=F)
	}


	# per-SNP heritability for D2 (generating the prior for joint modeling of two diseases)
	partf_mat2 = partf2[,c(2,3)]
	partf_mat2 = cbind(!(partf_mat2[,1]|partf_mat2[,2]),partf_mat2)
	partf_mat2 = as.matrix(partf_mat2)
	MM <- diag(colSums(partf_mat2))
	for(i in 1:(nrow(MM)-1)){
	  for(j in (i+1):ncol(MM)){
	    MM[i,j] <- sum(partf_mat2[,i]&partf_mat2[,j])
	    MM[j,i] <- MM[i,j]
	  }
	}
	tau2 = solve(MM,hgf2) 
	h_snp2 = partf_mat2%*%tau2
	h_snp2[h_snp2<0] = rep(min(h_snp2[h_snp2>0]),sum(h_snp2<0))
	h_snp_file12 = data.frame(bim_train_full[,1],bim_train_full[,2],h_snp1,h_snp2)
	write.table(h_snp_file12,paste0(out_full,'/','h_snp_D1_D2.txt'),quote=F,row.names=F,col.names=F)

	# dividing validation set to two random sets of individuals
	## to use plink to divide individuals, we will need two list of individual IDs
	n <- nrow(fam_test)
	n1 <- floor(n/2)
	s1 <- sample(1:n,n1)
	s2 <- setdiff(1:n,s1)
	sel1 <- fam_test[s1,c(1,2)]
	sel2 <- fam_test[s2,c(1,2)]
	## write two sets of individual IDs (to be used for plink, see dividing_testing() below)
	write.table(sel1,paste0(out_val,'/','sel1.txt'),quote=F,col.names=F,row.names=F)
	write.table(sel2,paste0(out_val,'/','sel2.txt'),quote=F,col.names=F,row.names=F)
}

divide_testing = function(simu_dir,val_dir,numr){
	seg1 = 'source ~/.bashrc; '
	seg2 = paste0('cd ',val_dir,numr,'; ')
	seg3 = paste0('plink --bfile ',simu_dir,numr,'/test --keep sel1.txt --make-bed --out test_sel1;')
	seg4 = paste0('plink --bfile ',simu_dir,numr,'/test --keep sel2.txt --make-bed --out test_sel2;')
	paste0(seg1,seg2,c(seg3,seg4))
}

get_coord = function(val_dir,input_dir,numr,N){
	seg1 = 'source ~/.bashrc; '
	seg2 = 'cd /net/zhao/yh367/prediction_codes/; '
	seg3 = paste0('python coord_genotypes.py --gf=',val_dir,numr,'/','test_sel1 --ssf=', input_dir,numr,'/Input_sumstats.txt --ssf_format=BASIC --N=',N,' --out=', input_dir,numr,'/coord_sel1 --vgf=', val_dir,numr,'/','test_sel1 --vbim=', val_dir,numr,'/','test_sel1.bim')
	seg4 = paste0('python coord_genotypes.py --gf=',val_dir,numr,'/','test_sel2 --ssf=', input_dir,numr,'/Input_sumstats.txt --ssf_format=BASIC --N=',N,' --out=', input_dir,numr,'/coord_sel2 --vgf=', val_dir,numr,'/','test_sel2 --vbim=', val_dir,numr,'/','test_sel2.bim')
	paste0(seg1,seg2,c(seg3,seg4))
}

get_ldpred_cmds = function(input_dir,numr,output_dir,r,h2,N){
	dir.create(paste0(output_dir,numr))
	dir.create(paste0(output_dir,numr,'/ldpred_sel1/'))
	dir.create(paste0(output_dir,numr,'/ldpred_sel2/'))
	pv <- c(1,0.3,0.1,0.03,0.01,0.003,0.001,0.0003,0.0001,0.00003,0.00001)
	ldpred_sel1 <- paste0('source ~/.bashrc; cd /net/zhao/yh367/prediction_codes; python LDpred_no_output.py --coord=',input_dir,numr,'/coord_sel1  --ld_radius=',r,' --PS=',pv,' --local_ld_file_prefix=',input_dir,'ldfile.',numr,'_sel1 --N=',N,'  --H2=',h2,'  --out=',output_dir,numr,'/ldpred_sel1/test_p_',pv)
	ldpred_sel2 <- paste0('source ~/.bashrc; cd /net/zhao/yh367/prediction_codes; python LDpred_no_output.py --coord=',input_dir,numr,'/coord_sel2  --ld_radius=',r,' --PS=',pv,' --local_ld_file_prefix=',input_dir,'ldfile.',numr,'_sel2 --N=',N,'  --H2=',h2,'  --out=',output_dir,numr,'/ldpred_sel2/test_p_',pv)
	c(ldpred_sel1,ldpred_sel2)
}

get_h2_cmds = function(input_dir,numr,output_dir,r,h2,N){
	dir.create(paste0(output_dir,numr,'/h2_sel1/'))
	dir.create(paste0(output_dir,numr,'/h2_sel2/'))
	pv <- c(1,0.3,0.1,0.03,0.01,0.003,0.001,0.0003,0.0001,0.00003,0.00001)
	h2_sel1 <- paste0('source ~/.bashrc; cd /net/zhao/yh367/prediction_codes; python AnnoPred_h2.py --coord=',input_dir,numr,'/coord_sel1  --ld_radius=',r,' --ps=',pv,'  --PS=',input_dir,numr,'/h_snp.txt  --local_ld_file_prefix=',input_dir,'ldfile.',numr,'_sel1 --N=',N,'  --H2=',h2,'  --out=',output_dir,numr,'/h2_sel1/test_p_',pv)
	h2_sel2 <- paste0('source ~/.bashrc; cd /net/zhao/yh367/prediction_codes; python AnnoPred_h2.py --coord=',input_dir,numr,'/coord_sel2  --ld_radius=',r,' --ps=',pv,'  --PS=',input_dir,numr,'/h_snp.txt  --local_ld_file_prefix=',input_dir,'ldfile.',numr,'_sel2 --N=',N,'  --H2=',h2,'  --out=',output_dir,numr,'/h2_sel2/test_p_',pv)
	c(h2_sel1,h2_sel2)
}

get_p_inf_cmds = function(input_dir,numr,output_dir,r,h2,N){
	dir.create(paste0(output_dir,numr,'/p_inf_sel1/'))
	dir.create(paste0(output_dir,numr,'/p_inf_sel2/'))
	pv <- c(1,0.3,0.1,0.03,0.01,0.003,0.001,0.0003,0.0001,0.00003,0.00001)
	p_inf_sel1 <- paste0('source ~/.bashrc; cd /net/zhao/yh367/prediction_codes; python AnnoPred_pT.py --coord=',input_dir,numr,'/coord_sel1  --ld_radius=',r,' --ps=',pv,'  --PS=',input_dir,numr,'/p_priors/p0_',pv,'.txt  --local_ld_file_prefix=',input_dir,'ldfile.',numr,'_sel1 --N=',N,'  --H2=',h2,'  --out=',output_dir,numr,'/p_inf_sel1/test_p_',pv)
	p_inf_sel2 <- paste0('source ~/.bashrc; cd /net/zhao/yh367/prediction_codes; python AnnoPred_pT.py --coord=',input_dir,numr,'/coord_sel2  --ld_radius=',r,' --ps=',pv,'  --PS=',input_dir,numr,'/p_priors/p0_',pv,'.txt  --local_ld_file_prefix=',input_dir,'ldfile.',numr,'_sel2 --N=',N,'  --H2=',h2,'  --out=',output_dir,numr,'/p_inf_sel2/test_p_',pv)
	c(p_inf_sel1,p_inf_sel2)
}

get_ldpred_cmds_no_h2 = function(input_dir,numr,output_dir,r,N){
	dir.create(paste0(output_dir,numr))
	dir.create(paste0(output_dir,numr,'/ldpred_sel1/'))
	dir.create(paste0(output_dir,numr,'/ldpred_sel2/'))
	pv <- c(1,0.3,0.1,0.03,0.01,0.003,0.001,0.0003,0.0001,0.00003,0.00001)
	ldpred_sel1 <- paste0('source ~/.bashrc; cd /net/zhao/yh367/prediction_codes; python LDpred_no_output.py --coord=',input_dir,numr,'/coord_sel1  --ld_radius=',r,' --PS=',pv,' --local_ld_file_prefix=',input_dir,'ldfile.',numr,'_sel1 --N=',N,'  --out=',output_dir,numr,'/ldpred_sel1/test_p_',pv)
	ldpred_sel2 <- paste0('source ~/.bashrc; cd /net/zhao/yh367/prediction_codes; python LDpred_no_output.py --coord=',input_dir,numr,'/coord_sel2  --ld_radius=',r,' --PS=',pv,' --local_ld_file_prefix=',input_dir,'ldfile.',numr,'_sel2 --N=',N,'  --out=',output_dir,numr,'/ldpred_sel2/test_p_',pv)
	c(ldpred_sel1,ldpred_sel2)
}

get_h2_cmds_no_h2 = function(input_dir,numr,output_dir,r,N){
	dir.create(paste0(output_dir,numr,'/h2_sel1/'))
	dir.create(paste0(output_dir,numr,'/h2_sel2/'))
	pv <- c(1,0.3,0.1,0.03,0.01,0.003,0.001,0.0003,0.0001,0.00003,0.00001)
	h2_sel1 <- paste0('source ~/.bashrc; cd /net/zhao/yh367/prediction_codes; python AnnoPred_h2.py --coord=',input_dir,numr,'/coord_sel1  --ld_radius=',r,' --ps=',pv,'  --PS=',input_dir,numr,'/h_snp.txt  --local_ld_file_prefix=',input_dir,'ldfile.',numr,'_sel1 --N=',N,'  --out=',output_dir,numr,'/h2_sel1/test_p_',pv)
	h2_sel2 <- paste0('source ~/.bashrc; cd /net/zhao/yh367/prediction_codes; python AnnoPred_h2.py --coord=',input_dir,numr,'/coord_sel2  --ld_radius=',r,' --ps=',pv,'  --PS=',input_dir,numr,'/h_snp.txt  --local_ld_file_prefix=',input_dir,'ldfile.',numr,'_sel2 --N=',N,'  --out=',output_dir,numr,'/h2_sel2/test_p_',pv)
	c(h2_sel1,h2_sel2)
}

get_pT_cmds_no_h2 = function(input_dir,numr,output_dir,r,N){
	dir.create(paste0(output_dir,numr,'/p_inf_sel1/'))
	dir.create(paste0(output_dir,numr,'/p_inf_sel2/'))
	pv <- c(1,0.3,0.1,0.03,0.01,0.003,0.001,0.0003,0.0001,0.00003,0.00001)
	p_inf_sel1 <- paste0('source ~/.bashrc; cd /net/zhao/yh367/prediction_codes; python AnnoPred_pT.py --coord=',input_dir,numr,'/coord_sel1  --ld_radius=',r,' --ps=',pv,'  --PS=',input_dir,numr,'/p_priors/p0_',pv,'.txt  --local_ld_file_prefix=',input_dir,'ldfile.',numr,'_sel1 --N=',N,'  --out=',output_dir,numr,'/p_inf_sel1/test_p_',pv)
	p_inf_sel2 <- paste0('source ~/.bashrc; cd /net/zhao/yh367/prediction_codes; python AnnoPred_pT.py --coord=',input_dir,numr,'/coord_sel2  --ld_radius=',r,' --ps=',pv,'  --PS=',input_dir,numr,'/p_priors/p0_',pv,'.txt  --local_ld_file_prefix=',input_dir,'ldfile.',numr,'_sel2 --N=',N,'  --out=',output_dir,numr,'/p_inf_sel2/test_p_',pv)
	c(p_inf_sel1,p_inf_sel2)
}

get_all_prs_cmds = function(input_dir,numr,output_dir,r,N){
	dir.create(paste0(output_dir,numr))
	dir.create(paste0(output_dir,numr,'/P_T_sel1/'))
	dir.create(paste0(output_dir,numr,'/P_T_sel2/'))
	rv <- c(1,0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1)
	P_T_sel1 <- paste0('source ~/.bashrc; cd /net/zhao/yh367/prediction_codes; python LD_pruning_thres_edit.py --coord=',input_dir,numr,'/coord_sel1  --ld_radius=',r,' --max_r2=',rv,'  --out=',output_dir,numr,'/P_T_sel1/test_sel1')
	P_T_sel2 <- paste0('source ~/.bashrc; cd /net/zhao/yh367/prediction_codes; python LD_pruning_thres_edit.py --coord=',input_dir,numr,'/coord_sel2  --ld_radius=',r,' --max_r2=',rv,'  --out=',output_dir,numr,'/P_T_sel2/test_sel2')
	c(P_T_sel1,P_T_sel2)
}

get_results = function(output_dir, dis, method1, method2){
  res_1_auc <- matrix(0,length(method1),11)
  res_2_auc <- matrix(0,length(method2),11)
  pv <- c(1,0.3,0.1,0.03,0.01,0.003,0.001,0.0003,0.0001,0.00003,0.00001)
  for(k in 1:length(method1)){
      setwd(paste0(output_dir,method1[k]))
      fl.list <- paste0(dis,'_p_',pv,'_LDpred-auc.txt')
      for(i in 1:11){
        fl <- readLines(fl.list[i])
        aucL <- unlist(strsplit(fl[length(fl)-1],': '))
        res_1_auc[k,i] <- as.numeric(aucL[2])
      }
      
      setwd(paste0(output_dir,method2[k]))
      fl.list <- paste0(dis,'_p_',pv,'_LDpred-auc.txt')
      for(i in 1:11){
        fl <- readLines(fl.list[i])
        aucL <- unlist(strsplit(fl[length(fl)-1],': '))
        res_2_auc[k,i] <- as.numeric(aucL[2])
      }  
  }
  avg_aucs <- rep(0,length(method1))
  for(i in 1:length(method1)){
    am1 <- which.max(res_2_auc[i,])
    am2 <- which.max(res_1_auc[i,])
    avg_aucs[i] <- (res_1_auc[i,am1]+res_2_auc[i,am2])/2
  }
  output1 <- rep(0,length(method1))
  output2 <- rep(0,length(method1))
  for(k in 1:length(method1)){
    
    fl <- readLines(paste0(output_dir,method1[k],'/',dis,'_p_1_LDpred-auc-inf.txt'))
    aucL <- unlist(strsplit(fl[length(fl)],': '))
    output1[k] <- as.numeric(aucL[2])
    
    
    fl <- readLines(paste0(output_dir,method2[k],'/',dis,'_p_1_LDpred-auc-inf.txt'))
    aucL <- unlist(strsplit(fl[length(fl)],': '))
    output2[k] <- as.numeric(aucL[2])
  }
  inf_output = (output1 + output2)/2

  list(results = round(avg_aucs,4),inf_output=round(inf_output,4))
}


get_Results_tuning = function(resultsDirectory, dis, method1, method2){
  l = length(method1)
  if(l==1){
    pv <- c(1,0.3,0.1,0.03,0.01,0.003,0.001,0.0003,0.0001,0.00003,0.00001)
    res_1_auc <- rep(0,11)
    res_2_auc <- rep(0,11)
    setwd(paste0(resultsDirectory,method1))
    fl.list <- paste0(dis,'_p_',pv,'_LDpred-auc.txt')
    for(i in 1:11){
      fl <- readLines(fl.list[i])
      aucL <- unlist(strsplit(fl[length(fl)-1],': '))
      res_1_auc[i] <- as.numeric(aucL[2])
    }
    
    setwd(paste0(resultsDirectory,method2))
    fl.list <- paste0(dis,'_p_',pv,'_LDpred-auc.txt')
    for(i in 1:11){
      fl <- readLines(fl.list[i])
      aucL <- unlist(strsplit(fl[length(fl)-1],': '))
      res_2_auc[i] <- as.numeric(aucL[2])
    }
    am1 <- which.max(res_2_auc)
    am2 <- which.max(res_1_auc)
    avg_auc <- (res_1_auc[am1]+res_2_auc[am2])/2

    fl <- readLines(paste0(resultsDirectory,method1,'/',dis,'_p_1_LDpred-auc-inf.txt'))
    aucL <- unlist(strsplit(fl[length(fl)],': '))
    output1 <- as.numeric(aucL[2])
    
    fl <- readLines(paste0(resultsDirectory,method2,'/',dis,'_p_1_LDpred-auc-inf.txt'))
    aucL <- unlist(strsplit(fl[length(fl)],': '))
    output2 <- as.numeric(aucL[2])

    inf_output = (output1 + output2)/2
    results = list(auc=round(avg_auc,4),inf_auc=round(inf_output,4),sel1=c(am1,pv[am1]),sel2=c(am2,pv[am2]))
    return(results)
  }
  if(l==2){
    res_1_auc <- matrix(0,2,11)
    res_2_auc <- matrix(0,2,11)
    #method1 <- c('h2_no_zero_1','p_inf1')
    #method2 <- c('h2_no_zero_2','p_inf2')
    pv <- c(1,0.3,0.1,0.03,0.01,0.003,0.001,0.0003,0.0001,0.00003,0.00001)
    for(k in 1:2){
        setwd(paste0(resultsDirectory,method1[k]))
        fl.list <- paste0(dis,'_p_',pv,'_LDpred-auc.txt')
        for(i in 1:11){
          fl <- readLines(fl.list[i])
          aucL <- unlist(strsplit(fl[length(fl)-1],': '))
          res_1_auc[k,i] <- as.numeric(aucL[2])
        }
        
        setwd(paste0(resultsDirectory,method2[k]))
        fl.list <- paste0(dis,'_p_',pv,'_LDpred-auc.txt')
        for(i in 1:11){
          fl <- readLines(fl.list[i])
          aucL <- unlist(strsplit(fl[length(fl)-1],': '))
          res_2_auc[k,i] <- as.numeric(aucL[2])
        }  
    }
    am1 <- which.max(res_2_auc)
    am2 <- which.max(res_1_auc)
    avg_auc <- (res_1_auc[am1]+res_2_auc[am2])/2
    sel1_p1 <- which.max(res_1_auc[1,])
    sel1_p2 <- which.max(res_1_auc[2,])
    if(res_1_auc[1,sel1_p1]>res_1_auc[2,sel1_p2]){
      sel1_method = 1
      sel1_p = sel1_p1
    }else{
      sel1_method = 2
      sel1_p = sel1_p2
    }
    sel2_p1 <- which.max(res_2_auc[1,])
    sel2_p2 <- which.max(res_2_auc[2,])
    if(res_2_auc[1,sel2_p1]>res_2_auc[2,sel2_p2]){
      sel2_method = 1
      sel2_p = sel2_p1
    }else{
      sel2_method = 2
      sel2_p = sel2_p2
    }
    avg_auc_another = (res_1_auc[sel2_method,sel2_p] + res_2_auc[sel1_method,sel1_p])/2
    if(avg_auc == avg_auc_another){
     print('we r cool')
    }else{
     print('shit happens')
     return('False results')
    }
    output1 <- rep(0,2)
    output2 <- rep(0,2)
    for(k in 1:2){
      fl <- readLines(paste0(resultsDirectory,method1[k],'/',dis,'_p_1_LDpred-auc-inf.txt'))
      aucL <- unlist(strsplit(fl[length(fl)],': '))
      output1[k] <- as.numeric(aucL[2])
      
      fl <- readLines(paste0(resultsDirectory,method2[k],'/',dis,'_p_1_LDpred-auc-inf.txt'))
      aucL <- unlist(strsplit(fl[length(fl)],': '))
      output2[k] <- as.numeric(aucL[2])
    }
    inf_output = (output1[sel2_method] + output2[sel1_method])/2
    results = list(auc=round(avg_auc_another,4),inf_auc=round(inf_output,4),sel1=c(sel1_method,pv[sel1_p]),sel2=c(sel2_method,pv[sel2_p]))
    return(results)
  }
}

get_P_T_PRS_all_PRS_significant = function(results_dir1,results_dir2,dis){
	pv <- c(1,0.3,0.1,0.03,0.01,0.003,0.001,0.0003,0.0001,3E-5,1E-5,1E-6,1E-7,5E-8,1E-8)
	rv <- c(1,0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1)
	pv <- paste(pv)
	pv[8] <- '0.0003'
	pv[9] <- '0.0001'
	rv <- paste(rv)
	rv[1] <- '1.0'
	
	output1 <- matrix(0,length(rv),length(pv))
	setwd(results_dir1)
	for(i in 1:length(rv)){
	  auc_fl <- paste0(dis,'_sel1_r_',rv[i],'_p_',pv,'_auc.txt')
	  for(j in 1:length(pv)){
	    fl <- readLines(auc_fl[j])
	    aucL <- unlist(strsplit(fl[length(fl)],': '))
	    output1[i,j] <- as.numeric(aucL[2])
	  }   
	}
	
	output2 <- matrix(0,length(rv),length(pv))
	setwd(results_dir2)
	for(i in 1:length(rv)){
	  auc_fl <- paste0(dis,'_sel2_r_',rv[i],'_p_',pv,'_auc.txt')
	  for(j in 1:length(pv)){
	    fl <- readLines(auc_fl[j])
	    aucL <- unlist(strsplit(fl[length(fl)],': '))
	    output2[i,j] <- as.numeric(aucL[2])
	  }   
	}
	
	am1 <- which.max(output2)
	am2 <- which.max(output1)
	P_T = round((output1[am1]+output2[am2])/2,4)
	PRS_all = round((output1[1,1]+output2[1,1])/2,4)
	PRS_significant = round((output1[1,length(pv)-1]+output2[1,length(pv)-1])/2,4)
	return(list(PRS_P_T=P_T,PRS_all=PRS_all,PRS_significant=PRS_significant))
}
