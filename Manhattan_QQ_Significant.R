rm(list = ls())
options(stringsAsFactors=F)

suppressMessages(library("qqman"))
suppressMessages(library("data.table"))

args       = commandArgs(trailingOnly=TRUE)
aname      = args[1]
adir       = args[2]
mfolder    = args[3]
qfolder    = args[4]
sfolder    = args[5]
lfolder    = args[6]
jfolder    = args[7]
population = args[8]
locuszoom_path = args[9]
#locuszoom_path='/ysm-gpfs/pi/zhao/from_louise/bl537/locuszoom/bin/locuszoom'
#### No Need to Change ####

#mfolder = "/ysm-gpfs/pi/zhao/from_louise/bl537/GWAS/VA/ManhattanPlot"
mdir = paste0(mfolder, "/", aname, ".png")

#qfolder = "/ysm-gpfs/pi/zhao/from_louise/bl537/GWAS/VA/qqPlot"
qdir = paste0(qfolder, "/", aname, ".png")


a = as.data.frame(fread(adir, header=T, sep='auto'))
colnames(a)

b = a[, c("CHROM", "POS", "P")]
colnames(b) = c("CHR", "BP", "P")
b = b[order(b$CHR, b$BP),]

## manhattan plot ##
m = b[complete.cases(b),]

# Genomic Inflation Factor 
chisq  <- qchisq(1 - m$P, 1)
l <- round(median(chisq)/qchisq(0.5,1), digits=3)

png(mdir, width=1000, height=600, res=110, type='cairo')
manhattan(m, chr = "CHR", bp = "BP", p = "P",
          col = c("#4682B4", "#B4464B"), chrlabs = as.character(c(1:22)), 
          main = bquote(.(aname) ~ " (SNP-level)   " ~ lambda == .(l)),
          suggestiveline = -log10(1e-5), genomewideline = -log10(5e-8),
          ylim = c(0, max(-log10(m$P), -log10(5e-8), na.rm=T)+2),
          logp = TRUE,
          las=2, cex = 0.5)
dev.off()


## qq plot ##
png(qdir, width=600, height=600, type='cairo', res=120)
qq(m$P, main = bquote(.(aname) ~ " " ~ lambda == .(l)))
dev.off()




## Part3. Significant ##

flag = function(dd){
  flag = rep(1, nrow(dd))
  for (i in 1:(nrow(dd)-1)){
    if ((dd$BP[i+1]-dd$BP[i])>1e6) {flag[i+1] = flag[i] + 1} 
  }
  return(flag)
}

oodir = paste0(sfolder, "/", aname, ".txt")
lodir = paste0(jfolder, "/step1b_", aname)

d = a[a$PVALUE<5e-8, c('CHROM', 'POS', 'P', 'SNP')]
colnames(d) = c('CHR', 'BP', 'P', 'SNP')
d = d[order(d$CHR, d$BP),]

if (nrow(d)>=1){
  chromList = unique(d$CHR)
  for (i in 1:length(chromList)){
    index = chromList[i]
    
    dd = d[d$CHR==index,]
    
    if (nrow(dd)==1){Flag = 1}else{Flag = flag(dd)}
    o = cbind(dd, Flag)
    
    for (j in 1:max(o$Flag)) {
      oo = o[o$Flag==j, ]
      line = oo[oo$P==min(oo$P),]
      ll = paste0('module load R; ', locuszoom_path, ' --metal ', adir, ' --refsnp "',line$SNP,'" --chr ',line$CHR,' --flank 500kb --markercol SNP --pvalcol P --plotonly --prefix ', lfolder,'/zoom --pop ', population, ' --build hg19 --source 1000G_March2012')
      print(line)
      if (max(i,j)==1){
        write.table(line[, c(1:4)], file=oodir, quote = F, sep ="\t", row.names = F, col.names = T)
        #write.table(mess, file=lodir, quote = F, row.names = F, col.names = F)
        write.table(ll, file=lodir, quote = F, append = T, row.names = F, col.names = F)
      }else{
        write.table(line[, c(1:4)], file=oodir, quote = F, sep ="\t", row.names = F, col.names = F, append = T)
        write.table(ll, file=lodir, quote = F, append = T, row.names = F, col.names = F)
      }
    }
  }
}else{
  write.table('No PVALUE<5e-8!', file=oodir, quote = F, sep ="\t", row.names = F, col.names = F)
  write.table('No PVALUE<5e-8!', file=lodir, quote = F, append = T, row.names = F, col.names = F)
}






