
library(GEOquery)
library(survival)


expr_raw = read.table(file="~/Desktop/Recomb_2019/GSE96058_gene_expression_3273_samples_and_136_replicates_transformed.csv",header=TRUE,sep=',',row.names=1)
clin_raw = getGEO(filename="~/Desktop/Recomb_2019/GSE96058-GPL11154_series_matrix.txt.gz")
clin_dat = clin_raw@phenoData@data

genes = c('SLITRK3','GALP','OR4C13','VN1R4','LCE3C','STXBP5','IYD','ARID1B','SMR3A','DIP2B','LOC100101266','C19orf66','PARP12','EXOC1','CEL','SLC9A1','CCDC28B','PCMT1','C10orf131','HARS2','CELP','C7orf53','OR52E6','HBS1L','MEMO1P3','HLA-K','GVINP2','RPS10P20','HSPA8P1','GAPDHP45','ANKRD30BP2','RPS20P25','KRT18P62','OR7E10P','SCML2P1','RPL5P28','RPS27P12','PRR13P5','IGKV2D-23','GPS2','GPS2P1','GPS2:GPS2P1')
types = c('Alltype','Alltype(non-senior)','Normal','Normal(non-senior)','LumA','LumA(non-senior','LumB','LumB(non-senior)','Basal','Basal(non-senior)','Her2','Her2(non-senior)','ERpos_PRpos_HER2pos','ERpos_PRpos_HER2pos(non-senior)','ERneg_PRpos_HER2pos','ERneg_PRpos_HER2pos(non-senior)','ERpos_PRneg_HER2pos','ERpos_PRneg_HER2pos(non-senior)','ERpos_PRpos_HER2neg','ERpos_PRpos_HER2neg(non-senior)','ERneg_PRneg_HER2pos','ERneg_PRneg_HER2pos(non-senior)','ERpos_PRneg_HER2neg','ERpos_PRneg_HER2neg(non-senior)','ERneg_PRpos_HER2neg','ERneg_PRpos_HER2neg(non-senior)','ERneg_PRneg_HER2neg','ERneg_PRneg_HER2neg(non-senior)')
wald_out = matrix(NA,42,28)
row.names(wald_out)=genes;
colnames(wald_out)=types;

#Overall
cnt=0;
sigcnt=0;
for (nonsen in 0:1){
  for (gene in genes){
    if (gene %in% row.names(expr_raw)){
      cnt = cnt+1;
      message(paste0(gene,' found:'))
      x = as.numeric(expr_raw[gene,])
      t = as.numeric(clin_dat[,'overall survival days:ch1'])
      e = as.numeric(clin_dat[,'overall survival event:ch1'])
      a = as.numeric(gsub(".*: ","",clin_dat[,'characteristics_ch1.2']))
      names(x) = colnames(expr_raw)
      names(t) = clin_dat[,'title']
      names(e) = clin_dat[,'title']
      names(a) = clin_dat[,'title']
      pats = intersect(names(t),names(x))
      if(nonsen == 1){
        pats2 = intersect(names(a)[(a>21) & (a<65)],pats)
      }else if(nonsen == 0){
        pats2 = intersect(names(a)[(a>0) & (a<Inf)],pats)
      }
      if (length(unique((x[pats2]>quantile(x[pats2], c(.5), na.rm=TRUE))))>1){
        #km_fit <- survfit(Surv(t[pats2], e[pats2]) ~ x[pats2]>quantile(x[pats2], c(.5), na.rm=TRUE))
        #plot(km_fit)
        #lines(km_fit[2], col=2)
        #survdiff(Surv(t[pats2], e[pats2]) ~ x[pats2]>quantile(x[pats2], c(.5), na.rm=TRUE))
        tmp = summary(coxph(Surv(t[pats2], e[pats2]) ~ x[pats2]>quantile(x[pats2], c(.5), na.rm=TRUE)))
        message(paste0(names(tmp$waldtest)[1],': ',tmp$waldtest[1],': ',names(tmp$waldtest)[2],': ',tmp$waldtest[2],': ',names(tmp$waldtest)[3],': ',tmp$waldtest[3]))
        if(tmp$waldtest[3]<0.05){
          sigcnt=sigcnt+1;
        }
        wald_out[gene,nonsen+1] = tmp$waldtest[3];
      }else{
        message('gene expression is homogeneous')
      }
      if(tmp$waldtest[3]<0.05){sigcnt=sigcnt+1;}
    }else{
      message(paste0(gene,' not in file'))
    }
  }
}


# Type
for(nonsen in 0:1){
  for (p50typ in 0:4){
    for (gene in genes){
      if (gene %in% row.names(expr_raw)){
        cnt=cnt+1;
        message(paste0(gene,' found:'))
        x = as.numeric(expr_raw[gene,])
        t = as.numeric(gsub(".*: ","",clin_dat[,'overall survival days:ch1']))
        e = as.numeric(gsub(".*: ","",clin_dat[,'overall survival event:ch1']))
        a = as.numeric(gsub(".*: ","",clin_dat[,'characteristics_ch1.2']))
        p50 = gsub(".*: ","",clin_dat[,'characteristics_ch1.20'])
        names(x) = colnames(expr_raw)
        names(t) = clin_dat[,'title']
        names(e) = clin_dat[,'title']
        names(a) = clin_dat[,'title']
        names(p50) = clin_dat[,'title']
        pats = intersect(names(t),names(x))
        if(nonsen == 1){
          pats2 = intersect(names(a)[(a>21) & (a<65)],pats)
        }else if(nonsen == 0){
          pats2 = intersect(names(a)[(a>0) & (a<Inf)],pats)
        }
        if(p50typ==0){pats3 = intersect(names(p50)[p50=="Normal"],pats2)}
        else if(p50typ==1){pats3 = intersect(names(p50)[p50=="LumA"],pats2)}
        else if(p50typ==2){pats3 = intersect(names(p50)[p50=="LumB"],pats2)}
        else if(p50typ==3){pats3 = intersect(names(p50)[p50=="Basal"],pats2)}
        else if(p50typ==4){pats3 = intersect(names(p50)[p50=="Her2"],pats2)}
        if (length(unique((x[pats3]>quantile(x[pats3], c(.5), na.rm=TRUE))))>1){
          #km_fit <- survfit(Surv(t[pats3], e[pats3]) ~ x[pats3]>quantile(x[pats3], c(.5), na.rm=TRUE))
          #plot(km_fit)
          #lines(km_fit[2], col=2)
          #survdiff(Surv(t[pats2], e[pats2]) ~ x[pats2]>quantile(x[pats2], c(.5), na.rm=TRUE))
          tmp = summary(coxph(Surv(t[pats3], e[pats3]) ~ x[pats3]>quantile(x[pats3], c(.5), na.rm=TRUE)))
          message(paste0(names(tmp$waldtest)[1],': ',tmp$waldtest[1],': ',names(tmp$waldtest)[2],': ',tmp$waldtest[2],': ',names(tmp$waldtest)[3],': ',tmp$waldtest[3]))
          if(tmp$waldtest[3]<0.05){
            sigcnt=sigcnt+1;
          }
          wald_out[gene,p50typ*2+3+nonsen] = tmp$waldtest[3];
        }else{
          message('gene expression is homogeneous')
        }
      }else{
        message(paste0(gene,' not in file'))
      }
    }
  }
}

# receptor status
for(nonsen in 0:1){
  for(stat in 0:7){
    for (gene in genes){
      if (gene %in% row.names(expr_raw)){
        cnt=cnt+1;
        message(paste0(gene,' found:'))
        x = as.numeric(expr_raw[gene,])
        t = as.numeric(gsub(".*: ","",clin_dat[,'overall survival days:ch1']))
        e = as.numeric(gsub(".*: ","",clin_dat[,'overall survival event:ch1']))
        a = as.numeric(gsub(".*: ","",clin_dat[,'characteristics_ch1.2']))
        er = as.numeric(gsub(".*: ","",clin_dat[,'characteristics_ch1.16']))
        pr = as.numeric(gsub(".*: ","",clin_dat[,'characteristics_ch1.17']))
        her2 = as.numeric(gsub(".*: ","",clin_dat[,'characteristics_ch1.18']))
        names(x) = colnames(expr_raw)
        names(t) = clin_dat[,'title']
        names(e) = clin_dat[,'title']
        names(a) = clin_dat[,'title']
        names(er) = clin_dat[,'title']
        names(pr) = clin_dat[,'title']
        names(her2) = clin_dat[,'title']
        pats = intersect(names(t),names(x))
        if(nonsen == 1){
          pats2 = intersect(names(a)[(a>21) & (a<65)],pats)
        }else if(nonsen == 0){
          pats2 = intersect(names(a)[(a>0) & (a<Inf)],pats)
        }
        if(stat == 0){pats3 = intersect(intersect(names(er)[er==1],intersect(names(pr==1),names(her2==1))),pats2)}
        else if(stat == 1){pats3 = intersect(intersect(names(er)[er==0],intersect(names(pr)[pr==1],names(her2)[her2==1])),pats2)}
        else if(stat == 2){pats3 = intersect(intersect(names(er)[er==1],intersect(names(pr)[pr==0],names(her2)[her2==1])),pats2)}
        else if(stat == 3){pats3 = intersect(intersect(names(er)[er==0],intersect(names(pr)[pr==0],names(her2)[her2==1])),pats2)}
        else if(stat == 4){pats3 = intersect(intersect(names(er)[er==0],intersect(names(pr)[pr==0],names(her2)[her2==1])),pats2)}
        else if(stat == 5){pats3 = intersect(intersect(names(er)[er==1],intersect(names(pr)[pr==0],names(her2)[her2==0])),pats2)}
        else if(stat == 6){pats3 = intersect(intersect(names(er)[er==0],intersect(names(pr)[pr==1],names(her2)[her2==0])),pats2)}
        else if(stat == 7){pats3 = intersect(intersect(names(er)[er==0],intersect(names(pr)[pr==0],names(her2)[her2==0])),pats2)}
        if (length(unique((x[pats3]>quantile(x[pats3], c(.5), na.rm=TRUE))))>1 && sum(e[pats3]) > 2 && length(pats3) > 2){
          #km_fit <- survfit(Surv(t[pats3], e[pats3]) ~ x[pats3]>quantile(x[pats3], c(.5), na.rm=TRUE))
          #plot(km_fit)
          #lines(km_fit[2], col=2)
          #survdiff(Surv(t[pats2], e[pats2]) ~ x[pats2]>quantile(x[pats2], c(.5), na.rm=TRUE))
          tmp = summary(coxph(Surv(t[pats3], e[pats3]) ~ x[pats3]>quantile(x[pats3], c(.5), na.rm=TRUE)))
          message(paste0(names(tmp$waldtest)[1],': ',tmp$waldtest[1],': ',names(tmp$waldtest)[2],': ',tmp$waldtest[2],': ',names(tmp$waldtest)[3],': ',tmp$waldtest[3]))
          if(tmp$waldtest[3]<0.05){
            sigcnt=sigcnt+1;
          }
          wald_out[gene,stat*2+13+nonsen] = tmp$waldtest[3];
        }else{
          message('gene expression is homogeneous')
        }
      }else{
        message(paste0(gene,' not in file'))
      }
    }
  }
}
print(paste0(sigcnt,' significant genes found from ',cnt,' genes subtype combinations in dataset'));
wald_out_adj = t(apply(wald_out,1,function(x) p.adjust(x, method = "BH")))
print(paste0(sum(wald_out_adj<0.05,na.rm=TRUE),' significant genes (with multiple testing correction) found from ',sum(!is.na(wald_out_adj)),' genes subtype combinations in dataset'));
write.table(formatC(wald_out,format='e',digits=5),file='~/Desktop/Recomb_2019/Supplementary_table_1v2.txt',quote=FALSE,sep='\t',col.names=TRUE,row.names=TRUE)
write.table(formatC(wald_out_adj,format='e',digits=5),file='~/Desktop/Recomb_2019/Supplementary_table_BH_1v2.txt',quote=FALSE,sep='\t',col.names=TRUE,row.names=TRUE)


