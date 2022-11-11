#' @Title HEART  find DE gene in single cell data.
#'
#' @param counts a data frame or a matrix for gene expression data.
#' the row is cell, the col is gene.
#' @param cell.data data frame containing informations of individual cells.
#'  It contains at least 2 columns, one for the cell name, which is the same as
#'  the row name of cell.data and counts, and one for the grouping variable of each cell,
#'  a dichotomous variable.
#' @param ln a logical value indicating whether to perform log transformation on counts.
#' Default value is FALSE.
#' @param mt  a character string. The adjustment methods include "holm", "hochberg", "hommel",
#' "bonferroni", "BH", "BY", "fdr", "none".
#' Default value is 'holm'.
#' @param thr  the threshold indication whether a gene is a DE gene.
#' Default value is 0.05.  If the adjusted p-value of a gene is less than 0.05,
#' it is considered a DE gene.
#' @param all.results  a logical value indicating whether to output all genes' results.
#' Default value is FALSE. If FALSE, only the results for the DE genes are output.

#' @return A data frame have gene name, p-value, adjusted p-value,
#'  the ratio of each part between groups.


heart=function(counts,cell.data,ln=F,
               mt='holm',thr=0.05,all.results=F){
  thr.mt='soft'
  counts=as.matrix(counts)
  if (ln==T) {new.counts=log2(1+counts)
  }else {new.counts=counts}

  # library(matrixTests)
  comb.lev<-function(data,cell.data,p.thr){
    data.g=rowSums(data!=0)

    x<-data[data.g!=0,]  #  filter zero expression
    group=cell.data$group
    g.f<-levels(as.factor(group))

    x1=as.matrix(x[,group==g.f[1]])
    x2=as.matrix(x[,group==g.f[2]])
    n1=ncol(x1)
    n2=ncol(x2)

    non0<-function(y){length(which(y!=0))}
    s.sqre<-function(y){
      sum((y-mean(y))^2)
    }

    # non-zero part mean
    new.m<-function(y){ifelse(non0(y)==0,return(0),return(mean(y[which(y!=0)])))}
    # non-zero sum square
    new.sqre<-function(y){ifelse(non0(y)==0,return(0),return(s.sqre(y[which(y!=0)])))}
    #non-zero????
    m1=apply(x1, 1, non0)  ## by row
    m2=apply(x2, 1, non0)
    #mean-x,mean-y
    x.mean<-apply(x1, 1, new.m)
    y.mean<-apply(x2, 1, new.m)
    #sum????-x,sum????-y
    x.sqre<-apply(x1, 1, new.sqre)
    y.sqre<-apply(x2, 1, new.sqre)

    p1=rowSums(x1!=0)/ncol(x1)
    p2=rowSums(x2!=0)/ncol(x2)
    p=rowSums(x!=0)/ncol(x)


    p_var<-sqrt(p*(1-p)*((1/n1)+(1/n2)))
    z=(p1-p2)/p_var

    p.sta<-z^2 # (chisq(1))
    t.sta<-(y.mean-x.mean)/sqrt((m1+m2)*(x.sqre+y.sqre)/((m1+m2-2)*m1*m2))
    non_0.lev=function(data,group){
      group=as.factor(group)
      rsl.1=apply(t(data), 2, function(y){
        matrixTests::col_brownforsythe(y[which(y!=0)],group[which(y!=0)])})

      rsl.2=do.call(rbind,rsl.1)
      rsl.2
    }
    lev=non_0.lev(x,cell.data$group)
    lev.sta=lev$statistic
    lev.p=lev$pvalue  ### warning与NA：sd=0导致的

    ########  #######

    t.matrix=data.frame(x.mean,y.mean,m1,m2,x.sqre,y.sqre,n1,n2,t.sta)
    t.matrix=within(t.matrix,{
      cond.sqre=1*I(exp(x.sqre+y.sqre)==1)  ##  x.sqre=0 or y.sqre=0
      cond.mean=1*I(exp(x.mean*y.mean)==1)   ##  x.mean=0 or y.mean=0
      # multi=fBasics::rowMaxs(cbind(x.mean/y.mean,y.mean/x.mean))
      multi=fBasics::rowMins(cbind(x.mean/y.mean,y.mean/x.mean))
      # logf=cond.sqre*((1-cond.mean)*(ifelse(multi==0,0,1/multi))
      #                 +cond.mean*abs(x.mean-y.mean))
      logf=(1-cond.mean)*(ifelse(multi==0,0,1/multi))+cond.mean*abs(x.mean-y.mean)   ## log fc

      t.end.sta=cond.sqre*logf+(1-cond.sqre)*(ifelse(0*t.sta=='NaN',0,t.sta))
      # t.df=m1+m2-2
      fisher.t=ifelse((m1+m2-2)>1,2-2*pt(abs(t.end.sta),df=m1+m2-2),
                      2-2*pt(abs(t.end.sta),df=2))

    })


    ##  2022_add_fold_

    t2=data.frame(p1,p2,
                  x.mean,y.mean,
                  x.sqre,y.sqre)

    max.ra=function(x,y){

      ll=cbind(ifelse(y!=0,x/y,NA),
               ifelse(x!=0,y/x,NA))

      # l.end=rowMaxs(ll,na.rm=T)
      l.end=apply(ll,2,function(z){max(z,na.rm = T)})

      l.end
    }


    ratio.mat=data.frame(gene.nm=rownames(t2),
                         prop.ratio=max.ra(t2$p1,t2$p2),
                         mean.ratio=max.ra(t2$x.mean,t2$y.mean),
                         var.ratio=max.ra(t2$x.sqre,t2$y.sqre))

    ##fisher ???\u{38fa}??????: Q(m,n),\u{4aa}?\u{23c}???
    fisher.p=2-2*pnorm(abs(z))    #li ????p
    fisher.t=2-2*pt(abs(t.sta),df=m1+m2-2) #  original version of Tp-value
    fisher.w=lev.p  #li ????var


    ## 20210615 add--p-value modify
    #   merge statistics
    mg.sta=data.frame(p1,p2,m1,m2,p,p_var,z,p.sta,fisher.p,x1.m=x.mean,x2.m=y.mean,
                      x.ss=x.sqre,y.ss=y.sqre,
                      t.s=t.matrix$t.end.sta,t.p=t.matrix$fisher.t,lev.sta,lev.p)
    ### mg.sta只是用来看结果，并未参与下一步的计算
    ## 基于mg.sta  属于相关的统计量
    var1=c('p1','p2','x1.m','x2.m')
    x.var=ifelse(mg.sta$m1<2,NA,mg.sta$x.ss/(mg.sta$m1-1))
    y.var=ifelse(mg.sta$m2<2,NA,mg.sta$y.ss/(mg.sta$m2-1))

    rsl.stas=mg.sta[var1]
    var.nm=c(paste0(rep(c('exp.prop','mean(non-zero)','var(non-zero)'),each=2),
                    '_G',rep(g.f,3)))
    rsl.stas=cbind(rsl.stas,x.var,y.var)
    colnames(rsl.stas)<-var.nm
    rsl.stas=data.frame(gene.nm=rownames(rsl.stas),rsl.stas)  ###输出的summary
    # 下一步的计算，只用到了sta中的信息

    sta=data.frame(fisher.p,change.t=t.matrix$fisher.t,fisher.w)

    l=data.frame(sta,df=apply(sta, 1, function(y){sum(abs(y)<100000,na.rm=T)}))
    cal=-2*log(sta)
    cal=data.frame(cal,cal.sum=apply(cal, 1, function(y){sum(y,na.rm=T)}))
    cal$change.pval=1-pchisq(cal$cal.sum,df=l$df)  ##.pvalue

    fisher.sta=-2*log(fisher.t)-2*log(fisher.w)-2*log(fisher.p) #chisq(6)
    fish.pval<-1-pchisq(fisher.sta,df=6)  #


    # fisher.sta=cal$cal.sum  #  20210615 new fisher combine statistics
    # fish.pvalue=cal$p_cal   ## 20210615 new p-calue with theoretical df

    ### 二分法求最佳自由度
    find.df=function(data,a,b,eps=1e-8){
      sta=na.omit(data)
      sta=sta[sta<1e10& sta>0]
      lnx=log(sta)

      fx=function(q){
        fx=digamma(q)-mean(lnx)+log(2)
        fx
      }
      ##

      if(fx(a)*fx(b)>0)
        list(fail='Wrong interval range!')
      else{
        repeat {
          if( abs(b-a)<=eps) break
          p=(a+b)/2
          if (fx(a)*fx(p)<0) b<-p else a<-p
        }
        p
      }
    }
    orig.df=find.df(fisher.sta,0.01,20)   ## original method's freedom
    change.df=find.df(cal$cal.sum,0.01,20) ## change method's freedom


    #######
    fish.padj=1-pchisq(fisher.sta,df=orig.df) ##
    change.padj=1-pchisq(cal$cal.sum,df=change.df) ##

    # rsl.fsh=cbind(fisher.sig,fisher.p.sig,fisher.t.sig,fisher.var.sig,fisher.sta,
    #               fish.pvalue,f.p.adj,fisher.p,fisher.t,fisher.w)

    rsl.fsh=cbind(cal,change.padj,fisher.sta,fish.pval,fish.padj)
    rsl.fsh=data.frame(gene.nm=rownames(rsl.fsh),rsl.fsh)
    p.var=c('gene.nm','fish.pval','fish.padj','change.pval','change.padj')
    rsl.fsh=rsl.fsh[p.var]

    ###????
    rslt<-list(test.p=rsl.fsh,
               gene.summary=rsl.stas,
               ratio=ratio.mat)
    rslt
  }
  # lev.sm1=rslt$fisher
  lev1<-comb.lev(new.counts,cell.data)
  p.thr=0.01
  end.output=function(list,var,mt,p.thr,thr){
    lev.sm1<-list$test.p
    gene.summary=list$gene.summary

    lev.sm1=lev.sm1[c('gene.nm',var)]
    lev.sm2=na.omit(lev.sm1)
    data.table::setorderv(lev.sm2,var)  ### 注意，此处改为setorderv

    lev.sm3=lev.sm2[lev.sm2[var]<p.thr,]   ###

    fdr.nm=paste0('FDR_',var)
    P=as.matrix(lev.sm3[var])

    FDR=p.adjust(P,method = mt)
    lev.sm3[fdr.nm]=FDR

    lev.fdr=lev.sm3[which(lev.sm3[fdr.nm]<=thr),]
    end.dt=merge(gene.summary,lev.fdr,by='gene.nm')
    data.table::setorderv(end.dt,fdr.nm)
    end.dt
  }

  ## 主要对比三个 fish.pval，fish.padj,change.padj
  pval.end=end.output(lev1,var='fish.pval',mt='holm',p.thr=p.thr,thr=thr)
  padj.end=end.output(lev1,var='fish.padj',mt='holm',p.thr=p.thr,thr=thr)
  change.end=end.output(lev1,var='change.padj',mt='holm',p.thr=p.thr,thr=thr)
  change.end2=end.output(lev1,var='change.pval',mt='holm',p.thr=p.thr,thr=thr)

  pv=nrow(pval.end)
  pa=nrow(padj.end)
  pc=nrow(change.end)

  rt1=pa/pv   # 通常>1   ## fish.pval和fish.padj的比值
  rt2=min(pv/pc,pv/pc)  ## fish.pval和change.padj的比值
  rt3=min(pv/pa,pa/pv)  ## fish.padj和change.padj的比值

  if (rt1 >=1.35){
    end.out=padj.end  ##  rt1>1.35,意味着pval过于严格，效果不好，应该输出padj
  }else if (pc>=pa){
    end.out=padj.end
  } else if (rt2 >0.8 & rt3 >0.8){
    if (pv>=pc ){
      end.out=pval.end
    }else{
      end.out=change.end
    }
    ########## pc和pv 输出大的,前提是change 不是最大的
  }
  colnames(end.out)[8:9]<-c('p_val','adjusted_p_val')
  # end.out

  end.out=merge(end.out,lev1$ratio,by='gene.nm')
  all.statistics=merge(lev1$test.p,lev1$gene.summary,by='gene.nm')
  # end.list=list(end.out=end.out,
  #               all.statistics=all.statistics)

  if(all.results==T){
    end=list(end.out=end.out,
             all.statistics=all.statistics)
  }
  if(all.results==F){end=end.out}
  end
}


