
eps<-10^(-5);
is.odd <- function(x) x %% 2 != 0
#' List of functions
#' freq.less function
#'
#' This function finds the sum of counts that the x-sample observations is greater than or less than the ones from the y-sample.
#' @param x,y x and y are numerical vectors of different subsamples. The length of the two vectors can vary.
#' @details When there is a tie between any pair of observations, 0.5 is added to the count. Missing value is allowed. Missing value is only added to the calculation when it is compared with another missing value from the other subsample.
#' @return Two values are returned: less.count and more.count. The first one is the total count that the observations in x-sample is less than the ones from the y-sample, and the second output is the total count that the observations in x-sample is more than the ones from the y-sample. When there is a tie, 0.5 is added to the count, instead of 1 or 0.
#' @examples
#' freq.less(x=c(1,2,4,9,0,0,NA),y=c(1,4,9,NA))
#' @export
freq.less<-function(x,y)
{
  if(all(is.na(x))|all(is.na(y))) return(c(NA,NA));

  n.fx<-length(x);n.fy<-length(y); n.all<-n.fx+n.fy;

  x[abs(x)<eps]<-0;y[abs(y)<eps]<-0;

  less.count<-sum(outer(x,y,'<')+outer(x,y,'==')*.5,na.rm=TRUE)
  more.count<-sum(outer(x,y,'>')+outer(x,y,'==')*.5,na.rm=TRUE)

  return(c(less.count, more.count))
}


#' multi.freq function
#'
#' This function find trend in a sample by comparing neighboring subsamples. The subsamples are stored in a list in R.
#' @param fsam a list in R. The order of the vectors in the list follows the order of the subsamples.
#' @details The first vector of data in the list will be compared with the second vector in the list by using function freq.less. Then the second vector will be compared with the 3rd vector if there is one. The statistics collected are based on computing: \deqn{\frac{1}{n_ln_{l+1}}\sum_{i=1}^{n_l}\sum_{j=1}^{n_{l+1}}1(x_{li}<x_{(l+1)j})}
#' @return count.vec it is a collection of a sequence less.count, more.count based on freq.less function.
#' @references Wang, Y., Stapleton, A. E., & Chen, C. (2018). Two-sample nonparametric stochastic order inference with an application in plant physiology. Journal of Statistical Computation and Simulation, 88(14), 2668-2683.
#' @examples
#' x1=c(1,2,4,9,0,0,NA);x2=c(1,4,9,NA);x3=c(2,5,10);
#' sam=list(x1,x2,x3); #
#' multi.freq(sam);
#' @export

multi.freq<-function(fsam)
{
  n.col<-length(fsam);
  count.vec<-NULL;
  for(fi in 2:n.col)
  {temp.f<-freq.less(fsam[[fi-1]],fsam[[fi]])
  count.vec<-c(count.vec,temp.f);}
  return(count.vec)
}

#' simu.ustat.pattern function
#'
#' This function create two independent subsamples of various subsample sizes, with a given probability vector.
#' @param mean.prob.vec a vector of length 2. Its first element represents the probability that a random observation from one subsample is less than the the one from another subsample..
#' @param effn.subs a vector contains two subsample sizes.
#' @param n.rep the total number of repetition.
#' @importFrom stats "rnorm"
#' @importFrom stats "qnorm"
#' @details each subsample is generated from a normal distribution, with an average generated from the mean.prob.vec.
#' @return simu.tab a list of length n.rep. Each element of the list is a 2 by 2 matrix, showing the comparison results from function multi.freq.
#' @references Wang, Y., Stapleton, A. E., & Chen, C. (2018). Two-sample nonparametric stochastic order inference with an application in plant physiology. Journal of Statistical Computation and Simulation, 88(14), 2668-2683.
#' @export
#' @examples
#' simu.ustat.pattern(c(0.8,0.2),c(5,8),n.rep=100)

simu.ustat.pattern<-function(mean.prob.vec, effn.subs,n.rep=10^2)
{
  mean.prob.vec<-replace(mean.prob.vec, mean.prob.vec<eps, eps);
  mean.prob.vec<-replace(mean.prob.vec, mean.prob.vec>1-eps, 1-eps);
  n.col<-length(mean.prob.vec);
  h<-cumsum(2^.5*qnorm(mean.prob.vec[1]));

  simu.norm.sam1.all<-matrix(rnorm(effn.subs[1]*n.rep, mean=0),nrow=n.rep);
  simu.norm.sam2.all<-matrix(rnorm(effn.subs[2]*n.rep, mean=h),nrow=n.rep);
  simu.tab<-lapply(c(1:n.rep), function(x) multi.freq(list(simu.norm.sam1.all[x,],simu.norm.sam2.all[x,])))
  return(simu.tab)

}


#' chi.stat function
#'
#' This function calculates the $M$ statistics value as defined in the reference paper.
#' @param ftab it is a matrix with dimension 2 by \eqn{K}.
#' @details The \eqn{M} statistics is defined as: \deqn{M=\sum_{l=1}^{K}\left(\frac{(O_{x,l}-E_{x,l})^2}{E_{x,l}}+\frac{(O_{x,l}-E_{x,l})^2}{\left(n_ln_{l+1}-E_{x,l}\right)}\right)+\sum_{l=1}^{K}\left(\frac{(O_{y,l}-E_{y,l})^2}{E_{y,l}}+\frac{(O_{y,l}-E_{y,l})^2}{\left(m_lm_{l+1}-E_{y,l}\right)}\right).}
#' @return chi.val, a chisuqre type of statistics value
#' @references Wang, Y., Stapleton, A. E., & Chen, C. (2018). Two-sample nonparametric stochastic order inference with an application in plant physiology. Journal of Statistical Computation and Simulation, 88(14), 2668-2683.
#' @export
#' @examples
#' chi.stat(ftab=rbind(c(20,10,20),c(15,15,20)))

chi.stat<-function(ftab)
{
  tot<-sum(ftab);
  expv<-outer(rowSums(ftab)/tot,  colSums(ftab)/tot, '*')*tot;
  signal<-(colMeans(expv)<eps)*(1:dim(ftab)[2]);
  indx<-setdiff(signal,0);
  ftemp<-((ftab-expv)^2/expv)
  chi.val<-ifelse(length(indx)==0, sum(ftemp), sum(ftemp[,-indx]));
  return(chi.val);
}


#' gen.decision function
#' This function compares frequency comparison counts of subsamples from two independent samples, and calculates a simulated p-value with a novel bootstrap method proposed in the reference paper.
#' @param est.prob a matrix of two rows, with each row represents the the sequential comparison results of subsamples from a sample.
#' @param effn.subsam1 the subsample sizes from sample 1.
#' @param effn.subsam2 the subsample sizes from sample 2.
#' @param fn.rep the total number of replications.
#' @param alpha the size of type I error.
#' @details The dimensions of est.prob, effn.subsam1 and effn.subsam2 need to match. For example, the first two entries of the first two rows from est.prob are pf comparison results from subsample1 and subsample2 of sample1. Thus the sum of the two entries is the product of the two subsample sizes.
#' @return critical.value the critical value of the test based on the alpha level provided
#' @return chi-stat the chisqure type test statistics value from the sample provided.
#' @return pvalue the simulated p-value.
#' @references Wang, Y., Stapleton, A. E., & Chen, C. (2018). Two-sample nonparametric stochastic order inference with an application in plant physiology. Journal of Statistical Computation and Simulation, 88(14), 2668-2683.
#' @importFrom magrittr "%>%"
#' @export
#' @examples
#' freq.mat<-rbind(c(20,5,10,15,20,5),c(15,10,15,10,20,5));
#' n.sam1<-rep(5,4);n.sam2<-rep(5,4); n.rep=1000;
#' gen.decision(freq.mat,n.sam1,n.sam2,n.rep);
#' ### This command will replicate the first p-value in Table 4 of the reference paper.
#' freq.mat<-rbind(c(40,10,20,30,40,10),c(30,20,30,20,40,10));
#' n.sam1<-c(5,10,5,10);n.sam2<-c(10,5,10,5); n.rep=1000;
#' gen.decision(freq.mat,n.sam1,n.sam2,n.rep)
#' ### This command will replicate the second p-value in Table 4 of the reference paper.

gen.decision<-function(est.prob, effn.subsam1,effn.subsam2,fn.rep=10^3,alpha=.05)
{
  n.col=length(effn.subsam1);
  base.n=rbind(rep(effn.subsam1[-n.col]*effn.subsam1[-1],rep(2,n.col-1)),
               rep(effn.subsam2[-n.col]*effn.subsam2[-1],rep(2,n.col-1)));

  mean.prob<-colSums(est.prob)/colSums(base.n);
  id.na<-((mean.prob%>%is.na)*c(1:length(mean.prob)))%>%setdiff(0);
  if(length(id.na)==0) id.na=NA;
  id.keep<-setdiff(c(1:length(mean.prob)),id.na)


  mean.prob.simu<-matrix(mean.prob[id.keep],ncol=2,byrow=TRUE);
  id.effb<-unique(round(id.keep/2+eps));



  effn.sub.simu1<-cbind(effn.subsam1[-n.col],effn.subsam1[-1])[id.effb,];
  effn.sub.simu2<-cbind(effn.subsam2[-n.col],effn.subsam2[-1])[id.effb,];

  num.fold<-dim(mean.prob.simu)[1];
  if(num.fold==1) { simu.sam1= simu.ustat.pattern(mean.prob.simu[1,], effn.subs=effn.sub.simu1,n.rep=fn.rep)
  simu.sam2= simu.ustat.pattern(mean.prob.simu[1,], effn.subs=effn.sub.simu2,n.rep=fn.rep)
  }   else {simu.sam1=list(rep(NULL,fn.rep));
  for(i in 1:num.fold) simu.sam1=mapply(c, simu.sam1, simu.ustat.pattern(mean.prob.simu[i,], effn.subs=effn.sub.simu1[i,],n.rep=fn.rep), SIMPLIFY=FALSE)
  simu.sam2=list(rep(NULL,fn.rep));
  for(i in 1:num.fold) simu.sam2=mapply(c, simu.sam2, simu.ustat.pattern(mean.prob.simu[i,], effn.subs=effn.sub.simu2[i,],n.rep=fn.rep), SIMPLIFY=FALSE)
  }

  simu.tab.list<-mapply(rbind,simu.sam1,simu.sam2,SIMPLIFY=FALSE)

  aa.temp2<-unlist(lapply(simu.tab.list, chi.stat));

  f.cri<-sort(aa.temp2)[(1-alpha)*fn.rep];
  f.stat<-chi.stat(est.prob[,id.keep]);
  pval<-sum(aa.temp2>f.stat)/fn.rep;
  result<-c(f.cri,f.stat,pval); names(result)<-c('critical.value', 'chi-stat','pvalue')
  return(result)

}


#' pow.ana.gen.decision function
#'
#' This function evaluates the type I error of the proposed test.
#' @param mean.prob1, the probability that observations of a subsample is less than the ones from another subsample, in sample #1.
#' @param mean.prob2, the probability that observations of a subsample is less than the ones from another subsample, in sample #2.
#' @param effn.subsam1 the subsample sizes from sample 1.
#' @param effn.subsam2 the subsample sizes from sample 2.
#' @param rseed a random seed.
#' @param boot.rep the number of repetitions needed to calculated simulated p-value,
#' @param N.rep the total number of bootstrap repetitions needed for calculating type I errors.
#' @param alpha.level the type I error level that will be assessed.
#' @return the simulated type I error.
#' @references Wang, Y., Stapleton, A. E., & Chen, C. (2018). Two-sample nonparametric stochastic order inference with an application in plant physiology. Journal of Statistical Computation and Simulation, 88(14), 2668-2683.
#' @importFrom stats "sd"
#' @importFrom magrittr "%>%"
#' @export
#' @examples
#' prob.vec<-c(.4,.2,.3,.6);
#' sub.sizes1<-c(2,4,3,5,3);sub.sizes2<-c(6,3,2,4,2)
#' pow.ana.gen.decision(prob.vec,prob.vec,sub.sizes1, sub.sizes1)
#' pow.ana.gen.decision(prob.vec,prob.vec,sub.sizes1, sub.sizes1,alpha.level=0.1)

pow.ana.gen.decision<-function(mean.prob1,mean.prob2,effn.subsam1,effn.subsam2,N.rep=10^1, boot.rep=10^1,rseed=1234,alpha.level=0.05)
{
  set.seed(rseed);
  n.col=length(mean.prob1);

  simu.sam1=list(rep(NULL,N.rep));mean.prob.simu=cbind(mean.prob1,1-mean.prob1)
  for(i in 1:n.col) simu.sam1=mapply(c, simu.sam1, simu.ustat.pattern(mean.prob.simu[i,], effn.subs=c(effn.subsam1[i],effn.subsam1[i+1]),n.rep=N.rep), SIMPLIFY=FALSE)
  simu.sam2=list(rep(NULL,N.rep));mean.prob.simu=cbind(mean.prob2,1-mean.prob2)
  for(i in 1:n.col) simu.sam2=mapply(c, simu.sam2, simu.ustat.pattern(mean.prob.simu[i,], effn.subs=c(effn.subsam2[i],effn.subsam2[i+1]),n.rep=N.rep), SIMPLIFY=FALSE)

  simu.tab.list<-mapply(rbind,simu.sam1,simu.sam2,SIMPLIFY=FALSE)
  dec<-lapply(c(1:N.rep),  function(x) {
    gen.decision(matrix(unlist(simu.tab.list[x]),nrow=2), effn.subsam1=effn.subsam1,effn.subsam2=effn.subsam2,fn.rep=boot.rep,alpha=.05);})
  out<-matrix(unlist(dec),nrow=N.rep, byrow=TRUE);
  return(sum(out[,3]<alpha.level)/N.rep);
}



#' sub.test function
#' This function calculates the simulated p-value of comparing the trend in subsamples from two independent samples.
#' @param sam1, the first sample.
#' @param sam2, the second sample
#' @param fn.rep2 the total number of bootstrap repetitions needed for calculating the simulated p-value.
#' @return critical.value the critical value of the test based on the alpha level provided
#' @return chi-stat the chisqure type test statistics value from the sample provided.
#' @return pvalue the simulated p-value.
#' @references Wang, Y., Stapleton, A. E., & Chen, C. (2018). Two-sample nonparametric stochastic order inference with an application in plant physiology. Journal of Statistical Computation and Simulation, 88(14), 2668-2683.
#' @export
#' @examples
#' attach(seedwt.multi.subsample)
#' Lev.TN<-levels(TreatmentName);
#' Lev.Line<-levels(Line);
#' n<-dim(seedwt.multi.subsample)[1];
#' level.show=c(1:8);fn.rep3=10^2;
#' line.name<-Lev.Line[1]; t1.name<-Lev.TN[1];t2.name<-Lev.TN[3];
#' ### To compare the GA treatment and the PACGA treatment from line B73
#' par(mfrow=c(1,2))
#' idx<-subset((TreatmentName==t1.name)*(Line==line.name)*(1:n),Env %in% level.show)
#' idx2<-subset((TreatmentName==t2.name)*(Line==line.name)*(1:n),Env %in% level.show)

#' boxplot(seedwt[idx]~Env[idx],xlab="ENV levels",ylab=paste('seedwt from',t1.name),
#'          ylim=c(0,12),cex.lab=1.5,cex.axis=1.8);
#' boxplot(seedwt[idx2]~Env[idx2], xlab="ENV levels",ylab=paste('seedwt from',t2.name),
#'          cex.lab=1.5,cex.axis=1.8);
#' mtext( paste ("Line Name:",line.name), side = 3,outer = TRUE, cex = 2.2,line = -3)
#'temp.sw1<-seedwt[idx];lab<-Env[idx]; uni.lab<-unique(lab)
#'sam.1<-lapply(1:length(uni.lab), function(x) temp.sw1[lab==uni.lab[x]])
#'temp.sw2<-seedwt[idx2];lab2<-Env[idx2]; uni.lab2<-unique(lab2)
#'sam.2<-lapply(1:length(uni.lab2), function(x) temp.sw2[lab2==uni.lab2[x]])
#'print(paste("working with line ",line.name,'and treatment',t1.name ,'vs',t2.name ))
#'resu<-sub.test(sam.1,sam.2,fn.rep2=fn.rep3);
#'## This will show a similar result as the first experiment of section 5 in the paper.

sub.test<-function(sam1,sam2,fn.rep2)## sam1 and sam2 are two lists ###
{
  est.prob.1<-multi.freq(sam1); est.prob.2<-multi.freq(sam2);
  est.prob<-rbind(est.prob.1,est.prob.2);
  effn.subsam1=unlist(lapply(c(1:length(sam1)),function(x) length(setdiff(sam1[[x]],NA))))
  effn.subsam2=unlist(lapply(c(1:length(sam2)),function(x) length(setdiff(sam2[[x]],NA))))
  result<- gen.decision(est.prob,effn.subsam1=effn.subsam1,effn.subsam2=effn.subsam2,fn.rep=fn.rep2,alpha=.05);
  return(result)

}


#' seedwt.multi.subsample dataset
#' @details multiple maize inbreds were exposed to all combinations of the following stressors: drought, nitrogen, and density stress. Plants were grown in an experimental plot divided into eight sections, and each of the sections received a combination of between zero and three of the stresses previously mentioned, so that all possible stress combinations were included. More details about the experiment can be found in the references
#' @importFrom usethis "use_data"
#' @references Wang, Y., Stapleton, A. E., & Chen, C. (2018). Two-sample nonparametric stochastic order inference with an application in plant physiology. Journal of Statistical Computation and Simulation, 88(14), 2668-2683.
#' @references Stutts, L., Wang, Y., & Stapleton, A. E. (2018). Plant growth regulators ameliorate or exacerbate abiotic, biotic and combined stress interaction effects on Zea mays kernel weight with inbred-specific patterns. Environmental and experimental botany, 147, 179-188.
#' @export

seedwt.multi.subsample<-read.csv("inst/extdata/NamedLevelsseedweightrans_periods_zerosV2.csv")

usethis::use_data(seedwt.multi.subsample, overwrite = T);


