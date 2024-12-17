

rm(list=ls())# clear the current environment
graphics.off()# clear the current plots
setwd("C:/Users/Etudiant/Documents/thèse/1ere année/1er papier/finals/1_st paper")
source("Mixing_models_functions.R")
library(deSolve)
library(philentropy)
library("viridisLite")




#EXPERIMENT TIME
period=2000
t <- seq(0,period+1,1) # output times for the ODE (d)

#CHANGE HERE IF THE SOURCE NUMBER OR ISOTOPE NUMBER CHANGE
TEF=list(data.frame("s1"=c(0),"s2"=c(0),"s3"=c(0)),data.frame("s1"=c(0),"s2"=c(0),"s3"=c(0))) #the TEF=0 for all sources and isotope
conc=list(data.frame("s1"=c(1),"s2"=c(1),"s3"=c(1)),data.frame("s1"=c(1),"s2"=c(1),"s3"=c(1))) #the conc=1 for all sources and isotope
iter=0.1

lambda=0.02
liste_lambda=list(data.frame(Date=t,l=rep(lambda,period+2)),data.frame(Date=t,l=rep(lambda,period+2)))
#WARNING IF T IS TOO SHORT COMPARED TO LAMBDA
if((2*log(2))/lambda >period)warning(paste("period should be at least",(2*log(2))/lambda,"days"))




#GENERATING REGULAR SWITCHS in diet
pB_func <- function(t,b,a,k) {a*sin(t*b)+k} 

b=0.01#switch freq
a=0.5 #amplitude
k=0.5 #oscillation axis
combis2=pB_func(t,b,a,k)
combis3=1-combis2

# windows()
# plot(NULL,xlim=c(0,period+2),ylim=c(0,2),xlab="Time (d)",ylab="Contribution")
# points(t,rep(0,period+2),type="l",lwd=1.4,col="black",lty=2)
# points(t,combis2,type="l",lwd=1.4,col="red")
# points(t,combis3,type="l",lwd=1.4,col="green")
# legend("topright",legend=c("Source1","Source2","Source3"),cex=0.8,col=c("black","red","green"),lty=c(2,1,1))
# 

#SOURCES VALUES
#all sources constant over time
# sources_list=list(data.frame(Date=t,s1=rep(10,period+2),s2=rep(0,period+2),s3=rep(5,period+2)),
#                   data.frame(Date=t,s1=rep(5,period+2),s2=rep(0,period+2),s3=rep(10,period+2)))

#SOURCES WITH RANDOM VARIATION
#GENERATING
# s3=c(10)
# for (i in 1:(period+1)) {
#   n=rnorm(1,0,0.1)
#   s3=c(s3,s3[i]+n)
# }
# plot(t,s3)
# write.csv2(s3,"s3_N_DMM_degradated.csv")
s2_C=read.csv2("s2_C_DMM_degradated.csv")[-1][[1]]
s3_C=read.csv2("s3_C_DMM_degradated.csv")[-1][[1]]
s2_N=read.csv2("s2_N_DMM_degradated.csv")[-1][[1]]
s3_N=read.csv2("s3_N_DMM_degradated.csv")[-1][[1]]

sources_list=list(data.frame(Date=t,s1=rep(10,period+2)
                             ,s2=s2_C
                             ,s3=s3_C),
                  data.frame(Date=t,s1=rep(-5,period+2),
                             s2=s2_N,
                             s3=s3_N))


sources_approxfun=list(c(approxfun(sources_list[[1]]$Date,sources_list[[1]]$s1,method='linear',rule=2)
                         ,approxfun(sources_list[[1]]$Date,sources_list[[1]]$s2,method='linear',rule=2)
                         ,approxfun(sources_list[[1]]$Date,sources_list[[1]]$s3,method='linear',rule=2)),
                       c(approxfun(sources_list[[2]]$Date,sources_list[[2]]$s1,method='linear',rule=2)
                         ,approxfun(sources_list[[2]]$Date,sources_list[[2]]$s2,method='linear',rule=2)
                         ,approxfun(sources_list[[2]]$Date,sources_list[[2]]$s3,method='linear',rule=2)))


data_feat=data_features(sources_list )
#CONSUMER INITAL VALUE
conso_ini=c(5,10)

#GENERATING δX and δXinf
res=NULL
for (n in 1:data_feat[1]){ #data_feat[1] is the number of isotopes
  res=c(res,list(RUN_TIMdyn2(t,#time
                             state0=c(X=conso_ini[n]),#initial conditions
                             par=c(#two parameters lambda and Xinf i.e. forcing functions 
                               lambda = approxfun(liste_lambda[[n]]$Date,liste_lambda[[n]]$l,method = 'const',rule=2), 
                               Xinf=approxfun(t, Mixing(t=t,sources_approx=sources_approxfun[[n]],p=data.frame(c(0),approxfun(t,combis2)(t),approxfun(t,combis3)(t))
                                                        ,TEF=TEF[[n]],Conc=conc[[n]])))))) 
}

#VISUALISATION OF THE EXPERIMENT
#######
nb_iterations=100
lambda_T=NULL
bias_DMM=NULL
bias_SMM=NULL
bias_SMM_delta=NULL
for (i in 1:nb_iterations){
#sampling times
nb_val=10
samp_time=c(sample(1:period,nb_val))
samp_time=as.numeric(samp_time)
samp_2=order(samp_time,decreasing = FALSE)
samp_3=samp_time[samp_2]
# write.csv2(samp_3,"sampled_values",row.names = F)
# samp_time=read.csv2("sampled_values")
samp_time=data.frame(samp_3)
#extracting the consumer values at the sampling times
samp_C=data.frame(row.names = NULL)
samp_N=data.frame(row.names = NULL)
for (i in 1:nrow(samp_time)){
  samp_C=rbind(samp_C,subset(res[[1]],res[[1]]$time==samp_time[i,]))
  samp_N=rbind(samp_N,subset(res[[2]],res[[2]]$time==samp_time[i,]))
}


#VISUALISATION OF THE EXPERIMENT bis
# windows()
# close.screen(all = TRUE)   
# split.screen(1:2)
# screen(1) ; plot(NA, NA, ylim = c(-10,15), xlim = c(0,period), 
#                  xlab = c("time(d) "), ylab ="" , las=1, cex.lab=1.4,font.lab=2, cex.axis=1.4,main="")
# mtext(expression(paste(" δ"^{13},"C (‰)")),side=2,line=4,padj=1.5,cex=1.3,font=2)
# points(res[[1]]$time,sources_list[[1]]$s1,type="l",col="red",lty=2)
# points(res[[1]]$time,res[[1]]$X,col="aquamarine",type="l",lwd=2)
# points(res[[1]]$time,sources_list[[1]]$s2,type="l",col="red")
# points(res[[1]]$time,sources_list[[1]]$s3,type="l",col="red")
# points(samp_C$time,samp_C$X,type="p",pch=16,col="black")
# screen(2);plot(NA, NA, ylim = c(-10,15), xlim = c(0,period), 
#                xlab = c("time(d) "), ylab ="" , las=1, cex.lab=1.4,font.lab=2, cex.axis=1.4,main="")
# mtext(expression(paste(" δ"^{15},"N (‰)")),side=2,line=4,padj=1.5,cex=1.3,font=2)
# # points(res[[2]]$time,res[[2]]$Xinf,lty=3,type="l",col="black"
# points(res[[2]]$time,sources_list[[2]]$s1,type="l",col="red",lty=2)
# points(res[[2]]$time,res[[2]]$X,col="aquamarine",type="l",lwd=2)
# points(res[[2]]$time,sources_list[[2]]$s2,type="l",col="red")
# points(res[[2]]$time,sources_list[[2]]$s3,type="l",col="red")
# points(samp_N$time,samp_N$X,type="p",pch=16,col="black")
# legend("bottomleft",legend=c("Consumer","Source1","Source2-3","Sampled consumer")
#        ,col=c("aquamarine","red","red","black"),lty=c(1,2,1,NA),pch=c(NA,NA,NA,16),cex=1.1,
#        lwd=2, inset = 0.03, bty = "n")

#GENERATING THE MEAN REAL DIET BETWEEN TWO SAMPLINGS

mean_real_diet_combi2=NULL
for (i in 1:(nrow(samp_time)-1)){
  mean_real_diet_combi2=c(mean_real_diet_combi2,mean(combis2[samp_time[i,]:samp_time[i+1,]]))
}
#data fram of the real mean diet between the sampling times
real_mean_diet=data.frame(time=samp_time,V1=rep(0,nrow(samp_time)),V2=c(NA,mean_real_diet_combi2),
                     V3=c(NA,1-mean_real_diet_combi2))


##RUNNING THE DEGARDED DMM
#formating data
list_sampled_conso=list(samp_C,samp_N)
TEF_deg=list(data.frame("s1"=rep(0,nrow(samp_time)),"s2"=rep(0,nrow(samp_time)),"s3"=rep(0,nrow(samp_time)))
             ,data.frame("s1"=rep(0,nrow(samp_time)),"s2"=rep(0,nrow(samp_time)),"s3"=rep(0,nrow(samp_time)))) #the TEF=0 for all sources and isotope
conc_deg=list(data.frame("s1"=rep(1,nrow(samp_time)),"s2"=rep(1,nrow(samp_time)),"s3"=rep(1,nrow(samp_time)))
              ,data.frame("s1"=rep(1,nrow(samp_time)),"s2"=rep(1,nrow(samp_time)),"s3"=rep(1,nrow(samp_time)))) #the conc=1 for all sources and isotope

# #Running the DMM keeping the best solution
# res_DMM=data.frame(row.names = NULL)
# for (i in 1:(length(samp_time)-1)) {
#   period_data=Choose_data(sources_list=sources_list,consu_list=list_sampled_conso
#                           ,conc_list=conc_deg,iter=iter,obs_nb=i) 
#   res_DMM=rbind(res_DMM,DMM_func(details_consu=period_data[[3]],details_consu_t1=period_data[[5]]
#                                  ,sources_approx=sources_approxfun,
#                                  list_TEF=TEF,details_conc=period_data[[4]]
#                                  ,Combi=period_data[[1]],X=1,lambda_list=liste_lambda
#                                  ,data_feat=data_features(sources_list))
#   )
# }

#Running all models for the best solution

res_mod=all_results(sources_list=sources_list,data_feat=c(2,3,nb_val),consu_list=list_sampled_conso
                ,conc_list=conc_deg,iter=iter,list_TEF=TEF,X=1,
                sources_approx=sources_approxfun,lambda_list=liste_lambda)

SMM=res_mod[[1]]
SMM_delta=res_mod[[2]]
DMM=res_mod[[3]]

####BIAS PART

#estimating the bias on source3
bias_DMM=c(bias_DMM,abs(mean_real_diet_combi2-DMM$Var2-DMM$Var1))
bias_SMM=c(bias_SMM,abs(mean_real_diet_combi2-SMM$Var2[-1]-SMM$Var1[-1]))
bias_SMM_delta=c(bias_SMM_delta,abs(mean_real_diet_combi2-SMM_delta$Var2-SMM_delta$Var1))
#T value 
T=NULL
for(i in 1:(nrow(samp_time)-1)) {T=c(T,samp_time[i+1,]-samp_time[i,])}
lambda_T=c(lambda_T,lambda*T)
}
##

windowsFonts(
  A=windowsFont("Arial Black")
)


#plot
# windows()
# plot(lambda_T,bias_SMM[-1],xlim=c(0,lambda_T[length(lambda_T)]),ylim=c(0,1),pch=16,col="black",cex=1.4)
# points(lambda_T,bias_SMM_delta,pch=16,col="pink")
# points(lambda_T,bias_DMM,pch=16,col="red")
# windows()
plot(NULL,xlim=c(0,10),ylim=c(0,1),xlab="λT",ylab="", las=1, cex.lab=1.4,font=1)
mtext("Bias (%)",side=2,line=4,padj=1.5,cex=1.3,font=1)
 # points(lambda_T,bias_SMM,type="p",pch=16,lty=2,col="blue")
# points(lambda_T,bias_SMM_delta,,type="p",pch=16,lty=2,col="red")
points(lambda_T,bias_DMM,type="p",pch=16,cex=1.1,lwd=1.4,col="black")
legend("topright",legend=c("DMM"),col=c("black")
       ,lty=c(1),cex=1.1,lwd=2, inset = 0.03, bty = "n")

 windows()
plot(NULL,xlim=c(0,10),ylim=c(0,1),xlab="λT",ylab="", las=1, cex.lab=1.4,font=1)
mtext("Bias (%)",side=2,line=4,padj=1.5,cex=1.3,font=1)
points(lambda_T,bias_SMM,type="p",pch=16,lty=2,col="grey50")
points(lambda_T,bias_SMM_delta,,type="p",pch=16,lty=2,col="grey")
# points(lambda_T,bias_DMM,type="p",pch=16,cex=1.1,lwd=1.4,col="black")
legend("topright",legend=c("SMM_delta","SMM"),col=c("grey50","grey"),
       lty=1,cex=1.1,lwd=2, inset = 0.03, bty = "n")



data_ggplot_essai= data.frame(lambda_T=lambda_T,SMM=bias_SMM,SMM_delta=bias_SMM_delta,DMM=bias_DMM)
data_ggplot_2=subset(data_ggplot_essai,data_ggplot_essai$lambda_T<10)
library(ggplot2)

ggplot(data_ggplot_2, aes(x=lambda_T, y=SMM)) + geom_point()+geom_smooth()+
  labs(x="λT", y="Bias")+ 
  ggtitle("SMM")+
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_y_continuous( limits=c(0, 1)) 
  
ggplot(data_ggplot_2, aes(x=lambda_T, y=SMM_delta)) + geom_point()+geom_smooth()+
  labs(x="λT", y="Bias")+ 
  ggtitle("SMMΔ")+
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_y_continuous( limits=c(0, 1)) 

ggplot(data_ggplot_2, aes(x=lambda_T, y=DMM)) + geom_point()+geom_smooth()+
  labs(x="λT", y="Bias")+ 
  ggtitle("DMM")+
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_y_continuous( limits=c(0, 1)) 

