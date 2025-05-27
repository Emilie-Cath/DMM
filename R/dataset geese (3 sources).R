rm(list=ls())# clear the current environment
graphics.off()# clear the current plots
setwd("C:/Users/Etudiant/Documents/codes_2")
source("fonctions_modeles_de_melange.R")
#install.packages("deSolve")
#install.packages("ggplot2")
library("deSolve")
library("ggplot2")
library("ggridges")
library("viridisLite")

#---------------DATA---------------------------
iter=0.01
lambda_plasma=0.33
lambda_bloodcell=0.03
Dates=c(0,92,123,396,426,457,488,547,730)
Sources_Carbon=data.frame(Date=Dates,
                          Zos=rep(-11.17,9),
                          Grass=rep(-30.88,9),
                          Ulvaxentero=rep(-12.88,9))
Sources_Nitrogen=data.frame(Date=Dates,
                            Zos=rep(6.49,9),
                            Grass=rep(4.43,9),
                            Ulvaxentero=rep(10.5,9))
liste_sources=list(Sources_Carbon,Sources_Nitrogen)

# Conc_Carbon=data.frame(Zos=rep(0.36,9),Grass=rep(0.4,9),Ulvaxentero=rep(0.21,9))
# Conc_Nitrogen=data.frame(Zos=rep(0.01,9),Grass=rep(0.35,9),Ulvaxentero=rep(0.015,9))
# liste_Conc=list(Conc_Carbon,Conc_Nitrogen)


Conc_Carbon=data.frame(Zos=rep(0.35,9),Grass=rep(0.4,9),Ulvaxentero=rep(0.2,9))
Conc_Nitrogen=data.frame(Zos=rep(0.01,9),Grass=rep(0.01,9),Ulvaxentero=rep(0.01,9))
liste_Conc=list(Conc_Carbon,Conc_Nitrogen)

TDF_Carbon=data.frame(Zos=1.63,Grass=1.63,Ulvaxentero=1.63)
TDF_Nitrogen=data.frame(Zos=3.54,Grass=3.54,Ulvaxentero=3.54)
liste_TDF=list(TDF_Carbon,TDF_Nitrogen)

# 
# #PLASMA
# Conso_plasma_Carbon=data.frame(Date=Dates,
#                                Iso_C=c(-12.09,-15.79,-21,-12.42,-14.84,-20.23,-25.7,-27.32,-10.8),
#                                lambda=rep(lambda_plasma,9))
# Conso_plasma_Nitrogen=data.frame(Date=Dates,
#                                  Iso_N=c(10.32,11.77,10.6,11,12.12,10.82,8.97,8.3,9.41),
#                                  lambda=rep(lambda_plasma,9))
# liste_Conso_plasma=list(Conso_plasma_Carbon,Conso_plasma_Nitrogen)

#blood cells

Conso_bloodcell_Carbon=data.frame(Date=Dates,
                                  Iso_C=c(-12.64,-14.53,-20.07,-13.99,-12.36,-21.28,-26.33,-27.56,-12.68),
                                  lambda=rep(lambda_bloodcell,9))
Conso_bloodcell_Nitrogen=data.frame(Date=Dates,
                                    Iso_C=c(9.09,10.30,10.51,9.55,10.30,9.17,8.84,8.09,8.49),
                                    lambda=rep(lambda_bloodcell,9))
liste_Conso_bloodcell=list(Conso_bloodcell_Carbon,Conso_bloodcell_Nitrogen)

Sources_C=c(approxfun(Sources_Carbon$Date,Sources_Carbon$Zos,method='linear',rule=2),
            approxfun(Sources_Carbon$Date,Sources_Carbon$Grass,method='linear',rule=2),
            approxfun(Sources_Carbon$Date,Sources_Carbon$Ulvaxentero,method='linear',rule=2))
Sources_N=c(approxfun(Sources_Nitrogen$Date,Sources_Nitrogen$Zos,method='linear',rule=2),
            approxfun(Sources_Nitrogen$Date,Sources_Nitrogen$Grass,method='linear',rule=2),
            approxfun(Sources_Nitrogen$Date,Sources_Nitrogen$Ulvaxentero,method='linear',rule=2))
liste_approx=list(Sources_C,Sources_N)
#plasma
# liste_lambda=list(data.frame(Conso_plasma_Carbon$Date,Conso_plasma_Carbon$lambda),data.frame(Conso_plasma_Nitrogen$Date,Conso_plasma_Nitrogen$lambda))
#blood cells
# 
liste_lambda=list(data.frame(Conso_bloodcell_Carbon$Date,Conso_bloodcell_Carbon$lambda),data.frame(Conso_bloodcell_Nitrogen$Date,Conso_bloodcell_Nitrogen$lambda))

#------------running the thre models---------------
#fait tourner
res=all_results(sources_list=liste_sources,data_feat=data_features(liste_sources),consu_list=liste_Conso_bloodcell
                ,conc_list=liste_Conc,iter=iter,list_TEF=liste_TDF,X=50,sources_approx=liste_approx,lambda_list=liste_lambda)

# res_1=all_results(sources_list=liste_sources,data_feat=data_features(liste_sources),consu_list=liste_Conso_bloodcell
#                   ,conc_list=liste_Conc,iter=iter,list_TEF=liste_TDF,X=1,sources_approx=liste_approx,lambda_list=liste_lambda)


# write.csv2(res[[1]],file="SMM_all_sol_ge3.csv")
# write.csv2(res[[2]],file="SMM_delta_50_ge3.csv")
# write.csv2(res[[3]],file="DMM_all_sol_ge3.csv")
SMM_50=read.csv2("SMM_50_ge3.csv")[-1]
SMM_delta_50=read.csv2("SMM_delta_50_ge3.csv")[-1]
DMM_50=read.csv2("DMM_50_ge3.csv")[-1]
res_50=c(res[[1]],res[[2]],res[[3]])

# res_50=c(list(SMM_50),list(SMM_delta_50),list(DMM_50))

# write.csv2(res_1[[1]],file="SMM_1_ge3.csv")
# write.csv2(res_1[[2]],file="SMM_delta_1_ge3.csv")
# write.csv2(res_1[[3]],file="DMM_1_ge3.csv")
# SMM_1=read.csv2("SMM_1_ge3.csv")[-1]
# SMM_delta_1=read.csv2("SMM_delta_1_ge3.csv")[-1]
# DMM_1=read.csv2("DMM_1_ge3.csv")[-1]
# res_1=c(list(SMM_1),list(SMM_delta_1),list(DMM_1))
#------------GRAPHS_bloodcells---------------

#biplot visuals
sign_biplot(consu_list=liste_Conso_bloodcell,sources_list=liste_sources,list_TEF=liste_TDF
            ,title="",range_C=c(-30,-8),range_N=c(7,15),source_names=c("Zos","grass","UlvaxEntero"))
poly_sources_evol(consu_list=liste_Conso_bloodcell,sources_list=liste_sources,list_TEF=liste_TDF,sources_name=c("Zos","grass","UlvaxEntero"),range_C=c(-30,-8),range_N=c(7,15))



#graphe_traj
traj_graph(sources_list=liste_sources,data_feat=data_features(liste_sources),sources_approx=liste_approx,consu_list=liste_Conso_bloodcell,iso_num=1
           ,lambda_list=liste_lambda,conc_list=liste_Conc,list_TEF=liste_TDF,iter=iter,X=50,time_range=c(0,730),range_y=c(-30,-8),title="C"
           ,y_label=expression(paste(" δ"^{13},"C (‰)")))
traj_graph(sources_list=liste_sources,data_feat=data_features(liste_sources),sources_approx=liste_approx,consu_list=liste_Conso_bloodcell,iso_num=2
           ,lambda_list=liste_lambda,conc_list=liste_Conc,list_TEF=liste_TDF,iter=iter,X=50,time_range=c(0,730),range_y=c(7,15),title="N"
           ,y_label=expression(paste(" δ"^{15},"N (‰)")))

#ridgeline
ridgeline_chart(res_model=SMM_50,model_name="SMM",sources_names=c("Zos","grass","UlvaxEntero"),X=50)
ridgeline_chart(res_model=SMM_delta_50,model_name="SMM_delta",sources_names=c("Zos","grass","UlvaxEntero"),X=50)
ridgeline_chart(res_model=DMM_50,model_name="DMM",sources_names=c("Zos","grass","UlvaxEntero"),X=50)

#distance 50 sol
dist_graph(obs_time=Dates,res_model=res_50,range_y=c(0,0.25))

#distance 1 sol
dist_graph(obs_time=Dates,res_model=res_1,range_y=c(0,0.25))
#source contrib
sources_contrib(mod_results = res_1,obs_time=Dates,title="Zoostera",source_num=1,"topleft")
sources_contrib(mod_results = res_1,obs_time=Dates,title="Grass",source_num=2,"topleft")
sources_contrib(mod_results = res_1,obs_time=Dates,title="Ulva",source_num=3,"topleft")
# sources_contrib(mod_results = res_1,obs_time=Dates,title="Entero",source_num=4,"topleft")

#contrib over time
contrib_time(res_model=res[[1]],model_name="SMM",date=Dates,data_feat=data_features(liste_sources)
             ,X=50,sources_name=c("Zos","grass","UlvaxEntero"),where_legend="topright")
contrib_time(res_model=SMM_delta_50,model_name="SMM_delta",date=Dates,data_feat=data_features(liste_sources)
             ,X=50,sources_name=c("Zos","grass","UlvaxEntero"),where_legend="topright")
contrib_time(res_model=res[[3]],model_name="DMM",date=Dates,data_feat=data_features(liste_sources)
             ,X=50,sources_name=c("Zos","grass","UlvaxEntero"),where_legend="topright")

#------------GRAPHS_PLASMA---------------
res=all_results(sources_list=liste_sources,data_feat=data_features(liste_sources),consu_list=liste_Conso_plasma
                ,conc_list=liste_Conc,iter=iter,list_TEF=liste_TDF,X=50,sources_approx=liste_approx,lambda_list=liste_lambda)
ridgeline_chart(res_model=res[[1]],model_name="SMM",sources_names=c("Zos","grass","UlvaxEntero"),X=50)
ridgeline_chart(res_model=res[[3]],model_name="DMM",sources_names=c("Zos","grass","UlvaxEntero"),X=50)
contrib_time(res_model=res[[1]],model_name="SMM",date=Dates,data_feat=data_features(liste_sources)
             ,X=50,sources_name=c("Zos","grass","UlvaxEntero"),where_legend="topright")
contrib_time(res_model=res[[3]],model_name="DMM",date=Dates,data_feat=data_features(liste_sources)
             ,X=50,sources_name=c("Zos","grass","UlvaxEntero"),where_legend="topright")

#------------BIAS ESTIMATION---------------
bias_SMM=data.frame(row.names=NULL)
for (i in 1:nrow(DMM_1)) {bias_SMM=rbind(bias_SMM,abs(DMM_1[i,3]+DMM_1[i,4]
                                       -SMM_1[i+1,3]-SMM_1[i+1,4]))}
colnames(bias_SMM)<-c("bias")
lambda_T <- c(3.285,1.107,9.750,1.071,1.107,1.107,2.107,6.535)
bias_SMM=cbind(Time=DMM_1[[1]],bias_SMM,lambda_T)

# plot(bias_SMM$lambda_T,bias_SMM$bias)

#------------BIAS ESTIMATION---------------
#We keep the 1% best solution 
# SMM_50=read.csv2("SMM_50_ge3.csv")[-1]
# DMM_50=read.csv2("DMM_50_ge3.csv")[-1]
SMM_50=res[[1]]
DMM_50=res[[3]]

moments_DMM=unique(DMM_50$time)

#first time for SMM
t_1=subset(SMM_50,SMM_50$time==0)
v1=quantile(t_1$Var1,probs = c(0.25,0.5,0.75))
v2=quantile(t_1$Var2,probs = c(0.25,0.5,0.75))
v3=quantile(t_1$Var3,probs = c(0.25,0.5,0.75))
df_stats_SMM_0=data.frame(time=0,v1_Q1=v1[1],v1_M=v1[2],v1_Q3=v1[3],
                      v2_Q1=v2[1],v2_M=v2[2],v2_Q3=v2[3],
                      v3_Q1=v3[1],v3_M=v3[2],v3_Q3=v3[3])

#loop 
df_stats_SMM_others=data.frame(row.names = F)
df_stats_DMM=data.frame(row.names = F)
for (i in 1:length(moments_DMM)){ #for every time
  val=moments_DMM[i]
  subs_SMM=subset(SMM_50,SMM_50$time==val)
  subs_DMM=subset(DMM_50,time==val)
  #SMM
  sv1=quantile(subs_SMM$Var1,probs = c(0.25,0.5,0.75))
  sv2=quantile(subs_SMM$Var2,probs = c(0.25,0.5,0.75))
  sv3=quantile(subs_SMM$Var3,probs = c(0.25,0.5,0.75))
  df_stats_SMM_others=rbind(df_stats_SMM_others
                            ,data.frame(time=val
                                        ,v1_Q1=sv1[1],v1_M=sv1[2],v1_Q3=sv1[3]
                                        ,v2_Q1=sv2[1],v2_M=sv2[2],v2_Q3=sv2[3]
                                        ,v3_Q1=sv3[1],v3_M=sv3[2],v3_Q3=sv3[3]))
  #DMM
  dv1=quantile(subs_DMM$Var1,probs = c(0.25,0.5,0.75))
  dv2=quantile(subs_DMM$Var2,probs = c(0.25,0.5,0.75))
  dv3=quantile(subs_DMM$Var3,probs = c(0.25,0.5,0.75))
  df_stats_DMM=rbind(df_stats_DMM
                            ,data.frame(time=val
                                        ,v1_Q1=dv1[1],v1_M=dv1[2],v1_Q3=dv1[3]
                                        ,v2_Q1=dv2[1],v2_M=dv2[2],v2_Q3=dv2[3]
                                        ,v3_Q1=dv3[1],v3_M=dv3[2],v3_Q3=dv3[3]))
}

df_stats_SMM=rbind(df_stats_SMM_0,df_stats_SMM_others)

#GRAPH DMM
Q1_col=c(2,5,8)
M_col=c(3,6,9)
Q3_col=c(4,7,10)

nb_sources=3
color_main=c(rgb(0,0.45,0.7),rgb(0.83,0.36,0), #up to 7 sources
             rgb(0,0.6,0.45),rgb(0.94,0.89,0.25)
             ,rgb(0.8,0.47,0.65),rgb(0.9,0.6,0),
             rgb(0.33,0.7,0.9))
color_poly=c(rgb(0,0.45,0.7,0.2),rgb(0.83,0.36,0,0.2),rgb(0,0.6,0.45,0.2),rgb(0.94,0.89,0.25,0.2)
             ,rgb(0.8,0.47,0.65,0.2),rgb(0.9,0.6,0,0.2),rgb(0.33,0.7,0.9,0.2) )
#Graph DMM
windows()
par(xpd=TRUE, mar=c(8,4,4,3))
plot(NULL,xlim=c(0,730),ylim=c(0,1),xlab="time(day)",ylab="",main="DMM",
     , las=1, cex.lab=1.3,font=1)
mtext("Contribution",side=2,line=4,padj=1.5,cex=1.3,font=1)
for (i in 1:nb_sources) {
  
  polygon(x=c(df_stats_DMM$time,rev(df_stats_DMM$time))
                                 ,y=c(df_stats_DMM[[Q3_col[i]]],rev(df_stats_DMM[[Q1_col[i]]]))
                                 ,col=color_poly[i],border=NA)}

for (i in 1:nb_sources){
  points(df_stats_DMM$time,df_stats_DMM[[M_col[i]]],type="l"
                               ,col=color_main[i],
                               pch=16,lwd=2)}
legend("topright",legend=c("Zos","Grass","UlvaxEntero"),cex=1,col=c(color_main[1:nb_sources]),lty="solid",box.lty=0)


#GRAPH SMM
windows()
par(xpd=TRUE, mar=c(8,4,4,3))
plot(NULL,xlim=c(0,730),ylim=c(0,1),xlab="time(day)",ylab="",main="SMM",
     , las=1, cex.lab=1.3,font=1)
mtext("Contribution",side=2,line=4,padj=1.5,cex=1.3,font=1)
for (i in 1:nb_sources) {
  
  polygon(x=c(df_stats_SMM$time,rev(df_stats_SMM$time))
          ,y=c(df_stats_SMM[[Q3_col[i]]],rev(df_stats_SMM[[Q1_col[i]]]))
          ,col=color_poly[i],border=NA)}

for (i in 1:nb_sources){
  points(df_stats_SMM$time,df_stats_SMM[[M_col[i]]],type="l"
         ,col=color_main[i],
         pch=16,lwd=2)}
# legend("topright",legend=c("Zos","Grass","UlvaxEntero"),cex=1,col=c(color_main[1:nb_sources]),lty="solid",box.lty=0)


#BIAS
Bias=data.frame(row.names = F)
lambda_T <- c(3.285,1.107,9.750,1.071,1.107,1.107,2.107,6.535)
for (i in 1:nrow(df_stats_DMM)){
  B1=(abs(df_stats_DMM$v1_M[i]-df_stats_SMM_others$v1_M[i])+abs(df_stats_DMM$v2_M[i]-df_stats_SMM_others$v2_M[i])+abs(df_stats_DMM$v3_M[i]-df_stats_SMM_others$v3_M[i]))*50
  Bias=rbind(Bias, data.frame(time=df_stats_DMM$time[i],lambda_T=lambda_T[i], B=B1))
}


plot(NULL,xlim=c(0,10),ylim=c(0,100),xlab="λT",ylab="β (%)"
     , las=1, cex.lab=1.4,font.lab=2, cex.axis=1.4,main="")
points(Bias$lambda_T,Bias$B,type="p",pch=16)
text(x=Bias$lambda_T+0.4,y=Bias$B+1,labels = Bias$time,cex=0.5,font=2,col="black")



###

dev.off()
#DMM
par(fig=c(0,0.4,0,1), new=TRUE,mar = c(4, 3.5, 4.1, 0.2))
plot(NULL,xlim=c(0,730),ylim=c(0,1),xlab="",ylab=""
     , las=1, cex.lab=1.4,font.lab=2, cex.axis=1.4,main="")
mtext("Contribution",side=2,line=4,padj=1.9,cex=1.2,font=1)
mtext("Time (d)",side=1,line=1,padj=1.5,at=10,cex=1.2,font=1)
for (i in 1:nb_sources) {
  
  polygon(x=c(df_stats_DMM$time,rev(df_stats_DMM$time))
          ,y=c(df_stats_DMM[[Q3_col[i]]],rev(df_stats_DMM[[Q1_col[i]]]))
          ,col=color_poly[i],border=NA)}

for (i in 1:nb_sources){
  points(df_stats_DMM$time,df_stats_DMM[[M_col[i]]],type="l"
         ,col=color_main[i],
         pch=16,lwd=2)}
# legend("topright",legend=c("Zos","Grass","UlvaxEntero"),cex=1,col=c(color_main[1:nb_sources]),lty="solid",box.lty=0)
#SMM
par(fig=c(0.4,0.8,0,1), new=TRUE,mar = c(4, 3.5, 4.1, 0.2))
plot(NULL,xlim=c(0,730),ylim=c(0,1),xlab="",ylab=""
     , las=1, cex.lab=1.4,font.lab=2, cex.axis=1.4,main="")
mtext("Contribution",side=2,line=4,padj=1.9,cex=1.2,font=1)
mtext("Time (d)",side=1,line=1,padj=1.5,at=10,cex=1.2,font=1)
for (i in 1:nb_sources) {
  
  polygon(x=c(df_stats_SMM$time,rev(df_stats_SMM$time))
          ,y=c(df_stats_SMM[[Q3_col[i]]],rev(df_stats_SMM[[Q1_col[i]]]))
          ,col=color_poly[i],border=NA)}

for (i in 1:nb_sources){
  points(df_stats_SMM$time,df_stats_SMM[[M_col[i]]],type="l"
         ,col=color_main[i],
         pch=16,lwd=2)}
# legend("topright",legend=c("Zos","Grass","UlvaxEntero"),cex=1,col=c(color_main[1:nb_sources]),lty="solid",box.lty=0)
#BIAS
# par(fig=c(0.8,1,0,1), new=TRUE)
plot(NULL,xlim=c(0,10),ylim=c(0,100),xlab="",ylab=""
     , las=1, cex.lab=1.4,font.lab=2, cex.axis=1.4,main="")
mtext("β (%)",side=2,line=4,padj=2.2,cex=1.2,font=1)
mtext("λT",side=1,line=1,padj=1.5,at=5,cex=1.2,font=1)
points(Bias$lambda_T,Bias$B,type="p",pch=16)
text(x=Bias$lambda_T-0.2,y=Bias$B+3,labels = Bias$time,cex=0.5,font=2,col="black")


plot(c(0,2),c(0,1),type = 'n', axes = F,xlab = '', ylab = '', main = '')
legend("center",legend=c("Zos","Grass","UlvaxEntero","Sampling date"),cex=1,col=c(color_main[1:nb_sources],"black","black"),lty=c(1,1,1,NA),pch=c(NA,NA,NA,16),box.lty=0,lwd=2)
