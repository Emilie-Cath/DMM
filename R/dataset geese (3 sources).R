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
                                  Iso_C=c(-27.63,-26.33,-21.26,-20.07,-14.50,-14.02,-12.65,-12.54,-12.32),
                                  lambda=rep(lambda_bloodcell,9))
Conso_bloodcell_Nitrogen=data.frame(Date=Dates,
                                    Iso_C=c(8.09,8.84,9.18,10.47,10.28,9.54,8.48,9.04,10.26),
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

#------------fait tourner---------------
#fait tourner
res=all_results(sources_list=liste_sources,data_feat=data_features(liste_sources),consu_list=liste_Conso_bloodcell
                ,conc_list=liste_Conc,iter=iter,list_TEF=liste_TDF,X=50,sources_approx=liste_approx,lambda_list=liste_lambda)

res_1=all_results(sources_list=liste_sources,data_feat=data_features(liste_sources),consu_list=liste_Conso_bloodcell
                  ,conc_list=liste_Conc,iter=iter,list_TEF=liste_TDF,X=1,sources_approx=liste_approx,lambda_list=liste_lambda)

# write.csv2(res[[1]],file="SMM_50_ge3.csv")
# write.csv2(res[[2]],file="SMM_delta_50_ge3.csv")
# write.csv2(res[[3]],file="DMM_50_ge3.csv")
SMM_50=read.csv2("SMM_50_ge3.csv")[-1]
SMM_delta_50=read.csv2("SMM_delta_50_ge3.csv")[-1]
DMM_50=read.csv2("DMM_50_ge3.csv")[-1]
res_50=c(list(SMM_50),list(SMM_delta_50),list(DMM_50))

# write.csv2(res_1[[1]],file="SMM_1_ge3.csv")
# write.csv2(res_1[[2]],file="SMM_delta_1_ge3.csv")
# write.csv2(res_1[[3]],file="DMM_1_ge3.csv")
SMM_1=read.csv2("SMM_1_ge3.csv")[-1]
SMM_delta_1=read.csv2("SMM_delta_1_ge3.csv")[-1]
DMM_1=read.csv2("DMM_1_ge3.csv")[-1]
res_1=c(list(SMM_1),list(SMM_delta_1),list(DMM_1))
#------------GRAPHs---------------

#biplot visuals
sign_biplot(consu_list=liste_Conso_bloodcell,sources_list=liste_sources,list_TEF=liste_TDF
            ,title="",range_C=c(-30,-8),range_N=c(7,15),source_names=c("Zos","grass","UlvaxEntero"))
poly_sources_evol(consu_list=liste_Conso_bloodcell,sources_list=liste_sources,list_TEF=liste_TDF,sources_name=c("Zos","grass","UlvaxEntero"),range_C=c(-30,-8),range_N=c(7,15))



#graphe_traj
traj_graph(sources_list=liste_sources,sources_approx=liste_approx,consu_list=liste_Conso_bloodcell,iso_num=1
           ,lambda_list=liste_lambda,conc_list=liste_Conc,list_TEF=liste_TDF,iter=iter,X=50,time_range=c(0,730),range_y=c(-30,-8),title="C"
           ,y_label=expression(paste(" δ"^{13},"C (‰)")))
traj_graph(sources_list=liste_sources,sources_approx=liste_approx,consu_list=liste_Conso_bloodcell,iso_num=2
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
sources_contrib(mod_results = res_1,obs_time=Dates,title="Entero",source_num=4,"topleft")

#contrib over time
contrib_time(res_model=SMM_50,model_name="SMM",date=Dates,data_feat=data_features(liste_sources)
             ,X=50,sources_name=c("Zos","grass","UlvaxEntero"),where_legend="topright")
contrib_time(res_model=SMM_delta_50,model_name="SMM_delta",date=Dates,data_feat=data_features(liste_sources)
             ,X=50,sources_name=c("Zos","grass","UlvaxEntero"),where_legend="topright")
contrib_time(res_model=DMM_50,model_name="DMM",date=Dates,data_feat=data_features(liste_sources)
             ,X=50,sources_name=c("Zos","grass","UlvaxEntero"),where_legend="topright")
#---------------------travail---------------------
# SIMM_3iso<-function(deltaN,deltaC,sA_C,sB_C,sC_C,sA_N,sB_N,sC_N,TDF_AN,TDF_BN,TDF_CN,TDF_AC,TDF_BC,TDF_CC){
#   pA<-(((sB_C+TDF_BC)-(sC_C+TDF_CC))*(deltaN-(sC_N+TDF_CN))-(deltaC-(sC_C+TDF_CC))*((sB_N+TDF_BN)-(sC_N+TDF_CN)))/(((sB_C+TDF_BC)-(sC_C+TDF_CC))*((sA_N+TDF_AN)-(sC_N+TDF_CN))+((sC_N+TDF_CN)-(sB_N+TDF_BN))*((sA_C+TDF_AC)-(sC_C+TDF_CC)))
#   pB<-(deltaC-(sC_C+TDF_CC)-pA*((sA_C+TDF_AC)-(sC_C+TDF_CC)))/((sB_C+TDF_BC)-(sC_C+TDF_CC))
#   pC<-1-pA-pB
#   return(data.frame(pA=pA,pB=pB,pC=pC))
# }
# 
# 
# B=Choose_data(sources_list=liste_sources,consu_list=liste_Conso_plasma,conc_list=liste_Conc,iter=iter,obs_nb=4)
# 
# SMM_sol_exacte=SIMM_3iso(deltaN=B[[3]][[2]][[2]],deltaC=B[[3]][[1]][[2]],sA_C=B[[2]][[1]][[1]],sB_C=B[[2]][[1]][[2]],sC_C=B[[2]][[1]][[3]],sA_N=B[[2]][[2]][[1]],sB_N=B[[2]][[2]][[2]],sC_N=B[[2]][[2]][[3]],TDF_AN=3.54,TDF_BN=3.54,TDF_CN=3.54,TDF_AC=1.63,TDF_BC=1.63,TDF_CC=1.63)
# 
# A=Choose_data(sources_list=liste_sources,consu_list=liste_Conso_plasma,conc_list=liste_Conc,iter=iter,obs_nb=3)
# DMM_sol_exacte=DMM_func(details_consu=A[[3]],details_consu_t1=A[[5]],sources_approx=liste_approx
#                         ,list_TEF=liste_TDF,details_conc=A[[4]],Combi=SMM_sol_exacte,X=1,lambda_list=liste_lambda,data_feat=data_features(liste_sources))
# DMM_tout=DMM_func(details_consu=A[[3]],details_consu_t1=A[[5]],sources_approx=liste_approx
#                         ,list_TEF=liste_TDF,details_conc=A[[4]],Combi=A[[1]],X=2,lambda_list=liste_lambda,data_feat=data_features(liste_sources))
#  Mel_C=data.frame(row.names=NULL)
#  Mel_N=data.frame(row.names=NULL)
# for (i in 1:nrow(A[[1]])) {
#   Mel_C=rbind(Mel_C,Mixing(t=60,sources_approx=liste_approx[[1]],p=A[[1]][i,],TEF=liste_TDF[[1]],Conc=A[[4]][[1]]))
#   Mel_N=rbind(Mel_N,Mixing(t=60,sources_approx=liste_approx[[2]],p=A[[1]][i,],TEF=liste_TDF[[2]],Conc=A[[4]][[2]]))}
#  Melange=data.frame(Mel_C,Mel_N)
# Mel_combi=cbind(Melange,A[[1]])
# 
# diff_C=diff_obs_pred(Mix=Mel_C,A[[3]][[1]][[2]])
# diff_N=diff_obs_pred(Mix=Mel_N,A[[3]][[2]][[2]])
# difference=cbind(diff_C,diff_N)
# 
# 
# # Mel_C_sol_ex=data.frame(Mixing(t=A[[1]][[1]][[1]],sources_approx=liste_approx[[1]],p=SMM_sol_exacte,TEF=liste_TDF[[1]],Conc=A[[4]][[1]]))
# # Mel_N_sol_ex=data.frame(Mixing(t=A[[1]][[1]][[1]],sources_approx=liste_approx[[2]],p=SMM_sol_exacte,TEF=liste_TDF[[2]],Conc=A[[4]][[2]]))
# # 
# # diff_C_sol_ex=diff_obs_pred(Mix=Mel_C_sol_ex,A[[3]][[1]][[2]])
# # diff_N_sol_ex=diff_obs_pred(Mix=Mel_N_sol_ex,A[[3]][[2]][[2]])
# # difference=cbind(diff_C_sol_ex,diff_N_sol_ex)
# 
# dist=distance(difference)
# 
# dist_combi=cbind(dist,A[[1]])
# dist_ordered=dist_combi[order(dist_combi[[1]],decreasing = TRUE),]
# sources_approx=liste_approx
# num_iso=1
# #dist_ordered[nrow(dist_ordered),-1]
# traj=RUN_TIMdyn2(t=seq(123,396,1)
#                  ,state0=c(X=A[[3]][[num_iso]][[2]]),
#                  par=c(lambda = approxfun(liste_lambda[[num_iso]][[1]], liste_lambda[[num_iso]][[2]],method = 'const',rule=2), 
#                        Xinf=approxfun(c(A[[3]][[num_iso]][[1]],A[[5]][[num_iso]][[1]]),
#                                       Mixing(t=c(A[[3]][[num_iso]][[1]],A[[5]][[num_iso]][[1]]),
#                                              sources_approx=sources_approx[[num_iso]],p=SMM_sol_exacte,TEF=liste_TDF[[num_iso]]
#                                              ,Conc=A[[4]][[num_iso]]),rule=2)))[1:2]
# 
# 
# plot(NULL,xlim=c(0,730),ylim=c(5,17),xlab="Time (d)",ylab="",main="title",cex.lab=1.3,las=1,cex.axis=1.1)
# mtext("delta 15 N",side=2,line=4,padj=1.5,cex=1.3,font=2)
# points(traj[[1]],traj[[2]])
# points(396,10.47,col="red")
# 
# didi=sqrt(((A[[3]][[1]][[2]]+12.83938)/A[[3]][[1]][[2]]+(A[[3]][[2]][[2]]-8.726203)/A[[3]][[2]][[2]])^2)
# 
#------------BIAS ESTIMATION---------------
bias_SMM=data.frame(row.names=NULL)
for (i in 1:nrow(DMM_1)) {bias_SMM=rbind(bias_SMM,abs(DMM_1[i,3]+DMM_1[i,4]
                                       -SMM_1[i+1,3]-SMM_1[i+1,4]))}
colnames(bias_SMM)<-c("bias")
lambda_T <- c(3.285,1.107,9.750,1.071,1.107,1.107,2.107,6.535)
bias_SMM=cbind(Time=DMM_1[[1]],bias_SMM,lambda_T)

# plot(bias_SMM$lambda_T,bias_SMM$bias)

#-----euclidian distance-----
Sources_Carbon=data.frame(Date=Dates,
                          Zos=rep(-11.17,9),
                          Grass=rep(-30.88,9),
                          Ulvaxentero=rep(-12.88,9))
Sources_Nitrogen=data.frame(Date=Dates,
                            Zos=rep(6.49,9),
                            Grass=rep(4.43,9),
                            Ulvaxentero=rep(10.5,9))

equ_val=data.frame(row.names=NULL)
for (i in 1:nrow(DMM_1)){equ_val=rbind(equ_val,
                                       cbind(
                                         equ_C=(Sources_Carbon[1,2]*DMM_1[i,3]
                                            +Sources_Carbon[1,3]*DMM_1[i,4]+
                                                Sources_Carbon[1,4]*DMM_1[i,5]),
                                       equ_N=(Sources_Nitrogen[1,2]*DMM_1[i,3]
                                          +Sources_Nitrogen[1,3]*DMM_1[i,4]+
                                            Sources_Nitrogen[1,4]*DMM_1[i,5])))}

val_equ_vs_conso=cbind(equ_val,Conso_bloodcell_Carbon[-1,2],Conso_bloodcell_Nitrogen[-1,2])
dist_euc=sqrt((val_equ_vs_conso$equ_C-val_equ_vs_conso$`Conso_bloodcell_Carbon[-1, 2]`)^2+(val_equ_vs_conso$equ_N-val_equ_vs_conso$`Conso_bloodcell_Nitrogen[-1, 2]`)^2)
bias_SMM=cbind(bias_SMM,dist_euc)

#setting colors +ordering col
deg_col=viridis(n=20, alpha = 1, begin = 0, end = 1, direction = -1, option = "D")
legend_image <- as.raster(matrix(rev(deg_col), ncol=1))
ord_col=order(bias_SMM[,4])
col_bias=NULL
for (i in 1:nrow(bias_SMM)){col_bias=c(col_bias,round(bias_SMM[i,4]*11/20))}
bias_SMM=cbind(bias_SMM,col_bias)
###ordering data
ord=order(bias_SMM$bias)
ord_bias_SMM=bias_SMM[ord,]
###graph
windows()
par(mfrow=c(1,2))
#legend
plot(c(0,2),c(0,1),type = 'n', axes = F,xlab = '', ylab = '', main = 'Euclidian distance')
# Tracer les graduations de la légende, ici le paramétrage en est strictement manuel 4 paliers entre 0 et 1
text(x=1.5, y = seq(0,1,l=6), labels = seq(1,11,l=6))
# Ajouter l'image de gradient
rasterImage(legend_image, 0, 0, 1,1)

plot(NA,NA,xlim=c(0,12),ylim=c(0,1),main="",xlab=c("λT"),ylab=c("Bias")
     , las=1, cex.lab=1.4,font.lab=2, cex.axis=1.4)
points(ord_bias_SMM$lambda_T,ord_bias_SMM$bias,type="p",col=deg_col[ord_bias_SMM$col_bias],pch=16,cex=1.4)
text(ord_bias_SMM$lambda_T+1.2,ord_bias_SMM$bias,labels = ord_bias_SMM$Time,font=2,cex=1)

##Big graph 
# windows()
# par(mfrow=c(2,2))
#SMM
ridgeline_chart(res_model=SMM_50,model_name="SMM",sources_names=c("Zos","grass","UlvaxEntero"),X=50)
#DMM
ridgeline_chart(res_model=DMM_50,model_name="DMM",sources_names=c("Zos","grass","UlvaxEntero"),X=50)
#bias
plot(NA,NA,xlim=c(0,10),ylim=c(0,1),main="",xlab=c("λT"),ylab=c("Bias")
     , las=1, cex.lab=1.4,font.lab=2, cex.axis=1.4)
points(ord_bias_SMM$lambda_T,ord_bias_SMM$bias,type="p",col=ordered_colors,pch=16,cex=1.4)
text(ord_bias_SMM$lambda_T+0.05,ord_bias_SMM$bias+0.01,labels = ord_bias_SMM$Time,font=2,cex=1)
#legend
plot(c(0,2),c(0,1),type = 'n', axes = F,xlab = '', ylab = '', main = 'Euclidian distance')
# Tracer les graduations de la légende, ici le paramétrage en est strictement manuel 4 paliers entre 0 et 1
text(x=1.5, y = seq(0,1,l=4), labels = seq(1,4.6,l=4))
# Ajouter l'image de gradient
rasterImage(legend_image, 0, 0, 1,1)
# install.packages("ggpubr")
library(ggpubr)
figure <-ggarrange(g_1,g_2,labels=c("a","b"),ncol=2,nrow=1)

