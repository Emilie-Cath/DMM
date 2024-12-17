rm(list=ls())# clear the current environment
graphics.off()# clear the current plots
setwd("C:/Users/Etudiant/Documents/codes_2")
source("fonctions_modeles_de_melange.R")
#install.packages("deSolve")
#install.packages("ggplot2")
library("deSolve")
library("ggplot2")
library("ggridges")

#---------------DATA---------------------------
iter=0.05
lambda_plasma=0.33
lambda_bloodcell=0.03
Dates=c(0,92,123,396,426,457,488,547,730)
Sources_Carbon=data.frame(Date=Dates,
                          Zos=rep(-11.17,9),
                          Grass=rep(-30.88,9),
                          Ulva=rep(-11.17,9),
                          Entero=rep(-14.06,9))
Sources_Nitrogen=data.frame(Date=Dates,
                            Zos=rep(6.49,9),
                            Grass=rep(4.43,9),
                            Ulva=rep(11.19,9),
                            Entero=rep(9.82,9))
liste_sources=list(Sources_Carbon,Sources_Nitrogen)

Conc_Carbon=data.frame(Zos=rep(0.35,9),Grass=rep(0.4,9),Ulva=rep(0.2,9),Entero=rep(0.18,9))
Conc_Nitrogen=data.frame(Zos=rep(0.03,9),Grass=rep(0.035,9),Ulva=rep(0.02,9),Entero=rep(0.01,9))
liste_Conc=list(Conc_Carbon,Conc_Nitrogen)

#Conc_Carbon=data.frame(Zos=rep(1,9),Grass=rep(1,9),Ulva=rep(1,9),Entero=rep(1,9))
#Conc_Nitrogen=data.frame(Zos=rep(0.01,9),Grass=rep(0.35,9),Ulva=rep(0.02,9),Entero=rep(0.01,9))
#liste_Conc=list(Conc_Carbon,Conc_Nitrogen)

TDF_Carbon=data.frame(Zos=1.63,Grass=1.63,Ulva=1.63,Entero=1.63)
TDF_Nitrogen=data.frame(Zos=3.54,Grass=3.54,Ulva=3.54,Entero=3.54)
liste_TDF=list(TDF_Carbon,TDF_Nitrogen)

#PLASMA
# Conso_plasma_Carbon=data.frame(Date=Dates,
#                                Iso_C=c(-12.09,-15.79,-21,-12.42,-14.84,-20.23,-25.7,-27.32,-10.8),
#                                lambda=rep(lambda_plasma,9))
# Conso_plasma_Nitrogen=data.frame(Date=Dates,
#                                  Iso_C=c(10.32,11.77,10.6,11,12.12,10.82,8.97,8.3,9.41),
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
            approxfun(Sources_Carbon$Date,Sources_Carbon$Ulva,method='linear',rule=2),
            approxfun(Sources_Carbon$Date,Sources_Carbon$Entero,method='linear',rule=2))
Sources_N=c(approxfun(Sources_Nitrogen$Date,Sources_Nitrogen$Zos,method='linear',rule=2),
            approxfun(Sources_Nitrogen$Date,Sources_Nitrogen$Grass,method='linear',rule=2),
            approxfun(Sources_Nitrogen$Date,Sources_Nitrogen$Ulva,method='linear',rule=2),
            approxfun(Sources_Nitrogen$Date,Sources_Nitrogen$Entero,method='linear',rule=2))
liste_approx=list(Sources_C,Sources_N)
#plasma
#liste_lambda=list(data.frame(Conso_plasma_Carbon$Date,Conso_plasma_Carbon$lambda),data.frame(Conso_plasma_Nitrogen$Date,Conso_plasma_Nitrogen$lambda))
#blood cells
liste_lambda=list(data.frame(Conso_bloodcell_Carbon$Date,Conso_bloodcell_Carbon$lambda),data.frame(Conso_bloodcell_Nitrogen$Date,Conso_bloodcell_Nitrogen$lambda))
#------------fait tourner---------------
#fait tourner
res=all_results(sources_list=liste_sources,data_feat=data_features(liste_sources),consu_list=liste_Conso_bloodcell
                ,conc_list=liste_Conc,iter=iter,list_TEF=liste_TDF,X=50,sources_approx=liste_approx,lambda_list=liste_lambda)

res_1=all_results(sources_list=liste_sources,data_feat=data_features(liste_sources),consu_list=liste_Conso_bloodcell
                  ,conc_list=liste_Conc,iter=iter,list_TEF=liste_TDF,X=1,sources_approx=liste_approx,lambda_list=liste_lambda)

# write.csv2(res[[1]],file="SMM_50_ge4.csv")
# write.csv2(res[[2]],file="SMM_delta_50_ge4.csv")
#write.csv2(res[[3]],file="DMM_50_ge4.csv")
SMM_50=read.csv2("SMM_50_ge4.csv")
SMM_50=SMM_50[-1]
SMM_delta_50=read.csv2("SMM_delta_50_ge4.csv")
SMM_delta_50=SMM_delta_50[-1]
DMM_50=read.csv2("DMM_50_ge4.csv")
DMM_50=DMM_50[-1]
res_50=c(list(SMM_50),list(SMM_delta_50),list(DMM_50))

# write.csv2(res_1[[1]],file="SMM_1_ge4.csv")
# write.csv2(res_1[[2]],file="SMM_delta_1_ge4.csv")
# write.csv2(res_1[[3]],file="DMM_1_ge4.csv")
SMM_1=read.csv2("SMM_1_ge4.csv")
SMM_1=SMM_1[-1]
SMM_delta_1=read.csv2("SMM_delta_1_ge4.csv")
SMM_delta_1=SMM_delta_1[-1]
DMM_1=read.csv2("DMM_1_ge4.csv")
DMM_1=DMM_1[-1]
res_1=c(list(SMM_1),list(SMM_delta_1),list(DMM_1))

#------------GRAPH---------------

#biplot visuals
sign_biplot(consu_list=liste_Conso_bloodcell,sources_list=liste_sources,list_TEF=liste_TDF
            ,title="",range_C=c(-30,-8),range_N=c(7,15),source_names=c("Zos","grass","Ulva","Entero"))
poly_sources_evol(consu_list=liste_Conso_bloodcell,sources_list=liste_sources,list_TEF=liste_TDF,sources_name=c("Zos","grass","Ulva","Entero"),range_C=c(-30,-8),range_N=c(7,15))



#graphe_traj
traj_graph(sources_list=liste_sources,sources_approx=liste_approx,consu_list=liste_Conso_bloodcell,iso_num=1
           ,lambda_list=liste_lambda,conc_list=liste_Conc,list_TEF=liste_TDF,iter=iter,X=50,time_range=c(0,730),range_y=c(-30,-8),title="C"
           ,y_label=expression(paste(" δ"^{13},"C (‰)")))
traj_graph(sources_list=liste_sources,sources_approx=liste_approx,consu_list=liste_Conso_bloodcell,iso_num=2
           ,lambda_list=liste_lambda,conc_list=liste_Conc,list_TEF=liste_TDF,iter=iter,X=50,time_range=c(0,730),range_y=c(7,15),title="N"
           ,y_label=expression(paste(" δ"^{15},"N (‰)")))

#ridgeline
ridgeline_chart(res_model=SMM_50,model_name="SMM",sources_names=c("Zos","grass","Ulva","Entero"),X=50)
ridgeline_chart(res_model=SMM_delta_50,model_name="SMM_delta",sources_names=c("Zos","grass","Ulva","Entero"),X=50)
ridgeline_chart(res_model=DMM_50,model_name="DMM",sources_names=c("Zos","grass","Ulva","Entero"),X=50)

#distance
dist_graph(obs_time=Dates,res_model=res_50,range_y=c(0,0.25))

#source contrib
sources_contrib(mod_results = res_1,obs_time=Dates,title="Zoostera",source_num=1,"topleft")
sources_contrib(mod_results = res_1,obs_time=Dates,title="Grass",source_num=2,"topleft")
sources_contrib(mod_results = res_1,obs_time=Dates,title="Ulva",source_num=3,"topleft")
sources_contrib(mod_results = res_1,obs_time=Dates,title="Entero",source_num=4,"topleft")


#aggregation

DMM_50_agg=data.frame(row.names = NULL)
for (i in 1:nrow(DMM_50)) {
  DMM_50_agg=rbind(DMM_50_agg,DMM_50[i,5]+DMM_50[i,6])
}
DMM_50_agg=cbind(DMM_50[-c(5,6)],DMM_50_agg)
colnames(DMM_50_agg)[5] <-c("Var3")
ridgeline_chart(res_model=DMM_50_agg,model_name="DMM",sources_names=c("Zos","grass","Ulva xEntero"),X=50)


SMM_50_agg=data.frame(row.names = NULL)
for (i in 1:nrow(SMM_50)) {
  SMM_50_agg=rbind(SMM_50_agg,SMM_50[i,5]+SMM_50[i,6])
}
SMM_50_agg=cbind(SMM_50[-c(5,6)],SMM_50_agg)
colnames(SMM_50_agg)[5] <-c("Var3")
ridgeline_chart(res_model=SMM_50_agg,model_name="SMM",sources_names=c("Zos","grass","Ulva xEntero"),X=50)


SMM_delta_50_agg=data.frame(row.names = NULL)
for (i in 1:nrow(SMM_delta_50)) {
  SMM_delta_50_agg=rbind(SMM_delta_50_agg,SMM_delta_50[i,6]+SMM_delta_50[i,7])
}
SMM_delta_50_agg=cbind(SMM_delta_50[-c(6,7)],SMM_delta_50_agg)
colnames(SMM_delta_50_agg)[6] <-c("Var3")
ridgeline_chart(res_model=SMM_delta_50_agg,model_name="SMM_delta",sources_names=c("Zos","grass","Ulva xEntero"),X=50)

#--------------------------
traj_graph_2 <- function(sources_list,sources_approx,consu_list,iso_num,lambda_list,conc_list,list_TEF,iter,X,time_range,range_y,title,y_label) {
  windows()
  plot(NULL,xlim=time_range,ylim=range_y,xlab="Time (d)",ylab="",main=title,cex.lab=1.3,las=1,cex.axis=1.1)
  mtext(y_label,side=2,line=4,padj=1.5,cex=1.3,font=2)
  data_feat=data_features(sources_list)
  for (j in 1:(data_feat[3]-1)){ #pour toutes les observations -1 
    A=Choose_data(sources_list=sources_list,consu_list=consu_list,conc_list=conc_list,iter=iter,obs_nb=j)
    
    for (k in 1:nrow(A[[1]])) {mix=RUN_TIMdyn2(t=seq(A[[3]][[iso_num]][[1]],A[[5]][[iso_num]][[1]],1)
                                            ,state0=c(X=A[[3]][[iso_num]][[2]]),
                                            par=c(lambda = approxfun(lambda_list[[iso_num]][[1]], lambda_list[[iso_num]][[2]],method = 'const',rule=2), 
                                                  Xinf=approxfun(c(A[[3]][[iso_num]][[1]],A[[5]][[iso_num]][[1]]),
                                                                 Mixing(t=c(A[[3]][[iso_num]][[1]],A[[5]][[iso_num]][[1]]),
                                                                        sources_approx=sources_approx[[iso_num]],p=A[[1]][k,],TEF=list_TEF[[iso_num]]
                                                                        ,Conc=A[[4]][[iso_num]]),rule=2)))[1:2]
    points(mix[[1]],mix[[2]],type="l",col="grey")}
  }
  points(consu_list[[iso_num]][[1]],consu_list[[iso_num]][[2]],type="p",
         col="red",pch=16) #on ajoute les valeurs observées
  
}

traj_graph_2(sources_list=liste_sources,sources_approx=liste_approx,consu_list=liste_Conso_bloodcell,iso_num=2
           ,lambda_list=liste_lambda,conc_list=liste_Conc,list_TEF=liste_TDF,iter=iter,X=50,time_range=c(0,730),range_y=c(7,15),title="N"
           ,y_label=expression(paste(" δ"^{15},"N (‰)")))

A=Choose_data(sources_list=liste_sources,consu_list=liste_Conso_bloodcell,conc_list=liste_Conc,iter=iter,obs_nb=3)


Mel_C=data.frame(row.names=NULL)
Mel_N=data.frame(row.names=NULL)
for (i in 1:nrow(A[[1]])) {
  Mel_C=rbind(Mel_C,Mixing(t=60,sources_approx=liste_approx[[1]],p=A[[1]][i,],TEF=liste_TDF[[1]],Conc=A[[4]][[1]]))
  Mel_N=rbind(Mel_N,Mixing(t=60,sources_approx=liste_approx[[2]],p=A[[1]][i,],TEF=liste_TDF[[2]],Conc=A[[4]][[2]]))}
Melange=data.frame(Mel_C,Mel_N)
Mel_combi=cbind(Melange,A[[1]])

diff_C=diff_obs_pred(Mix=Mel_C,-21.26)
diff_N=diff_obs_pred(Mix=Mel_N,9.18)
difference=cbind(diff_C,diff_N)

dist=distance(difference)

dist_combi=cbind(dist,A[[1]])
dist_ordered=dist_combi[order(dist_combi[[1]],decreasing = TRUE),]
sources_approx=liste_approx
num_iso=2
traj=RUN_TIMdyn2(t=seq(123,396,1)
                 ,state0=c(X=A[[3]][[num_iso]][[2]]),
                 par=c(lambda = approxfun(liste_lambda[[num_iso]][[1]], liste_lambda[[num_iso]][[2]],method = 'const',rule=2), 
                       Xinf=approxfun(c(A[[3]][[num_iso]][[1]],A[[5]][[num_iso]][[1]]),
                                      Mixing(t=c(A[[3]][[num_iso]][[1]],A[[5]][[num_iso]][[1]]),
                                             sources_approx=sources_approx[[num_iso]],p=dist_ordered[nrow(dist_ordered),-1],TEF=liste_TDF[[num_iso]]
                                             ,Conc=A[[4]][[num_iso]]),rule=2)))[1:2]


plot(NULL,xlim=c(0,730),ylim=c(5,17),xlab="Time (d)",ylab="",main="title",cex.lab=1.3,las=1,cex.axis=1.1)
mtext("delta 15 N",side=2,line=4,padj=1.5,cex=1.3,font=2)
points(traj[[1]],traj[[2]])
points(396,10.47,col="red")


