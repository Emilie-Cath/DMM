####VERSION GENERALISEE AVEC MES CODES###

rm(list=ls())# clear the current environment
graphics.off()# clear the current plots
setwd("C:/Users/Etudiant/Documents/thèse/1ere année/1er papier/finals/1_st paper")
source("C:/Users/Etudiant/Documents/thèse/1ere année/1er papier/finals/1_st paper/R/Mixing_models_functions.R")
library(deSolve)
library(philentropy)
library("viridisLite")



#-------------ISOTOPIC SPACE EFFECT--------------------

colors=c( "#CC79A7"  ,"#E69F00","#56B4E9","#009E73"  ,"#D55E00"  )
#CHANGE HERE IF THE SOURCE NUMBER OR ISOTOPE NUMBER CHANGE
TEF=list(data.frame("s1"=c(0),"s2"=c(0),"s3"=c(0)),data.frame("s1"=c(0),"s2"=c(0),"s3"=c(0))) #the TEF=0 for all sources and isotope
conc=list(data.frame("s1"=c(1),"s2"=c(1),"s3"=c(1)),data.frame("s1"=c(1),"s2"=c(1),"s3"=c(1))) #the conc=1 for all sources and isotope
iter=0.05

#ISOTOPIC SPACE
nb_iso_spaces=5
s3=data.frame(I1=rep(5,nb_iso_spaces),I2=rep(10,nb_iso_spaces))
s1=data.frame(I1=c(10,6,7,12,6),I2=c(5,9,6,2,9))
s2=data.frame(I1=c(0,5,2,-5,4),I2=c(0,8,2,-5,5))

plot(NULL,xlim=c(-10,15),ylim=c(-10,10),xlab="",ylab="",las=1)
mtext(expression(paste(" δ"^{15},"N (‰)")),side=2,line=4,padj=1.5,at=1,cex=1.2,font=2)
mtext(expression(paste(" δ"^{13},"C (‰)")),side=1,line=1,padj=1.5,at=5,cex=1.2,font=2)
for (i in 1:nrow(s3)){
  points(c(s1[i,1],s2[i,1],s3[i,1],s1[i,1]),c(s1[i,2],s2[i,2],s3[i,2],s1[i,2]),type="l",col=colors[i],lwd = 2)
}

#ESTIMATING EUCLIDIAN DISTANCE
#coord1=c(coordi1,coordi2)










euclidian_distance <-function(coord1,coord2){
  return(sqrt((coord1[1]-coord2[1])^2+(coord1[2]-coord2[2])^2))
}




#DIET
#the consumer eats 100% of source 3
combi=data.frame(s1=c(0),s2=c(0),s3=c(1))


#RESULTS
resum_data_SMM=data.frame(B=c(),l_T=c(),dist_eucli=c(),iso_space=c())

p=seq(2,300,length=4) #period
l=seq(0.002,0.02,length=5) #lambda

for (s in 1:length(l)){ #for every lambda
  # s=1
  lambda=l[s]
  print(c("l",s))
  
  for (f in 1:length(p)){ #for every period
    # f=1
    print(c("p",f))
    period=round(p[f])
    t <- seq(0,period,1) # output times for the ODE (d)
    liste_lambda=list(data.frame(Date=t,l=rep(lambda,length(t))),data.frame(Date=t,l=rep(lambda,length(t))))
    
    ##THE LOOP STARTS HERE
    for(m in 1:nrow(s3)){ #for each consumer/ isotopic space
      # m=1
      #SOURCES VALUES CONSTANT DEPENDING ON ISOTOPIC SPACE
      source1=s1[m,]
      source2=s2[m,]
      source3=s3[m,]
      sources_list=list(data.frame(Date=t,s1=rep(source1[[1]],period+1),s2=rep(source2[[1]],period+1),s3=rep(source3[[1]],period+1)), #pattern a
                            data.frame(Date=t,s1=rep(source1[[1]],period+1),s2=rep(source2[[2]],period+1),s3=rep(source3[[2]],period+1)))
      
      
      conso_ini=s2[m,]
      print(m)
 
        
        #CREATING THE SOURCE APPROXFUN
        sources_approxfun=list(c(approxfun(sources_list[[1]]$Date,sources_list[[1]]$s1,method='linear',rule=2)
                                 ,approxfun(sources_list[[1]]$Date,sources_list[[1]]$s2,method='linear',rule=2)
                                 ,approxfun(sources_list[[1]]$Date,sources_list[[1]]$s3,method='linear',rule=2)),
                               c(approxfun(sources_list[[2]]$Date,sources_list[[2]]$s1,method='linear',rule=2)
                                 ,approxfun(sources_list[[2]]$Date,sources_list[[2]]$s2,method='linear',rule=2)
                                 ,approxfun(sources_list[[2]]$Date,sources_list[[2]]$s3,method='linear',rule=2)))
        
        #EXPLAINS HOW MANY ISO,SOURCES AND OBS
        data_feat=data_features(sources_list )
        res=NULL
        ind=NULL
        new_conso_value=NULL
        #GENERATING DX AND DXinf VALUES
        for (k in 1:data_feat[1]) #FOR EACH ISOTOPE
        {res=c(res,list(RUN_TIMdyn2(t,#time
                                    state0=c(X=conso_ini[[k]]),#initial conditions
                                    par=c(#two parameters lambda and Xinf i.e. forcing functions 
                                      lambda = approxfun(liste_lambda[[k]]$Date,liste_lambda[[k]]$l,method = 'const',rule=2), 
                                      Xinf=approxfun(t, Mixing(t=t,sources_approx=sources_approxfun[[k]],p=combi,TEF=TEF[[k]]
                                                               ,Conc=conc[[k]])))))) #COMPUTING THE DX AND DXinf

        
        
        new_conso_value=c(new_conso_value,list(data.frame(t=seq(0,period,1),DX=res[[k]][[2]])))
        }

        
        #USING DX VALUES AS DATA, estimating the SMM and SMM delta
        
        
        #GENERATING THE POSSIBLE SOLUTIONS
        sol<- expand.grid(rep(list(seq(0,1,iter)),data_feat[2])) #data_feat[2] is the number of food sources
        sol <- subset(sol,rowSums(sol)==1)

        
        #SMM and SMM_delta estimation at the last "t"
        val_conso=NULL
        val_0=NULL
        for (j in 1:data_feat[[1]]) {val_conso=c(val_conso,list(new_conso_value[[j]][period,]))
        val_0=c(val_0,list(new_conso_value[[j]][1,]))}
        
        SMM=SMM_func(t=period,sources_approx=sources_approxfun,list_TEF=TEF,details_conc=conc
                     ,combi=sol,data_feat=data_feat,details_consu=val_conso,X=nrow(sol))
 
        
        #CALCULATING EUCLIDIAN DISTANCE BETWEEN DIST INI AND EQ
        
   
        
        resum_data_SMM=rbind(resum_data_SMM,data.frame(B=(abs(combi[[1]]-SMM[[2]])+abs(combi[[2]]-SMM[[3]])+abs(combi[[3]]-SMM[[4]]))*50
                                                       ,lambda_T=c(lambda*period)
                                                       ,dist_eucli=c(euclidian_distance(c(source2[[1]],source2[[2]]),c(source3[[1]],source3[[2]])))
                                                       ,iso_space=c(m),SMM,Dc_obs=val_conso[[1]]$DX,Dn_obs=val_conso[[2]]$DX
                                                       ,s1_C=source1[[1]],s2_C=source2[[1]],s3_C=source3[[1]]
                                                       ,s1_N=source1[[2]],s2_N=source2[[2]],s3_N=source3[[2]]))

    }#loop on ini value
    
  }#loop on period
}#loop on lambda

#ADDING DIFFERENCE PRED-OBS IN THE DATA FRAME

diff_C_SMM=NULL
diff_N_SMM=NULL
Dc_pred_SMM=NULL
Dn_pred_SMM=NULL
Dn_pred_SMM_delta=NULL
for (i in 1: nrow(resum_data_SMM)) {
  #SMM
  dc_SMM=resum_data_SMM$s1_C[i]*resum_data_SMM$Var1[i]+resum_data_SMM$s2_C[i]*resum_data_SMM$Var2[i]+resum_data_SMM$s3_C[i]*resum_data_SMM$Var3[i]
  dn_SMM=resum_data_SMM$s1_N[i]*resum_data_SMM$Var1[i]+resum_data_SMM$s2_N[i]*resum_data_SMM$Var2[i]+resum_data_SMM$s3_N[i]*resum_data_SMM$Var3[i]
  Dc_pred_SMM=c(Dc_pred_SMM,dc_SMM)
  Dn_pred_SMM=c(Dn_pred_SMM,dn_SMM)
  diff_C_SMM=c(diff_C_SMM,abs(resum_data_SMM$Dc_obs[i]-dc_SMM))
  diff_N_SMM=c(diff_N_SMM,abs(resum_data_SMM$Dn_obs[i]-dn_SMM))
}

#adding it to the dataset
data_SMM_dX_pred=cbind(resum_data_SMM,Dc_pred_SMM,Dn_pred_SMM,diff_C_SMM,diff_N_SMM)

##KEEPING A TOLERANCE
data_tol=data_SMM_dX_pred[data_SMM_dX_pred$diff_C_SMM<0.2&data_SMM_dX_pred$diff_N_SMM<0.2,]


#CHECKING THE NUMBER OF SOLUTIONS
lambda_T=unique(data_SMM_dX_pred$lambda_T)

df_nb_sol=data.frame(row.names = FALSE)
for (i in 1:nb_iso_spaces){#loop on isospaces
  df_sources=subset(data_tol,data_tol$iso_space==i)
    for (k in 1:length(lambda_T)){#loop on lambdaT
      val=lambda_T[k]
      df_lambda=subset(df_sources,lambda_T==val)
      df_nb_sol=rbind(df_nb_sol,data.frame(isospace=i,lambda_T=val,nb_sol=nrow(df_lambda)))
    }#loop on lambdaT
}#loop on isospaces

#stats
data_stats=data.frame(row.names = FALSE)

for (j in 1:nb_iso_spaces){ #loop on isospace
  data_tol_3=subset(data_tol,iso_space==j)
  for (i in 1:length(lambda_T)){#loop on lambda T
    val=lambda_T[i]
    subs_1=subset(data_tol_3,lambda_T==val)
    info=as.numeric(quantile(subs_1$B,probs=c(0.25,0.5,0.75)))
    data_stats=rbind(data_stats,data.frame(iso_space=j,lambda_T=val,Q1=info[1],Med=info[2],Q3=info[3]))
  }#loop on lambda T
}#loop on isospace


##KEEPING THE 5% best solutions
# ### 5% BEST SOLUTION
# nb_sol=11
# best_sol_data=data.frame(row.names = NULL)
# data_stats=data.frame(row.names = NULL)
# for (k in 1:nb_iso_spaces ){ #loop sources
#   source_data=subset(data_SMM_dX_pred,iso_space==k)
#   for (i in 1:length(lambda_T)){val=lambda_T[i] #loop λT
#   subs_2=subset(source_data,lambda_T==val)
#   best_sol=tail(subs_2,nb_sol)
#   stat=as.numeric(quantile(best_sol$B,probs=c(0.25,0.5,0.75)))
#   best_sol_data=rbind(best_sol_data,best_sol)
#   data_stats=rbind(data_stats,data.frame(iso_space=k,lambda_T=val,Q1=stat[1],M=stat[2],Q3=stat[3]))
#   }#loop λT
# }#loop sources



size=rev(seq(1,2,length=nb_iso_spaces))
#GRAPH
plot(NULL,xlim=c(0,6),ylim=c(0,100),xlab="λT",ylab="β(%)",las=1)
for (i in 1:nb_iso_spaces){
  subs=subset(data_stats,iso_space==i)
  points(subs$lambda_T,subs$Q1,col=colors[i],type="p",pch=3,cex=size[i])
  points(subs$lambda_T,subs$Q3,col=colors[i],type="p",pch=3,cex=size[i])
}
for (i in 1:nb_iso_spaces){
  subs=subset(data_stats,iso_space==i)
  points(subs$lambda_T,subs$M,col=colors[i],type="p",pch=16,cex=size[i])
}

#----in silico experiment (2I,3S)----

#CHANGE HERE IF THE SOURCE NUMBER OR ISOTOPE NUMBER CHANGE
TEF=list(data.frame("s1"=c(0),"s2"=c(0),"s3"=c(0)),data.frame("s1"=c(0),"s2"=c(0),"s3"=c(0))) #the TEF=0 for all sources and isotope
conc=list(data.frame("s1"=c(1),"s2"=c(1),"s3"=c(1)),data.frame("s1"=c(1),"s2"=c(1),"s3"=c(1))) #the conc=1 for all sources and isotope
iter=0.05


#CONSUMER INITAL VALUE
nb_ini_val=20

resum_data_SMM=data.frame(B=c(),l_T=c(),sources=c(),val_ini=c(),DC=c(),DN=c())
resum_data_SMM_delta=data.frame(B=c(),l_T=c(),sources=c(),val_ini=c())

initial_values=read.csv2("C:/Users/Etudiant/Documents/thèse/1ere année/1er papier/finals/1_st paper/Data/val_ini_2.csv",row.names = NULL)


#DIET

#the consumer eats 100% of source 3
# combi=data.frame(s1=c(0),s2=c(0),s3=c(1))
# conso=initial_values[c(12,14,20),]
#theconsumer eats a mix and the solution value is (3.5,4)
combi=data.frame(s1=c(0.15),s2=c(0.45),s3=c(0.4))
conso=initial_values[c(13,7,1),]


# ###ISOTOPIC SPACE
# #3 sources
# plot(NA, NA, ylim = c(-5,10), xlim = c(-5,10),
#      xlab = "", ylab = "" ,las=1,main="a")
# mtext(expression(paste(" δ"^{15},"N (‰)")),side=2,line=4,padj=1.5,at=5,cex=1.5,font=2)
# mtext(expression(paste(" δ"^{13},"C (‰)")),side=1,line=1,padj=1.5,at=5,cex=1.5,font=2)
# points(conso$I1,conso$I2,pch=16)
# points(c(0,5,10,0),c(0,10,5,0),type="l",col="red")
# points(c(0,5,10,0),c(0,10,5,0),type="p",pch=16,col="red")
# text(x=c(0.5,5.5,9.5),y=c(0,10,5),labels = c("2","3","1"),cex=1.4,font=2,col="red")
# 
# 
# 


#SETTING λT values
p=seq(2,300,length=4) #T
l=seq(0.002,0.02,length=5) #lambda

#keeping estimated DX values 
DX=data.frame(lambda_T=c(),val_ini=c(),t_f=c(),DC=c(),DN=c())

for (s in 1:length(l)){ #LOOP ON LAMBDA 
  # s=1
  lambda=l[s]
  print(c("l",s))
  
  for (f in 1:length(p)){ #LOOP ON T 
    # f=1
    print(c("p",f))
    period=round(p[f])
    t <- seq(0,period,1) # output times for the ODE (d)
    liste_lambda=list(data.frame(Date=t,l=rep(lambda,length(t))),data.frame(Date=t,l=rep(lambda,length(t))))
    #SOURCES VALUES 
    sources_list_var=c(list(list(data.frame(Date=t,s1=rep(10,period+1),s2=rep(0,period+1),s3=rep(5,period+1)), #pattern a
                                 data.frame(Date=t,s1=rep(5,period+1),s2=rep(0,period+1),s3=rep(10,period+1))))
                       ,
                       
                       list(list(data.frame(Date=t,s1=rep(10,period+1),s2=seq(-5,0,length=period+1),s3=seq(15,3,length=period+1)), #pattern b
                                 data.frame(Date=t,s1=rep(5,period+1),s2=rep(0,period+1),s3=rep(10,period+1)))),
                       
                       list(list(data.frame(Date=t,s1=rep(10,period+1),s2=seq(0,-5,length=period+1),s3=seq(0,-20,length=period+1)), #pattern c
                                 data.frame(Date=t,s1=rep(5,period+1),s2=rep(0,period+1),s3=rep(10,period+1)))),
                       
                       list(list(data.frame(Date=t,s1=rep(10,period+1),s2=seq(5,10,length=period+1),s3=seq(0,20,length=period+1)), #pattern d
                                 data.frame(Date=t,s1=rep(5,period+1),s2=rep(0,period+1),s3=rep(10,period+1)))),
                       
                       list(list(data.frame(Date=t,s1=rep(10,period+1),s2=seq(0,-5,length=period+1),s3=seq(15,0,length=period+1)), #pattern e
                                 data.frame(Date=t,s1=rep(5,period+1),s2=rep(0,period+1),s3=rep(10,period+1))))
    )
    
    
    
    
    
    for(m in 1:nrow(conso)){ #LOOP ON CONSUMER
    # m=20
    conso_ini=conso[m,]
    print(m)
    
    indicator_value=NULL
    bias_SMM=NULL
    bias_SMM_delta=NULL
    bias_DMM=NULL
    dist_KL_SMM=NULL
    dist_KL_SMM_delta=NULL
    SMM_arch=data.frame(row.names = NULL)
    SMM_delta_arch=data.frame(row.names = NULL)
    DMM_arch=data.frame(row.names = NULL)
    
    
    for (n in 1:length(sources_list_var)) { #LOOP ON SOURCE VARIATION
    # n=1
    
    sources_list=sources_list_var[[n]]
    
    #CREATING THE SOURCE APPROXFUN
    sources_approxfun=list(c(approxfun(sources_list[[1]]$Date,sources_list[[1]]$s1,method='linear',rule=2)
                             ,approxfun(sources_list[[1]]$Date,sources_list[[1]]$s2,method='linear',rule=2)
                             ,approxfun(sources_list[[1]]$Date,sources_list[[1]]$s3,method='linear',rule=2)),
                           c(approxfun(sources_list[[2]]$Date,sources_list[[2]]$s1,method='linear',rule=2)
                             ,approxfun(sources_list[[2]]$Date,sources_list[[2]]$s2,method='linear',rule=2)
                             ,approxfun(sources_list[[2]]$Date,sources_list[[2]]$s3,method='linear',rule=2)))
    
    #EXPLAINS HOW MANY ISO,SOURCES AND OBS
    data_feat=data_features(sources_list )
    res=NULL
    ind=NULL
    new_conso_value=NULL
    #GENERATING DX AND DXinf VALUES
    for (k in 1:data_feat[1]) #FOR EACH ISOTOPE
    {res=c(res,list(RUN_TIMdyn2(t,#time
                                state0=c(X=conso_ini[[k]]),#initial conditions
                                par=c(#two parameters lambda and Xinf i.e. forcing functions 
                                  lambda = approxfun(liste_lambda[[k]]$Date,liste_lambda[[k]]$l,method = 'const',rule=2), 
                                  Xinf=approxfun(t, Mixing(t=t,sources_approx=sources_approxfun[[k]],p=combi,TEF=TEF[[k]]
                                                           ,Conc=conc[[k]])))))) #COMPUTING THE DX AND DXinf
    
    
    
    new_conso_value=c(new_conso_value,list(data.frame(t=seq(0,period,1),DX=res[[k]][[2]])))
    }
    
    #USING DX VALUES AS DATA, estimating the SMM and SMM delta
    
    
    #GENERATING THE POSSIBLE SOLUTIONS
    sol<- expand.grid(rep(list(seq(0,1,iter)),data_feat[2])) #data_feat[2] is the number of food sources
    sol <- subset(sol,rowSums(sol)==1)
    
    ###SMM and SMM_delta estimation at the last "t"
    val_conso=NULL
    val_0=NULL
    for (j in 1:data_feat[[1]]) {val_conso=c(val_conso,list(new_conso_value[[j]][period,]))
    val_0=c(val_0,list(new_conso_value[[j]][1,]))}
    
    #KEEPING DX VALUES AT THE LAST "t"
    DX=rbind(DX,data.frame(lambda_T=lambda*period,val_ini=m,source=n,t_f=val_conso[[1]]$t,DC=val_conso[[1]]$DX,
                           DN=val_conso[[2]]$DX,s1_C=sources_approxfun[[1]][[1]](period)
                           ,s2_C=sources_approxfun[[1]][[2]](period),s3_C=sources_approxfun[[1]][[3]](period)
                           ,s1_N=sources_approxfun[[2]][[1]](period),s2_N=sources_approxfun[[2]][[2]](period)
                           ,s3_N=sources_approxfun[[2]][[3]](period)))
    
    SMM=SMM_func(t=period,sources_approx=sources_approxfun,list_TEF=TEF,details_conc=conc
                 ,combi=sol,data_feat=data_feat,details_consu=val_conso,X=nrow(sol))
    
    SMM_delta=SMM_delta_func(t=period,sources_approx=sources_approxfun,list_TEF=TEF,details_conc=conc,combi=sol
                             ,data_feat=data_feat,details_consu=val_conso,X=nrow(sol),lambda_list=liste_lambda)


    
    #BIAS TT SYSTEME
    resum_data_SMM=rbind(resum_data_SMM,data.frame(B=(abs(combi[[1]]-SMM[[2]])+abs(combi[[2]]-SMM[[3]])+abs(combi[[3]]-SMM[[4]]))*50,lambda_T=c(lambda*period)
                                                   ,sources=c(n),val_ini=c(m),SMM,DC=val_conso[[1]]$DX,DN=val_conso[[2]]$DX
                                                   ,s1_C=sources_approxfun[[1]][[1]](period),s2_C=sources_approxfun[[1]][[2]](period),s3_C=sources_approxfun[[1]][[3]](period)
                                                   ,s1_N=sources_approxfun[[2]][[1]](period),s2_N=sources_approxfun[[2]][[2]](period)
                                                   ,s3_N=sources_approxfun[[2]][[3]](period)))
    resum_data_SMM_delta=rbind(resum_data_SMM_delta,data.frame(B=(abs(combi[[1]]-SMM_delta[[2]])+abs(combi[[2]]-SMM_delta[[3]])+abs(combi[[3]]-SMM_delta[[4]]))*50,lambda_T=c(lambda*period)
                                                   ,sources=c(n),val_ini=c(m),SMM_delta,DC=val_conso[[1]]$DX,DN=val_conso[[2]]$DX
                                                   ,s1_C=sources_approxfun[[1]][[1]](period),s2_C=sources_approxfun[[1]][[2]](period),s3_C=sources_approxfun[[1]][[3]](period)
                                                   ,s1_N=sources_approxfun[[2]][[1]](period),s2_N=sources_approxfun[[2]][[2]](period)
                                                   ,s3_N=sources_approxfun[[2]][[3]](period)))

    }#loop on sources
    
    }#loop on ini value
    
  }#loop on period
}#loop on lambda




#ADDING DIFFERENCE PRED-OBS IN THE DATA FRAME

diff_C_SMM=NULL
diff_N_SMM=NULL
Dc_pred_SMM=NULL
Dn_pred_SMM=NULL
diff_C_SMM_delta=NULL
diff_N_SMM_delta=NULL
Dc_pred_SMM_delta=NULL
Dn_pred_SMM_delta=NULL
for (i in 1: nrow(resum_data_SMM)) {
  #SMM
  dc_SMM=resum_data_SMM$s1_C[i]*resum_data_SMM$Var1[i]+resum_data_SMM$s2_C[i]*resum_data_SMM$Var2[i]+resum_data_SMM$s3_C[i]*resum_data_SMM$Var3[i]
  dn_SMM=resum_data_SMM$s1_N[i]*resum_data_SMM$Var1[i]+resum_data_SMM$s2_N[i]*resum_data_SMM$Var2[i]+resum_data_SMM$s3_N[i]*resum_data_SMM$Var3[i]
  Dc_pred_SMM=c(Dc_pred_SMM,dc_SMM)
  Dn_pred_SMM=c(Dn_pred_SMM,dn_SMM)
  diff_C_SMM=c(diff_C_SMM,abs(resum_data_SMM$DC[i]-dc_SMM))
  diff_N_SMM=c(diff_N_SMM,abs(resum_data_SMM$DN[i]-dn_SMM))
  #SMM_delta
  dc_SMM_delta=resum_data_SMM_delta$s1_C[i]*resum_data_SMM_delta$Var1[i]+resum_data_SMM_delta$s2_C[i]*resum_data_SMM_delta$Var2[i]+resum_data_SMM_delta$s3_C[i]*resum_data_SMM_delta$Var3[i]
  dn_SMM_delta=resum_data_SMM_delta$s1_N[i]*resum_data_SMM_delta$Var1[i]+resum_data_SMM_delta$s2_N[i]*resum_data_SMM_delta$Var2[i]+resum_data_SMM_delta$s3_N[i]*resum_data_SMM_delta$Var3[i]
  Dc_pred_SMM_delta=c(Dc_pred_SMM_delta,dc_SMM_delta)
  Dn_pred_SMM_delta=c(Dn_pred_SMM_delta,dn_SMM_delta)
  diff_C_SMM_delta=c(diff_C_SMM_delta,abs(resum_data_SMM_delta$DC[i]-dc_SMM_delta))
  diff_N_SMM_delta=c(diff_N_SMM_delta,abs(resum_data_SMM_delta$DN[i]-dn_SMM_delta))
}

#adding it to the dataset
  data_SMM_dX_pred=cbind(resum_data_SMM,Dc_pred_SMM,Dn_pred_SMM,diff_C_SMM,diff_N_SMM)
  data_SMM_delta_dX_pred=cbind(resum_data_SMM_delta,Dc_pred_SMM_delta,Dn_pred_SMM_delta
                               ,diff_C_SMM_delta,diff_N_SMM_delta)
#--------------------------------------------------#
#SWITCH EFFECT
# write.csv2(data_dX_pred,file="total_dataset_diet1.csv",row.names = F)
# write.csv2(DX,file="diet_1_SMM_Values.csv",row.names = F)
# data_dX_pred=read.csv2("total_dataset_diet1.csv")
# DX=read.csv2("diet_1_SMM_Values.csv")


#CHECKING THE NUMBER OF SOLUTIONS
lambda_T=unique(data_SMM_dX_pred$lambda_T)
# data_tol=data_dX_pred[data_dX_pred$sources==5&data_dX_pred$diff_C<1&data_dX_pred$diff_N<1,]
data_tol=data_SMM_dX_pred[data_SMM_dX_pred$diff_C<0.5&data_SMM_dX_pred$diff_N<0.5,]
nb_sources=5
df_nb_sol=data.frame(row.names = FALSE)
for (i in 1:nb_sources){#loop on sources
  df_sources=subset(data_tol,data_tol$sources==i)
  for (j in 1:nrow(conso)){#loop on val ini
    df_val_ini=subset(df_sources,df_sources$val_ini==j)
    for (k in 1:length(lambda_T)){#loop on lambdaT
      val=lambda_T[k]
      df_lambda=subset(df_val_ini,df_val_ini$lambda_T==val)
      df_nb_sol=rbind(df_nb_sol,data.frame(source=i,val_ini=j,lambda_T=val,nb_sol=nrow(df_lambda)))
    }#loop on lambdaT
  }#loop on val ini
}#loop on sources

df_0s=subset(df_nb_sol,nb_sol==0)



# #SWITCH EFFECT GRAPH
# data_stats=data.frame(row.names = FALSE)
# for(j in 1:nrow(conso)){
#   data_val_3=subset(data_dX_pred,val_ini==j)
# 
# #TOLERANCE OF 1 per mill FOR EACH ISOTOPE
# data_tol=subset(data_val_3,diff_C<1)
# data_tol2=subset(data_tol,diff_N<1)
# 
# #stats for constant sources
data_stats=data.frame(row.names = FALSE)
data_tol_2=subset(data_tol,data_tol$sources==1)
for (j in 1:nrow(conso)){
  data_tol_3=subset(data_tol_2,data_tol_2$val_ini==j)
for (i in 1:length(lambda_T)){
  val=lambda_T[i]
  subs_1=subset(data_tol_3,lambda_T==val)
  info=as.numeric(quantile(subs_1$B,probs=c(0.25,0.5,0.75)))
  data_stats=rbind(data_stats,data.frame(val_ini=j,lambda_T=val,Q1=info[1],Med=info[2],Q3=info[3]))
}#loop on lambda
}#loop on initial value

# col=c("red","blue","green")
# plot(NULL,xlim=c(0,6),ylim=c(0,100),xlab="λT",ylab="β (%)")
# for(i in 1:nrow(conso)) {
#   subs_2=subset(data_stats,val_ini==i)
#   points(subs_2$lambda_T,subs_2$Q1,col=col[i],type="p")
#   points(subs_2$lambda_T,subs_2$Q3,col=col[i],type="p")
#   points(subs_2$lambda_T,subs_2$Med,col=col[i],type="p",pch=16)
# }
# legend("topright",legend=c("0% switch","50% switch","100% switch"),col=col,pch=16)

# 
# # ### 5% BEST SOLUTION
# nb_sol=11
# best_sol_data=data.frame(row.names = NULL)
# stat_data=data.frame(row.names = NULL)
# for (k in 1:length(sources_list_var) ){ #loop sources
#   source_data=subset(data_dX_pred,data_dX_pred$sources==k)
#   for (j in 1:nrow(conso)){ #loop val ini
# # j=1
#   subs_1=subset(source_data,source_data$val_ini==j)
#   for (i in 1:length(lambda_T)){val=lambda_T[i] #loop λT
#   subs_2=subset(subs_1,subs_1$lambda_T==val)
#   best_sol=tail(subs_2,nb_sol)
#   stat=as.numeric(quantile(best_sol$B,probs=c(0.25,0.5,0.75)))
#   best_sol_data=rbind(best_sol_data,best_sol)
#   stat_data=rbind(stat_data,data.frame(source=k,val_ini=j,lambda_T=val,Q1=stat[1],M=stat[2],Q3=stat[3]))
#   }#loop λT
# } #loop val ini
# }#loop sources
# data_qu_on_garde=subset(best_sol_data,best_sol_data$val_ini==1)
# boxplot(diff_C~lambda_T,data=data_qu_on_garde)
# boxplot(diff_N~lambda_T,data=data_qu_on_garde)

# col=c("red","blue","green")
# plot(NULL,xlim=c(0,6),ylim=c(0,100),xlab="λT",ylab="β (%)")
# for (i in 1:nrow(conso)){
#   subs_stat=subset(stat_data,val_ini==i)
#   points(subs_stat$lambda_T,subs_stat$Q1,type="p",col=col[i])
#   points(subs_stat$lambda_T,subs_stat$Q3,type="p",col=col[i])
#   points(subs_stat$lambda_T,subs_stat$M,type="p",col=col[i],pch=16,cex=1.2)
# }
# 
# val_1=subset(best_sol_data,val_ini==3)
# boxplot(diff_C~lambda_T,data=val_1)
# boxplot(diff_N~lambda_T,data=val_1)

##     GRAPH  
deg_col=viridis(n=10, alpha = 1, begin = 0, end = 1, direction = -1, option = "D")
#diet1
colors=deg_col[c(1,5,10)]
eq=data.frame(I1=c(5),I2=c(10))
#diet2
# colors=deg_col[c(1,4,6)]
# eq=data.frame(I1=c(3.5),I2=c(4))

dev.off() # ferme les anciennes fenêtes
dev.new(width =  2250, height = 2625, unit = "px")
par(fig=c(0,0.7,0,1), new=TRUE,mar = c(4, 3.5, 4.1, 0.2))
plot(NULL,xlim=c(0,6),ylim=c(0,100),xlab="",ylab=""
     , las=1, cex.lab=1.2,font.lab=1, cex.axis=1.2,main="")
mtext("β (%)",side=2,line=4,padj=2,cex=1.2,font=1)
mtext("λT",side=1,line=1,padj=1.5,at=3,cex=1.2,font=1)
for (i in 1:nrow(conso)){
  subs_stats=subset(data_stats,val_ini==i)
  points(subs_stats$lambda_T,subs_stats$Q1,type="p",pch=3,col=colors[i])
  points(subs_stats$lambda_T,subs_stats$Q3,type="p",pch=3,col=colors[i])
}
for (i in 1:nrow(conso)){
  subs_stats=subset(data_stats,val_ini==i)
  points(subs_stats$lambda_T,subs_stats$M,type="p",pch=16,cex=1.2,col=colors[i])
}
par(fig=c(0.7,1,0,1), new=TRUE,mar = c(4, 3.5, 4.1, 0.2))
#poly sources
plot(NA, NA, ylim = c(0,10), xlim = c(0,10),
     xlab = "", ylab = "" ,las=1,main="",cex.axis=1.2)
mtext(expression(paste(" δ"^{15},"N (‰)")),side=2,line=4,padj=1.5,at=5,cex=1.2,font=2)
mtext(expression(paste(" δ"^{13},"C (‰)")),side=1,line=1,padj=1,at=5,cex=1.2,font=2)
points(c(0,5,10,0),c(0,10,5,0),type="l",col="black")
points(c(0,5,10,0),c(0,10,5,0),type="p",pch=16,col="black")
points(eq$I1,eq$I2,pch=17,cex=1.5)
points(conso$I1,conso$I2,pch=16,col=colors,cex=1.3)


#-------------------------------------------------------------------------------#

####Final script source effect
lambda_T=unique(data_SMM_dX_pred$lambda_T)
#for constant sources we keep the 0.5 tolerance
#SMM
data_tol_SMM_s1=data_SMM_dX_pred[data_SMM_dX_pred$sources==1&data_SMM_dX_pred$diff_C<0.5&data_SMM_dX_pred$diff_N<0.5,]
stat_data_SMM_s1=data.frame(row.names = NULL)
#SMM_delta
data_tol_SMM_delta_s1=data_SMM_delta_dX_pred[data_SMM_delta_dX_pred$sources==1&data_SMM_delta_dX_pred$diff_C<0.5&data_SMM_delta_dX_pred$diff_N<0.5,]
stat_data_SMM_delta_s1=data.frame(row.names = NULL)


for (j in 1:nrow(conso)){ #loop val ini
  # j=1
  #SMM
  subs_1=subset(data_tol_SMM_s1,val_ini==j)
  #SMM delta
  subs_1_bis=subset(data_tol_SMM_delta_s1,val_ini==j)
  for (i in 1:length(lambda_T)){val=lambda_T[i] #loop λT
  subs_2=subset(subs_1,subs_1$lambda_T==val)#SMM
  subs_2_bis=subset(subs_1_bis,subs_1_bis$lambda_T==val)#SMM delta
  stat_SMM=as.numeric(quantile(subs_2$B,probs=c(0.25,0.5,0.75)))
  stat_SMM_delta=as.numeric(quantile(subs_2_bis$B,probs=c(0.25,0.5,0.75)))
  stat_data_SMM_s1=rbind(stat_data_SMM_s1,data.frame(sources=1,val_ini=j,lambda_T=val
                                                     ,Q1=stat_SMM[1],M=stat_SMM[2],Q3=stat_SMM[3]))
  stat_data_SMM_delta_s1=rbind(stat_data_SMM_delta_s1,data.frame(sources=1,val_ini=j,lambda_T=val
                                                                 ,Q1=stat_SMM_delta[1],M=stat_SMM_delta[2],Q3=stat_SMM_delta[3]))
  }#loop λT
} #loop val ini

#for non-constant sources we keep the 5% best solutions
# # ### 5% BEST SOLUTION
nb_sol=11
#SMM 
best_sol_SMM=data.frame(row.names = NULL)
stat_data_SMM=data.frame(row.names = NULL)
#SMM delta
best_sol_SMM_delta=data.frame(row.names = NULL)
stat_data_SMM_delta=data.frame(row.names = NULL)
for (k in 2:length(sources_list_var) ){ #loop sources from 2 to 5
  data_SMM=subset(data_SMM_dX_pred,sources==k)
  data_SMM_delta=subset(data_SMM_delta_dX_pred,sources==k)
  for (j in 1:nrow(conso)){ #loop val ini
# j=1
  subs_1=subset(data_SMM,data_SMM$val_ini==j) #SMM
  subs_1_bis=subset(data_SMM_delta,data_SMM$val_ini==j) #SMM delta
  
  for (i in 1:length(lambda_T)){val=lambda_T[i] #loop λT
  subs_2=subset(subs_1,subs_1$lambda_T==val)
  subs_2_bis=subset(subs_1_bis,subs_1_bis$lambda_T==val)
  best_sol_1=tail(subs_2,nb_sol)
  best_sol_2=tail(subs_2_bis,nb_sol)
  stat_1=as.numeric(quantile(best_sol_1$B,probs=c(0.25,0.5,0.75)))
  stat_2=as.numeric(quantile(best_sol_2$B,probs=c(0.25,0.5,0.75)))
  best_sol_SMM=rbind(best_sol_SMM,best_sol_1)
  best_sol_SMM_delta=rbind(best_sol_SMM_delta,best_sol_2)
  stat_data_SMM=rbind(stat_data_SMM,data.frame(sources=k,val_ini=j,lambda_T=val,Q1=stat_1[1],M=stat_1[2],Q3=stat_1[3]))
  stat_data_SMM_delta=rbind(stat_data_SMM_delta,data.frame(sources=k,val_ini=j,lambda_T=val,Q1=stat_2[1],M=stat_2[2],Q3=stat_2[3]))
  }#loop λT
} #loop val ini
}#loop sources

##JOINING THE TWO DATAFRAMES
final_data_set_SMM=rbind(stat_data_SMM_s1,stat_data_SMM)
final_data_set_SMM_delta=rbind(stat_data_SMM_delta_s1,stat_data_SMM_delta)
##GRAPH SWITCH EFFECT
# colors=c("#000000","#E69F00" ,"#56B4E9","#009E73","#0072B2" )
# plot(NULL,xlim=c(0,6),ylim=c(0,100),xlab="λT",ylab="β (%)"
#      , las=1, cex.lab=1.4,font.lab=2, cex.axis=1.4,main="")
# val_data_1=subset(final_data_set,final_data_set$val_ini==1)
# 
# for (i in 1:length(sources_list_var)){
#   subs_stats=subset(val_data_1,val_data_1$source==i)
#   points(subs_stats$lambda_T,subs_stats$Q1,type="p",pch=3,col=colors[i])
#   points(subs_stats$lambda_T,subs_stats$Q3,type="p",pch=3,col=colors[i])
# }
# for (i in 1:length(sources_list_var)){
#   subs_stats=subset(val_data_1,val_data_1$source==i)
#   points(subs_stats$lambda_T,subs_stats$M,type="p",pch=16,cex=1.2,col=colors[i])
# }



#####source effect
process_source_effect <- function(data) {
  extr_values_ordered=NULL
  for (i in 1: 3) {
    extr_data=subset(data,data$val_ini==i)
    val_cst=subset(extr_data,extr_data$sources==1)
    val_2=subset(extr_data,extr_data$sources==2)
    val_3=subset(extr_data,extr_data$sources==3)
    val_4=subset(extr_data,extr_data$sources==4)
    val_5=subset(extr_data,extr_data$sources==5)
    
    comp=data.frame("1"=val_cst$M,"2"=val_2$M,"3"=val_3$M,"4"=val_4$M,"5"=val_5$M)
    max=apply(comp,1,max)
    min=apply(comp,1,min)
    
    extr_values=data.frame(lambda_T=val_cst$lambda_T,min,max)
    in_order=order(extr_values$lambda_T)
    extr_values_ordered=c(extr_values_ordered,list(extr_values[in_order,]))
  }
  return(extr_values_ordered)
}

#-------------------------GRAPHS--------------------------------#
####graph source effect
# data_SMM=final_data_set_SMM
graph_source_effect <- function(data_SMM,data_SMM_delta,nb_values,col_val) {
  
  #GENERATING THE COLOR PALETTE
  deg_col=viridis(n=nb_values, alpha = 1, begin = 0, end = 1, direction = -1, option = "D")
  legend_image <- as.raster(matrix(rev(deg_col), ncol=1))
  
  #GENERATING THE COLOR RANK OF CHOOSEN INITIAL VALUES

  #process SMM
  process_SMM=process_source_effect(data_SMM)
  low_dist_SMM=process_SMM[[1]]
  med_low_dist_SMM=process_SMM[[2]]
  med_high_dist_SMM=process_SMM[[3]]
  #process SMM_delta
  process_SMM_delta=process_source_effect(data_SMM_delta)
  low_dist_SMM_delta=process_SMM_delta[[1]]
  med_low_dist_SMM_delta=process_SMM_delta[[2]]
  med_high_dist_SMM_delta=process_SMM_delta[[3]]
  
  
  #cstes source values
  data_cst=subset(data_SMM,data_SMM$sources==1)
  ref_val=NULL
  for (i in 1:length(unique(data_cst$val_ini))){
    subs=subset(data_cst,data_cst$val_ini==i)
    ref_val=c(ref_val,list(subs$M))}
  #GRAPH
  windows()
  #if switch legend
  # par(mfrow=c(2,3))
  
  # #generating the legend
  # plot(c(0,2),c(0,1),type = 'n', axes = F,xlab = '', ylab = '', main = 'Switch (%)')
  # # Tracer les graduations de la légende, ici le paramétrage en est strictement manuel 5 paliers entre 0 et 100
  # text(x=1.5, y = seq(0,1,l=5), labels = seq(0,100,l=5))
  # # ADD GRADIENT IMAGE
  # rasterImage(legend_image, 0, 0, 1,1)
  # 
  
  par(mfrow=c(2,2))
  #graph_1
  plot(NA, NA, ylim = c(0,100), xlim = c(0,6),
       xlab = c("λT"), ylab =c("β (%)") , las=1, cex.lab=1.4,font.lab=1, cex.axis=1.4,main="")
  polygon(c(low_dist_SMM$lambda_T,rev(low_dist_SMM$lambda_T)),c(low_dist_SMM$min,rev(low_dist_SMM$max)),col="grey",border=NA)
  polygon(c(low_dist_SMM_delta$lambda_T,rev(low_dist_SMM_delta$lambda_T)),c(low_dist_SMM_delta$min,rev(low_dist_SMM_delta$max)),col="grey50",border=NA)
  points(unique(data_cst$lambda_T),ref_val[[1]],type="p",col=deg_col[col_val[1]],pch=16,cex=1.4)
  legend("topright",legend=c("SMM","SMMΔ"),col=c("grey","grey50"),lty=1,cex=1.2)
  
  #graph_2
  plot(NA, NA, ylim = c(0,100), xlim = c(0,6),
       xlab = c("λT"), ylab =c("β (%)") , las=1, cex.lab=1.4,font.lab=1, cex.axis=1.4,main="")
  polygon(c(med_low_dist_SMM$lambda_T,rev(med_low_dist_SMM$lambda_T)),c(med_low_dist_SMM$min,rev(med_low_dist_SMM$max)),col="grey",border=NA)
  polygon(c(med_low_dist_SMM_delta$lambda_T,rev(med_low_dist_SMM_delta$lambda_T)),c(med_low_dist_SMM_delta$min,rev(med_low_dist_SMM_delta$max)),col="grey50",border=NA)
  points(unique(data_cst$lambda_T),ref_val[[2]],type="p",col=deg_col[col_val[2]],pch=16,cex=1.4)
  
  
  # #graph_3
  plot(NA, NA, ylim = c(0,100), xlim = c(0,6),
       xlab = c("λT"), ylab =c("β (%)") , las=1, cex.lab=1.4,font.lab=1, cex.axis=1.4,main="")
  polygon(c(med_high_dist_SMM$lambda_T,rev(med_high_dist_SMM$lambda_T)),c(med_high_dist_SMM$min,rev(med_high_dist_SMM$max)),col="grey",border=NA)
  polygon(c(med_high_dist_SMM_delta$lambda_T,rev(med_high_dist_SMM_delta$lambda_T)),c(med_high_dist_SMM_delta$min,rev(med_high_dist_SMM_delta$max)),col="grey50",border=NA)
  points(unique(data_cst$lambda_T),ref_val[[3]],type="p",col=deg_col[col_val[3]],pch=16,cex=1.4)
  

}

#diet 1
# col_val=c(1,5,10)
#diet 2
col_val=c(1,4,6)

graph_source_effect(data_SMM=final_data_set_SMM,data_SMM_delta=final_data_set_SMM_delta,nb_values=10,col_val=col_val)

#-----------------------------------In silico experiment (2I,4S)--------------------------------------

#CHANGE HERE IF THE SOURCE NUMBER OR ISOTOPE NUMBER CHANGE

TEF=list(data.frame("s1"=c(0),"s2"=c(0),"s3"=c(0),"s4"=c(0)),data.frame("s1"=c(0),"s2"=c(0),"s3"=c(0),"s4"=c(0))) #the TEF=0 for all sources and isotope
conc=list(data.frame("s1"=c(1),"s2"=c(1),"s3"=c(1),"s4"=c(1)),data.frame("s1"=c(1),"s2"=c(1),"s3"=c(1),"s4"=c(1))) #the conc=1 for all sources and isotope
iter=0.05


#CONSUMER INITAL VALUE
nb_ini_val=20

initial_values=read.csv2("C:/Users/Etudiant/Documents/thèse/1ere année/1er papier/finals/1_st paper/Data/val_ini_2.csv",row.names = NULL)
resum_data_SMM=data.frame(B=c(),l_T=c(),sources=c(),val_ini=c(),DC_obs=c(),DN_obs=c())
resum_data_SMM_delta=data.frame(B=c(),l_T=c(),sources=c(),val_ini=c())



#DIET
#the consumer eats 100% of source 3 diet 4
# combi=data.frame(s1=c(0),s2=c(0),s3=c(1),s4=c(0))
# conso=initial_values[c(12,14,20),]
#theconsumer eats a mix and the solution value is diet 3
combi=data.frame(s1=c(0.2),s2=c(0.2),s3=c(0.4),s4=c(0.2))
conso=initial_values[c(9,13,20),]
#keeping estimated DX values 
DX=data.frame(lambda_T=c(),t_f=c(),DC=c(),DN=c())

# #ISOTOPIC SPACE
# plot(NA, NA, ylim = c(0,10), xlim = c(0,10),
#      xlab = "", ylab = "" ,las=1,main="b")
# mtext(expression(paste(" δ"^{15},"N (‰)")),side=2,line=4,padj=1.5,at=5,cex=1.6,font=2)
# mtext(expression(paste(" δ"^{13},"C (‰)")),side=1,line=1,padj=1.5,at=5,cex=1.6,font=2)
# points(conso$I1,conso$I2,pch=16)
# points(c(0,2,5,10,0),c(0,6,10,5,0),type="l",col="red")
# points(c(0,2,5,10,0),c(0,6,10,5,0),type="p",pch=16,col="red")
# text(x=c(0.5,2.5,5.5,9.5),y=c(0,6,10,5),labels = c("2","4","3","1"),cex=1.4,font=2,col="red")
# 

p=seq(2,300,length=4) #period
l=seq(0.002,0.02,length=5) #lambda

for (s in 1:length(l)){ #for every lambda
   # s=1
  lambda=l[s]
  print(c("l",s))
  
  for (f in 1:length(p)){ #for every period
    # f=1
    print(c("p",f))
    period=round(p[f])
    t <- seq(0,period,1) # output times for the ODE (d)
    liste_lambda=list(data.frame(Date=t,l=rep(lambda,length(t))),data.frame(Date=t,l=rep(lambda,length(t))))
    #SOURCES VALUES WITH BIGGER VARIATION
    sources_list_var=c(list(list(data.frame(Date=t,s1=rep(10,period+1),s2=rep(0,period+1),s3=rep(5,period+1),s4=rep(2,period+1)), #pattern a
                                 data.frame(Date=t,s1=rep(5,period+1),s2=rep(0,period+1),s3=rep(10,period+1),s4=rep(6,period+1))))
                       ,
                       
                       list(list(data.frame(Date=t,s1=rep(10,period+1),s2=seq(-5,0,length=period+1),s3=seq(15,3,length=period+1),s4=rep(2,period+1)), #pattern b
                                 data.frame(Date=t,s1=rep(5,period+1),s2=rep(0,period+1),s3=rep(10,period+1),s4=rep(6,period+1)))),
                       
                       list(list(data.frame(Date=t,s1=rep(10,period+1),s2=seq(0,-5,length=period+1),s3=seq(0,-20,length=period+1),s4=rep(2,period+1)), #pattern c
                                 data.frame(Date=t,s1=rep(5,period+1),s2=rep(0,period+1),s3=rep(10,period+1),s4=rep(6,period+1)))),
                       
                       list(list(data.frame(Date=t,s1=rep(10,period+1),s2=seq(5,10,length=period+1),s3=seq(0,20,length=period+1),s4=rep(2,period+1)), #pattern d
                                 data.frame(Date=t,s1=rep(5,period+1),s2=rep(0,period+1),s3=rep(10,period+1),s4=rep(6,period+1)))),
                       
                       list(list(data.frame(Date=t,s1=rep(10,period+1),s2=seq(0,-5,length=period+1),s3=seq(15,0,length=period+1),s4=rep(2,period+1)), #pattern e
                                 data.frame(Date=t,s1=rep(5,period+1),s2=rep(0,period+1),s3=rep(10,period+1),s4=rep(6,period+1))))
    )
    
    
    
    
    
    for(m in 1:nrow(conso)){ #LOOP ON CONSUMER
     # m=3
    conso_ini=conso[m,]
    print(m)
    
    indicator_value=NULL
    bias_SMM=NULL
    bias_SMM_delta=NULL
    bias_DMM=NULL
    dist_KL_SMM=NULL
    dist_KL_SMM_delta=NULL
    SMM_arch=data.frame(row.names = NULL)
    SMM_delta_arch=data.frame(row.names = NULL)
    DMM_arch=data.frame(row.names = NULL)
    
    for (n in 1:length(sources_list_var)) { #LOOP ON SOURCES
    # n=1
    sources_list=sources_list_var[[n]]
    
    #CREATING THE SOURCE APPROXFUN
    sources_approxfun=list(c(approxfun(sources_list[[1]]$Date,sources_list[[1]]$s1,method='linear',rule=2)
                             ,approxfun(sources_list[[1]]$Date,sources_list[[1]]$s2,method='linear',rule=2)
                             ,approxfun(sources_list[[1]]$Date,sources_list[[1]]$s3,method='linear',rule=2)
                             ,approxfun(sources_list[[1]]$Date,sources_list[[1]]$s4,method='linear',rule=2)),
                           c(approxfun(sources_list[[2]]$Date,sources_list[[2]]$s1,method='linear',rule=2)
                             ,approxfun(sources_list[[2]]$Date,sources_list[[2]]$s2,method='linear',rule=2)
                             ,approxfun(sources_list[[2]]$Date,sources_list[[2]]$s3,method='linear',rule=2)
                             ,approxfun(sources_list[[2]]$Date,sources_list[[2]]$s4,method='linear',rule=2)))
    
    #EXPLAINS HOW MANY ISO,SOURCES AND OBS
    data_feat=data_features(sources_list )
    res=NULL
    new_conso_value=NULL
    #GENERATING DX AND DXinf VALUES
    for (k in 1:data_feat[1]) #FOR EACH ISOTOPE
    {res=c(res,list(RUN_TIMdyn2(t,#time
                                state0=c(X=conso_ini[[k]]),#initial conditions
                                par=c(#two parameters lambda and Xinf i.e. forcing functions 
                                  lambda = approxfun(liste_lambda[[k]]$Date,liste_lambda[[k]]$l,method = 'const',rule=2), 
                                  Xinf=approxfun(t, Mixing(t=t,sources_approx=sources_approxfun[[k]],p=combi,TEF=TEF[[k]]
                                                           ,Conc=conc[[k]])))))) #COMPUTING THE DX AND DXinf
    
    
    new_conso_value=c(new_conso_value,list(data.frame(t=seq(0,period,1),DX=res[[k]][[2]])))
    }
    
    
    #USING DX VALUES AS DATA, estimating the SMM and SMM delta
    
    
    #GENERATING THE POSSIBLE SOLUTIONS
    sol<- expand.grid(rep(list(seq(0,1,iter)),data_feat[2])) #data_feat[2] is the number of food sources
    sol <- subset(sol,rowSums(sol)==1)
    
    
    #SMM and SMM_delta estimation at the last "t"
    val_conso=NULL
    val_0=NULL
    for (j in 1:data_feat[[1]]) {val_conso=c(val_conso,list(new_conso_value[[j]][period,]))
    val_0=c(val_0,list(new_conso_value[[j]][1,]))}
    
    
    #KEEPING DX VALUES AT THE LAST "t"
    DX=rbind(DX,data.frame(lambda_T=lambda*period,val_ini=m,t_f=val_conso[[1]]$t,DC=val_conso[[1]]$DX,
                           DN=val_conso[[2]]$DX))
    
    
    SMM=SMM_func(t=period,sources_approx=sources_approxfun,list_TEF=TEF,details_conc=conc
                 ,combi=sol,data_feat=data_feat,details_consu=val_conso,X=nrow(sol))
    SMM_delta=SMM_delta_func(t=period,sources_approx=sources_approxfun,list_TEF=TEF,details_conc=conc,combi=sol
                             ,data_feat=data_feat,details_consu=val_conso,X=nrow(sol),lambda_list=liste_lambda)

    
    
    #BIAS TT SYSTEME
    resum_data_SMM=rbind(resum_data_SMM,data.frame(B=(abs(combi[[1]]-SMM[[2]])+abs(combi[[2]]-SMM[[3]])+abs(combi[[3]]-SMM[[4]])+abs(combi[[4]]-SMM[[5]]))*50
                                                   ,lambda_T=c(lambda*period)
                                                   ,sources=c(n),val_ini=c(m),SMM,DC=val_conso[[1]]$DX,DN=val_conso[[2]]$DX
                                                   ,s1_C=sources_approxfun[[1]][[1]](period),s2_C=sources_approxfun[[1]][[2]](period),s3_C=sources_approxfun[[1]][[3]](period)
                                                   ,s4_C=sources_approxfun[[1]][[4]](period),s1_N=sources_approxfun[[2]][[1]](period),s2_N=sources_approxfun[[2]][[2]](period)
                                                   ,s3_N=sources_approxfun[[2]][[3]](period),s4_N=sources_approxfun[[2]][[4]](period)))
    resum_data_SMM_delta=rbind(resum_data_SMM_delta,data.frame(B=(abs(combi[[1]]-SMM_delta[[2]])+abs(combi[[2]]-SMM_delta[[3]])+abs(combi[[3]]-SMM_delta[[4]])+abs(combi[[4]]-SMM_delta[[5]]))*50,lambda_T=c(lambda*period)
                                                               ,sources=c(n),val_ini=c(m),SMM_delta,DC=val_conso[[1]]$DX,DN=val_conso[[2]]$DX
                                                               ,s1_C=sources_approxfun[[1]][[1]](period),s2_C=sources_approxfun[[1]][[2]](period),s3_C=sources_approxfun[[1]][[3]](period)
                                                               ,s4_C=sources_approxfun[[1]][[4]](period),s1_N=sources_approxfun[[2]][[1]](period),s2_N=sources_approxfun[[2]][[2]](period)
                                                               ,s3_N=sources_approxfun[[2]][[3]](period),s4_N=sources_approxfun[[2]][[4]](period)))
    }#loop on sources
    }#loop on ini value
    
  }#loop on period
}#loop on lambda




lambda_T=unique(resum_data_SMM$lambda_T)


#Equilibrium is c'4.4,6.2)



diff_C_SMM=NULL
diff_N_SMM=NULL
Dc_pred_SMM=NULL
Dn_pred_SMM=NULL
diff_C_SMM_delta=NULL
diff_N_SMM_delta=NULL
Dc_pred_SMM_delta=NULL
Dn_pred_SMM_delta=NULL
for (i in 1: nrow(resum_data_SMM)) {
  #SMM
  dc_SMM=resum_data_SMM$s1_C[i]*resum_data_SMM$Var1[i]+resum_data_SMM$s2_C[i]*resum_data_SMM$Var2[i]+resum_data_SMM$s3_C[i]*resum_data_SMM$Var3[i]+resum_data_SMM$s4_C[i]*resum_data_SMM$Var4[i]
  dn_SMM=resum_data_SMM$s1_N[i]*resum_data_SMM$Var1[i]+resum_data_SMM$s2_N[i]*resum_data_SMM$Var2[i]+resum_data_SMM$s3_N[i]*resum_data_SMM$Var3[i]+resum_data_SMM$s4_N[i]*resum_data_SMM$Var4[i]
  Dc_pred_SMM=c(Dc_pred_SMM,dc_SMM)
  Dn_pred_SMM=c(Dn_pred_SMM,dn_SMM)
  diff_C_SMM=c(diff_C_SMM,abs(resum_data_SMM$DC[i]-dc_SMM))
  diff_N_SMM=c(diff_N_SMM,abs(resum_data_SMM$DN[i]-dn_SMM))
  #SMM_delta
  dc_SMM_delta=resum_data_SMM_delta$s1_C[i]*resum_data_SMM_delta$Var1[i]+resum_data_SMM_delta$s2_C[i]*resum_data_SMM_delta$Var2[i]+resum_data_SMM_delta$s3_C[i]*resum_data_SMM_delta$Var3[i]+resum_data_SMM_delta$s4_C[i]*resum_data_SMM_delta$Var4[i]
  dn_SMM_delta=resum_data_SMM_delta$s1_N[i]*resum_data_SMM_delta$Var1[i]+resum_data_SMM_delta$s2_N[i]*resum_data_SMM_delta$Var2[i]+resum_data_SMM_delta$s3_N[i]*resum_data_SMM_delta$Var3[i]+resum_data_SMM_delta$s4_N[i]*resum_data_SMM_delta$Var4[i]
  Dc_pred_SMM_delta=c(Dc_pred_SMM_delta,dc_SMM_delta)
  Dn_pred_SMM_delta=c(Dn_pred_SMM_delta,dn_SMM_delta)
  diff_C_SMM_delta=c(diff_C_SMM_delta,abs(resum_data_SMM_delta$DC[i]-dc_SMM_delta))
  diff_N_SMM_delta=c(diff_N_SMM_delta,abs(resum_data_SMM_delta$DN[i]-dn_SMM_delta))
}

#adding it to the dataset
data_SMM_dX_pred=cbind(resum_data_SMM,Dc_pred_SMM,Dn_pred_SMM,diff_C_SMM,diff_N_SMM)
data_SMM_delta_dX_pred=cbind(resum_data_SMM_delta,Dc_pred_SMM_delta,Dn_pred_SMM_delta
                             ,diff_C_SMM_delta,diff_N_SMM_delta)

# write.csv2(data_SMM_dX_pred,file="diet_3_SMM_values.csv",row.names = FALSE)
# write.csv2(data_SMM_delta_dX_pred,file="diet_3_SMM_delta_values.csv",row.names = FALSE)
data_SMM_dX_pred=read.csv2("diet_4_SMM_values.csv")
data_SMM_delta_dX_pred=read.csv2("diet_4_SMM_delta_values.csv")
#-----------------WORK----------------#
# # distance_bis=data.frame(row.names = NULL)
# # for (i in 1:nrow(data_dX_pred)){
# #     i=1
# #   diff=c((data_dX_pred$diff_C[i]/data_dX_pred$DC_obs[i])^2
# #          ,(data_dX_pred$diff_N[i]/data_dX_pred$DN_obs[i])^2)
# #   distance_bis=rbind(distance_bis,sum(diff))
# # }
# 
# 
# 
# # #SWITCH EFFECT GRAPH
# 
# ###KEEPING THE 1% BEST SOLUTIONS (17)
# #KEEPING THE 2% best solution (33)
# #KEEPING THE 5% best solution (83)
# # nb_sol=17
# # best_sol_data=data.frame(row.names = NULL)
# # stat_data=data.frame(row.names = NULL)
# # for (j in 1:nrow(conso)){ #loop val ini
# # # j=1
# #   subs_1=subset(data_dX_pred,data_dX_pred$val_ini==j)
# #   for (i in 1:length(lambda_T)){val=lambda_T[i] #loop λT
# #   subs_2=subset(subs_1,subs_1$lambda_T==val)
# #   best_sol=tail(subs_2,nb_sol)
# #   stat=as.numeric(quantile(best_sol$B,probs=c(0.25,0.5,0.75)))
# #   best_sol_data=rbind(best_sol_data,best_sol)
# #   stat_data=rbind(stat_data,data.frame(val_ini=j,lambda_T=val,Q1=stat[1],M=stat[2],Q3=stat[3]))
# #   }#loop λT
# # } #loop val ini
# # col=c("red","blue","green")
# # plot(NULL,xlim=c(0,6),ylim=c(0,100),xlab="λT",ylab="β (%)")
# # for (i in 1:nrow(conso)){
# #   subs_stat=subset(stat_data,val_ini==i)
# #   points(subs_stat$lambda_T,subs_stat$Q1,type="p",col=col[i])
# #   points(subs_stat$lambda_T,subs_stat$Q3,type="p",col=col[i])
# #   points(subs_stat$lambda_T,subs_stat$M,type="p",col=col[i],pch=16,cex=1.2)
# # }
# # 
# # val_1=subset(best_sol_data,val_ini==3)
# # boxplot(diff_C~lambda_T,data=val_1)
# # boxplot(diff_N~lambda_T,data=val_1)
# 
# ###KEEPING A TOLERANCE 
# 
# # tol_data=subset(data_dX_pred,data_dX_pred$diff_C<2)
# # tol_data_2=subset(tol_data,tol_data$diff_N<4)
# tol_data_3=data_dX_pred[data_dX_pred$diff_C<0.2&data_dX_pred$diff_N<0.2,]
# stat_data=data.frame(row.names = NULL)
# df_nb_sol=data.frame(row.names = NULL)
# for (i in 1:nrow(conso)) {#LOOP on val_ini
#   subs_1=subset(tol_data_3,tol_data_3$val_ini==i)
#   for (j in 1:length(lambda_T)) { #LOOP ON λT
#     val=lambda_T[j]
#     subs_2=subset(subs_1,subs_1$lambda_T==val)
#     df_nb_sol=rbind(df_nb_sol,data.frame(val_ini=i,λT=val,nb_sol=nrow(subs_2)))
#     stat=as.numeric(quantile(subs_2$B,probs=c(0.25,0.5,0.75)))
#     stat_data=rbind(stat_data,data.frame(val_ini=i,lambda_T=val,Q1=stat[1],M=stat[2],Q3=stat[3]))
#   }#LOOP ON λT
# }#LOOP on val_ini
# 
# col=c("red","blue","green")
# plot(NULL,xlim=c(0,6),ylim=c(0,100),xlab="λT",ylab="β (%)")
# for (i in 1:nrow(conso)){
#   subs_stat=subset(stat_data,val_ini==i)
#   points(subs_stat$lambda_T,subs_stat$Q1,type="p",col=col[i])
#   points(subs_stat$lambda_T,subs_stat$Q3,type="p",col=col[i])
# }
# for (i in 1:nrow(conso)){
#   subs_stat=subset(stat_data,val_ini==i)
#   points(subs_stat$lambda_T,subs_stat$M,type="p",col=col[i],pch=16,cex=1.2)
# }

#--------------------------REAL CODE-----------------------------#

#---------------SWITCH-----------------# 
#CHECKING THE NUMBER OF SOLUTIONS
lambda_T=unique(data_SMM_dX_pred$lambda_T)
# data_tol=data_dX_pred[data_dX_pred$sources==5&data_dX_pred$diff_C<1&data_dX_pred$diff_N<1,]
data_tol=data_SMM_dX_pred[data_SMM_dX_pred$diff_C<0.5&data_SMM_dX_pred$diff_N<0.5,]
nb_sources=5
df_nb_sol=data.frame(row.names = FALSE)
for (i in 1:nb_sources){#loop on sources
  df_sources=subset(data_tol,data_tol$sources==i)
  for (j in 1:nrow(conso)){#loop on val ini
    df_val_ini=subset(df_sources,df_sources$val_ini==j)
    for (k in 1:length(lambda_T)){#loop on lambdaT
      val=lambda_T[k]
      df_lambda=subset(df_val_ini,df_val_ini$lambda_T==val)
      df_nb_sol=rbind(df_nb_sol,data.frame(source=i,val_ini=j,lambda_T=val,nb_sol=nrow(df_lambda)))
    }#loop on lambdaT
  }#loop on val ini
}#loop on sources



# #stats for constant sources
data_stats=data.frame(row.names = FALSE)
data_tol_2=subset(data_tol,data_tol$sources==1)
for (j in 1:nrow(conso)){
  data_tol_3=subset(data_tol_2,data_tol_2$val_ini==j)
  for (i in 1:length(lambda_T)){
    val=lambda_T[i]
    subs_1=subset(data_tol_3,lambda_T==val)
    info=as.numeric(quantile(subs_1$B,probs=c(0.25,0.5,0.75)))
    data_stats=rbind(data_stats,data.frame(val_ini=j,lambda_T=val,Q1=info[1],Med=info[2],Q3=info[3]))
  }#loop on lambda
}#loop on initial value


##     GRAPH  
deg_col=viridis(n=10, alpha = 1, begin = 0, end = 1, direction = -1, option = "D")
legend_image <- as.raster(matrix(rev(deg_col), ncol=1))
rasterImage(legend_image, 0, 0, 0.8,1)
plot(c(0,2),c(0,1),type = 'n', axes = F,xlab = '', ylab = '', main = '')
# Tracer les graduations de la légende, ici le paramétrage en est strictement manuel 5 paliers entre 0 et 1
text(x=1, y = seq(0,1,l=5), labels = seq(0,100,l=5))
# Ajouter l'image de gradient
rasterImage(legend_image, 0, 0, 0.6,1)
#diet1
# colors=deg_col[c(1,5,10)]
#diet2
 colors=deg_col[c(2,3,8)]

# plot(NULL,xlim=c(0,6),ylim=c(0,100),xlab="λT",ylab="β (%)"
#      , las=1, cex.lab=1.4,font.lab=2, cex.axis=1.4,main="")
# for (i in 1:nrow(conso)){
#   subs_stats=subset(data_stats,val_ini==i)
#   points(subs_stats$lambda_T,subs_stats$Q1,type="p",pch=3,col=colors[i])
#   points(subs_stats$lambda_T,subs_stats$Q3,type="p",pch=3,col=colors[i])
# }
# for (i in 1:nrow(conso)){
#   subs_stats=subset(data_stats,val_ini==i)
#   points(subs_stats$lambda_T,subs_stats$M,type="p",pch=16,cex=1.2,col=colors[i])
# }

####Copy de au dessus


##     GRAPH  
deg_col=viridis(n=10, alpha = 1, begin = 0, end = 1, direction = -1, option = "D")
#diet3 eq 4 6.2
colors=deg_col[c(2,3,8)]
eq=data.frame(I1=c(4),I2=c(6.2))
# diet4
# colors=deg_col[c(1,5,10)]
# eq=data.frame(I1=c(5),I2=c(10))


dev.off() # ferme les anciennes fenêtes
dev.new(width =  2250, height = 2625, unit = "px")
par(fig=c(0,0.7,0,1), new=TRUE,mar = c(4, 3.5, 4.1, 0.2))
plot(NULL,xlim=c(0,6),ylim=c(0,100),xlab="",ylab=""
     , las=1, cex.lab=1.4,font.lab=2, cex.axis=1.4,main="")
mtext("β (%)",side=2,line=4,padj=2,cex=1.2,font=1)
mtext("λT",side=1,line=1,padj=1.5,at=3,cex=1.2,font=1)
for (i in 1:nrow(conso)){
  subs_stats=subset(data_stats,val_ini==i)
  points(subs_stats$lambda_T,subs_stats$Q1,type="p",pch=3,col=colors[i])
  points(subs_stats$lambda_T,subs_stats$Q3,type="p",pch=3,col=colors[i])
}
for (i in 1:nrow(conso)){
  subs_stats=subset(data_stats,val_ini==i)
  points(subs_stats$lambda_T,subs_stats$M,type="p",pch=16,cex=1.2,col=colors[i])
}

par(fig=c(0.7,1,0,1), new=TRUE,mar = c(4, 3.5, 4.1, 0.2))
#poly sources
plot(NA, NA, ylim = c(0,10), xlim = c(0,10),
     xlab = "", ylab = "" ,las=1,main="",cex.axis=1.2)
mtext(expression(paste(" δ"^{15},"N (‰)")),side=2,line=4,padj=1.5,at=5,cex=1.2,font=2)
mtext(expression(paste(" δ"^{13},"C (‰)")),side=1,line=1,padj=1,at=5,cex=1.2,font=2)
points(c(0,2,5,10,0),c(0,6,10,5,0),type="l")
points(c(0,2,5,10,0),c(0,6,10,5,0),type="p",pch=16,)
points(eq$I1,eq$I2,pch=17,cex=1.5)
points(conso$I1,conso$I2,pch=16,col=colors,cex=1.3)

#------------------SOURCE EFFECT-----------------#

####Final script source effect
lambda_T=unique(data_SMM_dX_pred$lambda_T)
#for constant sources we keep the 0.5 tolerance
#SMM
data_tol_SMM_s1=data_SMM_dX_pred[data_SMM_dX_pred$sources==1&data_SMM_dX_pred$diff_C<0.5&data_SMM_dX_pred$diff_N<0.5,]
stat_data_SMM_s1=data.frame(row.names = NULL)
#SMM_delta
data_tol_SMM_delta_s1=data_SMM_delta_dX_pred[data_SMM_delta_dX_pred$sources==1&data_SMM_delta_dX_pred$diff_C<0.5&data_SMM_delta_dX_pred$diff_N<0.5,]
stat_data_SMM_delta_s1=data.frame(row.names = NULL)


for (j in 1:nrow(conso)){ #loop val ini
  # j=1
  #SMM
  subs_1=subset(data_tol_SMM_s1,val_ini==j)
  #SMM delta
  subs_1_bis=subset(data_tol_SMM_delta_s1,val_ini==j)
  for (i in 1:length(lambda_T)){val=lambda_T[i] #loop λT
  subs_2=subset(subs_1,subs_1$lambda_T==val)#SMM
  subs_2_bis=subset(subs_1_bis,subs_1_bis$lambda_T==val)#SMM delta
  stat_SMM=as.numeric(quantile(subs_2$B,probs=c(0.25,0.5,0.75)))
  stat_SMM_delta=as.numeric(quantile(subs_2_bis$B,probs=c(0.25,0.5,0.75)))
  stat_data_SMM_s1=rbind(stat_data_SMM_s1,data.frame(sources=1,val_ini=j,lambda_T=val
                                                     ,Q1=stat_SMM[1],M=stat_SMM[2],Q3=stat_SMM[3]))
  stat_data_SMM_delta_s1=rbind(stat_data_SMM_delta_s1,data.frame(sources=1,val_ini=j,lambda_T=val
                                                                 ,Q1=stat_SMM_delta[1],M=stat_SMM_delta[2],Q3=stat_SMM_delta[3]))
  }#loop λT
} #loop val ini

#for non-constant sources we keep the 5% best solutions
# # ### 5% BEST SOLUTION
nb_sol=83
#SMM 
best_sol_SMM=data.frame(row.names = NULL)
stat_data_SMM=data.frame(row.names = NULL)
#SMM delta
best_sol_SMM_delta=data.frame(row.names = NULL)
stat_data_SMM_delta=data.frame(row.names = NULL)
for (k in 2:5 ){ #loop sources from 2 to 5
  data_SMM=subset(data_SMM_dX_pred,sources==k)
  data_SMM_delta=subset(data_SMM_delta_dX_pred,sources==k)
  for (j in 1:nrow(conso)){ #loop val ini
# j=1
  subs_1=subset(data_SMM,data_SMM$val_ini==j) #SMM
  subs_1_bis=subset(data_SMM_delta,data_SMM$val_ini==j) #SMM delta
  
  for (i in 1:length(lambda_T)){val=lambda_T[i] #loop λT
  subs_2=subset(subs_1,subs_1$lambda_T==val)
  subs_2_bis=subset(subs_1_bis,subs_1_bis$lambda_T==val)
  best_sol_1=tail(subs_2,nb_sol)
  best_sol_2=tail(subs_2_bis,nb_sol)
  stat_1=as.numeric(quantile(best_sol_1$B,probs=c(0.25,0.5,0.75)))
  stat_2=as.numeric(quantile(best_sol_2$B,probs=c(0.25,0.5,0.75)))
  best_sol_SMM=rbind(best_sol_SMM,best_sol_1)
  best_sol_SMM_delta=rbind(best_sol_SMM_delta,best_sol_2)
  stat_data_SMM=rbind(stat_data_SMM,data.frame(sources=k,val_ini=j,lambda_T=val,Q1=stat_1[1],M=stat_1[2],Q3=stat_1[3]))
  stat_data_SMM_delta=rbind(stat_data_SMM_delta,data.frame(sources=k,val_ini=j,lambda_T=val,Q1=stat_2[1],M=stat_2[2],Q3=stat_2[3]))
  }#loop λT
} #loop val ini
}#loop sources

##JOINING THE TWO DATAFRAMES
final_data_set_SMM=rbind(stat_data_SMM_s1,stat_data_SMM)
final_data_set_SMM_delta=rbind(stat_data_SMM_delta_s1,stat_data_SMM_delta)

#graph
#diet 4
col_val=c(1,5,10)
#diet 3
# col_val=c(2,3,8)
graph_source_effect(data_SMM=final_data_set_SMM,data_SMM_delta=final_data_set_SMM_delta,nb_values=10,col_val=col_val)



##-----------SOURCES POLYGONS------------------
deg_col=viridis(n=10, alpha = 1, begin = 0, end = 1, direction = -1, option = "D")
conso_bis=initial_values[c(13,7,1),]
colors=deg_col[c(1,4,6)]
# #4 sources
# #3 sources
plot(NA, NA, ylim = c(0,10), xlim = c(0,10),
     xlab = "", ylab = "" ,las=1,main="")
mtext(expression(paste(" δ"^{15},"N (‰)")),side=2,line=4,padj=1.5,at=5,cex=1.6,font=2)
mtext(expression(paste(" δ"^{13},"C (‰)")),side=1,line=1,padj=1.5,at=5,cex=1.6,font=2)
points(c(0,5,10,0),c(0,10,5,0),type="l",col="red")
points(c(0,5,10,0),c(0,10,5,0),type="p",pch=16,col="red")
points(conso_bis$I1,conso_bis$I2,pch=16,col=colors,cex=1.3)

# text(x=c(0.5,5.5,9.5),y=c(0,10,5),labels = c("2","3","1"),cex=1.4,font=2,col="red")

# 


deg_col=viridis(n=10, alpha = 1, begin = 0, end = 1, direction = -1, option = "D")
conso_bis=initial_values[c(12,14,20),]
colors=deg_col[c(1,5,10)]
# #4 sources

plot(NA, NA, ylim = c(0,10), xlim = c(0,10),
     xlab = "", ylab = "" ,las=1,main="")
mtext(expression(paste(" δ"^{15},"N (‰)")),side=2,line=4,padj=1.5,at=5,cex=1.6,font=2)
mtext(expression(paste(" δ"^{13},"C (‰)")),side=1,line=1,padj=1.5,at=5,cex=1.6,font=2)
points(c(0,2,5,10,0),c(0,6,10,5,0),type="l",col="red")
points(c(0,2,5,10,0),c(0,6,10,5,0),type="p",pch=16,col="red")
points(conso_bis$I1,conso_bis$I2,pch=16,col=colors,cex=1.3)
# text(x=c(0.5,2.5,5.5,9.5),y=c(0,6,10,5),labels = c("2","4","3","1"),cex=1.4,font=2,col="red")

