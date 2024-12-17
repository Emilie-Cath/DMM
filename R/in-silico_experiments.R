####VERSION GENERALISEE AVEC MES CODES###

rm(list=ls())# clear the current environment
graphics.off()# clear the current plots
setwd("C:/Users/Etudiant/Documents/thèse/1ere année/1er papier/finals/1_st paper")
source("Mixing_models_functions.R")
library(deSolve)
library(philentropy)
library("viridisLite")

#-------------FUNCTIONS--------------------


#data the results of the loop (data frame with the bias, lambda*t,soource var, initial value number and the model results)
#nb_values the number of switchs values tested
#ini_val the data frame with the initials values
#ref a vector with the ref values for each isotope ex:c(0,1)

process_switch_effect <- function(data,nb_values,ini_val,ref) {
  data_cste=subset(data,data$sources==1) #isolating the data where the source is constant
  
  
  ##euclidian distance
  ini_ref=cbind(val_nb=seq(1,nb_values),ini_val,ref_I1=rep(ref[1],nrow(ini_val)),ref_I2=rep(ref[2],nrow(ini_val)))
  dist_eucli=sqrt((ini_ref$ref_I1-ini_ref$I1)^2+(ini_ref$ref_I2-ini_ref$I2)^2)
  val_dist=cbind(val_ini=seq(1,nb_values),ini_val,dist_eucli)
  ord=order(val_dist$dist_eucli)
  val_dist_ordered=val_dist[ord,]
  
  
  
  
  return(val_dist_ordered)
}

#switch_effect the output of the process_switch_effect function
#scale_eucl_dist is the extremes values of the euclidian distance vector c(1)
#sources_I1 the sources values for the Iso 1 c(10,0,5)
#sources_I2 the sources values for the Iso 2 c(5,0,10)
#poly either T or F, to know if the polygon should be represented
graph_switch_effect <- function(switch_effect,scale_eucl_dist,nb_values,ref,sources_I1,sources_I2,poly,low_dist,med_low_dist,
                                med_high_dist,high_dist) {
  
  #colors
  deg_col=viridis(n=nb_values, alpha = 1, begin = 0, end = 1, direction = -1, option = "D")
  legend_image <- as.raster(matrix(rev(deg_col), ncol=1))
  nb_cat=nb_values/4
  low_col=deg_col[1:nb_cat]
  med_low_col=deg_col[(nb_cat+1):(2*nb_cat)]
  med_high_col=deg_col[(2*nb_cat+1):(3*nb_cat)]
  high_col=deg_col[(3*nb_cat+1):(4*nb_cat)]
  #data
  
  
  windows()
  par(mfrow=c(2,3))
  #legend
  plot(c(0,2),c(0,1),type = 'n', axes = F,xlab = '', ylab = '', main = 'Euclidian distance')
  # Tracer les graduations de la légende, ici le paramétrage en est strictement manuel 5 paliers entre 0 et 1
  text(x=1.5, y = seq(0,1,l=5), labels = seq(scale_eucl_dist[1],scale_eucl_dist[2],l=5))
  # Ajouter l'image de gradient
  rasterImage(legend_image, 0, 0, 1,1)
  #polygon
  plot(NA, NA, ylim = c(0,10), xlim = c(0,10),
       xlab = "", ylab = "" ,las=1)
  mtext(expression(paste(" δ"^{15},"N (‰)")),side=2,line=4,padj=1.5,at=5,cex=1,font=2)
  mtext(expression(paste(" δ"^{13},"C (‰)")),side=1,line=1,padj=1.5,at=5,cex=1,font=2)
  if (poly==T) {points(c(sources_I1,sources_I1[1]),c(sources_I2,sources_I2[1]),type="l",col="red")
    points(c(sources_I1,sources_I1[length(sources_I1)]),c(sources_I2,sources_I2[length(sources_I2)]),type="p",pch=16,col="red")
    text(x=sources_I1+0.5,y=sources_I2,labels = seq(1,length(sources_I1)),cex=1.4,font=2)}
  points(switch_effect$I1,switch_effect$I2,col=deg_col,type="p",pch=16)
  points(ref[1],ref[2],pch=17,col="black",cex=1.4)
  
  #1
  plot(NA, NA, ylim = c(0,1), xlim = c(0,6),
       xlab = c("λT"), ylab =c("Bias") , las=1, cex.lab=1.4,font.lab=2, cex.axis=1.4,main="Low distance ")
  points(low_dist$lambda_T,low_dist$B,type="p",pch=16,cex=1.4,col=low_col)
  #2
  plot(NA, NA, ylim = c(0,1), xlim = c(0,6),
       xlab = c("λT"), ylab =c("Bias") , las=1, cex.lab=1.4,font.lab=2, cex.axis=1.4,main="med Low distance ")
  points(med_low_dist$lambda_T,med_low_dist$B,type="p",pch=16,cex=1.4,col=med_low_col)
  #3
  plot(NA, NA, ylim = c(0,1), xlim = c(0,6),
       xlab = c("λT"), ylab =c("Bias") , las=1, cex.lab=1.4,font.lab=2, cex.axis=1.4,main="med high distance ")
  points(med_high_dist$lambda_T,med_high_dist$B,type="p",pch=16,cex=1.4,col=med_high_col)
  #4
  plot(NA, NA, ylim = c(0,1), xlim = c(0,6),
       xlab = c("λT"), ylab =c("Bias") , las=1, cex.lab=1.4,font.lab=2, cex.axis=1.4,main="High distance ")
  points(high_dist$lambda_T,high_dist$B,type="p",pch=16,cex=1.4,col=high_col)
}


#data_cst a subset of the data with only the constant sources
#ordered_data is the output of process_switch_effect function
#scale_eucl_dist is the extremes values of the euclidian distance vector c(1)
#sources_I1 the sources values for the Iso 1 c(10,0,5)
#sources_I2 the sources values for the Iso 2 c(5,0,10)

graph_switch_effect_2<- function(data_cst,ordered_data,nb_values
                                 ,scale_eucl_dist,ref,sources_I1,sources_I2) {
  data_merge=merge(data_cst,ordered_data,by=c("val_ini"))
  ord=order(data_merge$dist_eucli)
  merge_ordered=data_merge[ord,] #ordering data according to val_ini and its eucli dist
  
  order_value=ordered_data$val_ini
  
  deg_col=viridis(n=nb_values, alpha = 1, begin = 0,
                  end = 1, direction = -1, option = "D") #setting the palette
  legend_image <- as.raster(matrix(rev(deg_col), ncol=1))
  # #graph1
  windows()
  plot(NA, NA, ylim = c(0,1), xlim = c(0,6),
       xlab = c("λT"), ylab =c("Bias") , las=1, cex.lab=1.4,font.lab=2, cex.axis=1.4,main="")
  for (i in 1:nb_values) {subs=subset(merge_ordered,merge_ordered$val_ini==order_value[i])
  ord_2=order(subs$lambda_T)
  merge_ordered_2=subs[ord_2,]
  points(merge_ordered_2$lambda_T,merge_ordered_2$B,type="p",pch=16,cex=2,lwd=5,col=deg_col[i])}
  
  ##graph2
  windows()
  par(mfrow=c(1,2))
  #legend
  plot(c(0,2),c(0,1),type = 'n', axes = F,xlab = '', ylab = '', main = 'Euclidian distance')
  # Tracer les graduations de la légende, ici le paramétrage en est strictement manuel 5 paliers entre 0 et 1
  text(x=1.5, y = seq(0,1,l=5), labels = seq(scale_eucl_dist[1],scale_eucl_dist[2],l=5))
  # Ajouter l'image de gradient
  rasterImage(legend_image, 0, 0, 1,1)
  #polygon
  plot(NA, NA, ylim = c(0,10), xlim = c(0,10),
       xlab = "", ylab = "" ,las=1)
  mtext(expression(paste(" δ"^{15},"N (‰)")),side=2,line=4,padj=1.5,at=5,cex=1,font=2)
  mtext(expression(paste(" δ"^{13},"C (‰)")),side=1,line=1,padj=1.5,at=5,cex=1,font=2)
  points(c(sources_I1,sources_I1[1]),c(sources_I2,sources_I2[1]),type="l",col="red")
  points(c(sources_I1,sources_I1[length(sources_I1)]),c(sources_I2,sources_I2[length(sources_I2)]),type="p",pch=16,col="red")
  text(x=sources_I1+0.5,y=sources_I2,labels = seq(1,length(sources_I1)),cex=1.4,font=2)
  points(ordered_data$I1,ordered_data$I2,col=deg_col,type="p",pch=16)
  points(ref[1],ref[2],pch=17,col="black",cex=1.4)
  
}




#val a vector with the 4 studied initial values ex c(26,1,8,31)
#row_val a vector with the row position of each chosen value ex c(1,12,19,32)
process_source_effect <- function(data,val) {
  extr_values_ordered=NULL
  for (i in 1: length(val)) {
    extr_data=subset(data,data$val_ini==val[i])
    val_cst=subset(extr_data,extr_data$sources==1)
    val_2=subset(extr_data,extr_data$sources==2)
    val_3=subset(extr_data,extr_data$sources==3)
    val_4=subset(extr_data,extr_data$sources==4)
    val_5=subset(extr_data,extr_data$sources==5)
    
    comp=data.frame("1"=val_cst$B,"2"=val_2$B,"3"=val_3$B,"4"=val_4$B,"5"=val_5$B)
    max=apply(comp,1,max)
    min=apply(comp,1,min)
    
    extr_values=data.frame(lambda_T=val_cst$lambda_T,min,max)
    in_order=order(extr_values$lambda_T)
    extr_values_ordered=c(extr_values_ordered,list(extr_values[in_order,]))
  }
  return(extr_values_ordered)
}


#process is the output of process_source_effect
graph_source_effect <- function(data,val,nb_values,val_row,scale_eucl_dist) {
  #colors
  deg_col=viridis(n=nb_values, alpha = 1, begin = 0, end = 1, direction = -1, option = "D")
  legend_image <- as.raster(matrix(rev(deg_col), ncol=1))
  col_val=NULL
  for (j in 1:length(val_row)) {col_val=c(col_val,deg_col[val_row[j]])}
  #process
  process=process_source_effect(data,val)
  low_dist=process[[1]]
  med_low_dist=process[[2]]
  med_high_dist=process[[3]]
  high_dist=process[[4]]
  
  #cstes source values
  data_cst=subset(data,data$sources==1)
  ref_val=NULL
  for (i in 1:length(val)){ref_val=c(ref_val,list(subset(data_cst,data_cst$val_ini==val[i])))}
  
  
  windows()
  par(mfrow=c(2,3))
  
  #legend
  plot(c(0,2),c(0,1),type = 'n', axes = F,xlab = '', ylab = '', main = 'Euclidian distance')
  # Tracer les graduations de la légende, ici le paramétrage en est strictement manuel 5 paliers entre 0 et 1
  text(x=1.5, y = seq(0,1,l=5), labels = seq(scale_eucl_dist[1],scale_eucl_dist[2],l=5))
  # Ajouter l'image de gradient
  rasterImage(legend_image, 0, 0, 1,1)
  
  
  #graph_1
  plot(NA, NA, ylim = c(0,1), xlim = c(0,6),
       xlab = c("λT"), ylab =c("Bias") , las=1, cex.lab=1.4,font.lab=2, cex.axis=1.4,main="Low distance source variation effect")
  polygon(c(low_dist$lambda_T,rev(low_dist$lambda_T)),c(low_dist$min,rev(low_dist$max)),col="grey",border=NA)
  points(ref_val[[1]]$lambda_T,ref_val[[1]]$B,type="p",col=col_val[[1]],pch=16,cex=1.4)
  
  #graph_2
  plot(NA, NA, ylim = c(0,1), xlim = c(0,6),
       xlab = c("λT"), ylab =c("Bias") , las=1, cex.lab=1.4,font.lab=2, cex.axis=1.4,main="Med low")
  polygon(c(med_low_dist$lambda_T,rev(med_low_dist$lambda_T)),c(med_low_dist$min,rev(med_low_dist$max)),col="grey",border=NA)
  points(ref_val[[2]]$lambda_T,ref_val[[2]]$B,type="p",col=col_val[[2]],pch=16,cex=1.4)
  
  
  #graph_3
  plot(NA, NA, ylim = c(0,1), xlim = c(0,6),
       xlab = c("λT"), ylab =c("Bias") , las=1, cex.lab=1.4,font.lab=2, cex.axis=1.4,main="Med high")
  polygon(c(med_high_dist$lambda_T,rev(med_high_dist$lambda_T)),c(med_high_dist$min,rev(med_high_dist$max)),col="grey",border=NA)
  points(ref_val[[3]]$lambda_T,ref_val[[3]]$B,type="p",col=col_val[[3]],pch=16,cex=1.4)
  
  
  #graph_4
  plot(NA, NA, ylim = c(0,1), xlim = c(0,6),
       xlab = c("λT"), ylab =c("Bias") , las=1, cex.lab=1.4,font.lab=2, cex.axis=1.4,main="High")
  polygon(c(high_dist$lambda_T,rev(high_dist$lambda_T)),c(high_dist$min,rev(high_dist$max)),col="grey",border=NA)
  points(ref_val[[4]]$lambda_T,ref_val[[4]]$B,type="p",col=col_val[[4]],pch=16,cex=1.4)
  
}




#-------------Loop8: the graph maker with all the effects  (2I,3S)--------------------

#CHANGE HERE IF THE SOURCE NUMBER OR ISOTOPE NUMBER CHANGE
TEF=list(data.frame("s1"=c(0),"s2"=c(0),"s3"=c(0)),data.frame("s1"=c(0),"s2"=c(0),"s3"=c(0))) #the TEF=0 for all sources and isotope
conc=list(data.frame("s1"=c(1),"s2"=c(1),"s3"=c(1)),data.frame("s1"=c(1),"s2"=c(1),"s3"=c(1))) #the conc=1 for all sources and isotope
iter=0.05


#CONSUMER INITAL VALUE
nb_ini_val=20

resum_data_SMM=data.frame(B=c(),l_T=c(),sources=c(),val_ini=c())
resum_data_SMM_delta=data.frame(B=c(),l_T=c(),sources=c(),val_ini=c())

conso=read.csv2("val_ini_2.csv",row.names = NULL)


#DIET
#the consumer eats 100% of source 3
# combi=data.frame(s1=c(0),s2=c(0),s3=c(1))
#theconsumer eats a mix and the solution value is (3.5,4)
combi=data.frame(s1=c(0.15),s2=c(0.45),s3=c(0.4))


#SETTING λT values
p=seq(2,300,length=4) #T
l=seq(0.002,0.02,length=5) #lambda

for (s in 1:length(l)){ #for every lambda
  # s=1
  lambda=l[s]
  print(c("l",s))
  
  for (f in 1:length(p)){ #for everyT 
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
    
    
    
    
    ##THE LOOP STARTS HERE
    for(m in 1:nrow(conso)){ #for each consumer
      
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
      for (n in 1:length(sources_list_var)) { #for each source variation
        
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

        SMM=SMM_func(t=period,sources_approx=sources_approxfun,list_TEF=TEF,details_conc=conc
                     ,combi=sol,data_feat=data_feat,details_consu=val_conso,X=1)
  
        SMM_delta=SMM_delta_func(t=period,sources_approx=sources_approxfun,list_TEF=TEF,details_conc=conc,combi=sol
                                 ,data_feat=data_feat,details_consu=val_conso,X=1,lambda_list=liste_lambda)
        
        
        #estimating bias and archiving results
        resum_data_SMM=rbind(resum_data_SMM,data.frame(B=abs(combi[[1]]+combi[[2]]-SMM[[2]]-SMM[[3]]),lambda_T=c(lambda*period),sources=c(n),val_ini=c(m),SMM[-1]))
        resum_data_SMM_delta=rbind(resum_data_SMM_delta,data.frame(B=abs(combi[[1]]+combi[[2]]-SMM_delta[[2]]-SMM_delta[[3]]),lambda_T=c(lambda*period),sources=c(n),val_ini=c(m),SMM_delta[-1]))
        
      }#loop on sources

    }#loop on ini value
    
  }#loop on period
}#loop on lambda


 # write.csv2(resum_data_SMM,file="data_SMM_val_ini_2_diete2.csv")
 # write.csv2(resum_data_SMM_delta,file="data_SMM_delta_val_ini_2_diete2.csv")
#  # write.csv2(resum_data_SMM,file="data_SMM_val_X_sc_variables.csv")
 # write.csv2(resum_data_SMM_delta,file="data_SMM_delta_val_X_sc_variables.csv")






#-------------Loop9: the graph maker with all the effects  (2I,4S)--------------------
##SOURCE POLYGON

# #3 sources
# plot(NA, NA, ylim = c(-5,10), xlim = c(-5,10),
#      xlab = "", ylab = "" ,las=1,main="a")
# mtext(expression(paste(" δ"^{15},"N (‰)")),side=2,line=4,padj=1.5,at=5,cex=1.6,font=2)
# mtext(expression(paste(" δ"^{13},"C (‰)")),side=1,line=1,padj=1.5,at=5,cex=1.6,font=2)
# # points(val_ini$I1,val_ini$I2,pch=16)
# points(c(0,5,10,0),c(0,10,5,0),type="l",col="red")
# points(c(0,5,10,0),c(0,10,5,0),type="p",pch=16,col="red")
# text(x=c(0.5,5.5,9.5),y=c(0,10,5),labels = c("2","3","1"),cex=1.4,font=2,col="red")

# 
# #4 sources
# plot(NA, NA, ylim = c(0,10), xlim = c(0,10),
#      xlab = "", ylab = "" ,las=1,main="b")
# mtext(expression(paste(" δ"^{15},"N (‰)")),side=2,line=4,padj=1.5,at=5,cex=1.6,font=2)
# mtext(expression(paste(" δ"^{13},"C (‰)")),side=1,line=1,padj=1.5,at=5,cex=1.6,font=2)
# points(val_ini$I1,val_ini$I2,pch=16)
# points(c(0,2,5,10,0),c(0,6,10,5,0),type="l",col="red")
# points(c(0,2,5,10,0),c(0,6,10,5,0),type="p",pch=16,col="red")
# text(x=c(0.5,2.5,5.5,9.5),y=c(0,6,10,5),labels = c("2","4","3","1"),cex=1.4,font=2,col="red")





#CHANGE HERE IF THE SOURCE NUMBER OR ISOTOPE NUMBER CHANGE
TEF=list(data.frame("s1"=c(0),"s2"=c(0),"s3"=c(0),"s4"=c(0)),data.frame("s1"=c(0),"s2"=c(0),"s3"=c(0),"s4"=c(0))) #the TEF=0 for all sources and isotope
conc=list(data.frame("s1"=c(1),"s2"=c(1),"s3"=c(1),"s4"=c(1)),data.frame("s1"=c(1),"s2"=c(1),"s3"=c(1),"s4"=c(1))) #the conc=1 for all sources and isotope
iter=0.05


#CONSUMER INITAL VALUE
nb_ini_val=20

conso=read.csv2("val_ini_2.csv",row.names = NULL)
resum_data_SMM=data.frame(B=c(),l_T=c(),sources=c(),val_ini=c())
resum_data_SMM_delta=data.frame(B=c(),l_T=c(),sources=c(),val_ini=c())



#DIET
#the consumer eats 100% of source 3
# combi=data.frame(s1=c(0),s2=c(0),s3=c(1),s4=c(0))
#theconsumer eats a mix and the solution value is (3.5,4)
combi=data.frame(s1=c(0.2),s2=c(0.2),s3=c(0.4),s4=c(0.2))



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
    
    
    
    
    ##THE LOOP STARTS HERE
    for(m in 1:nrow(conso)){ #for each consumer
      # m=1
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
      for (n in 1:length(sources_list_var)) { #for each source variation
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
        
        SMM=SMM_func(t=period,sources_approx=sources_approxfun,list_TEF=TEF,details_conc=conc
                     ,combi=sol,data_feat=data_feat,details_consu=val_conso,X=1)
        SMM_delta=SMM_delta_func(t=period,sources_approx=sources_approxfun,list_TEF=TEF,details_conc=conc,combi=sol
                                 ,data_feat=data_feat,details_consu=val_conso,X=1,lambda_list=liste_lambda)
    
        
        #ESTIMATING SMM AND SMM-DELTA WITH THE EXACT SOL MIXING MODEL at the last t
     
        resum_data_SMM=rbind(resum_data_SMM,data.frame(B=abs(combi[[1]]+combi[[2]]+combi[[4]]-SMM[[2]]-SMM[[3]]-SMM[[5]]),lambda_T=c(lambda*period),sources=c(n),val_ini=c(m),SMM[-1]))
        resum_data_SMM_delta=rbind(resum_data_SMM_delta,data.frame(B=abs(combi[[1]]+combi[[2]]+combi[[4]]-SMM_delta[[2]]-SMM_delta[[3]]-SMM_delta[[5]]),lambda_T=c(lambda*period),sources=c(n),val_ini=c(m),SMM_delta[-1]))
        
      }#loop on sources
      # lines(rep(lambda*period,length(sources_list_var)),bias_SMM,pch=16,type="p",cex=1.5,col=col_conso[m])
      # lines(rep(lambda*period,length(sources_list_var)),bias_SMM_delta,pch=16,type="p",col="red",cex=1.5)
    }#loop on ini value
    
  }#loop on period
}#loop on lambda

# 
#  write.csv2(resum_data_SMM,file="data_SMM_val_ini_2_diete3.csv")
#  write.csv2(resum_data_SMM_delta,file="data_SMM_delta_val_ini_2_diete3.csv")
#  write.csv2(resum_data_SMM_delta,file="data_SMM_delta_32_val_MIX_4S.csv")
# 
 
#-------------GRAPHS 0.2,0.2,0.4,0.2--------------------

#sources cstes
data_3=read.csv2("data_SMM_val_ini_2_diete3.csv")
val_ini=read.csv2("val_ini_2.csv",row.names = NULL)

ordered_data=process_switch_effect(data=data_3,nb_values=20,ini_val=val_ini,ref=c(4.4,6.2))
data_3_cst=subset(data_3,data_3$sources==1)

low_dist=subset(data_3_cst,val_ini==9|val_ini==14|val_ini==6|val_ini==15|val_ini==7)
med_low_dist=subset(data_3_cst,val_ini==13|val_ini==8|val_ini==5|val_ini==10|val_ini==4)
med_high_dist=subset(data_3_cst,val_ini==11|val_ini==12|val_ini==3|val_ini==17|val_ini==16)
high_dist=subset(data_3_cst,val_ini==2|val_ini==1|val_ini==18|val_ini==19|val_ini==20)

graph=graph_switch_effect(switch_effect=ordered_data,scale_eucl_dist=c(1,8),nb_values=20
                          ,ref=c(4.4,6.2),sources_I1=c(10,0,2,5),sources_I2=c(5,0,6,10),poly=T,
                          low_dist,med_low_dist,
                          med_high_dist,high_dist)

graph_switch_effect_2(data_cst=data_3_cst,ordered_data=ordered_data,nb_values=20
                      ,scale_eucl_dist=c(1,8),ref=c(4.4,6.2)
                      ,sources_I1=c(10,0,2,5),sources_I2=c(5,0,6,10))
#process SMM
data_SMM=read.csv2("data_SMM_val_ini_2_diete3.csv")
graph_source_effect(data=data_SMM,val=c(9,5,16,20),nb_values=20,val_row=c(1,8,15,20),scale_eucl_dist=c(0.9,7.5))

#process_SMM_delta
data_SMM_delta=read.csv2("data_SMM_delta_val_ini_2_diete3.csv")
graph_source_effect(data=data_SMM_delta,val=c(9,5,16,20),nb_values=20,val_row=c(1,8,15,20),scale_eucl_dist=c(0.9,7.5))


# WORK IN PROGRESS GRAPH
data_3_cst=subset(data_3,data_3$sources==1)
deg_col=viridis(n=20, alpha = 1, begin = 0, end = 1, direction = -1, option = "D")
legend_image <- as.raster(matrix(rev(deg_col), ncol=1))
plot(NA, NA, ylim = c(0,1), xlim = c(0,6),
     xlab = c("λT"), ylab =c("Bias") , las=1, cex.lab=1.4,font.lab=2, cex.axis=1.4,main="Low distance ")
points(data_3_cst$lambda_T,data_3_cst$B,type="p",pch=16,cex=1.4)

essai_merge=merge(data_3_cst,ordered_data,by=c("val_ini"))
ord=order(essai_merge$dist_eucli)
merge_ordered=essai_merge[ord,]
nb_values=20
order_value=ordered_data$val_ini
plot(NA, NA, ylim = c(0,1), xlim = c(0,6),
     xlab = c("λT"), ylab =c("Bias") , las=1, cex.lab=1.4,font.lab=2, cex.axis=1.4,main="Low distance ")
for (i in 1:nb_values) {subs=subset(merge_ordered,merge_ordered$val_ini==order_value[i])
points(subs$lambda_T,subs$B,type="p",pch=16,cex=1.4,col=deg_col[i])}

#-------------GRAPHS 1,0,0 v2--------------------

data_1=read.csv2("data_SMM_val_ini_2_diete1.csv")
val_ini=read.csv2("val_ini_2.csv",row.names = NULL) #val ini 2
#ordering the initial value according to their euclidian distance
ordered_data=process_switch_effect(data=data_1,nb_values=20,ini_val=val_ini,ref=c(5,10))
#constant sources part
data_1_cst=subset(data_1,data_1$sources==1)
#splitting the initial values between four groups according to their eucl dist
low_dist=subset(data_1_cst,val_ini==12|val_ini==11|val_ini==8|val_ini==7|val_ini==9)
med_low_dist=subset(data_1_cst,val_ini==14|val_ini==4|val_ini==6|val_ini==5|val_ini==15)
med_high_dist=subset(data_1_cst,val_ini==13|val_ini==3|val_ini==2|val_ini==1|val_ini==10)
high_dist=subset(data_1_cst,val_ini==17|val_ini==16|val_ini==18|val_ini==19|val_ini==20)
#graphs for sctes sources
graph=graph_switch_effect(switch_effect=ordered_data,scale_eucl_dist=c(2,11),nb_values=20
                          ,ref=c(5,10),sources_I1=c(10,0,5),sources_I2=c(5,0,10),poly=T,
                          low_dist,med_low_dist,
                          med_high_dist,high_dist)


graph_switch_effect_2(data_cst=data_1_cst,ordered_data=ordered_data,nb_values=20
                      ,scale_eucl_dist=c(1,11),ref=c(5,10)
                      ,sources_I1=c(10,0,5),sources_I2=c(5,0,10))


#process SMM
data_SMM=read.csv2("data_SMM_val_ini_2_diete1.csv")
graph_source_effect(data=data_SMM,val=c(12,14,3,20),nb_values=20,val_row=c(1,6,12,20),scale_eucl_dist=c(0,11))

#process SMM_delta
data_SMM_delta=read.csv2("data_SMM_delta_val_ini_2_diete1.csv")
graph_source_effect(data=data_SMM_delta,val=c(12,14,3,20),nb_values=20,val_row=c(1,6,12,20),scale_eucl_dist=c(0,11))

#-------------work space--------------------
data_1=read.csv2("data_SMM_val_ini_2_diete1.csv")
val_ini=read.csv2("val_ini_2.csv",row.names = NULL) #val ini 2
#ordering the initial value according to their euclidian distance
ordered_data=process_switch_effect(data=data_1,nb_values=20,ini_val=val_ini,ref=c(5,10))
#on garde c(12,14,4,20)
#process SMM
data_SMM=read.csv2("data_SMM_val_ini_2_diete1.csv")
# graph_source_effect(data=data_SMM,val=c(12,14,3,20),nb_values=20,val_row=c(1,6,12,10),scale_eucl_dist=c(0,11))

s_ef=process_source_effect(data=data_SMM,val=c(12,14,4,20))

graph_source_effect <- function(data,val,nb_values,val_row,scale_eucl_dist) {
  #colors
  deg_col=viridis(n=nb_values, alpha = 1, begin = 0, end = 1, direction = -1, option = "D")
  legend_image <- as.raster(matrix(rev(deg_col), ncol=1))
  col_val=NULL
  for (j in 1:length(val_row)) {col_val=c(col_val,deg_col[val_row[j]])}
  #process
  process=process_source_effect(data,val)
  low_dist=process[[1]]
  med_low_dist=process[[2]]
  med_high_dist=process[[3]]
  high_dist=process[[4]]
  
  #cstes source values
  data_cst=subset(data,data$sources==1)
  ref_val=NULL
  for (i in 1:length(val)){ref_val=c(ref_val,list(subset(data_cst,data_cst$val_ini==val[i])))}
  
  
  windows()
  par(mfrow=c(2,3))
  
  #legend
  plot(c(0,2),c(0,1),type = 'n', axes = F,xlab = '', ylab = '', main = 'Euclidian distance')
  # Tracer les graduations de la légende, ici le paramétrage en est strictement manuel 5 paliers entre 0 et 1
  text(x=1.5, y = seq(0,1,l=5), labels = seq(scale_eucl_dist[1],scale_eucl_dist[2],l=5))
  # Ajouter l'image de gradient
  rasterImage(legend_image, 0, 0, 1,1)
  
  
  #graph_1
  plot(NA, NA, ylim = c(0,1), xlim = c(0,6),
       xlab = c("λT"), ylab =c("Bias") , las=1, cex.lab=1.4,font.lab=2, cex.axis=1.4,main="Low distance source variation effect")
  polygon(c(low_dist$lambda_T,rev(low_dist$lambda_T)),c(low_dist$min,rev(low_dist$max)),col="grey",border=NA)
  points(ref_val[[1]]$lambda_T,ref_val[[1]]$B,type="p",col=col_val[[1]],pch=16,cex=1.4)
  
  #graph_2
  plot(NA, NA, ylim = c(0,1), xlim = c(0,6),
       xlab = c("λT"), ylab =c("Bias") , las=1, cex.lab=1.4,font.lab=2, cex.axis=1.4,main="Med low")
  polygon(c(med_low_dist$lambda_T,rev(med_low_dist$lambda_T)),c(med_low_dist$min,rev(med_low_dist$max)),col="grey",border=NA)
  points(ref_val[[2]]$lambda_T,ref_val[[2]]$B,type="p",col=col_val[[2]],pch=16,cex=1.4)
  
  
  #graph_3
  plot(NA, NA, ylim = c(0,1), xlim = c(0,6),
       xlab = c("λT"), ylab =c("Bias") , las=1, cex.lab=1.4,font.lab=2, cex.axis=1.4,main="Med high")
  polygon(c(med_high_dist$lambda_T,rev(med_high_dist$lambda_T)),c(med_high_dist$min,rev(med_high_dist$max)),col="grey",border=NA)
  points(ref_val[[3]]$lambda_T,ref_val[[3]]$B,type="p",col=col_val[[3]],pch=16,cex=1.4)
  
  
  #graph_4
  plot(NA, NA, ylim = c(0,1), xlim = c(0,6),
       xlab = c("λT"), ylab =c("Bias") , las=1, cex.lab=1.4,font.lab=2, cex.axis=1.4,main="High")
  polygon(c(high_dist$lambda_T,rev(high_dist$lambda_T)),c(high_dist$min,rev(high_dist$max)),col="grey",border=NA)
  points(ref_val[[4]]$lambda_T,ref_val[[4]]$B,type="p",col=col_val[[4]],pch=16,cex=1.4)
  
}
#-------------GRAPHS 0.15,0.45,0.4 V2--------------------

data_2=read.csv2("data_SMM_val_ini_2_diete2.csv")
val_ini=read.csv2("val_ini_2.csv",row.names = NULL)
#ordering the initial value according to their euclidian distance
ordered_data=process_switch_effect(data=data_2,nb_values=20,ini_val=val_ini,ref=c(3.5,4))
#constant sources part
data_2_cst=subset(data_2,data_2$sources==1)
#splitting the initial values between four groups according to their eucl dist

low_dist=subset(data_2_cst,val_ini==13|val_ini==15|val_ini==10|val_ini==17|val_ini==14)
med_low_dist=subset(data_2_cst,val_ini==16|val_ini==6|val_ini==9|val_ini==18|val_ini==5)
med_high_dist=subset(data_2_cst,val_ini==3|val_ini==7|val_ini==19|val_ini==4|val_ini==20)
high_dist=subset(data_2_cst,val_ini==8|val_ini==2|val_ini==11|val_ini==12|val_ini==1)
#graphs for sctes sources
graph=graph_switch_effect(switch_effect=ordered_data,scale_eucl_dist=c(0.8,6.4),nb_values=20
                          ,ref=c(3.5,4),sources_I1=c(10,0,5),sources_I2=c(5,0,10),poly=T,
                          low_dist,med_low_dist,
                          med_high_dist,high_dist)

graph_switch_effect_2(data_cst=data_2_cst,ordered_data=ordered_data,nb_values=20
                      ,scale_eucl_dist=c(0.8,6.4),ref=c(3.5,4)
                      ,sources_I1=c(10,0,5),sources_I2=c(5,0,10))
#process SMM
data_SMM=read.csv2("data_SMM_val_ini_2_diete2.csv")
graph_source_effect(data=data_SMM,val=c(15,6,19,1),nb_values=20,val_row=c(2,7,12,20),scale_eucl_dist=c(0.8,6.4))

#process SMM_delta
data_SMM_delta=read.csv2("data_SMM_delta_val_ini_2_diete2.csv")
graph_source_effect(data=data_SMM_delta,val=c(15,6,19,1),nb_values=20,val_row=c(2,7,12,20),scale_eucl_dist=c(0.8,6.4))

#-------------ISOTOPIC SPACE EFFECT--------------------

#CHANGE HERE IF THE SOURCE NUMBER OR ISOTOPE NUMBER CHANGE
TEF=list(data.frame("s1"=c(0),"s2"=c(0),"s3"=c(0)),data.frame("s1"=c(0),"s2"=c(0),"s3"=c(0))) #the TEF=0 for all sources and isotope
conc=list(data.frame("s1"=c(1),"s2"=c(1),"s3"=c(1)),data.frame("s1"=c(1),"s2"=c(1),"s3"=c(1))) #the conc=1 for all sources and isotope
iter=0.1

#ISOTOPIC SPACE
nb_iso_spaces=6
s3=data.frame(I1=rep(5,nb_iso_spaces),I2=rep(10,nb_iso_spaces))
s1=data.frame(I1=c(10,6,7,12,6,14),I2=c(5,9,6,2,9,-2))
s2=data.frame(I1=c(0,5,2,-5,4,-10),I2=c(0,8,2,-5,5,-10))

plot(NULL,xlim=c(-10,15),ylim=c(-10,10),xlab="",ylab="",las=1)
mtext(expression(paste(" δ"^{15},"N (‰)")),side=2,line=4,padj=1.5,at=1,cex=1.2,font=2)
mtext(expression(paste(" δ"^{13},"C (‰)")),side=1,line=1,padj=1.5,at=5,cex=1.2,font=2)
for (i in 1:nrow(s3)){
  points(c(s1[i,1],s2[i,1],s3[i,1],s1[i,1]),c(s1[i,2],s2[i,2],s3[i,2],s1[i,2]),type="l")
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
      
      indicator_value=NULL
      bias_SMM=NULL
      bias_SMM_delta=NULL
      bias_DMM=NULL
      dist_KL_SMM=NULL
      dist_KL_SMM_delta=NULL
      SMM_arch=data.frame(row.names = NULL)
      SMM_delta_arch=data.frame(row.names = NULL)
      DMM_arch=data.frame(row.names = NULL)
        
        
        
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
        ind= c(ind,mean(indicator(res[[k]],sources_list =sources_list[[k]] ,lambda=lambda))) #COMPUTING THE INDICATOR FOR EACH ISO
        
        
        
        new_conso_value=c(new_conso_value,list(data.frame(t=seq(0,period,1),DX=res[[k]][[2]])))
        }
        #ESTIMATING THE INDICATOR
        indicator_value=c(indicator_value,sum(ind))
        
        
        
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
                     ,combi=sol,data_feat=data_feat,details_consu=val_conso,X=1)
 
        
        #CALCULATING EUCLIDIAN DISTANCE BETWEEN DIST INI AND EQ
        
   
        
        resum_data_SMM=rbind(resum_data_SMM,data.frame(B=abs(combi[[1]]+combi[[2]]-SMM[[2]]-SMM[[3]])
                                                       ,lambda_T=c(lambda*period)
                                                       ,dist_eucli=c(euclidian_distance(c(source2[[1]],source2[[2]]),c(source3[[1]],source3[[2]])))
                                                       ,iso_space=c(m),SMM[-1]))

    }#loop on ini value
    
  }#loop on period
}#loop on lambda

deg_col=viridis(n=nb_iso_spaces, alpha = 1, begin = 0, end = 1, direction = -1, option = "D")
legend_image <- as.raster(matrix(rev(deg_col), ncol=1))
dist_possibles=resum_data_SMM$dist_eucli[1:6]
dist_eucli_sorted=sort(dist_possibles)
colors <-deg_col[rank(dist_eucli_sorted)]
# colfunc<-colorRampPalette(c("red","#196F3D"))
# colors <- (colfunc(nb_iso_spaces))

plot(NULL,xlim=c(0,6),ylim=c(0,1),xlab="λT",ylab="Bias")
for (i in 1:nrow(s3)){
   i=6
  results=subset(resum_data_SMM,resum_data_SMM$iso_space==i)
  points(results$lambda_T,results$B,type="p",col=colors[i],pch=16)
}



# Tracer cette fois-ci le graphique avec sa légende
# 
# dev.off() # ferme les anciennes fenêtes
# 
# # Fenêtre du graphique (en bas à droite sur la partie x de 20% à 100% et la partie y de 0% à 90% ==> c(0.2,1,0,0.9) -  + Une marge importante en bas et à droite ==> mar = c(4,0,0.5,2)
# 
# par(fig=c(0.2,1,0,0.9), mar = c(6,4,0.5,2), new=TRUE)
# 
# plot(NULL,xlim=c(0,6),ylim=c(0,1))
# for (i in 1:nrow(s3)){
#   # i=5
#   results=subset(resum_data_SMM,resum_data_SMM$iso_space==i)
#   points(results$lambda_T,results$B,type="p",col=colors[i],pch=16)
# }
# 
# # Graphique 2 : fenêtre de la légende
# 
# par(fig=c(0,0.2,0,0.9), mar = c(2,0,0.5,2), new=TRUE)
# 
# # Préparer une image raster : gradient de couleur de la légende
# 
# legend_image <- as.raster(matrix(rev(deg_col), ncol=1))
# # Tracer un graphique blanc à la place de la légende
# 
# plot(c(0,2),c(0,1),type = 'n', axes = F,xlab = '', ylab = '', main = 'legend title',cex.main=0.5)
# 
# # Tracer les graduations de la légende, ici le paramétrage en est strictement manuel 5 paliers entre 0 et 1
# 
# text(x=1.5, y = seq(0,1,l=5), labels = seq(0,1,l=5))
# 
# # Ajouter l'image de gradient
# 
# rasterImage(legend_image, 0, 0, 1,1)
# 
# 
# # write.csv2(resum_data_SMM,file="data_SMM_val_ini_2_diete2.csv")
# # write.csv2(resum_data_SMM_delta,file="data_SMM_delta_val_ini_2_diete2.csv")
# #  # write.csv2(resum_data_SMM,file="data_SMM_val_X_sc_variables.csv")
# # write.csv2(resum_data_SMM_delta,file="data_SMM_delta_val_X_sc_variables.csv")