rm(list=ls())# clear the current environment
graphics.off()# clear the current plots
setwd("C:/Users/Etudiant/Documents/thèse/1ere année/1er papier/finals/1_st paper")
library(philentropy)
library("viridisLite")

#-------------FUNCTIONS--------------------

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

#data_cst a subset of the data with only the constant sources
#ordered_data is the output of process_switch_effect function
#scale_eucl_dist is the extremes values of the euclidian distance vector c(1)
#sources_I1 the sources values for the Iso 1 c(10,0,5)
#sources_I2 the sources values for the Iso 2 c(5,0,10)


graph_switch_effect_3<- function(data_cst,ordered_data,ref
                                 ,sources_I1,sources_I2) {
  deg_col=viridis(n=20, alpha = 1, begin = 0,
                  end = 1, direction = -1, option = "D") #setting the palette
  legend_image <- as.raster(matrix(rev(deg_col), ncol=1))
  
  
  data_merge=merge(data_cst,ordered_data,by=c("val_ini"))
  ord=order(data_merge$dist_eucli)
  merge_ordered=data_merge[ord,] #ordering data according to val_ini and its eucli dist
  
  order_value=ordered_data$val_ini
  
  # #graph1
  windows()
  plot(NA, NA, ylim = c(0,1), xlim = c(0,6),
       xlab = c("λT"), ylab =c("Bias") , las=1, cex.lab=1.4,font.lab=2, cex.axis=1.4,main="")
  for (i in 1:length(order_value)) {
    subs=subset(merge_ordered,merge_ordered$val_ini==order_value[i])
    ord_2=order(subs$lambda_T)
    merge_ordered_2=subs[ord_2,]
    col_val=round(merge_ordered_2[1,ncol(merge_ordered_2)]*20/11) #which color value
    points(merge_ordered_2$lambda_T,merge_ordered_2$B,type="p",pch=16,cex=2,lwd=5,col=deg_col[col_val])
  }
  
  ##graph2
  col_points=NULL
  for (i in 1:nrow(ordered_data)){
    col_points=c(col_points,deg_col[round(ordered_data[i,ncol(ordered_data)]*20/11)])
  }
  
  
  windows()
  par(mfrow=c(1,2))
  #legend
  plot(c(0,2),c(0,1),type = 'n', axes = F,xlab = '', ylab = '', main = 'Euclidian distance')
  # Tracer les graduations de la légende, ici le paramétrage en est strictement manuel 5 paliers entre 0 et 1
  text(x=1.5, y = seq(0,1,l=5), labels = seq(1,11,l=5))
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
  points(ordered_data$I1,ordered_data$I2,col=col_points,type="p",pch=16)
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

#-------------loading data--------------------
val_ini=read.csv2("val_ini_2.csv",row.names = NULL)

#diet 1
data_SMM_d1=read.csv2("data_SMM_val_ini_2_diete1.csv")[-1]
data_SMM_delta_d1=read.csv2("data_SMM_delta_val_ini_2_diete1.csv")[-1]

#diet2
data_SMM_d2=read.csv2("data_SMM_val_ini_2_diete2.csv")[-1]
data_SMM_delta_d2=read.csv2("data_SMM_delta_val_ini_2_diete2.csv")[-1]

#diet3
data_SMM_d3=read.csv2("data_SMM_val_ini_2_diete3.csv")[-1]
data_SMM_delta_d3=read.csv2("data_SMM_delta_val_ini_2_diete3.csv")[-1]


#-------------diet 1--------------------



#ordering the initial value according to their euclidian distance
ordered_data=process_switch_effect(data=data_SMM_d1,nb_values=20,ini_val=val_ini,ref=c(5,10))
#constant sources part
data_SMM_d1_cst=subset(data_SMM_d1,data_SMM_d1$sources==1)
#splitting the initial values between four groups according to their eucl dist
low_dist=subset(data_SMM_d1_cst,val_ini==12|val_ini==11|val_ini==8|val_ini==7|val_ini==9)
med_low_dist=subset(data_SMM_d1_cst,val_ini==14|val_ini==4|val_ini==6|val_ini==5|val_ini==15)
med_high_dist=subset(data_SMM_d1_cst,val_ini==13|val_ini==3|val_ini==2|val_ini==1|val_ini==10)
high_dist=subset(data_SMM_d1_cst,val_ini==17|val_ini==16|val_ini==18|val_ini==19|val_ini==20)
#graphs for sctes sources
graph=graph_switch_effect(switch_effect=ordered_data,scale_eucl_dist=c(2,11),nb_values=20
                          ,ref=c(5,10),sources_I1=c(10,0,5),sources_I2=c(5,0,10),poly=T,
                          low_dist,med_low_dist,
                          med_high_dist,high_dist)


graph_switch_effect_2(data_cst=data_SMM_d1_cst,ordered_data=ordered_data,nb_values=20
                      ,scale_eucl_dist=c(1,11),ref=c(5,10)
                      ,sources_I1=c(10,0,5),sources_I2=c(5,0,10))



graph_switch_effect_3(data_cst=data_SMM_d1_cst,ordered_data=ordered_data,ref=c(5,10)
                      ,sources_I1=c(10,0,5),sources_I2=c(5,0,10))

#process SMM

graph_source_effect(data=data_SMM_d1,val=c(12,14,3,20),nb_values=20,val_row=c(1,6,12,20),scale_eucl_dist=c(0,11))

#process SMM_delta

graph_source_effect(data=data_SMM_delta_d1,val=c(12,14,3,20),nb_values=20,val_row=c(1,6,12,20),scale_eucl_dist=c(0,11))

#-------------diet 2--------------------

#ordering the initial value according to their euclidian distance
ordered_data=process_switch_effect(data=data_SMM_d2,nb_values=20,ini_val=val_ini,ref=c(3.5,4))
#constant sources part
data_2_cst=subset(data_SMM_d2,data_SMM_d2$sources==1)
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

graph_switch_effect_3(data_cst=data_2_cst,ordered_data=ordered_data,ref=c(3.5,4)
                      ,sources_I1=c(10,0,5),sources_I2=c(5,0,10))

#process SMM

graph_source_effect(data=data_SMM_d2,val=c(15,6,19,1),nb_values=20,val_row=c(2,7,12,20),scale_eucl_dist=c(0.8,6.4))

#process SMM_delta
graph_source_effect(data=data_SMM_delta_d2,val=c(15,6,19,1),nb_values=20,val_row=c(2,7,12,20),scale_eucl_dist=c(0.8,6.4))

#-------------diet 3--------------------

#sources cstes
ordered_data=process_switch_effect(data=data_SMM_d3,nb_values=20,ini_val=val_ini,ref=c(4.4,6.2))
data_3_cst=subset(data_SMM_d3,data_SMM_d3$sources==1)

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

graph_switch_effect_3(data_cst=data_3_cst,ordered_data=ordered_data,ref=c(4.4,6.2)
                      ,sources_I1=c(10,0,2,5),sources_I2=c(5,0,6,10))

#process SMM
graph_source_effect(data=data_SMM_d3,val=c(9,5,16,20),nb_values=20,val_row=c(1,8,15,20),scale_eucl_dist=c(0.9,7.5))

#process_SMM_delta
graph_source_effect(data=data_SMM_delta_d3,val=c(9,5,16,20),nb_values=20,val_row=c(1,8,15,20),scale_eucl_dist=c(0.9,7.5))

#-------------source effect--------------------
data_SMM=data_SMM_d1
data_SMM_delta=data_SMM_delta_d1
#diet1
val=c(12,14,3,20)
ref=c(5,10)
# # #diet 2
# val=c(15,6,19,1)
# ref=c(3.5,4)
#diet 3
# ref=c(4.4,6.2)
# val=c(9,5,16,20)


# nb_values=20
val_row=c(1,6,12,20)
ordered_data=process_switch_effect(data=data_SMM_d1,nb_values=20,ini_val=val_ini,ref=ref)

# merged_data=merge(data,ordered_data_SMM,by=c("val_ini"))

#setting colors
deg_col=viridis(n=20, alpha = 1, begin = 0,
                end = 1, direction = -1, option = "D") #setting the palette
legend_image <- as.raster(matrix(rev(deg_col), ncol=1))

col_base=NULL
for (i in 1:length(val_row)) {
  col=round(ordered_data[val_row[i],ncol(ordered_data)]*20/11)
  if (col==0) {col=1}
  col_base=c(col_base,deg_col[col])
}

#process is the output of process_source_effect


  #process
  process_SMM=process_source_effect(data_SMM,val)
  low_dist_SMM=process_SMM[[1]]
  med_low_dist_SMM=process_SMM[[2]]
  med_high_dist_SMM=process_SMM[[3]]
  high_dist_SMM=process_SMM[[4]]
  
  process_SMM_delta=process_source_effect(data_SMM_delta,val)
  low_dist_SMM_delta=process_SMM_delta[[1]]
  med_low_dist_SMM_delta=process_SMM_delta[[2]]
  med_high_dist_SMM_delta=process_SMM_delta[[3]]
  high_dist_SMM_delta=process_SMM_delta[[4]]
  
  #cstes source values
  data_cst=subset(data_SMM,data_SMM$sources==1)
  ref_val=NULL
  for (i in 1:length(val)){ref_val=c(ref_val,list(subset(data_cst,data_cst$val_ini==val[i])))}
  
  
  windows()
  par(mfrow=c(2,2))
  
  # #legend
  # plot(c(0,2),c(0,1),type = 'n', axes = F,xlab = '', ylab = '', main = 'Euclidian distance')
  # # Tracer les graduations de la légende, ici le paramétrage en est strictement manuel 5 paliers entre 0 et 1
  # text(x=1.5, y = seq(0,1,l=5), labels = seq(1,11,l=5))
  # # Ajouter l'image de gradient
  # rasterImage(legend_image, 0, 0, 1,1)
  # 
  # 
  #graph_1
  plot(NA, NA, ylim = c(0,1), xlim = c(0,6),
       xlab = c("λT"), ylab =c("Bias") , las=1, cex.lab=1.4,font.lab=2, cex.axis=1.4,main="Low distance")
  polygon(c(low_dist_SMM$lambda_T,rev(low_dist_SMM$lambda_T)),c(low_dist_SMM$min,rev(low_dist_SMM$max)),col="grey",border=NA)
  polygon(c(low_dist_SMM_delta$lambda_T,rev(low_dist_SMM_delta$lambda_T)),c(low_dist_SMM_delta$min,rev(low_dist_SMM_delta$max)),col="grey50",border=NA)
  points(ref_val[[1]]$lambda_T,ref_val[[1]]$B,type="p",col=col_base[1],pch=16,cex=1.4)
  legend("topright",c("SMM","SMMΔ"),col=c("grey","grey50"),lty=1)
  
  #graph_2
  plot(NA, NA, ylim = c(0,1), xlim = c(0,6),
       xlab = c("λT"), ylab =c("Bias") , las=1, cex.lab=1.4,font.lab=2, cex.axis=1.4,main="Med low")
  polygon(c(med_low_dist_SMM$lambda_T,rev(med_low_dist_SMM$lambda_T)),c(med_low_dist_SMM$min,rev(med_low_dist_SMM$max)),col="grey",border=NA)
  polygon(c(med_low_dist_SMM_delta$lambda_T,rev(med_low_dist_SMM_delta$lambda_T)),c(med_low_dist_SMM_delta$min,rev(med_low_dist_SMM_delta$max)),col="grey50",border=NA)
  points(ref_val[[2]]$lambda_T,ref_val[[2]]$B,type="p",col=col_base[2],pch=16,cex=1.4)
  
  
  #graph_3
  plot(NA, NA, ylim = c(0,1), xlim = c(0,6),
       xlab = c("λT"), ylab =c("Bias") , las=1, cex.lab=1.4,font.lab=2, cex.axis=1.4,main="Med high")
  polygon(c(med_high_dist_SMM$lambda_T,rev(med_high_dist_SMM$lambda_T)),c(med_high_dist_SMM$min,rev(med_high_dist_SMM$max)),col="grey",border=NA)
  polygon(c(med_high_dist_SMM_delta$lambda_T,rev(med_high_dist_SMM_delta$lambda_T)),c(med_high_dist_SMM_delta$min,rev(med_high_dist_SMM_delta$max)),col="grey50",border=NA)
  points(ref_val[[3]]$lambda_T,ref_val[[3]]$B,type="p",col=col_base[3],pch=16,cex=1.4)
  
  
  #graph_4
  plot(NA, NA, ylim = c(0,1), xlim = c(0,6),
       xlab = c("λT"), ylab =c("Bias") , las=1, cex.lab=1.4,font.lab=2, cex.axis=1.4,main="High")
  polygon(c(high_dist_SMM$lambda_T,rev(high_dist_SMM$lambda_T)),c(high_dist_SMM$min,rev(high_dist_SMM$max)),col="grey",border=NA)
  polygon(c(high_dist_SMM_delta$lambda_T,rev(high_dist_SMM_delta$lambda_T)),c(high_dist_SMM_delta$min,rev(high_dist_SMM_delta$max)),col="grey50",border=NA)
  points(ref_val[[4]]$lambda_T,ref_val[[4]]$B,type="p",col=col_base[4],pch=16,cex=1.4)
  

