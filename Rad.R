# Responses: abdominal surgeries, ostomies, abscesses, hospitalizations, obstructions
#            (abdominal surgeries, ostomies, abscesses), (hospitalizations, obstructions)
# Interested in presences and counts
# Predictors: demographic, radiological subjective, labs, medications

library(tibble)
library(readxl)
library(writexl)
library(caret)
library(dplyr)
library(mice)
library(bnlearn)
library(Rgraphviz)

setwd("~/R/FistulaStudy")
source=data.frame(read_excel("RadData2.xlsx"))
df=data.frame(read_excel("RadData2.xlsx"))   

########## Data Entry ##########
# Getting original responses #####

# respData.original=source[,c(60,62,66,68,69,71)]
respData.original=source[,c(62,66,68,69,71)]
respData.original[,5]=factor(respData.original[,5])

#respData.original$BinProcedures=as.numeric(respData.original[,1]!=0)
respData.original$BinAbdSurg=factor(as.numeric(respData.original[,1]!=0))
respData.original$BinOstomy=factor(as.numeric(respData.original[,2]!=0))
respData.original$BinHospitlizations=factor(as.numeric(respData.original[,3]!=0))
respData.original$BinObstructions=factor(as.numeric(respData.original[,4]!=0))

respData.original=respData.original[,c(6:9,5,1:4)]
names(respData.original)[6:9]=
  c("#AbdSurg","#Ostomy","#Hospitalization","#Obstructions")

respData.original.types=c(1,1,1,1,1,0,0,0,0)

# Take response data out of df
df=source[,-c(60,61,62,66,68,69,71)]

# Getting combined responses #####

respData.combined=data.frame(
  factor(as.numeric(rowSums((sapply(respData.original[,3:4],as.numeric)-1))!=0)),
  factor(as.numeric(rowSums((sapply(respData.original[,c(1,2,5)],as.numeric)-1))!=0)),
  rowSums(sapply(respData.original[,8:9],as.numeric)),
  rowSums(sapply(respData.original[,5:7],as.numeric)))

names(respData.combined)=c("BinHospitalizations","BinSurgeries",
                           "#Hospitalizations","#Surgeries")

respData.combined.types=c(1,1,0,0)

# Coding variables #####
# Code new manifestations to binary

# Ensure numeric columns are numeric
df.numVars=c(2,6,12,24,26,28,30,32,57:60,70,71,74,78,81,85,86)
for(j in df.numVars)
  df[,j]=as.numeric(df[,j])

# Ensure categorical columns are categorical
df.catvars=c(3:5,8:11,14:16,19:22,34:55,66,67,69,75,83,87:99)
for(j in df.catvars)
  df[,j]=factor(df[,j])

# PredData definition #####
predData1=df[,c(2:6,8:10,12,14,19,20,24,26,28,30,32,
                34:36,38,41:44,46,49:51,53,
                40,48,55)] # <-- Jack-inhibitors
# Variable exclusions include:
# variables after medication variables
# date variables

names(predData1)[4]="Smoker"

standardize=function(predData){
  for(j in 1:ncol(predData)){
    if(!is.factor(predData[,j])){
      nonNA.rows=which(!is.na(predData[,j]))
      nonNA.sub = predData[nonNA.rows,j]
      predData[nonNA.rows,j]=(nonNA.sub-mean(nonNA.sub))/sd(nonNA.sub)
    }
  }
  return (predData)
}


predData1.stan=standardize(predData1)

# Imputation of predData1.stan #####

predData1.stan.imp1=mice(predData1.stan,maxit=0,seed=10+07+2021)
predData1.stan.predM=predData1.stan.imp1$predictorMatrix
predData1.stan.meth=predData1.stan.imp1$method

predData1.stan.imp=mice(predData1.stan,maxit=15,m=200,
                   predictorMatrix=predData1.stan.predM,
                   method=predData1.stan.meth,
                   print=F,seed=10+07+2021)

rm(list=c("predData1.stan.imp1","predData1.stan.predM","predData1.stan.meth"))
# Unitab Helper functions #####
getFactorUniLRT=function(glm){
  return(1-pchisq(glm$null.deviance-glm$deviance,glm$df.null-glm$df.residual))
}
getUniModels=function(predData, respData, types, maxP) {
  
  # Returns table (and list of models) for each response where numeric variables
  # or at least 1 level for factor variables are significant at confidence maxP
  
  sigPreds=sigResps=pvals=ORs=N=c()
  row=m=1
  Models=list()
  resps=c()
  for(i in 1:ncol(predData)) {
    for(j in 1:ncol(respData)) {
      if(types[j]==1)
        mod=glm(respData[,j]~predData[,i], family=binomial())
      else mod=glm(respData[,j]~predData[,i], family=poisson())
    
      if(!is.factor(predData[,j])){
        if(coef(summary(mod))[2,4]<=maxP) {
          sigPreds[row]=paste(names(predData)[i], " (", i, ")", sep="")
          sigResps[row]=paste(names(respData)[j], " (", j, ")", sep="")
          pvals[row]=coef(summary(mod))[2,4]
          ORs[row]=exp(coef(summary(mod))[2,1])
          N[row]=sum(complete.cases(predData[,i]))
          row=row+1
          Models[[m]]=mod
          m=m+1
        }
      }
      # If predictor is factor, return in table if at least one factor level has p <= maxP
      if(is.factor(predData[,j])){
        if((length(levels(predData[,j]))==2 & coef(summary(mod))[2,4]<=maxP) | 
           getFactorUniLRT(mod)<=maxP) {
          sigPreds[row]=paste(names(predData)[i], " (", i, ")", sep="")
          sigResps[row]=paste(names(respData)[j], " (", j, ")", sep="")
          pvals[row]=coef(summary(mod))[2,4]
          ORs[row]=exp(coef(summary(mod))[2,1])
          N[row]=sum(complete.cases(predData[,i]))
          row=row+1
          Models[[m]]=mod
          m=m+1
        }
      }
    }
  }
  List=list(data.frame(sigPreds,sigResps,pvals,ORs,N),Models,resps)
  return(List)
}
getIndices=function(strings){
  s2=read.table(text=strings,sep="(",as.is=T)$V2
  return(as.numeric(read.table(text=s2,sep=")",as.is=T)$V1))
}
getPredIndices=function(unitab,responseString){
  tab=unitab[which(unitab[,2]==responseString),]
  return(getIndices(tab[,1]))
}
getUniModels.imp=function(pred,pred.imp,resp,types,maxP){
  
  # Same as uniModels function above, but for pooled models on imputed dataset
  
  sigPreds=sigResps=pvals=ORs=c()
  row=m=1
  Models=list()
  resps=c()
  for(i in 1:ncol(pred)){
    for(j in 1:ncol(resp)){
      if(types[j]==1)
        mod=with(pred.imp, glm(resp[,j]~., data=data.frame(pred[,i]), family=binomial()))
      else mod=with(pred.imp, glm(resp[,j]~., data=data.frame(pred[,i]), family=poisson()))
      if(!is.factor(pred[,j])){
        if(summary(pool(mod))[2,6]<=maxP){
          sigPreds[row]=paste(names(pred)[i], " (", i, ")", sep="")
          sigResps[row]=paste(names(resp)[j], " (", j, ")", sep="")
          pvals[row]=(summary(pool(mod)))[2,6]
          ORs[row]=exp(summary(pool(mod))[2,2])
          row=row+1
          Models[[m]]=mod
          m=m+1
        }
      }
      if(is.factor(pred[,j])){
        
        pvals.mod=summary(pool(mod))[-1,6]
        condition=sum(as.numeric(pvals<=maxP))!=0
        
        # condition is equivalent to at least one factor level having p <= maxP
        
        if(condition){
          sigPreds[row]=paste(names(pred)[i], " (", i, ")", sep="")
          sigResps[row]=paste(names(resp)[j], " (", j, ")", sep="")
          pvals[row]=summary(pool(mod))[2,6]
          ORs[row]=exp(summary(pool(mod))[2,2])
          row=row+1
          Models[[m]]=mod
          m=m+1
        }
      }
    }
  }
  List=list(data.frame(sigPreds,sigResps,pvals,ORs),Models,resps)
  return(List)
}
getFistulaTable=function(fistula,resp){
  # Equivalent to confusion matrix using fistula as predictor and resp as reference
  return(table(fistula,resp))
}
# QUESTIONS #####

# Getting table of significant univariate models for predData1 on non-combined responses #####

unitab1=getUniModels(predData1,respData.original,respData.original.types,.05)[[1]]
unitab1.imp=getUniModels.imp(predData1,predData1.imp2,respData.original,respData.original.types,.05)[[1]]

# Predictors #####

predData1.demPreds=c(1:5)
predData1.radPreds=c(6:8)
predData1.subPreds=c(9:10)
predData1.labPreds=c(13:17)
predData1.medPreds=c(18:30)

preds.include=c(predData1.radPreds,
                predData1.labPreds)


# DAG Helper functions #####
DAG_by_imputation=function(imp, bl, resp, predNames, respName){
  
  # Get bootstrapped DAG for each version of imputed data, then average those DAGS
  
  dags=list()
  for(i in 1:(imp$m)){
    df=data.frame(mice::complete(imp,i),resp)
    names(df)=c(predNames,respName)
    dags[[i]]=averaged.network(boot.strength(df,R=30,algorithm="hc",algorithm.args=list(blacklist=bl)))
  }
  arcs=custom.strength(dags,nodes=names(df),cpdag=F)
  dag.avg=averaged.network(arcs)
  return(list(arcs, dag.avg))
}
get_blacklist=function(df, dem_indices){
  
  # Blacklist arcs for each dataset (e.g., race does not affect age)
  # Also, assume non-demographic variables do not affect non-BMI demographic variables 
  # (e.g., high ESR does not affect gender)
  
  bl=tiers2blacklist(list(names(df)[dem_indices[-length(dem_indices)]],
                          names(df)[-dem_indices[-length(dem_indices)]]))
  bl=rbind(bl,tiers2blacklist(list(c("Age"),
                                   c("Gender","Race","Smoker","BMI"))))
  bl=rbind(bl,tiers2blacklist(list(c("Gender"),
                                   c("Age","Race","Smoker","BMI"))))
  bl=rbind(bl,tiers2blacklist(list(c("Race"),
                                   c("Age","Gender","Smoker","BMI"))))
  bl=rbind(bl,tiers2blacklist(list(c("Smoker"),c("BMI"))))
  bl=rbind(bl,tiers2blacklist(list(names(df)[-ncol(df)],
                                   names(df)[ncol(df)])))
  return(bl)
  
}
get_DAGS=function(pred, resp, respName, HBI_Lab_cols, imp, seed){
  
  # Get DAGS from datasets
  # 1) Dataset of all variables (88 cases)
  # 2) Dataset with all variables except HBI
  # 3) Dataset with all variables except HBI and labs
  # 4) Dataset with all variables except labs
  # 5) Imputed dataset
  
  # For each dataset, return
  # 1) Bootstrap of DAGS (200 iterations for non-imputed dataset DAGS)
  # 2) Average of bootstrap networks
  
  
  set.seed(seed)
  
  dags=list()
  for(i in 1:5){
    
    if(i==1){
      cases=which(complete.cases(pred))
      df=data.frame(pred[cases,],resp[cases])
    }
    else if(i==2){
      cases=which(complete.cases(pred[,-HBI_Lab_cols[1]]))
      df=data.frame(pred[cases,-HBI_Lab_cols[1]],resp[cases])
    }
    else if(i==3){
      cases=which(complete.cases(pred[,-HBI_Lab_cols]))
      df=data.frame(pred[cases,-HBI_Lab_cols],resp[cases])
    }
    else if(i==4){
      cases=which(complete.cases(pred[,-HBI_Lab_cols[-1]]))
      df=data.frame(pred[cases,-HBI_Lab_cols[-1]],resp[cases])
    }
    else if(i==5){
      df=data.frame(pred, resp)
    }
    
    names(df)[ncol(df)]=respName
    bl=get_blacklist(df,1:5)
    
    if(i!=5){
      dags[[i]]=boot.strength(df,R=200,algorithm="hc",algorithm.args=list(blacklist=bl))
      dags[[i+5]]=averaged.network(dags[[i]])
    }
    else {
      temp=DAG_by_imputation(imp,bl,resp,names(df)[-ncol(df)],names(df)[ncol(df)])
      dags[[i]]=temp[[1]]
      dags[[i+5]]=temp[[2]]
    }
  }
  return(dags)
}
get_necessary_nodes=function(dag.strength, t){
  
  # Finds nodes in given DAG bootstrap with strength above threshold t
  # Returns
  # 1) Names of nodes with an arc of strength >= t
  # 2) Rows of dag.strength with strength < t
  
  nodes=c()
  indices_to_remove_from_boot=c()
  itrfb=1
  for(i in 1:nrow(dag.strength)){
    if(dag.strength[i,3]>=t)
      nodes=union(nodes,c(dag.strength[i,1],dag.strength[i,2]))
    else {
      indices_to_remove_from_boot[itrfb]=i
      itrfb=itrfb+1
    } 
  }
  return(list(nodes,indices_to_remove_from_boot))
}
lean_strength_plot=function(DAGS,type){
  
  # Prints only nodes in DAG with at least 1 sufficiently strong arc
  
  nn=get_necessary_nodes(DAGS[[type]],attr(DAGS[[type]],"threshold"))
  
  lean_strength=DAGS[[type]]
  lean_avg=DAGS[[type+5]]
  
  lean_strength=lean_strength[-nn[[2]],]
  lean_avg$nodes=lean_avg$nodes[nn[[1]]]
  
  return(strength.plot(lean_avg,lean_strength,shape="ellipse"))
}
# Original response DAGS #####

# BinAbdSurg DAGs
BinAbdSurg.DAGS=get_DAGS(predData1.stan,respData.original[,1],"BinAbdSurgeries",c(9,13:17),predData1.stan.imp,1)
lean_strength_plot(BinAbdSurg.DAGS,1) # Complete cases
lean_strength_plot(BinAbdSurg.DAGS,2) # Missing only HBI
lean_strength_plot(BinAbdSurg.DAGS,3) # Missing HBI and labs
lean_strength_plot(BinAbdSurg.DAGS,4) # Missing labs
lean_strength_plot(BinAbdSurg.DAGS,5) # Imputed dataset

# For each response, create multivariate model wtih variables connected to response in at least 1 DAG

# Fistula, Steroid34, CRP (causal to fistula)
BinAbdSurg.DAGS.glm=glm(respData.original[,1]~.,data=predData1.stan[,c(6,18)],family=binomial())
summary(BinAbdSurg.DAGS.glm)

# Fistula, Steroids34 still significant when considering stricture, inflammation, labs, abdPain
summary(glm(respData.original[,1]~.,data=predData1.stan[,c(6,18,13,7:8,14:17,10)],family=binomial()))

BinAbdSurg.DAGS.glm.imp=pool(with(predData1.stan.imp,glm(respData.original[,1]~.,
                                                         data=predData1.stan[,c(6,18,13)],
                                                         family=binomial())))
summary(BinAbdSurg.DAGS.glm.imp)

# BinOstomy DAGS
BinOstomy.DAGS=get_DAGS(predData1.stan,respData.original[,2],"BinOstomy",c(9,13:17),predData1.stan.imp,2)
lean_strength_plot(BinOstomy.DAGS,1)
lean_strength_plot(BinOstomy.DAGS,2)
lean_strength_plot(BinOstomy.DAGS,3)
lean_strength_plot(BinOstomy.DAGS,4)
lean_strength_plot(BinOstomy.DAGS,5)

# Variables taken from unitab (DAGs say no relation)
# CRP, ESR, Biologics46
BinOstomy.DAGS.glm=glm(respData.original[,2]~.,data=predData1.stan[,c(13,14,26)],family=binomial())
summary(BinOstomy.DAGS.glm)

# BinHospitalization DAGS
BinHospitlization.DAGS=get_DAGS(predData1.stan,respData.original[,3],"BinHospitalization",c(9,13:17),predData1.stan.imp,3)
lean_strength_plot(BinHospitlization.DAGS,1)
lean_strength_plot(BinHospitlization.DAGS,2)
lean_strength_plot(BinHospitlization.DAGS,3)
lean_strength_plot(BinHospitlization.DAGS,4)
lean_strength_plot(BinHospitlization.DAGS,5)

# Fistula, Steroids34, HBI
BinHospitlization.DAGS.glm=glm(respData.original[,3]~.,data=predData1.stan[,c(6,18,9)],family=binomial())
summary(BinHospitlization.DAGS.glm)

summary(glm(respData.original[,3]~.,data=predData1.stan[,c(6,18,9,13:17)],family=binomial()))

BinHospitlization.DAGS.glm.imp=pool(with(predData1.stan.imp,glm(respData.original[,3]~.,
                                                                data=predData1.stan[,c(6,18,9)],
                                                                family=binomial())))
summary(BinHospitlization.DAGS.glm.imp)

# BinObstructions DAGS
BinObstruction.DAGS=get_DAGS(predData1.stan,respData.original[,4],"BinObstruction",c(9,13:17),predData1.stan.imp,4)
lean_strength_plot(BinObstruction.DAGS,1)
lean_strength_plot(BinObstruction.DAGS,2)
lean_strength_plot(BinObstruction.DAGS,3)
lean_strength_plot(BinObstruction.DAGS,4)
lean_strength_plot(BinObstruction.DAGS,5)

# Stricture, Steroids34
BinObstruction.DAGS.glm=glm(respData.original[,4]~.,data=predData1.stan[,c(7,18)],family=binomial())
summary(BinObstruction.DAGS.glm)

BinObstruction.DAGS.glm.imp=pool(with(predData1.stan.imp,glm(respData.original[,4]~.,
                                                             data=predData1.stan[,c(7,18)],
                                                             family=binomial())))
summary(BinObstruction.DAGS.glm.imp)

# BinAbscess DAGS
BinAbscess.DAGS=get_DAGS(predData1.stan,respData.original[,5],"BinAbscess",c(9,13:17),predData1.stan.imp,5)
lean_strength_plot(BinAbscess.DAGS,1)
lean_strength_plot(BinAbscess.DAGS,2)
lean_strength_plot(BinAbscess.DAGS,3)
lean_strength_plot(BinAbscess.DAGS,4)
lean_strength_plot(BinAbscess.DAGS,5)

# Fistula, Abscess w/i 1 year prior
BinAbscess.DAGS.glm=glm(respData.original[,5]~.,data=predData1.stan[,c(6,12)],family=binomial())
summary(BinAbscess.DAGS.glm)

# NumAbdSurg DAGS
NumAbdSurg.DAGS=get_DAGS(predData1.stan,respData.original[,6],"NumAbdSurgeries",c(9,13:17),predData1.stan.imp,6)
lean_strength_plot(NumAbdSurg.DAGS,1)
lean_strength_plot(NumAbdSurg.DAGS,2)
lean_strength_plot(NumAbdSurg.DAGS,3)
lean_strength_plot(NumAbdSurg.DAGS,4)
lean_strength_plot(NumAbdSurg.DAGS,5)

# Fistula, Immunodmodulators 51
NumAbdSurg.DAGS.glm=glm(respData.original[,6]~.,data=predData1.stan[,c(6,29)],family=poisson())
summary(NumAbdSurg.DAGS.glm)

# Fistula is still significant after including radiological and lab predictors
summary(glm(respData.original[,6]~.,data=predData1.stan[,c(6,29,7:8,13:17)],family=poisson()))

# NumOstomy DAGS
NumOstomy.DAGS=get_DAGS(predData1.stan,respData.original[,7],"NumOstomy",c(9,13:17),predData1.stan.imp,7)

# NumHospitalization DAGS
NumHospitalization.DAGS=get_DAGS(predData1.stan,respData.original[,8],"NumHospitalization",c(9,13:17),predData1.stan.imp,8)
lean_strength_plot(NumHospitalization.DAGS,1)
lean_strength_plot(NumHospitalization.DAGS,2)
lean_strength_plot(NumHospitalization.DAGS,3)
lean_strength_plot(NumHospitalization.DAGS,4)
lean_strength_plot(NumHospitalization.DAGS,5)

# Fistula, HBI, Steroids 34
NumHospitalization.DAGS.glm=glm(respData.original[,8]~.,data=predData1.stan[,c(9,18,13)],family=poisson())
summary(NumHospitalization.DAGS.glm)

# NumObstruction DAGS
NumObstruction.DAGS=get_DAGS(predData1.stan,respData.original[,9],"NumObstruction",c(9,13:17),predData1.stan.imp,9)
lean_strength_plot(NumObstruction.DAGS,1)
lean_strength_plot(NumObstruction.DAGS,2)
lean_strength_plot(NumObstruction.DAGS,3)
lean_strength_plot(NumObstruction.DAGS,4)
lean_strength_plot(NumObstruction.DAGS,5)

# Steroids 34, Stricture
NumObstruction.DAGS.glm=glm(respData.original[,9]~.,data=predData1.stan[,c(18,7)],family=poisson())
summary(NumObstruction.DAGS.glm)

# Combined response DAGS #####

# Combined BinAbdSurg DAGS
BinAbdSurg.c.DAGS=get_DAGS(predData1.stan,respData.combined[,1],"BinAbdSurgeries",c(9,13:17),predData1.stan.imp,11)
lean_strength_plot(BinAbdSurg.c.DAGS,1)
lean_strength_plot(BinAbdSurg.c.DAGS,2)
lean_strength_plot(BinAbdSurg.c.DAGS,3)
lean_strength_plot(BinAbdSurg.c.DAGS,4)

BinAbdSurg.c.DAGS.glm=glm(respData.combined[,1]~.,data=predData1.stan[,c(6,18)],family=binomial())
summary(BinAbdSurg.c.DAGS.glm)

# Combined BinHospitalization DAGS
BinHospitlization.c.DAGS=get_DAGS(predData1.stan,respData.combined[,1],"BinHospitalization",c(9,13:17),predData1.stan.imp,12)
lean_strength_plot(BinHospitlization.c.DAGS,1,.2)
lean_strength_plot(BinHospitlization.c.DAGS,2,.2)
lean_strength_plot(BinHospitlization.c.DAGS,3,.2)
lean_strength_plot(BinHospitlization.c.DAGS,4,.2)
lean_strength_plot(BinHospitlization.c.DAGS,5,.2)

# Fistula tables #####

# Original responses
unitab1[which(unitab1[,1]=="Fistula. (6)"),]

getFistulaTable(predData1[,6],respData.original[,1])
getFistulaTable(predData1[,6],respData.original[,2])
getFistulaTable(predData1[,6],respData.original[,3])
getFistulaTable(predData1[,6],respData.original[,4])
getFistulaTable(predData1[,6],respData.original[,5])
getFistulaTable(predData1[,6],respData.original[,6])
getFistulaTable(predData1[,6],respData.original[,7])
getFistulaTable(predData1[,6],respData.original[,8])
getFistulaTable(predData1[,6],respData.original[,9])

# Combined responses

getFistulaTable(predData1[,6],respData.combined[,1])
getFistulaTable(predData1[,6],respData.combined[,2])
getFistulaTable(predData1[,6],respData.combined[,3])
getFistulaTable(predData1[,6],respData.combined[,4])

#####
source.dateVars=c(7,13,17,18,25,27,29,31,33,56,64,65,67,70,72,75,83,86,89,91)
rbind(source.dateVars,as.vector(colSums(!is.na(source[,source.dateVars]))))

predData=df[,c(2:6,8:12,14:16,18:20,22,24,26,28,30,32,34:56,)]

respCountIndices=c(62,66,68,69,71,74)

ostomyBin=factor(as.numeric(df[,66]!=0)) # Ostomy > 0 factor
hospitBin=factor(as.numeric(df[,68]!=0)) # Hospitalization > 0 factor
obstruBin=factor(as.numeric(df[,69]!=0)) # Obstruction > 0 factor

respData=cbind(df[,61],ostomyBin,hospitBin,obstruBin,df[,71],
               df[,c(74,62,66,68,69)])

names(respData)[c(1,5,7:10)]=c("AbdSurgBin","AbscessBin","#AbdSurg",
                               "#Ostomies","#Hospitalizations","#TreatedObstructions")

# Just putting responses to back of DF
df=df[,-c(61,respCountIndices)]
df=cbind(df, respData)

# Fixing numeric columns stored as chars
df[,67]=as.numeric(df[,67]!="0")
df[,87]=as.numeric(df[,87]!="0") # Original had "NO CHROHS" at r198

v=c(1,48,70,87,125,172,184)
wm2 = rep(NA, times=199)
wm2[v]=c(7,6,7,6,6,7,6)
df[v,23]=c("2","5","1","1","1","1","1")
df=cbind(df, wm2)
df=df[,c(1:23,110,24:109)]

# Variable indices with categorical (nominal) interpretations
catVars=c(3:5,8:11,14:16,19:24,35:56,68,70,72,76,84,88:106)

# Insures they are stored as nominal variables
for (j in catVars){
  temp=factor(df[,j])
  df[,j]=addNA(df[,j])
  levels(df[,j])=c(levels(temp),0)
}

# Variables with numeric interpretation
numVars=c(1:2,6,12,25,27,29,31,33,58:62,
          71:72,75,79,82,86:87,107:110)

for(j in numVars){
  df[,j]=as.numeric(df[,j])
}

# Get rid of deceased NA observation
df[199,106]=respData[199,6]=0

# # of negative outcomes for each patient
numOutcomes=rowSums(respData[6:10])
df=cbind(df, numOutcomes)
respData=cbind(respData, numOutcomes)
predData=df[,-c(7,73,78,81,101:112)]

# Combine surgeries and ostomies
surg.n = as.numeric(respData[,1])
ost.n  = as.numeric(respData[,2])-1

surgPost.Bin=as.factor(as.numeric(surg.n | ost.n))
surgPost.Num = rowSums(respData[,7:8])

respData=cbind(respData,as.factor(surgPost.Bin),surgPost.Num)
names(respData)[12]="surgPost.Bin"
