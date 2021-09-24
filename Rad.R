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
library(pROC)

setwd("C:/local/R/FistulaStudy")
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
                34:36,38,40:44,46,48:51,53,55)]
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

predData1.stan.imp=mice(predData1.stan,maxit=15,m=50,
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
DAG_by_imputation=function(imp, bl, resp, predNames, respNames){
  
  # Get bootstrapped DAG for each version of imputed data, then average those DAGS
  
  dags=list()
  strengths=matrix(NA, imp$m, 200)
  
  
  for(i in 1:(imp$m)){
    df=data.frame(mice::complete(imp,i),resp)
    names(df)=c(predNames,respNames)
  
    dags[[i]]=averaged.network(boot.strength(df,R=200,algorithm="hc",algorithm.args=list(blacklist=bl)))
    
    # dag=boot.strength(df,R=200,algorithm="hc",algorithm.args=list(blacklist=bl))
    # strengths[i,]=boot.strength(df,R=200,algorithm="hc",algorithm.args=list(blacklist=bl))$strength
    
    
  }
  arcs=custom.strength(dags,nodes=names(df),cpdag=F)
  
  
  
  dag.avg=averaged.network(arcs)
  return(list(arcs, dag.avg))
}
get_blacklist=function(df,resp_indices,meds1_indices,meds2_indices,meds3_indices){
  
  # Blacklist arcs for each dataset (e.g., race does not affect age)
  # Also, assume non-demographic variables do not affect non-BMI demographic variables 
  # (e.g., high ESR does not affect gender)
  
  bl=tiers2blacklist(list(c("Age"),
                          c("Gender","Race","Smoker","BMI")))
  bl=rbind(bl,tiers2blacklist(list(c("Gender"),
                                   c("Age","Race","Smoker","BMI"))))
  bl=rbind(bl,tiers2blacklist(list(c("Race"),
                                   c("Age","Gender","Smoker","BMI"))))
  bl=rbind(bl,tiers2blacklist(list(c("Smoker"),c("BMI"))))
  bl=rbind(bl,tiers2blacklist(list(c("Intra.abdominal.Abscess.history.","Abscess.w.in.1yr.prior...Baseline."),
                                   c("Fistula."))))
  
  # Blacklist arcs from med variables to variables recorded before they were recorded
  bl=rbind(bl,tiers2blacklist(list(names(df)[-c(meds3_indices,resp_indices)],
                                   names(df)[meds3_indices])))
  bl=rbind(bl,tiers2blacklist(list(names(df)[-c(meds2_indices,meds3_indices,resp_indices)],
                                   names(df)[meds2_indices])))
  bl=rbind(bl,tiers2blacklist(list(names(df)[-c(meds1_indices,meds2_indices,meds3_indices,resp_indices)],
                                   names(df)[meds1_indices])))
  
  for(r in resp_indices)
    bl=rbind(bl,tiers2blacklist(list(names(df)[-r],
                                     names(df)[r])))
  
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
get_necessary_nodes=function(dag.strength,t){
  
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
lean_strength_plot=function(DAGS,type,thresh){
  
  # Prints only nodes in DAG with at least 1 sufficiently strong arc
  
  nn=get_necessary_nodes(DAGS[[type]],attr(DAGS[[type]],"threshold"))
  
  lean_strength=DAGS[[type]]
  lean_avg=DAGS[[type+5]]
  
  lean_strength=lean_strength[-nn[[2]],]
  lean_avg$nodes=lean_avg$nodes[nn[[1]]]
  
  return(strength.plot(lean_avg,lean_strength,shape="ellipse",threshold=thresh))
}
get_connected_nodes=function(avg_model,node_name){
  
  connected_nodes_names=c(node_name)
  i=0
  
  repeat{
    
    i=i+1
    
    l1=union(avg_model$nodes[[connected_nodes_names[i]]]$mb,
             avg_model$nodes[[connected_nodes_names[i]]]$nbr)
    l2=union(avg_model$nodes[[connected_nodes_names[i]]]$parents,
             avg_model$nodes[[connected_nodes_names[i]]]$children)
    connected_nodes_names=union(connected_nodes_names,union(l1,l2))
    
    if(i==length(connected_nodes_names))
      break
    
  }
  
  return(connected_nodes_names)
  
}
essential_nodes_plot=function(DAGS,dataset_num,response_names){
  
  nn=c()
  for(response in response_names)
    nn=union(nn,get_connected_nodes(DAGS[[dataset_num+5]],response))
  
  # nn=get_connected_nodes(DAGS[[dataset_num+5]],response_name)
  
  essential_strength_rows=which(DAGS[[dataset_num]][,1]%in%nn & 
                           DAGS[[dataset_num]][,2]%in%nn)
  essential_strength=DAGS[[dataset_num]][essential_strength_rows,]
  
  essential_avg_model=DAGS[[dataset_num+5]]
  essential_avg_model_node_indices=which(names(essential_avg_model$nodes)%in%nn)
  essential_avg_model$nodes=essential_avg_model$nodes[essential_avg_model_node_indices]
  essential_arc_rows=which(essential_avg_model$arcs[,1]%in%nn)
  essential_avg_model$arcs=essential_avg_model$arcs[essential_arc_rows,]
  
  return(strength.plot(essential_avg_model,essential_strength,shape="ellipse"))
  
}
# Original response DAGS #####

# All response DAGS #####
meds1_indices=18:23
meds2_indices=24:28
meds3_indices=29:33

all_response_dags_data=data.frame(predData1,respData.original[,1:5])
all_response_dags=rep(list(0),10)
all_response_dags_cc=list(which(complete.cases(all_response_dags_data)),
                          which(complete.cases(all_response_dags_data[,-c(9)])),
                          which(complete.cases(all_response_dags_data[,-c(9,13:17)])),
                          which(complete.cases(all_response_dags_data[,-c(13:17)]))
)
all_response_dags_bl=list(get_blacklist(all_response_dags_data[all_response_dags_cc[[1]],],1:5,34:38),
                          get_blacklist(all_response_dags_data[all_response_dags_cc[[2]],-c(9)],1:5,33:37),
                          get_blacklist(all_response_dags_data[all_response_dags_cc[[3]],-c(9,13:17)],1:5,28:32),
                          get_blacklist(all_response_dags_data[all_response_dags_cc[[4]],-c(13:17)],1:5,29:33),
                          get_blacklist(all_response_dags_data[all_response_dags_cc[[1]],],1:5,34:38)
)

set.seed(0)

all_response_dags[[1]]=boot.strength(all_response_dags_data[all_response_dags_cc[[1]],],R=200,algorithm="hc",
                                     algorithm.args=list(blacklist=all_response_dags_bl[[1]]))
all_response_dags[[2]]=boot.strength(all_response_dags_data[all_response_dags_cc[[2]],-c(9)],R=200,algorithm="hc",
                                     algorithm.args=list(blacklist=all_response_dags_bl[[2]]))
all_response_dags[[3]]=boot.strength(all_response_dags_data[all_response_dags_cc[[3]],-c(9,13:17)],R=200,algorithm="hc",
                                     algorithm.args=list(blacklist=all_response_dags_bl[[3]]))
all_response_dags[[4]]=boot.strength(all_response_dags_data[all_response_dags_cc[[4]],-c(13:17)],R=200,algorithm="hc",
                                     algorithm.args=list(blacklist=all_response_dags_bl[[4]]))
all_response_dags[[6]]=averaged.network(all_response_dags[[1]])
all_response_dags[[7]]=averaged.network(all_response_dags[[2]])
all_response_dags[[8]]=averaged.network(all_response_dags[[3]])
all_response_dags[[9]]=averaged.network(all_response_dags[[4]])

temp=DAG_by_imputation(predData1.stan.imp,
                       all_response_dags_bl[[5]],
                       respData.original[,1:5],
                       names(predData1.stan),
                       names(respData.original)[1:5])
all_response_dags[[5]] =temp[[1]]
all_response_dags[[10]]=temp[[2]]
rm(list=c("temp"))

essential_nodes_plot(all_response_dags,1,names(all_response_dags_data)[34:38]) # All predictors (88)
essential_nodes_plot(all_response_dags,2,names(all_response_dags_data)[34:38]) # All predictors except HBI (137)
essential_nodes_plot(all_response_dags,3,names(all_response_dags_data)[34:38]) # All predictors except labs and HBI (199)
essential_nodes_plot(all_response_dags,4,names(all_response_dags_data)[34:38]) # All predictors except labs (120)
essential_nodes_plot(all_response_dags,5,names(all_response_dags_data)[34:38]) # Imputed dataset (199 x 50)

# Get strengths for arcs to BinAbdSurg
filter(all_response_dags[[1]], from=="Steroids...20mg.within.2wks...34" & to=="BinAbdSurg")
filter(all_response_dags[[2]], from=="Steroids...20mg.within.2wks...34" & to=="BinAbdSurg")
filter(all_response_dags[[2]], from=="Fistula." & to=="BinAbdSurg")
filter(all_response_dags[[3]], from=="Fistula." & to=="BinAbdSurg")
filter(all_response_dags[[5]], from=="Fistula." & to=="BinAbdSurg")

# Get strengths for arcs to BinAbscess
filter(all_response_dags[[1]], from=="Abscess.w.in.1yr.prior...Baseline." & to=="Intra.abd.Abscess.")
filter(all_response_dags[[3]], from=="Fistula." & to=="Intra.abd.Abscess.")
filter(all_response_dags[[4]], from=="Abscess.w.in.1yr.prior...Baseline." & to=="Intra.abd.Abscess.")
filter(all_response_dags[[5]], from=="Fistula." & to=="Intra.abd.Abscess.")

# Get strengths for arcs to BinObstruction
filter(all_response_dags[[1]], from=="Steroids...20mg.within.2wks...34" & to=="BinObstructions")
filter(all_response_dags[[2]], from=="Steroids...20mg.within.2wks...34" & to=="BinObstructions")
filter(all_response_dags[[3]], from=="Stricture." & to=="BinObstructions")
filter(all_response_dags[[4]], from=="Steroids...20mg.within.2wks...34" & to=="BinObstructions")
filter(all_response_dags[[5]], from=="Stricture." & to=="BinObstructions")

# DAGGITY Dag Code
dag {
  AW1YP [adjusted,pos="0.580,-1.322"]
  AbdSurg_E3 [outcome,pos="-1.352,-0.008"]
  CRP [pos="-1.098,-2.077"]
  Fistula [exposure,pos="-0.511,-1.891"]
  IAAH [pos="0.123,-2.302"]
  IAB_E2 [outcome,pos="0.233,-0.018"]
  Inflammation [pos="-0.502,-1.086"]
  Obstruction_E2 [outcome,pos="-0.497,0.502"]
  Steroid34 [pos="-1.070,1.659"]
  Stricture [pos="-0.497,-0.322"]
  AW1YP -> IAB_E2
  CRP -> Fistula
  Fistula -> AbdSurg_E3
  Fistula -> IAB_E2
  Fistula -> Inflammation
  IAAH -> AW1YP
  IAAH -> Fistula
  Inflammation -> Stricture
  Steroid34 -> AbdSurg_E3
  Steroid34 -> IAB_E2
  Steroid34 -> Obstruction_E2
  Stricture -> Obstruction_E2
}
