# Responses: abdominal surgeries, ostomies, abscesses, hospitalizations, obstructions
#            (abdominal surgeries, ostomies, abscesses), (hospitalizations, obstructions)
# Interested in presences and counts
# Predictors: demographic, radiological subjective, labs, medications

library(tibble)
library(readxl)
library(writexl)
library(caret)
library(dplyr) # filter
library(mice)
library(bnlearn)
library(Rgraphviz)
library(pROC)

df=data.frame(read_excel('C:/_local/fistula_study/RadData2.xlsx'))

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
df=source[,-c(61,62,66,68,69,71)]

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

# Imputation of predData1 #####
predData1.imp1=mice(predData1,maxit=0,seed=10+07+2021)
predData1.predM=predData1.imp1$predictorMatrix
predData1.meth=predData1.imp1$method

predData1.imp=mice(predData1,maxit=25,m=50,
                   predictorMatrix=predData1.predM,
                   method=predData1.meth,
                   print=F,seed=10+07+2021)

rm(list=c("predData1.imp1","predData1.predM","predData1.meth"))
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
  strengths=matrix(NA, nrow=702, ncol=imp$m)
  directions=matrix(NA, nrow=702, ncol=imp$m)
  
  for(i in 1:(imp$m)){
    df=data.frame(mice::complete(imp,i),resp)
    names(df)=c(predNames,respNames)
    
    bs_=boot.strength(df,R=200,algorithm="hc",algorithm.args=list(blacklist=bl))
  
    strengths[,i]=bs_[,3]
    directions[,i]=bs_[,4]
    
    #dags[[i]]=averaged.network(boot.strength(df,R=200,algorithm="hc",algorithm.args=list(blacklist=bl)))
    #dag=boot.strength(df,R=200,algorithm="hc",algorithm.args=list(blacklist=bl))

  }

  bs_$strength=rowMeans(strengths)
  bs_$direction=rowMeans(directions)
  
  return(bs_)  
  
  #arcs=custom.strength(dags,nodes=names(df),cpdag=F)
  #dag.avg=averaged.network(arcs)
  #return(list(arcs, dag.avg))
  
}
get_blacklist=function(df,resp_indices, dem_indices,
                       #meds1_indices,meds2_indices,meds3_indices
                       meds_indices){
  
  # Blacklist arcs for each dataset (e.g., race does not affect age)
  # Also, assume non-demographic variables do not affect non-BMI demographic variables 
  # (e.g., high ESR does not affect gender)
  
  bl=tiers2blacklist(list(c("Age"),
                          c("Gender","Race","Smoker","BMI")))
  bl=rbind(bl,tiers2blacklist(list(c("Gender"),
                                   c("Age","Race","Smoker","BMI"))))
  bl=rbind(bl,tiers2blacklist(list(c("Race"),
                                   c("Age","Gender","Smoker","BMI"))))
  bl=rbind(bl,tiers2blacklist(list(c("Smoker"),
                                   c("BMI"))))
  bl=rbind(bl,tiers2blacklist(list(c(names(df)[dem_indices]),
                                   c(names(df)[-dem_indices]))))
  
  bl=rbind(bl,tiers2blacklist(list(c("Intra.abdominal.Abscess.history.","Abscess.w.in.1yr.prior...Baseline."),
                                   c("Fistula."))))
  bl=rbind(bl,tiers2blacklist(list(c("Stricture.","Inflammation."),
                                   c("Fistula."))))
  bl=rbind(bl,(tiers2blacklist(list(c("Inflammation."),
                                    c("Stricture.")))))
  
  #bl=rbind(bl,tiers2blacklist(list(names(df)[-c(meds3_indices,resp_indices)],
  #                                 names(df)[meds3_indices])))
  #bl=rbind(bl,tiers2blacklist(list(names(df)[-c(meds2_indices,meds3_indices,resp_indices)],
  #                                 names(df)[meds2_indices])))
  #bl=rbind(bl,tiers2blacklist(list(names(df)[-c(meds1_indices,meds2_indices,meds3_indices,resp_indices)],
  #                                 names(df)[meds1_indices])))
  bl=rbind(bl,tiers2blacklist(list(c(names(df)[-c(meds_indices,resp_indices)]),
                                   c(names(df)[meds_indices]))))
  
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
    bl=get_blacklist(df,)
    
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
essential_nodes_plot_noT=function(DAGS,dataset_num,response_names){
  
  avg_dag=averaged.network(DAGS[[dataset_num]])
  nn=c()
  for(response in response_names)
    nn=union(nn,get_connected_nodes(avg_dag,response_names))
  
  essential_strength_rows=which(DAGS[[dataset_num]][,1]%in%nn & 
                                DAGS[[dataset_num]][,2]%in%nn
                                )
  essential_strength=DAGS[[dataset_num]][essential_strength_rows,]
  
  essential_avg_model=avg_dag
  print(paste('Threshold = ', avg_dag$learning$args[[1]]))
  
  essential_avg_model_node_indices=which(names(essential_avg_model$nodes)%in%nn)
  essential_avg_model$nodes=essential_avg_model$nodes[essential_avg_model_node_indices]
  essential_arc_rows=which(essential_avg_model$arcs[,1]%in%nn)
  essential_avg_model$arcs=essential_avg_model$arcs[essential_arc_rows,]
  
  
  illegal_df = data.frame(essential_avg_model$learning$illegal)
  essential_illegal_indices = which(illegal_df[,1]%in%nn & illegal_df[,2]%in%nn)

  essential_avg_model$learning$illegal = illegal_df[essential_illegal_indices,]
  return (list(essential_avg_model, essential_strength))
  
#  (strength.plot(essential_avg_model,
#                       essential_strength,
#                       shape="ellipse"
#                )
#   )
  
}
essential_nodes_plot_withT=function(DAGS,dataset_num,response_names,thresh){
  
  avg_dag=averaged.network(DAGS[[dataset_num]],threshold=thresh)
  nn=c()
  for (response in response_names)
    nn=union(nn,get_connected_nodes(avg_dag,response_names))
  
  essential_stength_rows=which(DAGS[[dataset_num]][,1]%in%nn & 
                               DAGS[[dataset_num]][,2]%in%nn &
                               DAGS[[dataset_num]][,3] >= thresh)
  essential_strength=DAGS[[dataset_num]][essential_stength_rows,]
  
  essential_avg_model=avg_dag
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

func = function(data, ind){
  return(1*(rowSums(apply(data[,ind],2,as.numeric))>=1))
}

# Combining meds for table
meds_combined=data.frame(Steorids=factor(func(predData1,c(18,24,29))),
                         Mesalamine=factor(func(predData1,c(19,25,30))),
                         Immunomodulators=factor(func(predData1,c(20,26,31))),
                         Biologics=factor(func(predData1,c(21,27,32))),
                         Jack_inhibitors=factor(func(predData1,c(22,28,33)))
)


all_response_dags_data=data.frame(predData1[,1:17],meds_combined,respData.original[,1:5])

# Recent-surgery-cases
rcs_cases = which(as.numeric(df[,16])-1==1)

all_response_dags_data = all_response_dags_data[-rcs_cases, ]

# Abscess at time of fistula
abscess_df = df[df[,71]==1,]
abscess_df = abscess_df[,c(7,72)]
abscess_at_imaging_cases = c(2,8,26,28,36,39,162)
all_response_dags_data = all_response_dags_data[-abscess_at_imaging_cases,]

# all_response_dags_data=all_response_dags_data[,c(1:22,31,23:26,32,27:30,33,34:38)]

# Imputation of data
all_response_dags_data.imp1=mice(all_response_dags_data[,1:22],maxit=0,seed=10+07+2021)
all_response_dags_data.predM=all_response_dags_data.imp1$predictorMatrix
all_response_dags_data.meth=all_response_dags_data$method

all_response_dags_data.imp=mice(all_response_dags_data[,1:22],maxit=25,m=50,
                                predictorMatrix=all_response_dags_data.predM,
                                method=all_response_dags_data.meth,
                                print=F,seed=10+07+2021)

rm(list=c('all_response_dags_data.imp1','all_response_dags_data.predM','all_response_dags_data.meth'))

all_response_dags=rep(list(0),5)
all_response_dags_cc=list(which(complete.cases(all_response_dags_data)),
                          which(complete.cases(all_response_dags_data[,-c(9)])),
                          which(complete.cases(all_response_dags_data[,-c(9,13:17)])),
                          which(complete.cases(all_response_dags_data[,-c(13:17)])),
                          which(complete.cases(all_response_dags_data))
                         )

all_response_dags_bl=list(get_blacklist(all_response_dags_data[all_response_dags_cc[[1]],],23:27,1:5,18:22),
                          get_blacklist(all_response_dags_data[all_response_dags_cc[[2]],-c(9)],22:26,1:5,17:21),
                          get_blacklist(all_response_dags_data[all_response_dags_cc[[3]],-c(9,13:17)],17:21,1:5,12:16),
                          get_blacklist(all_response_dags_data[all_response_dags_cc[[4]],-c(13:17)],18:22,1:5,13:17),
                          get_blacklist(all_response_dags_data[all_response_dags_cc[[1]],],23:27,1:5,18:22)
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
all_response_dags[[5]]=DAG_by_imputation(all_response_dags_data.imp,
                                         all_response_dags_bl[[5]],
                                         respData.original[-abscess_at_imaging_cases,1:5],
                                         names(all_response_dags_data)[1:22],
                                         names(respData.original)[1:5]
                                         )

all_response_dags[[6]]=averaged.network(all_response_dags[[1]])
all_response_dags[[7]]=averaged.network(all_response_dags[[2]])
all_response_dags[[8]]=averaged.network(all_response_dags[[3]])
all_response_dags[[9]]=averaged.network(all_response_dags[[4]])
all_response_dags[[10]]=averaged.network(all_response_dags[[5]])

# Get fistula table

table_preds=c(2:10,12,14:15,18:20,22,24,26,28,30,32,
              34:36,38,40:44,46,48:51,53,55,57:60,73,94)
stats_df = df[,table_preds]

stats_df[,6] = as.Date(stats_df[,6])
stats_df[,13] = as.Date(stats_df[,13])
stats_df[, 'years_since_diagnosis'] = as.numeric((stats_df[,6] - stats_df[,13])/365)

negative_years = which(stats_df[,'years_since_diagnosis']<=0)
stats_df[negative_years, 'years_since_diagnosis'] = 0

stats_df = stats_df[,-c(6,13)]
stats_df = stats_df[-abscess_at_imaging_cases, ]

stats_df[,40][stats_df[,40]!='0'] = '1'
stats_df[,40] = as.factor(stats_df[,40])

stats_df[,41] = as.character(stats_df[,41])
stats_df[,41][stats_df[,41]!='0'] = '1'
stats_df[,41] = as.factor(stats_df[,41])

numeric_vars=c(1,5,9,15:19,37:39,42)
for (n in numeric_vars)
   stats_df[,n] = as.numeric(stats_df[,n])

fistula_cases = which(stats_df[,6]==1)
pred = pop = fistu = nfist = P = c()

for (v in 1:ncol(stats_df)){
   
   if(names(stats_df)[v]=='Fistula.')
      next
   
   if(!is.factor(stats_df[,v])){
      
      pred[v] = paste(names(stats_df)[v], ' (SD)')
      
      mean = round(mean(stats_df[,v], na.rm=T), 3)
      sd = round(sd(stats_df[,v], na.rm=T), 3)
      
      f_mean = round(mean(stats_df[fistula_cases,v], na.rm=T),3)
      f_sd = round(sd(stats_df[fistula_cases, v], na.rm=T),3)
      
      n_mean = round(mean(stats_df[-fistula_cases, v], na.rm=T),3)
      n_sd = round(sd(stats_df[-fistula_cases, v], na.rm=T),3)
      
      pop[v] = paste(mean, ' (', sd, ')')
      fistu[v] = paste(f_mean, ' (', f_sd, ')')
      nfist[v] = paste(n_mean, ' (', n_sd, ')')
      
      mod = glm(stats_df[,'Fistula.']~stats_df[,v], family=binomial)
      P[v] = summary(mod)$coefficients[2, 4]
      
   }
   
   else {
      
      p_table = data.frame(table(stats_df[, v]))
      f_table = data.frame(table(stats_df[fistula_cases, v]))
      n_table = data.frame(table(stats_df[-fistula_cases, v]))
   
      p_order = order(p_table[, 'Var1'])
      f_order = order(f_table[, 'Var1'])
      n_order = order(n_table[, 'Var1'])
   
      f_table = f_table[f_order, ]
      n_table = n_table[n_order, ]
      
      pred_string = names(stats_df)[v]
      p_string = ''
      f_string = ''
      n_string = ''
      
      if (nrow(p_table)!=2) {
         
         for (r in 1:nrow(f_table)){
            
            pred_string = paste(pred_string, f_table[r, 'Var1'], '/')
            p_string = paste(p_string, p_table[r, 'Freq'], ':')
            f_string = paste(f_string, f_table[r, 'Freq'], ':')
            n_string = paste(n_string, n_table[r, 'Freq'], ':')
            
         }
         
      }
      
      else {
         
         p_string = paste(p_table[2, 2], ' (', round(100*p_table[2,2]/sum(p_table[,2]), 4), '%)')
         f_string = paste(f_table[2, 2], ' (', round(100*f_table[2,2]/sum(f_table[,2]), 4), '%)')
         n_string = paste(n_table[2, 2], ' (', round(100*n_table[2,2]/sum(n_table[,2]), 4), '%)')
         
      }
      
      
      
      cases = which(complete.cases(stats_df[,v]))
      mod = glm(stats_df[cases, 'Fistula.']~stats_df[cases, v], family=binomial)
      null = glm(stats_df[cases, 'Fistula.']~1, family=binomial)
      
      pred[v] = pred_string
      pop[v] = p_string
      fistu[v] = f_string
      nfist[v] = n_string
      P[v] = anova(mod, null, test='LRT')$`Pr(>Chi)`[2]
      
   }

}

new_stats_table = data.frame(pred, pop, fistu, nfist, P)
write_xlsx(new_stats_table, path='C:/_local/fistula_study/stats_no_abscess_at_imaging.xlsx')

essential_nodes_plot_noT(all_response_dags,1,names(all_response_dags_data)[23:27]) # All predictors (88)
#essential_nodes_plot_noT(all_response_dags,2,names(all_response_dags_data)[23:27]) # All predictors except HBI (137)
essential_nodes_plot_noT(all_response_dags,3,names(all_response_dags_data)[23:27]) # All predictors except labs and HBI (199)
#essential_nodes_plot_noT(all_response_dags,4,names(all_response_dags_data)[23:27]) # All predictors except labs (120)
essential_nodes_plot_noT(all_response_dags,5,names(all_response_dags_data)[23:27]) # Imputed dataset (199 x 50)

strength.plot(all_response_dags[[6]], all_response_dags[[1]], shape='ellipse')
strength.plot(all_response_dags[[8]], all_response_dags[[3]], shape='ellipse')
strength.plot(all_response_dags[[10]], all_response_dags[[5]], shape='ellipse')


# BinAbdSurg GLM
# Fistula
BinAbdSurg.DAGS.glm.preds=c(6)
BinAbdSurg.DAGS.glm=glm(respData.original[-rcs_cases,1]~.,data=all_response_dags_data[,BinAbdSurg.DAGS.glm.preds],family=binomial())
summary(BinAbdSurg.DAGS.glm)

# BinAbscess GLM
# Fistula, Abscess w/i 1 year prior
BinAbscess.DAGS.preds=c(6, 12)
BinAbscess.DAGS.glm=glm(respData.original[-rcs_cases,5]~.,data=all_response_dags_data[,BinAbscess.DAGS.preds],family=binomial())
summary(BinAbscess.DAGS.glm)

# BinObstruction
# Stricture
BinObstruction.DAGS.preds=c(7)
BinObstruction.DAGS.glm=glm(respData.original[-rcs_cases,4]~.,data=all_response_dags_data[,BinObstruction.DAGS.preds],family=binomial())
summary(BinObstruction.DAGS.glm)


dag1_hm=order(all_response_dags[[1]][,3], decreasing=T)==
all_response_dags[[1]][dag1_hm[1:30],]

dag2_hm=order(all_response_dags[[2]][,3], decreasing=T)
all_response_dags[[2]][dag2_hm[1:30],]

dag3_hm=order(all_response_dags[[3]][,3], decreasing=T)
all_response_dags[[3]][dag3_hm[1:30],]

# Get strengths for arcs in BinAbdSurg minimal adjustment set
filter(all_response_dags[[1]], from=="Fistula." & to=="BinAbdSurg")
filter(all_response_dags[[2]], from=="Fistula." & to=="BinAbdSurg")
filter(all_response_dags[[3]], from=="Fistula." & to=="BinAbdSurg")

# Get strengths for arcs in BinAbscess minimal adjustment set
filter(all_response_dags[[1]], from=="Fistula." & to=="Intra.abd.Abscess.")
filter(all_response_dags[[2]], from=="Fistula." & to=="Intra.abd.Abscess.")
filter(all_response_dags[[3]], from=="Fistula." & to=="Intra.abd.Abscess.")
filter(all_response_dags[[1]], from=="Abscess.w.in.1yr.prior...Baseline." & to=="Intra.abd.Abscess.")
filter(all_response_dags[[2]], from=="Abscess.w.in.1yr.prior...Baseline." & to=="Intra.abd.Abscess.")
filter(all_response_dags[[3]], from=="Abscess.w.in.1yr.prior...Baseline." & to=="Intra.abd.Abscess.")
filter(all_response_dags[[1]], from=='Intra.abdominal.Abscess.history.' & to=='Abscess.w.in.1yr.prior...Baseline.')
filter(all_response_dags[[1]], from=='Intra.abdominal.Abscess.history.' & to=='Abscess.w.in.1yr.prior...Baseline.')
filter(all_response_dags[[1]], from=='Intra.abdominal.Abscess.history.' & to=='Abscess.w.in.1yr.prior...Baseline.')


# Get strengths for arcs to BinObstruction
filter(all_response_dags[[1]], from=='Stricture.' & to=='BinObstructions')
filter(all_response_dags[[2]], from=='Stricture.' & to=='BinObstructions')
filter(all_response_dags[[3]], from=='Stricture.' & to=='BinObstructions')


# DAGGITY Dag Code
dag {
  Abd_Surg [pos="-0.709,1.207"]
  Abscess [outcome,pos="-0.024,1.306"]
  Abscess_1Yr_Prior [pos="0.036,0.128"]
  Abscess_History [pos="-0.180,-0.633"]
  Fistula [exposure,pos="-0.729,0.412"]
  Inflammation [pos="-1.007,-0.713"]
  Obstruction [pos="-1.753,1.234"]
  Stricture [pos="-1.753,0.323"]
  Abscess_1Yr_Prior -> Abscess
  Abscess_History -> Abscess_1Yr_Prior
  Abscess_History -> Fistula
  Fistula -> Abd_Surg
  Fistula -> Abscess
  Inflammation -> Fistula
  Inflammation -> Stricture
  Stricture -> Obstruction
}

get_new_unitab=function(pred, resp){
  
  P=c()
  for(j in 1:ncol(pred)){
    P[j]=chisq.test(pred[,j],resp)$p.value
    if(P[j]<.05)
      print(paste(names(pred)[j],': ',P[j],sep=''))

  }
  
  return(P)
  
}



# BinAbdSurg DAGs ######
BinAbdSurg.DAGS=get_DAGS(predData1.stan,respData.original[,1],"BinAbdSurgeries",c(9,13:17),predData1.stan.imp,1)
essential_nodes_plot(BinAbdSurg.DAGS,1,"BinAbdSurgeries") # Complete cases
essential_nodes_plot(BinAbdSurg.DAGS,2,"BinAbdSurgeries") # Missing only HBI
essential_nodes_plot(BinAbdSurg.DAGS,3,"BinAbdSurgeries") # Missing HBI and labs
essential_nodes_plot(BinAbdSurg.DAGS,4,"BinAbdSurgeries") # Missing labs
essential_nodes_plot(BinAbdSurg.DAGS,5,"BinAbdSurgeries") # Imputed dataset

# Fistula, Steroid34

BinAbdSurg.DAGS.glm.cc=which(complete.cases(predData1.stan[,BinAbdSurg.DAGS.glm.preds]))
BinAbdSurg.DAGS.glm.prob=predict(BinAbdSurg.DAGS.glm,type="response")
BinAbdSurg.DAGS.glm.ROC=roc(respData.original[BinAbdSurg.DAGS.glm.cc,1]~BinAbdSurg.DAGS.glm.prob,direction="<")
plot.roc(BinAbdSurg.DAGS.glm.ROC, main="Abdominal Surgeries",
         print.auc=T,print.auc.cex=1.1,print.thres=T,print.thres.pch=17,print.thres.cex=1.1)

# BinOstomy DAGS ######
BinOstomy.DAGS=get_DAGS(predData1.stan,respData.original[,2],"BinOstomy",c(9,13:17),predData1.stan.imp,2)
lean_strength_plot(BinOstomy.DAGS,1)
lean_strength_plot(BinOstomy.DAGS,2)
lean_strength_plot(BinOstomy.DAGS,3)
lean_strength_plot(BinOstomy.DAGS,4)
lean_strength_plot(BinOstomy.DAGS,5)

# Variables taken from unitab (DAGs say no relation)
# CRP, ESR, Biologics46
BinOstomy.DAGS.glm=glm(respData.original[,2]~.,data=predData.1.stan[,c(13,14,26)],family=binomial())
summary(BinOstomy.DAGS.glm)

# BinHospitalization DAGS #####
BinHospitlization.DAGS=get_DAGS(predData1.stan,respData.original[,3],"BinHospitalization",c(9,13:17),predData1.stan.imp,3)
essential_nodes_plot(BinHospitlization.DAGS,1,"BinHospitalization")
essential_nodes_plot(BinHospitlization.DAGS,2,"BinHospitalization")
essential_nodes_plot(BinHospitlization.DAGS,3,"BinHospitalization")
essential_nodes_plot(BinHospitlization.DAGS,4,"BinHospitalization")
essential_nodes_plot(BinHospitlization.DAGS,5,"BinHospitalization")

BinHospitlization.DAGS.glm=glm(respData.original[,3]~.,data=predData1.stan[,c(6,9,18)],family=binomial())
summary(BinHospitlization.DAGS.glm)

BinHospitlization.DAGS.glm.cc=which(complete.cases(predData1.stan[,c(6,9,18)]))
BinHospitlization.DAGS.glm.prob=predict(BinHospitlization.DAGS.glm,type="response")
BinHospitlization.DAGS.glm.ROC=roc(respData.original[BinHospitlization.DAGS.glm.cc,3]~
                                     BinHospitlization.DAGS.glm.prob,
                                   direction="<")
plot(BinHospitlization.DAGS.glm.ROC)

summary(glm(respData.original[,3]~.,data=predData1.stan[,c(6,18,9,13:17)],family=binomial()))

BinHospitlization.DAGS.glm.imp=pool(with(predData1.stan.imp,glm(respData.original[,3]~.,
                                                                data=predData1.stan[,c(6,18,9)],
                                                                family=binomial())))
summary(BinHospitlization.DAGS.glm.imp)

# BinObstructions DAGS #####
BinObstruction.DAGS=get_DAGS(predData1.stan,respData.original[,4],"BinObstruction",c(9,13:17),predData1.stan.imp,4)
essential_nodes_plot(BinObstruction.DAGS,1,"BinObstruction")
essential_nodes_plot(BinObstruction.DAGS,2,"BinObstruction")
essential_nodes_plot(BinObstruction.DAGS,3,"BinObstruction")
essential_nodes_plot(BinObstruction.DAGS,4,"BinObstruction")
essential_nodes_plot(BinObstruction.DAGS,5,"BinObstruction")


# Stricture, Steroids34


# BinAbscess DAGS #####
BinAbscess.DAGS=get_DAGS(predData1.stan,respData.original[,5],"BinAbscess",c(9,13:17),predData1.stan.imp,5)
essential_nodes_plot(BinAbscess.DAGS,1,"BinAbscess")
essential_nodes_plot(BinAbscess.DAGS,2,"BinAbscess")
essential_nodes_plot(BinAbscess.DAGS,3,"BinAbscess")
essential_nodes_plot(BinAbscess.DAGS,4,"BinAbscess")
essential_nodes_plot(BinAbscess.DAGS,5,"BinAbscess")

# Fistula, Abscess w/i 1 year prior
BinAbscess.DAGS.glm.preds=c(6,12)
BinAbscess.DAGS.glm=glm(respData.original[,5]~.,data=predData1.stan[,BinAbscess.DAGS.glm.preds],family=binomial())
summary(BinAbscess.DAGS.glm)

BinAbscess.DAGS.glm.cc=which(complete.cases(predData1.stan[,BinAbscess.DAGS.glm.preds]))
BinAbscess.DAGS.glm.prob=predict(BinAbscess.DAGS.glm,type="response")
BinAbcess.DAGS.glm.roc=roc(respData.original[BinAbscess.DAGS.glm.cc,5]~BinAbscess.DAGS.glm.prob,direction="<")
plot.roc(BinAbcess.DAGS.glm.roc, main="Abscess",
         print.auc=T,print.auc.cex=1.1,print.thres=T,print.thres.pch=17,print.thres.cex=1.1)


# BinObstruction_or_Abscess #####
BinObstruction_or_Abscess_y=as.numeric(as.numeric(respData.original[,4])-1 | as.numeric(respData.original[,5])-1)
BinObstruction_or_Abscess.DAGS=get_DAGS(predData1,BinObstruction_or_Abscess_y,"BinObstruction_or_Abscess",c(9,13:17),predData1.stan.imp,31)

essential_nodes_plot(BinObstruction_or_Abscess.DAGS,1,"BinObstruction_or_Abscess")
essential_nodes_plot(BinObstruction_or_Abscess.DAGS,2,"BinObstruction_or_Abscess")
essential_nodes_plot(BinObstruction_or_Abscess.DAGS,3,"BinObstruction_or_Abscess")
essential_nodes_plot(BinObstruction_or_Abscess.DAGS,4,"BinObstruction_or_Abscess")
essential_nodes_plot(BinObstruction_or_Abscess.DAGS,5,"BinObstruction_or_Abscess")

BinObstruction_or_Abscess.DAGS.glm.preds=c(6,18)
BinObstruction_or_Abscess.DAGS.glm=glm(BinObstruction_or_Abscess_y~.,data=predData1[,BinObstruction_or_Abscess.DAGS.glm.preds],family=binomial())
summary(BinObstruction_or_Abscess.DAGS.glm)

# BinSurgery_or_Obstruction_or_Abscess #####
BinSurgery_or_Obstruction_or_Abscess_y=as.numeric(as.numeric(respData.original[,1])-1 | as.numeric(respData.original[,4])-1 | as.numeric(respData.original[,5])-1)
BinSurgery_or_Obstruction_or_Abscess.DAGS=get_DAGS(predData1,BinSurgery_or_Obstruction_or_Abscess_y,"BinSurgery_or_Obstruction_or_Abscess",c(9,13:17),predData1.stan.imp,31)

essential_nodes_plot(BinSurgery_or_Obstruction_or_Abscess.DAGS,1,"BinSurgery_or_Obstruction_or_Abscess")
essential_nodes_plot(BinSurgery_or_Obstruction_or_Abscess.DAGS,2,"BinSurgery_or_Obstruction_or_Abscess")
essential_nodes_plot(BinSurgery_or_Obstruction_or_Abscess.DAGS,3,"BinSurgery_or_Obstruction_or_Abscess")
essential_nodes_plot(BinSurgery_or_Obstruction_or_Abscess.DAGS,4,"BinSurgery_or_Obstruction_or_Abscess")
essential_nodes_plot(BinSurgery_or_Obstruction_or_Abscess.DAGS,5,"BinSurgery_or_Obstruction_or_Abscess")

BinSurgery_or_Obstruction_or_Abscess.DAGS.glm.preds=c(13,18)
BinSurgery_or_Obstruction_or_Abscess.DAGS.glm=glm(BinSurgery_or_Obstruction_or_Abscess_y~.,data=predData1[,BinSurgery_or_Obstruction_or_Abscess.DAGS.glm.preds],family=binomial())
summary(BinSurgery_or_Obstruction_or_Abscess.DAGS.glm)


BinAbdSurg.all_dags.preds=c(6,18)
BinAbdSurg.all_dags.glm=glm(all_response_dags_data[,34]~.,data=all_response_dags_data[,BinAbdSurg.all_dags.preds],family=binomial())
summary(BinAbdSurg.all_dags.glm)
plot(BinAbdSurg.all_dags.glm$residuals)

BinStricture.all_dags.preds=c(7,18)
BinObstruction.all_dags.glm=glm(all_response_dags_data[,37]~.,data=all_response_dags_data[,BinStricture.all_dags.preds],family=binomial())
summary(BinObstruction.all_dags.glm)
plot(BinObstruction.all_dags.glm$residuals)

BinAbscess.all_dags.preds=c(6,11)
BinAbscess.all_dags.glm=glm(all_response_dags_data[,38]~.,data=all_response_dags_data[,BinAbscess.all_dags.preds],family=binomial())
summary(BinAbscess.all_dags.glm)
plot(BinAbscess.all_dags.glm$residuals)

#### Num DAGS #############################################################################################################
###########################################################################################################################
###########################################################################################################################
###########################################################################################################################

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
NumHospitalization.DAGS.glm=glm(respData.original[,8]~.,data=predData1.stan[,c(9,18)],family=poisson())
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


# Other things ############################################################################################################
###########################################################################################################################
###########################################################################################################################
###########################################################################################################################

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


# BinAbdSurg multivariate #####
# Regular multivariate
abdSurg.bin.o.preds=sort(union(preds.include,
                                  getPredIndices(unitab1,"BinAbdSurg (1)")))
abdSurg.bin.o.glm=glm(respData.original[,1]~.,
                      data=predData1[,abdSurg.bin.o.preds],
                      family=binomial())
summary(abdSurg.bin.o.glm)
# Imputed multivariate
abdSurg.bin.o.preds.imp=union(abdSurg.bin.o.preds,
                              getPredIndices(unitab1.imp,"BinAbdSurg (1)"))

# Compare predictor sets
abdSurg.bin.o.preds
abdSurg.bin.o.preds.imp

abdSurg.bin.o.glm.imp=with(predData1.imp2,
                           glm(respData.original[,1]~.,
                               data=predData1[,abdSurg.bin.o.preds.imp],
                               family=binomial()))
summary(pool(abdSurg.bin.o.glm.imp))
# BinOstomy #####
ostomy.bin.o.preds=sort(union(preds.include,
                              getPredIndices(unitab1,"BinOstomy (2)")))
ostomy.bin.o.glm=glm(respData.original[,2]~.,
                     data=predData1[,ostomy.bin.o.preds],
                     family=binomial())
summary(ostomy.bin.o.glm)
# Imputed multivariate
ostomy.bin.o.preds.imp=sort(union(preds.include,
                              getPredIndices(unitab1.imp,"BinOstomy (2)")))

# Compare predictor sets
ostomy.bin.o.preds
ostomy.bin.o.preds.imp

ostomy.bin.o.glm.imp=with(predData1.imp2,
                          glm(respData.original[,2]~.,
                              data=predData1[,ostomy.bin.o.preds.imp],
                              family=binomial()))
# PROBLEM: to many variables included which predicts ostomy at univariate level
summary(pool(ostomy.bin.o.glm.imp))
# BinHospitalization #####
hospitalization.bin.o.preds=sort(union(preds.include,
                                          getPredIndices(unitab1,"BinHospitlizations (3)")))
hospitalization.bin.o.glm=glm(respData.original[,3]~.,
                              data=predData1[,hospitalization.bin.o.preds],
                              family=binomial())
summary(hospitalization.bin.o.glm)
# Imputed multivariate
hospitalization.bin.o.preds.imp=sort(union(preds.include,
                                       getPredIndices(unitab1.imp,"BinHospitlizations (3)")))

# Compare predictor sets
hospitalization.bin.o.preds
hospitalization.bin.o.preds.imp

hospitalization.bin.o.glm.imp=with(predData1.imp2,
                                   glm(respData.original[,3]~.,
                                       data=predData1[,hospitalization.bin.o.preds.imp],
                                       family=binomial()))
summary(pool(hospitalization.bin.o.glm.imp))
# BinObstruction #####
obstruction.bin.o.preds=sort(union(preds.include,
                                   getPredIndices(unitab1,"BinObstructions (4)")))
obstruction.bin.o.glm=glm(respData.original[,4]~.,
                          data=predData1[,obstruction.bin.o.preds],
                          family=binomial())
summary(obstruction.bin.o.glm)
# Imputed multivariate
obstruction.bin.o.preds.imp=sort(union(preds.include,
                                           getPredIndices(unitab1.imp,"BinObstructions (4)")))

# Compare predictor sets
obstruction.bin.o.preds
obstruction.bin.o.preds.imp

obstruction.bin.o.glm.imp=with(predData1.imp2,
                                   glm(respData.original[,4]~.,
                                       data=predData1[,obstruction.bin.o.preds.imp],
                                       family=binomial()))
summary(pool(obstruction.bin.o.glm.imp))
#PROBLEM: nearly all imputed variables nonsignificant for obstructions at univariate level

# BinAbscess #####
abscess.bin.o.preds=sort(union(preds.include,
                                  getPredIndices(unitab1,"Intra.abd.Abscess. (5)")))
abscess.bin.o.glm=glm(respData.original[,5]~.,
                      data=predData1[,abscess.bin.o.preds],
                      family=binomial())
summary(abscess.bin.o.glm)
# Imputated multivariate
abscess.bin.o.preds.imp=sort(union(preds.include,
                              getPredIndices(unitab1.imp,"Intra.abd.Abscess. (5)")))

# Compare predictor sets
abscess.bin.o.preds
abscess.bin.o.preds.imp

abscess.bin.o.glm.imp=with(predData1.imp2,
                           glm(respData.original[,5]~.,
                               data=predData1[,abscess.bin.o.preds.imp],
                               family=binomial()))
summary(pool(abscess.bin.o.glm.imp))

# #AbdSurgeries #####
abdSurg.num.o.preds=sort(union(preds.include,
                               getPredIndices(unitab1,"#AbdSurg (6)")))
abdSurg.num.o.glm=glm(respData.original[,6]~.,
                      data=predData1[,abdSurg.num.o.preds],
                      family=poisson())
summary(abdSurg.num.o.glm)
# Imputed multivariate
abdSurg.num.o.preds.imp=sort(union(preds.include,
                                   getPredIndices(unitab1,"#AbdSurg (6)")))

# Compare predictor sets
abdSurg.num.o.preds
abdSurg.num.o.preds.imp

abdSurg.num.o.glm.imp=with(predData1.imp2,
                           glm(respData.original[,6]~.,
                               data=predData1[,abdSurg.num.o.preds.imp],
                               family=poisson()))
summary(pool(abdSurg.num.o.glm.imp))

# #Ostomies #####
ostomy.num.o.preds=sort(union(preds.include,
                              getPredIndices(unitab1,"#Ostomy (7)")))
ostomy.num.o.glm=glm(respData.original[,7]~.,
                     data=predData1[,ostomy.num.o.preds],
                     family=poisson())
summary(ostomy.num.o.glm)

# Imputed multivariate
ostomy.num.o.preds.imp=sort(union(ostomy.num.o.preds,
                                  getPredIndices(unitab1.imp,"#Ostomy (7)")))

# Compare predictor sets
ostomy.num.o.preds
ostomy.num.o.preds.imp

ostomy.num.o.glm.imp=with(predData1.imp2,
                          glm(respData.original[,7]~.,
                              data=predData1[,ostomy.num.o.preds.imp],
                              family=poisson()))
# PROBLEM: perfect separation on imputed ostomy glm
summary(pool(ostomy.num.o.glm.imp))

# #Hospitalizations #####
hospitalization.num.o.preds=sort(union(preds.include,
                                       getPredIndices(unitab1,"#Hospitalization (8)")))
hospitalization.num.o.glm=glm(respData.original[,8]~.,
                              data=predData1[,hospitalization.num.o.preds],
                              family=poisson())
summary(hospitalization.num.o.glm)
# Imputed multivariate
hospitalization.num.o.preds.imp=sort(union(hospitalization.num.o.preds,
                                           getPredIndices(unitab1.imp,"#Hospitalization (8)")))

# Compare predictor sets
hospitalization.num.o.preds
hospitalization.num.o.preds.imp

hospitalization.num.o.glm.imp=with(predData1.imp2,
                                   glm(respData.original[,8]~.,
                                       data=predData1[,hospitalization.num.o.preds.imp],
                                       family=poisson()))
summary(pool(hospitalization.num.o.glm.imp))

# #Obstructions #####
obstruction.num.o.preds=sort(union(preds.include,
                                   getPredIndices(unitab1,"#Obstructions (9)")))
obstruction.num.o.glm=glm(respData.original[,9]~.,
                          data=predData1[,obstruction.num.o.preds],
                          family=poisson())
summary(obstruction.num.o.glm)
# Imputed multivariate
obstruction.num.o.preds.imp=sort(union(obstruction.num.o.preds,
                                       getPredIndices(unitab1.imp,"#Obstructions (9)")))

# Compare predictor sets
obstruction.num.o.preds
obstruction.num.o.preds.imp

obstruction.num.o.glm.imp=with(predData1.imp2,
                               glm(respData.original[,9]~.,
                                   data=predData1[,obstruction.num.o.preds.imp],
                                   family=poisson()))
summary(pool(obstruction.num.o.glm.imp))

# Combined responses unitab #####
unitab2=getUniModels(predData1,respData.combined,respData.combined.types,.05)[[1]]
View(unitab2)

unitab2.imp=getUniModels.imp(predData1,predData1.imp2,respData.combined,respData.combined.types,.05)[[1]]

preds.include=c(predData1.radPreds, predData1.labPreds)

# BinHospitalizations combined #####
hospitalizations.bin.c.preds=sort(union(preds.include,
                                        getPredIndices(unitab2,"BinHospitalizations (1)")))
hospitalizations.bin.c.glm=glm(respData.combined[,1]~.,
                               data=predData1[,hospitalizations.bin.c.preds],
                               family=binomial())
summary(hospitalizations.bin.c.glm)
# Imputed multivariate
hospitalizations.bin.c.preds.imp=sort(union(hospitalizations.bin.c.preds,
                                            getPredIndices(unitab2.imp,"BinHospitalizations (1)")))

# Compare predictor sets
hospitalizations.bin.c.preds
hospitalizations.bin.c.preds.imp

hospitalizations.bin.c.glm.imp=with(predData1.imp2,
                                    glm(respData.combined[,1]~.,
                                        data=predData1[,hospitalizations.bin.c.preds.imp],
                                        family=binomial()))
summary(pool(hospitalizations.bin.c.glm.imp))
# BinSurgeries combined #####
surgeries.bin.c.preds=sort(union(preds.include,
                                 getPredIndices(unitab2,"BinSurgeries (2)")))
surgeries.bin.c.glm=glm(respData.combined[,2]~.,
                        data=predData1[,surgeries.bin.c.preds],
                        family=binomial())
summary(surgeries.bin.c.glm)
# Imputed multivariate
surgeries.bin.c.preds.imp=sort(union(surgeries.bin.c.preds,
                                     getPredIndices(unitab2,"BinSurgeries (2)")))

# Compare predictor sets
surgeries.bin.c.preds
surgeries.bin.c.preds.imp

surgeries.bin.c.glm.imp=with(predData1.imp2,
                             glm(respData.combined[,2]~.,
                                 data=predData1[,surgeries.bin.c.preds.imp],
                                 family=binomial()))
summary(pool(surgeries.bin.c.glm.imp))
# #Hospitalizations combined #####
hospitalization.num.c.preds=sort(union(preds.include,
                                      getPredIndices(unitab2,"#Hospitalizations (3)")))
hospitalization.num.c.glm=glm(respData.combined[,3]~.,
                              data=predData1[,hospitalization.num.c.preds],
                              family=poisson())
summary(hospitalization.num.c.glm)
# Imputed multivariate
hospitalization.num.c.preds.imp=sort(union(hospitalization.num.c.preds,
                                           getPredIndices(unitab2.imp,"#Hospitalizations (3)")))

# Compare predictor sets
hospitalization.num.c.preds
hospitalization.num.c.preds.imp

hospitalization.num.c.glm.imp=with(predData1.imp2,
                                   glm(respData.combined[,3]~.,
                                       data=predData1[,hospitalization.num.c.preds.imp],
                                       family=poisson()))
summary(pool(hospitalization.num.c.glm.imp))

# #Surgeries combined #####
surgeries.num.c.preds=sort(union(preds.include,
                                 getPredIndices(unitab2,"#Surgeries (4)")))
surgeries.num.c.glm=glm(respData.combined[,4]~.,
                        data=predData1[,surgeries.num.c.preds],
                        family=poisson())
summary(surgeries.num.c.glm)
# Imputed multivariate
surgeries.num.c.preds.imp=sort(union(surgeries.num.c.preds,
                                     getPredIndices(unitab2.imp,"#Surgeries (4)")))

# Compare predictor sets
surgeries.num.c.preds
surgeries.num.c.preds.imp

surgeries.num.c.glm.imp=with(predData1.imp2,
                             glm(respData.combined[,4]~.,
                                 data=predData1[,surgeries.num.c.preds.imp],
                                 family=poisson()))
summary(pool(surgeries.num.c.glm.imp))
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
