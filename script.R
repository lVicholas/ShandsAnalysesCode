library(readxl)
library(writexl)
library(MASS)

setwd("C:/local/R/ThomasStudy")
data=data.frame(read_xlsx("Encounters in IBD With Race and Ethnicity.xlsx"))

# Extract each of 8 tables
ambulatory1=data[3:10,3:6]
ambulatory2=data[3:10,10:13]
ed1=data[13:20,3:6]
ed2=data[13:20,10:13]
ed_to_inpatient1=data[23:30,3:6]
ed_to_inpatient2=data[23:30,10:13]
inpatient1=data[33:40,3:6]
inpatient2=data[33:40,10:13]

table_list=list(ambulatory1,ambulatory2,ed1,ed2,ed_to_inpatient1,ed_to_inpatient2,inpatient1,inpatient2)

appendTotals=function(table){
  nt=table
  nt=rbind(nt,colSums(nt))
  nt=cbind(nt,rowSums(nt))
  return(data.frame(nt))
}

# For each of 8 tables, ensure entries are numeric, add row/col sums, give proper row/col names
for(i in 1:length(table_list)){
  
  tab=data.frame(table_list[[i]])
  tab=apply(tab,2,as.numeric)
  tab=appendTotals(tab)
  
  if(i %% 2 == 1)
    names(tab)=c("IBD_in_Caucasians_(NONHISPANIC_+_WHITE)","IBD_in_Hispanics",
                             "IBD_in_African_Americans","IBD_in_Asian_Americans",
                             "Row_Totals")
  else names(tab)=c("IBD_+anti-TNF_in_Caucasians_(NON-HISPANIC_+_WHITE)","IBD_+_anti-TNF_in_Hispanics",
                                "IBD_+anti-TNF_in_African_Americans","IBD_+_anti-TNF_in_Asian_Americans",
                                "Row_Totals")
  
  row.names(tab)=c("0-9yo","10-17yo","18-34yo","35-44yo","45-54yo","55-64yo","65-74yo","75-84yo",
                               "Column_Totals")
  
  table_list[[i]]=tab #/tab[nrow(tab),ncol(tab)]
  
}

# Most counts (by race) suffer from overdispersion
cbind(apply(table_list[[1]],2,sd),apply(table_list[[1]],2,mean))
cbind(apply(table_list[[2]],2,sd),apply(table_list[[2]],2,mean))
cbind(apply(table_list[[3]],2,sd),apply(table_list[[3]],2,mean))
cbind(apply(table_list[[4]],2,sd),apply(table_list[[4]],2,mean))
cbind(apply(table_list[[5]],2,sd),apply(table_list[[5]],2,mean))
cbind(apply(table_list[[6]],2,sd),apply(table_list[[6]],2,mean))
cbind(apply(table_list[[7]],2,sd),apply(table_list[[7]],2,mean))
cbind(apply(table_list[[8]],2,sd),apply(table_list[[8]],2,mean))

# Most counts (by age group) suffer from overdispersion
cbind(apply(table_list[[1]],1,sd),apply(table_list[[1]],1,mean))
cbind(apply(table_list[[2]],1,sd),apply(table_list[[2]],1,mean))
cbind(apply(table_list[[3]],1,sd),apply(table_list[[3]],1,mean))
cbind(apply(table_list[[4]],1,sd),apply(table_list[[4]],1,mean))
cbind(apply(table_list[[5]],1,sd),apply(table_list[[5]],1,mean))
cbind(apply(table_list[[6]],1,sd),apply(table_list[[6]],1,mean))
cbind(apply(table_list[[7]],1,sd),apply(table_list[[7]],1,mean))
cbind(apply(table_list[[8]],1,sd),apply(table_list[[8]],1,mean))

melt_by_race=function(tab){
  
  return=data.frame(matrix(NA,nrow=nrow(tab)*2),ncol=2)
  names(return)=c("Race","Count")
  
  return[1:(nrow(return)/2),1:2]=cbind(rep(names(tab)[1],nrow(tab)),tab[1:nrow(tab),1])
  return[(1+nrow(return)/2):nrow(return),1:2]=cbind(rep(names(tab)[2],nrow(tab)),tab[1:nrow(tab),2])
  
  return=data.frame(return)
  return[,2]=as.numeric(return[,2])
  
  return(data.frame(return))
  
}

# Since overdispersion is present, fit negative-binomial models to race variables
get_race_mods=function(table){
  mods=list()
  for(j in 2:(ncol(table)-1)){
    data=melt_by_race(table[1:8,c(1,j)])
    mods[[j-1]]=glm.nb(Count~as.factor(Race),data=data)
  }
  return(mods)
}

# Get models for each response comparing Causassians to other minorities
list_of_list_of_race_mods=list()
for(i in 1:length(table_list))
  list_of_list_of_race_mods[[i]]=get_race_mods(table_list[[i]])

print_list_mods=function(mods){
  for(i in 1:length(mods))
    print(summary(mods[[i]]))
}
write_mod_list_to_excel=function(mods, workbook_name){
  sheets=list()
  for(i in 1:length(mods)){
    sheets[[paste("sheet",i,"Name",sep="")]]=data.frame(attr(coef(summary(mods[[i]])),"dimnames")[[1]],
                                                        coef(summary(mods[[i]])))
  }
  write_xlsx(sheets,workbook_name)
}

# Print and Write_to_Excel block #####
print_list_mods(list_of_list_of_race_mods[[1]])
print_list_mods(list_of_list_of_race_mods[[2]])
print_list_mods(list_of_list_of_race_mods[[3]])
print_list_mods(list_of_list_of_race_mods[[4]])
print_list_mods(list_of_list_of_race_mods[[5]])
print_list_mods(list_of_list_of_race_mods[[6]])
print_list_mods(list_of_list_of_race_mods[[7]])
print_list_mods(list_of_list_of_race_mods[[8]])

write_mod_list_to_excel(list_of_list_of_race_mods[[1]], "Ambulatory1.xlsx")
write_mod_list_to_excel(list_of_list_of_race_mods[[2]], "Ambulatory2.xlsx")
write_mod_list_to_excel(list_of_list_of_race_mods[[3]], "Ed1.xlsx")
write_mod_list_to_excel(list_of_list_of_race_mods[[4]], "Ed2.xlsx")
write_mod_list_to_excel(list_of_list_of_race_mods[[5]], "Ed_to_inpatient1.xlsx")
write_mod_list_to_excel(list_of_list_of_race_mods[[6]], "Ed_to_inpatient2.xlsx")
write_mod_list_to_excel(list_of_list_of_race_mods[[7]], "Inpatient1.xlsx")
write_mod_list_to_excel(list_of_list_of_race_mods[[8]], "Inpatient2.xlsx")

# Models by age group #####