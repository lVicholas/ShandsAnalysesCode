library(readxl)
library(writexl)
library(MASS)

# Data not included due to sensitivity
setwd("C:/local/R/IBD_database_study")
ibd_endounters=data.frame(read_xlsx("Encounters in IBD With Race and Ethnicity.xlsx"))

response_names=c("Ambulatory_total","Ed_total","Ed_to_inpatient_total","Inpatient_total",
                 "Ambulatory_+_anti-TNF","Ed_+_anti-TNF","Ed_to_inpatient_+_anti-TNF","Inpatient_+_anti-TNF",
                 "Ambulatory_-_anti-TNF","Ed_-_anti-TNF","Ed_to_inpatient_-_anti-TNF","Inpatient_-_anti-TNF")

# Total # of people in database for each race
# 1st 4 elements for Caucasians, Hispanics, African-Americans, Asian-Americans
# 2nd 4 elements for # of each race taking an anti-TNF
# 3rd 4 elements for # of each race non taking anti-TNF
race_totals=c(32138,6690,4835,370,3855,657,489,31)
race_totals=t(data.frame(c(race_totals,race_totals[1:4]-race_totals[5:8])))
names(race_totals)=c("Caucassians_total","Hispanics_total","African_Americans_total","Asian_Americans_total",
                     "Caucassians_+_anti-TNF","Hispanics_+_anti-TNF","African_Americans_+_anti-TNF","Asian_Americans_+_anti-TNF",
                     "Caucassians_-_anti-TNF","Hispanics_-_anti-TNF","African_Americans_-_anti-TNF","Asian_Americans_-_anti-TNF")

# Numbers of people who did not use a response (e.g. ambulatory visit)
# Also separated by anti-TNF status
race_zero_counts=list()
race_zero_counts[["Ambulatory_total"]]=c(2001,527,272,33)
race_zero_counts[["ED_total"]]=c(12438,1525,985,167)
race_zero_counts[["ED_to_Inpatient_total"]]=c(18522,3383,2175,243)
race_zero_counts[["Inpatient_total"]]=c(18253,3628,2524,219)

race_zero_counts[["Ambulatory_+_anti-TNF"]]=c(39,8,0,0)
race_zero_counts[["ED_+_anti-TNF"]]=c(1985,161,138,23)
race_zero_counts[["ED_to_inpatient_+_anti-TNF"]]=c(2583,403,235,27)
race_zero_counts[["Inpatient_+_anti-TNF"]]=c(2256,196,250,17)

race_zero_counts[["Ambulatory_-_anti-TNF"]]=race_zero_counts[[1]]-race_zero_counts[[5]]
race_zero_counts[["ED_-_anti-TNF"]]=race_zero_counts[[2]]-race_zero_counts[[6]]
race_zero_counts[["ED_to_inpatient_-_anti-TNF"]]=race_zero_counts[[3]]-race_zero_counts[[7]]
race_zero_counts[["Inpatient_-_anti-TNF"]]=race_zero_counts[[4]]-race_zero_counts[[8]]

# Get tables for comparisons across race #####
race_t_test=function(r1,r2,response,TNF){
  
  # For races r1 and r2, perform Welch t-test for difference in proportions of
  # people in each race responsible for at least 1 occurence of response
  
  n1=as.numeric(race_totals[4*TNF+r1])
  c1=n1-race_zero_counts[[4*TNF+response]][r1]
  n2=as.numeric(race_totals[4*TNF+r2])
  c2=n2-race_zero_counts[[4*TNF+response]][r2]
  
  x=rep(c(1,0),c(c1,n1-c1))
  y=rep(c(1,0),c(c2,n2-c2))
  return(t.test(x,y,var.equal=F))
  
}
get_table_race_t_tests=function(r1,r2){
  
  # Gets table of Welch t-tests for differnce in proportions between races of number of people
  # responsible for at least 1 occurence of response
  # E.g., test whether Caucasians are more likely to have at least 1 ambulatory visit than Hispanics
  # Also tests based on anti-TNF status
  
  rows=response_names
  N1=N2=c()
  C1=C2=c()
  r1.p=r2.p=c()
  T.vec=P=c()
  CI=c()
  
  for(i in 1:4){
    for(j in 0:2){
      
      N1[4*j+i]=as.numeric(race_totals[4*j+r1])
      N2[4*j+i]=as.numeric(race_totals[4*j+r2])
      C1[4*j+i]=N1[4*j+i]-race_zero_counts[[4*j+i]][r1]
      C2[4*j+i]=N2[4*j+i]-race_zero_counts[[4*j+i]][r2]
      
      x=rep(c(1,0),c(C1[4*j+i],N1[4*j+i]-C1[4*j+i]))
      y=rep(c(1,0),c(C2[4*j+i],N2[4*j+i]-C2[4*j+i]))
      test=t.test(x,y,var.equal=F)
      
      r1.p[4*j+i]=round(test$estimate[1],4)
      r2.p[4*j+i]=round(test$estimate[2],4)
      
      T.vec[4*j+i]=round(test$statistic,4)
      P[4*j+i]=round(test$p.value,4)
      
      CI[4*j+i]=paste("(",round(test$conf.int[1],3),", ",round(test$conf.int[2],3),")",sep="")
      
    }
  }
  
  D=data.frame(rows,N1,C1,N2,C2,r1.p,r2.p,T.vec,P,CI)
  names(D)=c("Response",
             paste(names(race_totals)[r1],"n"), paste("# of successes",names(race_totals)[r1]),
             paste(names(race_totals)[r2],"n"), paste("# of successes",names(race_totals)[r2]),
             paste(names(race_totals)[r1],"proportion successes"), 
             paste(names(race_totals)[r2],"proportion successes"),
             "Welch test t-stat","Welch test p-val","95% Conf.Int")
  return(D)
  
}

# Writing tables to files
caucasian_vs_other_races_test=list(Caucasians_VS_Hispanics=get_table_race_t_tests(1,2),
                                   Caucasians_VS_African_Americans=get_table_race_t_tests(1,3),
                                   Caucasians_VS_Asian_Americans=get_table_race_t_tests(1,4))
write_xlsx(caucasian_vs_other_races_test,"Caucasians_vs_Other_Races.xlsx")

hispanics_vs_other_races=list(Hispanics_VS_Caucasians=get_table_race_t_tests(2,1),
                              Hispanics_VS_African_Americans=get_table_race_t_tests(2,3),
                              Hispanics_VS_Asian_Americans=get_table_race_t_tests(2,4))
write_xlsx(hispanics_vs_other_races,"Hispanics_vs_Other_Races.xlsx")

# Get tables for comparisons across TNF status #####

get_TNF_comparisons=function(race){
  
  # Compares proportions of people with at least 1 use of response for each race with and without anti-TNF
  # E.g., whether Hispanics on anti-TNF are more likely to have at least 1 ED visit than Hispanics not on anti-TNF
  
  resp_names=c("Ambulatory","ED","ED_to_inpatient","Inpatient")
  N1=N2=c()
  C1=C2=c()
  TNF.p=non_TNF.p=c()
  T.vec=P=c()
  CI=c()
  
  for(i in 1:4){
    
    N1[i]=race_totals[race+4]
    N2[i]=race_totals[race+8]
    
    C1[i]=N1[i]-race_zero_counts[[i+4]][race]
    C2[i]=N2[i]-race_zero_counts[[i+8]][race]
    
    x=rep(c(1,0),c(C1[i],N1[i]-C1[i]))
    y=rep(c(1,0),c(C2[i],N2[i]-C2[i]))
    test=t.test(x,y,var.equal=F)
    
    TNF.p[i]=round(test$estimate[1],4)
    non_TNF.p[i]=round(test$estimate[2],4)
    
    T.vec[i]=round(test$statistic,4)
    P[i]=round(test$p.value,4)
    
    CI[i]=paste("(",round(test$conf.int[1],3),", ",round(test$conf.int[2],3),")",sep="")
    
  }
  
  D=data.frame(resp_names,N1,C1,N2,C2,TNF.p,non_TNF.p,T.vec,P,CI)
  names(D)=c("Response",
             "n for anti-TNF", "# of successes for anti-TNF",
             "n for non-anti-TNF","# of successes for non-anti-TNF",
             "anti-TNF proportion of successes","non-anti-TNF proportion of successes",
             "Welch t-statistic","Welch p-value","95% Conf.int")
  return(D)
  
}

TNF_comparisons=list(Caucasians=get_TNF_comparisons(1),
                     Hispanics=get_TNF_comparisons(2),
                     African_Americans=get_TNF_comparisons(3),
                     Asian_Americans=get_TNF_comparisons(4))
write_xlsx(TNF_comparisons,"TNF_comparisons.xlsx")

###### Non-IBD encounters ##################################################################################################
############################################################################################################################

non_ibd_race_totals=c(11308102,4015296,3056464,238151)
non_ibd_race_totals.ambulance=c(8253663,2825098,2294987,163346)
non_ibd_race_totals.ed_visit=c(5585354,2246372,1865123,88837)

non_ibd_encounters.ambulance.race_zeroes=c(3056437,1190732,761751,73835)
non_ibd_encounters.ed_visit.race_zeroes=c(5735179,1770452,1192327,149483)

ambulance_and_CRP_2_zeroes=c(497,155,64,10)
ambulance_and_hgb_12_zeroes=c(965,251,153,20)

my_t_test=function(total1,zero1,total2,zero2){
  
  x=rep(c(0,1),times=c(zero1,total1-zero1))
  y=rep(c(0,1),times=c(zero2,total2-zero2))
  return(t.test(x,y,var.equal=F))
  
}
get_race_t_tests=function(totals,zeroes,r1,r2){
  
  return(
    my_t_test(totals[r1],zeroes[r1],totals[r2],zeroes[r2])
  )
  
}

# Non-IBD Encounters
# Ambulance
white_vs_hispanic_non_ibd.ambulance=get_race_t_tests(non_ibd_race_totals.ambulance,non_ibd_encounters.ambulance.race_zeroes,1,2)
white_vs_african_american_non_ibd.ambulance=get_race_t_tests(non_ibd_race_totals.ambulance,non_ibd_encounters.ambulance.race_zeroes,1,3)
white_vs_asian_american_non_ibd.ambulance=get_race_t_tests(non_ibd_race_totals.ambulance,non_ibd_encounters.ambulance.race_zeroes,1,4)

hispanic_vs_african_american_non_ibd.ambulance=get_race_t_tests(non_ibd_race_totals.ambulance,non_ibd_encounters.ambulance.race_zeroes,2,3)
hispanic_vs_asian_american_non_ibd.ambulance=get_race_t_tests(non_ibd_race_totals.ambulance,non_ibd_encounters.ambulance.race_zeroes,2,4)

african_american_vs_asian_american_non_ibd.ambulance=get_race_t_tests(non_ibd_race_totals.ambulance,non_ibd_encounters.ambulance.race_zeroes,3,4)

white_vs_hispanic_non_ibd.ambulance
white_vs_african_american_non_ibd.ambulance
white_vs_asian_american_non_ibd.ambulance

hispanic_vs_african_american_non_ibd.ambulance
hispanic_vs_asian_american_non_ibd.ambulance

african_american_vs_asian_american_non_ibd.ambulance

# ED vistic
white_vs_hispanic_non_ibd.ed=get_race_t_tests(non_ibd_race_totals.ed_visit,non_ibd_encounters.ed_visit.race_zeroes,1,2)
white_vs_african_american_non_ibd.ed=get_race_t_tests(non_ibd_race_totals.ed_visit,non_ibd_encounters.ed_visit.race_zeroes,1,3)
white_vs_asian_american_non_ibd.ed=get_race_t_tests(non_ibd_race_totals.ed_visit,non_ibd_encounters.ed_visit.race_zeroes,1,4)

hispanic_vs_african_american_non_ibd.ed=get_race_t_tests(non_ibd_race_totals.ed_visit,non_ibd_encounters.ed_visit.race_zeroes,2,3)
hispanic_vs_asian_american_non_ibd.ed=get_race_t_tests(non_ibd_race_totals.ed_visit,non_ibd_encounters.ed_visit.race_zeroes,2,4)

african_american_vs_asian_american_non_ibd.ed=get_race_t_tests(non_ibd_race_totals.ed_visit,non_ibd_encounters.ed_visit.race_zeroes,3,4)

# Ambulance + CRP > 2
white_vs_hispanic.ambulance_crp=get_race_t_tests(race_totals,ambulance_and_CRP_2_zeroes,1,2)
white_vs_african_american.ambulance_crp=get_race_t_tests(race_totals,ambulance_and_CRP_2_zeroes,1,3)
white_vs_asian_american.ambulance_crp=get_race_t_tests(race_totals,ambulance_and_CRP_2_zeroes,1,4)

hispanic_vs_african_american.ambulance_crp=get_race_t_tests(race_totals,ambulance_and_CRP_2_zeroes,2,3)
hispanic_vs_asian_american.ambulance_crp=get_race_t_tests(race_totals,ambulance_and_CRP_2_zeroes,2,4)

african_american_vs_asian_american.ambulance_crp=get_race_t_tests(race_totals,ambulance_and_CRP_2_zeroes,3,4)

white_vs_hispanic.ambulance_crp
white_vs_african_american.ambulance_crp
white_vs_asian_american.ambulance_crp

hispanic_vs_african_american.ambulance_crp
hispanic_vs_asian_american.ambulance_crp

african_american_vs_asian_american.ambulance_crp

# Ambulance + HGB < 12
white_vs_hispanic.ambulance_hgb=get_race_t_tests(race_totals,ambulance_and_hgb_12_zeroes,1,2)
white_vs_african_american.ambulance_hgb=get_race_t_tests(race_totals,ambulance_and_hgb_12_zeroes,1,3)
white_vs_asian_american.ambulance_hgb=get_race_t_tests(race_totals,ambulance_and_hgb_12_zeroes,1,4)

hispanic_vs_african_american.ambulance_hgb=get_race_t_tests(race_totals,ambulance_and_hgb_12_zeroes,2,3)
hispanic_vs_asian_american.ambulance_hgb=get_race_t_tests(race_totals,ambulance_and_hgb_12_zeroes,2,4)

african_american_vs_asian_american.ambulance_hgb=get_race_t_tests(race_totals,ambulance_and_hgb_12_zeroes,3,4)

white_vs_hispanic.ambulance_hgb
white_vs_african_american.ambulance_hgb
white_vs_asian_american.ambulance_hgb

hispanic_vs_african_american.ambulance_hgb
hispanic_vs_asian_american.ambulance_hgb

african_american_vs_asian_american.ambulance_hgb
