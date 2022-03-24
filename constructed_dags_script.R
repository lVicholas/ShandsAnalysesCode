# Constructed DAGS

library(readxl)
library(writexl)
source = data.frame(read_excel('C:/_local/fistula_study/RadData2.xlsx'))

df = source

# Data cleaning #####

# Ensure numeric columns are numeric
df_numVars=c(2,6,12,24,26,28,30,32,57:60,70,71,74,78,81,85,86)
for(j in df_numVars)
   df[,j]=as.numeric(df[,j])

# Ensure categorical columns are categorical
df_catvars=c(3:5,8:11,14:16,19:22,34:55,66,67,69,75,83,87:99)
for(j in df_catvars)
   df[,j]=factor(df[,j])


y_vars = c(62,66,68,69,71)
junk = df[,y_vars]

y = data.frame(1:199)

y$abdSurg=factor(as.numeric(junk[,1]!=0))
y$ostomy=factor(as.numeric(junk[,2]!=0))
y$hospitlizations=factor(as.numeric(junk[,3]!=0))
y$obstruction=factor(as.numeric(junk[,4]!=0))
y$abscess=factor(junk[,5])

y=y[,-1]

x_vars = c(2:6,8:10,12,14,19,20,24,26,28,30,32,
           34:36,38,40:44,46,48:51,53,55)
X = df[,x_vars]

levels(X[,3]) = c('White', 'Black', 'Other')
levels(X[,4]) = c('Non-smoker', 'Former smoker', 'Smoker')

rm(list=c('junk'))

# End data cleaning #####

# DAG Code #####
'
dag {
Fistula [exposure,pos="-0.525,-0.539"]
abdominal_pain [pos="-0.877,-1.051"]
abscess [outcome,pos="-0.318,0.093"]
abscess_history [pos="-0.353,-1.099"]
abscess_wi_1yr [pos="-0.535,-1.000"]
age [pos="-1.353,-1.189"]
albumin [pos="-1.053,-0.523"]
biologics1 [pos="-1.310,-0.326"]
biologics2 [pos="-1.306,-0.149"]
biologics3 [pos="-1.310,0.088"]
bmi [pos="-1.090,-1.467"]
crp [pos="-1.335,-0.515"]
esr [pos="-1.240,-0.520"]
gender [pos="-1.233,-1.394"]
hbi [pos="-0.840,-0.598"]
hgb [pos="-1.147,-0.515"]
immunomodulators1 [pos="-0.983,-0.318"]
immunomodulators2 [pos="-0.983,-0.152"]
immunomodulators3 [pos="-0.983,0.082"]
inflammation_read [pos="-0.782,-0.905"]
jack_inhibitor1 [pos="-1.177,-0.329"]
jack_inhibitor2 [pos="-1.171,-0.152"]
jack_inhibitor3 [pos="-1.174,0.085"]
mesalamine1 [pos="-0.808,-0.332"]
mesalamine2 [pos="-0.798,-0.152"]
mesalamine3 [pos="-0.801,0.074"]
race [pos="-0.699,-1.514"]
smoker [pos="-0.961,-1.476"]
steroid1 [pos="-0.689,-0.343"]
steroid2 [pos="-0.675,-0.149"]
steroid3 [pos="-0.678,0.072"]
stricture_read [pos="-1.063,-0.924"]
wbc [pos="-0.948,-0.523"]
Fistula -> biologics2
Fistula -> biologics3
Fistula -> hbi
Fistula -> immunomodulators2
Fistula -> immunomodulators3
Fistula -> jack_inhibitor2
Fistula -> jack_inhibitor3
Fistula -> mesalamine2
Fistula -> mesalamine3
Fistula -> steroid2
Fistula -> steroid3
abdominal_pain -> Fistula
abdominal_pain -> biologics2
abdominal_pain -> biologics3
abdominal_pain -> hbi
abdominal_pain -> immunomodulators2
abdominal_pain -> immunomodulators3
abdominal_pain -> jack_inhibitor2
abdominal_pain -> jack_inhibitor3
abdominal_pain -> mesalamine2
abdominal_pain -> mesalamine3
abdominal_pain -> steroid2
abdominal_pain -> steroid3
abscess_history -> Fistula
abscess_history -> abscess
abscess_history -> abscess_wi_1yr
abscess_history -> hbi
abscess_wi_1yr -> Fistula
abscess_wi_1yr -> abscess
abscess_wi_1yr -> hbi
abscess_wi_1yr -> mesalamine2
abscess_wi_1yr -> steroid2
abscess_wi_1yr -> steroid3
age -> Fistula
age -> abdominal_pain
age -> abscess_history
age -> abscess_wi_1yr
age -> albumin
age -> biologics1
age -> biologics2
age -> biologics3
age -> bmi
age -> crp
age -> esr
age -> hbi
age -> hgb
age -> immunomodulators1
age -> immunomodulators2
age -> immunomodulators3
age -> jack_inhibitor1
age -> jack_inhibitor2
age -> jack_inhibitor3
age -> mesalamine1
age -> mesalamine2
age -> mesalamine3
age -> steroid1
age -> steroid2
age -> steroid3
age -> stricture_read
age -> wbc
albumin -> abdominal_pain
albumin -> biologics2
albumin -> biologics3
albumin -> immunomodulators2
albumin -> immunomodulators3
albumin -> inflammation_read
albumin -> jack_inhibitor2
albumin -> jack_inhibitor3
albumin -> mesalamine2
albumin -> mesalamine3
albumin -> steroid3
biologics1 -> Fistula
biologics1 -> abscess
biologics1 -> biologics2
biologics1 -> biologics3
biologics1 -> immunomodulators3
biologics1 -> jack_inhibitor3
biologics1 -> mesalamine3
biologics1 -> steroid3
biologics2 -> abscess
biologics2 -> biologics3
biologics2 -> immunomodulators3
biologics2 -> jack_inhibitor3
biologics2 -> mesalamine3
biologics2 -> steroid3
biologics3 -> abscess
bmi -> Fistula
bmi -> abdominal_pain
bmi -> biologics1
bmi -> biologics2
bmi -> biologics3
bmi -> hbi
bmi -> immunomodulators1
bmi -> immunomodulators2
bmi -> immunomodulators3
bmi -> jack_inhibitor1
bmi -> jack_inhibitor2
bmi -> jack_inhibitor3
bmi -> mesalamine1
bmi -> mesalamine2
bmi -> mesalamine3
bmi -> steroid1
bmi -> steroid2
bmi -> stricture_read
crp -> abdominal_pain
crp -> biologics2
crp -> biologics3
crp -> immunomodulators2
crp -> immunomodulators3
crp -> inflammation_read
crp -> jack_inhibitor2
crp -> jack_inhibitor3
crp -> mesalamine2
crp -> mesalamine3
crp -> steroid2
crp -> steroid3
crp -> stricture_read
esr -> abdominal_pain
esr -> biologics2
esr -> biologics3
esr -> immunomodulators2
esr -> immunomodulators3
esr -> inflammation_read
esr -> jack_inhibitor2
esr -> jack_inhibitor3
esr -> mesalamine3
esr -> steroid3
gender -> Fistula
gender -> abdominal_pain
gender -> albumin
gender -> bmi
gender -> crp
gender -> esr
gender -> hbi
gender -> hgb
gender -> inflammation_read
gender -> stricture_read
gender -> wbc
hgb -> abdominal_pain
hgb -> biologics2
hgb -> biologics3
hgb -> immunomodulators2
hgb -> immunomodulators3
hgb -> inflammation_read
hgb -> jack_inhibitor2
hgb -> jack_inhibitor3
hgb -> mesalamine2
hgb -> mesalamine3
hgb -> steroid2
hgb -> steroid3
immunomodulators1 -> Fistula
immunomodulators1 -> abscess
immunomodulators1 -> immunomodulators2
immunomodulators2 -> abscess
immunomodulators2 -> immunomodulators3
immunomodulators3 -> abscess
inflammation_read -> Fistula
inflammation_read -> abdominal_pain
inflammation_read -> biologics2
inflammation_read -> biologics3
inflammation_read -> hbi
inflammation_read -> immunomodulators2
inflammation_read -> immunomodulators3
inflammation_read -> jack_inhibitor2
inflammation_read -> jack_inhibitor3 [pos="-0.858,-0.712"]
inflammation_read -> mesalamine2
inflammation_read -> mesalamine3
inflammation_read -> steroid2
inflammation_read -> steroid3
inflammation_read -> stricture_read
jack_inhibitor1 -> Fistula
jack_inhibitor1 -> abscess
jack_inhibitor1 -> biologics3
jack_inhibitor1 -> immunomodulators3
jack_inhibitor1 -> jack_inhibitor3
jack_inhibitor2 -> abscess
jack_inhibitor2 -> jack_inhibitor3
jack_inhibitor3 -> abscess
mesalamine1 -> Fistula
mesalamine1 -> abscess
mesalamine1 -> mesalamine2
mesalamine2 -> abscess
mesalamine2 -> mesalamine3
mesalamine3 -> abscess
race -> Fistula
race -> abdominal_pain
race -> abscess_history
race -> abscess_wi_1yr
race -> albumin
race -> crp
race -> esr
race -> hgb
race -> inflammation_read
race -> stricture_read
race -> wbc
smoker -> Fistula
smoker -> abscess_history
smoker -> abscess_wi_1yr
smoker -> biologics1
smoker -> biologics2
smoker -> biologics3
smoker -> bmi [pos="-0.957,-1.186"]
smoker -> hbi
smoker -> immunomodulators1
smoker -> immunomodulators2
smoker -> immunomodulators3
smoker -> inflammation_read
smoker -> jack_inhibitor1
smoker -> jack_inhibitor2
smoker -> jack_inhibitor3
smoker -> mesalamine1
smoker -> mesalamine2
smoker -> mesalamine3
smoker -> steroid1
smoker -> steroid2
smoker -> steroid3
smoker -> stricture_read
steroid1 -> Fistula
steroid1 -> abscess
steroid1 -> steroid2
steroid2 -> abscess
steroid2 -> steroid3
steroid3 -> abscess
stricture_read -> abscess [pos="-0.803,-0.597"]
stricture_read -> biologics2
stricture_read -> biologics3
stricture_read -> hbi
stricture_read -> immunomodulators2
stricture_read -> immunomodulators3
stricture_read -> jack_inhibitor2
stricture_read -> jack_inhibitor3
stricture_read -> mesalamine2
stricture_read -> mesalamine3
stricture_read -> steroid2
stricture_read -> steroid3
wbc -> abdominal_pain
wbc -> biologics2
wbc -> biologics3
wbc -> immunomodulators2
wbc -> immunomodulators3
wbc -> inflammation_read
wbc -> jack_inhibitor2
wbc -> jack_inhibitor3
wbc -> mesalamine2
wbc -> mesalamine3
wbc -> steroid2
wbc -> steroid3
}
'
# End DAG Code #####

glm_X = X
glm_y = y

# No patients with prior surgery in 6 months
surgery_6_months_cases = which(as.numeric(df[,16])-1 == 1)
glm_X = glm_X[-surgery_6_months_cases, ]
glm_y = glm_y[-surgery_6_months_cases, ]

# No patients with abscess at imaging
abscess_at_imaging_cases = c(2, 8, 26, 28, 36, 39, 162)
glm_X = glm_X[-abscess_at_imaging_cases, ]
glm_y = glm_y[-abscess_at_imaging_cases, ]

# Subgroup stats #####
get_stats_table = function(df, fistula_col, stat_col){
   
   t = table(df[,fistula_col], df[,stat_col], 
             dnn=c('Fistula', names(df)[stat_col]))
   
   
   return(t)
   
}

df_no_recent_surgery = df[-surgery_6_months_cases, c(8, 61, 66, 68, 69, 71, 94)]
df_no_recent_surgery[,2] = as.factor(df_no_recent_surgery[,2])
df_no_recent_surgery[,3] = as.factor(as.numeric(as.numeric(df_no_recent_surgery[,3])-1!=0))
df_no_recent_surgery[,4] = as.factor(as.numeric(df_no_recent_surgery[,4]!=0))
df_no_recent_surgery[,5] = as.factor(as.numeric(as.numeric(df_no_recent_surgery[,5])-1!=0))
df_no_recent_surgery[,6] = as.factor(df_no_recent_surgery[,6])
df_no_recent_surgery[,7] = as.factor(as.numeric(as.numeric(df_no_recent_surgery[,7])-1!=0))

get_stats_table(df_no_recent_surgery, 1, 2)
get_stats_table(df_no_recent_surgery, 1, 3)
get_stats_table(df_no_recent_surgery, 1, 4)
get_stats_table(df_no_recent_surgery, 1, 5)
get_stats_table(df_no_recent_surgery, 1, 6)
get_stats_table(df_no_recent_surgery, 1, 7)

get_stats_table_p_value = function(df, fistula_col, stat_col){
   
   glm = glm(df[,fistula_col] ~ df[,stat_col], family=binomial)
   p = summary(glm)$coefficients[2, 4]
   return(p)
   
}

get_stats_table_p_value(df_no_recent_surgery, 1, 2)
get_stats_table_p_value(df_no_recent_surgery, 1, 3)
get_stats_table_p_value(df_no_recent_surgery, 1, 4)
get_stats_table_p_value(df_no_recent_surgery, 1, 5)
get_stats_table_p_value(df_no_recent_surgery, 1, 6)
get_stats_table_p_value(df_no_recent_surgery, 1, 7)

# END subgroup stats #####

# Predictors in Minimal Adjustment Set
# Age, Gender, Race, Smoker, BMI,
# Fistula, Sctricture, Inflammation, Abdominal pain,
# Abscess_history, Abscess_within_1yr_prior,
# Steorids_1, Mesalamine_1, Immunomodulators_1, Biologics_1, Jak_Inhibitor_1

preds = c(1, 2, 3, 4, 5, 6, 7, 8, 10, 11, 12, 18, 19, 20, 21)

abdSurg_preds = preds
ostomy_preds = preds
hospitalization_preds = preds
obstruction_preds = preds
abscess_preds = preds

abdSurg_glm = glm(glm_y[, 1] ~ ., data = glm_X[, abdSurg_preds],family=binomial())
summary(abdSurg_glm)

ostomy_glm = glm(glm_y[, 2] ~ ., data=glm_X[, ostomy_preds], family=binomial)
summary(ostomy_glm)

hospitilizations_glm = glm(glm_y[, 3] ~ ., data=glm_X[, hospitalization_preds], family=binomial)
summary(hospitilizations_glm)$coefficients[,1]

obstruction_glm = glm(glm_y[, 4] ~ ., data=glm_X[, obstruction_preds], family=binomial)
summary(obstruction_glm)

abscess_glm = glm(glm_y[, 5] ~ ., data = glm_X[, abscess_preds],family=binomial())
summary(abscess_glm)

library(stringr)

split = function(s){
   
   r = c()
   s_split = str_split(s, '\n')
   for(i in 1:length(s_split[[1]]))
      r[i] = as.numeric(s_split[[1]][i])
   
   return(r[-length(r)])
   
}

junk = '-1.8533
-0.018
0.4258
-0.4785
0.9377
0.0844
0.028
0.0051
1.4414
0.8303
0.0435
0.1953
-0.3781
1.5568
1.0385
1.0654
-0.5175
0.4277
'
split(junk)
round(exp(split(junk)), 4)

# Write glms to excel helper functions ######
get_significances = function(P){
   
   sig = c()
   for (i in 1:length(P)){
      
      sig[i] = ''
      if (P[i] <= .05 & P[i] > .01)
         sig[i] = '*'
      else if (P[i] <= .01 & P[i] > .001)
         sig[i] = '**'
      else if (P[i] <= .001)
         sig[i] = '***'
      
   }
   
   return(sig)
   
}
glm_summary_to_df = function(glm){
   
   names = c('Intercept',
             'Age', 'Gender: Male', 'Race: Black', 'Race: Other', 
             'Smoker: Former smoker', 'Smoker: Smoker', 'BMI',
             'Fistula at imaging', 'Stricture at imaging', 'Inflammation at imaging',
             'Abdominal pain at imaging', 
             'Intra abdominal abscess history', 'Abscess within 1 year prior',
             'Steroids within 2 weeks prior imaging',
             'Mesalamine within 2 weeks prior imaging',
             'Immunomodulators within 2 weeks prior imaging',
             'Biologics within 8 weeks prior imaging',
             'Jak inhibitor within 2 weeks prior imaging')
   
   coef_df = data.frame(summary(glm)$coefficients)
   coef_df = round(coef_df, 4)
   coef_df['Predictor'] = names
   coef_df = data.frame(coef_df[, c(5, 1:4)])

   coef_df['Significance'] = get_significances(coef_df[,5])
   
   return(coef_df)
   
}
# End glms to excel helper functions ######

# Write glms to excel
full_dfs = list()
full_dfs[[1]] = glm_summary_to_df(abdSurg_glm)
full_dfs[[2]] = glm_summary_to_df(ostomy_glm)
full_dfs[[3]] = glm_summary_to_df(hospitilizations_glm)
full_dfs[[4]] = glm_summary_to_df(obstruction_glm)
full_dfs[[5]] = glm_summary_to_df(abscess_glm)

# write_xlsx(full_dfs, 'C:/_local/fistula_study/full_data_glms.xlsx')
# write_xlsx(full_dfs, 'C:/_local/fistula_study/subgroup_data_glms.xlsx')

# End write glms to excel
