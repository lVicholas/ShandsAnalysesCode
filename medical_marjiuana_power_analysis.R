# Medical Marijuana Power Analysis

#
# 1.	What is your age? 
# a)	17-25 years <-- (.35, .65)
# b)	26-60 years 
# c)	>60 years
# I assume b) and c) <-- (.35, .65)
# 
# 3.	Are you currently (choose all that apply)
# ***** All approximately equal *****
# a)	A college student <-- .2
# b)	Employed <-- .2
# c)	Self-employed <-- .2
# d)	Not working or studying at this time (includes those on disability) <-- .2
# e)	Retired <-- .2
# 
# 4.	What is your diagnosis? 
# a)	Ulcerative colitis (UC) <-- (0)
# b)	Crohns disease (CD) <-- (.8, 1)
# c)	Indeterminate Inflammatory Bowel Disease (IBD) <-- (0, .2)
# 
# 8.	My use of marijuana is
# a)	Recreational <-- (.5, 1)
# b)	Medical <-- (.5, 1)
# c)  Combination ** I am assuming medical >= .5 <---> medical
#

library(pwr)
library(ggplot2)

get_power_from_props = function(n, pop_p, trt_p1, trt_p2){
   
   h = ES.h(trt_p1, trt_p2)
   n1 = n * pop_p
   n2 = n * (1 - pop_p)
   
   test = pwr.2p2n.test(h, n1, n2)
   return(test)
   
} 
power_curve_from_props = function(pop_p, trt_p1, trt_p2, n_space=seq(30, 500, 1)){
   
   h = ES.h(trt_p1, trt_p2)
   
   powers = c()
   for (n_id in 1:length(n_space)){
      
      n = n_space[n_id]
      n1 = floor(n * pop_p)
      n2 = floor(n * (1 - pop_p))
      powers[n_id] = pwr.2p2n.test(h, n1, n2)$power
      
   }
   
   frame = data.frame(n = n_space, power = powers)
   return(frame)
   
} # Get sequence of (n, power)
get_min_n_for_pwr = function(pop_p, trt_p1, trt_p2, target_power=.8){
   
   h = ES.h(trt_p1, trt_p2)
   
   for(n in 35:2000){
      
      n1 = floor(n * pop_p)
      n2 = floor( n * (1 - pop_p) )
      current_power = pwr.2p2n.test(h, n1, n2)$power
      
      if (current_power >= target_power)
         return(n)
      
   }
   
   return(-1)
   
} # Min n for power p
power_table = function(pop_p, trt_p1_list, trt_p2_list, target_power=.8){
   
   row_names = c()
   for(i in 1:length(trt_p1_list))
      row_names[i] = paste('trt_p1_', trt_p1_list[i], '')
   
   col_names = c()
   for(i in 1:length(trt_p2_list))
      col_names[i] = paste('trt_p2_', trt_p2_list[i], '')
   
   n_mat = matrix(0, nrow=length(trt_p1_list), ncol=length(trt_p2_list))
   for(i in 1:nrow(n_mat)){
      for(j in 1:ncol(n_mat)){
         trt_p1 = trt_p1_list[i]
         trt_p2 = trt_p2_list[j]
         p = get_min_n_for_pwr(pop_p, trt_p1, trt_p2)
         n_mat[i, j] = p
      }
   }
   
   df = data.frame(n_mat)
   colnames(df) = col_names
   row.names(df) = row_names
   
   return(df)
}
get_power_curves_given_pop_prop = function(pop_prop, trt_p1_list, trt_p2_list){
   
   if(length(trt_p1_list) != length(trt_p2_list))
      print('proportion lists have different lengths')
   
   curves = list()
   
   for(i in 1:length(trt_p1_list)){
      
      p1 = trt_p1_list[i]
      p2 = trt_p2_list[i]
      
      curve = power_curve_from_props(pop_prop, p1, p2)
      curve['power'] = curve['power']
      curves[[i]] = curve
      
   }
   
   return(curves)
   
}
plot_power_curves_given_pop_prop = function(pop_prop, trt_p_df, cartesian=F){
   
   trt_p1_list = trt_p_df[, 1]
   trt_p2_list = trt_p_df[, 2]
   
   if(cartesian){
      
      cartesian_product = expand.grid(trt_p1_list, trt_p2_list)
      trt_p1_list = cartesian_product[, 1]
      trt_p2_list = cartesian_product[, 2]
      
   }
   
   ns = 500 - 30 + 1
   
   curve_names = c()
   for(i in 1:length(trt_p1_list)){
      
      cns = (i-1)*ns + 1
      cne = i*ns
      
      p1 = trt_p1_list[i]
      p2 = trt_p2_list[i]
      curve_name = paste('(p1=', p1, ', p2=', p2, ')', sep='')  
      curve_names[ cns:cne ] = curve_name
      
   }
   
   p1_name = colnames(trt_p_df)[1]
   p2_name = colnames(trt_p_df)[2]
   legend_name = paste(p1_name, ', ', p2_name, ', pop_split=', pop_prop, sep='')
   
   
   dfs = get_power_curves_given_pop_prop(pop_prop, trt_p1_list, trt_p2_list)
   DF = do.call('rbind', dfs)
   DF = cbind(DF, curve_names)
   colnames(DF) = c('n', 'power', 'group')

   plot = (ggplot(DF, aes(x=n, y=power, color=factor(group))) 
           + geom_point() 
           + geom_hline(yintercept=.8, col='black', size=1.3)
           + scale_color_discrete(name=legend_name)
   )
   return(plot)
                 
}

### ### ### ### ### ### ### ### Age vs. marijuana use ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ######### ### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ######### ### ### ### ### ### ### ### ### ### ###

### ### ### Age vs. recreational marijuana use ### ### ###

age_17_25_population_proportion = .5

use_recreational_17_25 = seq(.25, 1.0, .1)
use_recreational_26_60 = seq(.05, .35, .1)

age_vs_recreational_df = data.frame(recr_17_25=use_recreational_17_25,
                                    recr_26_60=use_recreational_26_60
                                    )

plot_power_curves_given_pop_prop(age_17_25_population_proportion, 
                                 age_vs_recreational_df
                                 #cartesian=T
                                )

# Rows --> 17-25 proportion recreational use
# Columns --> 26+ proportion recreational use
# Power = .8
df = power_table(age_17_25_population_proportion, use_recreational_17_25, use_recreational_26_60)
View(df)

### ### ### Age vs. medical marijuana use ### ### ###

age_17_25_population_proportion = .5

use_medical_17_25 = seq(.10, .45, .05)
use_medical_26_60 = seq(.10, .55, .05)

# Rows --> 17-25 proportion recreational use
# Columns --> 26+ proportion recreational use
# Power = .8
df = power_table(age_17_25_population_proportion, use_medical_17_25, use_medical_26_60)
View(df)

### ### ### Type of marijuana use in young people ### ### ###

age_17_25_population_proportion = .5

use_recreational_17_25 = seq(.10, .45, .05)
use_medical_17_25 = seq(.10, .55, .05)

# Rows --> 17-25 proportion recreational use
# Columns --> 26+ proportion recreational use
# Power = .8
df = power_table(age_17_25_population_proportion, use_recreational_17_25, use_medical_17_25)
View(df)

### ### ### Type of marijuana use in old people ### ### ###

age_17_25_population_proportion = .5

use_recreational_26_60 = seq(.05, .35, .05)
use_medical_26_60 = seq(.10, .55, .05)

# Rows --> 17-25 proportion recreational use
# Columns --> 26+ proportion recreational use
# Power = .8
df = power_table(age_17_25_population_proportion, use_recreational_17_25, use_medical_17_25)
View(df)

### ### ### ### ### ### ### ### Employment vs. marijuana use ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ######## ### ### ### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ######## ### ### ### ### ### ### ### ### ### ### ### ###

# Diagnosis vs. recreational marijuana use

crohns_population_proportion = .8

use_recreational_crohns = seq(.10, .65, .05)
use_recreational_ibd = seq(.10, .65, .05)

df = power_table(crohns_population_proportion, use_recreational_crohns, use_recreational_ibd)
View(df)

# Diagnosis vs medical marijuana use

crohns_population_proportion = .8

use_medical_crohns = seq(.15, .65, .05)
use_medical_ibd = seq(.10, .65, .05)

df = power_table(crohns_population_proportion, use_medical_crohns, use_medical_ibd)
View(df)

