library(data.table)
library(dplyr)
library(ggplot2)
library(lavaan)
library(VIM)
library(mice)

data=read.csv("data.csv")
mydata=as.data.table(copy(data))

#### data preparation



## replace 888 as missing value(NA)
mydata[mydata==888]=NA
typeof(mydata$vocabprocess_processing_speed_target)
mydata$pvt_mean_rt= sapply(mydata$pvt_mean_rt, as.numeric)
typeof(mydata$pvt_mean_rt)
## replace 0 as missing value in wasi_sum_rawscores
mydata$wasi_sum_rawscores[mydata$wasi_sum_rawscores==0]=NA

### catergorical replacement
## gender
gender=function(x){
  if(x==0){
    x="male"
  }else if(x==1){
    x="female"
  }
}

mydata$gender=sapply(mydata$gender, gender)

## diagnosis
diagn=function(x){
  if(x==0){
    x="Non-autistic"
  }else if(x==1){
    x="Autistic"
  }
}
mydata$diagnosis=sapply(mydata$diagnosis, diagn)

## where_english
where_e=function(x){
  if(x==1){
    x="Home"
  }else if(x==2){
    x="nursery"
  }else if(x==3){
    x="Playground"
  }else if(x==4){
    x="School"
  }
}




mydata$where_english=sapply(mydata$where_english, where_e)

mydata[, age_y:=mydata$age_m/12]



##################################################################################
nrow(na.omit(mydata))
## only 29 rows without any missing data, so we can't use likewise deletion (complete-case analysis)

## number of missing value for each column
sapply(mydata, function(x) sum(is.na(x)))
## number of missing value in SCQ is 38 because autistic group doens't have it


## k-Nearst Neighbour Imputation
library(MASS)
library(VIM)
# install.packages("MASS")
##k=5 for default
## k~sqrt(N)
mydata_knn = kNN(mydata[, -c(5, 12, 18)], k=9)
mydata_knn=subset(mydata_knn, select = part_no:pvt_count_falsestarts )
mydata_knn[, tomi_compmean:=mydata$tomi_compmean]



mydata_knn$tomi_compmean=ifelse(is.na(mydata_knn$tomi_compmean), 
                                round((mydata_knn$tomi_basic+
                                         mydata_knn$tomi_advanced +
                                         mydata_knn$tomi_early)/3, 1),
                                mydata_knn$tomi_compmean)
#detach("package:MASS", TRUE)


mydata_knn[, et_falsebelief_testtrial_preference_score:=
             mydata$et_falsebelief_testtrial_preference_score]
mydata_knn$et_falsebelief_testtrial_preference_score=ifelse(
  is.na(mydata_knn$et_falsebelief_testtrial_preference_score), 
  round((mydata_knn$et_falsebelief_Testtrial_dwell_time_to_correct-
           mydata_knn$et_falsebelief_testtrial_dwell_time_to_incorrect),1),
  mydata_knn$et_falsebelief_testtrial_preference_score)


mydata_knn[, SCQ:=mydata$SCQ]


### create a new column for year of learning second language
mydata_knn[, year_biligual := round(mydata_knn$age_m /12 - mydata_knn$age_acquisition, 3)]
mydata_knn=as.data.frame(mydata_knn)
## reorder the column 
mydata_knn=mydata_knn[, colnames(mydata_knn)[c(1:10, 41, 11:15, 42, 16:40, 43, 44)]]


#########################################
## data processing for where_english


## since participants may have broad defintion for the choice in where_english, 
## they may define playground and nursery as childcare or environment where parents are not presented. 
## Thus, I group nursery and playgroud together as a new category, which called "Other"

## function to group nursery and playground as new category, "Other"
where_e_combine =function(x){
  if(x== "nursery" || x== "Playground"){
    x= "Other"
  }else{
    x=x
  }
}

mydata_knn$where_english=sapply(mydata_knn$where_english, where_e_combine)




###############################################
library(dplyr)
detach("package:MASS", TRUE)

mydata_knn$pvt_mean_rt=as.numeric(mydata_knn$pvt_mean_rt)

## not need to assume the normal distribution 
# independent 2-group Mann-Whitney U Test

## contain diagnosis, but no SCQ
library(MASS)
mydata_knn_numeric=mydata_knn[, -c(1,2,25, 43)]
wilcox_result=lapply(mydata_knn_numeric[, -2], function(x) wilcox.test(x~ diagnosis, data = mydata_knn_numeric))
pval =sapply(wilcox_result, `[[`, "p.value")
pval =sapply(pval, round, 4)
mydata_knn_numeric=as.data.table(mydata_knn_numeric)
mydata_knn_au_num=mydata_knn_numeric[diagnosis=="Autistic"]
mydata_knn_nau_num=mydata_knn_numeric[diagnosis=="Non-autistic"]

mean_au=lapply(mydata_knn_au_num[, -2], mean)
mean_au=sapply(mean_au, round, 2)
mean_nau=lapply(mydata_knn_nau_num[, -2], mean)
mean_nau=sapply(mean_nau, round, 2)

sign_diff=function(x){
  if(x>0.05){
    return("Not significantly different")
  }else if(x<0.05){
    return("significantly different")
  }
}

sign_diff_2=sapply(pval, sign_diff)

diag_stat=as.data.frame(cbind( mean_au, mean_nau, pval, sign_diff_2))


diag_stat[sign_diff_2=="significantly different", ]
nrow(diag_stat[sign_diff_2=="significantly different", ])


## tranform dataframe to latex output
library(knitr)
kable(diag_stat[, -4], "latex")



###############################################
## see the difference between where_english
## Kruskal-Wallis rank sum test 
## non-parametric test (doesn't assume normal distribution )
## extenstion of  Mann–Whitney U test, it compares more than two groups at the same time
mydata_knn_numeric_we=mydata_knn[, -c(1,2, 4,43)]

wilcox_result2=lapply(mydata_knn_numeric_we[, -22], 
                      function(x) kruskal.test(x~ where_english,
                                               data = mydata_knn_numeric_we))
pval2 =sapply(wilcox_result2, `[[`, "p.value")
pval2=sapply(pval2, round, 4)

mydata_knn_numeric_we=as.data.table(mydata_knn_numeric_we)
mydata_knn_home=mydata_knn_numeric_we[where_english=="Home"]
mydata_knn_school=mydata_knn_numeric_we[where_english=="School"]
mydata_knn_other=mydata_knn_numeric_we[where_english=="Other"]

mean_home=lapply(mydata_knn_home[, -22], mean)
mean_home=sapply(mean_home, round, 2)

mean_school=lapply(mydata_knn_school[, -22], mean)
mean_school=sapply(mean_school, round, 2)

mean_other=lapply(mydata_knn_other[, -22], mean)
mean_other=sapply(mean_other, round, 2)

sign_diff_we=sapply(pval2, sign_diff)

we_stat=as.data.frame(cbind( mean_home, mean_school, mean_other, pval2, sign_diff_we))

we_stat[sign_diff_we=="significantly different", ]
nrow(we_stat[sign_diff_we=="significantly different", ])

## only tomi_basic and pvt_count_falsestarts in exam score are significantly different for three groups in where english

## delete the bi_exposure variable
we_stat=we_stat[-c(1, 4, 15:21,39 ) ,]
nrow(we_stat[we_stat$sign_diff_we=="significantly different", ])
we_stat[we_stat$sign_diff_we=="significantly different", ]

## tranform dataframe to latex output
library(knitr)
kable(we_stat[, -5], "latex")


#########################################
## split into autistic and non-autistic group
au=mydata[mydata$diagnosis=="Autistic"]

nau=mydata[mydata$diagnosis=="Non-autistic"]



###########################################################
## Factor Analysis

library(psych)
library(GPArotation)



social=mydata_knn[, c(8:16)]

cor_social=as.data.frame(cor(social))
cor_social[cor_social < 0.9 | cor_social ==1] <- ""
## cancel the ef_preference score, then we can use factor analysis


social=mydata_knn[, c(8:10, 12:16)]
social_scale=as.data.frame(scale(social))
## KMO measure (if MSA<0.5, factor analysis is not appropriate)
KMO(cor(social))

ev=eigen(cor(social))
ev$values
factor=c(1:8)
scree=data.frame(factor, ev$values)
plot(scree, ylab="Eigen value", type="b", main="screeplot for social cognition")
abline(h=1, col="blue")


social_fa2=fa(social, nfactors =2 , rotate="varimax", fm="ml")
# factanal(social, 2, rotation="promax")
## cummulative variation 0.59
social_fa2
print(social_fa2$loadings, cutoff = 0.3)

social_fa2_loading =as.data.frame(social_fa2$loadings[, 1:2])
social_fa2_loading$ML2 = sapply(social_fa2_loading$ML2, round, 3)
social_fa2_loading$ML1 = sapply(social_fa2_loading$ML1, round, 3)

library(knitr)
kable(social_fa2_loading, "latex")

#########################################################

executive=mydata_knn[, c(26:42)]
## 17 variables

cor_exe=as.data.frame(cor(executive))
cor_exe[abs(cor_exe) < 0.8 | cor_exe ==1] <- ""
### pvt_mean_rt and pvt_mean_lapse_rt has cor=0.9
## delete pvt_mean_lapse_rt since it has 19 missing value
## flanker_mean_rt and flanker_mean_rt_inconverge has cor=0.82, 
## delete flanker_mean_rt_inconverge (13 missing value)
## cor(brief_raw_working_memory, brief_raw_task_monitor)=0.809
## cor(brief_raw_plan_organise, brief_raw_task_monitor)=0.8256
## delete brief_raw_task_monitor

executive=mydata_knn[, c(26:32, 34,35, 37:40, 42)]

KMO(cor(executive))

## delete pvt_count_falsestarts since MSA=0.44
executive=mydata_knn[, c(26:32, 34,35, 37:40)]
KMO(cor(executive))
exe_scale=as.data.frame(scale(executive))

ev=eigen(cor(executive))
ev$values
factor=c(1:13)
scree=data.frame(factor, ev$values)
plot(scree, ylab="Eigen value", type="b", main="screeplot for executive")
abline(h=1, col="blue")
fa.parallel(executive,  fm="ml", fa="fa")
## choose 2

exe_fa2=fa(executive, nfactors = 2, rotate="varimax", fm="ml")
exe_fa2
print(exe_fa2$loadings, cutoff = 0.3)


exe_fa2_loading =as.data.frame(exe_fa2$loadings[, 1:2])
exe_fa2_loading$ML2 = sapply(exe_fa2_loading$ML2, round, 3)
exe_fa2_loading$ML1 = sapply(exe_fa2_loading$ML1, round, 3)

library(knitr)
kable(exe_fa2_loading, "latex")

####################################################################
### comfirmatory factor analysis (CFA)
library(lavaan)
library(semPlot)

## standardize data first

knn_std=as.data.frame(scale(mydata_knn[, -c(1,2,4,25,43)]))

where_e_dummy=model.matrix(~where_english, data = mydata_knn)
where_e_dummy=as.data.frame(where_e_dummy[,-1])
# where_e_dummy=as.data.frame(lapply(where_e_dummy, as.factor))


knn_std = as.data.frame (cbind(knn_std, data$diagnosis, 
                               where_e_dummy[, 1], where_e_dummy[, 2]))
colnames(knn_std) = c(colnames(mydata_knn[, -c(1,2,4,25,43)]), 
                      "diagnosis", colnames(where_e_dummy))


################################################################

knn_std=as.data.table(knn_std)
non_au_std=knn_std[knn_std$diagnosis==0]
au_std = knn_std[knn_std$diagnosis==1]


#### the terrible cfa model without doing EFA first
cfa_model=' 
social_cognition =~ tomi_early + tomi_basic + tomi_advanced + tomi_compmean + 
tom_tb_totalscore + et_figurestask_dwell_time_interacting + 
et_figurestask_dwell_time_not_interacting + 
et_falsebelief_Testtrial_dwell_time_to_correct + 
et_falsebelief_testtrial_dwell_time_to_incorrect + 
et_falsebelief_testtrial_preference_score

executive_function =~ brief_raw_inhibit + brief_raw_self.monitor +
brief_raw_shift + brief_raw_emotional_control + 
brief_raw_initiate + brief_raw_working_memory + brief_raw_plan_organise +
brief_raw_task_monitor + brief_raw_organisation_of_materials +
flanker_percenterrors_congruent + flanker_percenterrors_incongruent +
flanker_mean_rt_congruent + flanker_mean_rt_incongruent +
pvt_mean_rt + pvt_number_of_lapses + pvt_mean_lapse_rt +
pvt_count_falsestarts

language =~ vocabprocess_processing_speed_target + bpvs_raw

'
fit_cfa=cfa(cfa_model, data=knn_std)
summary(fit_cfa, fit.measures=TRUE, standardized=TRUE)

######################################################
#### for three cognition outcome (see the covariance)

cfa_fa_sep='

social_et =~ et_falsebelief_Testtrial_dwell_time_to_correct +
et_falsebelief_testtrial_dwell_time_to_incorrect +
et_figurestask_dwell_time_interacting +
et_figurestask_dwell_time_not_interacting 

social_tom =~ tomi_early + tomi_basic + tomi_advanced + tom_tb_totalscore

exe_brief =~brief_raw_inhibit + brief_raw_shift + brief_raw_emotional_control 
+ brief_raw_initiate + brief_raw_working_memory + brief_raw_plan_organise  
 +brief_raw_self.monitor +
brief_raw_organisation_of_materials 

exe_flanker_pvt =~ flanker_mean_rt_congruent +
flanker_percenterrors_congruent + flanker_percenterrors_incongruent +
pvt_number_of_lapses + pvt_mean_rt  

executive =~ exe_brief + exe_flanker_pvt

social_cognition =~ social_tom + social_et 

language =~ vocabprocess_processing_speed_target + bpvs_raw
'

fit_cfa_sep=cfa(cfa_fa_sep, data=knn_std)
summary(fit_cfa_sep, fit.measures=TRUE, standardized=TRUE)
semPaths(fit_cfa_sep,"std", edge.label.cex = 1 , color = "lightblue", nCharNodes = 6)


###########################################################
## social

social_cfa_all='
social_tom =~ tomi_early + tomi_basic + tomi_advanced + 
tom_tb_totalscore   

social_et =~ et_falsebelief_Testtrial_dwell_time_to_correct +
et_falsebelief_testtrial_dwell_time_to_incorrect +
et_figurestask_dwell_time_interacting +
et_figurestask_dwell_time_not_interacting 

social=~ social_tom +social_et

social ~  bilec_total_input +  bilec_total_output  + age_acquisition+
wasi_sum_rawscores +age_m +diagnosis
'

social_cfa_all_model=cfa(social_cfa_all, data=knn_std, ordered="diagnosis")
summary(social_cfa_all_model, fit.measures=TRUE, standardized=TRUE)

semPaths(social_cfa_all_model,"std"
         , edge.label.cex = 1 , color = "lightblue", nCharNodes = 7)



###########################################################
## executive 

exe_cfa_all='
exe_brief =~brief_raw_inhibit + brief_raw_shift + brief_raw_emotional_control 
+ brief_raw_initiate + brief_raw_working_memory + brief_raw_plan_organise  
 +brief_raw_self.monitor +
brief_raw_organisation_of_materials 

exe_flanker_pvt =~ flanker_mean_rt_congruent  +
flanker_percenterrors_congruent + flanker_percenterrors_incongruent +
pvt_number_of_lapses + pvt_mean_rt 

executive =~ exe_brief + exe_flanker_pvt

executive ~   bilec_total_input + bilec_total_output  + age_acquisition+
wasi_sum_rawscores +age_m +diagnosis
'
##  + brief_raw_task_monitor
## 
exe_cfa_all_fit = cfa(exe_cfa_all, data=knn_std, ordered="diagnosis")
summary(exe_cfa_all_fit, fit.measures=TRUE, standardized=TRUE)


#####################################################
## language 

language_cfa_all='
language =~ vocabprocess_processing_speed_target + bpvs_raw

language ~  bilec_total_input + bilec_total_output  + age_acquisition+
wasi_sum_rawscores +age_m +diagnosis
' 

fit_language_cfa_all=cfa(language_cfa_all, data=knn_std, ordered = "diagnosis")
summary(fit_language_cfa_all, fit.measures=TRUE, standardized=TRUE)


##################################################################
### language group by diagnosis
language_cfa='
language =~ vocabprocess_processing_speed_target + bpvs_raw
language ~ bilec_total_input + bilec_total_output  + age_acquisition
+wasi_sum_rawscores +age_m 
' 

fit_language_cfa_gp=cfa(language_cfa, data=knn_std, group = "diagnosis")
summary(fit_language_cfa_gp, fit.measures=TRUE, standardized=TRUE)





##################################################################
### social cognition  group by diagnosis
social_cfa='
social_tom =~ tomi_early + tomi_basic + tomi_advanced +tom_tb_totalscore   

social_et =~ et_falsebelief_Testtrial_dwell_time_to_correct +
et_figurestask_dwell_time_interacting +
et_figurestask_dwell_time_not_interacting 
+et_falsebelief_testtrial_dwell_time_to_incorrect 

social=~ social_tom +social_et   

social ~ bilec_total_input + bilec_total_output  + age_acquisition+
wasi_sum_rawscores +age_m 
'

social_cfa_fit_group = cfa(social_cfa, data=knn_std, group = "diagnosis")
summary(social_cfa_fit_group, fit.measures=TRUE, standardized=TRUE)


##################################################################
### Executive function  group by diagnosis


exe_cfa_2='
exe_brief =~brief_raw_inhibit + brief_raw_shift + brief_raw_emotional_control 
+ brief_raw_initiate + brief_raw_working_memory + brief_raw_plan_organise  
 +brief_raw_self.monitor + brief_raw_task_monitor +
brief_raw_organisation_of_materials 

exe_flanker_pvt =~ flanker_mean_rt_incongruent  +
flanker_percenterrors_congruent + flanker_percenterrors_incongruent +
 pvt_mean_rt 

executive =~ pvt_number_of_lapses +exe_brief + exe_flanker_pvt


executive ~ bilec_total_input + bilec_total_output  + age_acquisition+
wasi_sum_rawscores +age_m 
'


exe_cfa_fit_gp = cfa(exe_cfa_2, data=knn_std, group = "diagnosis", check.gradient = FALSE)
summary(exe_cfa_fit_gp, fit.measures=TRUE, standardized=TRUE)



##########################################################
## pca

## pca for language

pca_language = prcomp(formula= ~ vocabprocess_processing_speed_target + bpvs_raw,
                      data = mydata_knn, scale=TRUE)


vars_language=pca_language$sdev^2
prop_language=vars_language/sum(vars_language)

cum_prop_language=cumsum(prop_language)
plot(cum_prop_language, main="cummulative proportion for language")

cum_prop_language


####################################################
### language vs bilingualism
language_pca_x1=pca_language$x[, 1]
data_l=as.data.frame(cbind(language_pca_x1, knn_std[, c(1,4,15:21)], 
                           factor(knn_std$diagnosis)))
colnames(data_l)=c("language_pca_x1", colnames(knn_std[, c(1,4,15:21)]), "diagnosis")
fit_pca_l=lm(language_pca_x1~., data=data_l)
summary(fit_pca_l)

fit_pca_l_step=step(fit_pca_l, direction = "both", trace=FALSE)
summary(fit_pca_l_step)

par(mfrow=c(2,2))
plot(fit_pca_l_step)
par(mfrow=c(1,1))



##### autistic group 
data_l=as.data.table(data_l)
data_l_au=data_l[diagnosis==1]
data_l_au=data_l_au[, -11]

fit_pca_l_au=lm(language_pca_x1~., data=data_l_au)
summary(fit_pca_l_au)

fit_pca_l_step_au=step(fit_pca_l_au, direction = "both", trace=FALSE)
summary(fit_pca_l_step_au)

par(mfrow=c(2,2))
plot(fit_pca_l_step_au)
par(mfrow=c(1,1))



### Non-Autistic group
##### autistic group 
data_l=as.data.table(data_l)
data_l_nau=data_l[diagnosis==0]
data_l_nau=data_l_nau[, -11]

fit_pca_l_nau=lm(language_pca_x1~., data=data_l_nau)
summary(fit_pca_l_nau)

fit_pca_l_step_nau=step(fit_pca_l_nau, direction = "both", trace=FALSE)
summary(fit_pca_l_step_nau)

par(mfrow=c(2,2))
plot(fit_pca_l_step_nau)
par(mfrow=c(1,1))



