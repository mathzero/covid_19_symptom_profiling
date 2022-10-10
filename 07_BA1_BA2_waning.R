
#' First clear the environment of variables
rm(list=ls(all=TRUE))
# get root director of project
root.dir <- getwd()
setwd(dir = "/Group/react2_study5/report_phases_combined/projects/omicron_symptom_profiling/")
outpath <- paste0(root.dir,"/output/")
figpath <-  paste0(root.dir,"/plots/")



source("E:/Group/functions/load_packages.R", local = T)
source("E:/Group/functions/full_join_multiple_datasets.R", local = T)
source("E:/Group/functions/wrangle_cleaner_functions.R", local = T)
source("E:/Group/functions/cats_and_covs.R", local = T)
source("E:/Group/react2_study5/report_phases_combined/projects/function_scripts/create_subfolder.R", local = T)
source("E:/Group/react2_study5/report_phases_combined/projects/function_scripts/forest_plot.R", local = T)
source("E:/Group/react2_study5/report_phases_combined/projects/function_scripts/save_styled_table.R", local = T)
source("E:/Group/react2_study5/report_phases_combined/projects/function_scripts/stability_selection.R", local = T)

# 
# 
# source("E:/Group/react2_study5/report_phases_combined/projects/delta_symptom_prediction/code/00_bits_and_pieces.R", local = T)
# source("E:/Group/react2_study5/report_phases_combined/projects/delta_symptom_prediction/code/00_functions.R", local = T)
# source("E:/Group/react2_study5/report_phases_combined/projects/symptom_prediction_children/code/00_bits_and_pieces.R", local = T)


#' Pull in packages needed
package.list <- c("prevalence","mgcv","knitr","MASS","kableExtra","table1","dplyr",
                  "tidyr", "pheatmap","scales","OverReact","ggstance",
                  "ggplot2","ggsci", "RColorBrewer", "tidyverse", "lubridate", 
                  "readr","ggthemes", "questionr", "foreach", "doParallel","withr",
                  "patchwork","randomcoloR","focus","ggtext","showtext"
)

load_packages(package.list)


# Import REACT-1 data -----------------------------------------------------
source("E:/Group/react2_study5/report_phases_combined/projects/omicron_symptom_profiling/code/00_bits_and_pieces.R", local = T)
source("E:/Group/react2_study5/report_phases_combined/projects/omicron_symptom_profiling//code/00_data_prep.R", local = T)

# create subfolder
createMySubfolder(subfolderName = "waning_analysis")


# Date edit
dfRes <- dfRes %>% mutate(vaccinethird=as.Date(vaccinethird, format = "%m/%d/%Y"),
                          vaccine_second=as.Date(vaccine_second, format = "%m/%d/%Y",origin = "01/01/1970"),
                          vaccine_second_sym=as.Date(vaccine_second_sym, format = "%m/%d/%Y",origin = "01/01/1970"),
                          vaccinesecondsym=as.Date(vaccinesecondsym, format = "%m/%d/%Y",origin = "01/01/1970"),
                          days_since_boost = as.numeric(d_comb-vaccinethird),
                          days_since_second_vax = as.numeric(case_when(!is.na(vaccine_second) ~ d_comb-vaccine_second,
                                                            !is.na(vaccine_second_sym) ~ d_comb-vaccine_second_sym,
                                                            !is.na(vaccinesecondsym) ~ d_comb-vaccinesecondsym)),
                          
                          calendar_time = as.numeric(d_comb-as.Date("2022-01-01")),
                          impact_a_lot=case_when(covidability7==1 ~ 1,
                                                  T ~ 0),
                          impact_a_lot_or_a_litte=case_when(covidability7==1 ~ 1,
                                                            covidability7==2 ~ 1,
                                                 T ~ 0),
                          BA2 = case_when(variant_inferred_detail=="BA.2 (Omicron)" ~ 1,
                                          variant_inferred_detail=="BA.1 (Omicron)" ~ 0,
                                          T ~ NA_real_)
                          )

table(dfRes$days_since_second_vax)
# create symptom count variable
dfRes$symptom_count_26 <- rowSums(dfRes[,covid_yesnos], na.rm=T)
dfRes$symptom_count_4 <- rowSums(dfRes[,sympnames_type_df$symptom_code[1:4]], na.rm=T)

dfRes$vaccine_second_sym %>% table()

# subset data
dfRes_pos=dfRes %>% filter(variant_inferred_detail%in% c("BA.1 (Omicron)","BA.2 (Omicron)"), round %in%c(17:19))
dfRes_pos_boost=dfRes %>% filter(variant_inferred_detail%in% c("BA.1 (Omicron)","BA.2 (Omicron)"), 
                                 round %in%c(17:19), !covidability7%in% c(4:5),
                                 vaccdose==3,days_since_boost>=14)
dfRes_pos_2_or_3=dfRes %>% filter(variant_inferred_detail%in% c("BA.1 (Omicron)","BA.2 (Omicron)"), 
                                 round %in%c(17:19), !covidability7%in% c(4:5),
                                 vaccdose%in%c(2,3)) %>% 
  mutate(boosted = factor(case_when(vaccdose==2~ "No",
                             T ~ "Yes"), levels = c("No","Yes")),
         time_since_last_vax= as.numeric(case_when(!is.na(days_since_boost) ~ days_since_boost,
                                        !is.na(days_since_second_vax) ~ days_since_second_vax,
                                        # T ~ NA_real_
                                        ))) %>% 
  filter(time_since_last_vax>=14)
table(dfRes_pos_2_or_3$age_group_named)
# 
# 
# # Waning analysis for infection in r17-19, BA1/BA2 ----------------------
# dfRes_pos_boost$age_group_named <- factor(dfRes_pos_boost$age_group_named,
#                                           levels = c("18-24" ,"25-34","35-44","45-54","55-64","65-74", "74+"))
# 
# mod_glm_severe = glm(formula = as.formula("impact_a_lot ~ age_group_named + sex  + BA2 + days_since_boost +I(days_since_boost^2)+ 
#                                           round"),
#               data = dfRes_pos_boost,family = "binomial")
# 
# 
# mod_glm_severe_calendar_time = glm(formula = as.formula("impact_a_lot ~ age_group_named + sex  +
#                                   prior_covid_28plus+ BA2 + days_since_boost +I(days_since_boost^2)+ 
#                                           round + calendar_time"),
#                      data = dfRes_pos_boost,family = "binomial")
# 
# mod_glm_a_lot_or_a_little = glm(formula = 
#                                   as.formula("impact_a_lot_or_a_litte ~ age_group_named + sex  + BA2 + days_since_boost+
#                                              I(days_since_boost^2)+ round"),
#                      data = dfRes_pos_boost,family = "binomial")
# 
# 
# jtools::summ(mod_glm_severe,exp=T)
# jtools::summ(mod_glm_severe_calendar_time,exp=T)
# jtools::summ(mod_glm_a_lot_or_a_little,exp=T)
# 
# # output severe table
# tab_severity=jtools::summ(mod_glm_severe,exp=T)$coeftable %>% as.data.frame() %>% 
#   select(-`z val.`)
# tab_severity$OR_with_CI=paste0(round(tab_severity$`exp(Est.)`,2), ", [", round(tab_severity$`2.5%`,2),",",
#                                round(tab_severity$`97.5%`,2),"]")
# tab_severity$variable = rownames(tab_severity)
# tab_severity <- tab_severity %>% select(variable,OR_with_CI,p)
# colnames(tab_severity)= c("Independent variable","Odds ratio with CI","P-value")
# # tab_severity <- tab_severity %>% select(variable, everything())
# tab_severity$`Independent variable` <- c("Intercept","Age group: 25-34",
#                                          "Age group: 35-44",
#                                          "Age group: 45-54",
#                                          "Age group: 55-64",
#                                          "Age group: 65-74",
#                                          "Age group: 74+",
#                                          "Sex: male",
#                                          # "Prior COVID-19: Yes",
#                                          "BA.2",
#                                          "Days since booster vaccine",
#                                          "[Days since booster vaccine]^2",
#                                          "Round 18",
#                                          "Round 19"
#                                          )
# 
# 
# 
# # output severe tablewith calendar time
# tab_severity_calendar_time=jtools::summ(mod_glm_severe_calendar_time,exp=T)$coeftable %>% as.data.frame() %>% 
#   select(-`z val.`)
# tab_severity_calendar_time$OR_with_CI=paste0(round(tab_severity_calendar_time$`exp(Est.)`,2), ", [", round(tab_severity_calendar_time$`2.5%`,2),",",
#                                round(tab_severity_calendar_time$`97.5%`,2),"]")
# tab_severity_calendar_time$variable = rownames(tab_severity_calendar_time)
# tab_severity_calendar_time <- tab_severity_calendar_time %>% select(variable,OR_with_CI,p)
# colnames(tab_severity_calendar_time)= c("Independent variable","Odds ratio with CI","P-value")
# # tab_severity_calendar_time <- tab_severity_calendar_time %>% select(variable, everything())
# tab_severity_calendar_time$`Independent variable` <- c("Intercept","Age group: 25-34",
#                                          "Age group: 35-44",
#                                          "Age group: 45-54",
#                                          "Age group: 55-64",
#                                          "Age group: 65-74",
#                                          "Age group: 74+",
#                                          "Sex: male",
#                                          "Prior COVID-19 (28+ days ago): Yes",
#                                          "BA.2",
#                                          "Days since booster vaccine",
#                                          "[Days since booster vaccine]^2",
#                                          "Round 18",
#                                          "Round 19",
#                                          "Calendar time (days since 1 Jan 2022)"
# )
# 
# 
# 
# # output severe table 2
# tab_severity_2=jtools::summ(mod_glm_a_lot_or_a_little,exp=T)$coeftable %>% as.data.frame() %>% 
#   select(-`z val.`)
# colnames(tab_severity_2)= c("Odds ratio","Lower (95% CI)","Upper (95% CI)","P-value")
# tab_severity_2$variable = rownames(tab_severity_2)
# tab_severity_2 <- tab_severity_2 %>% select(variable, everything())
# 
# tab_severity$pval <- OverReact::pValEpiConverter(tab_severity$`P-value`)
# tab_severity_2$pval <- OverReact::pValEpiConverter(tab_severity_2$`P-value`)
# 
# 
# 
# 
# # 2/3 vax Waning analysis for infection in r17-19, BA1/BA2 ----------------------
# 
# dfRes_pos_2_or_3$age_group_named <- factor(dfRes_pos_2_or_3$age_group_named,
#                                           levels = c("18-24" ,"25-34","35-44","45-54","55-64","65-74", "74+"))
# 
# 
# mod_glm_severe_2_3 = glm(formula = as.formula("impact_a_lot ~ age_group_named + sex  + BA2 + boosted + time_since_last_vax +I(time_since_last_vax^2)+ 
#                                           round"),
#                      data = dfRes_pos_2_or_3,family = "binomial")
# 
# # 
# # mod_glm_severe_calendar_time = glm(formula = as.formula("impact_a_lot ~ age_group_named + sex  +
# #                                   prior_covid_28plus+ BA2 + days_since_boost +I(days_since_boost^2)+ 
# #                                           round + calendar_time"),
# #                                    data = dfRes_pos_2_or_3,family = "binomial")
# # 
# # mod_glm_a_lot_or_a_little = glm(formula = 
# #                                   as.formula("impact_a_lot_or_a_litte ~ age_group_named + sex  + BA2 +boosted + days_since_boost+
# #                                              I(days_since_boost^2)+ round"),
# #                                 data = dfRes_pos_2_or_3,family = "binomial")
# 
# 
# jtools::summ(mod_glm_severe_2_3,exp=T)
# # jtools::summ(mod_glm_severe_calendar_time,exp=T)
# # jtools::summ(mod_glm_a_lot_or_a_little,exp=T)
# 
# # output severe table
# tab_severity_2_3=jtools::summ(mod_glm_severe_2_3,exp=T)$coeftable %>% as.data.frame() %>% 
#   select(-`z val.`)
# tab_severity_2_3$OR_with_CI=paste0(round(tab_severity_2_3$`exp(Est.)`,2), ", [", round(tab_severity_2_3$`2.5%`,2),",",
#                                round(tab_severity_2_3$`97.5%`,2),"]")
# tab_severity_2_3$variable = rownames(tab_severity_2_3)
# tab_severity_2_3 <- tab_severity_2_3 %>% select(variable,OR_with_CI,p)
# colnames(tab_severity_2_3)= c("Independent variable","Odds ratio with CI","P-value")
# # tab_severity_2_3 <- tab_severity_2_3 %>% select(variable, everything())
# tab_severity_2_3$`Independent variable` <- c("Intercept",
#                                              "Age group: 25-34",
#                                          "Age group: 35-44",
#                                          "Age group: 45-54",
#                                          "Age group: 55-64",
#                                          "Age group: 65-74",
#                                          "Age group: 74+",
#                                          "Sex: male",
#                                          # "Prior COVID-19: Yes",
#                                          "BA.2",
#                                          "Boosted: yes",
#                                          "Days since last vaccine",
#                                          "[Days since last vaccine]^2",
#                                          "Round 18",
#                                          "Round 19"
# )
# 
# tab_severity_2_3$pval <- OverReact::pValEpiConverter(tab_severity_2_3$`P-value`)
# 
# 
# 
# 
# # Loop over all symptoms --------------------------------------------------
# 
# 
# sympFunc <- function(symp,varnum=9,dat,
#                      vars = c("age_group_named","sex", "BA2", "prior_covid", 
#                               "days_since_boost", "round")){
#   f=as.formula(paste0(symp," ~ ",paste0(vars,collapse="+")))
#   mod = glm(formula = f,data = dat,family = "binomial")
#   tab=jtools::summ(mod,exp=T)$coeftable %>% as.data.frame()
#   tab$symptom = sympnames_type_df$symptom[sympnames_type_df$symptom_code==symp]
#   return(tab[varnum,])
# }
# 
# # function to apply function above and clean df
# sympFuncTabGenerator <- function(varnum = 9,dat=dfRes_pos_boost,
#                                  vars = c("age_group_named","sex", "BA2", 
#                                           "prior_covid", "days_since_boost", "round")){
#   allmods=lapply(sympnames_type_df$symptom_code,sympFunc,varnum = varnum,dat=dat,vars=vars)
#   names(allmods)=sympnames_type_df$symptom
#   allmods=bind_rows(allmods) %>% select(-`z val.`)
#   colnames(allmods)= c("Odds ratio","Lower (95% CI)","Upper (95% CI)","P-value","Symptom")
#   allmods <- allmods %>% select(Symptom, everything())
#   allmods$`P-value` <- OverReact::pValEpiConverter(allmods$`P-value`)
#   return(allmods)
# }
# 
# # run over all symptoms
# allmods_ba2 <- sympFuncTabGenerator(9)
# allmods_priorinfect <- sympFuncTabGenerator(10)
# 
# savePrettyExcelWorkbook(listOfTables = list(tab_a_lot=tab_severity, 
#                                             tab_a_lot_plus_time=tab_severity_calendar_time,
#                                             tab_severity_2_3=tab_severity_2_3,
#                                             ors_ba2=allmods_ba2,
#                                             ors_priorinfection=allmods_priorinfect),
#                         workbookName = "severity_modelling",outpath = outpath)
# 
# 
# 
# # Model among all positives -----------------------------------------------
# 
# allmods_priorinfect_allpos <- sympFuncTabGenerator(varnum = 10,dat=dfRes %>% filter(estbinres==1,
#                                                                                    round %in%17:19),
#                                             vars = c("age_group_named","sex", "BA2", 
#                                                      "prior_covid",  "round",
#                                                      "boosted")
#                                             )
# allmods_boosted_allpos <- sympFuncTabGenerator(varnum = 9,dat=dfRes %>% filter(estbinres==1,
#                                                                                    round%in% 17:19),
#                                                    vars = c("age_group_named","sex",  "boosted",
#                                                             "prior_covid","round"
#                                                    ))
# 
# 
# dat$prior_covid %>% table(dat$symptomatic_26) %>% prop.table(1)
# 
# 
# 
# 
# # Analyse time since prior infection --------------------------------------
# 
# ### This is the rationale for not using prior infection as a raw variable in the models -
# ### especially in the later rounds, it captures a lot of people who are still testing positive from
# ### an infection in the past 28 days (we assume)
# 
# dfRes %>% 
#   filter(!infection_to_swab>600,!infection_to_swab<0) %>%
#   ggplot(aes(x=infection_to_swab, fill=round)) +
#   geom_histogram(bins = 100) + 
#   facet_wrap(.~round) +
#   OverReact::theme_react(strip_text_size = 9)
# 
# 
# dfRes %>% 
#   filter(!infection_to_swab>600,!infection_to_swab<0, estbinres==1) %>%
#   ggplot(aes(x=infection_to_swab, fill=round)) +
#   geom_histogram(bins = 100) + 
#   facet_wrap(.~round) +
#   OverReact::theme_react(strip_text_size = 9)



# Incremental model maker -------------------------------------------------

dfRes_pos_2_or_3$days_since_symptom_onset <- dfRes_pos_2_or_3$days_since_symptom_onset/7
dfRes_pos_2_or_3$calendar_time <- dfRes_pos_2_or_3$calendar_time/7
dfRes_pos_2_or_3$time_since_last_vax <- dfRes_pos_2_or_3$time_since_last_vax/7





myvars=c("BA2_bin","age_group_named","sex",
         "boosted","time_since_last_vax", "prior_covid_28plus","days_since_symptom_onset",
         "calendar_time")



cov_name_list_2 <- c("Omicron variant","Age","Sex","Boosted (Yes)","Weeks since last vaccination",
                            "Prior COVID-19 (28+ days ago)","Weeks since symptom onset",
                     "Calendar time"
                     )

dfRes_pos_2_or_3$days_since_symptom_onset %>% table(exclude="none")
names(cov_name_list_2) <- myvars
dfRes_pos_2_or_3$sex <- as.factor(dfRes_pos_2_or_3$sex)
class(dfRes_pos_2_or_3$BA2)
dfRes_pos_2_or_3$BA2_bin <- factor(case_when(dfRes_pos_2_or_3$BA2==1 ~ "BA.2",
                                      T~"BA.1"),levels=c("BA.1","BA.2"))
mod_inc=OverReact::ModelMakerMulti(dat = dfRes_pos_2_or_3,add_askerisks_for_pvals = T,
                                   list_of_variables_of_interest = myvars,outcome="impact_a_lot",
                                   sf=2,simpleround =T,joint_adjustment_vars=myvars,
                                   cov_name_list=cov_name_list_2,family = "binomial"
                                   )
mod_inc_df=mod_inc$df_output
mod_inc_df <- mod_inc_df[!grepl("ntercept",mod_inc_df$Category),]
mod_inc_df <- mod_inc_df %>% select(-plus_BA2_bin)
mod_inc_df$crude_mod_OR[3:nrow(mod_inc_df)] <- ""
mod_inc_df$plus_age_group_named[9:nrow(mod_inc_df)] <- ""
mod_inc_df$plus_sex[11:nrow(mod_inc_df)] <- ""
mod_inc_df$plus_boosted[13:nrow(mod_inc_df)] <- ""
mod_inc_df$plus_time_since_last_vax[14:nrow(mod_inc_df)] <- ""
mod_inc_df$plus_prior_covid_28plus[16:nrow(mod_inc_df)] <- ""
mod_inc_df$plus_days_since_symptom_onset[17:nrow(mod_inc_df)] <- ""
# mod_inc_df$plus_round[20:nrow(mod_inc_df)] <- ""
mod_inc_df[mod_inc_df=="NA (NA,NA) "] <- "-"
mod_inc_df[mod_inc_df=="NA (NA,NA)"] <- "-"


### With all vaccine numbers ###
dfRes_pos_anyvax=dfRes %>% 
  filter(variant_inferred_detail%in% c("BA.1 (Omicron)","BA.2 (Omicron)"), 
                                  round %in%c(17:19), !covidability7%in% c(4:5)) %>% 
  mutate(boosted = factor(case_when(vaccdose==2~ "No",
                                    T ~ "Yes"), levels = c("No","Yes")),
         vax_count=factor(case_when(vaccdose==0 ~ "Unvaccinated/one vaccine",
                             vaccdose==1 ~ "Unvaccinated/one vaccine",
                             vaccdose==2 ~ "Two vaccines",
                             vaccdose==3 ~ "Three vaccines",
                             T~ NA_character_), levels=c(
                               "Unvaccinated/one vaccine","Two vaccines","Three vaccines"
                             )),
         time_since_last_vax= as.numeric(case_when(!is.na(days_since_boost) ~ days_since_boost,
                                                   !is.na(days_since_second_vax) ~ days_since_second_vax,
                                                   # T ~ NA_real_
         )))
dfRes_pos_anyvax$BA2_bin <- factor(case_when(dfRes_pos_anyvax$BA2==1 ~ "BA.2",
                                             T~"BA.1"),levels=c("BA.1","BA.2"))
dfRes_pos_anyvax$days_since_symptom_onset <- dfRes_pos_anyvax$days_since_symptom_onset/7
dfRes_pos_anyvax$calendar_time <- dfRes_pos_anyvax$calendar_time/7
dfRes_pos_anyvax$time_since_last_vax <- dfRes_pos_anyvax$time_since_last_vax/7

myvars2=c("BA2_bin","age_group_named","sex",
         "vax_count","prior_covid_28plus","days_since_symptom_onset","calendar_time")
# cov_name_list_2=as.list(myvars)
cov_name_list_3 <- c("Omicron variant","Age","Sex","Number of vaccines",
                     "Prior COVID-19 (28+ days ago)","Weeks since symptom onset","Calendar time"
)

names(cov_name_list_3) <- myvars2
dfRes_pos_anyvax$sex <- as.factor(dfRes_pos_anyvax$sex)

mod_inc_2=OverReact::ModelMakerMulti(dat = dfRes_pos_anyvax,add_askerisks_for_pvals = T,
                                   list_of_variables_of_interest = myvars2,outcome="impact_a_lot",
                                   sf=2,simpleround =T,joint_adjustment_vars=myvars2,
                                   cov_name_list=cov_name_list_3
)
table(dfRes_pos_anyvax$vax_count)

mod_inc_2_df=mod_inc_2$df_output
mod_inc_2_df <- mod_inc_2_df[!grepl("ntercept",mod_inc_2_df$Category),]
mod_inc_2_df <- mod_inc_2_df %>% select(-plus_BA2_bin)
mod_inc_2_df$crude_mod_OR[3:nrow(mod_inc_2_df)] <- ""
mod_inc_2_df$plus_age_group_named[10:nrow(mod_inc_2_df)] <- ""
mod_inc_2_df$plus_sex[12:nrow(mod_inc_2_df)] <- ""
mod_inc_2_df$plus_vax_count[15:nrow(mod_inc_2_df)] <- ""
mod_inc_2_df$plus_prior_covid_28plus[17:nrow(mod_inc_2_df)] <- ""
mod_inc_2_df$plus_days_since_symptom_onset[18:nrow(mod_inc_2_df)] <- ""
# mod_inc_2_df$plus_round[21:nrow(mod_inc_2_df)] <- ""
mod_inc_2_df[mod_inc_2_df=="NA (NA,NA) "] <- "-"
mod_inc_2_df[mod_inc_2_df=="NA (NA,NA)"] <- "-"


# Symptom count as outcome ------------------------------------------------


mod_inc_sxcount=OverReact::ModelMakerMulti(dat = dfRes_pos_2_or_3,add_askerisks_for_pvals = T,
                                   list_of_variables_of_interest = myvars,outcome="symptom_count_26",
                                   sf=2,simpleround =T,joint_adjustment_vars=myvars,
                                   cov_name_list=cov_name_list_2,family = "poisson"
)
mod_inc_sxcount_df=mod_inc_sxcount$df_output
mod_inc_sxcount_df <- mod_inc_sxcount_df[!grepl("ntercept",mod_inc_sxcount_df$Category),]
mod_inc_sxcount_df <- mod_inc_sxcount_df %>% select(-plus_BA2_bin)
mod_inc_sxcount_df$crude_mod_OR[3:nrow(mod_inc_sxcount_df)] <- ""
mod_inc_sxcount_df$plus_age_group_named[9:nrow(mod_inc_sxcount_df)] <- ""
mod_inc_sxcount_df$plus_sex[11:nrow(mod_inc_sxcount_df)] <- ""
mod_inc_sxcount_df$plus_boosted[13:nrow(mod_inc_sxcount_df)] <- ""
mod_inc_sxcount_df$plus_time_since_last_vax[14:nrow(mod_inc_sxcount_df)] <- ""
mod_inc_sxcount_df$plus_prior_covid_28plus[16:nrow(mod_inc_sxcount_df)] <- ""
mod_inc_sxcount_df$plus_days_since_symptom_onset[17:nrow(mod_inc_sxcount_df)] <- ""
# mod_inc_sxcount_df$plus_round[20:nrow(mod_inc_sxcount_df)] <- ""
mod_inc_sxcount_df[mod_inc_sxcount_df=="NA (NA,NA) "] <- "-"
mod_inc_sxcount_df[mod_inc_sxcount_df=="NA (NA,NA)"] <- "-"




### SAVE ###

savePrettyExcelWorkbook(listOfTables = list(incremental_models=mod_inc_df,
                                            incremental_models_allvax=mod_inc_2_df,
                                            incremental_models_sxcount=mod_inc_sxcount_df),
                        workbookName = "severity_mods_incremental",outpath = outpath)
