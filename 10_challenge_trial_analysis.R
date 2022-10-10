
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
                  "patchwork","randomcoloR","focus","car"
)

load_packages(package.list)


# Import REACT-1 data -----------------------------------------------------
source("E:/Group/react2_study5/report_phases_combined/projects/omicron_symptom_profiling/code/00_bits_and_pieces.R", local = T)
source("E:/Group/react2_study5/report_phases_combined/projects/omicron_symptom_profiling//code/00_data_prep.R", local = T)

# create subfolder
createMySubfolder(subfolderName = "challenge_trial_analysis")


# filter to ages 18-30 and not shielding
dfRes <- dfRes %>%   
  filter(round %in% c(2:7,8:10,13:15,16,17,18,19)) %>% 
  filter(age <=30, age >=18,shield2!=1) %>% 
  mutate(round = as.numeric(round),
         wave = factor(case_when(round <8 ~ "Rounds 2-7 (Wild type)",
                                 round <11 ~ "Rounds 8-10 (Alpha)",
                                 round %in% c(13:15) ~ "Rounds 13-15 (Delta)",
                                 round%in% c(17) ~ "Round 17 (BA.1 Omicron)",
                                 round%in% c(19) ~ "Rounds 19 (BA.2 Omicron)"
         ),
         levels = c("Rounds 2-7 (Wild type)",
                    "Rounds 8-10 (Alpha)",
                    "Rounds 13-15 (Delta)",
                    "Round 17 (BA.1 Omicron)",
                    "Rounds 19 (BA.2 Omicron)")),
         wave_omi = factor(case_when(round <8 ~ "Rounds 2-7 (Wild type)",
                                     round <11 ~ "Rounds 8-10 (Alpha)",
                                     round %in% c(13:15) ~ "Rounds 13-15 (Delta)",
                                     round%in% c(17:19) ~ "Rounds 17-19 (Omicron)",
         ),
         levels = c("Rounds 2-7 (Wild type)",
                    "Rounds 8-10 (Alpha)",
                    "Rounds 13-15 (Delta)",
                    "Rounds 17-19 (Omicron)")),
         
         brediff_cat = case_when(brediff==2 ~ "Yes, Even when sitting/lying down",
                                 brediff==1 ~ "Yes, affecting normal activities",
                                 brediff==3 ~ "No, affecting normal activities",
                                 T ~ NA_character_
         ),
         covidability_cat = factor(case_when(symptomatic_26==0 ~ "No symptoms reported",
                                             covidability7 == 1 ~ "A lot",
                                             covidability7 == 2 ~ "A little",
                                             covidability7 == 3 ~ "Not at all",
                                             covidability7 %in% c(4,5) ~ "Don't know / PNA / Non-reponse",
                                             T ~ "Don't know / PNA / Non-reponse",
         ), 
         levels = c("A lot","A little","Not at all",
                    "Don't know / PNA / Non-reponse","No symptoms reported")))

dfRes$one_of_four <- as.factor(dfRes$one_of_four)
dfRes$symptomatic <- as.factor(dfRes$symptomatic)
dfRes$estbinres_char <- ifelse(dfRes$estbinres==1,"Yes","No")
dfRes$round <- as.character(dfRes$round)
dfRes$ethnic_new %>% table(dfRes$round)
table(dfRes$round,dfRes$covida_cat)
cov_name_list$prior_inf = "Prior COVID-19 infection"
cov_name_list$symptom_count_26 = "Number of reported symptoms"
cov_name_list$symptom_count_4 = "Number of reported symptoms (classic symptoms only)"
cov_name_list$vax_status_cat = "Vaccination status"
cov_name_list$covidability_cat="Symptoms affecting day-to-day activities"
cov_name_list$seekmed_cat="Sought medical attention for symptoms"
cov_name_list$brediff_cat="Breathing difficulties"

# define list of rowvars
rowvar_list=c("all_participants", "sex","age_group_named","ethnic_new","prior_inf",
              "vax_status_cat",
              "symptomatic", "symptom_count_26","seekmed_cat","covidability_cat",
              "days_since_symptom_onset",covid_yesnos)

dfRes$vax_status_cat[is.na(dfRes$vax_status_cat)] <- "NA"
dfRes$vax_status_cat %>% table(exclude="none")
dfRes$symptomatic_26 = as.numeric(rowSums(dfRes[,covid_yesnos],na.rm = T)>0)
dfRes$symptom_count_26 <- rowSums(dfRes[,covid_yesnos], na.rm=T)


# Plot symptom prevalence per wave ----------------------------------------

createTableOne <- function(dat,colvar = "wave_omi"){
  
  # Run on only positives
  table_one_pos <- OverReact::crossTabMulti(dat = dat,rowvar_list =  
                                              rowvar_list,
                                            colvar = colvar,
                                            cov_names = cov_name_list,confint = T,include_percentages = T,
                                            rowwise_precentages = F,statistical_test = F) 
  
  # table_one_pos_plot <- OverReact::makeXtabPlottable(myxtab = table_one_pos_plot)
  
  # Tinker
  table_one_pos <- table_one_pos %>% filter(Category!=0)
  switchindex <- table_one_pos$Category==1
  table_one_pos$Category[switchindex] <- table_one_pos$Variable[switchindex]
  table_one_pos$Variable[switchindex] <- "Symptoms"
  ### Replace NaNs
  table_one_pos[table_one_pos=="0 (NaN%, [0-NaN])"] <- "0"
  return(table_one_pos)
}

### run function over subgroups
table_one_all <- createTableOne(dat = dfRes)
table_one_vax <- createTableOne(dat = dfRes%>% filter(vax_status_cat %in% c(
                                                         "Three does","Two does")))
table_one_unvax <- createTableOne(dat = dfRes%>% filter((vax_status_cat %in% c(
                                                           "NA", "Not vaccinated", "One does"))))

table_one_pos <- createTableOne(dat = dfRes %>% filter( estbinres==1))
table_one_pos_vax <- createTableOne(dat = dfRes %>% filter( estbinres==1,
                                                            vax_status_cat %in% c(
                                                              "Three does","Two does")))
table_one_pos_unvax <- createTableOne(dat = dfRes %>% filter( estbinres==1,
                                                            (vax_status_cat %in% c(
                                                              "NA", "Not vaccinated", "One does"))))




# UNivariate analysis -----------------------------------------------------

dfRes <- dfRes %>% 
  mutate_at(covid_yesnos_firstsymp,binaryCleaner_1_0) %>% 
  mutate_at(covid_yesnos_month,binaryCleaner_1_0)  %>% 
  mutate(doublevaxxed_or_boosted=case_when(vax_status_number>=2 ~1,
                                           TRUE ~0),
         sex=factor(sex)) %>% 
  mutate_at(c(sympnames_type_df$symptom_code,covid_yesnos_firstsymp,covid_yesnos_month,
              "sex","doublevaxxed_or_boosted"),as.factor) 



# Run models --------------------------------------------------------------



runUnivariate <- function(dat, variant="Omicron",rounds =c(17), adj_level =4,
                          joint_adjustment_vars = c("age", "sex")){
  
  # OMICRON #
  
  uiv_omicron =OverReact::ModelMakerMulti(dat = dat %>% filter(round %in% rounds),
                                          list_of_variables_of_interest = sympnames_type_df$symptom_code,
                                          outcome = "estbinres",
                                          cov_name_list = cov_name_list,
                                          joint_adjustment_vars = joint_adjustment_vars)
  df=uiv_omicron$df_output
  df <- df[df$Category!="0 [reference]",]
  df[df$Category!="0 [reference]",]
  df <- df[grepl("change",df$Variable),]
  df <- df %>% select(-Category)
  df$Variant = variant
  
  return(df)
  
}

# function for all variants
runMultipleUnivariate <- function(dat,runwild=T){
  ### Delta
  univ_delta=runUnivariate(dat, variant="Delta",rounds =c(12:15))
  
  
  if(runwild){
    ### Alpha
    univ_alpha=runUnivariate(dat, variant="Alpha",rounds =c(8:10))
    
    ### WIld type
    univ_wildtype=runUnivariate(dat, variant="Wildtype",rounds =c(2:7))
    # bind
    res=rbind(univ_wildtype,univ_alpha,univ_delta)
    
  }else{
    res=univ_delta
    
  }
  
 
  return(res)
  
  }


univ_res_all=runMultipleUnivariate(dfRes)
univ_res_unvax=runMultipleUnivariate(dfRes%>% filter((vax_status_cat %in% c(
  "NA", "Not vaccinated", "One does"))))
univ_res_vax=runMultipleUnivariate(dat = dfRes%>% filter(vax_status_cat %in% c(
                                                      "Three does","Two does")),
                                   runwild = F)



# Export ------------------------------------------------------------------

savePrettyExcelWorkbook(listOfTables = list(table_one_all =table_one_all,
                                            table_one_vax = table_one_vax,
                                            table_one_unvax=table_one_unvax,
                                            table_one_pos = table_one_pos,
                                            table_one_pos_vax=table_one_pos_vax,
                                            table_one_pos_unvax=table_one_pos_unvax,
                                            univ_res_all=univ_res_all,
                                            univ_res_vax=univ_res_vax,
                                            univ_res_unvax=univ_res_unvax),
                        workbookName = "react1_results_for_challenge_trial",outpath = outpath)


