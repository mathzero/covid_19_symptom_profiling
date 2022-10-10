
#' First clear the environment of variables
rm(list=ls(all=TRUE))
# get root director of project
root.dir <- getwd()
# setwd(dir = "/Group/react2_study5/report_phases_combined/projects/omicron_symptom_profiling/")
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




#' Pull in packages needed
package.list <- c("prevalence","mgcv","knitr","MASS","kableExtra","table1","dplyr",
                  "tidyr", "pheatmap","OverReact",
                  "ggplot2","ggsci", "RColorBrewer", "tidyverse", "lubridate", 
                  "readr","ggthemes", "questionr", "foreach", "doParallel","withr",
                  "patchwork","randomcoloR","focus"
)

load_packages(package.list)


# Import REACT-1 data -----------------------------------------------------
source("E:/Group/react2_study5/report_phases_combined/projects/omicron_symptom_profiling/code/00_bits_and_pieces.R",
       local = T)
source("E:/Group/react2_study5/report_phases_combined/projects/omicron_symptom_profiling/code/00_data_prep.R", 
       local = T)


# create subfolder
createMySubfolder(subfolderName = "eda")

# Table one ---------------------------------------------------------------
dfRes$symptomatic_26

dfRes_tab <- dfRes %>% 
  filter(round %in% c(2:7,8:10,13:15,16,17,18,19)) %>% 
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

# create symptom count variable
dfRes_tab$symptom_count_26 <- rowSums(dfRes_tab[,covid_yesnos], na.rm=T)
dfRes_tab$symptom_count_4 <- rowSums(dfRes_tab[,sympnames_type_df$symptom_code[1:4]], na.rm=T)

dfRes_tab$seekmed_cat %>% table(dfRes_tab$wave, exclude="none") %>% prop.table() *100
dfRes_tab$covidability7[dfRes_tab$estbinres==1] %>% 
  table(dfRes_tab$symptomatic_26[dfRes_tab$estbinres==1], exclude="none") %>% 
  prop.table(2) *100



dfRes_tab$one_of_four <- as.factor(dfRes_tab$one_of_four)
dfRes_tab$symptomatic <- as.factor(dfRes_tab$symptomatic)
dfRes_tab$estbinres_char <- ifelse(dfRes_tab$estbinres==1,"Yes","No")
dfRes_tab$round <- as.character(dfRes_tab$round)
dfRes_tab$ethnic_new %>% table(dfRes_tab$round)
table(dfRes$round,dfRes$covida_cat)
cov_name_list$prior_inf = "Prior COVID-19 infection"
cov_name_list$symptom_count_26 = "Number of reported symptoms"
cov_name_list$symptom_count_4 = "Number of reported symptoms (classic symptoms only)"
cov_name_list$vax_status_noDate_v2 = "Vaccination status"
cov_name_list$covidability_cat="Symptoms affecting day-to-day activities"
cov_name_list$seekmed_cat="Sought medical attention for symptoms"
cov_name_list$brediff_cat="Breathing difficulties"

# define list of rowvars
rowvar_list=c("all_participants", "sex","age_group_named","ethnic_new","estbinres_char","prior_inf",
              "vax_status_noDate_v2",
              "symptomatic", "symptom_count_26","seekmed_cat","covidability_cat",
              "days_since_symptom_onset",covid_yesnos,sympnames_type_df$first_symptom_code[1:26])

# Run on all participants
table_one <- OverReact::crossTabMulti(dat = dfRes_tab,rowvar_list =  rowvar_list,
                                      colvar = "wave",cov_names = cov_name_list,
                                      confint = T,include_percentages = T,comma_thousands = T,
                                      rowwise_precentages = F,statistical_test = F)





# Tinker
table_one <- table_one %>% filter(Category!=0)
switchindex <- table_one$Category==1
table_one$Category[switchindex] <- table_one$Variable[switchindex]
table_one$Variable[switchindex] <- "Symptoms"

### Replace NaNs
table_one[table_one=="0 (NaN%)"] <- "0"
table_one$Variable[(nrow(table_one)-25):nrow(table_one)] <- "First reported symptoms"



# Table one just omicron --------------------------------------------------

# Run on all participants
table_one_omicron <- OverReact::crossTabMulti(dat = dfRes_tab,rowvar_list =  rowvar_list,
                                      colvar = "wave_omi",cov_names = cov_name_list,
                                      confint = T,include_percentages = T,comma_thousands = T,
                                      rowwise_precentages = F,statistical_test = F)
# Tinker
table_one_omicron <- table_one_omicron %>% filter(Category!=0)
switchindex <- table_one_omicron$Category==1
table_one_omicron$Category[switchindex] <- table_one_omicron$Variable[switchindex]
table_one_omicron$Variable[switchindex] <- "Symptoms"

### Replace NaNs
table_one_omicron[table_one_omicron=="0 (NaN%)"] <- "0"
table_one_omicron$Variable[(nrow(table_one_omicron)-25):nrow(table_one_omicron)] <- "First reported symptoms"




# Run on positives --------------------------------------------------------

dfRes %>% filter(estbinres==1) %>% group_by(variant_inferred_detail) %>% 
  summarise(mean_days=mean(days_since_symptom_onset,na.rm=T))

# Run on only positives
table_one_pos <- OverReact::crossTabMulti(dat = dfRes_tab %>% filter( estbinres==1),rowvar_list =  rowvar_list,
                                      colvar = "variant_inferred_detail",cov_names = cov_name_list,
                                      confint = T,include_percentages = T,comma_thousands = T,
                                      rowwise_precentages = F,statistical_test = F)
# Tinker
table_one_pos <- table_one_pos %>% filter(Category!=0)
switchindex <- table_one_pos$Category==1
table_one_pos$Category[switchindex] <- table_one_pos$Variable[switchindex]
table_one_pos$Variable[switchindex] <- "Symptoms"
### Replace NaNs
table_one_pos[table_one_pos=="0 (NaN%)"] <- "0"

# Add first symptom name
table_one_pos$Variable[(nrow(table_one_pos)-25):nrow(table_one_pos)] <- "First reported symptoms"

# Run on only positives
table_one_pos_60_plus <- OverReact::crossTabMulti(dat = dfRes_tab %>% filter( estbinres==1, age>=60),
                                                  rowvar_list =  rowvar_list,
                                          colvar = "variant_inferred_detail",cov_names = cov_name_list,
                                          confint = F,include_percentages = T,
                                          rowwise_precentages = F,statistical_test = F)
# Tinker
table_one_pos_60_plus <- table_one_pos_60_plus %>% filter(Category!=0)
switchindex <- table_one_pos_60_plus$Category==1
table_one_pos_60_plus$Category[switchindex] <- table_one_pos_60_plus$Variable[switchindex]
table_one_pos_60_plus$Variable[switchindex] <- "Symptoms"
### Replace NaNs
table_one_pos_60_plus[table_one_pos_60_plus=="0 (NaN%)"] <- "0"



# FIrst symptoms ----------------------------------------------------------

# add first symptoms to covnamelist
cov_name_list[sympnames_type_df$first_symptom_code[1:26]] <- sympnames_type_df$symptom[1:26]

# Sort out symptoms
dfRes_tab <- dfRes_tab %>% 
  mutate_at(all_of(sympnames_type_df$first_symptom_code[1:26]),binaryCleaner_1_0) 


### Replace NA with 0 in symptoms
dfRes_tab[sympnames_type_df$first_symptom_code[1:26]][is.na(dfRes_tab[sympnames_type_df$first_symptom_code[1:26]])] <- 0

# Run on only positives
table_one_pos_first_symp <- OverReact::crossTabMulti(dat = dfRes_tab %>% 
                                                       filter( estbinres==1, symptomatic==1),
                                                     rowvar_list =  
                                                       sympnames_type_df$first_symptom_code[1:26],
                                          colvar = "variant_inferred_detail",
                                          cov_names = cov_name_list,
                                          confint = F,
                                          include_percentages = T,
                                          rowwise_precentages = F,statistical_test = F) 

table_one_pos_first_symp <- table_one_pos_first_symp %>% filter(Category!=0)


### save workbooks
savePrettyExcelWorkbook(listOfTables = list(tab1_all=table_one, table_one_omicron=table_one_omicron,
                                            tab1_pos=table_one_pos,
                                            tab1_pos_60_plus=table_one_pos_60_plus,
                                            tab1_first_symptom=table_one_pos_first_symp),
                        workbookName = "table_one",outpath = outpath)




# Plot symptom prevalence per wave ----------------------------------------

dfRes_tab$symptomatic_26 = as.numeric(rowSums(dfRes_tab[,covid_yesnos],na.rm = T)>0)

# Run on only positives
table_one_pos <- OverReact::crossTabMulti(dat = dfRes_tab %>% filter( estbinres==1),rowvar_list =  
                                                 c("symptomatic_26", covid_yesnos),
                                          colvar = "variant_inferred_detail",
                                          cov_names = cov_name_list,confint = T,include_percentages = T,
                                          rowwise_precentages = F,statistical_test = F) 

# table_one_pos_plot <- OverReact::makeXtabPlottable(myxtab = table_one_pos_plot)

table_one_pos$`Sum / mean(SD)`


table_one_pos <- table_one_pos%>% 
  filter(Category!=0) %>% 
  select(-Category, -`Sum / mean(SD)`) %>% 
  pivot_longer(cols = -Variable)

table_one_pos_plot <- table_one_pos

# convert % to numeric
table_one_pos_plot$lower <- as.numeric(OverReact::xtabPercentageExtractor(mystring = 
                                                                                 table_one_pos_plot$value,
                                                                               lookbehind = "\\[",
                                                                          lookahead = "[-]"))
table_one_pos_plot$upper <- as.numeric(OverReact::xtabPercentageExtractor(mystring = 
                                                                            table_one_pos_plot$value,
                                                                          lookbehind = "[-]",
                                                                          lookahead = "\\]"))

table_one_pos_plot$percentage <- as.numeric(OverReact::xtabPercentageExtractor(mystring = 
                                                                                 table_one_pos_plot$value,
                                                                               lookbehind = "\\(",
                                                                               lookahead = "\\%"))



dodge_mult=0.8
dodger <- position_dodge2(width = dodge_mult,reverse = F)

stripes_df <- table_one_pos_plot %>% arrange(-percentage) %>% 
  filter(name=="Rounds 2-7 (Wild type)") %>% 
  left_join(sympnames_type_df, by = c("Variable"= "symptom"))%>% 
  group_by(symptom_type) %>% 
  mutate(stripe_col = case_when(row_number() %%2 ==0 ~ "NA", TRUE ~ "grey90")
  )


### Get summary df to fix column widths
summ <- table_one_pos_plot %>% 
  left_join(sympnames_type_df, by = c("Variable"= "symptom")) %>% 
  group_by(symptom_type) %>% 
  summarise(n=n())

# join
sympnames_type_df <- sympnames_type_df %>% left_join(summ)
# set col width
sympnames_type_df$col_width=dodge_mult*sympnames_type_df$n/max(sympnames_type_df$n)

# PLot
p_prevs=table_one_pos_plot %>% 
  left_join(sympnames_type_df, by = c("Variable"= "symptom")) %>% 
  mutate(name = factor(name, levels = unique(name)),
         symptom_type = factor(symptom_type, levels = c("Overall",  "Smell/taste","Respiratory", 
                                                "Coryzal", "Gastrointestinal","Fatigue", "Other"))) %>% 
  ggplot() +
  # geom_tile(data = stripes_df , 
  #           aes(x=reorder(Variable,percentage), y =1, height = Inf,
  #                                  fill = stripe_col),
  #           alpha=0.8,
  #           col = "grey70",
  #           linetype = "dashed",
  #           show.legend = F) +
  # scale_fill_manual(values = c("white","grey96")) +
  # ggnewscale::new_scale_fill() +
  geom_bar(aes(x=reorder(Variable,percentage), y= percentage, fill = name, width=col_width),
           stat = "identity",
           position = dodger, col = "black", size = 0.01) +
  geom_errorbar(position = dodger,
                aes(x=reorder(Variable,percentage),ymin= lower, ymax=upper, width=col_width),
                size=0.3) +
  theme_react(strip_text_size = 10) +
  scale_x_discrete(labels = function(x) str_wrap(x, width =17)) +
  scale_y_continuous(breaks = scales::breaks_width(10)) +
  # scale_fill_brewer(palette = "Reds") +
  OverReact::scale_fill_imperial(palette = "default") +
  # coord_flip() +
  # facet_grid(scales = "free",rows = "symptom_type",space = "fixed",shrink = T) +
  ggforce::facet_col(facets = "symptom_type",scales = "free",space = "free") +
  labs(x="", y="% of PCR positive respondents with symptom in past week", fill = "") +
  theme(legend.position = "bottom",
        panel.spacing = unit(1, "cm"))
p_prevs


# save
OverReact::saveREACTplot(p = p_prevs,figpath = figpath,filename = "symptom_prevalence_by_variant",
                         width = 7.5,height = 10)




# As above for negatives --------------------------------------------------


# Run on only negitives
table_one_neg_plot <- OverReact::crossTabMulti(dat = dfRes_tab %>% filter( estbinres==0),rowvar_list =  
                                                 c("symptomatic_26", covid_yesnos),
                                               colvar = "wave",
                                               cov_names = cov_name_list,confint = T,include_percentages = T,
                                               rowwise_precentages = F,statistical_test = F) 

# table_one_neg_plot <- OverReact::makeXtabPlottable(myxtab = table_one_neg_plot)

table_one_neg_plot


table_one_neg_plot <- table_one_neg_plot%>% 
  filter(Category!=0) %>% 
  select(-Category, -`Sum / mean(SD)`) %>% 
  pivot_longer(cols = -Variable)

# convert % to numeric
table_one_neg_plot$lower <- as.numeric(OverReact::xtabPercentageExtractor(mystring = 
                                                                            table_one_neg_plot$value,
                                                                          lookbehind = "\\[",
                                                                          lookahead = "[-]"))
table_one_neg_plot$upper <- as.numeric(OverReact::xtabPercentageExtractor(mystring = 
                                                                            table_one_neg_plot$value,
                                                                          lookbehind = "[-]",
                                                                          lookahead = "\\]"))

table_one_neg_plot$percentage <- as.numeric(OverReact::xtabPercentageExtractor(mystring = 
                                                                                 table_one_neg_plot$value,
                                                                               lookbehind = "\\(",
                                                                               lookahead = "\\%"))



dodge_mult=0.8
dodger <- position_dodge2(width = dodge_mult,reverse = F)

stripes_df <- table_one_neg_plot %>% arrange(-percentage) %>% 
  filter(name=="Rounds 2-7 (Wild type)") %>% 
  left_join(sympnames_type_df, by = c("Variable"= "symptom"))%>% 
  group_by(symptom_type) %>% 
  mutate(stripe_col = case_when(row_number() %%2 ==0 ~ "NA", TRUE ~ "grey90")
  )


### Get summary df to fix column widths
summ <- table_one_neg_plot %>% 
  left_join(sympnames_type_df, by = c("Variable"= "symptom")) %>% 
  group_by(symptom_type) %>% 
  summarise(n=n())

# join
sympnames_type_df <- sympnames_type_df %>% left_join(summ)
# set col width
sympnames_type_df$col_width=dodge_mult*sympnames_type_df$n/max(sympnames_type_df$n)

# PLot
p_prevs_neg=table_one_neg_plot %>% 
  left_join(sympnames_type_df, by = c("Variable"= "symptom")) %>% 
  mutate(name = factor(name, levels = unique(name)),
         symptom_type = factor(symptom_type, levels = c("Overall",  "Smell/taste","Respiratory", 
                                                        "Coryzal", "Gastrointestinal","Fatigue", "Other"))) %>% 
  ggplot() +
  # geom_tile(data = stripes_df , 
  #           aes(x=reorder(Variable,percentage), y =1, height = Inf,
  #                                  fill = stripe_col),
  #           alpha=0.8,
  #           col = "grey70",
  #           linetype = "dashed",
  #           show.legend = F) +
  # scale_fill_manual(values = c("white","grey96")) +
  # ggnewscale::new_scale_fill() +
  geom_bar(aes(x=reorder(Variable,percentage), y= percentage, fill = name, width=col_width),
           stat = "identity",
           position = dodger, col = "black", size = 0.01) +
  geom_errorbar(position = dodger,
                aes(x=reorder(Variable,percentage),ymin= lower, ymax=upper, width=col_width),
                size=0.3) +
  theme_react(strip_text_size = 10) +
  scale_x_discrete(labels = function(x) str_wrap(x, width =17)) +
  scale_y_continuous(breaks = scales::breaks_width(10)) +
  # scale_fill_brewer(palette = "Reds") +
  OverReact::scale_fill_imperial(palette = "default") +
  # coord_flip() +
  # facet_grid(scales = "free",rows = "symptom_type",space = "fixed",shrink = T) +
  ggforce::facet_col(facets = "symptom_type",scales = "free",space = "free") +
  labs(x="", y="% of PCR negative respondents with symptom in past week", fill = "") +
  theme(legend.position = "bottom",
        panel.spacing = unit(1, "cm"))
p_prevs_neg


# save
OverReact::saveREACTplot(p = p_prevs_neg,figpath = figpath,
                         filename = "symptom_prevalence_by_variant_pcr_neg",
                         width = 7.5,height = 10)




# Plot comparing severity -------------------------------------------------

# get plot dat to save compute
plotdat_severity=dfRes_tab %>% filter(!is.na(covidability_cat),
                                      !is.na(variant_inferred_detail),
                                      estbinres==1)

p_sev <-  plotdat_severity %>% 
  group_by(variant_inferred_detail,covidability_cat) %>% 
  summarise(n=n()) %>% 
  group_by(variant_inferred_detail) %>% 
  mutate(percent=100*n/sum(n)) %>% 
    ggplot(aes(x=variant_inferred_detail, y=percent, fill=covidability_cat)) +
  geom_col(position="dodge") +
  OverReact::theme_react() +
  scale_fill_brewer(palette = "Reds",direction = -1) +
  labs(x="", "%")
p_sev






# Time since symptom onset plot by variant -----------------------------------------------
library(ggbeeswarm)
class(dfRes$variant_inferred)
p_symptom_onset <- dfRes %>% 
  filter(round!=1, estbinres==1, !is.na(days_since_symptom_onset),
         !is.na(variant_inferred_detail)) %>% 
  mutate(variant_inferred=case_when(is.na(variant_inferred) ~ "Transition phase",
                                      round == 19 ~ "Round 19 (Omicron BA.2)",
                                      round == 17 ~ "Round 17 (Omicron BA.1)",
                                      round == 18 ~ "Transition phase",
                                    T ~ as.character(variant_inferred)),
         days_since_symptom_onset_cat=factor(case_when(days_since_symptom_onset==11 ~ "11+",
                                                       T ~ as.character(days_since_symptom_onset)),
                                             levels = c(1:10,"11+"))) %>% 
  ggplot(aes(y=days_since_symptom_onset, x=variant_inferred_detail,col=variant_inferred_detail,
             fill=variant_inferred_detail)) +
  OverReact::scale_color_imperial(palette = "default") +
  OverReact::scale_fill_imperial(palette = "default") +
  geom_violin(alpha=0.5)+
  geom_boxplot(col="black", fill=NA) +
  # geom_point(position = position_jitter(width=0.1, height=0.1), alpha=0.1)+
  # ggbeeswarm::geom_quasirandom(alpha=0.1, groupOnX = F,position = position_jitter(width=10)) +
  OverReact::theme_react() +
  theme(legend.position = "none") +
  scale_y_continuous(breaks=seq(0,11,1)) +
  labs(x="",y="Days since symptom onset",
       title = "B")
  # coord_flip()
  # facet_wrap(.~variant_inferred_detail, nrow=1)
p_symptom_onset

# save
OverReact::saveREACTplot(p = p_symptom_onset,figpath = figpath,
                         filename = "days_since_sx_onset",
                         width = 7,height = 3.8)



# Number of sx reported by days since sx onset ----------------------------
dfRes$symptom_count_26 = rowSums(dfRes[,sympnames_type_df$symptom_code[1:26]], na.rm=T)
hist(dfRes$symptom_count_26)

p_symptom_count_onset <- dfRes %>% 
  filter(round!=1, estbinres==1, !is.na(variant_inferred_detail), 
         !is.na(days_since_symptom_onset)) %>% 
  mutate(variant_inferred=case_when(is.na(variant_inferred) ~ "Transition phase",
                                    round == 19 ~ "Round 19 (Omicron BA.2)",
                                    round == 17 ~ "Round 17 (Omicron BA.1)",
                                    round == 18 ~ "Transition phase",
                                    T ~ as.character(variant_inferred))) %>% 
  group_by(days_since_symptom_onset,variant_inferred_detail,
           days_since_symptom_onset_cat=factor(case_when(days_since_symptom_onset==11 ~ "11+",
                                                   T ~ as.character(days_since_symptom_onset)),
                                               levels = c(1:10,"11+"))) %>% 
  summarise(n=n(),
            mean_sx=mean(symptom_count_26, na.rm=T),
            sd_sx=sd(symptom_count_26, na.rm=T),
            se_sx=sd_sx/sqrt(n),
            mean_sx_lower=mean_sx-(1.96*se_sx),
            mean_sx_upper=mean_sx+(1.96*se_sx)
            ) %>% 
  ggplot(aes(x=days_since_symptom_onset_cat, y=mean_sx, col=variant_inferred_detail,
             fill=variant_inferred_detail)) +
  OverReact::scale_color_imperial(palette = "default") +
  OverReact::scale_fill_imperial(palette = "default") +
  geom_point() +
  geom_line() +
  geom_errorbar(aes(ymin=mean_sx_lower, ymax=mean_sx_upper))+
  # geom_point(position = position_jitter(width=0.1, height=0.1), alpha=0.1)+
  # ggbeeswarm::geom_quasirandom(alpha=0.1, groupOnX = F,position = position_jitter(width=10)) +
  OverReact::theme_react(strip_text_size = 9) +
  theme(legend.position = "none") +
  facet_wrap(.~variant_inferred_detail, nrow=1)+
  scale_y_continuous(breaks=seq(0,12,1)) +
  # scale_x_continuous(breaks=seq(0,11,1)) +
  labs(y="Number of symptoms reported",x="Days since symptom onset")
p_symptom_count_onset


p_symp_comb=p_symptom_count_onset/p_symptom_onset


# save
OverReact::saveREACTplot(p = p_symp_comb,figpath = figpath,
                         filename = "sx_panel_plot_days_sx_onset",
                         width = 9,height = 6)



# Symptom distribution ----------------------------------------------------
dfRes$random_number=as.numeric(sample(11:14,nrow(dfRes),replace = T))
dfRes$random_number=as.numeric(rep(11:14,(nrow(dfRes)/4)))

### Filter down data to get df for means
dat_for_plot=dfRes %>% 
  filter(round!=1, estbinres==1, !is.na(variant_inferred_detail)) %>% 
  mutate(variant_inferred=case_when(is.na(variant_inferred) ~ "Transition phase",
                                    round == 19 ~ "Round 19 (Omicron BA.2)",
                                    round == 17 ~ "Round 17 (Omicron BA.1)",
                                    round == 18 ~ "Transition phase",
                                    T ~ as.character(variant_inferred)),
         days_since_symptom_onset_cat=(case_when(days_since_symptom_onset==11 ~ "11+",
                                                 T ~ as.character(days_since_symptom_onset))),
         
         days_since_symptom_onset_inferred = (case_when(days_since_symptom_onset<11 ~ days_since_symptom_onset,
                                                        days_since_symptom_onset==11 ~ random_number,
                                                        T ~NA_real_)))


### Plot of distributions
p_symptom_count_onset_dist <- dat_for_plot %>% 
  group_by(days_since_symptom_onset_inferred,variant_inferred_detail) %>% 
  summarise(n=n(),
            mean_sx=mean(symptom_count_26, na.rm=T),
            sd_sx=sd(symptom_count_26, na.rm=T),
            se_sx=sd_sx/sqrt(n),
            mean_sx_lower=mean_sx-(1.96*se_sx),
            mean_sx_upper=mean_sx+(1.96*se_sx)
  ) %>% 
  ggplot(aes(x=days_since_symptom_onset_inferred, y=n, col=variant_inferred_detail,
             fill=variant_inferred_detail)) +
  OverReact::scale_color_imperial(palette = "default") +
  OverReact::scale_fill_imperial(palette = "default") +
  geom_col()+
  # geom_point(position = position_jitter(width=0.1, height=0.1), alpha=0.1)+
  # ggbeeswarm::geom_quasirandom(alpha=0.1, groupOnX = F,position = position_jitter(width=10)) +
  # OverReact::theme_react(strip_text_size = 9) +
  OverReact::theme_react(strip_text_size = 9) +
  theme(legend.position = "none") +
  facet_wrap(.~variant_inferred_detail, nrow=5)+
  # scale_y_continuous(breaks=seq(0,12,1)) +
  scale_x_continuous(breaks=seq(0,14,1)) +
  labs(title = "A", x="Days since symptom onset*", y="N")
p_symptom_count_onset_dist



class(means)
means=dat_for_plot %>% group_by(variant_inferred_detail) %>% 
  summarise(n=n(),
    mean_dur=round(mean(days_since_symptom_onset_inferred,na.rm=T),2),
            sd_dur=sd(days_since_symptom_onset_inferred,na.rm=T),
    se_dur=sd_dur/sqrt(n),
    lower=round(mean_dur-1.96*se_dur,2),
    upper=round(mean_dur+1.96*se_dur,2)) %>% 
  as.data.frame() %>% 
  mutate(mean_dur_comb=paste0(mean_dur," [",lower,"-",upper,"]"))

p_symptom_count_onset_boxplot <- dat_for_plot %>% 
  # left_join(means) %>% 
  ggplot(aes(y=days_since_symptom_onset_inferred, x=variant_inferred_detail,
             col=variant_inferred_detail)) +
  OverReact::scale_color_imperial(palette = "default") +
  geom_boxplot()+
  ggpubr::stat_compare_means(comparisons = list(c(1,2),c(2,3),c(3,4),c(4,5)
                                                # c(1,4),c(2,4),c(3,4),c(1,5),c(2,5),c(3,5)
                                                )
                             ,size=3,label = "p.signif") +
  geom_text(data=means,aes(label=paste0("Mean = ",mean_dur_comb),
            x=variant_inferred_detail,
                           y=mean_dur),inherit.aes = F, size=3,nudge_y = 0.4) +
  OverReact::theme_react(strip_text_size = 9) +
  theme(legend.position = "none") +
  labs(title = "A", y="Days since symptom onset*", x="")
p_symptom_count_onset_boxplot


## add title label to p_symptom_count_onset
p_symptom_count_onset <- p_symptom_count_onset+labs(title="B")

### Create panel layouts
p_symp_comb_2=(p_symptom_count_onset_dist)/p_symptom_count_onset +
  plot_layout(heights=c(4,2))
p_symp_comb_2

p_symp_comb_3=p_symptom_count_onset_boxplot/p_symptom_count_onset +
  plot_layout(heights=c(5,3))
p_symp_comb_3




# save
OverReact::saveREACTplot(p = p_symptom_count_onset_dist,figpath = figpath,
                         filename = "days_sx_onset_dist",
                         width = 8,height = 8)

# save
OverReact::saveREACTplot(p = p_symptom_count_onset_boxplot,figpath = figpath,
                         filename = "days_sx_onset_box",
                         width = 10,height = 5)
# save
OverReact::saveREACTplot(p = p_symp_comb_2,figpath = figpath,
                         filename = "sx_2_panel_plot_days_sx_onset_dist",
                         width = 10,height = 10)
# save
OverReact::saveREACTplot(p = p_symp_comb_3,figpath = figpath,
                         filename = "sx_2_panel_plot_days_sx_onset_boxplot",
                         width = 10,height = 10)



# Distribution of symptom counts in 2vax + 3 vax --------------------------

# create dodger
dodger=position_dodge2(width = 0.8)
# add age group
dfRes <- dfRes %>% mutate(age_group_broad=case_when(age<=40 ~ "Age 18-40",
                                                    age<=60 ~ "Age 41-60",
                                                    T ~ "Age 61+"),
                          age_group_broader=case_when(age<=40 ~ "Age 18-40",
                                                    T ~ "Age 41+"))
dfRes$vaccine_status_named %>% table()


p_symptom_count_onset_omicron_vax <- dfRes %>% 
  filter(round!=1, estbinres==1, !is.na(variant_inferred_detail), 
         !is.na(days_since_symptom_onset), round %in% c(17:19),
         vaccine_status_named %in% c("Three doses","Two doses")) %>% 
  mutate(variant_inferred=case_when(is.na(variant_inferred) ~ "Transition phase",
                                    round == 19 ~ "Round 19 (Omicron BA.2)",
                                    round == 17 ~ "Round 17 (Omicron BA.1)",
                                    round == 18 ~ "Transition phase",
                                    T ~ as.character(variant_inferred)),
         
         days_since_symptom_onset_cat=factor(case_when(days_since_symptom_onset==11 ~ "11+",
                                                       T ~ as.character(days_since_symptom_onset)),
                                             levels = c(1:10,"11+"))) %>% 
  group_by(days_since_symptom_onset_cat,variant_inferred_detail,
           age_group_broader,vaccine_status_named) %>% 
  summarise(n=n(),
            mean_sx=mean(symptom_count_26, na.rm=T),
            sd_sx=sd(symptom_count_26, na.rm=T),
            se_sx=sd_sx/sqrt(n),
            mean_sx_lower=max(0,mean_sx-(1.96*se_sx)),
            mean_sx_upper=mean_sx+(1.96*se_sx)
  ) %>% 
  mutate(vaccine_status_named=factor(vaccine_status_named, levels=c("Two doses",
                                                                    "Three doses"))) %>% 
  ggplot(aes(x=days_since_symptom_onset_cat, y=mean_sx, col=vaccine_status_named,
             fill=vaccine_status_named, group=vaccine_status_named)) +
  # OverReact::scale_color_imperial(palette = "default",reverse =T) +
  # OverReact::scale_fill_imperial(palette = "default",reverse =T) +
  scale_fill_manual(values = c(OverReact::imperial_palettes$cool[[1]],
                               OverReact::imperial_palettes$cool[[7]])) +
  scale_colour_manual(values = c(OverReact::imperial_palettes$cool[[1]],
                               OverReact::imperial_palettes$cool[[7]])) +
  geom_point() +
  geom_line() +
  # geom_errorbar(aes(ymin=mean_sx_lower, ymax=mean_sx_upper),
  #               position=dodger, alpha=0.3)+
  geom_ribbon(aes(ymin=mean_sx_lower, ymax=mean_sx_upper),
              position=dodger, alpha=0.1, col=NA)+
  OverReact::theme_react(strip_text_size = 9) +
  theme(legend.position = "bottom") +
  facet_wrap(.~variant_inferred_detail*age_group_broader, nrow=2)+
  labs(y="Number of symptoms reported",x="Days since symptom onset",
       fill="Vaccination status",col="Vaccination status")

p_symptom_count_onset_omicron_vax

# save
OverReact::saveREACTplot(p = p_symptom_count_onset_omicron_vax,figpath = figpath,
                         filename = "sx_by_onset_time_vax_age",
                         width = 7,height = 5.5)




# Test negatives symptom prevalence ---------------------------------------

cov_name_list$round="Round"
dfRes$dummy="All"
dfRes$round <- as.character(dfRes$round)
xtab_symptomatic_adults_negs=OverReact::makeTablesNew(dat = dfRes %>% filter(age>=18, estbinres==0),
                                                      result_var = "symptomatic_26",covariates = "round",
                                                      cov_name_list =cov_name_list,weights = NULL,
                                                      output_list = F,sens = 1,spec = 1)
xtab_symptomatic_adults_pos=OverReact::makeTablesNew(dat = dfRes %>% filter(age>=18, estbinres==1),
                                                      result_var = "symptomatic_26",covariates = "round",
                                                      cov_name_list =cov_name_list,weights = NULL,
                                                      output_list = F,sens = 1,spec = 1)


### Create round date pmidpoints
round_date_midpoints_df <- dfRes %>% 
  group_by(round) %>% 
  summarise(n=n(),
            round_date=median(d_comb, na.rm=T))
round_date_midpoints_df$round <- as.character((round_date_midpoints_df$round))

# Add ages
xtab_symptomatic_adults_negs <- xtab_symptomatic_adults_negs %>% 
  left_join(round_date_midpoints_df, by = c("Category"="round"))
xtab_symptomatic_adults_pos <- xtab_symptomatic_adults_pos %>% 
  left_join(round_date_midpoints_df, by = c("Category"="round"))
xtab_symptomatic_adults_pos$PCR="PCR positive"
xtab_symptomatic_adults_negs$PCR="PCR negative"


# combine and rename
xtab_symptomatic_comb=rbind(xtab_symptomatic_adults_negs,xtab_symptomatic_adults_pos) %>% 
  rename(Positive = Total,
         Symptomatic = Positive,
         Prop_symptomatic = Prevalence,
         Round = Category)

# Plot
pd=position_dodge(width=20)
maxdate=max(dfRes$d_comb,na.rm=T)

p_prop_negatives=xtab_symptomatic_comb %>% 
  filter(Round!=1) %>% 
  mutate(Round=factor(Round, levels=2:19),
         date=as.Date(round_date,format = "%d-%m-%Y")) %>% 
  ggplot(aes(x=date, y= Prop_symptomatic)) + 
  geom_point(position = pd) +
  geom_errorbar(aes(x=date, y= Prop_symptomatic,ymin =Lower, ymax = Upper), 
                width=10,
                position = pd) +
  scale_x_date(date_breaks = "2 month",date_labels = "%B \n%Y") +
  theme_bw() + 
  theme_adjust +
  scale_color_manual(values = myCols[c(5,1)]) +
  # annotate("rect",xmin=0.5, xmax=6.5, ymin=0, ymax=Inf, alpha=0.1) +
  annotate("rect",xmin=as.Date("23-12-2020",format = "%d-%m-%Y"), 
           xmax=as.Date("9-5-2021",format = "%d-%m-%Y"),
           ymin=0, ymax=Inf, alpha=0.09) +
  annotate("rect",xmin=as.Date("9-5-2021",format = "%d-%m-%Y"), 
           xmax=as.Date(as.Date("01-01-2022",format = "%d-%m-%Y"),format = "%d-%m-%Y"), 
           ymin=0, ymax=Inf, alpha=0.19)+
  annotate("rect",xmin=as.Date("01-01-2022",format = "%d-%m-%Y"), 
           xmax=as.Date(maxdate,format = "%d-%m-%Y"), 
           ymin=0, ymax=Inf, alpha=0.29)+
  annotate("text",x=as.Date("23-08-2020",format = "%d-%m-%Y"), y=80, label = "Wild type") +
  annotate("text",x=as.Date("28-02-2021",format = "%d-%m-%Y"), y=80, label = "Alpha")  +
  annotate("text",x=as.Date("13-09-2021",format = "%d-%m-%Y"), y=80, label = "Delta")  +
  annotate("text",x=as.Date("15-02-2022",format = "%d-%m-%Y"), y=80, label = "Omicron")  +
  facet_wrap(.~PCR,ncol=1)+
  labs(x="Date", y="Prevalence of any of 26 symptoms \n in REACT-1", col = "")
p_prop_negatives

OverReact::saveREACTplot(p = p_prop_negatives,figpath = figpath,filename = "background_symptom_prevalence",
                         width = 10, height =6
)

