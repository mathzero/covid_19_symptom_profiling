
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
                  "patchwork","randomcoloR","car","ggnewscale"
)

load_packages(package.list)


# Import REACT-1 data -----------------------------------------------------
source("E:/Group/react2_study5/report_phases_combined/projects/omicron_symptom_profiling/code/00_bits_and_pieces.R", local = T)
source("E:/Group/react2_study5/report_phases_combined/projects/omicron_symptom_profiling//code/00_data_prep.R", local = T)

# create subfolder
createMySubfolder(subfolderName = "univariate_time_since_sx_onset")
suffix="withoutpriorinfection"

# Omicron = 17
# Delta = 12-15
# Alpha = 8-10
# Wildtype = 2-7
# clean up key symptom variables
dfRes <- dfRes %>% 
  mutate_at(covid_yesnos_firstsymp,binaryCleaner_1_0) %>% 
  mutate_at(covid_yesnos_month,binaryCleaner_1_0)  %>% 
  mutate(doublevaxxed_or_boosted=case_when(vax_status_number>=2 ~1,
                                           TRUE ~0),
         sex=factor(sex),
         days_since_symptom_onset_cat = case_when(days_since_symptom_onset >10 ~"11 days or more",
                                                  days_since_symptom_onset >5~"6-10 days",
                                                  days_since_symptom_onset >0~"5 days or less",
                                                  T ~ NA_character_)) %>% 
  mutate_at(c(sympnames_type_df$symptom_code,covid_yesnos_firstsymp,covid_yesnos_month,
              "sex","doublevaxxed_or_boosted"),as.factor) 

dfRes$prior_covid <- as.factor(dfRes$prior_covid)
dfRes$age_group_named <- as.factor(dfRes$age_group_named)
dfRes$days_since_symptom_onset_cat <- as.factor(dfRes$days_since_symptom_onset_cat)

# Run models --------------------------------------------------------------



runUnivariate <- function(dat, variant="Omicron",rounds =c(17), adj_level =5,
                          joint_adjustment_vars = c("age_group_named", "sex", 
                                                    "days_since_symptom_onset_cat",
                                                    # "prior_covid",
                                                    "vax_status_number")){
  
  # OMICRON #
  uiv_omicron =OverReact::ModelMakerMulti(dat = dat %>% filter(round %in% rounds),
                                          list_of_variables_of_interest = sympnames_type_df$symptom_code,
                                          outcome = "estbinres",
                                          cov_name_list = cov_name_list,
                                          joint_adjustment_vars = joint_adjustment_vars)
  # extract plottable dat
  plot_df <- uiv_omicron$plot_output
  plot_df$Category <- plot_df$Variable
  plot_df$Variable <- "Symptom in past week"
  
  
  
  # Get summary df with average vals by group
  mean_ors <- plot_df %>% group_by(Category) %>% 
    filter(adjustment==adj_level, !is.na(OR)) %>% 
    summarise(mean_or=mean(OR, na.rm=T)) %>% 
    arrange(-mean_or) %>% 
    ungroup()
  
  bold_vect=rev(case_when(mean_ors$Category %in% sympnames_type_df$symptom[c(27:28)] ~ "bold",
                      T ~ "plain"))
  
  # join with original data
  plot_df_withmean <- plot_df %>% 
    filter(adjustment==adj_level, !is.na(OR)) %>% 
    left_join(mean_ors) %>% 
    arrange(-mean_or) %>% 
    filter(!is.na(OR)) %>% 
    mutate(Category =factor(Category, levels=unique(Category)))
  plot_df_withmean$variant <- variant
  p_compare <- plotReactForest(univ_df_plot = plot_df_withmean,adjustment_numbers = adj_level,
                               adjustment_descriptions = c("Age, sex, vaccination, prior infection status adjusted"),
                               insignificant_results_greyed_out = T,
                               legend.position = "bottom", 
                               palette = "cool")
  p_compare <- p_compare + labs(title = paste0(variant, " variant"), 
                   subtitle = "Odds ratios from logistic regression, adjusted for age, sex, vaccination and prior infection status") +
    theme(legend.position = "none",
          axis.text.y = element_text(face = bold_vect))
  
  return(list(plot=p_compare,
              results=plot_df_withmean))

}

# 
# dat$symptomatic_26 %>% table(dat$days_since_symptom_onset_cat,exclude="none")
# dat$symptnowaw_1 %>% table(dat$days_since_symptom_onset_cat,exclude="none")
# dat$symptnowaw_1 %>% table(dat$symptomatic_26,exclude="none")
# dat$symptnowaw_1 %>% table(dat$estbinres,exclude="none")

### BA.1
univ_ba1=runUnivariate(dat=dfRes%>% filter(round %in% c(17:19) & 
                                      (variant_inferred_detail ==  "BA.1 (Omicron)" | 
                                         estbinres==0) # filter for BA1/negative
                             ), 
                       variant="BA1",rounds =c(17,18,19))
univ_ba1$plot


### BA.2
univ_ba2=runUnivariate(dfRes%>%  filter(round %in% c(17:19) & 
                                      (variant_inferred_detail ==  "BA.2 (Omicron)" | estbinres==0) # filter for BA1/negative
                             ), variant="BA2",rounds =c(17,18,19))
univ_ba2$plot

### Delta
univ_delta=runUnivariate(dfRes, variant="Delta",rounds =c(12:15))
### Alpha
univ_alpha=runUnivariate(dfRes, variant="Alpha",rounds =c(8:10))
### WIld type
univ_wildtype=runUnivariate(dfRes, variant="Wildtype",rounds =c(2:7),adj_level = 4)

### output filtered for plots
u_w <- univ_wildtype$results %>% filter(adjustment==4) %>% 
  mutate(OR_concat = paste0(OverReact::specifyDecimal(OR,3,format = "g"), " (", OverReact::specifyDecimal(Lower,3,format = "g"),",",OverReact::specifyDecimal(Upper,3,format = "g"),")")) %>% 
  select(-mean_or)%>% 
  arrange(Category)
u_a <- univ_alpha$results %>% filter(adjustment==5) %>% 
  mutate(OR_concat = paste0(OverReact::specifyDecimal(OR,3,format = "g"), " (", OverReact::specifyDecimal(Lower,3,format = "g"),",",OverReact::specifyDecimal(Upper,3,format = "g"),")")) %>% 
  select(-mean_or)%>% 
  arrange(Category)
u_d <- univ_delta$results %>% filter(adjustment==5) %>% 
  mutate(OR_concat = paste0(OverReact::specifyDecimal(OR,3,format = "g"), " (", OverReact::specifyDecimal(Lower,3,format = "g"),",",OverReact::specifyDecimal(Upper,3,format = "g"),")")) %>% 
  select(-mean_or)%>% 
  arrange(Category)
u_ba1 <- univ_ba1$results %>% filter(adjustment==5) %>% 
  mutate(OR_concat = paste0(OverReact::specifyDecimal(OR,3,format = "g"), " (", OverReact::specifyDecimal(Lower,3,format = "g"),",",OverReact::specifyDecimal(Upper,3,format = "g"),")")) %>% 
  select(-mean_or) %>% 
  arrange(Category)
u_ba2 <- univ_ba2$results %>% filter(adjustment==5) %>% 
  mutate(OR_concat = paste0(OverReact::specifyDecimal(OR,3,format = "g"), " (", OverReact::specifyDecimal(Lower,3,format = "g"),",",OverReact::specifyDecimal(Upper,3,format = "g"),")")) %>% 
  select(-mean_or) %>% 
  arrange(Category)

# combine all univariate data
u_comb <- bind_rows(u_w,u_a,u_d,u_ba1,u_ba2) %>% 
  select(Category,variant,OR_concat) %>% 
  pivot_wider(id_cols = Category,names_from = variant,values_from = OR_concat)

### Reorder to save faff in the doc
u_comb <- u_comb[c(15,1:14,16:nrow(u_comb)),]




# save nice looking workbook
savePrettyExcelWorkbook(listOfTables = list(wildtype=u_w,
                                            alpha=u_a,
                                            delta=u_d,
                                            ba1=u_ba1,ba2=u_ba2,combined=u_comb),
                        workbookName = paste0("univariate_results_symp_onset_adj",suffix),
                        outpath = outpath)


# Bind all results and combine --------------------------------------------

df_all <- rbind(univ_ba1$results,univ_ba2$results,univ_delta$results,
                univ_alpha$results,univ_wildtype$results) 
# %>% 
#   filter(adjustment%in%c(3,4))
df_all <- df_all[!(df_all$variant != "Wildtype" & df_all$adjustment ==3),]

# Pivot wider
df_all_wide=df_all %>% pivot_wider(id_cols = Category,names_from = variant, values_from = c(OR,Lower, Upper))

df_all_wide <- df_all_wide %>% left_join(sympnames_type_df, by = c("Category" = "symptom"))


# Regular all-in-one forest plot ------------------------------------------

df_all_plot <- df_all %>% 
  mutate(adjustment = (variant))

bold_vect=rev(c(rep("plain",12),"bold", rep("plain",14)))


# do plot
p_forest <- plotReactForest(univ_df_plot = df_all_plot %>%  filter(Category!="Any of 26 symptoms"),
                            adjustment_numbers = c("Wildtype","Alpha","Delta","BA1","BA2"),
                            adjustment_descriptions = c("Wildtype","Alpha","Delta","BA1","BA2"),
                            insignificant_results_greyed_out = T,
                            palette = "default")+
  theme(axis.text.y = element_text(face = bold_vect))


p_forest

OverReact::saveREACTplot(p = p_forest,figpath = figpath,
                         filename = paste0("univ_forest_symp_onset_adjust",suffix),
                         width = 5,height = 8.5,savePDF = T)





# Forest with columns split by symptom type -------------------------------

df_all_plot_2 <- df_all_plot %>% left_join(sympnames_type_df, by = c("Category"= "symptom"))
bold_vect_2=rev(c(rep("plain",9),"bold", rep("plain",17)))


# do plot
p_forest_2 <- plotReactForest(univ_df_plot = df_all_plot_2%>%  
                                filter(Category!="Any of 26 symptoms") %>% 
                                filter(adjustment%in%c("BA1","BA2")),
                            adjustment_numbers = c("Wildtype","Alpha","Delta","BA1","BA2"),
                            adjustment_descriptions = c("Wildtype","Alpha","Delta","BA1","BA2"),
                            insignificant_results_greyed_out = T,
                            palette = "default")+
  theme(axis.text.y = element_text(face = bold_vect_2))+
  scale_colour_manual(values = c( "#EC7300", "#960078" ))


p_forest_2

OverReact::saveREACTplot(p = p_forest_2,figpath = figpath,
                         filename = paste0("univ_forest_ba1_ba2_compare",suffix),
                         width = 5,height = 8.5,savePDF = T)


# Forest without Omicron --------------------------------------------------


# Forest with columns split by symptom type -------------------------------

df_all_plot_3 <- df_all_plot %>% left_join(sympnames_type_df, by = c("Category"= "symptom"))
bold_vect_3=rev(c(rep("plain",12),"plain", rep("plain",14)))


# do plot
p_forest_3 <- plotReactForest(univ_df_plot = df_all_plot_3 %>%  
                                filter(Category!="Any of 26 symptoms")%>%
                                filter(!adjustment%in%c("BA1","BA2")),
                              adjustment_numbers = c("Wildtype","Alpha","Delta","BA1","BA2"),
                              adjustment_descriptions = c("Wildtype","Alpha","Delta","BA1","BA2"),
                              insignificant_results_greyed_out = T,
                              palette = "default")+
  theme(axis.text.y = element_text(face = bold_vect_3)) +
  scale_colour_manual(values = c("#002147","#00ACD7","#BBCE00"))

p_forest_3


### Combine and save
p_comb=p_forest_3+p_forest_2

OverReact::saveREACTplot(p = p_comb,figpath = figpath,
                         filename = paste0("univ_forest_panel",suffix),
                         width = 8,height = 8.5,savePDF = T)



# Apply eda plot layout ---------------------------------------------------

# 
# 
# 
# stripes_df <- df_all_plot_3 %>% 
#   arrange(-percentage) %>% 
#   filter(name=="Rounds 2-7 (Wild type)") %>% 
#   left_join(sympnames_type_df, by = c("Variable"= "symptom"))%>% 
#   group_by(symptom_type) %>% 
#   mutate(stripe_col = case_when(row_number() %%2 ==0 ~ "NA", TRUE ~ "grey90")
#   )


### Get summary df to fix column widths
summ <- df_all_plot_3 %>% 
  group_by(symptom_type) %>% 
  summarise(n=n())

# dodge
# dodger <- ggstance::position_dodge2v(height =  dodge_mult,preserve = "single",reverse = F)
dodge_mult=0.9
dodger <- position_dodge2(width = dodge_mult,reverse = F)


# join
sympnames_type_df <- sympnames_type_df %>% left_join(summ)

# set col width
sympnames_type_df$col_width=dodge_mult*sympnames_type_df$n/max(sympnames_type_df$n)


# Wrap --------------------------------------------------------------------


# dodge
# dodger <- ggstance::position_dodge2v(height =  dodge_mult,preserve = "single",reverse = F)
dodge_mult=0.5
dodger <- position_dodge2(width = dodge_mult,reverse = F)


df_all_plot_3$Category
df_all_plot_3 <- df_all_plot_3 %>% filter(Category!="Any of 26 symptoms")
plots=list()
x=1
for(symp in c("Smell/taste","Respiratory/cardiac", 
                "Coryzal (cold-like)", "Gastrointestinal","Fatigue", "Other","Influenza-like")){
  # PLot
  p_prevs=df_all_plot_3 %>% 
    left_join(sympnames_type_df) %>% 
    mutate(
      variant = factor(variant, levels = c("Wildtype","Alpha","Delta","BA1","BA2")),
      symptom_type = factor(symptom_type, levels = c("Overall",  "Smell/taste","Respiratory/cardiac", 
                                                     "Coryzal (cold-like)", "Gastrointestinal","Fatigue", "Other","Influenza-like"))) %>% 
    filter(symptom_type==symp) %>% 
    ggplot(aes(group=variant)) +
    geom_hline(yintercept = 1, size=0.3, linetype="dashed", col ="grey30") +
    geom_hline(yintercept = 4, size=0.2, linetype="dashed", col ="grey50") +
    geom_hline(yintercept = 16, size=0.2, linetype="dashed", col ="grey50") +
    geom_hline(yintercept = 64, size=0.2, linetype="dashed", col ="grey50") +
    geom_point(aes(x=reorder(Category,OR), y= OR, col = variant,width=col_width),
               size=1.3,position = dodger, shape = 15) +
    geom_errorbar(aes(x=reorder(Category,OR),y=OR,ymin= Lower, ymax=Upper,col = variant),
                  size=0.6,position = dodger, width=0.5) +
    scale_x_discrete(labels = function(x) str_wrap(x, width =17)) +
    scale_y_continuous(trans = "log2", #breaks = x_seq, 
                       labels =function(x) MASS::fractions(x),
                       limits = c(0.5,
                                  max(df_all_plot_3$Upper)))+     
    OverReact::scale_fill_imperial(palette = "default") +
    OverReact::scale_color_imperial(palette = "default") +
    theme_react(strip_text_size = 10) +
    # coord_flip() +
    # facet_grid(scales = "free",rows = "symptom_type",space = "fixed",shrink = T) +
    ggforce::facet_col(facets = "symptom_type",scales = "free",space = "free") +
    labs(x="", y="Odds ratio", fill = "", col ="", fill = "") +
    theme(legend.position = "bottom",
          panel.spacing = unit(1, "cm"),
          panel.grid.minor.x = element_blank(),
          panel.grid.major.x = element_blank())
  
  if(!x%in%c(1,2,4)){
    p_prevs <- p_prevs+theme(axis.text.y = element_blank(),
                             axis.title.y = element_blank())
  }
  
  plots[[x]]=p_prevs
  x=x+1
}


# Create design
layout="
CCCCCCCCCCCCCCCCHHHHHHHHHHHHHHHH
BBBBBBBBBBBBBBDDDDDDDDDDDDDDDDDD
EEEEEEEEEEEEEEGGGGGGGGGGGGGGGGGG
"



p_comb <- patchwork::wrap_plots(plots,guides = "collect", 
                      design = layout
                      ) & 
  theme(legend.position = "bottom")
p_comb

# save
OverReact::saveREACTplot(p = p_comb,figpath = figpath,
                         filename = paste0("symptom_ors_by_variant_adjusted_time_since_sx_onset_",suffix),
                         width = 11,height = 7.5)







       