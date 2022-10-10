
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
                  "patchwork","randomcoloR"
)

load_packages(package.list)


# Import REACT-1 data -----------------------------------------------------
source("E:/Group/react2_study5/report_phases_combined/projects/omicron_symptom_profiling/code/00_bits_and_pieces.R",
       local = T)
source("E:/Group/react2_study5/report_phases_combined/projects/omicron_symptom_profiling/code/00_data_prep.R", 
       local = T)


# create subfolder
createMySubfolder(subfolderName = "pooled_analysis")




# Run analysis ------------------------------------------------------------

symp=sympnames_type_df$symptom_code[[1]]
dfRes$variant_inferred %>% class()
dfRes$vax_status_number %>% table(dfRes$round, exclude="none")
dfRes <- dfRes %>% mutate(vax_count_imputed=case_when(round %in% 1:8 ~ 0,
                                                      T ~ vax_status_number),
                          variant_inferred_new=factor(case_when(round==17 ~ "Round 17 (Omicron BA.1)",
                                                         round==19 ~ "Round 19 (Omicron BA.2)",
                                                         T ~ as.character(variant_inferred)),
                                                      levels= c("Rounds 2-7 (Wild type)",
                                                                "Rounds 8-10 (Alpha)",
                                                              "Rounds 12-15 (Delta)",
                                                              "Round 17 (Omicron BA.1)",
                                                              "Round 19 (Omicron BA.2)")))
dfRes$variant_inferred_new %>% table(dfRes$round)

runFunctions <- function(symp){
  
  ### Define formula
  interaction=paste(c(symp,"variant_inferred_new"),collapse = ":")
  f=paste0("estbinres ~ ",paste(c("age_group_named", "sex", 
                                  "vax_count_imputed",
                                  "variant_inferred_new",
                                  symp),collapse = "+"),"+",
           interaction
  )
  f
  
  # Run models
  mymod_wild_alpha=glm(formula = f,family = "binomial",data = dfRes %>% filter(round%in%c(2:10)))
  mymod_alpha_delta=glm(formula = f,family = "binomial",data = dfRes %>% filter(round%in%c(8:15)))
  print("Halfway!")
  mymod_delta_BA1=glm(formula = f,family = "binomial",data = dfRes %>% filter(round%in%c(12:15,17)))
  mymod_BA1_BA2=glm(formula = f,family = "binomial",data = dfRes %>% filter(round%in%c(17,19)))
  
  
  ### get model summaries
  summ_wild_alpha=jtools::summ(mymod_wild_alpha,exp=T)[[1]] %>%  as.data.frame()
  summ_alpha_delta=jtools::summ(mymod_alpha_delta,exp=T)[[1]] %>%  as.data.frame()
  summ_delta_ba1=jtools::summ(mymod_delta_BA1,exp=T)[[1]] %>%  as.data.frame()
  summ_ba1_ba2=jtools::summ(mymod_BA1_BA2,exp=T)[[1]] %>%  as.data.frame()
  
  ### add interaction column description
  summ_wild_alpha$variant_interact="Wild-type:Alpha"
  summ_alpha_delta$variant_interact="Alpha:Delta"
  summ_delta_ba1$variant_interact="Delta:BA.1"
  summ_ba1_ba2$variant_interact="BA.1:BA.2"
  summ_wild_alpha$variable=rownames(summ_wild_alpha)
  summ_alpha_delta$variable=rownames(summ_alpha_delta)
  summ_delta_ba1$variable=rownames(summ_delta_ba1)
  summ_ba1_ba2$variable=rownames(summ_ba1_ba2)
  
  
  res=rbind(summ_wild_alpha,summ_alpha_delta,summ_delta_ba1,summ_ba1_ba2)
  res$symptom = sympnames_type_df$symptom[sympnames_type_df$symptom_code==symp]
  res <- res %>% select(symptom,variable, everything())
  return(res)
}


### Run over all symptoms

myresults_list=list()
for(symp in sympnames_type_df$symptom_code){
  print(symp)
  myresults_list[[symp]] <- runFunctions(symp)
}
names(myresults_df)
myresults_df=bind_rows(myresults_list)
myresults_df <- myresults_df %>% filter(grepl("variant_inferred_new",variable))
myresults_df <- myresults_df %>% filter(grepl(":",variable))
myresults_df <- myresults_df %>% janitor::clean_names()
myresults_df <- myresults_df %>% left_join(sympnames_type_df)
# sympnames_type_df$symptom_type
# 
# dodge_mult=0.8
# dodger <- position_dodge2(width = dodge_mult,reverse = F)
# 
# myresults_df %>% 
#   ggplot(aes(x=symptom,y = exp_est, 
#              col=variant_interact, fill=variant_interact)) +
#   geom_point(position = dodger) + 
#   geom_errorbar(aes(ymin=x2_5_percent,ymax=x97_5_percent),position = dodger) +
#   OverReact::theme_react() +
#   ggforce::facet_col(.~symptom_type,scales = "free_x")
# 
# 

# Plot layout grid --------------------------------------------------------


### Get summary df to fix column widths
summ <- myresults_df %>% 
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

# 
# # dodge
# # dodger <- ggstance::position_dodge2v(height =  dodge_mult,preserve = "single",reverse = F)
dodge_mult=0.5
dodger <- position_dodge2(width = dodge_mult,reverse = F)

symptype="Smell/taste"
plots=list()
x=1
for(symptype in c("Overall","Smell/taste","Respiratory/cardiac", 
              "Coryzal (cold-like)", "Gastrointestinal","Fatigue", "Other","Influenza-like")){
  # PLot
  p_prevs=myresults_df %>% 
    mutate(
      variant_interact = factor(variant_interact, levels = c("Wild-type:Alpha","Alpha:Delta","Delta:BA.1","BA.1:BA.2")),
      symptom_type = factor(symptom_type, levels = c("Overall",  "Smell/taste","Respiratory/cardiac", 
                                                     "Coryzal (cold-like)", "Gastrointestinal","Fatigue", "Other","Influenza-like"))) %>% 
    filter(symptom_type==symptype) %>% 
    ggplot(aes(group=variant_interact)) +
    geom_hline(yintercept = 1, size=0.3, linetype="dashed", col ="grey30") +
    geom_hline(yintercept = 4, size=0.2, linetype="dashed", col ="grey50") +
    geom_hline(yintercept = 16, size=0.2, linetype="dashed", col ="grey50") +
    # geom_hline(yintercept = 64, size=0.2, linetype="dashed", col ="grey50") +
    geom_point(aes(x=reorder(symptom,exp_est), y= exp_est, col = variant_interact,width=col_width),
               size=1.3,position = dodger, shape = 15) +
    geom_errorbar(aes(x=reorder(symptom,exp_est),y=exp_est,ymin= x2_5_percent, ymax=x97_5_percent,col = variant_interact),
                  size=0.6,position = dodger, width=0.5) +
    scale_x_discrete(labels = function(x) str_wrap(x, width =17)) +
    scale_y_continuous(trans = "log2", #breaks = x_seq, 
                       labels =function(x) MASS::fractions(x),
                       limits = c(min(myresults_df$x2_5_percent),
                                  max(myresults_df$x97_5_percent)))+     
    # OverReact::scale_fill_imperial(palette = "default") +
    # OverReact::scale_color_imperial(palette = "default") +
    scale_color_manual(values = as.character(OverReact::imperial_palettes$default[c(2:5)])) +
    scale_fill_manual(values = as.character(OverReact::imperial_palettes$default[c(2:5)])) +
    theme_react(strip_text_size = 10) +
    # coord_flip() +
    # facet_grid(scales = "free",rows = "symptom_type",space = "fixed",shrink = T) +
    ggforce::facet_col(facets = "symptom_type",scales = "free",space = "free") +
    labs(x="", y="Odds ratio", fill = "", col ="", fill = "") +
    theme(legend.position = "bottom",
          panel.spacing = unit(1, "cm"),
          panel.grid.minor.x = element_blank(),
          panel.grid.major.x = element_blank())
  p_prevs
  if(!x%in%c(1,2,4)){
    p_prevs <- p_prevs+theme(axis.text.y = element_blank(),
                             axis.title.y = element_blank())
  }
  
  plots[[x]]=p_prevs
  x=x+1
}

# Create design
layout="
AAAACCCCCCCCCCCCCCHHHHHHHHHHHHHH
DDDDDDDDDDDDDDDDDDDDDFFFFFFFFFFF
BBBBBBBEEEEEEEEEEEEEEGGGGGGGGGGG
"



p_comb <- patchwork::wrap_plots(plots,guides = "collect", 
                                design = layout
) & 
  theme(legend.position = "bottom")
p_comb


# save
OverReact::saveREACTplot(p = p_comb,figpath = figpath,
                         filename = paste0("pooled_analysis_interaction_effects"),
                         width = 11,height = 6.5)


  
