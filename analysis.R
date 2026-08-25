# Load packages
library(tidyverse)
setwd("/data/Wolfson-PNU-dementia/UKB_Projects/Piranavi/scripts/")

# read in 
dat = readRDS("/data/Wolfson-PNU-dementia/UKB_Projects/Piranavi/UKB_Piranavi_78867r675516_010726.rds") %>% tibble()
disease_codes = readRDS("../../../UKB_BJ/LegacyC100091/UKB_Category100091_78867r677320_Legacy_Raw.rds")
disease_codes = tibble(disease_codes)
proteomics = readRDS("/data/Wolfson-PNU-dementia/UKB_BJ/proteomics_MRI/outputs/proteomics_wide_qc.rds")

# curate code list
codes = read_lines("/data/Wolfson-PNU-dementia/UKB_Projects/Piranavi/icd_infection_codes.txt") %>% 
  str_split(",") %>% 
  unlist() %>% 
  str_remove_all("\\.") %>% 
  data.frame() %>% 
  dplyr::rename("code" = ".") %>% 
  tibble() %>% 
  mutate(code = str_remove_all(code," ")) %>%
  filter(code != "") 

# expand codes with a - 
extra_codes = list()
codes_with_hypen = codes %>% filter(grepl("-",code))  
for(i in c(1:nrow(codes_with_hypen))){
  message(i)
  split_code = codes_with_hypen[i,]$code %>% str_split_fixed("-",n=2)
  letter = substr(split_code[1],0,1)
  start = substr(split_code[1],2,nchar(split_code[1]))
  end = substr(split_code[2],2,nchar(split_code[2]))
  if(substr(start,1,2)==00){
    extra_codes[[length(extra_codes)+1]] = paste0(letter,"00",seq(start,end,by=1) )
  } else if(substr(start,1,1)==0){
    extra_codes[[length(extra_codes)+1]] = paste0(letter,"0",seq(start,end,by=1) )
  } else  {
        extra_codes[[length(extra_codes)+1]] = paste0(letter,seq(start,end,by=1) )
  }
}

codes_without_hypen = codes %>% filter(!grepl("-",code))  
all_codes = c(codes_without_hypen$code,unlist(extra_codes))
write_csv(data.frame(all_codes),"../outputs/infection_death_codelist.csv")

# get deaths 

deaths_due_to_infection = disease_codes %>% filter(!is.na(date_of_death_f40000_0_0)) %>% 
  dplyr::select(EID78867,age_at_death_f40007_0_0,contains("underlying_primary_cause_of_death_icd10_f40001"),contains("contributory_secondary_ca_f40002")) %>% 
  pivot_longer(cols = -c(1,2)) %>% 
  filter(!is.na(value)) %>%
  mutate(expanded_code = ifelse(nchar(value)==3,paste0(value,"X"),value)) %>%
  mutate(is_infection = ifelse(expanded_code %in% all_codes,"yes","no")) %>% 
  filter(is_infection=="yes")


# demographics 
## clean ethnicity data 
dat = dat %>%
  mutate(ethnicity_clean = case_when(
    ethnic_background_f21000_0_0 %in% c("White","British","Irish") ~ "White",
    ethnic_background_f21000_0_0 %in% c("Asian or Asian British","Chinese","Pakistani","Bangladeshi","Any other Asian background") ~ "Asian",
    ethnic_background_f21000_0_0 %in% c("Black or Black British","Caribbean","African","Any other Black background") ~ "Black",
    .default = as.character("Mixed/Other/Missing")
  ))
    
## clean smoking data 
dat = dat %>%
  mutate(smoking_status_f20116_0_0 = case_when(
    smoking_status_f20116_0_0 == "Previous" ~ "Ever",
    smoking_status_f20116_0_0 == "Current" ~ "Ever",
    smoking_status_f20116_0_0 == "Never" ~ "Never"    
  ))

# define censoring details
dat = dat %>% 
  mutate(death = ifelse(!is.na(age_at_death_f40007_0_0),"dead","censored")) %>%
  mutate(data_extraction_date = as.Date("2023-12-15",format="%Y-%m-%d")) %>%
  mutate(pseudo_dob = as.Date(paste0(year_of_birth_f34_0_0,"-01-01"),format="%Y-%m-%d")) %>%
  mutate(age_at_data_extraction = as.numeric(data_extraction_date - pseudo_dob ) / 365.25 ) %>%
  mutate(age_at_loss_to_fu = as.numeric(date_lost_to_followup_f191_0_0 - pseudo_dob ) / 365.25 ) %>%
  mutate(censor_age = case_when(
    !is.na(age_at_death_f40007_0_0) ~ age_at_death_f40007_0_0,
    !is.na(date_lost_to_followup_f191_0_0) ~ age_at_loss_to_fu,
    is.na(date_lost_to_followup_f191_0_0) & is.na(age_at_death_f40007_0_0) ~ age_at_data_extraction
  )) %>% 
  mutate(years_of_fu = censor_age - age_at_recruitment_f21022_0_0)

# define MS status
dat = dat %>% mutate(MS_status = ifelse(!is.na(date_g35_first_reported_multiple_sclerosis_f131042_0_0),"MS","Control"))

# filter to those with serology
dat = dat %>% filter(!is.na(ebv_seropositivity_for_epsteinbarr_virus_f23053_0_0))
table(dat$ebv_seropositivity_for_epsteinbarr_virus_f23053_0_0,dat$MS_status)
dat %>%
  group_by(ebv_seropositivity_for_epsteinbarr_virus_f23053_0_0) %>%
  dplyr::count(MS_status) %>%
  mutate(prop = n /sum(n))
dat %>%
  group_by(MS_status) %>%
  dplyr::count(ebv_seropositivity_for_epsteinbarr_virus_f23053_0_0) %>%
  mutate(prop = n /sum(n), total =sum(n))

fisher.test(dat$ebv_seropositivity_for_epsteinbarr_virus_f23053_0_0,dat$MS_status)

# define infection death
dat = dat %>% 
  mutate(death_due_to_infection = ifelse(EID78867 %in% deaths_due_to_infection$EID78867,"Yes","No"))
dat$death_due_to_infection %>% table


library(compareGroups)

# demographics 
dat %>% nrow()
dat %>% dplyr::count(ebv_seropositivity_for_epsteinbarr_virus_f23053_0_0) %>% mutate(prop = n/sum(n))
tbl = compareGroups::compareGroups(
  data = dat,
  ebv_seropositivity_for_epsteinbarr_virus_f23053_0_0 ~ 
    age_at_recruitment_f21022_0_0 + 
    year_of_birth_f34_0_0 + 
    sex_f31_0_0 + 
    townsend_deprivation_index_at_recruitment_f22189_0_0 + 
    body_mass_index_bmi_f21001_0_0 +
    smoking_status_f20116_0_0 +
    ethnicity_clean + 
    death + 
    censor_age + 
    years_of_fu + 
    death_due_to_infection +  
    white_blood_cell_leukocyte_count_f30000_0_0 +
    lymphocyte_count_f30120_0_0 +
    monocyte_count_f30130_0_0 +
    neutrophill_count_f30140_0_0 +
    eosinophill_count_f30150_0_0 +
    basophill_count_f30160_0_0 + 
    ebna1_antigen_for_epsteinbarr_virus_f23004_0_0 +
    vca_p18_antigen_for_epsteinbarr_virus_f23003_0_0 +
    ead_antigen_for_epsteinbarr_virus_f23006_0_0 +
    zebra_antigen_for_epsteinbarr_virus_f23005_0_0,
  method = c(2,2,3,2,2,3,3,3,2,2,3,2,2,2,2,2,2,3,2,2,2,2)
    )
compareGroups::createTable(tbl) %>% export2csv(file = "../outputs/demographics.csv")

# cross-sectional association with all disease codes 
disease_codes = disease_codes %>% filter(EID78867 %in% dat$EID78867 )
just_source_of_report_codes = disease_codes %>% dplyr::select(EID78867,contains("first_reported")) 

# loop through each outcome 
res = list()
for(i in c(2:ncol(just_source_of_report_codes))){
  message(i)
  this_outcome = just_source_of_report_codes %>% dplyr::select(1,all_of(i))
  this_disease = colnames(just_source_of_report_codes)[i]
  colnames(this_outcome)[2] = "disease"
  this_outcome = this_outcome %>% mutate(disease_binary = ifelse(is.na(disease),0,1))
  counts = this_outcome %>% filter(disease_binary == 1) %>% nrow()
  if(counts > 100){

    model_dat = dat %>% left_join(this_outcome,by="EID78867")

    # define prevalent disease 
    model_dat = model_dat %>%
      mutate(date_of_recruitment = as.Date(paste0("01-01-",year_of_birth_f34_0_0+age_at_recruitment_f21022_0_0),format="%d-%m-%Y")) %>%
      mutate(prevalent_disease = ifelse(!is.na(disease) & disease < (date_of_recruitment + 5*365.25),1,0))

    cc_counts = model_dat %>% 
      group_by(ebv_seropositivity_for_epsteinbarr_virus_f23053_0_0) %>% dplyr::count(prevalent_disease,.drop = F) %>%
      filter(!is.na(prevalent_disease)) %>% 
      pivot_wider(names_from = c(1,2),values_from = n)
    cc_counts$False_1 = ifelse(is.null(cc_counts$False_1),0,cc_counts$False_1)
    cc_counts$False_0 = ifelse(is.null(cc_counts$False_0),0,cc_counts$False_0)
    cc_counts$True_0 = ifelse(is.null(cc_counts$True_0),0,cc_counts$True_0)
    cc_counts$True_1 = ifelse(is.null(cc_counts$True_1),0,cc_counts$True_1)
    
    cc_counts = cc_counts %>% dplyr::rename("EBV -ve controls" = False_0,
                    "EBV -ve cases" = False_1,
                    "EBV +ve controls" = True_0,
                    "EBV +ve cases" = True_1)
    
    res[[length(res)+1]] = glm(data = model_dat, prevalent_disease ~ 
    townsend_deprivation_index_at_recruitment_f22189_0_0 + 
    age_at_recruitment_f21022_0_0 + as.character(sex_f31_0_0) + 
    as.character(ebv_seropositivity_for_epsteinbarr_virus_f23053_0_0), family = binomial(link="logit")) %>% 
      broom::tidy() %>% 
      mutate(model_outcome = this_disease) %>% 
      mutate(total_n = counts,
             OR = exp(estimate),
             lower_ci = exp(estimate - 1.96*std.error),
             upper_ci = exp(estimate + 1.96*std.error),
             ) %>% 
      bind_cols(cc_counts) %>% filter(grepl("ebv",term))%>% 
      mutate(crude_OR = (
        `EBV +ve cases` / `EBV +ve controls`
        ) / 
          (`EBV -ve cases` / `EBV -ve controls`)
      )
  }
}
res = do.call("bind_rows",res)
res = res %>% dplyr::select(-term)

# remove sex-specific outcomes 

res$model_outcome = str_remove_all(res$model_outcome,"date_")

res = res %>% filter(
  !grepl("^o",model_outcome) &
    !grepl("^n4",model_outcome) &
    !grepl("^n5",model_outcome) &
    !grepl("^n7",model_outcome) &
    !grepl("^n8",model_outcome) &
    !grepl("^n90",model_outcome) &
    !grepl("^n91",model_outcome) &
    !grepl("^n92",model_outcome) &
    !grepl("^n93",model_outcome) &
    !grepl("^n94",model_outcome) &
    !grepl("^n95",model_outcome) &
    !grepl("^n96",model_outcome) &
    !grepl("^n97",model_outcome) &
    !grepl("^n98",model_outcome)
)
res$model_outcome = str_remove_all(res$model_outcome,"first_reported_")

res = res %>% arrange(p.value)
n_tests = res %>% nrow()
bonf = 0.05/n_tests
res = res %>% mutate(fdr = p.adjust(p.value,method="fdr"))
write_csv(res,"../outputs/cross_sectional_phewas.csv")

# build matched dataset 
dat$townsend_decile = Hmisc::cut2(dat$townsend_deprivation_index_at_recruitment_f22189_0_0,g=10)
control_df = dat %>% filter(ebv_seropositivity_for_epsteinbarr_virus_f23053_0_0 == "True")
case_df = dat %>% filter(ebv_seropositivity_for_epsteinbarr_virus_f23053_0_0 != "True")

set.seed(123)
matched_list = list()
for(i in c(1:nrow(case_df))){
  message(i)
  # get thie case
  this_case = case_df[i,]

  # match on age & sex & Townsend
  this_matched_control = control_df %>%
    filter(!EID %in% matched_list) %>%
    filter(sex_f31_0_0 == this_case$sex_f31_0_0) %>%
    filter(townsend_decile == this_case$townsend_decile) %>%
    mutate(delta_age = abs(age_at_recruitment_f21022_0_0 - this_case$age_at_recruitment_f21022_0_0)) %>% 
    slice_min(delta_age,with_ties=F,n=1)
    

  matched_list[[i]] = this_matched_control$EID78867
  control_df <<- control_df %>% filter(!EID78867 %in% matched_list)
}
matched_list = unlist(matched_list)
overall_matched_df = dat %>% filter(EID78867 %in% c(matched_list,case_df$EID78867))

# save 
saveRDS(overall_matched_df,"../outputs/matched_dataset.rds")
saveRDS(dat,"../outputs/main_dataset.rds")

# read back in 
library(tidyverse)
setwd("/data/Wolfson-PNU-dementia/UKB_Projects/Piranavi/scripts/")

overall_matched_df = readRDS("../outputs/matched_dataset.rds")
dat = readRDS("../outputs/main_dataset.rds")

# check serostatus looks correct
dat$ebna1_antigen_for_epsteinbarr_virus_f23004_0_0
dat$vca_p18_antigen_for_epsteinbarr_virus_f23003_0_0
dat$zebra_antigen_for_epsteinbarr_virus_f23005_0_0
dat$zebra_antigen_for_epsteinbarr_virus_f23005_0_0
dat$ead_antigen_for_epsteinbarr_virus_f23006_0_0
ebv_serology_only = dat %>% dplyr::select(contains("antigen_for_epstein"))
GGally::ggpairs(ebv_serology_only)

# strict seronegativity definition
dat = dat %>% 
  mutate(strict_ebv_seronegativity = ifelse(
    ead_antigen_for_epsteinbarr_virus_f23006_0_0 <100 & 
    zebra_antigen_for_epsteinbarr_virus_f23005_0_0 <100 & 
    vca_p18_antigen_for_epsteinbarr_virus_f23003_0_0 < 250 & 
    ebna1_antigen_for_epsteinbarr_virus_f23004_0_0 < 250,
    "seronegative","seropositive"
  ))
dat %>% dplyr::count(strict_ebv_seronegativity) %>% mutate(prop = n/sum(n))
dat %>% group_by(MS_status) %>% dplyr::count(strict_ebv_seronegativity) %>% mutate(prop = n/sum(n))

# survival analysis 

library(survival)
library(survminer)
dat$years_of_fu %>% median
dat$years_of_fu %>% sum
dat %>% nrow()
dat %>% dplyr::count(death) %>% mutate(prop = n/sum(n))
n_deaths = dat %>% filter(death=="dead") %>% nrow()
n_infection_deaths = dat %>% filter(death_due_to_infection=="Yes") %>% nrow()
n_person_years = dat$years_of_fu %>% sum
message("deaths per 1,000 person-years:", n_deaths/n_person_years * 1000)
message("infection deaths per 1,000 person-years:", n_infection_deaths/n_person_years * 1000)

dat %>% dplyr::count(death_due_to_infection) %>% mutate(prop = n/sum(n))
dat %>% filter(death=="dead") %>% dplyr::count(death_due_to_infection) %>% mutate(prop = n/sum(n))
dat %>% group_by(ebv_seropositivity_for_epsteinbarr_virus_f23053_0_0) %>% dplyr::count(death_due_to_infection) %>% mutate(prop = n/sum(n))
counts = dat %>% group_by(ebv_seropositivity_for_epsteinbarr_virus_f23053_0_0) %>% dplyr::count(death_due_to_infection) %>% filter(death_due_to_infection=="Yes")
fu_time = dat %>% group_by(ebv_seropositivity_for_epsteinbarr_virus_f23053_0_0) %>% summarise(fu_time = sum(years_of_fu))
fu_time = fu_time %>% inner_join(counts,by="ebv_seropositivity_for_epsteinbarr_virus_f23053_0_0")
fu_time %>% mutate(mortality_rate = n/fu_time*1000)

# curate outcomes
dat$death_outcome <- ifelse(dat$death_due_to_infection == "Yes", 1, 0)
dat$survtime <- dat$years_of_fu

# build survival object
dat$ebv_simple = ifelse(dat$ebv_seropositivity_for_epsteinbarr_virus_f23053_0_0=="True","EBV +ve","EBV -ve")
fit <- survfit(Surv(survtime, death_outcome) ~ ebv_simple, data = dat)


# plot with ggsurvplot
p = ggsurvplot(
  fit,
  data = dat,
  risk.table = T,
  pval = F,
  fun = "pct",
  conf.int = TRUE,
  ylim = c(95,100),
  xlab = "Time from recruitment (years)",
  ylab = "Survival probability \n(% free from infection-related death)",
  legend.title = "EBV serostatus",
  conf.int.alpha = 0.25,             # lighter confidence bands
  legend.labs = c("EBV +ve", "EBV -ve"),
  palette = "Accent"
)
p = p$plot 
p = p + theme(
    text = element_text(size = 10),
    axis.title = element_text(size = 10),
    axis.text = element_text(size = 9),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9),
    plot.title = element_text(size = 11),
    plot.subtitle = element_text(size = 10),
    plot.caption = element_text(size = 8)
  )

# cox models

## strict definition 
strict_def = coxph(Surv(survtime, death_outcome) ~ 
        as.character(strict_ebv_seronegativity) + 
        as.character(sex_f31_0_0) + age_at_recruitment_f21022_0_0 + 
        townsend_deprivation_index_at_recruitment_f22189_0_0, data = dat) %>% broom::tidy() %>% 
  mutate(HR = exp(estimate), lower_ci = exp(estimate - 1.96*std.error), upper_ci = exp(estimate + 1.96*std.error)) %>%
  mutate(model = "Strict def")


## primary model
all_res = coxph(Surv(survtime, death_outcome) ~ 
        as.character(ebv_seropositivity_for_epsteinbarr_virus_f23053_0_0) + 
        as.character(sex_f31_0_0) + age_at_recruitment_f21022_0_0 + 
        townsend_deprivation_index_at_recruitment_f22189_0_0, data = dat) %>% broom::tidy() %>% 
  mutate(HR = exp(estimate), lower_ci = exp(estimate - 1.96*std.error), upper_ci = exp(estimate + 1.96*std.error)) %>%
  mutate(model = "Primary analysis")

## model w bmi + smoking
dat$smoking_status_f20116_0_0 = relevel(factor(dat$smoking_status_f20116_0_0),ref="Never")
all_res = all_res %>% 
bind_rows(
coxph(Surv(survtime, death_outcome) ~ 
        as.character(ebv_seropositivity_for_epsteinbarr_virus_f23053_0_0) + 
        as.character(sex_f31_0_0) + age_at_recruitment_f21022_0_0 + 
        townsend_deprivation_index_at_recruitment_f22189_0_0 +
        body_mass_index_bmi_f21001_0_0 +
        smoking_status_f20116_0_0, data = dat) %>% broom::tidy() %>% 
  mutate(HR = exp(estimate), lower_ci = exp(estimate - 1.96*std.error), upper_ci = exp(estimate + 1.96*std.error))%>%
  mutate(model = "Full model (+ Smoking + BMI)")
)

## matched 
overall_matched_df$death_outcome <- ifelse(overall_matched_df$death_due_to_infection == "Yes", 1, 0)
overall_matched_df$survtime <- overall_matched_df$years_of_fu

# build survival object
overall_matched_df$ebv_simple = ifelse(overall_matched_df$ebv_seropositivity_for_epsteinbarr_virus_f23053_0_0=="True","EBV +ve","EBV -ve")
all_res = all_res %>% 
bind_rows(
  coxph(Surv(survtime, death_outcome) ~ 
        as.character(ebv_seropositivity_for_epsteinbarr_virus_f23053_0_0) + 
        as.character(sex_f31_0_0) + age_at_recruitment_f21022_0_0 + 
        townsend_deprivation_index_at_recruitment_f22189_0_0, data = overall_matched_df) %>% broom::tidy() %>% 
  mutate(HR = exp(estimate), lower_ci = exp(estimate - 1.96*std.error), upper_ci = exp(estimate + 1.96*std.error))%>%
  mutate(model = "Matched cohort")
)

## ebna1
all_res = all_res %>% 
bind_rows(
coxph(Surv(survtime, death_outcome) ~ 
        ebna1_antigen_for_epsteinbarr_virus_f23004_0_0 + 
        as.character(sex_f31_0_0) + age_at_recruitment_f21022_0_0 + 
        townsend_deprivation_index_at_recruitment_f22189_0_0, data = dat) %>% broom::tidy() %>% 
  mutate(HR = exp(estimate), lower_ci = exp(estimate - 1.96*std.error), upper_ci = exp(estimate + 1.96*std.error))%>%
  mutate(model = "EBNA1 titre")
)

## vca
all_res = all_res %>% 
bind_rows(
coxph(Surv(survtime, death_outcome) ~ 
        vca_p18_antigen_for_epsteinbarr_virus_f23003_0_0 + 
        as.character(sex_f31_0_0) + age_at_recruitment_f21022_0_0 + 
        townsend_deprivation_index_at_recruitment_f22189_0_0, data = dat) %>% broom::tidy() %>% 
  mutate(HR = exp(estimate), lower_ci = exp(estimate - 1.96*std.error), upper_ci = exp(estimate + 1.96*std.error))%>%
  mutate(model = "VCA titre")
)

write_csv(all_res,"../outputs/all_cox_res.csv")



# forest 
primary = all_res %>% filter(model=="Full model (+ Smoking + BMI)") %>%
  mutate(term = case_when(
    term == "as.character(sex_f31_0_0)Male" ~ "Male gender",
    term == "as.character(ebv_seropositivity_for_epsteinbarr_virus_f23053_0_0)True" ~ "EBV +",
    term == "age_at_recruitment_f21022_0_0" ~ "Age",
    term == "townsend_deprivation_index_at_recruitment_f22189_0_0" ~ "Deprivation",
    term == "body_mass_index_bmi_f21001_0_0" ~ "BMI",
    term == "smoking_status_f20116_0_0Ever" ~ "Smoking"
  ))
primary = primary %>% arrange(desc(estimate))
primary$term = factor(primary$term,levels = primary$term,ordered=T) 
primary$p_sig = ifelse(primary$p.value<0.05,"*","")
p1 = ggplot(primary,aes(HR,term,fill=term,label=p_sig))+
  geom_errorbarh(mapping = aes(xmin=lower_ci,xmax=upper_ci,y=term),height=0.1)+
  geom_point(shape=21,size=3)+
  theme_classic()+
  geom_text(mapping = aes(x=upper_ci*1.1))+
  labs(x="Hazard ratio for\n infection-related death",y="Variable")+
  geom_vline(xintercept=1,linetype="dashed",color="black",alpha=0.5)+
  scale_x_log10(breaks = c(0.5,1,2,4),limits = c(0.4,5))+
  scale_fill_brewer(palette="Blues")+
  theme(legend.position="none", 
    text = element_text(size = 10),
    axis.title = element_text(size = 10),
    axis.text = element_text(size = 9),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9),
    plot.title = element_text(size = 11),
    plot.subtitle = element_text(size = 10),
    plot.caption = element_text(size = 8)
  )


png("../outputs/survival_curves.png",res=900,units="in",width=12,height=4)
cowplot::plot_grid(p,p1,align = "h")
dev.off()


