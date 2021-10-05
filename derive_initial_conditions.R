library(data.table)
library(qs)
library(ggplot2)

source("./backcalculation_functions.R")


# 
# SETUP
# 

# Load backcalculation output
load("./output/2021-09-30/backcalculation_output.RData")

# # Source script with delay functions
# cm_path = "./covidm_for_fitting/"
# source(paste0(cm_path,"/R/shared/cmS_misc.R"))

theme_set(cowplot::theme_cowplot(font_size = 10) + theme(strip.background = element_blank()))

# Set age groups to be the same as in deconvolution for the time being
agegroups_model = factor(agegroups,levels = agegroups)
min_ages_model = get_min_age(agegroups_model)

# Read in UN population data
pop = qread("../un_data/unwpp_data.qs")

#
# DISAGGREGATION OF INFECTIONS
#


# Make data table for incidence
inc_dt = CJ(country=backcalc_dt[,unique(country)],age=0:100,date=backcalc_dt[,unique(date)])
inc_dt = merge(inc_dt,pop[,.(country,age,population)],by=c("country","age"),all.x=T)
inc_dt[,age_group_model:=cut(age,c(min_ages_model,Inf),labels=agegroups_model,right=F)]
inc_dt = inc_dt[,.(population=sum(population)),by=.(country,age_group_model,date)]

# Add deconvolution age groups to data table
inc_dt[,age_group:=cut(get_min_age(age_group_model),c(min_ages,Inf),labels=agegroups,right=F)]

# Merge with deconvolution output
# N.B. duplicates exposures within same age_group
inc_dt = merge(inc_dt,backcalc_dt[,.(country,age_group,date,exposures)],by=c("country","age_group","date"))

# Split exposures by population fraction
inc_dt[,exposures:=exposures*population/sum(population),by=.(country,age_group,date)]

# Plot
ggplot(inc_dt,aes(x=date,y=exposures,group=age_group_model,color=age_group_model)) +
    geom_line() +
    facet_wrap(~country)

# Read in age-dependent symptomatic fraction
covid_scenario = qread(datapath("2-linelist_both_fit_fIa0.5-rbzvih.qs"));
colsy = names(covid_scenario)[grep("y_",names(covid_scenario))]
covy = colMeans(covid_scenario[,..colsy])

min_ages_y = as.integer(sub("y_","",names(covy)))
covy = data.table(age_group=c(paste0(min_ages_y[1:(length(min_ages_y)-1)],"-",min_ages_y[2:length(min_ages_y)]-1),paste0(min_ages_y[length(min_ages_y)],"+")),y=covy)
# Make data table for age-dependent symptomatic fraction
covy_dt = CJ(country=backcalc_dt[,unique(country)],age=0:100)
# Merge with population data
covy_dt = merge(covy_dt,pop[,.(country,age,population)],by=c("country","age"),all.x=T)
covy_dt[,age_group:=cut(age,c(min_ages_y,Inf),labels=covy[,age_group],right=F)]
# Merge with symptomatic fraction
covy_dt = merge(covy_dt,covy,by="age_group",all.x=T)
# Change age groups
covy_dt[,age_group_model:=cut(age,c(min_ages_model,Inf),labels=agegroups_model,right=F)]
# Calculate population-weighted average symptomatic fraction
covy_dt = covy_dt[,.(y=sum(y*population)/sum(population)),by=.(country,age_group_model)]

# Plot
ggplot(covy_dt,aes(x=age_group_model,y=y,group=country,color=country)) + geom_line()

# ADD VACCINATION NUMBERS
vax_dt = construct_vax_data_table(vax,inc_dt[,unique(date)]-Ab_delay,agegroups_model,pop)
vax_dt[,date_v:=date]
vax_dt[,date:=date+Ab_delay]
vax_dt_wide = dcast(vax_dt,country + date + age_group + population ~ type,value.var = c("first","second","prop","cum_prop"))

setnames(vax_dt_wide,c("age_group","first_va","first_vb","second_va","second_vb"),c("age_group_model","nS_Va1","nS_Vb1","nVa1_Va2","nVb1_Vb2"))


#
# INITIAL CONDITIONS CALCULATION
#


prev_dt = calc_init_condns(inc_dt,vax_dt_wide,agegroups_model,covy_dt,vrnt_prop,ve_params,dE,dIp,dIs,dIa)

# Plot to check
ggplot(prev_dt,aes(x=date,y=nS_E+nV_E+nV_L,color=age_group_model)) +
    geom_line() +
    facet_wrap(~country)

rm(backcalc_dt,backcalc_samps,backcalc_dt_non_ltc,backcalc_samps_non_ltc,
   backcalc_dt_ltc,backcalc_samps_ltc)

save.image(paste0(dir_out,"init_cond_output.RData"))
