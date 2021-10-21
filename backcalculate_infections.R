library(data.table)
library(qs)
library(ggplot2)
library(covidregionaldata)
library(ISOweek)
library(lubridate)
library(surveillance)
if(!require(incidental)){
    install.packages("incidental")
}
library(incidental)
library(mgcv)
library(nnet)
library(splines)
library(effects)
library(doParallel)
library(cowplot)
library(stringr)

source("./backcalculation_functions.R")
# source("./initial_conditions_functions.R")
# source("./simulation_functions.R")

# TODO:
# - sort out alignment of dates for deaths data and vaccination data - currently
# only works because max date in vax data + Ab_delay is after max date in deaths
# data - DONE
# - update to latest data - DOING
# - switch delay distributions to those in covidm - DONE
# - use country names or two-letter country iso codes throughout - DONE (using country names)
# - use age-dep vaccination data from OWID for countries missing age breakdown
# - resolve remaining issues with death time series, e.g. where data is missing
# scale cumulative deaths by distribution of total deaths over that time period - DONE
# - combine data for England, Wales and Scotland - INCLUDED ENGLAND
# - separate into LTC and non-LTC deaths - DONE, NEED TO CORRECT
#   - don't apply split to France data as it already excludes care home deaths - DONE
#   - use mean proportion of deaths that are among LTC residents for countries
#   missing LTC death data - DONE
# - use sampling approach to get integer numbers of deaths rather than 
# multiplying by proportions and ending up with decimals (so that data is in 
# correct format for RIDE)?
# - make proportion of vaccine type time dependent - DONE
# - include time-dependent variant proportions in vaccine efficacy - DONE
# - use country-specific IFRs?
# - correct dates for vaccination data to be end of ISO week - no, start of ISO week is fine
# - sort out Finland data - interpolate with WHO data? - DONE
# - review/correct proportions for LTC deaths, esp. France - DONE
# - include multiple vaccine doses
# - fix issue with early vaccination data for 70-79yo for Norway?
# - shift infection time series when plotting against seroprevalence data to account for delay from infection to seroconversion (~20-21 days)
# - download and save data automatically in script rather than manually
# - fix issue with disaggregation of deaths in raw data for Ireland if it's included
# - fix issue with IFR for Iceland (from jump in ei due to missing vaccine type data for 80+ year olds) if it's included
# - rerun to correct mistake with England vax schedule data


# 
# SETUP
# 


# Register parallel backend
# registerDoParallel(cores = detectCores()-1)
registerDoParallel(cores = 4)

date_fitting = today()

# Create directory to save output into
dir_out = paste0("./output/",date_fitting,"/")
dir.create(dir_out,recursive = T)

# Set plot theme
theme_set(cowplot::theme_cowplot(font_size = 10) + theme(strip.background = element_blank()))

# Set age groups
agegroups = c("0-39","40-49","50-59","60-69","70-79","80+") #c("0-9","10-19","20-29","30-39","40-49","50-59","60-69","70-79","80+")
min_ages = get_min_age(agegroups)
max_ages = get_max_age(agegroups)
max_ages[is.na(max_ages)] = Inf

# Get country iso codes
covid_data_path = "./fitting_data/"
datapath = function(x) paste0(covid_data_path, x)
country_iso_codes = readRDS(datapath("country_iso_codes.rds"))

# Get age groups from covidm
cm_path = "./covidm_for_fitting"
cm_build = F
cm_version = 2
source(paste0(cm_path,"/R/covidm_new.R"))
popUK = readRDS(datapath("popNHS.rds"))
matricesUK = readRDS(datapath("matricesNHS.rds"))
cm_populations = rbind(cm_populations[name != "United Kingdom"], popUK)
cm_matrices = c(cm_matrices, matricesUK)

# Get names of countries that are in cm_matrices
regions = cm_populations[location_type==4 & name %in% names(cm_matrices),as.character(unique(name))]

# Set delay distributions:
# latent period
dE  = cm_delay_gamma(2.5, 2.5, t_max = 15, t_step = 1)$p
# presymptomatic period
dIp = cm_delay_gamma(2.5, 4.0, t_max = 15, t_step = 1)$p
# symptomatic period
dIs = cm_delay_gamma(2.5, 4.0, t_max = 15, t_step = 1)$p
# asymptomatic period
dIa = cm_delay_gamma(5.0, 4.0, t_max = 15, t_step = 1)$p
# infection-to-hospitalisation delay
dHosp = cm_delay_gamma(6.0 + 2.5, 0.71, t_max = 60, t_step = 1)$p
# infection-to-death delay
dDeath = cm_delay_lnorm(15, 0.9, t_max = 60, t_step = 1)$p

# Build parameters for different regions ###
params = cm_parameters_SEI3R(regions, deterministic = T, 
                             date_start = "2020-01-01", 
                             date_end = date_fitting,
                             dE  = dE,
                             dIp = dIp,
                             dIs = dIs,
                             dIa = dIa)
params = cm_split_matrices_ex_in(params, 15)

# Get age-dependent susceptibility and symptomatic fraction
covid_scenario = qread(datapath("2-linelist_both_fit_fIa0.5-rbzvih.qs"));
colsu = names(covid_scenario)[grep("u_",names(covid_scenario))]
covu = unname(rep(colMeans(covid_scenario[,..colsu]), each = 2))
colsy = names(covid_scenario)[grep("y_",names(covid_scenario))]
covy = unname(rep(colMeans(covid_scenario[,..colsy]), each = 2))

for (i in seq_along(params$pop)) {
    params$pop[[i]]$u = 0.08*covu / mean(covu);
    params$pop[[i]]$u2 = 0.08*covu / mean(covu);
    params$pop[[i]]$u3 = 0.08*covu / mean(covu);
    params$pop[[i]]$y = covy;
    params$pop[[i]]$y2 = covy;
    params$pop[[i]]$y3 = covy;
}

# get list of age groups in covidm
agegroups_model = params$pop[[1]]$group_names
agegroups_model = factor(agegroups_model,levels = agegroups_model)
min_ages_model = get_min_age(agegroups_model)


# Set vaccine efficacy parameters
ve_params = list()

# Infection:
# eiX_vYZ = efficacy of dose Z of vaccine Y against infection with strain X
# AstraZeneca
# Dose 1
ve_params$ei_va1  = VE(0.70)
ve_params$ei2_va1 = VE(0.70)
ve_params$ei3_va1 = VE(0.43)
# Dose 2
ve_params$ei_va2  = VE(0.75)
ve_params$ei2_va2 = VE(0.75)
ve_params$ei3_va2 = VE(0.63)

# Pfizer / Moderna
# Dose 1
ve_params$ei_vb1  = VE(0.70)
ve_params$ei2_vb1 = VE(0.70)
ve_params$ei3_vb1 = VE(0.62)
# Dose 2
ve_params$ei_vb2  = VE(0.85)
ve_params$ei2_vb2 = VE(0.85)
ve_params$ei3_vb2 = VE(0.80)

# Symptomatic disease:
# ed_vYZiX = efficacy of dose Z of vaccine Y against disease given infection with strain X
# AstraZeneca
# Dose 1
ve_params$ed_va1i  = VE(0.70,   ve_params$ei_va1)
ve_params$ed_va1i2 = VE(0.70,   ve_params$ei2_va1)
ve_params$ed_va1i3 = VE(0.48,   ve_params$ei3_va1)
# Dose 2
ve_params$ed_va2i  = VE(0.80,   ve_params$ei_va2)
ve_params$ed_va2i2 = VE(0.80,   ve_params$ei2_va2)
ve_params$ed_va2i3 = VE(0.65,   ve_params$ei3_va2)

# Pfizer / Moderna
# Dose 1
ve_params$ed_vb1i  = VE(0.70,   ve_params$ei_vb1)
ve_params$ed_vb1i2 = VE(0.70,   ve_params$ei2_vb1)
ve_params$ed_vb1i3 = VE(0.62,   ve_params$ei3_vb1)
# Dose 2
ve_params$ed_vb2i  = VE(0.90,   ve_params$ei_vb2)
ve_params$ed_vb2i2 = VE(0.90,   ve_params$ei2_vb2)
ve_params$ed_vb2i3 = VE(0.81,   ve_params$ei3_vb2)

# Hospitalisation:
# eh_vYZiX = efficacy of dose Z of vaccine Y against hospitalisation given disease from strain X
# AstraZeneca
# Dose 1
ve_params$eh_va1d  = VE(0.85, 1 - (1 - ve_params$ei_va1)  * (1 - ve_params$ed_va1i))
ve_params$eh_va1d2 = VE(0.85, 1 - (1 - ve_params$ei2_va1) * (1 - ve_params$ed_va1i2))
ve_params$eh_va1d3 = VE(0.83, 1 - (1 - ve_params$ei3_va1) * (1 - ve_params$ed_va1i3))
# Dose 2
ve_params$eh_va2d  = VE(0.90, 1 - (1 - ve_params$ei_va2)  * (1 - ve_params$ed_va2i))
ve_params$eh_va2d2 = VE(0.90, 1 - (1 - ve_params$ei2_va2) * (1 - ve_params$ed_va2i2))
ve_params$eh_va2d3 = VE(0.95, 1 - (1 - ve_params$ei3_va2) * (1 - ve_params$ed_va2i3))

# Pfizer / Moderna
# Dose 1
ve_params$eh_vb1d  = VE(0.85, 1 - (1 - ve_params$ei_vb1)  * (1 - ve_params$ed_vb1i))
ve_params$eh_vb1d2 = VE(0.85, 1 - (1 - ve_params$ei2_vb1) * (1 - ve_params$ed_vb1i2))
ve_params$eh_vb1d3 = VE(0.92, 1 - (1 - ve_params$ei3_vb1) * (1 - ve_params$ed_vb1i3))
# Dose 2
ve_params$eh_vb2d  = VE(0.95, 1 - (1 - ve_params$ei_vb2)  * (1 - ve_params$ed_vb2i))
ve_params$eh_vb2d2 = VE(0.95, 1 - (1 - ve_params$ei2_vb2) * (1 - ve_params$ed_vb2i2))
ve_params$eh_vb2d3 = VE(0.97, 1 - (1 - ve_params$ei3_vb2) * (1 - ve_params$ed_vb2i3))

# Mortality:
# em_vYZdX = efficacy of dose Z of vaccine Y against death given disease from strain X
# AstraZeneca
# Dose 1
ve_params$em_va1d  = VE(0.85, 1 - (1 - ve_params$ei_va1)  * (1 - ve_params$ed_va1i))
ve_params$em_va1d2 = VE(0.85, 1 - (1 - ve_params$ei2_va1) * (1 - ve_params$ed_va1i2))
ve_params$em_va1d3 = VE(0.83, 1 - (1 - ve_params$ei3_va1) * (1 - ve_params$ed_va1i3))
# Dose 2
ve_params$em_va2d  = VE(0.95, 1 - (1 - ve_params$ei_va2)  * (1 - ve_params$ed_va2i))
ve_params$em_va2d2 = VE(0.95, 1 - (1 - ve_params$ei2_va2) * (1 - ve_params$ed_va2i2))
ve_params$em_va2d3 = VE(0.95, 1 - (1 - ve_params$ei3_va2) * (1 - ve_params$ed_va2i3))

# Pfizer / Moderna
# Dose 1
ve_params$em_vb1d  = VE(0.85, 1 - (1 - ve_params$ei_vb1)  * (1 - ve_params$ed_vb1i))
ve_params$em_vb1d2 = VE(0.85, 1 - (1 - ve_params$ei2_vb1) * (1 - ve_params$ed_vb1i2))
ve_params$em_vb1d3 = VE(0.92, 1 - (1 - ve_params$ei3_vb1) * (1 - ve_params$ed_vb1i3))
# Dose 2
ve_params$em_vb2d  = VE(0.95, 1 - (1 - ve_params$ei_vb2)  * (1 - ve_params$ed_vb2i))
ve_params$em_vb2d2 = VE(0.95, 1 - (1 - ve_params$ei2_vb2) * (1 - ve_params$ed_vb2i2))
ve_params$em_vb2d3 = VE(0.97, 1 - (1 - ve_params$ei3_vb2) * (1 - ve_params$ed_vb2i3))

# Onward transmission:
# et_vYZdX = efficacy of dose Z of vaccine Y against onward transmission given infection with strain X
# AstraZeneca
# Dose 1
ve_params$et_va1i  = VE(0.47)
ve_params$et_va1i2 = VE(0.47)
ve_params$et_va1i3 = VE(0.42)
# Dose 2
ve_params$et_va2i  = VE(0.47)
ve_params$et_va2i2 = VE(0.47)
ve_params$et_va2i3 = VE(0.42)

# Pfizer / Moderna
# Dose 1
ve_params$et_vb1i  = VE(0.47)
ve_params$et_vb1i2 = VE(0.47)
ve_params$et_vb1i3 = VE(0.42)
# Dose 2
ve_params$et_vb2i  = VE(0.47)
ve_params$et_vb2i2 = VE(0.47)
ve_params$et_vb2i3 = VE(0.42)

# Set delays from vaccination to development of Ab
Ab_delay1 = 28 # delay after 1st dose in days
Ab_delay2 = 14 # delay after 2nd dose in days


# 
# DATA 
# 


# POPULATION DATA

# Read in UN population data
pop = qread("../un_data/unwpp_data.qs")
# Filter to 2020 estimates
pop = pop[year==2020]
setnames(pop,"total","population")

# Change some country names to match deaths data
pop[country=="Bolivia (Plurinational State of)",country:="Bolivia"]
pop[country=="China, Taiwan Province of China",country:="Taiwan"]
pop[country=="Dem. People's Republic of Korea",country:="South Korea"]
pop[country=="Republic of Moldova",country:="Moldova"]
pop[country=="State of Palestine",country:="Palestine"]
pop[country=="United States of America",country:="United States"]
pop[country=="Venezuela (Bolivarian Republic of)",country:="Venezuela"]
pop[country=="Viet Nam",country:="Vietnam"]

# DEATH DATA

# Read in WHO death data
who_deaths = get_national_data(source="who")
setDT(who_deaths)
dir.create(paste0("../who_data/",date_fitting),recursive = T)
write.csv(who_deaths,paste0("../who_data/",date_fitting,"/who_data.csv"),row.names = F)

# Get England data
uk_data = get_regional_data("United Kingdom")
setDT(uk_data)
cols = c("region",intersect(names(who_deaths),names(uk_data)))
phe_deaths = uk_data[region=="England",..cols]
setnames(phe_deaths,"region","country")
dir.create(paste0("../phe_data/",date_fitting),recursive = T)
write.csv(phe_deaths,paste0("../phe_data/",date_fitting,"/phe_data.csv"),row.names = F)

# Bind to WHO data
who_deaths = rbind(who_deaths,phe_deaths,fill=T)

# Read in age-stratified death data
source_deaths = "coverage"
deaths_raw = read_death_data(source_deaths)

# Clean age-stratified death data
deaths = clean_death_data(source_deaths,deaths_raw,who_deaths)

# Plot to check
ggplot(deaths[,.(date,deaths_both=c(0,diff(cum_deaths_both))),by=.(country,age_group)],aes(x=date,y=deaths_both,group=age_group,color=age_group)) + 
    geom_line() + 
    facet_wrap(~country)

# Read in LTC death data
ltc_deaths = fread("../ltccovid_data/ltc_deaths.csv")

# VACCINATION DATA

# Read in ECDC vaccination data
# ecdc_vax = fread("../ecdc_data/ecdc_vaccination_data.csv")
ecdc_vax = read.csv("https://opendata.ecdc.europa.eu/covid19/vaccine_tracker/csv", na.strings = "", fileEncoding = "UTF-8-BOM")
setDT(ecdc_vax)
dir.create(paste0("../ecdc_data/",date_fitting),recursive = T)
write.csv(ecdc_vax,paste0("../ecdc_data/",date_fitting,"/ecdc_vaccination_data.csv"),row.names=F)

# Clean ECDC vaccination data
# out = clean_ecdc_vaccination_data(ecdc_vax,country_iso_codes)
# vax = out$vax
# num_type = out$num_type
vax = clean_ecdc_vaccination_data(ecdc_vax,country_iso_codes)

# Read in processed Public Health England (PHE) vaccination data
vaccPHE = readRDS("../phe_data/vax-covidm20211002121137.rds")

# Process to same format as cleaned ECDC data
vaxENG = process_phe_vaccination_data(vaccPHE,agegroups_model)

# Make population data table for England
popENG = CJ(country="England",age=0:100)
agegroups_pop = popUK[,unique(age)]
min_ages_pop = get_min_age(agegroups_pop)
popENG[,age_group:=cut(age,c(min_ages_pop,Inf),labels=agegroups_pop,right=F)]
popENG = merge(popENG,popUK,by.x=c("country","age_group"),by.y=c("name","age"))
popENG = popENG[,.(country,age,female=1000*f/5,male=1000*m/5)]
popENG[,population:=female+male]

pop = rbind(pop,popENG,fill=T)

# IFR

# Read in ensemble IFR estimate from O'Driscoll et al Nature 2020
ifr = fread(datapath("IFR_by_age_ODriscoll.csv"))

# VARIANT DATA

# Read in ECDC variant data
# ecdc_vrnt_data = fread("../ecdc_data/ecdc_variant_data.csv")
ecdc_vrnt_data = read.csv("https://opendata.ecdc.europa.eu/covid19/virusvariant/csv")
setDT(ecdc_vrnt_data)
dir.create(paste0("../ecdc_data/",date_fitting),recursive = T)
write.csv(ecdc_vrnt_data,paste0("../ecdc_data/",date_fitting,"/ecdc_variant_data.csv"),row.names=F)

# Convert ISO weeks to dates
ecdc_vrnt_data[,date:=as.IDate(ISOweek2date(paste0(sub("-","-W",year_week),"-7")))]

# Plot data for different sources
plot_variant_data(ecdc_vrnt_data,"GISAID")
plot_variant_data(ecdc_vrnt_data,"TESSy")

# Process ECDC variant data
vrnt_data = process_variant_data(ecdc_vrnt_data)

# Read in COG England data
# cog_vrnt_data = fread("../cog_data/lineages_by_ltla_and_week.tsv")
cog_vrnt_data = read.delim("https://covid-surveillance-data.cog.sanger.ac.uk/download/lineages_by_ltla_and_week.tsv", na.strings = "")
setDT(cog_vrnt_data)
dir.create(paste0("../cog_data/",date_fitting),recursive = T)
write.table(cog_vrnt_data,paste0("../cog_data/",date_fitting,"/lineages_by_ltla_and_week.tsv"),sep="\t",row.names = F)

# Convert week end date from Date to IDate
cog_vrnt_data[,WeekEndDate:=as.IDate(WeekEndDate)]

# Process COG data
vrnt_dataENG = process_cog_variant_data(cog_vrnt_data)

# Bind to ECDC data
vrnt_data = rbind(vrnt_data,vrnt_dataENG,fill=T)

# Plot non-Alpha-Delta/Alpha/Delta proportions over time
ggplot(vrnt_data,aes(x=date,y=prop_vrnt,group=vrnt,color=vrnt)) +
    geom_line() +
    facet_wrap(~country)

# Estimate variant proportions by multinomial regression
res = estimate_variant_proportions(vrnt_data,deaths[,max(date)])

# Plot estimates for both models
plot_variant_proportions(res$vrnt_prop1,vrnt_data)
plot_variant_proportions(res$vrnt_prop2,vrnt_data)

# Pick model 2 (model with a natural cubic spline function of sample date) as it 
# has lower AIC
vrnt_prop = res$vrnt_prop2

# Construct data tables of deaths, vaccinations and IFR for deconvolution age
# groups
cols = c("cum_deaths_male","cum_deaths_female","cum_deaths_both")
dt = construct_data_table(agegroups,deaths,pop,cols,ltc_deaths,vax,Ab_delay1,Ab_delay2,vaxENG,ifr,vrnt_prop,ve_params)

# Plot data to check
# by age
ggplot(dt,aes(x=date,y=deaths_i_both,group=age_group,color=age_group)) +
    geom_line() +
    facet_wrap(~country)

# overall against WHO data
ggplot() +
    geom_line(aes(x = date, y = deaths_new),who_deaths[country %in% dt[,unique(country)]]) +
    geom_line(aes(x=date,y=deaths_i_both),dt[,.(deaths_i_both=sum(deaths_i_both)),by=.(country,date)],size=0.2,color="red") +
    facet_wrap(~country)

# # overall INED against COVerAGE data
# ggplot() +
#     geom_line(aes(x=date,y=deaths_i_both),dt_ined[,.(deaths_i_both=sum(deaths_i_both)),by=.(country,date)]) +
#     geom_line(aes(x=date,y=deaths_i_both),dt[,.(deaths_i_both=sum(deaths_i_both)),by=.(country,date)],size=0.2,color="red") +
#     facet_wrap(~country)

# Drop data for Iceland as there are too few deaths for deconvolution and for 
# Ireland, Latvia and Romania as it is very incomplete
dt = dt[!(country %in% c("Iceland","Ireland","Latvia","Romania"))]

# Plot deaths
plot_deaths(dt)
ggsave(paste0(dir_out,"deaths_by_age_",source_deaths,".png"),width = 6,height = 4.8)

# Plot vaccine coverage
# country-level
start_date = as.Date("2020-12-01")
p1 = plot_ovrl_vax_cov(dt[date>=start_date])
# ggsave(paste0(dir_out,"vax_cov.png"),width = 5,height = 4)

# by country and age
ps1 = plot_vax_cov(dt[date>=start_date])
# ggsave(paste0(dir_out,"vax_cov_by_age.png"),width = 10,height = 8)

# Calculate and plot IFR with different assumptions:
# with scaling for Alpha and Delta severity
dt[,ifr_t:=pmax(((1-(1-ei)*cum_prop_v/(1-ei*cum_prop_v))+cum_prop_v/(1-ei*cum_prop_v)*(1-ei)*(1-ed)*(1-em))*(prop_vrnt+prop_vrnt2*1.5+prop_vrnt3*1.5*1.8)*ifr,0)]
plot_ifr(dt[date>=start_date])
ggsave(paste0(dir_out,"ifr_over_time_sev_scld.png"),width = 10,height = 8)
plot_avg_ifr(dt[date>=start_date])
ggsave(paste0(dir_out,"avg_ifr_over_time_sev_scld.png"),width = 5,height = 4)

# without scaling for Alpha and Delta severity
dt[,ifr_t:=pmax(((1-(1-ei)*cum_prop_v/(1-ei*cum_prop_v))+cum_prop_v/(1-ei*cum_prop_v)*(1-ei)*(1-ed)*(1-em))*ifr,0)]
ps2 = plot_ifr(dt[date>=start_date])
# ggsave(paste0(dir_out,"ifr_over_time_not_sev_scld.png"),width = 10,height = 8)
p2 = plot_avg_ifr(dt[date>=start_date])
# ggsave(paste0(dir_out,"avg_ifr_over_time_not_sev_scld.png"),width = 5,height = 4)
p = plot_grid(p1+theme(legend.position="none"),
              p2+theme(legend.position="none"),
              labels = c("A","B"))
l = get_legend(p1 + theme(legend.box.margin = margin(0,0,0,12)))
ggsave(paste0(dir_out,"vax_cov_and_IFR.png"),plot_grid(p,l,rel_widths = c(3,.4)),width = 10,height = 4)

ps = plot_grid(ps1+theme(legend.position="none"),
               ps2+theme(legend.position="none"),
               labels = c("A","B"))
ls = get_legend(ps1 + theme(legend.box.margin = margin(0,0,0,12)))
ggsave(paste0(dir_out,"vax_cov_and_IFR_by_age.png"),plot_grid(ps,ls,rel_widths = c(3,.4)),width = 10,height = 4)

# Plot variant proportions
plot_variant_proportions(vrnt_prop[country %in% dt[,unique(country)]],vrnt_data[country %in% dt[,unique(country)]])
ggsave(paste0(dir_out,"vrnt_prop_over_time.png"),width = 7.5,height = 6)

# 
# BACKCALCULATION
# 


# Convolve distributions to get incubation period and exposure-to-death delay distribution
dIncub = disc_conv(cm_delay_gamma(2.5, 2.5, t_max = 30, t_step = 1)$p,cm_delay_gamma(2.5, 4.0, t_max = 30, t_step = 1)$p)
dDeath = disc_conv(cm_delay_gamma(2.5, 2.5, t_max = 60, t_step = 1)$p,cm_delay_gamma(15, 2.2, t_max = 60, t_step = 1)$p)
# Normalise to ensure distributions sum to 1
dIncub = dIncub/sum(dIncub)
dDeath = dDeath/sum(dDeath)

# Frailty index for relative frailty of LTC residents compared to general population
frlty_idx = 3.8

## Backcalculate IFR-scaled infections
method = "ride"
out = run_backcalculation(dt,dDeath,dIncub,method = method)
backcalc_dt = out$backcalc_dt
backcalc_samps = out$backcalc_samps
backcalc_dt_non_ltc = out$backcalc_dt_non_ltc
backcalc_samps_non_ltc = out$backcalc_samps_non_ltc
backcalc_dt_ltc = out$backcalc_dt_ltc
backcalc_samps_ltc = out$backcalc_samps_ltc

rm(out)

# Calculate 95% credible intervals
out = calc_exposures_and_infections_CI(backcalc_dt_non_ltc,backcalc_samps_non_ltc,dIncub)
backcalc_dt_non_ltc = out$dt
exposures_samps_non_ltc = out$exposures_samps
infections_samps_non_ltc = out$infections_samps
out = calc_exposures_and_infections_CI(backcalc_dt_ltc,backcalc_samps_ltc,dIncub)
backcalc_dt_ltc = out$dt
exposures_samps_ltc = out$exposures_samps
infections_samps_ltc = out$infections_samps

rm(out)

# Calculate overall exposures 
exposures_samps = vector("list",length(exposures_samps_non_ltc))
infections_samps = vector("list",length(exposures_samps_non_ltc))

countries = backcalc_dt[,unique(country)]
agegroups_ltc = backcalc_dt_ltc[,unique(age_group)]
for (i in 1:length(countries)){
    print(i)
    for (j in 1:length(agegroups)){
        k = (i-1)*length(agegroups)+j
        cntry = countries[i]
        age_grp = agegroups[j]
        exposures_samps[[k]] = exposures_samps_non_ltc[[k]]
        infections_samps[[k]] = infections_samps_non_ltc[[k]]
        if (age_grp %in% agegroups_ltc){
            idx = (i-1)*length(agegroups_ltc)+which(agegroups_ltc==age_grp)
            exposures_samps[[k]] = exposures_samps[[k]] + exposures_samps_ltc[[idx]]
            infections_samps[[k]] = infections_samps[[k]] + infections_samps_ltc[[idx]]
        }
        backcalc_dt[country==cntry & age_group==age_grp,`:=`(exposures_dead_q95l=apply(backcalc_samps[[k]],2,function(x) quantile(x,probs = 0.025)),exposures_dead_q95u=apply(backcalc_samps[[k]],2,function(x) quantile(x,probs = 0.975)))]
        backcalc_dt[country==cntry & age_group==age_grp,`:=`(exposures_q95l=apply(exposures_samps[[k]],2,function(x) quantile(x,probs = 0.025)),exposures_q95u=apply(exposures_samps[[k]],2,function(x) quantile(x,probs = 0.975)))]
        backcalc_dt[country==cntry & age_group==age_grp,`:=`(infections_q95l=apply(infections_samps[[k]],2,function(x) quantile(x,probs = 0.025)),infections_q95u=apply(infections_samps[[k]],2,function(x) quantile(x,probs = 0.975)))]
    }
}

save.image(paste0(dir_out,"backcalculation_output.RData"))


#
# PLOTTING
#

# Read in ECDC age-stratified case data for comparison with inferred infection time series
ecdc_cases_by_age = read.csv("https://opendata.ecdc.europa.eu/covid19/agecasesnational/csv", na.strings = "", fileEncoding = "UTF-8-BOM")
setDT(ecdc_cases_by_age)
write.csv(ecdc_cases_by_age,paste0("../ecdc_data/",date_fitting,"/ecdc_cases_by_age.csv"),row.names=F)

# Read in and process seroprevalence data from SeroTracker for comparison with estimated cumulative proportion infected
serotracker_data = fread("../serotracker_data/SeroTracker Serosurveys Reporting Prevalence.csv")
sero_data = process_seroprevalence_data(serotracker_data)

# Plot backcalculation output
plot_output(paste0(dir_out,"backcalculation_output.RData"),ecdc_cases_by_age,sero_data)

