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
source("./initial_conditions_functions.R")

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
# - combine data for England, Wales and Scotland
# - separate into LTC and non-LTC deaths - DONE, NEED TO CORRECT
#   - don't apply split to France data as it already excludes care home deaths
#   - use mean proportion of deaths that are among LTC residents for countries
#   missing LTC death data - DONE
# - use sampling approach to get integer numbers of deaths rather than 
# multiplying by proportions and ending up with decimals (so that data is in 
# correct format for RIDE)?
# - make proportion of vaccine type time dependent
# - include time-dependent variant proportions in vaccine efficacy
# - use country-specific IFRs?
# - correct dates for vaccination data to be end of ISO week - no, start of ISO week is fine
# - sort out Finland data - interpolate with WHO data?


# 
# SETUP
# 

date_fitting = today()

# Create directory to save output into
dir_out = paste0("./output/",date_fitting,"/")
dir.create(dir_out,recursive = T)

# Set plot theme
theme_set(cowplot::theme_cowplot(font_size = 10) + theme(strip.background = element_blank()))

# Register parallel backend
registerDoParallel(cores = detectCores()-1)
# registerDoParallel(cores = 4)

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
regions = cm_populations[location_type==4 & name %in% names(cm_matrices),unique(name)]

# Set delay distributions
dE  = cm_delay_gamma(2.5, 2.5, t_max = 15, t_step = 1)$p
dIp = cm_delay_gamma(2.5, 4.0, t_max = 15, t_step = 1)$p
dIs = cm_delay_gamma(2.5, 4.0, t_max = 15, t_step = 1)$p
dIa = cm_delay_gamma(5.0, 4.0, t_max = 15, t_step = 1)$p
dHosp = cm_delay_gamma(6.0 + 2.5, 0.71, 60, 1)$p
dDeath = cm_delay_lnorm(15, 0.9, 60, 1)$p

# Build parameters for different regions ###
params = cm_parameters_SEI3R(regions, deterministic = T, 
                             date_start = "2020-01-01", 
                             date_end = date_fitting,
                             dE  = dE,
                             dIp = dIp,
                             dIs = dIs,
                             dIa = dIa)
params = cm_split_matrices_ex_in(params, 15)

# get list of age groups in covidm
agegroups_model <- params$pop[[1]]$group_names
min_ages_model = get_min_age(agegroups_model)
# max_ages_model = get_max_age(agegroups_model)
# max_ages_model[is.na(max_ages_model)] = Inf

# Set vaccine efficacy parameters
ve_params = list()
# eiX_vYZ = efficacy of dose Z of vaccine Y against infection with strain X
ve_params$ei_va2 = 0.68
ve_params$ei2_va2 = 0.68
ve_params$ei3_va2 = 0.6154
ve_params$ei_vb2 = 0.85
ve_params$ei2_vb2 = 0.85
ve_params$ei3_vb2 = 0.7999
# ed_vYZiX = efficacy of dose Z of vaccine Y against disease given infection with strain X
ve_params$ed_va2i = 0.3125
ve_params$ed_va2i2 = 0.3125
ve_params$ed_va2i3 = 0.2353094
ve_params$ed_vb2i = 0.2667
ve_params$ed_vb2i2 = 0.2667
ve_params$ed_vb2i3 = 0.187906
# eh_vYZiX = efficacy of dose Z of vaccine Y against hospitalisation given disease from strain X
ve_params$eh_va2d = 0.55
ve_params$eh_va2d2 = 0.55
ve_params$eh_va2d3 = 0.5147909
ve_params$eh_vb2d = 0.09
ve_params$eh_vb2d2 = 0.09
ve_params$eh_vb2d3 = 0.2215385
# em_vYZdX = efficacy of dose Z of vaccine Y against death given disease from strain X
ve_params$em_va2d = 0.77
ve_params$em_va2d2 = 0.77
ve_params$em_va2d3 = 0.6766406
ve_params$em_vb2d = 0.55
ve_params$em_vb2d2 = 0.55
ve_params$em_vb2d3 = 0.52
# et_vYZdX = efficacy of dose Z of vaccine Y against onward transmission given infection with strain X
ve_params$et_va2i = 0.5
ve_params$et_va2i2 = 0.5
ve_params$et_va2i3 = 0.4525
ve_params$et_vb2i = 0.6
ve_params$et_vb2i2 = 0.6
ve_params$et_vb2i3 = 0.5646

# Set delay from vaccination to development of Ab
Ab_delay = 14 # days

# 
# DATA 
# 


# POPULATION DATA

# Read in UN population data
pop = qread("../../OneDrive - London School of Hygiene and Tropical Medicine/LSHTM_RF/COVID/comix/data/unwpp_data.qs")
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
write.csv(who_deaths,"../who_data/who_data.csv",row.names = F)

# Read in INED age-stratified death data
ined_deaths = fread("../ined_data/AgeSex/Cum_deaths_by_age_sex.csv")
# pop = fread("../ined_data/AgeSex/covid_pooled_POP_2021-04-01.csv")
#
# names(pop) = tolower(names(pop))
# pop[,age_group:=sub("'","",age_group)]
# deaths = merge(deaths,pop,by=c("region","country","country_code","age_group"))

# Set variables to be included in cleaned data
# cols = c("cum_deaths_male","cum_deaths_female","cum_deaths_unknown","cum_deaths_both")
cols = c("cum_deaths_male","cum_deaths_female","cum_deaths_both")

# Clean INED death data
deaths = clean_ined_death_data(ined_deaths,who_deaths,cols)

# Plot to check
ggplot(deaths[country %in% c("Austria","Denmark","France","Spain"),.(date,deaths_both=c(0,diff(cum_deaths_both))),by=.(country,age_group)],aes(x=date,y=deaths_both,group=age_group,color=age_group)) + 
    geom_line() + 
    facet_wrap(~country)

# Read in COVerAGE death data
coverage_deaths = fread(cmd = "unzip -cq ../coverage_data/Output_10.zip", skip = 3)

# Clean coverage death data
out10 = clean_coverage_death_data(coverage_deaths,who_deaths,"deaths")

# Plot to check
ggplot(out10[country %in% c("Austria","Denmark","Finland","France","Spain"),.(date,deaths_both=c(0,diff(cum_deaths_both))),by=.(country,age_group)],aes(x=date,y=deaths_both,group=age_group,color=age_group)) + 
    geom_line() + 
    facet_wrap(~country)

# Read in LTC death data
ltc_deaths = fread("../ltccovid_data/ltc_deaths.csv")

# Read in ECDC vaccination data
ecdc_vax = fread("../ecdc_data/ecdc_vaccination_data.csv")

# Clean ECDC vaccination data
# out = clean_ecdc_vaccination_data(ecdc_vax,country_iso_codes)
# vax = out$vax
# num_type = out$num_type
vax = clean_ecdc_vaccination_data(ecdc_vax,country_iso_codes)

# Read in ensemble IFR estimate from O'Driscoll et al Nature 2020
ifr = fread(datapath("IFR_by_age_ODriscoll.csv"))

# Read in ECDC variant data
vrnt_data = fread("../ecdc_data/ecdc_variant_data.csv")

vrnt_data[,date:=ISOweek2date(paste0(sub("-","-W",year_week),"-7"))]

# Classify variants into non-Alpha-Delta, Alpha and Delta
vrnt_data[,vrnt:="Other"]
vrnt_data[variant=="B.1.1.7",vrnt:="Alpha"]
vrnt_data[variant=="B.1.617.2",vrnt:="Delta"]

# Plot variant proportions over time
# GISAID data
ggplot(vrnt_data[source=="GISAID"],aes(x=date,y=percent_variant,group=variant,color=variant)) +
    geom_line() +
    facet_wrap(~country)
# TESSy data
ggplot(vrnt_data[source=="TESSy"],aes(x=date,y=percent_variant,group=variant,color=variant)) +
    geom_line() +
    facet_wrap(~country)

# Aggregate data by Alpha/Delta status
vrnt_data = vrnt_data[source=="GISAID",.(number_detections_variant=sum(number_detections_variant)),by=.(country,country_code,year_week,date,new_cases,number_sequenced,percent_cases_sequenced,vrnt)]
vrnt_data[,prop_vrnt:=number_detections_variant/number_sequenced]

# Plot non-Alpha-Delta/Alpha/Delta proportions over time
ggplot(vrnt_data,aes(x=date,y=prop_vrnt,group=vrnt,color=vrnt)) +
    geom_line() +
    facet_wrap(~country)

# # Shift dates back by 4 days to account for reporting delay (as in covidm code)/
# # proportion corresponding to middle of ISO week
# vrnt_data[,date:=date-4]

# # Construct data table of variant proportions since start of pandemic
# vrnt_prop = CJ(country=vrnt_data[,unique(country)],date=dates2,vrnt=vrnt_data[,unique(vrnt)])
# vrnt_prop = merge(vrnt_prop,vrnt_data[,.(country,date,vrnt,prop_vrnt)],by=c("country","date","vrnt"),all.x=T)
#
# # Interpolate missing values
# vrnt_prop[date==min(date),prop_vrnt:=fcase(vrnt=="Alpha",0,
#                                            vrnt=="Delta",0,
#                                            vrnt=="Other",1)]
# vrnt_prop[,prop_vrnt:=approx(date,prop_vrnt,date)$y,by=.(country,vrnt)]
#
# # Plot to check
# ggplot(vrnt_prop,aes(x=date,y=prop_vrnt,group=vrnt,color=vrnt)) +
#     geom_line() +
#     facet_wrap(~country)

# Fit multinomial logistic model to estimate variant proportions over time
# Cast to wide format
vrnt_data_wide = dcast(vrnt_data,country + date ~ vrnt,value.var = "number_detections_variant")
# Exclude rows without any observations
vrnt_data_wide = vrnt_data_wide[!(Alpha==0 & Delta==0 & Other==0)]

countries1 = vrnt_data_wide[,unique(country)]
ncountries = length(countries1)
m1 = vector("list",ncountries)
m2 = vector("list",ncountries)
for (i in 1:ncountries){
    cntry = countries1[i]
    m1[[i]] = multinom(as.matrix(vrnt_data_wide[country==cntry,.(Other,Alpha,Delta)]) ~ date,vrnt_data_wide[country==cntry])
    m2[[i]] = multinom(as.matrix(vrnt_data_wide[country==cntry,.(Other,Alpha,Delta)]) ~ ns(date,df=2),vrnt_data_wide[country==cntry])
}

# Compare AIC and BIC for fitted models
AIC_BIC = data.table(model=c("m1","m2"),
AIC=c(sum(sapply(m1,AIC)),sum(sapply(m2,AIC))),
BIC=c(sum(sapply(m1,BIC)),sum(sapply(m2,BIC))))
print(AIC_BIC)
#    model      AIC      BIC
# 1:    m1 711620.8 712474.5
# 2:    m2 667952.1 669232.7
# Use model m2 as it has lower AIC and BIC

vrnt_prop1 = vector("list",ncountries)
vrnt_prop2 = vector("list",ncountries)
new_dt = data.table(date=seq.Date(vrnt_data_wide[,min(date)],vrnt_data_wide[,max(date)],by=1)) #data.table(date = seq.Date(vrnt_data_wide[,min(date)],dates2[length(dates2)],by=1))
for (i in 1:ncountries){
    cntry = countries1[i]
    # new_dt = data.table(date=seq.Date(vrnt_data_wide[country==cntry,min(date)],vrnt_data_wide[country==cntry,max(date)],by=1))
    vrnt_prop1[[i]] = data.table(new_dt,predict(m1[[i]],newdata = new_dt,type = "probs"))
    vrnt_prop1[[i]][,country:=cntry]
    vrnt_prop2[[i]] = data.table(new_dt,predict(m2[[i]],newdata = new_dt,type = "probs"))
    vrnt_prop2[[i]][,country:=cntry]
}
vrnt_prop1 = rbindlist(vrnt_prop1)
vrnt_prop2 = rbindlist(vrnt_prop2)
vrnt_prop1 = melt(vrnt_prop1,id.vars = c("country","date"),variable.name = "vrnt",value.name = "prop_vrnt")
vrnt_prop2 = melt(vrnt_prop2,id.vars = c("country","date"),variable.name = "vrnt",value.name = "prop_vrnt")

ggplot() +
    geom_line(aes(x=date,y=prop_vrnt,color=vrnt),vrnt_prop1) +
    geom_point(aes(x=date,y=prop_vrnt,color=vrnt),vrnt_data) +
    facet_wrap(~country)
ggplot() +
    geom_line(aes(x=date,y=prop_vrnt,color=vrnt),vrnt_prop2) +
    geom_point(aes(x=date,y=prop_vrnt,color=vrnt),vrnt_data) +
    facet_wrap(~country)

# TODO - work out how to calculate CIs with effects package. Might need to
# reformat input data to get it to work

# Construct data tables of deaths, vaccinations and IFR for deconvolution age
# groups
# dt_ined = construct_data_table(agegroups,deaths,pop,cols,ltc_deaths,vax,num_type,ifr)
# dt = construct_data_table(agegroups,out10,pop,c("cum_deaths_male","cum_deaths_female","cum_deaths_both"),ltc_deaths,vax,num_type,ifr)
dt_ined = construct_data_table(agegroups,deaths,pop,cols,ltc_deaths,vax,Ab_delay,ifr,vrnt_prop2,ve_params)
dt = construct_data_table(agegroups,out10,pop,c("cum_deaths_male","cum_deaths_female","cum_deaths_both"),ltc_deaths,vax,Ab_delay,ifr,vrnt_prop2,ve_params)

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

# overall INED against COVerAGE data
ggplot() +
    geom_line(aes(x=date,y=deaths_i_both),dt_ined[,.(deaths_i_both=sum(deaths_i_both)),by=.(country,date)]) +
    geom_line(aes(x=date,y=deaths_i_both),dt[,.(deaths_i_both=sum(deaths_i_both)),by=.(country,date)],size=0.2,color="red") +
    facet_wrap(~country)

# Drop data for Iceland as there are too few deaths for deconvolution and for 
# Ireland, Latvia and Romania as it is very incomplete
dt = dt[!(country %in% c("Iceland","Ireland","Latvia","Romania"))]


# 
# BACKCALCULATION
# 


# # Read in incubation period and onset-to-death delay data
# ip = readRDS("incubation_period.rds")
# odd = readRDS("onset_to_death_delay.rds")
# 
# # Convolve distributions to get infection-to-death delay distribution
# ip_pdf = function(x) dlnorm(x,ip$mean,ip$sd)
# odd_pdf = function(y) dlnorm(y,odd$mean,odd$sd)
# idd_pdf = function(z) integrate(function(x,z) odd_pdf(z-x)*ip_pdf(x),-Inf,Inf,z)$value
# ip_pdf = Vectorize(ip_pdf)
# idd_pdf = Vectorize(idd_pdf)
# 
# # Discretised incubation period and infection-to-death delay distributions
# ip_max = ip$max
# ip_pmf = c(0,sapply(0:(ip_max-1),function(x){integrate(function(z) ip_pdf(z),x,x+1)$value}))/integrate(function(z) ip_pdf(z),0,ip_max)$value
# idd_max = ip$max + odd$max
# idd_pmf = c(0,sapply(0.5:(idd_max-0.5),function(x){integrate(function(z) idd_pdf(z),x,x+1)$value}))/integrate(function(z) idd_pdf(z),0.5,idd_max+0.5)$value
# # Renormalise to correct rounding error
# idd_pmf = idd_pmf/sum(idd_pmf)
# 
# # mean_ip = exp(ip$mean+ip$sd^2/2)
# # mean_odd = exp(odd$mean+odd$sd^2/2)
# # mean_idd = round(mean_ip + mean_odd)

# Convolve distributions to get incubation period and exposure-to-death delay distribution
dIncub = disc_conv(cm_delay_gamma(2.5, 2.5, t_max = 30, t_step = 1)$p,cm_delay_gamma(2.5, 4.0, t_max = 30, t_step = 1)$p)
dDeath = disc_conv(cm_delay_gamma(2.5, 2.5, t_max = 60, t_step = 1)$p,cm_delay_gamma(15, 2.2, t_max = 60, t_step = 1)$p)
# Normalise to ensure distributions sum to 1
dIncub = dIncub/sum(dIncub)
dDeath = dDeath/sum(dDeath)

# Frailty index for relative frailty of LTC residents compared to general population
frlty_idx = 3.8

## Backcalculate IFR-scaled infections
# Extract non-LTC and LTC data
dt_non_ltc = dt[,.(country,age_group,date,deaths_i_both=deaths_i_both_non_ltc,ifr_t)]
dt_ltc = dt[get_min_age(age_group)>=60,.(country,age_group,date,deaths_i_both=deaths_i_both_ltc,ifr_t=frlty_idx*ifr_t)]

# Deconvolve deaths to get IFR-scaled exposures
tstart = Sys.time()
# dt_backproj = backcalc(dt,idd_pmf,ip_pmf,method = "backproj")
out = backcalc(dt_non_ltc,dDeath,dE,method = "backproj")
dt_non_ltc_backproj = out[[1]]
samps_non_ltc_backproj = out[[2]]
out = backcalc(dt_ltc,dDeath,dE,method = "backproj")
dt_ltc_backproj = out[[1]]
samps_ltc_backproj = out[[2]]
tend = Sys.time()
print(tend-tstart)

tstart = Sys.time()
# dt_ride = backcalc(dt,idd_pmf,ip_pmf,method = "ride") # takes about 30 mins
out = backcalc(dt_non_ltc,dDeath,dIncub,method = "ride")
dt_non_ltc_ride = out[[1]]
samps_non_ltc_ride = out[[2]]
out = backcalc(dt_ltc,dDeath,dIncub,method = "ride")
dt_ltc_ride = out[[1]]
samps_ltc_ride = out[[2]]
tend = Sys.time()
print(tend-tstart) # takes about 40 mins

# Merge non-LTC and LTC estimates
dt_ride = merge(dt_non_ltc_ride,dt_ltc_ride,by=c("country","age_group","date"),all.x=T,suffixes=c("_non_ltc","_ltc"))
setnafill(dt_ride,fill=0,cols=c("deaths_i_both_ltc","exposures_dead_ltc","exposures_ltc","infections_ltc"))
dt_ride[,`:=`(deaths_i_both=deaths_i_both_non_ltc+deaths_i_both_ltc,
              exposures_dead=exposures_dead_non_ltc+exposures_dead_ltc,
              exposures=exposures_non_ltc+exposures_ltc,
              infections=infections_non_ltc+infections_ltc)]
dt_ride = merge(dt_ride,unique(dt[,.(country,age_group,population)]),by=c("country","age_group"),all.x=T)

# Calculate cumulative infections
dt_ride = calc_cum_exposures_and_infections(dt_ride)

save.image(paste0(dir_out,"backcalculation_output.RData"))

#
# PLOTTING
#

load("./output/2021-09-15/backcalculation_output.RData")
# p_backproj = plot_infections(dt_backproj)
# ggsave("./output/infections_by_age_backproj.pdf",plot=p_backproj,width = 10,height = 8)
# p_ride = plot_infections(dt_ride)
# ggsave("./output/infections_by_age_ride.pdf",plot=p_ride,width = 10,height = 8)
# p_backproj_DNK_NOR = plot_infections(dt_backproj[country_code %in% c("DNK","NOR")])
# ggsave("./output/infections_by_age_backproj_DNK_NOR.pdf",plot=p_backproj_DNK_NOR,width = 6,height = 3)
# p_ride_DNK_NOR = plot_infections(dt_ride[country_code %in% c("DNK","NOR")])
# ggsave("./output/infections_by_age_ride_DNK_NOR.pdf",plot=p_ride_DNK_NOR,width = 6,height = 3)
p_non_ltc_backproj = plot_infections(dt_non_ltc_backproj)
ggsave(paste0(dir_out,"infections_by_age_non_LTC_backproj_COVerAGE.pdf"),plot=p_non_ltc_backproj,width = 10,height = 8)
p_ltc_backproj = plot_infections(dt_ltc_backproj)
ggsave(paste0(dir_out,"infections_by_age_LTC_backproj_COVerAGE.pdf"),plot=p_ltc_backproj,width = 10,height = 8)
p_non_ltc_ride = plot_infections(dt_non_ltc_ride) #[country %in% c("Denmark","Norway")])
p_non_ltc_ride = p_non_ltc_ride + ylim(c(0,5e4)) # restrict y-axis due to issues with Finland estimates
ggsave(paste0(dir_out,"infections_by_age_non_LTC_ride_COVerAGE.pdf"),plot=p_non_ltc_ride,width = 10,height = 8)
p_ltc_ride = plot_infections(dt_ltc_ride)
ggsave(paste0(dir_out,"infections_by_age_LTC_ride_COVerAGE.pdf"),plot=p_ltc_ride,width = 10,height = 8)

# Compare with ECDC age-stratified case data
ecdc_cases_by_age = fread("../ecdc_data/ecdc_case_data_by_age.csv")
# Should be end of ISO week - TODO correct in vaccination data!
ecdc_cases_by_age[,date:=ISOweek2date(paste0(sub("-","-W",year_week),"-7"))]
ecdc_cases_by_age[age_group=="<15yr",age_group:="0-14yr"]
ggplot(ecdc_cases_by_age[country %in% dt[,unique(country)]],aes(x=date,y=new_cases,group=age_group,color=age_group)) +
    geom_line() +
    facet_wrap(~country)
ggsave(paste0(dir_out,"ecdc_cases_by_age.pdf"),width = 10,height = 8)
ggplot(ecdc_cases_by_age[country_code %in% c("DK","NO")],aes(x=date,y=new_cases,group=age_group,color=age_group)) +
    geom_line() +
    facet_wrap(~country)
ggsave(paste0(dir_out,"ecdc_cases_by_age_DNK_NOR.pdf"),width = 6,height = 3)

# Plot estimated infections vs reported cases
# overall
ggplot() +
    geom_line(aes(x=date,y=infections),dt_ride[,.(infections=sum(infections)),by=.(country,date)]) +
    geom_line(aes(x=date,y=new_cases/7),ecdc_cases_by_age[country %in% dt_ride[,unique(country)],.(new_cases=sum(new_cases)),by=.(country,date)],linetype="dashed") +
    facet_wrap(~country) + ylim(c(0,6e4))
ggsave(paste0(dir_out,"infections_vs_obs_cases_ride_COVerAGE.pdf"),width = 10,height = 8)

# by age
agegroups_comp = c("0-49","50-79","80+")
min_ages_comp = get_min_age(agegroups_comp)
dt_ride[,age_group_comp:=cut(get_min_age(age_group),c(min_ages_comp,Inf),labels=agegroups_comp,right=F)]
ecdc_cases_by_age[,age_group_comp:=cut(get_min_age(age_group),c(min_ages_comp,Inf),labels=agegroups_comp,right=F)]
ggplot() +
    geom_line(aes(x=date,y=infections,group=age_group_comp,color=age_group_comp),dt_ride[,.(infections=sum(infections)),by=.(country,date,age_group_comp)]) +
    geom_line(aes(x=date,y=new_cases/7,group=age_group_comp,color=age_group_comp),ecdc_cases_by_age[country %in% dt_ride[,unique(country)],.(new_cases=sum(new_cases)),by=.(country,date,age_group_comp)],linetype="dashed") +
    facet_wrap(~country) +  ylim(c(0,5e4))
ggsave(paste0(dir_out,"infections_vs_obs_cases_by_age_ride_COVerAGE.pdf"),width = 10,height = 8)

# Denmark and Norway
ggplot() +
    geom_line(aes(x=date,y=infections,group=age_group_comp,color=age_group_comp),dt_ride[country %in% c("Denmark","Norway"),.(infections=sum(infections)),by=.(country,date,age_group_comp)]) +
    geom_line(aes(x=date,y=new_cases/7,group=age_group_comp,color=age_group_comp),ecdc_cases_by_age[country %in% c("Denmark","Norway"),.(new_cases=sum(new_cases)),by=.(country,date,age_group_comp)],linetype="dashed") +
    facet_wrap(~country)
ggsave(paste0(dir_out,"infections_vs_obs_cases_by_age_ride_DNK_NOR_COVerAGE.pdf"),width = 6,height = 3)

# Compare with SeroTracker seroprevalence data
sero_data = fread("../serotracker_data/SeroTracker Serosurveys Reporting Prevalence.csv")
# sero_data[,country_code:=dt[match(sero_data[,Country],country),country_code]]
# Filter to only include national and regional surveys with household and community and residual sera samples
sero_data = sero_data[`Grade of Estimate Scope` %in% c("National","Regional") &
                          `Sample Frame (groups of interest)` %in% c("Household and community samples","Residual sera"),
                      .(#country_code,
                        scope=`Grade of Estimate Scope`,
                        country=Country,
                        region=`Specific Geography`,
                        Start.date=dmy(`Sampling Start Date`),
                        End.date=dmy(`Sampling End Date`),
                        sample_frame=`Sample Frame (groups of interest)`,
                        Central.estimate=`Serum positive prevalence`,
                        Lower.bound=`Serum pos prevalence, 95pct CI Lower`,
                        Upper.bound=`Serum pos prevalence, 95pct CI Upper`,
                        Test=`Test Type`,
                        Institution=`Lead Institution`,
                        Data.source=URL,
                        N.tests=`Denominator Value`
                      )]
# Remove social housing serosurvey data for Denmark and early regional estimates for Italy
sero_data = sero_data[!((country=="Denmark" & Central.estimate=="17.30%")|(country=="Italy" & Start.date<=as.Date("2020-05-05")))]
sero_data[,date := Start.date + (End.date - Start.date)/2]

pct = function(x) as.numeric(str_replace_all(x, "%", "")) / 100

# Plot cumulative proportion infected vs seroprevalence
ggplot() +
    geom_line(aes(x=date,y=cum_prop_inf,group=age_group,color=age_group),dt_ride) +
    geom_point(aes(x=date,y=pct(Central.estimate),shape=scope),sero_data[country %in% dt_ride[,unique(country)]]) +
    facet_wrap(~country) + ylim(c(0,1))
ggsave(paste0(dir_out,"cum_prop_infected_vs_seroprev_COVerAGE.pdf"),width = 10,height = 8)

p_deaths = plot_deaths(dt)
ggsave(paste0(dir_out,"deaths_by_age_COVerAGE.pdf"),plot = p_deaths,width = 10,height = 8)

ggplot(dt,aes(x=date,y=cum_prop_va+cum_prop_vb,group=age_group,color=age_group)) +
    geom_line() +
    facet_wrap(~country)
ggsave(paste0(dir_out,"vax_cov_by_age_COVerAGE.pdf"),width = 10,height = 8)

# save.image("../backcalculation_output.RData")


#
# INCIDENCE CALCULATION
#


# DISAGGREGATE INFECTIONS INTO COVIDM AGE GROUPS
# # # Create variable with deaths shifted back by mean of infection-to-death distribution
# # dt_model[,death_i_both_shft:=shift(death_i_both,18,type="lead"),by=.(country_code,age_group)]
# # Take rolling mean of deaths in each age group with large averaging window
# dt_model[,death_i_both_ma:=frollmean(death_i_both,90)]
# # # Calculate proportions of deaths in each age group with rolling mean deaths
# # dt_model[,prop_deaths_ma:=death_i_both_ma/sum(death_i_both_ma),by=.(country_code,date)]
#
# # Calculate proportions of infections in each age group from rolling mean deaths
# # N.B. Doesn't account for delay between infection and death!
# # TODO - revisit this!!!
# dt_model[,inf_ma:=death_i_both_ma/ifr_t]
# dt_model[,age_group2:=cut(get_min_age(age_group),c(min_ages,Inf),labels=agegroups,right=F)]
# dt_model[,prop_inf_ma:=inf_ma/sum(inf_ma),by=.(country_code,age_group2,date)]
#
# # Plot proportion of deaths by age over time
# ggplot(dt_model,aes(x=date,y=death_i_both_ma,group=age_group,color=age_group)) +
#     geom_line() +
#     facet_wrap(~country)
# # Plot proportion of infections by age over time
# ggplot(dt_model[date>as.Date("2020-04-17") & age_group2=="0-39"],aes(x=date,y=prop_inf_ma,group=age_group,color=age_group)) +
#     geom_line() +
#     facet_wrap(~country)

# Make data table for incidence
inc_dt = CJ(country=dt_ride[,unique(country)],age=0:100,date=dt_ride[,unique(date)])
inc_dt = merge(inc_dt,pop[,.(country,age,population)],by=c("country","age"),all.x=T)
inc_dt[,age_group_model:=cut(age,c(min_ages_model,Inf),labels=agegroups_model,right=F)]
inc_dt = inc_dt[,.(population=sum(population)),by=.(country,age_group_model,date)]

# Add deconvolution age groups to data table
inc_dt[,age_group:=cut(get_min_age(age_group_model),c(min_ages,Inf),labels=agegroups,right=F)]

# Merge with deconvolution output
# N.B. duplicates exposures within same age_group
inc_dt = merge(inc_dt,dt_ride[,.(country,age_group,date,exposures)],by=c("country","age_group","date"))

# Split exposures by population fraction
inc_dt[,exposures:=exposures*population/sum(population),by=.(country,age_group,date)]

# Plot
ggplot(inc_dt,aes(x=date,y=exposures,group=age_group_model,color=age_group_model)) +
    geom_line() +
    facet_wrap(~country)

# # DISAGGREGATE INFECTIONS BY VARIANT
# 
# # Merge into incidence data table
# vrnt_prop_wide = dcast(vrnt_prop2,country+date ~ vrnt,value.var = "prop_vrnt")
# cols6 = c("Other","Alpha","Delta") #paste0("prop_vrnt",c("","2","3"))
# setnames(vrnt_prop_wide,c("Other","Alpha","Delta"),cols6)
# inc_dt = merge(inc_dt,vrnt_prop_wide,by=c("country","date"),all.x=T)
# inc_dt[,(cols6):=nafill(.SD,"nocb"),.SDcols=cols6,by=.(country)]
# inc_dt = melt(inc_dt,setdiff(names(inc_dt),cols6),variable.name = "vrnt",value.name = "prop_vrnt")
# 
# # Calculate numbers of exposures for each variant
# inc_dt[,`:=`(nS_E=prop_vrnt*exposures)]
# # inc_dt[,nSV_EL:=prop_vrnt*exposures]
# 
# # # Convolve exposures with latent period and infectious period distribution to
# # # get recoveries
# # inc_dt[,nEL_I:=disc_conv(nSV_EL,dE),by=.()]
# # inc_dt[,nI_R:=disc_conv(nEL_I,dIa),by=.()]


# DISAGGREGATE INFECTIONS INTO SX AND ASX INFECTIONS
# # Make data table of unique combinations of country and age
# countries = deaths[,unique(country)]
# country_codes = deaths[,unique(country_code)]
# base_dt = CJ(country_code = country_codes,age = 0:100)
#
# # Add population data
# base_dt = merge(base_dt,pop[,.(iso3,age,population)],by.x=c("country_code","age"),by.y=c("iso3","age"),all.x=T)

# Read in age-dependent symptomatic fraction
covid_scenario = qread(datapath("2-linelist_both_fit_fIa0.5-rbzvih.qs"));
cols4 = names(covid_scenario)[grep("y_",names(covid_scenario))]
# covy = colMeans(covid_scenario[,..cols4])
#
# min_ages_y = as.integer(sub("y_","",names(covy)))
# covy = data.table(age_group=c(paste0(min_ages_y[1:(length(min_ages_y)-1)],"-",min_ages_y[2:length(min_ages_y)]-1),paste0(min_ages_y[length(min_ages_y)],"+")),y=covy)
# covy_dt = copy(base_dt)
# covy_dt[,age_group:=cut(age,c(min_ages_y,Inf),labels=covy[,age_group],right=F)]
# covy_dt = merge(covy_dt,covy,by="age_group",all.x=T)
#
# # Change age groups
# covy_dt[,age_group:=cut(age,c(min_ages,Inf),labels=agegroups,right=F)]
# # Calculate population-weighted average symptomatic fraction
# covy_dt = covy_dt[,.(y=sum(y*population)/sum(population)),by=.(country,age_group)]
#
# # Plot
# ggplot(covy_dt,aes(x=age_group,y=y,group=country,color=country)) + geom_line()
#
# # Add symptomatic proportion
# inc_dt = merge(inc_dt,covy_dt,by=c("country","age_group"),all.x=T)

covy = unname(rep(colMeans(covid_scenario[,..cols4]), each = 2))
covy_dt = data.table(age_group_model = agegroups_model,y = covy)

# # Add symptomatic proportion to incidence data table
# inc_dt = merge(inc_dt,covy_dt,by="age_group_model")
# 
# # Convolve exposures with latent period distribution to get infections
# setorder(inc_dt,country,age_group_model,vrnt,date)
# inc_dt[,nE_I:=disc_conv(nS_E,dE),by=.(country,age_group_model,vrnt)]
# 
# # Calculate symptomatic and asymptomatic infections
# inc_dt[,`:=`(nE_Ip=y*nE_I,nE_Ia=(1-y)*nE_I)]
# 
# # CONVOLVE INFECTIONS TO GET NUMBERS ENTERING DIFFERENT INFECTION STATES
# # Convolve presymptomatic infections to get symptomatic infections
# inc_dt[,nIp_Is:=disc_conv(nE_Ip,dIp),by=.(country,age_group_model,vrnt)]
# 
# # Convolve infections with infectious period distribution to get recoveries
# inc_dt[,nIs_R:=disc_conv(nIp_Is,dIs),by=.(country,age_group_model,vrnt)]
# inc_dt[,nIa_R:=disc_conv(nE_Ia,dIa),by=.(country,age_group_model,vrnt)]
# 
# ggplot(melt(inc_dt[country=="Austria" & vrnt=="Other"],measure.vars = c("nS_E","nE_Ip","nE_Ia","nIp_Is","nIs_R","nIa_R")),aes(x=date,y=value,group=age_group_model,color=age_group_model)) +
#     geom_line() +
#     facet_wrap(~variable)

# ADD VACCINATION NUMBERS
vax_dt = construct_vax_data_table(vax,inc_dt[,unique(date)]-Ab_delay,agegroups_model)
vax_dt[,date_v:=date]
vax_dt[,date:=date+Ab_delay]
vax_dt_wide = dcast(vax_dt,country + date + age_group + population ~ type,value.var = c("first","second","prop","cum_prop"))

setnames(vax_dt_wide,c("age_group","first_va","first_vb","second_va","second_vb"),c("age_group_model","nS_Va1","nS_Vb1","nVa1_Va2","nVb1_Vb2"))

# inc_dt = merge(inc_dt,vax_dt_wide[,!c("population","prop_va","prop_vb","cum_prop_va","cum_prop_vb")],by=c("country","age_group_model","date"))


#
# INITIAL CONDITIONS CALCULATION
#


# # Calculate prevalences over time
# prev_dt = calc_init_condns(inc_dt,vax_dt_wide)
# 
# ggplot(melt(prev_dt[country=="Austria"],measure.vars = c("E","Ia","Ip","Is")),aes(x=date,y=value,group=age_group_model,color=age_group_model)) +
#     geom_line() +
#     facet_wrap(~variable)

prev_dt = calc_init_condns(inc_dt,vax_dt_wide,agegroups_model,covy,vrnt_prop2,ve_params)

ggplot(prev_dt[country!="Finland"],aes(x=date,y=nS_E+nV_E+nV_L+nV_R,color=age_group_model)) + geom_line() +
    facet_wrap(~country)


#
# REMAINING BURDEN CALCULATION
#


# # Extract current numbers of susceptibles and incidence of exposures and vaccinations
# # TODO - average incidence of exposures and vaccinations over previous few weeks
# # instead of using just last week?
# max_dates = prev_dt[,.(date=max(date)),by=.(country)]
# curr_prev_dt = merge(max_dates,prev_dt[,.(country,age_group_model,date,population,S)],by=c("country","date"),all.x=T)
# curr_inc_dt = merge(max_dates,inc_dt[,.(country,age_group_model,date,vrnt,nS_E)],by=c("country","date"),all.x=T)
# curr_vax_dt = merge(max_dates,vax_dt_wide[,.(country,age_group_model,date,nS_Va1,nS_Vb1,cum_prop_va,cum_prop_vb)],by=c("country","date"),all.x=T)
# 
# # Cast incidence data table to wide format for different variants
# curr_inc_dt_wide = dcast(curr_inc_dt,country + date + age_group_model ~ vrnt, value.var = "nS_E")
# setnames(curr_inc_dt_wide,c("Other","Alpha","Delta"),paste0("nS_E",c("","2","3")))
# 
# # Merge prevalence, incidence and vaccination data tables
# curr_prev_inc_dt = merge(curr_prev_dt,curr_inc_dt_wide,by=c("country","date","age_group_model"))
# curr_dt = merge(curr_prev_inc_dt,curr_vax_dt,by=c("country","date","age_group_model"),all.x=T)
# 
# # Calculate maximum number of days on which exposures can occur if incidence and
# # vaccinations remain at their current levels
# curr_dt[,d:=S/(nS_E+nS_E2+nS_E3+nS_Va1+nS_Vb1)]
# # Calculate maximum remaining cumulative number of exposures
# curr_dt[,`:=`(cum_nS_E=d*nS_E,cum_nS_E2=d*nS_E2,cum_nS_E3=d*nS_E3)]
# 
# # Add symptomatic proportion to data table
# curr_dt = merge(curr_dt,covy_dt,by="age_group_model")

# Infection hospitalisation rate (derived from Salje et al., Science)
ihr = data.table(age = 0:85,ihr = exp(-7.37 + 0.068 * 0:85) / (1 + exp(-7.37 + 0.068 * 0:85)))
ihr[,age_group_model:=cut(age,c(min_ages_model,Inf),labels=agegroups_model,right=F)]
ihr = ihr[,.(ihr=mean(ihr)),by=.(age_group_model)]

# # Add infection-hospitalisation rate to data table
# curr_dt = merge(curr_dt,ihr,by="age_group_model")
# 
# # Add IFR to data table
# base_dt = merge(CJ(country=curr_dt[,unique(country)],age=0:100),pop[,.(country,age,population)],by=c("country","age"),all.x=T)
# ifr_dt = construct_ifr_data_table(ifr,base_dt,min_ages_model,agegroups_model)
# setnames(ifr_dt,"age_group","age_group_model")
# curr_dt = merge(curr_dt,ifr_dt,by=c("country","age_group_model"))
# 
# # Calculate maximum remaining numbers of symptomatic infections, hospitalisations and deaths for each variant
# curr_dt[,`:=`(cum_nE_Ip = y*cum_nS_E,
#               cum_nE2_Ip2 = y*cum_nS_E2,
#               cum_nE3_Ip3 = y*cum_nS_E3)]
# curr_dt[,`:=`(cum_hosp = ihr*cum_nE_Ip,
#               cum_hosp2 = ihr*cum_nE2_Ip2,
#               cum_hosp3 = ihr*cum_nE3_Ip3)]
# curr_dt[,`:=`(cum_deaths = ifr*cum_nE_Ip,
#               cum_deaths2 = ifr*cum_nE2_Ip2,
#               cum_deaths3 = ifr*cum_nE3_Ip3)]
# 
# # Calculate overall symptomatic infections, hospitalisations and deaths
# curr_dt[,`:=`(ovrl_cum_nE_Ip = cum_nE_Ip+cum_nE2_Ip2+cum_nE3_Ip3,
#               ovrl_cum_hosp = cum_hosp+cum_hosp2+cum_hosp3,
#               ovrl_cum_deaths = cum_deaths+cum_deaths2+cum_deaths3)]
# cols5 = c("ovrl_cum_nE_Ip","ovrl_cum_hosp","ovrl_cum_deaths")
# curr_dt[,sub("cum","cum_inc",cols5):=lapply(.SD,function(x) x/population),.SDcols=cols5]
# 
# curr_dt[,cum_prop_v := cum_prop_va+cum_prop_vb]
# curr_dt[,`:=`(cum_v=cum_prop_v*population,cum_u=(1-cum_prop_v)*population)]
# curr_dt[,`:=`(pop_prop_v=cum_v/sum(population),pop_prop_u=cum_u/sum(population)),by=.(country)]
# 
# curr_dt[,age_group_model:=factor(age_group_model,levels=agegroups_model)]
# 
# ggplot(curr_dt[country!="Greece"],aes(x=cum_prop_v,y=ovrl_cum_inc_hosp*1e5,color=age_group_model)) +
#     geom_point() +
#     labs(x="Proportion fully vaccinated",y="Maximum remaining hospitalisations/100,000 population") +
#     scale_y_log10() +
#     facet_wrap(~country)
# 
# # Plot vaccination coverage pyramids
# # ggplot(melt(curr_dt,measure.vars = c("cum_u","cum_v")),aes(x=value,y=age_group_model,alpha=variable)) +
# #     geom_col(position = "stack") +
# #     scale_alpha_manual(name="",values=c(0.4,1),labels=c("Unvaccinated/\npartially vaccinated","Fully vaccinated")) +
# #     labs(x="Population",y="Age group") +
# #     facet_wrap(~country)
# ggplot(melt(curr_dt[country!="Greece"],measure.vars = c("pop_prop_u","pop_prop_v")),aes(x=value,y=age_group_model,alpha=variable)) +
#     geom_col(position = "stack") +
#     scale_alpha_manual(name="",values=c(0.4,1),labels=c("Unvaccinated/\npartially vaccinated","Fully vaccinated")) +
#     labs(x="Proportion of population",y="Age group") +
#     facet_wrap(~country)
# ggsave(paste0(dir_out,"vax_cov_pop_pyramids.pdf"),width = 10,height = 8)
# 
# # Plot maximum remaining hospitalisations by age
# ggplot(curr_dt[country!="Greece"],aes(x=ovrl_cum_inc_hosp*1e5,y=age_group_model)) +
#     geom_col() +
#     labs(x="Maximum remaining hospitalisations/100,000 population",y="Age group") +
#     facet_wrap(~country)
# ggsave(paste0(dir_out,"rem_hosps_by_age.pdf"),width = 10,height = 8)
# 
# ovrl_curr_dt = curr_dt[,.(population=sum(population),
#                           cum_va=sum(cum_prop_va*population),
#                           cum_vb=sum(cum_prop_vb*population),
#                           ovrl_cum_nE_Ip=sum(ovrl_cum_nE_Ip),
#                           ovrl_cum_hosp=sum(ovrl_cum_hosp),
#                           ovrl_cum_deaths=sum(ovrl_cum_deaths)),
#                        by=.(country,date)]
# ovrl_curr_dt[,`:=`(cum_prop_va=cum_va/population,
#                   cum_prop_vb=cum_vb/population)]
# ovrl_curr_dt[,sub("cum","cum_inc",cols5):=lapply(.SD,function(x) x/population),.SDcols=cols5]
# ovrl_curr_dt[,cum_prop_v := cum_prop_va+cum_prop_vb]
# 
# ggplot(ovrl_curr_dt[country!="Greece"],aes(x=cum_prop_v,y=ovrl_cum_inc_hosp*1e5)) +
#     geom_point() +
#     geom_text(aes(label=country),hjust=-0.1,vjust=0) +
#     xlim(c(0.1,0.7)) +
#     labs(x = "Proportion fully vaccinated",y = "Maximum remaining hospitalisations/100,000 population") +
#     scale_y_log10()
# ggsave(paste0(dir_out,"rem_hosps_vs_prop_full_vax.pdf"),width = 10,height = 8)
# # ggplot(melt(ovrl_curr_dt,measure.vars = c("ovrl_cum_inc_hosp","ovrl_cum_inc_deaths")),aes(x=cum_prop_v,y=value*1e5)) +
# #     geom_point() +
# #     geom_text(aes(label=country),hjust=-0.1,vjust=0) +
# #     xlim(c(0.1,0.7)) +
# #     labs(x = "Proportion fully vaccinated",y = "Maximum remaining per 100,000 population") +
# #     scale_y_log10() +
# #     facet_wrap(~variable)


base_dt = merge(CJ(country=prev_dt[,unique(country)],age=0:100),pop[,.(country,age,population)],by=c("country","age"),all.x=T)
ifr_dt = construct_ifr_data_table(ifr,base_dt,min_ages_model,agegroups_model)
setnames(ifr_dt,"age_group","age_group_model")

res = calc_rem_burden(prev_dt,agegroups_model,ihr,ifr_dt,dHosp,dDeath)
proj_prev_dt = res$proj_prev_dt
rem_burden_dt = res$rem_burden_dt

# Plot projected infection incidence (should be constant)
ggplot(proj_prev_dt[date<as.Date("2023-01-01") & !(country %in% c("Finland","France","Greece"))],aes(x=date,y=nS_E+nV_E+nV_L+nV_R,color=age_group_model)) +
    geom_line() +
    facet_wrap(~country)

# Plot maximum remaining hospitalisations and deaths by age
ggplot(rem_burden_dt[!(country %in% c("Finland","France","Greece"))],aes(x=cum_inc_hosp*1e5,y=age_group_model)) +
    geom_col() +
    labs(x="Maximum remaining hospitalisations/100,000 population",y="Age group") +
    facet_wrap(~country)
ggsave(paste0(dir_out,"rem_hosps_by_age2.pdf"),width = 10,height = 8)

max_dates = prev_dt[,.(date=max(date)),by=.(country)]
curr_prev_dt = merge(max_dates,prev_dt,by=c("country","date"))
rem_burden_dt = merge(rem_burden_dt,curr_prev_dt[,.(country,date,age_group_model,cum_prop_va,cum_prop_vb)],by=c("country","age_group_model"))

# Plot maximum overall remaining hospitalisations and deaths
cols5 = c("cum_hosp","cum_deaths")
ovrl_rem_burden_dt = rem_burden_dt[,.(population=sum(population),
                          cum_va=sum(cum_prop_va*population),
                          cum_vb=sum(cum_prop_vb*population),
                          cum_hosp=sum(cum_hosp),
                          cum_deaths=sum(cum_deaths)),
                       by=.(country,date)]
ovrl_rem_burden_dt[,`:=`(cum_prop_va=cum_va/population,
                  cum_prop_vb=cum_vb/population)]
ovrl_rem_burden_dt[,sub("cum","cum_inc",cols5):=lapply(.SD,function(x) x/population),.SDcols=cols5]
ovrl_rem_burden_dt[,cum_prop_v := cum_prop_va+cum_prop_vb]

ggplot(ovrl_rem_burden_dt[!(country %in% c("Finland","France","Greece"))],aes(x=cum_prop_v,y=cum_inc_hosp*1e5)) +
    geom_point() +
    geom_text(aes(label=country),hjust=-0.1,vjust=0) +
    xlim(c(0.1,0.7)) +
    labs(x = "Proportion fully vaccinated",y = "Maximum remaining hospitalisations/100,000 population") +
    scale_y_log10()
ggsave(paste0(dir_out,"rem_hosps_vs_prop_full_vax2.pdf"),width = 10,height = 8)

# Plot vaccination coverage pyramids
rem_burden_dt[,cum_prop_v := cum_prop_va+cum_prop_vb]
rem_burden_dt[,`:=`(cum_v=cum_prop_v*population,cum_u=(1-cum_prop_v)*population)]
rem_burden_dt[,`:=`(pop_prop_v=cum_v/sum(population),pop_prop_u=cum_u/sum(population)),by=.(country)]
ggplot(melt(rem_burden_dt[!(country %in% c("Finland","France","Greece"))],measure.vars = c("pop_prop_u","pop_prop_v")),aes(x=value,y=age_group_model,alpha=variable)) +
    geom_col(position = "stack") +
    scale_alpha_manual(name="",values=c(0.4,1),labels=c("Unvaccinated/\npartially vaccinated","Fully vaccinated")) +
    labs(x="Proportion of population",y="Age group") +
    facet_wrap(~country)
ggsave(paste0(dir_out,"vax_cov_pop_pyramids2.pdf"),width = 10,height = 8)



# # # Make data tables of rates of waning immunity for different variants and vaccines
# # wn_dt = data.table(vrnt = inc_dt[,unique(vrnt)], wn = rep(log(0.85)/-182.5,inc_dt[,length(unique(vrnt))]))
# # wv_dt = CJ(dose = c(1,2),type = vax[,unique(type)])
# # wv_dt[,wv:=c(rep(0,vax[,length(unique(type))]),rep(log(0.85)/-182.5,vax[,length(unique(type))]))]
# 
# 
# 
# # dt_model = construct_data_table(agegroups_model,deaths,pop,cols,ltc_deaths,vax,num_type,ifr)
# 
# # # Make data table of unique combinations of country and age
# # countries = deaths[,unique(country)]
# # country_codes = deaths[,unique(country_code)]
# # base_dt = CJ(country_code = country_codes,age = 0:100)
# #
# # # Add population data
# # base_dt = merge(base_dt,pop[,.(iso3,age,population)],by.x=c("country_code","age"),by.y=c("iso3","age"),all.x=T)
# #
# # # Make a copy
# # base_deaths_dt = copy(base_dt)
# #
# # # Add age groups from deaths data
# # for (i in seq_along(country_codes)){
# #     agegroups_cntry = sort(deaths[country_code==country_codes[i],unique(age_group)])
# #     min_ages_cntry = get_min_age(agegroups_cntry)
# #     min_ages_cntry[is.na(min_ages_cntry)] = 0
# #
# #     base_deaths_dt[country_code==country_codes[i],age_group:=cut(age,c(min_ages_cntry,Inf),labels=agegroups_cntry,right=F)]
# # }
# #
# # # Merge with deaths data table
# # # N.B. This duplicates deaths for the same age group, so we divide deaths between
# # # ages according to population fraction
# # # TODO - revisit this and divide deaths between age groups by relative IFR rather
# # # than population fraction?
# # base_deaths_dt = merge(base_deaths_dt,deaths,by=c("country_code","age_group"),all=T,allow.cartesian = T)
# # base_deaths_dt[,(cols):=lapply(.SD,as.numeric),.SDcols=cols]
# # base_deaths_dt[,(cols):=lapply(.SD,function(x) x*population/sum(population)),.SDcols=cols,by=.(country_code,date,age_group)]
# #
# # # Check death totals match
# # print(deaths[!(country_code %in% c("ENW","SCO")),lapply(.SD,function(x) sum(x,na.rm=T)),.SDcols=cols])
# # print(base_deaths_dt[!(country_code %in% c("ENW","SCO")),lapply(.SD,function(x) sum(x,na.rm=T)),.SDcols=cols])
# #
# # # Construct data tables of deaths for deconvolution age groups and covidm age groups
# # deaths_dt_model = construct_deaths_data_table(base_deaths_dt,agegroups_model,min_date)
# # deaths_dt = construct_deaths_data_table(base_deaths_dt,agegroups,min_date)
# #
# # # # Change age groups
# # # base_deaths_dt[,age_group:=cut(age,c(min_ages,Inf),labels=agegroups,right=F)]
# # # # Sum deaths in each age group
# # # base_deaths_dt = base_deaths_dt[,lapply(.SD,sum),.SDcols=c("population",cols),by=.(region,country,country_code,date,age_group)]
# # #
# # # dates = seq.Date(min_date,base_deaths_dt[,max(date),by=.(country)][,min(V1)],by=1)
# # # deaths_dt = CJ(country_code = base_deaths_dt[,unique(country_code)], date = dates, age_group = agegroups)
# # # deaths_dt = merge(deaths_dt,base_deaths_dt[,!c("region","country","population")],by=c("country_code","date","age_group"),all.x=T)
# # #
# # # # Exclude country and variable combinations with all missing data
# # # deaths_dt_long = melt(deaths_dt,id.vars = c("country_code","date","age_group"),variable.name="sex",value.name="cum_death")
# # # deaths_dt_long[,sex:=sub("cum_death_","",sex)]
# # # country_codes_notallNA = deaths_dt_long[,!all(is.na(cum_death)),by=.(country_code,sex)][V1==T,.(country_code,sex)]
# # # deaths_dt_long = deaths_dt_long[country_codes_notallNA,on=c("country_code","sex")]
# # #
# # # # Fill in cumulative deaths for earliest date so there is a value to interpolate from
# # # deaths_dt_long[date == min_date,cum_death:=0]
# # #
# # # # Interpolate cumulative deaths
# # # deaths_dt_long[,cum_death_i := approx(date,cum_death,dates)$y,by=.(country_code,sex,age_group)]
# # # # deaths_dt[,(paste0(cols,"_i")) := lapply(.SD,function(yi) approx(date,yi,dates)$y),.SDcols=cols,by=.(country_code,age_group)]
# # #
# # # # Calculate new daily deaths
# # # # FOR NOW - as large negative counts have been removed, treat remaining small
# # # # negative counts from differencing as 0s
# # # deaths_dt_long[,death_i := c(0,pmax(diff(cum_death_i),0)),by=.(country_code,sex,age_group)]
# # #
# # # # Cast to wide format
# # # deaths_dt = dcast(deaths_dt_long,country_code + date + age_group ~ sex, value.var = c("cum_death","cum_death_i","death_i"))
# # # # Add region, country and country code
# # # deaths_dt = merge(deaths_dt,unique(deaths[,.(region,country,country_code)]),by="country_code")
# #
# # # Plot data
# # # cumulative deaths
# # ggplot(deaths_dt,aes(x=date,y=cum_death_i_both,group=age_group,color=age_group)) +
# #     geom_line() +
# #     scale_y_log10() +
# #     facet_wrap(~country_code)
# # # new deaths
# # ggplot(deaths_dt[country_code != "USA"],aes(x=date,y=death_i_both,group=age_group,color=age_group)) +
# #     geom_line() +
# #     facet_wrap(~country_code)
# #
# # # LTC death data
# # ltc_deaths = fread("../ltccovid_data/ltc_deaths.csv")
# # names(ltc_deaths) = tolower(names(ltc_deaths))
# # names(ltc_deaths)[3:8] = c("approach","total_deaths","total_ltc_res_deaths","total_ltc_deaths","perc_ltc_res_deaths","perc_ltc_deaths")
# # cols5 = c("total_deaths","total_ltc_res_deaths","total_ltc_deaths")
# # ltc_deaths[,(cols5):=lapply(.SD,function(x) as.numeric(sub(",","",x))),.SDcols = cols5]
# # ltc_deaths[,`:=`(prop_ltc_res_deaths=total_ltc_res_deaths/total_deaths,prop_ltc_deaths=total_ltc_deaths/total_deaths)]
# #
# # # Merge with deaths data table
# # deaths_dt = merge(deaths_dt,ltc_deaths[,.(country,prop_ltc_res_deaths,prop_ltc_deaths)],by="country")
# #
# # # Calculate LTC deaths among those 60+ by multiplying total deaths by LTC proportion
# # # N.B. relies on having proportion of deaths in LTCs or LTC residents
# # deaths_dt[age_group %in% agegroups[max_ages<60],death_i_both_ltc := 0] # assume no deaths in under-60s
# # deaths_dt[age_group %in% agegroups[max_ages>60],death_i_both_ltc := death_i_both * fifelse(!is.na(prop_ltc_res_deaths),prop_ltc_res_deaths,prop_ltc_deaths)]
# # deaths_dt[,death_i_both_non_ltc := death_i_both - death_i_both_ltc]
# #
# # # VACCINATION DATA
# #
# # # Read in ECDC vaccination data
# # vax = fread("../ecdc_data/ecdc_vaccination_data.csv")
# #
# # # Remove non-country-level data
# # vax = vax[Region==ReportingCountry,]
# #
# # # Change country code for Greece
# # vax[ReportingCountry=="EL", ReportingCountry:="GR"]
# #
# # # Add country names
# # vax[,country:=country_iso_codes[match(vax[,ReportingCountry],iso_code),country]]
# #
# # # Convert ISO week to date
# # vax[,date:=ISOweek2date(paste0(YearWeekISO,"-1"))]
# # # Lose last 3 weeks of data as incomplete at time of download (current week may
# # # not be finished and delay in reporting of 1-2 weeks)
# # max_date = vax[,max(date)]
# # vax = vax[!(date %in% (max_date-7*(0:2)))]
# #
# # # Drop data for HCWs and LTCF residents as this is contained in age group data
# # vax = vax[!(TargetGroup %in% c("HCW","LTCF"))]
# # vax[,age_group:=sub("_","-",sub("Age|1_Age","",TargetGroup))]
# #
# # # Drop data disaggregated by age for those <18 for now and data aggregated for
# # # <60 and 60+
# # # TODO - revisit this
# # # print(vax[age_group=="UNK",sum(FirstDose+SecondDose+UnknownDose),by=.(ReportingCountry)])
# # vax = vax[!(age_group %in% c("0-4","5-9","10-14","15-17","<60","60+","ALL"))]
# # vax[age_group=="<18", age_group := "0-17"]
# #
# # # Lose unneeded variables
# # setnames(vax,c("Vaccine","ReportingCountry"),c("vaccine","country_code2"))
# # # TODO - check population denominators against UN estimates
# # vax = vax[,.(country,country_code2,year_week_iso=YearWeekISO,date,age_group,vaccine,first=FirstDose,second=SecondDose)] #,population=Denominator
# #
# # # Treat Pfizer and Moderna as vaccine B, and all others (AstraZeneca, Janssen,
# # # Sputnik, Beijing CBNG) as vaccine A
# # vax[,type:=fifelse(vaccine %in% c("COM","MOD"),"vb",fifelse(vaccine!="UNK","va","vu"))]
# #
# # # Calculate proportions of each vaccine type
# # # TODO - calculate over time
# # num_type = dcast(vax[,.(first=sum(first)),by=.(country,type)],country~type)
# # # Split vaccines of unknown type according to average proportion of each type in
# # # other countries FOR NOW - can correct with OWID data
# # cols3 = c("va","vb","vu")
# # num_type[,(cols3):=lapply(.SD,as.numeric),.SDcols=cols3]
# # num_type[!is.na(vu),`:=`(va=fifelse(is.na(va),vu*num_type[is.na(vu),mean(va/(va+vb))],va+vu*num_type[is.na(vu),mean(va/(va+vb))]),
# #                          vb=fifelse(is.na(va),vu*num_type[is.na(vu),mean(vb/(va+vb))],vb+vu*num_type[is.na(vu),mean(vb/(va+vb))]))]
# # num_type[,vu:=NULL]
# # num_type[,`:=`(prop_va=va/(va+vb),prop_vb=vb/(va+vb))]
# #
# # # FOR NOW - sum over vaccine type
# # vax = vax[,.(first = sum(first), second = sum(second)), by = .(country,country_code2,year_week_iso,date,age_group)] #,population
# #
# # # # Calculate proportion vaccinated
# # # vax[,prop_v := second/population]
# #
# # # See how many doses have unknown age group
# # print(vax[,sum(first+second)]) # 209715607
# # print(vax[age_group!="UNK",sum(first+second)]) #209476181
# # vax[,.(first=sum(first),second=sum(second)), by = .(age_group)] # only about 240,000 out of about 210 million, so ignore FOR NOW
# #
# # # Get vaccination age groups from vax
# # agegroups_vax = sort(vax[,unique(age_group)])
# # agegroups_vax = agegroups_vax[agegroups_vax!="UNK"]
# # min_ages_vax = get_min_age(agegroups_vax)
# # min_ages_vax[is.na(min_ages_vax)] = 0
# # max_ages_vax = get_max_age(agegroups_vax)
# # max_ages_vax[is.na(max_ages_vax)] = Inf
# #
# # # Make daily date sequence from earliest death date to latest date for which
# # # vaccine data is available - length should be a multiple of 7
# # Ab_delay = 14 # days
# # dates1 = seq.Date(min_date-ceiling(Ab_delay/7)*7,vax[,max(date)+6],by=1) # -Ab_delay to account for Ab_delay-day delay to Ab development, and max(date) + 6 to get end of last ISO week
# #
# # # Create data table with all combinations of countries, dates, vaccines and ages
# # # to store vaccine schedule
# # vax_dt = CJ(country=vax[,unique(country)],date=dates1,age=0:100)
# # vax_dt = merge(vax_dt,pop[,.(country,age,population)],by=c("country","age"),all.x=T)
# # vax_dt[,year_week_iso:=ISOweek(date)]
# # vax_dt[,age_group:=cut(age,c(min_ages_vax,Inf),labels=agegroups_vax,right=F)]
# # vax_dt[,country_code2:=country_iso_codes[match(vax_dt[,country],country),iso_code]]
# #
# # # Merge with vax data table
# # # N.B. This duplicates numbers of doses for the same ISO week and age group, so
# # # we then divide by 7 to get the average doses per day and divide doses between
# # # age groups according to population fraction
# # vax_dt = merge(vax_dt,vax[,!"date"],by=c("country","country_code2","year_week_iso","age_group"),all.x=T)
# # cols1 = c("first","second")
# # vax_dt[,(cols1):=lapply(.SD,as.numeric),.SDcols=cols1]
# # vax_dt[,(cols1):=lapply(.SD,function(x) x/7),.SDcols=cols1]
# # vax_dt[,(cols1):=lapply(.SD,function(x) x*population/sum(population)),.SDcols=cols1,by=.(country,date,age_group)]
# # print(vax_dt[,sum(first+second,na.rm=T)]) #209476181
# #
# # # Change age groups
# # vax_dt[,age_group:=cut(age,c(min_ages,Inf),agegroups,right=F)]
# # # Sum vaccinations in each age group
# # vax_dt = vax_dt[,lapply(.SD,function(x) sum(x,na.rm=T)),.SDcols=c("population",cols1),by=.(country,date,age_group)]
# # # Calculate proportion fully vaccinated
# # vax_dt[,prop_v:=second/population]
# # # Calculate cumulative proportion fully vaccinated ensuring data is in chronological date order
# # setorder(vax_dt,country,date,age_group)
# # vax_dt[,cum_prop_v:=cumsum(prop_v),by=.(country,age_group)]
# #
# # # Plot data
# # # cumulative proportion fully vaccinated
# # ggplot(vax_dt[date>=vax[,min(date)]],aes(x=date,y=cum_prop_v,group=age_group,color=age_group)) +
# #     geom_line() +
# #     facet_wrap(~country)
# #
# # # proportion newly fully vaccinated
# # ggplot(vax_dt[date>=vax[,min(date)]],aes(x=date,y=prop_v,group=age_group,color=age_group)) +
# #     geom_line() +
# #     facet_wrap(~country)
# #
# #
# # # Merge death and vaccination data
# # # See which countries are in both datasets
# # print(intersect(deaths_dt[,country],vax_dt[,country]))
# #
# # # Add Ab_delay days to date in vaccination data for development of Ab
# # vax_dt[,date_v:=date]
# # vax_dt[,date:=date+Ab_delay]
# # dates2 = as.Date(intersect(deaths_dt[,date],vax_dt[,date]),origin="1970-01-01")
# # deaths_vax_dt = merge(deaths_dt[date %in% dates2],vax_dt[date %in% dates2],by=c("country","date","age_group"),all.x=T)
# #
# # # IFR ESTIMATES
# # ifr_dt = copy(base_dt)
# #
# # # # Add IFR estimate from Levin et al Eur Jrnl Epi 2020
# # # ifr_dt[,`:=`(ifr=10^(-3.27+0.0524*age)/100)]
# # # # Calculate 95% CI for IFR assuming FOR NOW no covariance between intercept and
# # # # slope estimates
# # # # TODO - estimate intercept and slope covariance from lower and upper bounds in
# # # # supplementary spreadsheet in Levin et al 2020
# # # se_intcpt = 0.07*log(10)
# # # se_slope = 0.0013*log(10)
# # # ifr_dt[,se_ifr:=sqrt(se_intcpt^2+se_slope^2*age^2)*ifr]
# # # ifr_dt[,`:=`(ifr_lb=pmax(0,ifr-qnorm(0.975)*se_ifr),ifr_ub=pmin(1,ifr+qnorm(0.975)*se_ifr))]
# # #
# # # ifr_dt[,age_group:=cut(age,c(min_ages,Inf),labels=agegroups,right=F)]
# # # cols2 = c("ifr","se_ifr","ifr_lb","ifr_ub")
# # # ifr_dt=ifr_dt[,lapply(.SD,function(x) sum(x*population)/sum(population)),.SDcols=cols2,by=.(country_code,age_group)]
# #
# # # IFR estimate from O'Driscoll et al Nature 2020
# # ifr = fread(datapath("IFR_by_age_ODriscoll.csv"))
# # names(ifr) = tolower(names(ifr))
# # min_ages_ifr = get_min_age(ifr[,age_group])
# #
# # # Add IFR age groups
# # ifr_dt[,age_group:=cut(age,c(min_ages_ifr,Inf),labels=ifr[,age_group],right=F)]
# #
# # # Merge with IFR data table
# # ifr_dt = merge(ifr_dt,ifr[,.(age_group,ifr=median_perc_mean/100,ifr_lb=ci_95_lb_mean/100,ifr_ub=ci_95_ub_mean/100)],by="age_group")
# #
# # # Change age groups
# # ifr_dt[,age_group := cut(age,c(min_ages,Inf),labels=agegroups,right=F)]
# #
# # # Calculate population-weighted average for each age group
# # cols2 = c("ifr","ifr_lb","ifr_ub")
# # ifr_dt = ifr_dt[,lapply(.SD,function(x) sum(x*population)/sum(population)),.SDcols=cols2,by=c("country_code","age_group")]
# #
# # # Plot country IFRs
# # ggplot(ifr_dt,aes(x=age_group,y=ifr,group=country_code,color=country_code)) +
# #     geom_line() #+scale_y_log10()
# #
# # # Merge IFR with death and vaccination data
# # dt = merge(deaths_vax_dt,ifr_dt,by=c("country_code","age_group"))
# # dt = merge(dt,num_type[,.(country,prop_va,prop_vb)],by="country")
# #
# # # Vaccine efficacy parameters
# # # # eiX_vYZ = efficacy of dose Z of vaccine Y against infection with strain X
# # # ei_va2 = 0.68
# # # ei2_va2 = 0.68
# # # ei3_va2 = 0.6154
# # # ei_vb2 = 0.85
# # # ei2_vb2 = 0.85
# # # ei3_vb2 = 0.7999
# # # ed_vYZiX = efficacy of dose Z of vaccine Y against disease given infection with strain X
# # ed_va2i = 0.3125
# # ed_va2i2 = 0.3125
# # ed_va2i3 = 0.2353094
# # ed_vb2i = 0.2667
# # ed_vb2i2 = 0.2667
# # ed_vb2i3 = 0.187906
# # # em_vYZdX = efficacy of dose Z of vaccine Y against death given disease from strain X
# # em_va2d = 0.77
# # em_va2d2 = 0.77
# # em_va2d3 = 0.6766406
# # em_vb2d = 0.55
# # em_vb2d2 = 0.55
# # em_vb2d3 = 0.52
# #
# # # FOR NOW - use efficacy values for wild-type variant
# # # TODO - use proportion of each vaccine by age group and time and proportion of
# # # each variant over time to calculate more detailed age- and time-dependent IFR
# # # dt[,`:=`(ei=prop_va*ei_va2+prop_vb*ei_vb2,ed=prop_va*ed_va2i+prop_vb*ed_vb2i,em=prop_va*em_va2d+prop_vb*em_vb2d)]
# # dt[,`:=`(ed=prop_va*ed_va2i+prop_vb*ed_vb2i,em=prop_va*em_va2d+prop_vb*em_vb2d)]
# #
# # # Calculate IFR under vaccination
# # dt[,ifr_v:=(1-ed)*(1-em)*ifr]
# #
# # # Calculate IFR over time
# # dt[,ifr_t:=(1-cum_prop_v)*ifr+cum_prop_v*ifr_v]
# #
# # # Plot population-weighted average IFR over time for all countries
# # ggplot(dt_out[,.(ifr_t=sum(ifr_t*population)/sum(population)),by=.(country,date)],
# #        aes(x=date,y=ifr_t,group=country,color=country)) + geom_line()
# # ggsave("./output/avg_ifr_over_time.pdf",width = 5,height = 4)
# 
# # tstart = Sys.time()
# # for (i in seq_along(country_codes1)){
# #     print(i)
# #     deaths_i = dt[country_code==country_codes1[i]]
# #
# #     deaths_i_wide = dcast(deaths_i,date ~ age_group,value.var = "death_i_both")
# #     dates = deaths_i_wide[,date] # rename dates
# #     deaths_i_wide = as.matrix(deaths_i_wide[,!"date"],rownames.value=as.character(dates))
# #
# #     # BACK PROJECTION
# #     # Convert to sts object
# #     deaths_i_sts = sts(deaths_i_wide,start=c(lubridate::year(min(dates)),yday(min(dates))),frequency = 365)
# #     # # Plot death time series
# #     # plot(deaths_i_sts,type=observed~time)
# #
# #     # Deconvolve death curve to infection curve of those that died
# #     control = list(k = 0,eps = c(0.01,2),Tmark = nrow(deaths_i_sts) - mean_idd,alpha=0.05,eq3a.method = "C")
# #     idd_pmf_mat = matrix(idd_pmf,nrow = length(idd_pmf),ncol = ncol(deaths_i_sts))
# #     deaths_i_stsBP = backprojNP(deaths_i_sts,idd_pmf_mat,control = control)
# #     # plot(deaths_i_stsBP,xaxis.labelFormat=NULL,legend=NULL,lwd=c(1,1,2),lty=c(1,1,1),ylim=c(0,100),main="")
# #
# #     exposures_dead = data.table(deaths_i_stsBP@upperbound)
# #     names(exposures_dead) = agegroups
# #     exposures_dead[,date:=dates]
# #     exposures_dead_long = melt(exposures_dead,id.vars = "date",variable.name = "age_group",value.name = "exposures_dead")
# #     exposures_dead_long[,country_code:=country_codes1[i]]
# #
# #     exposures_dead_list[[i]] = exposures_dead_long
# #
# #     # # RIDE ALGORITHM
# #     # exposures_dead_ride_mat = matrix(nrow = nrow(deaths_i_wide),ncol = ncol(deaths_i_wide))
# #     # for (j in 1:ncol(deaths_i_wide)){
# #     #     # First input to fit_incidence needs to be an integer vector
# #     #     exposures_model = fit_incidence(as.integer(round(deaths_i_sts@observed[,j])),idd_pmf[2:(idd_max+1)])
# #     #     exposures_dead_ride_mat[,j] = exposures_model$Ihat
# #     #
# #     #     # TODO - get 95% CrI from exposures_model$Isamps
# #     #
# #     # }
# #     # exposures_dead_ride = data.table(exposures_dead_ride_mat)
# #     # names(exposures_dead_ride) = agegroups
# #     # exposures_dead_ride[,date:=dates]
# #     # exposures_dead_ride_long = melt(exposures_dead_ride,id.vars = "date",variable.name = "age_group",value.name = "exposures_dead")
# #     #
# #     # exposures_dead_ride_list[[i]] = exposures_dead_ride_long
# # }
# # tend = Sys.time()
# # print(tend - tstart)
