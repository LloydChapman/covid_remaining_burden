library(data.table)
library(qs)
library(ggplot2)
library(osfr)
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
# library(truncnorm)
library(cowplot)
library(patchwork)
library(stringr)

source("./backcalculation_functions.R")
# source("./initial_conditions_functions.R")
# source("./simulation_functions.R")
source("./process_data.R")
source("./derive_initial_conditions.R")
source("./calculate_remaining_burden.R")


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
theme_set(cowplot::theme_cowplot(font_size = 12) + theme(strip.background = element_blank()))

# Download death data
download_death_data("coverage")
download_death_data("ined")

# Set source of age-stratified death data
source_deaths = "coverage"

# Set deconvolution method
method = "ride"

# Set age groups for deconvolution
agegroups = c("0-39","40-49","50-59","60-69","70-79","80+") #c("0-9","10-19","20-29","30-39","40-49","50-59","60-69","70-79","80+")
min_ages = get_min_age(agegroups)
max_ages = get_max_age(agegroups)
max_ages[is.na(max_ages)] = Inf

# Set age groups to disaggregate deconvolution output over
agegroups_model = factor(agegroups,levels = agegroups)

# Get country iso codes
covid_data_path = "./fitting_data/"
datapath = function(x) paste0(covid_data_path, x)
country_iso_codes = readRDS(datapath("country_iso_codes.rds"))

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

# Get UK population
popUK = readRDS(datapath("popNHS.rds"))

# Make population data table for England
popENG = CJ(country="England",age=0:100)
agegroups_pop = popUK[,unique(age)]
min_ages_pop = get_min_age(agegroups_pop)
popENG[,age_group:=cut(age,c(min_ages_pop,Inf),labels=agegroups_pop,right=F)]
popENG = merge(popENG,popUK,by.x=c("country","age_group"),by.y=c("name","age"))
popENG = popENG[,.(country,age,
                   female=1000*fifelse(age<90,f/5,f/10),
                   male=1000*fifelse(age<90,m/5,m/10))]
popENG[,population:=female+male]

# Bind to UN population data table
pop = rbind(pop,popENG,fill=T)

# Set delay distributions:
# latent period
dE  = cm_delay_gamma(2.5, 2.5, t_max = 15, t_step = 1)$p
# presymptomatic period
dIp = cm_delay_gamma(2.5, 4.0, t_max = 15, t_step = 1)$p
# symptomatic period
dIs = cm_delay_gamma(2.5, 4.0, t_max = 15, t_step = 1)$p
# asymptomatic period
dIa = cm_delay_gamma(5.0, 4.0, t_max = 15, t_step = 1)$p

# Get age-dependent susceptibility and symptomatic fraction
covid_scenario = qread(datapath("2-linelist_both_fit_fIa0.5-rbzvih.qs"));
colsu = names(covid_scenario)[grep("u_",names(covid_scenario))]
covu = colMeans(covid_scenario[,..colsu])
colsy = names(covid_scenario)[grep("y_",names(covid_scenario))]
covy = colMeans(covid_scenario[,..colsy])

# Get vaccine efficacies
ve_params = get_vaccine_efficacies()

# Set delays from vaccination to development of Ab
Ab_delay1 = 28 # delay after 1st dose in days
Ab_delay2 = 14 # delay after 2nd dose in days

# Read in ensemble IFR estimate from O'Driscoll et al Nature 2020
ifr = fread(datapath("IFR_by_age_ODriscoll.csv"))

# Frailty index for relative frailty of LTC residents compared to general population
frlty_idx = 3.8

# 
# DATA PROCESSING
# 


out = process_data(source_deaths,country_iso_codes,agegroups,pop,Ab_delay1,Ab_delay2,ifr,ve_params,dir_out)
dt = out$dt
vax = out$vax
vaxENG = out$vaxENG
vrnt_prop = out$vrnt_prop
rm(out)

# See how many doses have unknown age group
print(vax[,sum(count)])
print(vax[age_group!="UNK",sum(count)])
print(vax[country %in% dt[,unique(country)],.(count=sum(count)), by = .(age_group)][,.(age_group,p_UNK=count/sum(count))]) # only small % so ignore FOR NOW
print(vax[country %in% dt[,unique(country)],.(count=sum(count)),by=.(country,age_group)][,.(age_group,p_UNK=count/sum(count)),by=.(country)][age_group=="UNK"])


# 
# BACKCALCULATION
# 


# Convolve distributions to get incubation period and exposure-to-death delay distribution
# ip_params = readRDS("incubation_period.rds")
# ip_mean = convert_to_nat_mean(ip_params$mean,ip_params$sd)
# ip_sd = convert_to_nat_sd(ip_params$mean,ip_params$sd)
# dIncub = cm_delay_lnorm(ip_mean,ip_sd/ip_mean,t_max = 30,t_step = 1)$p
# 
# od_params = readRDS("onset_to_death_delay.rds")
# od_mean = convert_to_nat_mean(od_params$mean,od_params$sd)
# od_sd = convert_to_nat_sd(od_params$mean,od_params$sd)
# dDeath = disc_conv(cm_delay_lnorm(ip_mean,ip_sd/ip_mean,t_max = 60,t_step = 1)$p,cm_delay_lnorm(od_mean,od_sd/od_mean,t_max = 60,t_step = 1)$p)
dIncub = disc_conv(cm_delay_gamma(2.5, 2.5, t_max = 30, t_step = 1)$p,cm_delay_gamma(2.5, 4.0, t_max = 30, t_step = 1)$p)
dDeath = disc_conv(cm_delay_gamma(2.5, 2.5, t_max = 60, t_step = 1)$p,cm_delay_gamma(15, 2.2, t_max = 60, t_step = 1)$p)

# Normalise to ensure distributions sum to 1
dIncub = dIncub/sum(dIncub)
dDeath = dDeath/sum(dDeath)

## Backcalculate IFR-scaled infections
out = run_backcalculation(dt,dDeath,dIncub,frlty_idx,method = method)
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


# Plot output against age-stratified case data and seroprevalence data for validation
# Read in ECDC age-stratified case data for comparison with inferred infection time series
ecdc_cases_by_age = get_data("https://opendata.ecdc.europa.eu/covid19/agecasesnational/csv","csv",
                             paste0("../ecdc_data/",date_fitting,"/"),"ecdc_cases_by_age.csv")

# Read in and process seroprevalence data from SeroTracker for comparison with estimated cumulative proportion infected
serotracker_data = fread("../serotracker_data/SeroTracker Serosurveys Reporting Prevalence.csv")
sero_data = process_seroprevalence_data(serotracker_data)

# Remove functions
rm(list=lsf.str())

# Save output
save.image(paste0(dir_out,"backcalculation_output.RData"))

# Plot backcalculation output
plot_output(paste0(dir_out,"backcalculation_output.RData"),ecdc_cases_by_age,sero_data)


#
# INITIAL CONDITIONS CALCULATION
#


# Calculate current numbers of susceptible and uninfected vaccinated individuals in each 
# age group in each country
prev_dt = derive_initial_conditions(paste0(dir_out,"backcalculation_output.RData"))

# Save
saveRDS(prev_dt,paste0(dir_out,"init_cond_output.RDS"))


#
# REMAINING BURDEN CALCULATION
#


# Calculate remaining burden of hospitalisations and deaths assuming that entire
# population is exposed now
out = calculate_remaining_burden(paste0(dir_out,"init_cond_output.RDS"),agegroups_model,pop,ifr,frlty_idx)
rem_burden_dt = out$rem_burden_dt
ovrl_rem_burden_dt = out$ovrl_rem_burden_dt

# Save
saveRDS(rem_burden_dt,paste0(dir_out,"rem_burden_output.RDS"))
saveRDS(ovrl_rem_burden_dt,paste0(dir_out,"ovrl_rem_burden_output.RDS"))


# Plot
plot_remaining_burden(rem_burden_dt,ovrl_rem_burden_dt)
