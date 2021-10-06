library(data.table)
library(qs)
library(ggplot2)

source("./backcalculation_functions.R")

# 
# SETUP
# 

# Load initial conditions calculation output
load("./output/2021-10-05/init_cond_output.RData")


#
# REMAINING BURDEN CALCULATION
#


# CALCULATIONS ASSUMING ALL REMAINING SUSCEPTIBLES INFECTED NOW AND CORRESPONDING PROPORTION OF BREAKTHROUGH INFECTIONS AMONG VACCINATED INDIVIDUALS

# Infection hospitalisation rate (derived from Salje et al., Science)
ihr = data.table(age = 0:85,ihr = exp(-7.37 + 0.068 * 0:85) / (1 + exp(-7.37 + 0.068 * 0:85)))
ihr[,age_group_model:=cut(age,c(min_ages_model,Inf),labels=agegroups_model,right=F)]
ihr = ihr[,.(ihr=mean(ihr)),by=.(age_group_model)]

# Add IFR to data table
base_dt = merge(CJ(country=prev_dt[,unique(country)],age=0:100),pop[,.(country,age,population)],by=c("country","age"),all.x=T)
ifr_dt = construct_ifr_data_table(ifr,base_dt,min_ages_model,agegroups_model)
setnames(ifr_dt,"age_group","age_group_model")

# Calculate remaining burden of hospitalisations and deaths
rem_burden_dt = calc_rem_burden(prev_dt,ihr,ifr_dt)


# 
# PLOTTING 
# 


# Plot vaccine coverage pyramids
rem_burden_dt[,`:=`(cum_v=cum_prop_v*population,cum_u=(1-cum_prop_v)*population)]
rem_burden_dt[,`:=`(pop_prop_v=cum_v/sum(population),pop_prop_u=cum_u/sum(population)),by=.(country)]

ggplot(melt(rem_burden_dt[country!="Greece"],measure.vars = c("pop_prop_u","pop_prop_v")),aes(x=value,y=age_group_model,alpha=variable)) +
    geom_col(position = "stack") +
    scale_alpha_manual(name="",values=c(0.4,1),labels=c("Unvaccinated/\npartially vaccinated","Fully vaccinated")) +
    labs(x="Proportion of population",y="Age group") +
    facet_wrap(~country)
ggsave(paste0(dir_out,"vax_cov_pop_pyramids.png"),width = 10,height = 8)

# Plot maximum remaining hospitalisations by age
ggplot(rem_burden_dt[country!="Greece"],aes(x=cum_inc_hosp*1e5,y=age_group_model)) +
    geom_col() +
    labs(x="Maximum remaining hospitalisations/100,000 population",y="Age group") +
    facet_wrap(~country)
ggsave(paste0(dir_out,"rem_hosps_by_age.png"),width = 10,height = 8)

# Plot maximum overall remaining hospitalisations against vaccine coverage
ovrl_rem_burden_dt = rem_burden_dt[,.(population=sum(population),
                          cum_v=sum(cum_prop_v*population,na.rm=T),
                          cum_hosp=sum(cum_hosp,na.rm=T),
                          cum_deaths=sum(cum_deaths,na.rm=T)),
                       by=.(country,date)]
ovrl_rem_burden_dt[,cum_prop_v:=cum_v/population]
cols = c("cum_hosp","cum_deaths")
ovrl_rem_burden_dt[,sub("cum","cum_inc",cols):=lapply(.SD,function(x) x/population),.SDcols=cols]

ggplot(ovrl_rem_burden_dt[country!="Greece"],aes(x=cum_prop_v,y=cum_inc_hosp*1e5)) +
    geom_point() +
    geom_text(aes(label=country),hjust=-0.1,vjust=0) +
    xlim(c(0.1,0.7)) +
    labs(x = "Proportion fully vaccinated",y = "Maximum remaining hospitalisations/100,000 population") +
    scale_y_log10()
ggsave(paste0(dir_out,"rem_hosps_vs_prop_full_vax.png"),width = 10,height = 8)


