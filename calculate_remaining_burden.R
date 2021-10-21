library(data.table)
library(qs)
library(ggplot2)
library(cowplot)
library(scales)

# source("./backcalculation_functions.R")

# 
# SETUP
# 


# source("./derive_initial_conditions.R")

# Load initial conditions calculation output
load("./output/2021-10-20/init_cond_output.RData")

source("./backcalculation_functions.R")

# Set plot theme
theme_set(cowplot::theme_cowplot(font_size = 10) + theme(strip.background = element_blank()))

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
res = calc_rem_burden(prev_dt,ihr,ifr_dt)

# Calculate 
cols = c("population","cum_prop_v",
         "cum_hosp","cum_deaths","cum_inc_hosp","cum_inc_deaths",
         "cum_hosp_u","cum_deaths_u","cum_inc_hosp_u","cum_inc_deaths_u",
         "cum_hosp_v","cum_deaths_v","cum_inc_hosp_v","cum_inc_deaths_v","cum_prop_exp")
med = res[,lapply(.SD,median),.SDcols=cols,by=.(country,age_group_model,date)]
q95l = res[,lapply(.SD,function(x) quantile(x,probs=0.025)),.SDcols=cols,by=.(country,age_group_model,date)]
q95u = res[,lapply(.SD,function(x) quantile(x,probs=0.975)),.SDcols=cols,by=.(country,age_group_model,date)]
rem_burden_dt = merge(med,q95l,by=c("country","age_group_model","date"),suffixes = c("","_q95l"))
rem_burden_dt = merge(rem_burden_dt,q95u,by=c("country","age_group_model","date"),suffixes = c("","_q95u"))


# 
# PLOTTING 
# 


# Calculate cumulative numbers vaccinated and unvaccinated, and corresponding population proportions
rem_burden_dt[,`:=`(cum_v=cum_prop_v*population,cum_u=(1-cum_prop_v)*population)]
rem_burden_dt[,`:=`(pop_prop_v=cum_v/sum(population),pop_prop_u=cum_u/sum(population)),by=.(country)]

# Plot vaccine coverage pyramids
p1 = ggplot(melt(rem_burden_dt,measure.vars = c("pop_prop_u","pop_prop_v")),aes(x=value,y=age_group_model,alpha=variable)) +
    geom_col(position = "stack") +
    scale_alpha_manual(name="",values=c(0.4,1),labels=c("Unvaccinated","Partially/fully\nvaccinated")) +
    labs(x="Proportion of population",y="Age group") +
    facet_wrap(~country)
# ggsave(paste0(dir_out,"vax_cov_pop_pyramids.png"),width = 10,height = 8)

# Plot maximum remaining hospitalisations and deaths by age
# p2 = ggplot(rem_burden_dt,aes(y=age_group_model)) +
#     geom_col(aes(x=cum_inc_hosp*1e5),alpha=0.5) +
#     geom_errorbarh(aes(xmin=cum_inc_hosp_q95l*1e5,xmax=cum_inc_hosp_q95u*1e5),height=0.5) +
#     labs(x="Maximum remaining hospitalisations/100,000 population",y="Age group") +
#     scale_x_log10() +
#     facet_wrap(~country)
p2 = ggplot(melt(rem_burden_dt,measure.vars = c("cum_inc_hosp_u","cum_inc_hosp_v")),aes(y=age_group_model)) +
    geom_col(aes(x=value*1e5,alpha=variable),position="stack") +
    geom_errorbarh(aes(xmin=cum_inc_hosp_q95l*1e5,xmax=cum_inc_hosp_q95u*1e5),height=0.5) +
    scale_alpha_manual(name="",values=c(0.4,1),labels=c("Unvaccinated","Partially/fully\nvaccinated")) +
    labs(x="Maximum remaining hospitalisations/100,000 population",y="Age group") +
    facet_wrap(~country)
# ggsave(paste0(dir_out,"rem_hosps_by_age.png"),width = 10,height = 8)

# p3 = ggplot(rem_burden_dt,aes(y=age_group_model)) +
#     geom_col(aes(x=cum_inc_deaths*1e5),alpha=0.5) +
#     geom_errorbarh(aes(xmin=cum_inc_deaths_q95l*1e5,xmax=cum_inc_deaths_q95u*1e5),height=0.5) +
#     labs(x="Maximum remaining deaths/100,000 population",y="Age group") +
#     scale_x_log10(breaks=c(0.1,1,10,100,1000),labels=c("0.1","1","10","100","1000")) +
#     facet_wrap(~country)
p3 = ggplot(melt(rem_burden_dt,measure.vars = c("cum_inc_deaths_u","cum_inc_deaths_v")),aes(y=age_group_model)) +
    geom_col(aes(x=value*1e5,alpha=variable),position="stack") +
    geom_errorbarh(aes(xmin=cum_inc_deaths_q95l*1e5,xmax=cum_inc_deaths_q95u*1e5),height=0.5) +
    scale_alpha_manual(name="",values=c(0.4,1),labels=c("Unvaccinated","Partially/fully\nvaccinated")) +
    labs(x="Maximum remaining deaths/100,000 population",y="Age group") +
    facet_wrap(~country)
# ggsave(paste0(dir_out,"rem_deaths_by_age.png"),width = 10,height = 8)

# Calculate maximum overall remaining hospitalisations and deaths
max_ages_model = get_max_age(agegroups_model)
max_ages_model[is.na(max_ages_model)] = Inf
ovrl_res = res[,.(population=sum(population),
                  cum_v=sum(cum_prop_v*population,na.rm=T),
                  cum_hosp=sum(cum_hosp,na.rm=T),
                  cum_deaths=sum(cum_deaths,na.rm=T),
                  cum_exp=sum(cum_prop_exp*population)),by=.(sample,country,date)]
ovrl_res[,`:=`(cum_prop_v=cum_v/population,cum_prop_exp=cum_exp/population)]
ovrl_res[,`:=`(cum_inc_hosp=cum_hosp/population,cum_inc_deaths=cum_deaths/population)]

cols = setdiff(names(ovrl_res),c("sample","country","date"))
ovrl_med = ovrl_res[,lapply(.SD,median),.SDcols=cols,by=.(country,date)]
ovrl_q95l = ovrl_res[,lapply(.SD,function(x) quantile(x,probs=0.025)),.SDcols=cols,by=.(country,date)]
ovrl_q95u = ovrl_res[,lapply(.SD,function(x) quantile(x,probs=0.975)),.SDcols=cols,by=.(country,date)]
ovrl_rem_burden_dt = merge(ovrl_med,ovrl_q95l,by=c("country","date"),suffixes = c("","_q95l"))
ovrl_rem_burden_dt = merge(ovrl_rem_burden_dt,ovrl_q95u,by=c("country","date"),suffixes = c("","_q95u"))
ovrl_rem_burden_dt = merge(ovrl_rem_burden_dt,res[,.(prop_pop_60plus=sum(population[age_group_model %in% agegroups_model[max_ages_model>60]])/sum(population)),by=.(country)],by="country",all.x=T)

# Burden of remaining hospitalisations
print(ovrl_rem_burden_dt[c(which.min(cum_inc_hosp),which.max(cum_inc_hosp)),.(country,cum_inc_hosp=1e5*cum_inc_hosp,cum_inc_hosp_q95l=1e5*cum_inc_hosp_q95l,cum_inc_hosp_q95u=1e5*cum_inc_hosp_q95u)])
# Absolute number of remaining hospitalisations
print(ovrl_rem_burden_dt[c(which.min(cum_hosp),which.max(cum_hosp)),.(country,cum_hosp,cum_hosp_q95l,cum_hosp_q95u)])

# Burden of remaining deaths
print(ovrl_rem_burden_dt[c(which.min(cum_inc_deaths),which.max(cum_inc_deaths)),.(country,cum_inc_deaths=1e5*cum_inc_deaths,cum_inc_deaths_q95l=1e5*cum_inc_deaths_q95l,cum_inc_deaths_q95u=1e5*cum_inc_deaths_q95u)])
# Absolute number of remaining deaths
print(ovrl_rem_burden_dt[c(which.min(cum_deaths),which.max(cum_deaths)),.(country,cum_deaths,cum_deaths_q95l,cum_deaths_q95u)])

# Total absolute number of remaining hospitalisations and deaths across all countries
print(ovrl_rem_burden_dt[,.(cum_hosp=sum(cum_hosp),cum_deaths=sum(cum_deaths))])

# Cumulative proportion infected
print(ovrl_rem_burden_dt[c(which.min(cum_prop_exp),which.max(cum_prop_exp)),.(country,cum_prop_exp,cum_prop_exp_q95l,cum_prop_exp_q95u)])

# Plot maximum overall remaining burden:
# hospitalisations against vaccine coverage
p4 = ggplot(ovrl_rem_burden_dt,aes(x=cum_prop_v)) +
    geom_point(aes(y=cum_inc_hosp*1e5,size=cum_prop_exp),alpha=0.5,stroke=0) +
    geom_errorbar(aes(ymin=cum_inc_hosp_q95l*1e5,ymax=cum_inc_hosp_q95u*1e5),width = 0.01,alpha=0.5) +
    geom_text(aes(y=cum_inc_hosp*1e5,label=country),hjust=-0.2,vjust=0.4) +
    xlim(c(NA,1)) +
    labs(x="Proportion who have received at least one dose",
         y="Maximum remaining hospitalisations/100,000 population",
         size="Cumulative\nproportion\ninfected") +
    scale_y_log10()
# ggsave(paste0(dir_out,"rem_hosps_vs_prop_vax.png"),width = 8,height = 6.4)

# p = plot_grid(p1+theme(legend.position="none"),p2,p4,labels = c("A","B","C"),nrow=2,ncol=2)
# ggsave(paste0(dir_out,"vax_cov_and_rem_hosps.png"),plot = p,width = 12,height = 9.6)
p = plot_grid(p1,p2+theme(legend.position="none"),p4,labels = c("A","B","C"),nrow=2,ncol=2,rel_widths=c(1.1,1))
ggsave(paste0(dir_out,"vax_cov_and_rem_hosps.png"),plot=p,width = 13,height = 10.4)

# deaths against vaccine coverage
p5 = ggplot(ovrl_rem_burden_dt,aes(x=cum_prop_v)) +
    geom_point(aes(y=cum_inc_deaths*1e5,size=cum_prop_exp),alpha=0.5,stroke=0) +
    geom_errorbar(aes(ymin=cum_inc_deaths_q95l*1e5,ymax=cum_inc_deaths_q95u*1e5),width = 0.01,alpha=0.5) +
    geom_text(aes(y=cum_inc_deaths*1e5,label=country),hjust=-0.2,vjust=0.4) +
    xlim(c(NA,1)) +
    labs(x="Proportion who have received at least one dose",
         y="Maximum remaining deaths/100,000 population",
         size="Cumulative\nproportion\ninfected") +
    scale_y_log10()
# ggsave(paste0(dir_out,"rem_deaths_vs_prop_vax.png"),width = 8,height = 6.4)
ggsave(paste0(dir_out,"vax_cov_and_rem_deaths.png"),plot_grid(p3,p5,rel_widths=c(1,1),labels=c("A","B")),width = 14.4,height = 6)

# hospitalisations against proportion aged 60+
p6 = ggplot(ovrl_rem_burden_dt,aes(x=prop_pop_60plus)) +
    geom_point(aes(y=cum_inc_hosp*1e5,size=cum_prop_exp),alpha=0.5,stroke=0) +
    geom_errorbar(aes(ymin=cum_inc_hosp_q95l*1e5,ymax=cum_inc_hosp_q95u*1e5),width = 0.002,alpha=0.5) +
    geom_text(aes(y=cum_inc_hosp*1e5,label=country),hjust=-0.2,vjust=0.4) +
    xlim(c(0.23,0.31)) +
    labs(x="Proportion aged 60+ years",
         y="Maximum remaining hospitalisations/100,000 population",
         size="Cumulative\nproportion\ninfected") +
    scale_y_log10()
# ggsave(paste0(dir_out,"rem_hosps_vs_prop_60plus.png"),width = 8,height = 6.4)

# deaths against proportion aged 60+
p7 = ggplot(ovrl_rem_burden_dt,aes(x=prop_pop_60plus)) +
    geom_point(aes(y=cum_inc_deaths*1e5,size=cum_prop_exp),alpha=0.5,stroke=0) +
    geom_errorbar(aes(ymin=cum_inc_deaths_q95l*1e5,ymax=cum_inc_deaths_q95u*1e5),width = 0.002,alpha=0.5) +
    geom_text(aes(y=cum_inc_deaths*1e5,label=country),hjust=-0.2,vjust=0.4) +
    xlim(c(0.23,0.31)) +
    labs(x="Proportion aged 60+ years",
         y="Maximum remaining deaths/100,000 population",
         size="Cumulative\nproportion\ninfected") +
    scale_y_log10()
# ggsave(paste0(dir_out,"rem_deaths_vs_prop_60plus.png"),width = 8,height = 6.4)

p8 = plot_grid(p6+theme(legend.position="none"),
               p7+theme(legend.position="none"),
               labels=c("A","B"))
l = get_legend(p6 + theme(legend.box.margin = margin(0,0,0,12)))
ggsave(paste0(dir_out,"rem_hosps_and_deaths_vs_prop_60plus.png"),plot_grid(p8,l,rel_widths = c(3,.4)),width = 15,height = 6)

