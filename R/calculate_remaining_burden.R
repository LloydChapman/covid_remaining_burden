calculate_remaining_burden = function(fnm,agegroups_model,pop,ihr,ifr,frlty_idx,ve_params){
    # Load initial conditions calculation output
    prev_dt = readRDS(fnm)
    
    # Get minimum ages of age groups
    min_ages_model = get_min_age(agegroups_model)
    
    # Construct infection hospitalisation rate data table
    base_dt = merge(CJ(country=prev_dt[,unique(country)],age=0:100),pop[,.(country,age,population)],by=c("country","age"),all.x=T)
    # ihr = data.table(age = 0:85,ihr = exp(-7.37 + 0.068 * 0:85) / (1 + exp(-7.37 + 0.068 * 0:85)))
    # ihr[,age_group_model:=cut(age,c(min_ages_model,Inf),labels=agegroups_model,right=F)]
    # ihr = ihr[!is.na(age_group_model),.(ihr=mean(ihr)),by=.(age_group_model)] # exclude values where age group is missing as ages are not included in the backcalculation output
    ihr_dt = construct_ihr_data_table(ihr,base_dt,min_ages_model,agegroups_model)
    setnames(ihr_dt,"age_group","age_group_model")
    
    # Construct IFR data table
    ifr_dt = construct_ifr_data_table(ifr,base_dt,min_ages_model,agegroups_model)
    setnames(ifr_dt,"age_group","age_group_model")
    
    # Calculate remaining burden of hospitalisations and deaths
    res = calc_rem_burden(prev_dt,ihr_dt,ifr_dt,frlty_idx,ve_params)
    
    # Calculate 
    # cols = c("population","cum_prop_v",
    #          "cum_hosp","cum_deaths","cum_inc_hosp","cum_inc_deaths",
    #          "cum_hosp_u","cum_deaths_u","cum_inc_hosp_u","cum_inc_deaths_u",
    #          "cum_hosp_v","cum_deaths_v","cum_inc_hosp_v","cum_inc_deaths_v",
    #          "cum_hosp_i","cum_deaths_i","cum_inc_hosp_i","cum_inc_deaths_i",
    #          "cum_prop_exp")
    cols = c("population","S","V","R",
             "cum_hosp","cum_deaths","cum_inc_hosp","cum_inc_deaths",
             "cum_hosp_u","cum_deaths_u","cum_inc_hosp_u","cum_inc_deaths_u",
             "cum_hosp_v","cum_deaths_v","cum_inc_hosp_v","cum_inc_deaths_v",
             "cum_hosp_i","cum_deaths_i","cum_inc_hosp_i","cum_inc_deaths_i",
             "cum_prop_exp")
    med = res[,lapply(.SD,median),.SDcols=cols,by=.(country,age_group_model,date)]
    q95l = res[,lapply(.SD,function(x) quantile(x,probs=0.025)),.SDcols=cols,by=.(country,age_group_model,date)]
    q95u = res[,lapply(.SD,function(x) quantile(x,probs=0.975)),.SDcols=cols,by=.(country,age_group_model,date)]
    rem_burden_dt = merge(med,q95l,by=c("country","age_group_model","date"),suffixes = c("","_q95l"))
    rem_burden_dt = merge(rem_burden_dt,q95u,by=c("country","age_group_model","date"),suffixes = c("","_q95u"))
    
    # Calculate cumulative numbers vaccinated and unvaccinated, and corresponding population proportions
    # rem_burden_dt[,`:=`(cum_v=cum_prop_v*population,cum_u=(1-cum_prop_v)*population)]
    # rem_burden_dt[,`:=`(pop_prop_v=cum_v/sum(population),pop_prop_u=cum_u/sum(population)),by=.(country)]
    rem_burden_dt[,`:=`(pop_prop_u=S/sum(population),pop_prop_v=V/sum(population),pop_prop_i=R/sum(population)),by=.(country)]
    
    # Calculate maximum overall remaining hospitalisations and deaths
    max_ages_model = get_max_age(agegroups_model)
    max_ages_model[is.na(max_ages_model)] = Inf
    ovrl_res = res[,.(population=sum(population),
                      S=sum(S),
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
    
    # Remaining susceptible proportion
    print(ovrl_rem_burden_dt[c(which.min(S/population),which.max(S/population)),.(country,prevS=S/population,prevS_q95l=S_q95l/population,prevS_q95u=S_q95u/population)])
    
    return(list(rem_burden_dt=rem_burden_dt,ovrl_rem_burden_dt=ovrl_rem_burden_dt))
}
