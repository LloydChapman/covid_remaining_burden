library(data.table)
library(qs)
library(ggplot2)
library(ISOweek)
library(doParallel)

source("./backcalculation_functions.R")

# Register parallel backend
# registerDoParallel(cores = detectCores()-1)
registerDoParallel(cores = 4)

derive_initial_conditions = function(fnm,agegroups_model){
    # Load backcalculation output
    load(fnm)
    
    # # Source script with delay functions
    # cm_path = "./covidm_for_fitting/"
    # source(paste0(cm_path,"/R/shared/cmS_misc.R"))
    
    # Set plot theme
    theme_set(cowplot::theme_cowplot(font_size = 10) + theme(strip.background = element_blank()))
    
    # Get minimum ages of age groups
    min_ages_model = get_min_age(agegroups_model)
    
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
    # vax = clean_ecdc_vaccination_data(ecdc_vax,country_iso_codes)
    # vaxENG = process_phe_vaccination_data(vaccPHE,factor(params$pop[[1]]$group_names,levels=params$pop[[1]]$group_names))
    
    dates1 = seq.Date(backcalc_dt[,min(date)]-ceiling(max(Ab_delay1,Ab_delay2)/7)*7,backcalc_dt[,max(date)],by=1)
    vax_dt = construct_vax_data_table(vax,dates1,agegroups_model,pop)
    vaxENG_dt = construct_vax_data_table(vaxENG,dates1,agegroups_model,pop)
    vax_dt = rbind(vax_dt,vaxENG_dt)
    
    # Add lags for Ab development to date
    vax_dt[dose==1,date:=date+Ab_delay1]
    vax_dt[dose==2,date:=date+Ab_delay2]
    
    # Cast vaccination data to wide format
    vax_dt_wide = dcast(vax_dt,country + date + age_group + population ~ type + dose,value.var = c("count","prop","cum_prop"))
    
    # Restrict to dates for which there is data for both first and second doses given Ab development delays
    vax_dt_wide = vax_dt_wide[between(date,backcalc_dt[,min(date)],backcalc_dt[,max(date)])]
    
    setnames(vax_dt_wide,c("age_group","count_va_1","count_vb_1","count_va_2","count_vb_2"),c("age_group_model","nS_Va1","nS_Vb1","nVa1_Va2","nVb1_Vb2"))
    
    
    #
    # DISAGGREGATION OF INFECTIONS AND INITIAL CONDITIONS CALCULATION
    #
    
    
    # Make data table for incidence
    base_dt = CJ(country=backcalc_dt[,unique(country)],age=0:100,date=backcalc_dt[,unique(date)])
    base_dt = merge(base_dt,pop[,.(country,age,population)],by=c("country","age"),all.x=T)
    base_dt[,age_group_model:=cut(age,c(min_ages_model,Inf),labels=agegroups_model,right=F)]
    base_dt = base_dt[,.(population=sum(population)),by=.(country,age_group_model,date)]
    
    # Add deconvolution age groups to data table
    base_dt[,age_group:=cut(get_min_age(age_group_model),c(min_ages,Inf),labels=agegroups,right=F)]
    
    # prev_list = vector("list",nrow(exposures_samps[[1]]))
    # TODO - split this into non-LTC and LTC exposures (as they have different IFRs)
    tstart = Sys.time()
    prev_list = foreach(i=1:100) %dopar% {#nrow(exposures_samps[[1]])) %dopar% {
        exposures_vec = unlist(lapply(exposures_samps,"[",i,))
        exposures_ltc_vec = unlist(lapply(exposures_samps_ltc,"[",i,))
        
        # Merge with deconvolution output
        # N.B. duplicates exposures within same age_group
        backcalc_dt[,exposures:=exposures_vec]
        backcalc_dt[age_group %in% agegroups[max_ages>60],exposures_ltc:=exposures_ltc_vec]
        inc_dt = merge(base_dt,backcalc_dt[,.(country,age_group,date,exposures,exposures_ltc)],by=c("country","age_group","date"))
        
        # Split exposures by population fraction
        inc_dt[,`:=`(exposures=exposures*population/sum(population),
                     exposures_ltc=exposures_ltc*population/sum(population)),
               by=.(country,age_group,date)]
        
        # # Plot
        # ggplot(inc_dt,aes(x=date,y=exposures,group=age_group_model,color=age_group_model)) +
        #     geom_line() +
        #     facet_wrap(~country)    
        
        # Calculate initial conditions
        # prev_list[[i]] = calc_init_condns(inc_dt,vax_dt_wide,agegroups_model,covy_dt,vrnt_prop,ve_params,dE,dIp,dIs,dIa)
        calc_init_condns(inc_dt,vax_dt_wide,agegroups_model,covy_dt,vrnt_prop,ve_params,dE,dIp,dIs,dIa)
    }
    tend = Sys.time()
    print(tend-tstart)
    
    # Bind list into data table
    prev_dt = rbindlist(prev_list,idcol="sample")
    
    # # Remove backcalculation output
    # rm(backcalc_dt,backcalc_dt_non_ltc,backcalc_dt_ltc,
    #    backcalc_samps,backcalc_samps_non_ltc,backcalc_samps_ltc,
    #    exposures_samps,exposures_samps_non_ltc,exposures_samps_ltc,
    #    infections_samps,infections_samps_non_ltc,infections_samps_ltc)
    # 
    # # Remove functions
    # rm(list=lsf.str())
    # 
    # # Save
    # save(list=ls(all.names=T),file=paste0(dir_out,"init_cond_output.RData"),envir=environment())
    return(prev_dt)
}


