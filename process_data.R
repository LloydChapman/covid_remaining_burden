library(data.table)
library(qs)
library(ggplot2)
library(covidregionaldata)
library(ISOweek)
library(lubridate)
library(mgcv)
library(nnet)
library(splines)
library(effects)
library(cowplot)
library(stringr)

source("backcalculation_functions.R")

process_data = function(source_deaths,country_iso_codes,pop,Ab_delay1,Ab_delay2,vrnt_prop,ve_params,dir_out){
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
    ecdc_vax = get_data("https://opendata.ecdc.europa.eu/covid19/vaccine_tracker/csv","csv",
                        paste0("../ecdc_data/",date_fitting,"/"),"ecdc_vaccination_data.csv")
    
    # Clean ECDC vaccination data
    # out = clean_ecdc_vaccination_data(ecdc_vax,country_iso_codes)
    # vax = out$vax
    # num_type = out$num_type
    vax = clean_ecdc_vaccination_data(ecdc_vax,country_iso_codes)
    
    # Read in processed Public Health England (PHE) vaccination data
    vaccPHE = readRDS("../phe_data/vax-covidm20211002121137.rds")
    
    # List of age groups in covidm
    agegroups_model = c("0-4","5-9","10-14","15-19","20-24","25-29","30-34","35-39",
                        "40-44","45-49","50-54","55-59","60-64","65-69","70-74","75+")
    agegroups_model = factor(agegroups_model,levels = agegroups_model)
    min_ages_model = get_min_age(agegroups_model)
    
    # Process to same format as cleaned ECDC data
    vaxENG = process_phe_vaccination_data(vaccPHE,agegroups_model)
    
    # IFR
    
    # Read in ensemble IFR estimate from O'Driscoll et al Nature 2020
    ifr = fread(datapath("IFR_by_age_ODriscoll.csv"))
    
    # VARIANT DATA

    # Read in ECDC variant data
    # ecdc_vrnt_data = fread("../ecdc_data/ecdc_variant_data.csv")
    ecdc_vrnt_data = get_data("https://opendata.ecdc.europa.eu/covid19/virusvariant/csv","csv",
                              paste0("../ecdc_data/",date_fitting,"/"),"ecdc_variant_data.csv")
    
    # Convert ISO weeks to dates
    ecdc_vrnt_data[,date:=as.IDate(ISOweek2date(paste0(sub("-","-W",year_week),"-7")))]
    
    # Plot data for different sources
    plot_variant_data(ecdc_vrnt_data,"GISAID")
    plot_variant_data(ecdc_vrnt_data,"TESSy")
    
    # Process ECDC variant data
    vrnt_data = process_variant_data(ecdc_vrnt_data)
    
    # Read in COG England data
    # cog_vrnt_data = fread("../cog_data/lineages_by_ltla_and_week.tsv")
    cog_vrnt_data = get_data("https://covid-surveillance-data.cog.sanger.ac.uk/download/lineages_by_ltla_and_week.tsv",
                             "tsv",paste0("../cog_data/",date_fitting,"/"),"lineages_by_ltla_and_week.tsv")
    
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
    
    return(list(dt=dt,vrnt_prop=vrnt_prop))
}