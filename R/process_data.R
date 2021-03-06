process_data = function(source_deaths,country_iso_codes,agegroups,pop,Ab_delay1,Ab_delay2,ifr,ve_params,dir_fig,date_fitting){
    # DEATH DATA
    
    # Read in WHO death data
    # who_deaths = get_national_data(source="who")
    # setDT(who_deaths)
    # dir.create(paste0("./data/who_data/",date_fitting),recursive = T)
    # write.csv(who_deaths,paste0("./data/who_data/",date_fitting,"/who_data.csv"),row.names = F)
    who_deaths = fread(paste0("./data/who_data/",date_fitting,"/who_data.csv"))
    who_deaths[,date:=as.Date(date)]
    
    # Get England data
    # uk_data = get_regional_data("United Kingdom")
    # setDT(uk_data)
    # cols = c("region",intersect(names(who_deaths),names(uk_data)))
    # phe_deaths = uk_data[region=="England",..cols]
    # setnames(phe_deaths,"region","country")
    # dir.create(paste0("./data/phe_data/",date_fitting),recursive = T)
    # write.csv(phe_deaths,paste0("./data/phe_data/",date_fitting,"/phe_data.csv"),row.names = F)
    phe_deaths = fread(paste0("./data/phe_data/",date_fitting,"/phe_data.csv"))
    phe_deaths[,date:=as.Date(date)]
    
    # Bind to WHO data
    who_deaths = rbind(who_deaths,phe_deaths,fill=T)
    
    # Read in age-stratified death data
    deaths_raw = read_death_data(source_deaths)
    
    # Clean age-stratified death data
    deaths = clean_death_data(source_deaths,deaths_raw,who_deaths)
    
    # If using the COVerAGE data, then also read in and clean INED data to add 
    # countries missing from COVerAGE data
    if (source_deaths=="coverage"){
        deaths_raw_ined = read_death_data("ined")
        deaths_ined = clean_death_data("ined",deaths_raw_ined,who_deaths)
    }
    
    # Use INED data for the Netherlands and Romania as they are missing from 
    # COVerAGE data and for Denmark and Italy as it seems more consistent with WHO data
    countries_ined = c("Denmark","Italy","Netherlands","Romania")
    deaths = rbind(deaths[!(country %in% countries_ined),!"code"],
                   deaths_ined[country %in% countries_ined,
                               .(country,date,age_group,cum_deaths_both,cum_deaths_female,cum_deaths_male)])
    
    # # Plot to check
    # ggplot(deaths[,.(date,deaths_both=c(0,diff(cum_deaths_both))),by=.(country,age_group)],aes(x=date,y=deaths_both,group=age_group,color=age_group)) + 
    #     geom_line() + 
    #     facet_wrap(~country)
    
    # Read in LTC death data
    ltc_deaths = fread("./data/ltccovid_data/ltc_deaths.csv")
    
    # VACCINATION DATA
    
    # Read in ECDC vaccination data
    # # ecdc_vax = fread("./data/ecdc_data/ecdc_vaccination_data.csv")
    # ecdc_vax = get_data("https://opendata.ecdc.europa.eu/covid19/vaccine_tracker/csv","csv",
    #                     paste0("./data/ecdc_data/",date_fitting,"/"),"ecdc_vaccination_data.csv")
    ecdc_vax = fread(paste0("./data/ecdc_data/",date_fitting,"/ecdc_vaccination_data.csv"))
    
    # Clean ECDC vaccination data
    # out = clean_ecdc_vaccination_data(ecdc_vax,country_iso_codes)
    # vax = out$vax
    # num_type = out$num_type
    out = clean_ecdc_vaccination_data(ecdc_vax,country_iso_codes)
    vax = out$vax
    vax_type = out$vax_type
    agg_vax = out$agg_vax
    
    # # Checks on coverage in Germany and the Netherlands
    # ggplot(vax[country=="Netherlands" & dose==1],aes(x=date,y=count,color=age_group)) + geom_line() + facet_wrap(~type)
    # dates = seq.Date(vax[,min(date)],vax[,max(date)],by=1)
    # agegroups_vax = vax[,unique(age_group)]
    # agegroups_vax = agegroups_vax[agegroups_vax!="UNK"]
    # vax_dt = construct_vax_data_table(vax,dates,agegroups_vax,pop)
    # ggplot(vax_dt[,.(cum_prop=sum(cum_prop)),by=.(country,date,age_group,dose)][dose==1],aes(x=date,y=cum_prop,color=age_group)) +
    #     geom_line() +
    #     ylim(c(0,1)) +
    #     facet_wrap(~country)
    # ecdc_vax[ReportingCountry=="NL" & YearWeekISO<"2021-W41",.(first=sum(FirstDose),second=sum(SecondDose))]
    # vax[country=="Netherlands",.(count=sum(count)),by=.(dose)]
    # vax_dt[country=="Germany" & date==max(date),.(cum_prop=sum(cum_prop)),by=.(date,age_group,dose)]
    # vax_dt[country=="Germany" & date==max(date) & (age_group %in% c("18-24","25-49","50-59")) & dose==1,.(cum_prop=2*sum(cum_prop*population)/sum(population))]
    # vax_dt[country=="Germany" & date==max(date) & (age_group %in% c("60-69","70-79","80+")) & dose==1,.(cum_prop=2*sum(cum_prop*population)/sum(population))]
    
    # # Read in processed Public Health England (PHE) vaccination data
    # vaccPHE = readRDS("./data/phe_data/vax-covidm20211002121137.rds")
    # 
    # # List of age groups in covidm
    # agegroups_model = c("0-4","5-9","10-14","15-19","20-24","25-29","30-34","35-39",
    #                     "40-44","45-49","50-54","55-59","60-64","65-69","70-74","75+")
    # agegroups_model = factor(agegroups_model,levels = agegroups_model)
    # min_ages_model = get_min_age(agegroups_model)
    # 
    # # Process to same format as cleaned ECDC data
    # vaxENG = process_phe_vaccination_data(vaccPHE,agegroups_model)
    
    # Read in UK government vaccination data for England
    # vaccENG = get_data("https://api.coronavirus.data.gov.uk/v2/data?areaType=nation&areaCode=E92000001&metric=vaccinationsAgeDemographics&format=csv",
    #                    "csv","./data/gov_uk_data/",paste0("nation_",date_fitting-1,".csv")) # date of download file is previous day
    vaccENG = fread(paste0("./data/gov_uk_data/nation_",date_fitting-1,".csv"))
    
    # Process to same format as cleaned ECDC vaccination data
    vaxENG = process_gov_uk_vaccination_data(vaccENG,vax_type)
    
    # Get RKI vaccination data
    # rki_vax = get_data(paste0("https://raw.githubusercontent.com/robert-koch-institut/COVID-19-Impfungen_in_Deutschland/master/Archiv/",
    #                           date_fitting-1,"_Deutschland_Landkreise_COVID-19-Impfungen.csv"),"csv",paste0("./data/rki_data/",date_fitting,"/"),"rki_vaccination_data.csv")
    rki_vax = fread(paste0("./data/rki_data/",date_fitting,"/rki_vaccination_data.csv"))
        
    # Process to same format as cleaned ECDC vaccination data
    vaxDEU = process_rki_vaccination_data(rki_vax,ecdc_vax,pop)
    
    # # Plot
    # ggplot(vaxENG[,.(type,prop=count/sum(count,na.rm=T)),by=.(date,age_group,dose)][type=="vb"],
    #        aes(x=date,y=prop,group=factor(dose),color=factor(dose))) +
    #     geom_line() + facet_wrap(~age_group)
    # ggplot(vaxENG[type=="vb" & dose==1],aes(x=date,y=count,color=age_group)) + 
    #     geom_line()
    
    # VARIANT DATA

    # Read in ECDC variant data
    # # ecdc_vrnt_data = fread("./data/ecdc_data/ecdc_variant_data.csv")
    # ecdc_vrnt_data = get_data("https://opendata.ecdc.europa.eu/covid19/virusvariant/csv","csv",
    #                           paste0("./data/ecdc_data/",date_fitting,"/"),"ecdc_variant_data.csv")
    ecdc_vrnt_data = fread(paste0("./data/ecdc_data/",date_fitting,"/ecdc_variant_data.csv"))
    
    # Convert ISO weeks to dates
    ecdc_vrnt_data[,date:=as.IDate(ISOweek2date(paste0(sub("-","-W",year_week),"-7")))]
    
    # # Plot data for different sources
    # plot_variant_data(ecdc_vrnt_data,"GISAID")
    # plot_variant_data(ecdc_vrnt_data,"TESSy")
    
    # Process ECDC variant data
    vrnt_data = process_variant_data(ecdc_vrnt_data)
    
    # Read in COG England data
    # # cog_vrnt_data = fread("./data/cog_data/lineages_by_ltla_and_week.tsv")
    # cog_vrnt_data = get_data("https://covid-surveillance-data.cog.sanger.ac.uk/download/lineages_by_ltla_and_week.tsv",
    #                          "tsv",paste0("./data/cog_data/",date_fitting,"/"),"lineages_by_ltla_and_week.tsv")
    cog_vrnt_data = fread(paste0("./data/cog_data/",date_fitting,"/lineages_by_ltla_and_week.tsv"))
    
    # Convert week end date from Date to IDate
    cog_vrnt_data[,WeekEndDate:=as.IDate(WeekEndDate)]
    
    # Process COG data
    vrnt_dataENG = process_cog_variant_data(cog_vrnt_data)
    
    # Bind to ECDC data
    vrnt_data = rbind(vrnt_data,vrnt_dataENG,fill=T)
    
    # # Plot non-Alpha-Delta/Alpha/Delta proportions over time
    # ggplot(vrnt_data,aes(x=date,y=prop_vrnt,group=vrnt,color=vrnt)) +
    #     geom_line() +
    #     facet_wrap(~country)
    
    # Estimate variant proportions by multinomial regression
    res = estimate_variant_proportions(vrnt_data,deaths[,max(date)])
    
    # # Plot estimates for both models
    # plot_variant_proportions(res$vrnt_prop1,vrnt_data)
    # plot_variant_proportions(res$vrnt_prop2,vrnt_data)
    
    # Pick model 2 (model with a natural cubic spline function of sample date) as it 
    # has lower AIC
    vrnt_prop = res$vrnt_prop2
    
    # Construct data tables of deaths, vaccinations and IFR for deconvolution age
    # groups
    cols = c("cum_deaths_male","cum_deaths_female","cum_deaths_both")
    dt = construct_data_table(agegroups,deaths,pop,cols,ltc_deaths,vax,Ab_delay1,Ab_delay2,vaxENG,vaxDEU,ifr,vrnt_prop,ve_params)
    
    # # Plot data to check
    # # by age
    # ggplot(dt,aes(x=date,y=deaths_i_both,group=age_group,color=age_group)) +
    #     geom_line() +
    #     facet_wrap(~country)
    # 
    # # overall against WHO data
    # ggplot() +
    #     geom_line(aes(x = date, y = deaths_new),who_deaths[country %in% dt[,unique(country)]]) +
    #     geom_line(aes(x=date,y=deaths_i_both),dt[,.(deaths_i_both=sum(deaths_i_both)),by=.(country,date)],size=0.2,color="red") +
    #     facet_wrap(~country)
    
    # # overall INED against COVerAGE data
    # ggplot() +
    #     geom_line(aes(x=date,y=deaths_i_both),dt_ined[,.(deaths_i_both=sum(deaths_i_both)),by=.(country,date)]) +
    #     geom_line(aes(x=date,y=deaths_i_both),dt[,.(deaths_i_both=sum(deaths_i_both)),by=.(country,date)],size=0.2,color="red") +
    #     facet_wrap(~country)
    
    # Drop data for Iceland as there are too few deaths for deconvolution and for 
    # Ireland, Latvia and Romania as it is very incomplete
    dt = dt[!(country %in% c("Iceland","Ireland","Latvia"))]
    
    # Plot deaths
    plot_deaths(dt)
    ggsave(paste0(dir_fig,"deaths_by_age_",source_deaths,".png"),width = 6,height = 4.8)
    
    # Plot vaccine coverage
    # country-level
    start_date = as.Date("2020-12-01")
    line_labels = T
    p1 = plot_ovrl_vax_cov(dt[date>=start_date],line_labels)
    # ggsave(paste0(dir_fig,"vax_cov.png"),width = 5,height = 4)
    
    # by country and age
    ps1 = plot_vax_cov(dt[date>=start_date])
    # ggsave(paste0(dir_fig,"vax_cov_by_age.png"),width = 10,height = 8)
    
    # # Calculate and plot IFR with different assumptions:
    # # with scaling for Alpha and Delta severity
    # dt[,ifr_t:=((1-(1-ei)*cum_prop_v/(1-ei*cum_prop_v))+cum_prop_v/(1-ei*cum_prop_v)*(1-ei)*(1-ed)*(1-em))*(prop_vrnt+prop_vrnt2*1.5+prop_vrnt3*1.5*1.8)*ifr]
    # plot_ifr(dt[date>=start_date])
    # ggsave(paste0(dir_fig,"ifr_over_time_sev_scld.png"),width = 10,height = 8)
    # plot_avg_ifr(dt[date>=start_date])
    # ggsave(paste0(dir_fig,"avg_ifr_over_time_sev_scld.png"),width = 5,height = 4)
    # 
    # # without scaling for Alpha and Delta severity
    # dt[,ifr_t:=((1-(1-ei)*cum_prop_v/(1-ei*cum_prop_v))+cum_prop_v/(1-ei*cum_prop_v)*(1-ei)*(1-ed)*(1-em))*ifr]
    ps2 = plot_ifr(dt[date>=start_date])
    # ggsave(paste0(dir_fig,"ifr_over_time_not_sev_scld.png"),width = 10,height = 8)
    p2 = plot_avg_ifr(dt[date>=start_date],line_labels)
    # ggsave(paste0(dir_fig,"avg_ifr_over_time_not_sev_scld.png"),width = 5,height = 4)
    if (line_labels){
        p = p1 + p2 + plot_annotation(tag_levels = "A")
        ggsave(paste0(dir_fig,"vax_cov_and_IFR.png"),p,width = 10,height = 4.75,bg = "white")
        ggsave(paste0(dir_fig,"vax_cov_and_IFR.pdf"),p,width = 10,height = 4.75)
    } else {
        p = plot_grid(p1+theme(legend.position="none"),
                      p2+theme(legend.position="none"),
                      labels = c("A","B"))
        l = get_legend(p1 + theme(legend.box.margin = margin(0,0,0,12)))
        ggsave(paste0(dir_fig,"vax_cov_and_IFR.png"),plot_grid(p,l,rel_widths = c(3,.5)),width = 10,height = 4.75,bg = "white")
        ggsave(paste0(dir_fig,"vax_cov_and_IFR.pdf"),plot_grid(p,l,rel_widths = c(3,.5)),width = 10,height = 4.75)          
    }
    
    ps = plot_grid(ps1+theme(legend.position="none",axis.text.x=element_text(angle=90,hjust=1),strip.text=element_text(size=9)),
                   ps2+theme(legend.position="none",axis.text.x=element_text(angle=90,hjust=1),strip.text=element_text(size=9)),
                   labels = c("A","B"))
    ls = get_legend(ps1 + theme(legend.box.margin = margin(0,0,0,12)))
    ggsave(paste0(dir_fig,"vax_cov_and_IFR_by_age.png"),plot_grid(ps,ls,rel_widths = c(3,.4)),width = 10,height = 4.75,bg = "white")
    
    # Plot variant proportions
    plot_variant_proportions(vrnt_prop[country %in% dt[,unique(country)]],vrnt_data[country %in% dt[,unique(country)]])
    ggsave(paste0(dir_fig,"vrnt_prop_over_time.png"),width = 7.5,height = 6)
    
    return(list(dt=dt,vax=vax,vaxENG=vaxENG,vaxDEU=vaxDEU,vrnt_prop=vrnt_prop))
}