# Functions for getting min and max age from age group string
get_min_age = function(x) as.numeric(sub("-.*","",sub("\\+|<","-",x)))
get_max_age = function(x) as.numeric(sub(".*-","",sub("\\+|<","-",x)))

split_deaths = function(deaths,cntry,start_date,end_date,agg_age_group,cols,dataset="ined"){
    deaths_cntry = deaths[country==cntry & date>=start_date & date<=end_date & age_group=="Total",.SD,.SDcols = c("country","date","age_group",cols)]
    if (dataset=="ined"){
        cum_deaths_cntry = deaths[country==cntry & date==min(date[country==cntry & date>end_date]) & age_group=="Total",.SD,.SDcols = c("country","date","age_group",cols)]    
    } else {
        cum_deaths_cntry = who_deaths[country==cntry & date==min(date[country==cntry & date>end_date],na.rm=T),.(country,date,age_group="Total",cum_deaths_both=deaths_total,cum_deaths_male=NA,cum_deaths_female=NA)]
    }
    
    deaths_cntry[,`:=`(prop_cum_deaths_male=cum_deaths_male/cum_deaths_cntry[,cum_deaths_male],
                       prop_cum_deaths_female=cum_deaths_female/cum_deaths_cntry[,cum_deaths_female],
                       # prop_cum_deaths_unknown=cum_deaths_unknown/cum_deaths_cntry[,cum_deaths_unknown],
                       prop_cum_deaths_both=cum_deaths_both/cum_deaths_cntry[,cum_deaths_both])]
    
    agegroups = deaths[country==cntry & date==min(date[country==cntry & date>end_date]) & !grepl("Total",age_group),age_group]
    agegroups = agegroups[get_min_age(agegroups)<get_max_age(agg_age_group)]
    deaths_cntry_dt = CJ(country=cntry,age_group=agegroups,date=deaths_cntry[,unique(date)])
    # deaths_cntry_dt = merge(deaths_cntry_dt,deaths_cntry[,.(date,prop_cum_deaths_male,prop_cum_deaths_female,prop_cum_deaths_unknown,prop_cum_deaths_both)],by="date",all.x=T)
    deaths_cntry_dt = merge(deaths_cntry_dt,deaths_cntry[,.(date,prop_cum_deaths_male,prop_cum_deaths_female,prop_cum_deaths_both)],by="date",all.x=T)
    setnafill(deaths_cntry_dt,fill=0,cols=paste0("prop_",cols))
    
    # deaths_cntry_dt = merge(deaths_cntry_dt,deaths[country==cntry & date==min(date[country==cntry & date>end_date]) & age_group %in% agegroups,.SD,.SDcols=c("region","country_no","age_group",cols)],by=c("age_group"))
    next_cum_deaths_cntry = deaths[country==cntry & (age_group %in% agegroups) & date==min(date[country==cntry & date>end_date]),.SD,.SDcols=c("age_group",cols)]
    deaths_cntry_dt = merge(deaths_cntry_dt,next_cum_deaths_cntry,by="age_group")
    deaths_cntry_dt[,`:=`(cum_deaths_male=prop_cum_deaths_male*cum_deaths_male,
                      cum_deaths_female=prop_cum_deaths_female*cum_deaths_female,
                      # cum_deaths_unknown=prop_cum_deaths_unknown*cum_deaths_unknown,
                      cum_deaths_both=prop_cum_deaths_both*cum_deaths_both)]
    deaths_cntry_dt[,`:=`(prop_cum_deaths_male=NULL,
                      prop_cum_deaths_female=NULL,
                      # prop_cum_deaths_unknown=NULL,
                      prop_cum_deaths_both=NULL)] 
    
    deaths1 = copy(deaths)
    if (is.infinite(agg_age_group)){
        deaths1 = deaths1[!(country==cntry & date>=start_date & date<=end_date)]
    } else {
        deaths1 = deaths1[!(country==cntry & date>=start_date & date<=end_date & age_group==agg_age_group)]
    }
    deaths1 = rbind(deaths1,deaths_cntry_dt,fill=T)
    
    return(deaths1)
}

clean_ined_death_data = function(ined_deaths,cols){
    deaths = copy(ined_deaths) 
        
    # # Number of datasets per country
    # deaths[,length(unique(excelsheet)),by=.(country)]
    # # Countries with more than 1 dataset
    # deaths[,length(unique(excelsheet)),by=.(country)][V1>1,country]
    # 
    # # See what the unique age groups are
    # deaths[,unique(age_group)]
    # # by country
    # deaths[,table(country,age_group)]
    
    # Change names
    setnames(deaths,sub("cum_deaths","cum_death",cols),cols)
    
    # Convert dates
    deaths[,date:=as.Date(death_reference_date,format="%d.%m.%Y")]
    
    # Calculate totals with and without age group by country
    print(dcast(deaths[age_group %in% c("Total known","Total unknown"),max(cum_deaths_both),by=.(country,excelsheet,age_group)],country+excelsheet~age_group))
    # Very small proportion with missing age so ignore FOR NOW
    
    # Drop data for Canada as it doesn't cover first wave of pandemic
    deaths = deaths[country!="Canada"]
    
    # Exclude dataset with lower cumulative death total for countries with 2 datasets
    # or that doesn't cover up to the present in the case of France, or the start of 
    # the pandemic in the case of Canada.
    # TODO - correct data for France to account for only including hospital deaths
    # find another source (Google?)
    deaths = deaths[!(#(country=="Canada" & grepl("GC",excelsheet)) |
        (country=="France" & grepl("Cepi",excelsheet)) |
            (country=="Sweden" & grepl("NBHW",excelsheet)) |
            (country=="England & Wales" & grepl("NHS",excelsheet)))]
    
    # Remove what appear to be erroneous entries for Switzerland, USA, France, Italy and Portugal
    deaths = deaths[!((country=="Switzerland" & date %in% as.Date(c("2020-11-08","2020-11-14"))) |
                          (country=="United States" & date %in% as.Date(c("2021-05-02","2021-05-09"))) |
                          (country=="France" & date %in% as.Date(c("2020-12-03","2020-12-07","2020-12-08","2021-06-04"))) |
                          (country=="Italy" & between(date,as.Date("2020-08-11"),as.Date("2020-09-07"))) |
                          (country=="Portugal" & date %in% c(as.Date("2021-01-24"),as.Date("2021-01-25"))))]
    
    # Split deaths in aggregated age groups for Denmark by scaling the 
    # cumulative number of deaths up to the first date with disaggregated age 
    # groups for each age group by the distribution of the overall number of 
    # deaths in the period with aggregated age groups
    deaths = split_deaths(deaths,"Denmark",as.Date("2020-04-06"),as.Date("2020-08-10"),"<60",cols)
    deaths = split_deaths(deaths,"Denmark",as.Date("2020-04-03"),as.Date("2020-04-05"),"<70",cols)
    
    # Split early deaths in aggregated age groups for Germany
    deaths[country=="Germany" & age_group=="0-59",age_group:="<60"]
    deaths = split_deaths(deaths,"Germany",as.Date("2020-03-30"),as.Date("2020-04-24"),"<60",cols)
    
    # Aggregate 90-99 and 100+ age group data for Germany so for loop below works
    deaths[country=="Germany" & age_group %in% c("90-99","100+"),age_group:="90+"]
    # TODO - make use of cols robust here, so that only those columns are selected
    # deaths = deaths[,lapply(.SD,sum),.SDcols=cols,by=setdiff(names(deaths),cols)]
    deaths = deaths[,lapply(.SD,sum),.SDcols=cols,by=.(region,country,country_no,excelsource,excelsheet,age_group,death_reference_date,death_reference_date_type,date)]
    
    # Exclude totals
    # TODO - update this to incorporate deaths missing ages
    deaths = deaths[!(age_group %in% c("Total","Total known","Total unknown"))]
    
    # Reorder
    setorder(deaths,country,date,age_group)
}

clean_COVerAGE_death_data = function(COVerAGE_deaths,cols){
    out10 = copy(COVerAGE_deaths)
    
    # Change variable names to lower case
    names(out10) = tolower(names(out10))
    
    # Change country name for US
    out10[country=="USA",country:="United States"]
    
    # Select country-level data
    out10 = out10[region=="All"]
    
    # Rename variables
    cols1 = paste0("cum_",cols)
    setnames(out10,cols,cols1)
    
    # Convert dates
    out10[,date:=as.Date(date,format = "%d.%m.%Y")]
    
    # Add age group
    out10[,age_group:=paste0(age,"-",age+ageint-1)]
    
    # Change sex coding to full words
    out10[,sex:=fcase(sex=="b","both",
                      sex=="m","male",
                      sex=="f","female")]
    
    # Remove erroneous entries
    # TODO - correct data for Ireland and Scotland
    out10 = out10[!(country=="Canada" & date==as.Date("2020-04-30"))] # values appear to be incorrect due to missing "TOT" value for 2020-04-30 in inputDB
    out10 = out10[!(country=="Denmark" & date==as.Date("2021-01-05"))] # values appear to be incorrect (lower than previous date)
    out10 = out10[!(country=="France" & (between(date,as.Date("2020-12-06"),as.Date("2020-12-11"))) | (date %in% as.Date(c("2020-12-03","2021-06-04","2021-06-21"))))] # corrections to published records lead to negative new deaths
    out10 = out10[!(country=="Greece" & (between(date,as.Date("2021-02-03"),as.Date("2021-02-05")) | (date %in% as.Date(c("2021-03-01")))))]
    # out10 = out10[!(country=="Ireland" & )]
    out10 = out10[!(country=="Italy" & (between(date,as.Date("2020-08-11"),as.Date("2020-09-01")) | date %in% as.Date(c("2021-07-15","2021-07-16"))))]
    out10 = out10[country=="Italy" & date==as.Date("2020-10-20") & age_group=="60-69" & sex=="both",cum_deaths:=3671] # typo in total in inputDB 2671 should be 3671
    out10 = out10[country=="Italy" & date==as.Date("2020-10-20") & age_group=="60-69" & sex=="male",cum_deaths:=2789] # typo in total in inputDB leads to error in cleaned data
    out10 = out10[country=="Italy" & date==as.Date("2020-10-20") & age_group=="60-69" & sex=="female",cum_deaths:=882] # typo in total in inputDB leads to error in cleaned data
    out10 = out10[!(country=="Portugal" & date %in% as.Date(c("2020-10-04","2021-04-24")))] # incorrect all zero or partial zero entries
    out10 = out10[!(country=="Switzerland" & date==as.Date("2021-05-12"))] # values for 70+ seem to be aggregated
    # out10 = out10[!(country=="Scotland" & )]
    out10 = out10[!(country=="Ukraine" & date %in% as.Date(c("2021-05-12","2021-05-13","2021-05-18")))] # missing values for 80+ age groups & lower value for 70-80 age group than previous date
    
    # Calculate total deaths for merging total death counts from WHO time series
    # N.B. Needs some modification to make it work for multiple variables (cases, deaths, tests)
    notallNA = apply(out10[,mget(cols1)],1,function(x) !all(is.na(x)))
    agg_out10 = out10[notallNA,lapply(.SD,sum),.SDcols=cols1,by=.(country,region,code,date,sex)]
    agg_out10[,age_group:="Total"]
    out10 = rbind(out10[notallNA],agg_out10,fill=T)
    
    # notallNA1 = apply(out10[,mget(cols1)],1,function(x) !all(is.na(x)))
    # min_dates = rbind(out10[notallNA1,.(date=min(date)),by=.(country)],out10[!(country %in% unique(country[notallNA1])),.(date=min(date)),by=.(country)])
    min_dates = out10[,.(date=min(date)),by=.(country)]
    countries = out10[,unique(country)]
    total_deaths_list = vector("list",length(countries))
    cols2 = paste0(sub("tests","tested",cols),"_total")
    cols3 = c("country","date",cols2)
    for (i in 1:length(countries)){
        cntry = countries[i]
        tmp = who_deaths[country==cntry & date<min_dates[country==cntry,date],..cols3]
        setnames(tmp,cols2,cols1)
        tmp[,`:=`(region="All",sex="both",age_group="Total")]
        total_deaths_list[[i]] = tmp
    }
    total_deaths = rbindlist(total_deaths_list)
    out10 = rbind(out10,total_deaths,fill=T)
        
    # Cast to wide format for compatibility with construct_data_table function
    if (length(cols1)==1){
        out10[,sex:=paste0(cols1,"_",sex)]
    }
    out10 = dcast(out10,country + code + date + age_group ~ sex, value.var = cols1)
    
    for (i in 1:length(countries)){
        cntry = countries[i]
        out10 = split_deaths(out10,cntry,who_deaths[country==cntry,min(date)],min_dates[country==cntry,date]-1,Inf,paste0(cols1,"_",c("both","male","female")),dataset="COVerAGE")
    }
    
    # Reorder
    setorder(out10,country,date,age_group)
    
    # Drop rows with overall cumulative deaths
    out10 = out10[age_group!="Total"]
    
    return(out10)
}

clean_ecdc_vaccination_data = function(ecdc_vax,country_iso_codes){
    vax = copy(ecdc_vax)
    
    # Remove non-country-level data
    vax = vax[Region==ReportingCountry,]
    
    # Change country code for Greece
    vax[ReportingCountry=="EL", ReportingCountry:="GR"]
    
    # Add country names
    vax[,country:=country_iso_codes[match(vax[,ReportingCountry],iso_code),country]]
    
    # Convert ISO week to date
    vax[,date:=ISOweek2date(paste0(YearWeekISO,"-1"))] # start of ISO week as data is reported by ISO week and counts get split across ISO week in construct_data_table
    # Lose last 3 weeks of data as incomplete at time of download (current week may
    # not be finished and delay in reporting of 1-2 weeks)
    max_date = vax[,max(date)]
    vax = vax[!(date %in% (max_date-7*(0:2)))]
    
    # Drop data for HCWs and LTCF residents as this is contained in age group data
    vax = vax[!(TargetGroup %in% c("HCW","LTCF"))]
    vax[,age_group:=sub("_","-",sub("Age|1_Age","",TargetGroup))]
    
    # Drop data disaggregated by age for those <18 for now and data aggregated for 
    # <60 and 60+
    # TODO - revisit this
    # print(vax[age_group=="UNK",sum(FirstDose+SecondDose+UnknownDose),by=.(ReportingCountry)])
    vax = vax[!(age_group %in% c("0-4","5-9","10-14","15-17","<60","60+","ALL"))]
    vax[age_group=="<18", age_group := "0-17"]
    
    # Lose unneeded variables
    setnames(vax,c("Vaccine","ReportingCountry"),c("vaccine","country_code2"))
    # TODO - check population denominators against UN estimates
    vax = vax[,.(country,country_code2,year_week_iso=YearWeekISO,date,age_group,vaccine,first=FirstDose,second=SecondDose)] #,population=Denominator
    
    # Treat Pfizer and Moderna as vaccine B, and all others (AstraZeneca, Janssen, 
    # Sputnik, Beijing CBNG) as vaccine A
    vax[,type:=fifelse(vaccine %in% c("COM","MOD"),"vb",fifelse(vaccine!="UNK","va","vu"))]
    
    # # Calculate proportions of each vaccine type
    # # TODO - calculate over time
    # num_type = dcast(vax[,.(first=sum(first)),by=.(country,type)],country~type)
    # # Split vaccines of unknown type according to average proportion of each type in
    # # other countries FOR NOW - can correct with OWID data
    # cols3 = c("va","vb","vu")
    # num_type[,(cols3):=lapply(.SD,as.numeric),.SDcols=cols3]
    # num_type[!is.na(vu),`:=`(va=fifelse(is.na(va),vu*num_type[is.na(vu),mean(va/(va+vb))],va+vu*num_type[is.na(vu),mean(va/(va+vb))]),
    #                          vb=fifelse(is.na(va),vu*num_type[is.na(vu),mean(vb/(va+vb))],vb+vu*num_type[is.na(vu),mean(vb/(va+vb))]))]
    # num_type[,vu:=NULL]
    # num_type[,`:=`(prop_va=va/(va+vb),prop_vb=vb/(va+vb))]
    # 
    # # FOR NOW - sum over vaccine type
    # vax = vax[,.(first = sum(first), second = sum(second)), by = .(country,country_code2,year_week_iso,date,age_group)] #,population
    
    # Sum by vaccine type
    vax = vax[,.(first = sum(first), second = sum(second)), by = setdiff(names(vax),c("vaccine","first","second"))] #,population
    
    # Split vaccines of unknown type according to average proportion of each type in
    # other countries (should really also be by age group) FOR NOW - can correct with OWID data
    vax_wide = dcast(vax,... ~ type,value.var = c("first","second"))
    cols3 = paste0(rep(c("first","second"),each=3),"_",c("va","vb","vu"))
    vax_wide[,(cols3):=lapply(.SD,as.numeric),.SDcols=cols3]
    vax_wide[!is.na(first_vu),`:=`(first_va=fifelse(is.na(first_va),first_vu*vax_wide[is.na(first_vu),mean(first_va/(first_va+first_vb),na.rm=T)],first_va+first_vu*vax_wide[is.na(first_vu),mean(first_va/(first_va+first_vb),na.rm=T)]),
                             first_vb=fifelse(is.na(first_vb),first_vu*vax_wide[is.na(first_vu),mean(first_vb/(first_va+first_vb),na.rm=T)],first_vb+first_vu*vax_wide[is.na(first_vu),mean(first_vb/(first_va+first_vb),na.rm=T)]))]
    vax_wide[!is.na(second_vu),`:=`(second_va=fifelse(is.na(second_va),second_vu*vax_wide[is.na(second_vu),mean(second_va/(second_va+second_vb),na.rm=T)],second_va+second_vu*vax_wide[is.na(second_vu),mean(second_va/(second_va+second_vb),na.rm=T)]),
                             second_vb=fifelse(is.na(second_vb),second_vu*vax_wide[is.na(second_vu),mean(second_vb/(second_va+second_vb),na.rm=T)],second_vb+second_vu*vax_wide[is.na(second_vu),mean(second_vb/(second_va+second_vb),na.rm=T)]))]
    vax_wide[,`:=`(first_vu=NULL,second_vu=NULL)]
    
    # Melt to long format
    vax = melt(vax_wide,measure.vars = patterns(first="first",second="second"),variable.name = "type")
    vax[,type:=fcase(type==1,"va",type==2,"vb")]
    
    # Remove rows with missing entries for both first and second doses
    vax = vax[!(is.na(first) & is.na(second))]
    
    # # Calculate proportions of each vaccine type over time
    # vax[,`:=`(prop_first=first/sum(first),prop_second=second/sum(second)),by=.(country,date,age_group)]
    
    # See how many doses have unknown age group
    print(vax[,sum(first+second)]) # 209715607
    print(vax[age_group!="UNK",sum(first+second)]) #209476181
    print(vax[,.(first=sum(first),second=sum(second)), by = .(age_group)]) # only about 240,000 out of about 210 million, so ignore FOR NOW
    
    # return(list(vax=vax,num_type=num_type))
    return(vax)
}

construct_vax_data_table = function(vax,dates,agegroups){
    # Get vaccination age groups from vax
    agegroups_vax = sort(vax[,unique(age_group)])
    agegroups_vax = agegroups_vax[agegroups_vax!="UNK"]
    min_ages_vax = get_min_age(agegroups_vax)
    min_ages_vax[is.na(min_ages_vax)] = 0
    
    # Get minimum ages of age groups
    min_ages = get_min_age(agegroups)
    
    # Create data table with all combinations of countries, dates, vaccines and ages
    # to store vaccine schedule
    vax_dt = CJ(country=vax[,unique(country)],date=dates,age=0:100,type=vax[,unique(type)])
    vax_dt = merge(vax_dt,pop[,.(country,age,population)],by=c("country","age"),all.x=T)
    vax_dt[,year_week_iso:=ISOweek(date)]
    vax_dt[,age_group:=cut(age,c(min_ages_vax,Inf),labels=agegroups_vax,right=F)]
    
    # Merge with vax data table
    # N.B. This duplicates numbers of doses for the same ISO week and age group, so 
    # we then divide by 7 to get the average doses per day and divide doses between 
    # age groups according to population fraction
    vax_dt = merge(vax_dt,vax[,!"date"],by=c("country","year_week_iso","age_group","type"),all.x=T)
    cols = c("first","second")
    vax_dt[,(cols):=lapply(.SD,as.numeric),.SDcols=cols]
    vax_dt[,(cols):=lapply(.SD,function(x) x/7),.SDcols=cols]
    vax_dt[,(cols):=lapply(.SD,function(x) x*population/sum(population)),.SDcols=cols,by=.(country,date,age_group,type)]
    print(vax_dt[,sum(first+second,na.rm=T)]) #209476181
    
    # Change age groups
    vax_dt[,age_group:=cut(age,c(min_ages,Inf),agegroups,right=F)]
    # Sum vaccinations in each age group by type
    vax_dt = vax_dt[,lapply(.SD,function(x) sum(x,na.rm=T)),.SDcols=c("population",cols),by=.(country,date,age_group,type)]
    # Calculate proportion fully vaccinated
    vax_dt[,prop:=second/population]
    # Calculate cumulative proportion fully vaccinated ensuring data is in chronological date order
    setorder(vax_dt,country,date,age_group,type)
    vax_dt[,cum_prop:=cumsum(prop),by=.(country,age_group,type)]
    
    return(vax_dt)
}

# construct_deaths_data_table = function(base_deaths_dt,agegroups,min_date){
#     base_deaths_dt1 = copy(base_deaths_dt)
#     
#     # Change age groups
#     min_ages = get_min_age(agegroups)
#     base_deaths_dt1[,age_group:=cut(age,c(min_ages,Inf),labels=agegroups,right=F)]
#     # Sum deaths in each age group
#     base_deaths_dt1 = base_deaths_dt1[,lapply(.SD,sum),.SDcols=c("population",cols),by=.(region,country,country_code,date,age_group)]
#     
#     dates = seq.Date(min_date,base_deaths_dt1[,max(date),by=.(country)][,min(V1)],by=1)
#     deaths_dt = CJ(country_code = base_deaths_dt1[,unique(country_code)], date = dates, age_group = agegroups)
#     deaths_dt = merge(deaths_dt,base_deaths_dt1[,!c("region","country","population")],by=c("country_code","date","age_group"),all.x=T)
#     
#     # Exclude country and variable combinations with all missing data
#     deaths_dt_long = melt(deaths_dt,id.vars = c("country_code","date","age_group"),variable.name="sex",value.name="cum_death")
#     deaths_dt_long[,sex:=sub("cum_death_","",sex)]
#     country_codes_notallNA = deaths_dt_long[,!all(is.na(cum_death)),by=.(country_code,sex)][V1==T,.(country_code,sex)]
#     deaths_dt_long = deaths_dt_long[country_codes_notallNA,on=c("country_code","sex")]
#     
#     # Fill in cumulative deaths for earliest date so there is a value to interpolate from
#     deaths_dt_long[date == min_date,cum_death:=0]
#     
#     # Interpolate cumulative deaths
#     deaths_dt_long[,cum_death_i := approx(date,cum_death,dates)$y,by=.(country_code,sex,age_group)]
#     # deaths_dt[,(paste0(cols,"_i")) := lapply(.SD,function(yi) approx(date,yi,dates)$y),.SDcols=cols,by=.(country_code,age_group)]
#     
#     # Calculate new daily deaths
#     # FOR NOW - as large negative counts have been removed, treat remaining small 
#     # negative counts from differencing as 0s
#     deaths_dt_long[,death_i := c(0,pmax(diff(cum_death_i),0)),by=.(country_code,sex,age_group)]
#     
#     # Cast to wide format
#     deaths_dt = dcast(deaths_dt_long,country_code + date + age_group ~ sex, value.var = c("cum_death","cum_death_i","death_i"))
#     # Add region, country and country code
#     deaths_dt = merge(deaths_dt,unique(deaths[,.(region,country,country_code)]),by="country_code")
#     
#     return(deaths_dt)
# }

construct_ifr_data_table = function(ifr,base_dt,min_ages,agegroups){
    ifr_dt = copy(base_dt)
    
    names(ifr) = tolower(names(ifr))
    min_ages_ifr = ifr[,get_min_age(age_group)]
    
    # Add IFR age groups
    ifr_dt[,age_group:=cut(age,c(min_ages_ifr,Inf),labels=ifr[,age_group],right=F)]
    
    # Merge with IFR data table
    ifr_dt = merge(ifr_dt,ifr[,.(age_group,ifr=median_perc_mean/100,ifr_lb=ci_95_lb_mean/100,ifr_ub=ci_95_ub_mean/100)],by="age_group")
    
    # Change age groups
    ifr_dt[,age_group := cut(age,c(min_ages,Inf),labels=agegroups,right=F)]
    
    # Calculate population-weighted average for each age group
    cols2 = c("ifr","ifr_lb","ifr_ub")
    ifr_dt = ifr_dt[,lapply(.SD,function(x) sum(x*population)/sum(population)),.SDcols=cols2,by=c("country","age_group")]    
    
}

# construct_data_table = function(agegroups,deaths,pop,cols,ltc_deaths,vax,num_type,ifr){
construct_data_table = function(agegroups,deaths,pop,cols,ltc_deaths,vax,Ab_delay,ifr){
    # Get minimum ages of age groups
    min_ages = get_min_age(agegroups)
    max_ages = get_max_age(agegroups)
    max_ages[is.na(max_ages)] = Inf
    
    # Make data table of unique combinations of country and age
    countries = deaths[,unique(country)]
    # country_codes = deaths[,unique(country_code)]
    base_dt = CJ(country = countries,age = 0:100)
    
    # Add population data
    base_dt = merge(base_dt,pop[,.(country,age,population)],by=c("country","age"),all.x=T)
    # base_dt = pop[country %in% countries,.(country,age,population)]
    
    # Make a copy
    base_deaths_dt = copy(base_dt)
    
    # Add age groups from deaths data
    for (i in 1:length(countries)){
        agegroups_cntry = sort(deaths[country==countries[i],unique(age_group)])
        min_ages_cntry = get_min_age(agegroups_cntry)
        min_ages_cntry[is.na(min_ages_cntry)] = 0
        agegroups_cntry = agegroups_cntry[order(min_ages_cntry)]
        
        base_deaths_dt[country==countries[i],age_group:=cut(age,c(min_ages_cntry,Inf),labels=agegroups_cntry,right=F)]
    }
    
    # Merge with deaths data table
    # N.B. This duplicates deaths for the same age group, so we divide deaths between 
    # ages according to population fraction
    # TODO - revisit this and divide deaths between age groups by relative IFR rather
    # than population fraction?
    base_deaths_dt = merge(base_deaths_dt,deaths,by=c("country","age_group"),allow.cartesian = T)
    base_deaths_dt[,(cols):=lapply(.SD,as.numeric),.SDcols=cols]
    base_deaths_dt[,(cols):=lapply(.SD,function(x) 
        if(sum(population,na.rm=T)!=0){
            x*population/sum(population)
        } else {
            x/.N
        }),
        .SDcols=cols,by=.(country,date,age_group)]
    
    # Check death totals match
    print(deaths[!(country %in% base_deaths_dt[is.na(population),country]),lapply(.SD,function(x) sum(x,na.rm=T)),.SDcols=cols])
    print(base_deaths_dt[!is.na(population),lapply(.SD,function(x) sum(x,na.rm=T)),.SDcols=cols])
    
    # Change age groups
    base_deaths_dt[,age_group:=cut(age,c(min_ages,Inf),labels=agegroups,right=F)]
    # Sum deaths in each age group
    base_deaths_dt = base_deaths_dt[,lapply(.SD,sum),.SDcols=c("population",cols),by=.(country,date,age_group)]
    
    # dates = seq.Date(min_date,base_deaths_dt[,max(date),by=.(country)][,min(V1)],by=1)
    # deaths_dt = CJ(country = base_deaths_dt[,unique(country)], date = dates, age_group = agegroups)
    # deaths_dt = merge(deaths_dt,base_deaths_dt[,!"population"],by=c("country","date","age_group"),all.x=T)
    deaths_dt_list = vector("list",length(countries))
    for (i in 1:length(countries)){
        # Make data table for each country with all dates and age groups in
        cntry = countries[i]
        date_range = base_deaths_dt[country==cntry,range(date),]
        # Start dates from 60 days before first death record to allow for maximum possible infection-to-death time in deconvolution below
        tmp = CJ(country = cntry, date = seq.Date(date_range[1]-60,date_range[2],by=1), age_group = agegroups)
        deaths_dt_list[[i]] = merge(tmp,base_deaths_dt[country==cntry,!"population"],by=c("country","date","age_group"),all.x=T)
    }
    deaths_dt = rbindlist(deaths_dt_list)

    # Exclude country and variable combinations with all missing data
    deaths_dt_long = melt(deaths_dt,id.vars = c("country","date","age_group"),variable.name="sex",value.name="cum_deaths")
    deaths_dt_long[,sex:=sub("cum_deaths_","",sex)]
    # countries_notallNA = deaths_dt_long[,!all(is.na(cum_deaths)),by=.(country,sex)][V1==T,.(country,sex)]
    # deaths_dt_long = deaths_dt_long[countries_notallNA,on=c("country","sex")]
    morethan1obs = deaths_dt_long[!is.na(cum_deaths),.(length(unique(date))),by=.(country,sex)][V1>1,.(country,sex)]
    deaths_dt_long = deaths_dt_long[morethan1obs,on=c("country","sex")]

    # Fill in cumulative deaths for earliest date so there is a value to interpolate from
    countries1 = deaths_dt_long[,unique(country)]
    for (i in 1:length(countries1)){
        deaths_dt_long[country==countries1[i] & date==min(date[country==countries1[i]]),cum_deaths:=0]    
    }

    # Interpolate cumulative deaths
    deaths_dt_long[,cum_deaths_i := approx(date,cum_deaths,date)$y,by=.(country,sex,age_group)]
    # deaths_dt[,(paste0(cols,"_i")) := lapply(.SD,function(yi) approx(date,yi,dates)$y),.SDcols=cols,by=.(country,age_group)]

    # Calculate new daily deaths
    # FOR NOW - as large negative counts have been removed, treat remaining small
    # negative counts from differencing as 0s
    # TODO - sort out whether diff should be padded with 0 (correct only if 
    # there have been no deaths by date of first record) or cum_deaths[1]
    # (correct if death data is complete) or NA
    deaths_dt_long[,deaths_i := c(0,pmax(diff(cum_deaths_i),0)),by=.(country,sex,age_group)]

    # Cast to wide format
    deaths_dt = dcast(deaths_dt_long,country + date + age_group ~ sex, value.var = c("cum_deaths","cum_deaths_i","deaths_i"))
    # # Add region, country and country code
    # deaths_dt = merge(deaths_dt,unique(deaths[,.(region,country,country_code)]),by="country")
    
    # Plot data
    # cumulative deaths
    print(ggplot(deaths_dt,aes(x=date,y=cum_deaths_i_both,group=age_group,color=age_group)) +
        geom_line() +
        scale_y_log10() +
        facet_wrap(~country))
    # new deaths
    print(ggplot(deaths_dt[country!="United States"],aes(x=date,y=deaths_i_both,group=age_group,color=age_group)) +
        geom_line() +
        facet_wrap(~country))
    
    # LTC death data
    names(ltc_deaths) = tolower(names(ltc_deaths))
    names(ltc_deaths)[3:8] = c("approach","total_deaths","total_ltc_res_deaths","total_ltc_deaths","perc_ltc_res_deaths","perc_ltc_deaths")
    cols5 = c("total_deaths","total_ltc_res_deaths","total_ltc_deaths")
    ltc_deaths[,(cols5):=lapply(.SD,function(x) as.numeric(sub(",","",x))),.SDcols = cols5]
    ltc_deaths[,`:=`(prop_ltc_res_deaths=total_ltc_res_deaths/total_deaths,prop_ltc_deaths=total_ltc_deaths/total_deaths)]
    
    # Use mean proportion of deaths in LTCs for countries without LTC death data
    ltc_deaths_mssng = data.table(country=setdiff(deaths_dt[,unique(country)],ltc_deaths[,country]),
                            prop_ltc_res_deaths=ltc_deaths[,mean(prop_ltc_res_deaths,na.rm=T)],
                            prop_ltc_deaths=ltc_deaths[,mean(prop_ltc_deaths,na.rm=T)])
    ltc_deaths = rbind(ltc_deaths,ltc_deaths_mssng,fill=T)
    
    # Merge with deaths data table
    deaths_dt = merge(deaths_dt,ltc_deaths[,.(country,prop_ltc_res_deaths,prop_ltc_deaths)],by="country")
    
    # Calculate LTC deaths among those 60+ by multiplying total deaths by LTC proportion
    # N.B. relies on having proportion of deaths in LTCs or LTC residents
    deaths_dt[age_group %in% agegroups[max_ages<60],deaths_i_both_ltc := 0] # assume no deaths in under-60s
    deaths_dt[age_group %in% agegroups[max_ages>60],deaths_i_both_ltc := deaths_i_both * fifelse(!is.na(prop_ltc_res_deaths),prop_ltc_res_deaths,prop_ltc_deaths)]
    deaths_dt[,deaths_i_both_non_ltc := deaths_i_both - deaths_i_both_ltc]
    
    # VACCINATION DATA
    
    # Get ISO week of earliest date with death record from deaths data
    min_iso_week = ISOweek(deaths_dt[,min(date)])
    # Get date of start of ISO week to ensure length of date vector is a multiple of 7
    min_date = ISOweek2date(paste0(min_iso_week,"-1")) 
    # Make daily date sequence from earliest death date to latest date for which 
    # vaccine data is available - length should be a multiple of 7
    dates1 = seq.Date(min_date-ceiling(Ab_delay/7)*7,vax[,max(date)+6],by=1) # -Ab_delay to account for Ab_delay-day delay to Ab development, and max(date) + 6 to get end of last ISO week
    
    # # Get vaccination age groups from vax
    # agegroups_vax = sort(vax[,unique(age_group)])
    # agegroups_vax = agegroups_vax[agegroups_vax!="UNK"]
    # min_ages_vax = get_min_age(agegroups_vax)
    # min_ages_vax[is.na(min_ages_vax)] = 0
    # max_ages_vax = get_max_age(agegroups_vax)
    # max_ages_vax[is.na(max_ages_vax)] = Inf
    # 
    # # Create data table with all combinations of countries, dates, vaccines and ages
    # # to store vaccine schedule
    # vax_dt = CJ(country=vax[,unique(country)],date=dates1,age=0:100,type=vax[,unique(type)])
    # vax_dt = merge(vax_dt,pop[,.(country,age,population)],by=c("country","age"),all.x=T)
    # vax_dt[,year_week_iso:=ISOweek(date)]
    # vax_dt[,age_group:=cut(age,c(min_ages_vax,Inf),labels=agegroups_vax,right=F)]
    # vax_dt[,country_code2:=country_iso_codes[match(vax_dt[,country],country),iso_code]]
    # 
    # # Merge with vax data table
    # # N.B. This duplicates numbers of doses for the same ISO week and age group, so 
    # # we then divide by 7 to get the average doses per day and divide doses between 
    # # age groups according to population fraction
    # vax_dt = merge(vax_dt,vax[,!"date"],by=c("country","country_code2","year_week_iso","age_group","type"),all.x=T)
    # cols1 = c("first","second")
    # vax_dt[,(cols1):=lapply(.SD,as.numeric),.SDcols=cols1]
    # vax_dt[,(cols1):=lapply(.SD,function(x) x/7),.SDcols=cols1]
    # vax_dt[,(cols1):=lapply(.SD,function(x) x*population/sum(population)),.SDcols=cols1,by=.(country,date,age_group,type)]
    # print(vax_dt[,sum(first+second,na.rm=T)]) #209476181
    # 
    # # Change age groups
    # vax_dt[,age_group:=cut(age,c(min_ages,Inf),agegroups,right=F)]
    # # Sum vaccinations in each age group by type
    # vax_dt = vax_dt[,lapply(.SD,function(x) sum(x,na.rm=T)),.SDcols=c("population",cols1),by=.(country,date,age_group,type)]
    # # Calculate proportion fully vaccinated
    # vax_dt[,prop:=second/population]
    # # Calculate cumulative proportion fully vaccinated ensuring data is in chronological date order
    # setorder(vax_dt,country,date,age_group,type)
    # vax_dt[,cum_prop:=cumsum(prop),by=.(country,age_group,type)]
    
    # Construct vaccination data table
    vax_dt = construct_vax_data_table(vax,dates1,agegroups)
    
    # Plot data
    # cumulative proportion fully vaccinated
    print(ggplot(vax_dt[date>=vax[,min(date)],.(cum_prop=sum(cum_prop)),by=.(country,age_group,date)],aes(x=date,y=cum_prop,group=age_group,color=age_group)) + 
        geom_line() +
        facet_wrap(~country))
    
    # proportion newly fully vaccinated
    print(ggplot(vax_dt[date>=vax[,min(date)],.(prop=sum(prop)),by=.(country,age_group,date)],aes(x=date,y=prop,group=age_group,color=age_group)) + 
        geom_line() +
        facet_wrap(~country))
    
    
    # Merge death and vaccination data
    # See which countries are in both datasets
    print(intersect(deaths_dt[,country],vax_dt[,country]))
    
    # Cast vaccination data to wide format
    vax_dt_wide = dcast(vax_dt,country + date + age_group + population ~ type,value.var = c("first","second","prop","cum_prop"))
    
    # Add Ab_delay days to date in vaccination data for development of Ab
    vax_dt_wide[,date_v:=date]
    vax_dt_wide[,date:=date+Ab_delay]
    # Make vector of dates that are covered by both the death data and the 
    # vaccination data
    # Use the last vaccination date as the cut-off for the vaccination data 
    # rather than the last Ab development date as we need to use the vaccination
    # data later for the numbers entering the vaccination compartments
    dates2 = as.Date(intersect(deaths_dt[,date],vax_dt_wide[,date_v]),origin="1970-01-01")
    deaths_vax_dt = merge(deaths_dt[date %in% dates2],vax_dt_wide[date %in% dates2],by=c("country","date","age_group"))
    
    # IFR ESTIMATES
    ifr_dt = copy(base_dt)
    
    # # Add IFR estimate from Levin et al Eur Jrnl Epi 2020
    # ifr_dt[,`:=`(ifr=10^(-3.27+0.0524*age)/100)]
    # # Calculate 95% CI for IFR assuming FOR NOW no covariance between intercept and 
    # # slope estimates
    # # TODO - estimate intercept and slope covariance from lower and upper bounds in 
    # # supplementary spreadsheet in Levin et al 2020
    # se_intcpt = 0.07*log(10)
    # se_slope = 0.0013*log(10)
    # ifr_dt[,se_ifr:=sqrt(se_intcpt^2+se_slope^2*age^2)*ifr]
    # ifr_dt[,`:=`(ifr_lb=pmax(0,ifr-qnorm(0.975)*se_ifr),ifr_ub=pmin(1,ifr+qnorm(0.975)*se_ifr))]
    # 
    # ifr_dt[,age_group:=cut(age,c(min_ages,Inf),labels=agegroups,right=F)]
    # cols2 = c("ifr","se_ifr","ifr_lb","ifr_ub")
    # ifr_dt=ifr_dt[,lapply(.SD,function(x) sum(x*population)/sum(population)),.SDcols=cols2,by=.(country,age_group)]
    
    names(ifr) = tolower(names(ifr))
    min_ages_ifr = get_min_age(ifr[,age_group])
    
    # Add IFR age groups
    ifr_dt[,age_group:=cut(age,c(min_ages_ifr,Inf),labels=ifr[,age_group],right=F)]
    
    # Merge with IFR data table
    ifr_dt = merge(ifr_dt,ifr[,.(age_group,ifr=median_perc_mean/100,ifr_lb=ci_95_lb_mean/100,ifr_ub=ci_95_ub_mean/100)],by="age_group")
    
    # Change age groups
    ifr_dt[,age_group := cut(age,c(min_ages,Inf),labels=agegroups,right=F)]
    
    # Calculate population-weighted average for each age group
    cols2 = c("ifr","ifr_lb","ifr_ub")
    ifr_dt = ifr_dt[,lapply(.SD,function(x) sum(x*population)/sum(population)),.SDcols=cols2,by=c("country","age_group")]    
    
    # Plot country IFRs
    ggplot(ifr_dt,aes(x=age_group,y=ifr,group=country,color=country)) + 
        geom_line() +
        # scale_y_log10() +
        theme(legend.position = "none")
        
    # Merge IFR with death and vaccination data
    dt = merge(deaths_vax_dt,ifr_dt,by=c("country","age_group"))
    # # N.B. Countries without vax data dropped here - check this!
    # dt = merge(dt,num_type[,.(country,prop_va,prop_vb)],by="country")
    
    # Vaccine efficacy parameters
    # # eiX_vYZ = efficacy of dose Z of vaccine Y against infection with strain X
    # ei_va2 = 0.68
    # ei2_va2 = 0.68
    # ei3_va2 = 0.6154
    # ei_vb2 = 0.85
    # ei2_vb2 = 0.85
    # ei3_vb2 = 0.7999
    # ed_vYZiX = efficacy of dose Z of vaccine Y against disease given infection with strain X
    ed_va2i = 0.3125
    ed_va2i2 = 0.3125
    ed_va2i3 = 0.2353094
    ed_vb2i = 0.2667
    ed_vb2i2 = 0.2667
    ed_vb2i3 = 0.187906
    # em_vYZdX = efficacy of dose Z of vaccine Y against death given disease from strain X
    em_va2d = 0.77
    em_va2d2 = 0.77
    em_va2d3 = 0.6766406
    em_vb2d = 0.55
    em_vb2d2 = 0.55
    em_vb2d3 = 0.52
    
    # # FOR NOW - use efficacy values for wild-type variant (N.B. same as for Delta)
    # # TODO - use proportion of each vaccine by age group and time and proportion of 
    # # each variant over time to calculate more detailed age- and time-dependent IFR
    # # dt[,`:=`(ei=prop_va*ei_va2+prop_vb*ei_vb2,ed=prop_va*ed_va2i+prop_vb*ed_vb2i,em=prop_va*em_va2d+prop_vb*em_vb2d)]
    # dt[,`:=`(ed=prop_va*ed_va2i+prop_vb*ed_vb2i,em=prop_va*em_va2d+prop_vb*em_vb2d)]
    
    # Calculate IFRs under vaccination with different vaccine types
    # dt[,ifr_v:=(1-ed)*(1-em)*ifr]
    dt[,`:=`(ifr_va=(1-ed_va2i)*(1-em_va2d)*ifr,ifr_vb=(1-ed_vb2i)*(1-em_vb2d)*ifr)]
    
    # Calculate IFR over time
    # dt[,ifr_t:=(1-cum_prop_v)*ifr+cum_prop_v*ifr_v]
    dt[,ifr_t := (1-cum_prop_va-cum_prop_vb)*ifr + cum_prop_va*ifr_va + cum_prop_vb*ifr_vb]
    
    # Plot population-weighted average IFR over time for all countries
    ggplot(dt[,.(ifr_t=sum(ifr_t*population)/sum(population)),by=.(country,date)],
           aes(x=date,y=ifr_t,group=country,color=country)) + geom_line()
    # ggsave("./output/avg_ifr_over_time.pdf",width = 5,height = 4)
    
    return(dt)
}

# Deconvolution function
deconv = function(dt,idd_pmf,method = "ride"){
    dt1 = copy(dt)
    
    countries = dt1[,unique(country)]
    
    dt_list = split(dt1[,.(country,age_group,date,deaths_i_both)],by=c("country","age_group"))
    
    mean_idd = sum((0:(length(idd_pmf)-1))*idd_pmf)
    
    if (method == "backproj"){ # BACK PROJECTION
        # exposures_dead_list = vector("list",length(countries))
        # for (i in seq_along(countries)){
        #     print(i)
        #     deaths_i = dt1[country==countries[i]]
        # 
        #     deaths_i_wide = dcast(deaths_i,date ~ age_group,value.var = "deaths_i_both")
        #     dates = deaths_i_wide[,date] # rename dates
        #     deaths_i_wide = as.matrix(deaths_i_wide[,!"date"],rownames.value=as.character(dates))
        # 
        #     # Convert to sts object
        #     deaths_i_sts = sts(deaths_i_wide,start=c(lubridate::year(min(dates)),yday(min(dates))),frequency = 365)
        #     # # Plot death time series
        #     # plot(deaths_i_sts,type=observed~time)
        # 
        #     # Deconvolve death curve to infection curve of those that died
        #     control = list(k = 0,eps = c(0.01,2),Tmark = nrow(deaths_i_sts) - mean_idd,alpha=0.05,eq3a.method = "C")
        #     idd_pmf_mat = matrix(idd_pmf,nrow = length(idd_pmf),ncol = ncol(deaths_i_sts))
        #     deaths_i_stsBP = backprojNP(deaths_i_sts,idd_pmf_mat,control = control)
        #     # plot(deaths_i_stsBP,xaxis.labelFormat=NULL,legend=NULL,lwd=c(1,1,2),lty=c(1,1,1),ylim=c(0,100),main="")
        # 
        #     exposures_dead = data.table(deaths_i_stsBP@upperbound)
        #     names(exposures_dead) = agegroups
        #     exposures_dead[,date:=dates]
        #     exposures_dead_long = melt(exposures_dead,id.vars = "date",variable.name = "age_group",value.name = "exposures_dead")
        #     exposures_dead_long[,country:=countries[i]]
        # 
        #     exposures_dead_list[[i]] = exposures_dead_long
        # }
        # 
        # exposures_dead = do.call(rbind,lapply(exposures_dead_list,"[[",1))
        
        exposures_dead_list = vector("list",length(dt_list))
        for (i in 1:length(dt_list)){
        # exposures_dead_list = foreach(i=1:length(dt_list)) %dopar% {
            dates = dt_list[[i]][,date]
            
            # Convert to sts object
            deaths_i_sts = sts(dt_list[[i]][,deaths_i_both],start=c(lubridate::year(min(dates)),yday(min(dates))),frequency = 365)
            
            # Deconvolve death curve to infection curve of those that died
            control = list(k = 2,eps = c(0.01,2),Tmark = nrow(deaths_i_sts) - mean_idd,B=10,alpha=0.05,eq3a.method = "C")
            deaths_i_stsBP = backprojNP(deaths_i_sts,idd_pmf,control = control)
            
            # list(deaths_i_stsBP@upperbound,deaths_i_stsBP@lambda)
            exposures_dead_list[[i]] = list(deaths_i_stsBP@upperbound,deaths_i_stsBP@lambda)
        }
        
    } else if (method == "ride"){ # RIDE ALGORITHM
        exposures_dead_list = vector("list",length(dt_list))
        for (i in 1:length(dt_list)){
        # exposures_dead_list = foreach(i=1:length(dt_list)) %dopar% {
            # TODO - check if delay distribution is only defined from day 1 onwards in incidental
            exposures_model = fit_incidence(dt_list[[i]][,as.integer(round(deaths_i_both))],idd_pmf[2:length(idd_pmf)]/sum(idd_pmf[2:length(idd_pmf)]))
            
            # list(exposures_model$Ihat,exposures_model$Isamps)
            exposures_dead_list[[i]] = list(exposures_model$Ihat,exposures_model$Isamps)
        }
        
    } else {
        stop("method must be either 'backproj' or 'ride'.")
    }

    exposures_dead = do.call(c,lapply(exposures_dead_list,"[[",1))
    exposures_dead_samps_list = lapply(exposures_dead_list,"[[",2)
    
    # Add exposures that died to data table
    dt1[,exposures_dead:=exposures_dead]

    return(list(dt1,exposures_dead_samps_list))
}

# Discrete convolution function
# matrix version
# disc_conv = function(x,delay_pmf){
#     dmax = length(delay_pmf) - 1
#     y = matrix(nrow = nrow(x),ncol = ncol(x))
#     for (j in 1:ncol(x)){
#         for (i in 1:nrow(x)){
#             k = 1:i
#             y[i,j] = sum(x[i-k+1] * delay_pmf[k],na.rm=T)
#         }
#     }
#     return(y)
# }
# vector version
disc_conv = function(x,delay_pmf){
    dmax = length(delay_pmf) - 1
    y = numeric(length(x))
    for (i in 1:length(x)){
        k=1:i
        y[i] = sum(x[i-k+1] * delay_pmf[k],na.rm=T)
    }
    return(y)
}

calc_exposures_and_infections = function(dt,ip_pmf){
    dt1 = copy(dt)
    
    # Divide by age- and time-dependent IFR to get exposures
    dt1[,exposures := exposures_dead/ifr_t]
    # Make sure data table is in right order for convolution
    setorder(dt1,country,age_group,date)
    # Convolve exposures with incubation period to get infections
    dt1[,infections := disc_conv(exposures,ip_pmf),by=.(country,age_group)]
    
    return(dt1)
}

calc_cum_exposures_and_infections = function(dt){
    dt1 = copy(dt)
    
    setorder(dt1,country,age_group,date)
    dt1[,`:=`(cum_exp = cumsum(exposures),cum_inf = cumsum(infections)),by=.(country,age_group)]
    dt1[,`:=`(cum_prop_exp = cum_exp/population,cum_prop_inf = cum_inf/population)]
    
    return(dt1)
}


backcalc = function(dt,idd_pmf,ip_pmf,method = "ride"){
    dt1 = copy(dt)
    
    # Deconvolve deaths to get IFR-scaled exposures
    out = deconv(dt1,idd_pmf,method = method)
    backcalc_dt = out[[1]]
    backcalc_samps = out[[2]]
    
    # Calculate exposures and infections
    backcalc_dt = calc_exposures_and_infections(backcalc_dt,ip_pmf)
    
    return(list(backcalc_dt,backcalc_samps))
}

plot_infections = function(dt){
    # plot inferred exposures and infections
    p = ggplot(dt,aes(x=date,group=age_group,color=age_group)) + 
        geom_line(aes(y=infections)) +
        geom_line(aes(y=exposures),linetype="dashed") +
        facet_wrap(~country)
    
    return(p)
}

plot_deaths = function(dt){
    # plot interpolated deaths
    p = ggplot(dt,aes(x=date,y=deaths_i_both,group=age_group,color=age_group)) + 
        ylab("deaths") +
        geom_line() +
        facet_wrap(~country)
    
    return(p)
}

calc_init_condns = function(inc_dt,vax_dt_wide){
    # Assume initially no infected individuals, so numbers in compartments are 
    # all 0 apart from S
    # FOR NOW - assume no waning
    # TODO - include L
    # TODO - include waning. N.B. will also need to update model for IFR
    
    inc_dt1 = copy(inc_dt)
    
    # Calculate differences in numbers entering and leaving compartments at each
    # time point 
    inc_dt1[,`:=`(diffE=nS_E-nE_I,
                 diffIp=nE_Ip-nIp_Is,
                 diffIs=nIp_Is-nIs_R,
                 diffIa=nE_Ia-nIa_R,
                 diffR=nIs_R+nIa_R)]
    # Calculate numbers in each compartment at each time point
    cols = c("E","Ip","Is","Ia","R")
    setorder(inc_dt1,country,vrnt,age_group_model,date)
    inc_dt1[,(cols):=lapply(.SD,cumsum),.SDcols=paste0("diff",cols),by=.(country,vrnt,age_group_model)]
    # # Calculate compartment prevalences
    # inc_dt1[,paste0("prev",cols):=lapply(.SD,function(x) x/population),.SDcols=cols]
    
    incS_dt = inc_dt1[,.(nS_E=sum(nS_E)),by=.(country,age_group_model,date)]
    incS_dt = merge(incS_dt,vax_dt_wide[,!c("prop_va","prop_vb","cum_prop_va","cum_prop_vb")],by=c("country","date","age_group_model"))
    # TODO - update when L is included
    incS_dt[,`:=`(diffS=-nS_E-nS_Va1-nS_Vb1,
                  diffVa1=nS_Va1-nVa1_Va2,
                  diffVa2=nVa1_Va2,
                  diffVb1=nS_Vb1-nVb1_Vb2,
                  diffVb2=nVb1_Vb2)]
    setorder(incS_dt,country,age_group_model,date)
    cols1 = c("S","Va1","Va2","Vb1","Vb2")
    incS_dt[,(cols1):=lapply(.SD,cumsum),.SDcols=paste0("diff",cols1),by=.(country,age_group_model)]
    # N.B. S in line above is still cumulative number of susceptibles that have
    # have lost susceptibility so need to add number initially susceptible 
    # (population)
    incS_dt[,S:=pmax(0,population+S)] # minumum of 0 once everybody has been infected/vaccinated
    # incS_dt[,paste0("prev",cols1):=lapply(.SD,function(x) x/population),.SDcols=cols1]
    
    inc_dt1[,vrnt:=fcase(vrnt=="Other","",
                        vrnt=="Alpha","2",
                        vrnt=="Delta","3")]
    inc_dt_wide = dcast(inc_dt1,country+age_group_model+date+population ~ vrnt,sep = "",value.var = cols)
    
    # Extract current prevalences
    cols2 = c("country","age_group_model","date",cols1)
    prev_dt = merge(inc_dt_wide,incS_dt[,..cols2],by=c("country","age_group_model","date"))
    
    # Calculate compartment prevalences
    cols3 = setdiff(names(prev_dt),c("country","age_group_model","date","population"))
    prev_dt[,paste0("prev",cols3):=lapply(.SD,function(x) x/population),.SDcols=cols3]
    
    ggplot(melt(prev_dt[age_group_model=="50-54"],measure.vars=patterns("prevE"),value.name = "prevE"),aes(x=date,y=prevE,color=variable)) + 
        geom_line() + 
        facet_wrap(~country)
    ggplot(melt(prev_dt[age_group_model=="50-54"],measure.vars=patterns("prevR"),value.name = "prevR"),aes(x=date,y=prevR,color=variable)) + 
        geom_line() + 
        facet_wrap(~country)
    
    return(prev_dt)
}