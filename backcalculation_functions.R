VE = function(x, given = 0){
    (x - given) / (1 - given)
}

get_vaccine_efficacies = function(){
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
    
    return(ve_params)

}

read_death_data = function(source_deaths="coverage"){
    if (source_deaths=="coverage"){
        deaths_raw = fread(cmd = "unzip -cq ../coverage_data/Output_10.zip", skip = 3)
    } else if (source_deaths=="ined"){
        deaths_raw = fread("../ined_data/AgeSex/Cum_deaths_by_age_sex.csv")
    } else {
        stop("Source of age-stratified death data currently unsupported. Please choose 'coverage' or 'ined'.")
    }
}

# Functions for getting min and max age from age group string
get_min_age = function(x) as.numeric(sub("-.*","",sub("\\+|<","-",x)))
get_max_age = function(x) as.numeric(sub(".*-","",sub("\\+|<","-",x)))

divide_na = function(x,y){
    if (!is.na(y)){
        if (y!=0){
            x/y
        } else {
            rep(0,length(x))
        }        
    } else {
        rep(NA,length(x))
    }
}

split_deaths = function(deaths,who_deaths,cntry,start_date,end_date,agg_age_group,cols,direction="back"){
    deaths_cntry = deaths[country==cntry & date>=start_date & date<=end_date & age_group=="Total",.SD,.SDcols = c("country","date","age_group",cols)]
    if (direction=="back"){
        cum_deaths_cntry = who_deaths[country==cntry & date==min(date[country==cntry & date>end_date],na.rm=T),.(country,date,age_group="Total",cum_deaths_both=deaths_total,cum_deaths_male=NA,cum_deaths_female=NA)]    
    } else if (direction=="forward"){
        cum_deaths_cntry = who_deaths[country==cntry & date==max(date[country==cntry & date<start_date],na.rm=T),.(country,date,age_group="Total",cum_deaths_both=deaths_total,cum_deaths_male=NA,cum_deaths_female=NA)]    
    }
    
    deaths_cntry[,`:=`(prop_cum_deaths_both=divide_na(cum_deaths_both,cum_deaths_cntry[,cum_deaths_both]),
                       prop_cum_deaths_male=divide_na(cum_deaths_male,cum_deaths_cntry[,cum_deaths_male]),
                       prop_cum_deaths_female=divide_na(cum_deaths_female,cum_deaths_cntry[,cum_deaths_female])#,
                       # prop_cum_deaths_unknown=divide_na(cum_deaths_unknown,cum_deaths_cntry[,cum_deaths_unknown])
                       )]
    
    if (direction=="back"){
        agegroups = deaths[country==cntry & date==min(date[country==cntry & date>end_date]) & !grepl("Total",age_group),age_group]    
    } else if (direction=="forward") {
        agegroups = deaths[country==cntry & date==start_date-1 & !grepl("Total",age_group),age_group]
    }
    agegroups = agegroups[get_min_age(agegroups)<get_max_age(agg_age_group)]
    deaths_cntry_dt = CJ(country=cntry,age_group=agegroups,date=deaths_cntry[,unique(date)])
    # deaths_cntry_dt = merge(deaths_cntry_dt,deaths_cntry[,.(date,prop_cum_deaths_male,prop_cum_deaths_female,prop_cum_deaths_unknown,prop_cum_deaths_both)],by="date",all.x=T)
    deaths_cntry_dt = merge(deaths_cntry_dt,deaths_cntry[,.(date,prop_cum_deaths_male,prop_cum_deaths_female,prop_cum_deaths_both)],by="date",all.x=T)
    # setnafill(deaths_cntry_dt,fill=0,cols=paste0("prop_",cols))
    
    # deaths_cntry_dt = merge(deaths_cntry_dt,deaths[country==cntry & date==min(date[country==cntry & date>end_date]) & age_group %in% agegroups,.SD,.SDcols=c("region","country_no","age_group",cols)],by=c("age_group"))
    if (direction=="back"){
        next_cum_deaths_cntry = deaths[country==cntry & (age_group %in% agegroups) & date==min(date[country==cntry & date>end_date]),.SD,.SDcols=c("age_group",cols)]    
    } else if (direction=="forward"){
        next_cum_deaths_cntry = deaths[country==cntry & (age_group %in% agegroups) & date==max(date[country==cntry & date<start_date]),.SD,.SDcols=c("age_group",cols)]
    }
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

clean_ined_death_data = function(ined_deaths,who_deaths,cols){
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
    deaths = split_deaths(deaths,who_deaths,"Denmark",as.Date("2020-04-06"),as.Date("2020-08-10"),"<60",cols)
    deaths = split_deaths(deaths,who_deaths,"Denmark",as.Date("2020-04-03"),as.Date("2020-04-05"),"<70",cols)
    
    # Split early deaths in aggregated age groups for Germany
    deaths[country=="Germany" & age_group=="0-59",age_group:="<60"]
    deaths = split_deaths(deaths,who_deaths,"Germany",as.Date("2020-03-30"),as.Date("2020-04-24"),"<60",cols)
    
    # Aggregate 90-99 and 100+ age group data for Germany so for loop in construct_data_table() works
    deaths[country=="Germany" & age_group %in% c("90-99","100+"),age_group:="90+"]
    # TODO - make use of cols robust here, so that only those columns are selected
    # deaths = deaths[,lapply(.SD,sum),.SDcols=cols,by=setdiff(names(deaths),cols)]
    deaths = deaths[,lapply(.SD,sum),.SDcols=cols,by=.(region,country,country_no,excelsource,excelsheet,age_group,death_reference_date,death_reference_date_type,date)]
    
    min_dates = deaths[,.(date=min(date)),by=.(country)]
    countries = deaths[,unique(country)]
    total_deaths_list = vector("list",length(countries))
    cols3 = c("country","date","deaths_total")
    for (i in 1:length(countries)){
        cntry = countries[i]
        tmp = who_deaths[country==cntry & date<min_dates[country==cntry,date],..cols3]
        setnames(tmp,"deaths_total","cum_deaths_both")
        tmp[,age_group:="Total"]
        total_deaths_list[[i]] = tmp
    }
    total_deaths = rbindlist(total_deaths_list)
    deaths = rbind(deaths,total_deaths,fill=T)
    
    for (i in 1:length(countries)){
        cntry = countries[i]
        deaths = split_deaths(deaths,who_deaths,cntry,who_deaths[country==cntry,min(date)],min_dates[country==cntry,date]-1,Inf,cols)
    }
    
    # Exclude totals
    # TODO - update this to incorporate deaths missing ages
    deaths = deaths[!(age_group %in% c("Total","Total known","Total unknown"))]
    
    # Reorder
    setorder(deaths,country,date,age_group)
}

clean_coverage_death_data = function(coverage_deaths,who_deaths,cols){
    deaths = copy(coverage_deaths)
    
    # Change variable names to lower case
    names(deaths) = tolower(names(deaths))
    
    # Change country name for US
    deaths[country=="USA",country:="United States"]
    
    # Select country-level data
    deaths = deaths[region=="All"]
    
    # Rename variables
    cols1 = paste0("cum_",cols)
    setnames(deaths,cols,cols1)
    
    # Convert dates
    deaths[,date:=as.Date(date,format = "%d.%m.%Y")]
    
    # Add age group
    deaths[,age_group:=paste0(age,"-",age+ageint-1)]
    
    # Change sex coding to full words
    deaths[,sex:=fcase(sex=="b","both",
                      sex=="m","male",
                      sex=="f","female")]
    
    # Check data by sex is not available beyond end of data for both sexes 
    # If it is, remove it
    max_dates_sex = dcast(deaths[,.(date=max(date)),by=.(country,sex)],country ~ sex,value.var="date")
    max_dates_sex[,beyond:=(female>both | male>both)]
    countries_beyond = max_dates_sex[beyond==T,.(country,date=both)]
    print(countries_beyond)
    for (i in 1:nrow(countries_beyond)){
        deaths = deaths[!(country==countries_beyond[i,country] & date>countries_beyond[i,date])]
    }
    
    # Remove erroneous entries
    # TODO - correct data for Ireland and Scotland
    deaths = deaths[!(country=="Canada" & date==as.Date("2020-04-30"))] # values appear to be incorrect due to missing "TOT" value for 2020-04-30 in inputDB
    deaths = deaths[!(country=="Denmark" & date==as.Date("2021-01-05"))] # values appear to be incorrect (lower than previous date)
    deaths = deaths[!(country=="Denmark" & date==as.Date("2021-06-03"))] # values are doubled for some reason
    deaths = deaths[!(country=="Finland" & between(date,as.Date("2020-04-10"),as.Date("2020-05-21")))] # aggregated values in raw data for under-60s are split across 10-year age groups
    deaths = deaths[!(country=="France" & (between(date,as.Date("2020-12-06"),as.Date("2020-12-11"))) | (date %in% as.Date(c("2020-12-03","2021-06-04","2021-06-21"))))] # corrections to published records lead to negative new deaths
    deaths = deaths[!(country=="Greece" & (between(date,as.Date("2021-02-03"),as.Date("2021-02-05")) | (date %in% as.Date(c("2021-03-01")))))]
    deaths = deaths[!(country=="Ireland" & date %in% as.Date(c("2020-04-11","2020-04-15","2020-04-16")))] # values are incorrect (lower than previous date)
    deaths = deaths[!(country=="Italy" & (between(date,as.Date("2020-08-11"),as.Date("2020-09-01")) | date %in% as.Date(c("2021-07-15","2021-07-16"))))]
    deaths = deaths[country=="Italy" & date==as.Date("2020-10-20") & age_group=="60-69" & sex=="both",cum_deaths:=3671] # typo in total in inputDB 2671 should be 3671
    deaths = deaths[country=="Italy" & date==as.Date("2020-10-20") & age_group=="60-69" & sex=="male",cum_deaths:=2789] # typo in total in inputDB leads to error in cleaned data
    deaths = deaths[country=="Italy" & date==as.Date("2020-10-20") & age_group=="60-69" & sex=="female",cum_deaths:=882] # typo in total in inputDB leads to error in cleaned data
    deaths = deaths[!(country=="Latvia" & date %in% as.Date(c("2021-07-20","2021-08-17")))] # values are incorrect (lower than previous date)
    deaths = deaths[!(country=="Portugal" & date %in% as.Date(c("2020-10-04","2021-04-24")))] # incorrect all zero or partial zero entries
    deaths = deaths[!(country=="Switzerland" & date==as.Date("2021-05-12"))] # values for 70+ seem to be aggregated
    # deaths = deaths[!(country=="Scotland" & )]
    deaths = deaths[!(country=="Ukraine" & date %in% as.Date(c("2021-05-12","2021-05-13","2021-05-18")))] # missing values for 80+ age groups & lower value for 70-80 age group than previous date
    
    # Calculate total deaths for merging total death counts from WHO time series
    # N.B. Needs some modification to make it work for multiple variables (cases, deaths, tests)
    notallNA = apply(deaths[,mget(cols1)],1,function(x) !all(is.na(x)))
    agg_deaths = deaths[notallNA,lapply(.SD,sum),.SDcols=cols1,by=.(country,region,code,date,sex)]
    agg_deaths[,age_group:="Total"]
    deaths = rbind(deaths[notallNA],agg_deaths,fill=T)
    
    # notallNA1 = apply(deaths[,mget(cols1)],1,function(x) !all(is.na(x)))
    # min_dates = rbind(deaths[notallNA1,.(date=min(date)),by=.(country)],deaths[!(country %in% unique(country[notallNA1])),.(date=min(date)),by=.(country)])
    min_dates = deaths[,.(date=min(date)),by=.(country)]
    max_dates = deaths[,.(date=max(date)),by=.(country)]
    countries = deaths[,unique(country)]
    total_deaths_list = vector("list",length(countries))
    cols2 = paste0(sub("tests","tested",cols),"_total")
    cols3 = c("country","date",cols2)
    for (i in 1:length(countries)){
        cntry = countries[i]
        tmp = who_deaths[country==cntry & (date<min_dates[country==cntry,date] | date>max_dates[country==cntry,date]),..cols3]
        total_deaths_list[[i]] = tmp
    }
    total_deaths = rbindlist(total_deaths_list)
    # Add WHO data for time periods for Finland and Italy with missing data
    tmp = who_deaths[(country=="Finland" & between(date,as.Date("2020-04-10"),as.Date("2020-05-21")))|
                          (country=="Italy" & between(date,as.Date("2021-01-12"),as.Date("2021-02-24"))),..cols3]
    total_deaths = rbind(total_deaths,tmp)
    setnames(total_deaths,cols2,cols1)
    total_deaths[,`:=`(region="All",sex="both",age_group="Total")]
    deaths = rbind(deaths,total_deaths,fill=T)
        
    # Cast to wide format for compatibility with construct_data_table function
    if (length(cols1)==1){
        deaths[,sex:=paste0(cols1,"_",sex)]
    }
    deaths = dcast(deaths,country + code + date + age_group ~ sex, value.var = cols1)
    
    for (i in 1:length(countries)){
        cntry = countries[i]
        if (cntry %in% who_deaths[,unique(country)]){
            # Fill in early deaths for each country by scaling the cumulative number of 
            # deaths up to the first date with disaggregated age groups for each age 
            # group by the distribution of the overall number of daily deaths in the WHO
            # data over time up to that point
            deaths = split_deaths(deaths,who_deaths,cntry,who_deaths[country==cntry,min(date)],min_dates[country==cntry,date]-1,Inf,paste0(cols1,"_",c("both","male","female")))
            # Fill in recent deaths for each country by scaling the cumulative number of 
            # deaths up to the last date with disaggregated age groups for each age 
            # group by the distribution of the overall number of daily deaths in the WHO
            # data over time since then
            deaths = split_deaths(deaths,who_deaths,cntry,max_dates[country==cntry,date]+1,who_deaths[country==cntry,max(date)],Inf,paste0(cols1,"_",c("both","male","female")),direction = "forward")
        }
    }
    # Do the same for periods with missing data for Finland and Italy
    deaths = split_deaths(deaths,who_deaths,"Finland",as.Date("2020-04-10"),as.Date("2020-05-21"),Inf,paste0(cols1,"_",c("both","male","female")))
    deaths = split_deaths(deaths,who_deaths,"Italy",as.Date("2021-01-12"),as.Date("2021-02-24"),Inf,paste0(cols1,"_",c("both","male","female")))
    
    # Reorder
    setorder(deaths,country,date,age_group)
    
    # Drop rows with overall cumulative deaths
    deaths = deaths[age_group!="Total"]
    
    return(deaths)
}

clean_death_data = function(source_deaths,deaths_raw,who_deaths){
    if (source_deaths=="coverage"){
        cols = "deaths" # variables to be included in cleaned data
        deaths = clean_coverage_death_data(deaths_raw,who_deaths,cols)
    } else if (source_deaths=="ined"){
        cols = c("cum_deaths_male","cum_deaths_female","cum_deaths_both") # variables to be included in cleaned data
        deaths = clean_ined_death_data(deaths_raw,who_deaths,cols)
    } else {
        stop("Source of age-stratified death data currently unsupported. Please choose 'coverage' or 'ined'.")
    }
}

process_variant_data = function(vrnt_data){
    vrnt_data1 = copy(vrnt_data)
    
    # Remove unknown variant entries
    vrnt_data1 = vrnt_data[variant!="UNK"]
    
    # Classify variants into non-Alpha-Delta, Alpha and Delta
    vrnt_data1[,vrnt:="Other"]
    vrnt_data1[variant=="B.1.1.7",vrnt:="Alpha"]
    vrnt_data1[variant=="B.1.617.2",vrnt:="Delta"]
    
    # Aggregate data by Alpha/Delta status
    # vrnt_data1 = vrnt_data1[source=="GISAID",.(number_detections_variant=sum(number_detections_variant)),by=.(country,country_code,year_week,date,new_cases,number_sequenced,percent_cases_sequenced,vrnt)]
    vrnt_data1 = vrnt_data1[(country!="Hungary" & source=="GISAID")|(country=="Hungary" & ((source=="GISAID" & year_week<"2021-24")|(source=="TESSy" & year_week>="2021-24"))),
                            .(number_detections_variant=sum(number_detections_variant,na.rm=T)),
                            by=.(country,country_code,year_week,date,new_cases,number_sequenced,number_sequenced_known_variant,percent_cases_sequenced,vrnt)]
    # vrnt_data1 = vrnt_data1[(!(country %in% c("Cyprus","Malta")) & (country!="Hungary" & source=="GISAID")|(country=="Hungary" & ((source=="GISAID" & year_week<"2021-24")|(source=="TESSy" & year_week>="2021-24"))))|
    #                             ((country %in% c("Cyprus","Malta")) & source=="TESSy"),
    #                         .(number_detections_variant=sum(number_detections_variant,na.rm=T)),
    #                         by=.(country,country_code,year_week,date,new_cases,number_sequenced,number_sequenced_known_variant,percent_cases_sequenced,vrnt)]
    vrnt_data1[,prop_vrnt:=number_detections_variant/number_sequenced_known_variant]
    
    # FOR NOW - drop data for Slovakia for dates with low numbers of sequences showing low proportion of WT
    # TODO - work out how to do this in a more principled way
    vrnt_data1 = vrnt_data1[!(country=="Slovakia" & year_week %in% c("2020-53","2021-01","2021-02","2021-06","2021-10"))]
    
    return(vrnt_data1)
}

process_cog_variant_data = function(cog_vrnt_data){
    cog_vrnt_data1 = copy(cog_vrnt_data)
    cog_vrnt_data1[,vrnt:="Other"]
    cog_vrnt_data1[Lineage=="B.1.1.7",vrnt:="Alpha"]
    cog_vrnt_data1[Lineage %in% c("B.1.617.2","AY.4","AY.5","AY.6","AY.9"),vrnt:="Delta"]
    # Exclude sequences with unknown lineage from 2021-08-14 counts, as number is very high
    cog_vrnt_data1 = cog_vrnt_data1[!(WeekEndDate==as.Date("2021-08-14") & Lineage=="None")]
    vrnt_dataENG = CJ(date = cog_vrnt_data1[,unique(WeekEndDate)], vrnt = cog_vrnt_data1[,unique(vrnt)])
    vrnt_dataENG = merge(vrnt_dataENG,cog_vrnt_data1[,.(number_detections_variant=sum(Count)),by=.(date=WeekEndDate,vrnt)],all.x=T)
    vrnt_dataENG[is.na(number_detections_variant),number_detections_variant:=0]
    vrnt_dataENG[,prop_vrnt:=number_detections_variant/sum(number_detections_variant),by=.(date)]
    vrnt_dataENG[,country:="England"]
    
    return(vrnt_dataENG)
}

plot_variant_data = function(vrnt_data,src){
    # Plot variant proportions over time
    p = ggplot(vrnt_data[source==src],aes(x=date,y=percent_variant,group=variant,color=variant)) +
        geom_line() +
        labs(title = src) +
        facet_wrap(~country)
    return(p)
}

estimate_variant_proportions = function(vrnt_data,max_date){
    # Fit multinomial logistic model to estimate variant proportions over time
    # Cast to wide format
    vrnt_data_wide = dcast(vrnt_data,country + date ~ vrnt,value.var = "number_detections_variant")
    # Exclude rows without any observations
    vrnt_data_wide = vrnt_data_wide[!(Alpha==0 & Delta==0 & Other==0)]
    
    countries = vrnt_data_wide[,unique(country)]
    ncountries = length(countries)
    m1 = vector("list",ncountries)
    m2 = vector("list",ncountries)
    for (i in 1:ncountries){
        cntry = countries[i]
        print(cntry)
        m1[[i]] = multinom(as.matrix(vrnt_data_wide[country==cntry,.(Other,Alpha,Delta)]) ~ date,vrnt_data_wide[country==cntry],maxit=10000)
        m2[[i]] = multinom(as.matrix(vrnt_data_wide[country==cntry,.(Other,Alpha,Delta)]) ~ ns(date,df=2),vrnt_data_wide[country==cntry],maxit=10000)
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
    new_dt = data.table(date=seq.Date(vrnt_data_wide[,min(date)],max_date,by=1)) #data.table(date = seq.Date(vrnt_data_wide[,min(date)],dates2[length(dates2)],by=1))
    for (i in 1:ncountries){
        cntry = countries[i]
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
    
    return(list(vrnt_prop1=vrnt_prop1,vrnt_prop2=vrnt_prop2))
    # TODO - work out how to calculate CIs with effects package. Might need to
    # reformat input data to get it to work
}

plot_variant_proportions = function(vrnt_prop,vrnt_data){
    p = ggplot() +
        geom_line(aes(x=date,y=prop_vrnt,color=vrnt),vrnt_prop) +
        geom_point(aes(x=date,y=prop_vrnt,color=vrnt),vrnt_data) +
        labs(x="Date",y="Proportion of sequences",color="Variant") +
        theme(axis.text.x=element_text(angle=45,hjust=1)) + 
        facet_wrap(~country)
    return(p)
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
    # vax = melt(vax_wide,measure.vars = patterns(first="first",second="second"),variable.name = "type")
    # vax[,type:=fcase(type==1,"va",type==2,"vb")]
    vax = melt(vax_wide,measure.vars = c("first_va","first_vb","second_va","second_vb"),variable.name = "dose",value.name = "count")
    vax[,`:=`(type=sub("[a-z]+_","",dose),dose=sub("_[a-z]+","",dose))]
    vax[,dose:=fifelse(dose=="first",1,2)]
    
    # Remove rows with missing entries for both first and second doses
    # vax = vax[!(is.na(first) & is.na(second))]
    vax = vax[!is.na(count)]
    
    # # Calculate proportions of each vaccine type over time
    # vax[,`:=`(prop_first=first/sum(first),prop_second=second/sum(second)),by=.(country,date,age_group)]
    
    # See how many doses have unknown age group
    print(vax[,sum(count)])
    print(vax[age_group!="UNK",sum(count)])
    print(vax[,.(count=sum(count)), by = .(age_group)]) # only small % so ignore FOR NOW
    
    # return(list(vax=vax,num_type=num_type))
    return(vax)
}

process_phe_vaccination_data = function(vaccPHE,agegroups_model){
    vaxENG = vector("list",length(vaccPHE))
    for (i in 1:length(vaxENG)){
        vaxENG[[i]] = data.table(date=rep(vaccPHE[[i]]$vt,each=length(agegroups_model)),
                                 age_group=rep(agegroups_model,length(vaccPHE[[i]]$vt)),
                                 va1=unlist(vaccPHE[[i]]$va1),
                                 va2=unlist(vaccPHE[[i]]$va2),
                                 vb1=unlist(vaccPHE[[i]]$vb1),
                                 vb2=unlist(vaccPHE[[i]]$vb2),
                                 pop=i)
    }
    vaxENG = rbindlist(vaxENG,fill=T)
    vaxENG = vaxENG[,lapply(.SD,sum),.SDcols=c("va1","va2","vb1","vb2"),by=.(date,age_group)]
    vaxENG[,`:=`(country="England",year_week_iso=date2ISOweek(date))]
    vaxENG_long = melt(vaxENG,measure.vars = c("va1","va2","vb1","vb2"),variable.name = "type",value.name = "count")
    vaxENG_long[,`:=`(type=substr(type,1,2),dose=substr(type,3,3))]
    # vaxENG_long[,dose:=fifelse(dose==1,"first","second")]
    # vaxENG = dcast(vaxENG_long,country+year_week_iso+date+age_group+type ~ dose)
    
    # # Plot cumulative number vaccinated to check
    # ggplot(vaxENG[type=="va",.(date,second=cumsum(second)),by=.(age_group)],aes(x=date,y=second,color=age_group)) + geom_line()
    # ggplot(vaxENG[type=="vb",.(date,second=cumsum(second)),by=.(age_group)],aes(x=date,y=second,color=age_group)) + geom_line()
    
    # Remove delays for Ab development
    vaxENG_long[dose==1,date:=date-28]
    vaxENG_long[dose==2,date:=date-14]
    
    # Restrict to dates for which there is data for both first and second doses given Ab development delays
    vaxENG_long = vaxENG_long[between(date,vaxENG[,min(date)]-14,vaxENG[,max(date)]-28)]
    
    return(vaxENG_long)
}

construct_vax_data_table = function(vax,dates,agegroups,pop){
    # Get vaccination age groups from vax
    agegroups_vax = sort(vax[,unique(age_group)])
    agegroups_vax = agegroups_vax[agegroups_vax!="UNK"]
    min_ages_vax = get_min_age(agegroups_vax)
    min_ages_vax[is.na(min_ages_vax)] = 0
    
    # Get minimum ages of age groups
    min_ages = get_min_age(agegroups)
    
    # Create data table with all combinations of countries, dates, vaccines and ages
    # to store vaccine schedule
    vax_dt = CJ(country=vax[,unique(country)],date=dates,age=0:100,type=vax[,unique(type)],dose=vax[,unique(dose)])
    vax_dt = merge(vax_dt,pop[,.(country,age,population)],by=c("country","age"),all.x=T)
    vax_dt[,year_week_iso:=ISOweek(date)]
    vax_dt[,age_group:=cut(age,c(min_ages_vax,Inf),labels=agegroups_vax,right=F)]
    
    # Merge with vax data table
    # N.B. This duplicates numbers of doses for the same ISO week and age group
    # for the ECDC data, so we then divide by 7 to get the average doses per day 
    # and divide doses between age groups according to population fraction
    if (vax[,length(unique(country))]==1 && vax[,unique(country)]=="England"){
        vax_dt = merge(vax_dt,vax[,!"year_week_iso"],by=c("country","date","age_group","type","dose"),all.x=T)
        vax_dt[,count:=as.numeric(count)]
    } else {
        vax_dt = merge(vax_dt,vax[,!"date"],by=c("country","year_week_iso","age_group","type","dose"),all.x=T)
        vax_dt[,count:=as.numeric(count)]
        vax_dt[,count:=count/7]
    }
    vax_dt[,count:=count*population/sum(population),by=.(country,date,age_group,type,dose)]
    print(vax_dt[,sum(count,na.rm=T)])
    
    # Change age groups
    vax_dt[,age_group:=cut(age,c(min_ages,Inf),agegroups,right=F)]
    # Sum vaccinations in each age group by type and dose
    vax_dt = vax_dt[,lapply(.SD,function(x) sum(x,na.rm=T)),.SDcols=c("population","count"),by=.(country,date,age_group,type,dose)]
    # Calculate proportion vaccinated
    vax_dt[,`:=`(prop=count/population)]
    # Calculate cumulative proportion vaccinated ensuring data is in chronological date order
    setorder(vax_dt,country,date,age_group,type,dose)
    vax_dt[,`:=`(cum_prop=cumsum(prop)),by=.(country,age_group,type,dose)]
    
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
    
    # Calculate standard deviation of posterior distribution of IFR assuming it is normal
    ifr_dt[,sigma:=(ifr-ifr_lb)/qnorm(0.975)]
}

# construct_data_table = function(agegroups,deaths,pop,cols,ltc_deaths,vax,num_type,ifr){
construct_data_table = function(agegroups,deaths,pop,cols,ltc_deaths,vax,Ab_delay1,Ab_delay2,vaxENG,ifr,vrnt_prop,ve_params){
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
        date_range = base_deaths_dt[country==cntry,range(date)]
        # # Start dates from 60 days before first death record to allow for maximum possible infection-to-death time in deconvolution below
        # tmp = CJ(country = cntry, date = seq.Date(date_range[1]-60,date_range[2],by=1), age_group = agegroups)
        # Start dates from first date in WHO data (as this is well before first recorded death for all European countries, so sufficiently early to account for infection-to-death delay)
        if (cntry=="England"){ # use minimum of dates for other countries for England, as it does not appear in WHO data
            tmp = CJ(country = cntry, date = seq.Date(base_deaths_dt[,min(date)],date_range[2],by=1), age_group = agegroups)
            deaths_dt_list[[i]] = merge(tmp,base_deaths_dt[country==cntry,!"population"],by=c("country","date","age_group"),all.x=T)
            deaths_dt_list[[i]][date<date_range[1],(cols):=0]
        } else {
            tmp = CJ(country = cntry, date = seq.Date(date_range[1],date_range[2],by=1), age_group = agegroups)
            deaths_dt_list[[i]] = merge(tmp,base_deaths_dt[country==cntry,!"population"],by=c("country","date","age_group"),all.x=T)
        }
    }
    deaths_dt = rbindlist(deaths_dt_list)

    # Exclude country and variable combinations with all missing data
    deaths_dt_long = melt(deaths_dt,id.vars = c("country","date","age_group"),variable.name="sex",value.name="cum_deaths")
    deaths_dt_long[,sex:=sub("cum_deaths_","",sex)]
    # countries_notallNA = deaths_dt_long[,!all(is.na(cum_deaths)),by=.(country,sex)][V1==T,.(country,sex)]
    # deaths_dt_long = deaths_dt_long[countries_notallNA,on=c("country","sex")]
    morethan1obs = deaths_dt_long[!is.na(cum_deaths),.(length(unique(date))),by=.(country,sex)][V1>1,.(country,sex)]
    deaths_dt_long = deaths_dt_long[morethan1obs,on=c("country","sex")]

    # # Fill in cumulative deaths for earliest date so there is a value to interpolate from
    # countries1 = deaths_dt_long[,unique(country)]
    # for (i in 1:length(countries1)){
    #     deaths_dt_long[country==countries1[i] & date==min(date[country==countries1[i]]),cum_deaths:=0]    
    # }

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
    
    # Calculate proportion of deaths among those 60+
    # Exclude France as deaths time series is for hospital deaths only
    deaths_dt[country!="France",prop_60plus := sum(deaths_i_both[age_group %in% agegroups[max_ages>60]])/sum(deaths_i_both),by=.(country)]
    # Calculate proportion of deaths among those 60+ in France from non-LTC deaths and proportion of deaths in LTCs
    # N.B. This overestimates total deaths in France by ~10,000 as the proportion of deaths in LTCs has decreased over time
    # TODO - make proportion of deaths in LTCs time dependent e.g. using INED data
    deaths_dt[country=="France",prop_60plus := prop_ltc_deaths + (1-prop_ltc_deaths)*sum(deaths_i_both[age_group %in% agegroups[max_ages>60]])/sum(deaths_i_both)]
    
    # Calculate LTC deaths among 60+ year-olds by:
    # 1. assuming all LTC deaths occur among those 60+
    # 2. dividing deaths among 60+ by proportion of total deaths among 60+ (to get total deaths), and 
    # 3. multiplying by proportion of total deaths that have occurred among LTC residents
    # (N.B. relies on having proportion of deaths in LTCs or LTC residents)
    # TODO - rethink/modify this calculation for non-European countries w/o care homes
    deaths_dt[age_group %in% agegroups[max_ages<60],deaths_i_both_ltc := 0] # assume no deaths in under-60s
    # Calculate deaths among those 60+ for France by scaling up hospital (non-LTC) deaths
    deaths_dt[country=="France" & age_group %in% agegroups[max_ages>60],deaths_i_both := deaths_i_both/(1-prop_ltc_deaths/prop_60plus)]
    deaths_dt[age_group %in% agegroups[max_ages>60],deaths_i_both_ltc := deaths_i_both * fifelse(!is.na(prop_ltc_res_deaths),prop_ltc_res_deaths,prop_ltc_deaths)/prop_60plus]
    # Calculate non-LTC deaths
    deaths_dt[,deaths_i_both_non_ltc := deaths_i_both - deaths_i_both_ltc]
    
    # VACCINATION DATA
    
    # Get ISO week of earliest date with death record from deaths data
    min_iso_week = ISOweek(deaths_dt[,min(date)])
    # Get date of start of ISO week to ensure length of date vector is a multiple of 7
    min_date = ISOweek2date(paste0(min_iso_week,"-1"))
    # Get last date for which vaccine data is available
    max_date = vax[,max(date)+6]
    # Make daily date sequence from earliest death date to latest date for which 
    # vaccine data is available - length should be a multiple of 7
    dates1 = seq.Date(min_date-ceiling(max(Ab_delay1,Ab_delay2)/7)*7,max_date,by=1) # -Ab_delay to account for Ab_delay-day delay to Ab development, and max(date) + 6 to get end of last ISO week
    
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
    vax_dt = construct_vax_data_table(vax,dates1,agegroups,pop)
    
    # Construct vaccination data table for England
    vaxENG_dt = construct_vax_data_table(vaxENG,dates1,agegroups,pop)
    
    # Bind together
    vax_dt = rbind(vax_dt,vaxENG_dt)
    
    # Plot data
    # cumulative proportion fully vaccinated
    print(ggplot(vax_dt[date>=vax[,min(date)] & dose==2,.(cum_prop=sum(cum_prop)),by=.(country,age_group,date)],aes(x=date,y=cum_prop,group=age_group,color=age_group)) + 
        geom_line() +
        facet_wrap(~country))
    
    # proportion newly fully vaccinated
    print(ggplot(vax_dt[date>=vax[,min(date)] & dose==2,.(prop=sum(prop)),by=.(country,age_group,date)],aes(x=date,y=prop,group=age_group,color=age_group)) + 
        geom_line() +
        facet_wrap(~country))
    
    
    # Merge death and vaccination data
    # See which countries are in both datasets
    print(intersect(deaths_dt[,country],vax_dt[,country]))
    
    # Add lags for Ab development to date
    vax_dt[dose==1,date:=date+Ab_delay1]
    vax_dt[dose==2,date:=date+Ab_delay2]
    
    # Cast vaccination data to wide format
    vax_dt_wide = dcast(vax_dt,country + date + age_group + population ~ type + dose,value.var = c("count","prop","cum_prop"))
    
    # Restrict to dates for which there is data for both first and second doses given Ab development delays
    vax_dt_wide = vax_dt_wide[between(date,min_date,max_date+min(Ab_delay1,Ab_delay2))]
    
    # # Add Ab_delay days to date in vaccination data for development of Ab
    # vax_dt_wide[,date_v:=date]
    # vax_dt_wide[,date:=date+Ab_delay]
    # Make vector of dates that are covered by both the death data and the 
    # vaccination data
    # Use the last vaccination date as the cut-off for the vaccination data 
    # rather than the last Ab development date as we need to use the vaccination
    # data later for the numbers entering the vaccination compartments
    dates2 = as.Date(intersect(deaths_dt[,date],vax_dt_wide[,date]),origin="1970-01-01")
    deaths_vax_dt = merge(deaths_dt[date %in% dates2],vax_dt_wide[date %in% dates2],by=c("country","date","age_group"))
    
    # IFR ESTIMATES
    # ifr_dt = copy(base_dt)
    # 
    # # # Add IFR estimate from Levin et al Eur Jrnl Epi 2020
    # # ifr_dt[,`:=`(ifr=10^(-3.27+0.0524*age)/100)]
    # # # Calculate 95% CI for IFR assuming FOR NOW no covariance between intercept and 
    # # # slope estimates
    # # # TODO - estimate intercept and slope covariance from lower and upper bounds in 
    # # # supplementary spreadsheet in Levin et al 2020
    # # se_intcpt = 0.07*log(10)
    # # se_slope = 0.0013*log(10)
    # # ifr_dt[,se_ifr:=sqrt(se_intcpt^2+se_slope^2*age^2)*ifr]
    # # ifr_dt[,`:=`(ifr_lb=pmax(0,ifr-qnorm(0.975)*se_ifr),ifr_ub=pmin(1,ifr+qnorm(0.975)*se_ifr))]
    # # 
    # # ifr_dt[,age_group:=cut(age,c(min_ages,Inf),labels=agegroups,right=F)]
    # # cols2 = c("ifr","se_ifr","ifr_lb","ifr_ub")
    # # ifr_dt=ifr_dt[,lapply(.SD,function(x) sum(x*population)/sum(population)),.SDcols=cols2,by=.(country,age_group)]
    # 
    # names(ifr) = tolower(names(ifr))
    # min_ages_ifr = get_min_age(ifr[,age_group])
    # 
    # # Add IFR age groups
    # ifr_dt[,age_group:=cut(age,c(min_ages_ifr,Inf),labels=ifr[,age_group],right=F)]
    # 
    # # Merge with IFR data table
    # ifr_dt = merge(ifr_dt,ifr[,.(age_group,ifr=median_perc_mean/100,ifr_lb=ci_95_lb_mean/100,ifr_ub=ci_95_ub_mean/100)],by="age_group")
    # 
    # # Change age groups
    # ifr_dt[,age_group := cut(age,c(min_ages,Inf),labels=agegroups,right=F)]
    # 
    # # Calculate population-weighted average for each age group
    # cols2 = c("ifr","ifr_lb","ifr_ub")
    # ifr_dt = ifr_dt[,lapply(.SD,function(x) sum(x*population)/sum(population)),.SDcols=cols2,by=c("country","age_group")]    
    
    ifr_dt = construct_ifr_data_table(ifr,base_dt,min_ages,agegroups)
    
    # Plot country IFRs
    ggplot(ifr_dt,aes(x=age_group,y=ifr,group=country,color=country)) + 
        geom_line() +
        # scale_y_log10() +
        theme(legend.position = "none")
        
    # Merge IFR with death and vaccination data
    dt = merge(deaths_vax_dt,ifr_dt,by=c("country","age_group"))
    # # N.B. Countries without vax data dropped here - check this!
    # dt = merge(dt,num_type[,.(country,prop_va,prop_vb)],by="country")
    
    # ESTIMATED VARIANT PROPORTIONS
    # Cast variant proportion data table to wide format
    vrnt_prop_wide = dcast(vrnt_prop,country+date~vrnt,value.var = "prop_vrnt")
    cols6 = paste0("prop_vrnt",c("","2","3"))
    setnames(vrnt_prop_wide,c("Other","Alpha","Delta"),cols6)
    
    # Merge with overall data table
    dt = merge(dt,vrnt_prop_wide,by=c("country","date"),all.x=T)
    # Backfill variant proportions with first non-NA observation
    dt[,(cols6):=nafill(.SD,type="nocb"),.SDcols=cols6,by=.(country)]
    
    # # FOR NOW - use efficacy values for wild-type variant (N.B. same as for Delta)
    # # TODO - use proportion of each vaccine by age group and time and proportion of 
    # # each variant over time to calculate more detailed age- and time-dependent IFR
    # # dt[,`:=`(ei=prop_va*ei_va2+prop_vb*ei_vb2,ed=prop_va*ed_va2i+prop_vb*ed_vb2i,em=prop_va*em_va2d+prop_vb*em_vb2d)]
    # dt[,`:=`(ed=prop_va*ed_va2i+prop_vb*ed_vb2i,em=prop_va*em_va2d+prop_vb*em_vb2d)]
    # dt[,`:=`(ei=fifelse(!(prop_va==0 & prop_vb==0),prop_va/(prop_va+prop_vb)*(prop_vrnt*ve_params$ei_va2+prop_vrnt2*ve_params$ei2_va2+prop_vrnt3*ve_params$ei3_va2) +
    #                         prop_vb/(prop_va+prop_vb)*(prop_vrnt*ve_params$ei_vb2+prop_vrnt2*ve_params$ei2_vb2+prop_vrnt3*ve_params$ei3_vb2),0),
    #          ed=fifelse(!(prop_va==0 & prop_vb==0),prop_va/(prop_va+prop_vb)*(prop_vrnt*ve_params$ed_va2i+prop_vrnt2*ve_params$ed_va2i2+prop_vrnt3*ve_params$ed_va2i3) +
    #                         prop_vb/(prop_va+prop_vb)*(prop_vrnt*ve_params$ed_vb2i+prop_vrnt2*ve_params$ed_vb2i2+prop_vrnt3*ve_params$ed_vb2i3),0),
    #          em=fifelse(!(prop_va==0 & prop_vb==0),prop_va/(prop_va+prop_vb)*(prop_vrnt*ve_params$em_va2d+prop_vrnt2*ve_params$em_va2d2+prop_vrnt3*ve_params$em_va2d3) +
    #                         prop_vb/(prop_va+prop_vb)*(prop_vrnt*ve_params$em_vb2d+prop_vrnt2*ve_params$em_vb2d2+prop_vrnt3*ve_params$em_vb2d3),0),
    #          et=fifelse(!(prop_va==0 & prop_vb==0),prop_va/(prop_va+prop_vb)*(prop_vrnt*ve_params$et_va2i+prop_vrnt2*ve_params$et_va2i2+prop_vrnt3*ve_params$et_va2i3) +
    #                         prop_vb/(prop_va+prop_vb)*(prop_vrnt*ve_params$et_vb2i+prop_vrnt2*ve_params$et_vb2i2+prop_vrnt3*ve_params$et_vb2i3),0))]
    dt[,`:=`(p_va1=pmax((cum_prop_va_1-cum_prop_va_2)/(cum_prop_va_1+cum_prop_vb_1),0),
             p_va2=pmin(cum_prop_va_2/(cum_prop_va_1+cum_prop_vb_1),1),
             p_vb1=pmax((cum_prop_vb_1-cum_prop_vb_2)/(cum_prop_va_1+cum_prop_vb_1),0),
             p_vb2=pmin(cum_prop_vb_2/(cum_prop_va_1+cum_prop_vb_1),1))]
    dt[,`:=`(ei=fifelse(!(cum_prop_va_1==0 & cum_prop_vb_1==0),
                            p_va1*(prop_vrnt*ve_params$ei_va1+prop_vrnt2*ve_params$ei2_va1+prop_vrnt3*ve_params$ei3_va1) +
                            p_va2*(prop_vrnt*ve_params$ei_va2+prop_vrnt2*ve_params$ei2_va2+prop_vrnt3*ve_params$ei3_va2) +
                            p_vb1*(prop_vrnt*ve_params$ei_vb1+prop_vrnt2*ve_params$ei2_vb1+prop_vrnt3*ve_params$ei3_vb1) +
                            p_vb2*(prop_vrnt*ve_params$ei_vb2+prop_vrnt2*ve_params$ei2_vb2+prop_vrnt3*ve_params$ei3_vb2),0),
             ed=fifelse(!(cum_prop_va_1==0 & cum_prop_vb_1==0),
                            p_va1*(prop_vrnt*ve_params$ed_va1i+prop_vrnt2*ve_params$ed_va1i2+prop_vrnt3*ve_params$ed_va1i3) +
                            p_va2*(prop_vrnt*ve_params$ed_va2i+prop_vrnt2*ve_params$ed_va2i2+prop_vrnt3*ve_params$ed_va2i3) +
                            p_vb1*(prop_vrnt*ve_params$ed_vb1i+prop_vrnt2*ve_params$ed_vb1i2+prop_vrnt3*ve_params$ed_vb1i3) +
                            p_vb2*(prop_vrnt*ve_params$ed_vb2i+prop_vrnt2*ve_params$ed_vb2i2+prop_vrnt3*ve_params$ed_vb2i3),0),
             em=fifelse(!(cum_prop_va_1==0 & cum_prop_vb_1==0),
                            p_va1*(prop_vrnt*ve_params$em_va1d+prop_vrnt2*ve_params$em_va1d2+prop_vrnt3*ve_params$em_va1d3) +
                            p_va2*(prop_vrnt*ve_params$em_va2d+prop_vrnt2*ve_params$em_va2d2+prop_vrnt3*ve_params$em_va2d3) +
                            p_vb1*(prop_vrnt*ve_params$em_vb1d+prop_vrnt2*ve_params$em_vb1d2+prop_vrnt3*ve_params$em_vb1d3) +
                            p_vb2*(prop_vrnt*ve_params$em_vb2d+prop_vrnt2*ve_params$em_vb2d2+prop_vrnt3*ve_params$em_vb2d3),0),
             et=fifelse(!(cum_prop_va_1==0 & cum_prop_vb_1==0),
                            p_va1*(prop_vrnt*ve_params$et_va1i+prop_vrnt2*ve_params$et_va1i2+prop_vrnt3*ve_params$et_va1i3) +
                            p_va2*(prop_vrnt*ve_params$et_va2i+prop_vrnt2*ve_params$et_va2i2+prop_vrnt3*ve_params$et_va2i3) +
                            p_vb1*(prop_vrnt*ve_params$et_vb1i+prop_vrnt2*ve_params$et_vb1i2+prop_vrnt3*ve_params$et_vb1i3) +
                            p_vb2*(prop_vrnt*ve_params$et_vb2i+prop_vrnt2*ve_params$et_vb2i2+prop_vrnt3*ve_params$et_vb2i3),0))]
    
    # Calculate IFRs under vaccination with different vaccine types
    # dt[,ifr_v:=(1-ed)*(1-em)*ifr]
    # dt[,`:=`(ifr_va=(prop_vrnt*(1-ed_va2i)*(1-em_va2d) +
    #                      prop_vrnt2*(1-ed_va2i2)*(1-em_va2d2) + 
    #                      prop_vrnt3*(1-ed_va2i3)*(1-em_va2d3))*ifr,
    #          ifr_vb=(prop_vrnt*(1-ed_vb2i)*(1-em_vb2d) +
    #                      prop_vrnt2*(1-ed_vb2i2)*(1-em_vb2d2) +
    #                      prop_vrnt3*(1-ed_vb2i3)*(1-em_vb2d3))*ifr)]
    
    # Calculate IFR over time
    # dt[,ifr_t := ifr]
    dt[,cum_prop_v:=cum_prop_va_1+cum_prop_vb_1]
    dt[,ifr_t:=pmax((1-(1-ei)*cum_prop_v/(1-ei*cum_prop_v)+(1-ei)*cum_prop_v/(1-ei*cum_prop_v)*(1-ed)*(1-em))*ifr,0)] # set to 0 where IFR is negative for older age groups for Iceland (due to cum_prop_v>1), as IFR can't be negative
    # dt[,ifr_t := (1-cum_prop_va-cum_prop_vb)*ifr + cum_prop_va*ifr_va + cum_prop_vb*ifr_vb]
    # dt[,ifr_t := (1-cum_prop_va-cum_prop_vb)*(prop_vrnt+prop_vrnt2+prop_vrnt3)*ifr + 
    #        cum_prop_va*ifr_va + 
    #        cum_prop_vb*ifr_vb]
    
    # # Calculate uncertainty bounds on time-varying IFR
    # dt[,ifr_t_lb:=pmax((1-(1-ei)*cum_prop_v/(1-ei*cum_prop_v)+(1-ei)*cum_prop_v/(1-ei*cum_prop_v)*(1-ed)*(1-em))*ifr_lb,0)]
    # dt[,ifr_t_ub:=pmax((1-(1-ei)*cum_prop_v/(1-ei*cum_prop_v)+(1-ei)*cum_prop_v/(1-ei*cum_prop_v)*(1-ed)*(1-em))*ifr_ub,0)]
    
    # Plot population-weighted average IFR over time for all countries
    print(ggplot(dt[,.(ifr_t=sum(ifr_t*population)/sum(population)),by=.(country,date)],
           aes(x=date,y=ifr_t,group=country,color=country)) + geom_line())
    # ggsave("./output/avg_ifr_over_time.png",width = 5,height = 4)
    
    return(dt)
}


plot_deaths = function(dt){
    # Plot interpolated deaths
    p = ggplot(dt,aes(x=date,y=deaths_i_both,group=age_group,color=age_group)) + 
        labs(x="Date",y="Deaths",color="Age group") +
        geom_line() +
        theme(axis.text.x=element_text(angle=45,hjust=1)) +
        facet_wrap(~country)
    return(p)
}

plot_ovrl_vax_cov = function(dt){
    p = ggplot(dt[,.(cum_prop_v=sum(cum_prop_v*population)/sum(population)),by=.(country,date)],
               aes(x=date,y=cum_prop_v,group=country,color=country)) +
        geom_line() + 
        labs(x="Date",y="Proportion who have received at least one dose",color="Country")
    return(p)
}

plot_vax_cov = function(dt){
    p = ggplot(dt,aes(x=date,y=cum_prop_v,group=age_group,color=age_group)) +
        geom_line() + 
        labs(x="Date",y="Proportion who have received at least one dose",color="Age group") +
        facet_wrap(~country)
    return(p)
}

plot_ifr = function(dt){
    dt[,`:=`(ifr_t_lb=ifr_t-qnorm(0.975)*sigma*ifr_t/ifr,
             ifr_t_ub=ifr_t+qnorm(0.975)*sigma*ifr_t/ifr)]
    p = ggplot(dt,aes(x=date)) + 
        geom_line(aes(y=ifr_t,group=age_group,color=age_group)) +
        geom_ribbon(aes(ymin=ifr_t_lb,ymax=ifr_t_ub,fill=age_group),linetype=0,alpha=0.5) +
        labs(x="Date",y="IFR",color="Age group",fill="Age group") + 
        # scale_y_log10() +
        facet_wrap(~country)
    return(p)    
}

plot_avg_ifr = function(dt){
    dt[,`:=`(ifr_t_lb=ifr_t-qnorm(0.975)*sigma*ifr_t/ifr,
             ifr_t_ub=ifr_t+qnorm(0.975)*sigma*ifr_t/ifr)]
    cols = c("ifr_t","ifr_t_lb","ifr_t_ub")
    p = ggplot(dt[,lapply(.SD,function(x) sum(x*population)/sum(population)),.SDcols=cols,by=.(country,date)],aes(x=date)) + 
        geom_line(aes(y=ifr_t,group=country,color=country)) +
        geom_ribbon(aes(ymin=ifr_t_lb,ymax=ifr_t_ub,fill=country),linetype=0,alpha=0.1) + 
        labs(x="Date",y="IFR",color="Country",fill="Country")
    print(dt[,.(ifr_t=sum(ifr_t*population)/sum(population)),by=.(country,date)][,.(ifr=max(ifr_t)),by=.(country)])
    return(p)
}

convert_to_nat_mean = function(mu,sigma){
    exp(mu+sigma^2/2)
}

convert_to_nat_sd = function(mu,sigma){
    sqrt((exp(sigma^2)-1)*exp(2*mu+sigma^2))
}

# Deconvolution function
deconv = function(dt,dDeath,method = "ride"){
    dt1 = copy(dt)
    
    # Make sure that rows of data table are in the correct order so that 
    # deconvolution output can just be added to table as a column
    # N.B. This is essential!
    setorder(dt1,country,age_group,date)
    
    countries = dt1[,unique(country)]
    
    dt_list = split(dt1[,.(country,age_group,date,deaths_i_both)],by=c("country","age_group"))
    
    mean_dDeath = sum((0:(length(dDeath)-1))*dDeath)
    
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
            control = list(k = 2,eps = c(0.01,2),Tmark = nrow(deaths_i_sts) - mean_dDeath,B=-1,alpha=0.05,eq3a.method = "C")
            deaths_i_stsBP = backprojNP(deaths_i_sts,dDeath,control = control)
            
            # list(deaths_i_stsBP@upperbound,deaths_i_stsBP@lambda)
            exposures_dead_list[[i]] = list(deaths_i_stsBP@upperbound,deaths_i_stsBP@lambda)
        }
        
    } else if (method == "ride"){ # RIDE ALGORITHM
        exposures_dead_list = vector("list",length(dt_list))
        # for (i in 1:length(dt_list)){
        exposures_dead_list = foreach(i=1:length(dt_list)) %dopar% {
            # TODO - check if delay distribution is only defined from day 1 onwards in incidental
            exposures_model = fit_incidence(dt_list[[i]][,as.integer(round(deaths_i_both))],dDeath[2:length(dDeath)]/sum(dDeath[2:length(dDeath)]),dof_grid = seq(6,30,by=2),linear_tail = 28,extrapolation_prior_precision = 100)
            
            list(exposures_model$Ihat,exposures_model$Isamps)
            # exposures_dead_list[[i]] = list(exposures_model$Ihat,exposures_model$Isamps)
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

calc_exposures_and_infections = function(dt,dIncub){
    dt1 = copy(dt)
    
    # Divide by age- and time-dependent IFR to get exposures
    dt1[,exposures := exposures_dead/ifr_t]
    # Make sure data table is in right order for convolution
    setorder(dt1,country,age_group,date)
    # Convolve exposures with incubation period to get infections
    dt1[,infections := disc_conv(exposures,dIncub),by=.(country,age_group)]
    
    return(dt1)
}

calc_cum_exposures_and_infections = function(dt){
    dt1 = copy(dt)
    
    setorder(dt1,country,age_group,date)
    dt1[,`:=`(cum_exp = cumsum(exposures),cum_inf = cumsum(infections)),by=.(country,age_group)]
    dt1[,`:=`(cum_prop_exp = cum_exp/population,cum_prop_inf = cum_inf/population)]
    
    return(dt1)
}


backcalc = function(dt,dDeath,dIncub,method = "ride"){
    dt1 = copy(dt)
    
    # Deconvolve deaths to get IFR-scaled exposures
    out = deconv(dt1,dDeath,method = method)
    backcalc_dt = out[[1]]
    backcalc_samps = out[[2]]
    
    # Calculate exposures and infections
    backcalc_dt = calc_exposures_and_infections(backcalc_dt,dIncub)
    
    return(list(backcalc_dt,backcalc_samps))
}

run_backcalculation = function(dt,dDeath,dIncub,frlty_idx,method = "ride"){
    # Extract non-LTC and LTC data
    dt_non_ltc = dt[,.(country,age_group,date,deaths_i_both=deaths_i_both_non_ltc,ifr,ifr_t,sigma)]
    dt_ltc = dt[get_min_age(age_group)>=60,.(country,age_group,date,deaths_i_both=deaths_i_both_ltc,ifr=frlty_idx*ifr,ifr_t=frlty_idx*ifr_t,sigma=frlty_idx*sigma)]
    
    tstart = Sys.time()
    # non-LTC
    out = backcalc(dt_non_ltc,dDeath,dIncub,method = method)
    backcalc_dt_non_ltc = out[[1]]
    backcalc_samps_non_ltc = out[[2]]
    # LTC
    out = backcalc(dt_ltc,dDeath,dIncub,method = method)
    backcalc_dt_ltc = out[[1]]
    backcalc_samps_ltc = out[[2]]
    tend = Sys.time()
    print(tend-tstart)
    
    # Merge non-LTC and LTC estimates
    backcalc_dt = merge(backcalc_dt_non_ltc,backcalc_dt_ltc,by=c("country","age_group","date"),all.x=T,suffixes=c("_non_ltc","_ltc"))
    setnafill(backcalc_dt,fill=0,cols=c("deaths_i_both_ltc","exposures_dead_ltc","exposures_ltc","infections_ltc"))
    backcalc_dt[,`:=`(deaths_i_both=deaths_i_both_non_ltc+deaths_i_both_ltc,
                  exposures_dead=exposures_dead_non_ltc+exposures_dead_ltc,
                  exposures=exposures_non_ltc+exposures_ltc,
                  infections=infections_non_ltc+infections_ltc)]
    backcalc_dt = merge(backcalc_dt,unique(dt[,.(country,age_group,population)]),by=c("country","age_group"),all.x=T)
    
    # backcalc_samps = mapply("+",backcalc_samps_non_ltc,backcalc_samps_ltc)
    countries = backcalc_dt[,unique(country)]
    agegroups = backcalc_dt[,unique(age_group)]
    agegroups_ltc = backcalc_dt_ltc[,unique(age_group)]
    # backcalc_samps = vector("list",length(backcalc_samps_non_ltc))
    backcalc_samps = backcalc_samps_non_ltc
    for (i in 1:length(countries)){
        for (j in 1:length(agegroups)){
            k = (i-1)*length(agegroups)+j
            # backcalc_samps[[k]] = backcalc_samps_non_ltc[[k]]
            if (agegroups[j] %in% agegroups_ltc){
                backcalc_samps[[k]] = backcalc_samps[[k]] + backcalc_samps_ltc[[(i-1)*length(agegroups_ltc)+which(agegroups_ltc==agegroups[j])]]
            }
        }
    }
    
    # Calculate cumulative infections
    backcalc_dt = calc_cum_exposures_and_infections(backcalc_dt)
    
    return(list(backcalc_dt=backcalc_dt,backcalc_samps=backcalc_samps,
                backcalc_dt_non_ltc=backcalc_dt_non_ltc,backcalc_samps_non_ltc=backcalc_samps_non_ltc,
                backcalc_dt_ltc=backcalc_dt_ltc,backcalc_samps_ltc=backcalc_samps_ltc))
}

calc_exposures_and_infections_CI = function(dt,samps,dIncub){
    # dt[,sigma:=(ifr-ifr_lb)/qnorm(0.975)]
    # dt[,sigma_t:=(1-(1-ei)*cum_prop_v/(1-ei*cum_prop_v)+(1-ei)*cum_prop_v/(1-ei*cum_prop_v)*(1-ed)*(1-em))*sigma]
    countries = dt[,unique(country)]
    agegroups = dt[,unique(age_group)]
    exposures_samps = vector("list",length(samps))
    infections_samps = vector("list",length(samps))
    dt1 = copy(dt)
    for (i in 1:length(countries)){
        print(i)
        for (j in 1:length(agegroups)){
            k = (i-1)*length(agegroups)+j
            cntry = countries[i]
            age_grp = agegroups[j]
            # ifr_vec = dt1[country==cntry & age_group==age_grp,ifr_t]
            # exposures_samps[[k]] = t(t(samps[[k]])/ifr_vec)
            # ifr_samps = matrix(rtruncnorm(nrow(samps[[k]])*ncol(samps[[k]]),a=0,b=Inf,mean=dt1[country==cntry & age_group==age_grp,ifr_t],sd=dt1[country==cntry & age_group==age_grp,sigma]),
            #                    nrow=nrow(samps[[k]]),ncol=ncol(samps[[k]]),byrow=T)
            nsamps = nrow(samps[[k]])
            nt = ncol(samps[[k]])
            ifr_samps = matrix(0,nrow=nsamps,ncol=nt)
            while (sum(ifr_samps<=0)>0){ # resample if there are any negative or zero IFRs
                ifr_samps = matrix(dt1[country==cntry & age_group==age_grp,ifr_t],nrow=nsamps,ncol=nt,byrow=T) + 
                    t(t(matrix(dt1[country==cntry & age_group==age_grp,ifr_t/ifr],
                               nrow=nsamps,ncol=nt,byrow=T)*rnorm(nsamps,0,dt1[country==cntry & age_group==age_grp,sigma][1])))
            }
            exposures_samps[[k]] = samps[[k]]/ifr_samps
            infections_samps[[k]] = apply(exposures_samps[[k]],2,function(x) disc_conv(x,dIncub))
            dt1[country==cntry & age_group==age_grp,`:=`(exposures_dead_q95l=apply(samps[[k]],2,function(x) quantile(x,probs = 0.025)),exposures_dead_q95u=apply(samps[[k]],2,function(x) quantile(x,probs = 0.975)))]
            dt1[country==cntry & age_group==age_grp,`:=`(exposures_q95l=apply(exposures_samps[[k]],2,function(x) quantile(x,probs = 0.025)),exposures_q95u=apply(exposures_samps[[k]],2,function(x) quantile(x,probs = 0.975)))]
            dt1[country==cntry & age_group==age_grp,`:=`(infections_q95l=apply(infections_samps[[k]],2,function(x) quantile(x,probs = 0.025)),infections_q95u=apply(infections_samps[[k]],2,function(x) quantile(x,probs = 0.975)))]
        }
    }
    return(list(dt=dt1,exposures_samps=exposures_samps,infections_samps=infections_samps))
}

plot_infections = function(dt,pop,yl=NULL){
    countries = dt[,unique(country)]
    agegroups = dt[,unique(age_group)]
    min_ages = get_min_age(agegroups)    
    # Add populations if they are not already in data table
    if (!("population" %in% names(dt))){
        base_dt = pop[country %in% countries,.(country,age,population)]
        base_dt[,age_group:=cut(age,c(min_ages,Inf),labels=agegroups,right=F)]
        base_dt = base_dt[,.(population=sum(population)),by=.(country,age_group)]
        dt = merge(dt,base_dt,by=c("country","age_group"),all.x=T)
    }

    # Plot inferred exposures
    p = ggplot(dt,aes(x=date,group=age_group,color=age_group)) +
        geom_line(aes(y=exposures)) +
        geom_ribbon(aes(ymin=exposures_q95l,ymax=exposures_q95u,fill=age_group),linetype=0,alpha=0.5) +
        labs(x="Date",y="Infections",color="Age group",fill="Age group") +
        facet_wrap(~country)
    
    if (!is.null(yl)){
       p = p + coord_cartesian(ylim=c(0,yl)) #dt[!(country %in% c("Denmark","Norway")),max(exposures_q95u)] 
    }

    # # Plot incidence of inferred exposures
    # p = ggplot(dt,aes(x=date,group=age_group,color=age_group)) +
    #     geom_line(aes(y=exposures/population*1e5)) +
    #     geom_ribbon(aes(ymin=exposures_q95l/population*1e5,ymax=exposures_q95u/population*1e5,fill=age_group),linetype=0,alpha=0.5) +
    #     ylim(c(0,dt[country!="Denmark",max(exposures_q95u/population*1e5)])) +
    #     facet_wrap(~country)
    
    # # Plot inferred infections
    # p = ggplot(dt,aes(x=date,group=age_group,color=age_group)) + 
    #     geom_line(aes(y=infections)) +
    #     geom_ribbon(aes(ymin=infections_q95l,ymax=infections_q95u,fill=age_group),linetype=0,alpha=0.5) +
    #     ylim(c(0,dt[!(country %in% c("Denmark","Norway")),max(exposures_q95u)])) +
    #     facet_wrap(~country)
    
    return(p)
}

process_seroprevalence_data = function(sero_data){
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
    # Remove social housing serosurvey data for Denmark
    sero_data = sero_data[!(country=="Denmark" & Central.estimate=="17.30%")]
    # Remove early regional estimates for Italy
    sero_data = sero_data[!(country=="Italy" & Start.date<=as.Date("2020-05-05"))]
    # Remove data for Finland with <20 tests
    sero_data = sero_data[!(country=="Finland" & N.tests<20)]
    
    # Change country name for surveys in England so that England data is plotted
    sero_data[region=="United Kingdom of Great Britain and Northern Ireland; England",country:="England"]
    sero_data[,date := Start.date + (End.date - Start.date)/2]
    
    return(sero_data)
}

pct = function(x) as.numeric(str_replace_all(x, "%", "")) / 100

plot_output = function(fnm,ecdc_cases_by_age,sero_data){
    load(fnm)
    p1_non_ltc = plot_infections(backcalc_dt_non_ltc,pop,5e4)
    # ggsave(paste0(dir_out,"infections_by_age_non_ltc_",method,"_",source_deaths,".png"),width = 10,height = 8)
    p1_ltc = plot_infections(backcalc_dt_ltc,pop,backcalc_dt_ltc[date<max(date)-28,max(exposures_q95u)])
    # ggsave(paste0(dir_out,"infections_by_age_ltc_",method,"_",source_deaths,".png"),width = 10,height = 8)
    p1 = plot_infections(backcalc_dt,pop,5e4)
    # ggsave(paste0(dir_out,"infections_by_age_",method,"_",source_deaths,".png"),width = 10,height = 8)
    
    # Compare with ECDC age-stratified case data
    # ecdc_cases_by_age = fread("../ecdc_data/ecdc_case_data_by_age.csv")
    ecdc_cases_by_age[,date:=ISOweek2date(paste0(sub("-","-W",year_week),"-7"))]
    ecdc_cases_by_age[age_group=="<15yr",age_group:="0-14yr"]
    ggplot(ecdc_cases_by_age[country %in% dt[,unique(country)]],aes(x=date,y=new_cases,group=age_group,color=age_group)) +
        geom_line() +
        facet_wrap(~country)
    ggsave(paste0(dir_out,"ecdc_cases_by_age.png"),width = 10,height = 8)
    ggplot(ecdc_cases_by_age[country_code %in% c("DK","NO")],aes(x=date,y=new_cases,group=age_group,color=age_group)) +
        geom_line() +
        facet_wrap(~country)
    ggsave(paste0(dir_out,"ecdc_cases_by_age_DNK_NOR.png"),width = 6,height = 3)
    
    # Plot estimated infections vs reported cases
    # overall
    ggplot() +
        geom_line(aes(x=date,y=infections),backcalc_dt[,.(infections=sum(infections)),by=.(country,date)]) +
        # geom_ribbon(aes(x=date,ymin=infections_q95l,ymax=infections_q95u),backcalc_dt[,.(infections_q95l=sum(infections_q95l),infections_q95u=sum(infections_q95u)),by=.(country,date)],linetype=0,alpha=0.5) +
        geom_line(aes(x=date,y=new_cases/7),ecdc_cases_by_age[country %in% backcalc_dt[,unique(country)],.(new_cases=sum(new_cases)),by=.(country,date)],linetype="dashed") +
        # coord_cartesian(ylim=c(0,1e5)) +
        facet_wrap(~country)
    ggsave(paste0(dir_out,"infections_vs_obs_cases_",method,"_",source_deaths,".png"),width = 6,height = 4.8)
    
    # by age
    agegroups_comp = c("0-49","50-79","80+")
    min_ages_comp = get_min_age(agegroups_comp)
    
    backcalc_dt[,age_group_comp:=cut(get_min_age(age_group),c(min_ages_comp,Inf),labels=agegroups_comp,right=F)]
    ecdc_cases_by_age[,age_group_comp:=cut(get_min_age(age_group),c(min_ages_comp,Inf),labels=agegroups_comp,right=F)]
    ggplot() +
        geom_line(aes(x=date,y=infections/population,group=age_group_comp,color=age_group_comp),backcalc_dt[,.(population=sum(population),infections=sum(infections)),by=.(country,date,age_group_comp)]) +
        geom_line(aes(x=date,y=(new_cases/7)/population,group=age_group_comp,color=age_group_comp),ecdc_cases_by_age[country %in% backcalc_dt[,unique(country)],.(population=sum(population),new_cases=sum(new_cases)),by=.(country,date,age_group_comp)],linetype="dashed") +
        # geom_line(aes(x=date,y=infections/population,group=age_group_comp,color=age_group_comp),backcalc_dt[,.(population=sum(population),infections=sum(infections)),by=.(country,date,age_group_comp)]) +
        # geom_line(aes(x=date,y=(new_cases/7)/population,group=age_group_comp,color=age_group_comp),ecdc_cases_by_age[country %in% backcalc_dt[,unique(country)],.(population=sum(population),new_cases=sum(new_cases)),by=.(country,date,age_group_comp)],linetype="dashed") +
        labs(x="Date",y="Infections/cases",color="Age group") +
        theme(axis.text.x=element_text(angle=45,hjust=1)) +
        ylim(0,0.007) +
        facet_wrap(~country)
    ggsave(paste0(dir_out,"infections_vs_obs_cases_by_age_",method,"_",source_deaths,".png"),width = 6,height = 4.8)
    
    # Denmark and Norway
    ggplot() +
        geom_line(aes(x=date,y=infections,group=age_group_comp,color=age_group_comp),backcalc_dt[country %in% c("Denmark","Norway"),.(infections=sum(infections)),by=.(country,date,age_group_comp)]) +
        geom_line(aes(x=date,y=new_cases/7,group=age_group_comp,color=age_group_comp),ecdc_cases_by_age[country %in% c("Denmark","Norway"),.(new_cases=sum(new_cases)),by=.(country,date,age_group_comp)],linetype="dashed") +
        facet_wrap(~country)
    ggsave(paste0(dir_out,"infections_vs_obs_cases_by_age_",method,"_",source_deaths,"_DNK_NOR.png"),width = 6,height = 3)
    
    # Compare with SeroTracker seroprevalence data
    # Plot cumulative proportion exposed vs seroprevalence
    # print(merge(backcalc_dt[,.(date=max(date)),by=.(country)],backcalc_dt,by=c("country","date"),all.x=T)[,.(cum_prop_exp=sum(cum_exp)/sum(population)),by=.(country)])
    # cum_exp_samps = vector("list",length(exposures_samps))
    cum_exp_u_samps = vector("list",length(exposures_samps))
    for (i in 1:length(countries)){
        for (j in 1:length(agegroups)){
            k = (i-1)*length(agegroups) + j
            cntry = countries[i]
            age_grp = agegroups[j]
            # cum_exp_samps[[k]] = t(apply(exposures_samps[[k]],1,cumsum))
            # Calculate proportions of infections that are among unvaccinated individuals
            prop_exp_u = dt[country==cntry & age_group==age_grp,1-(1-ei)*cum_prop_v/(1-ei*cum_prop_v)]
            cum_exp_u_samps[[k]] = t(apply(t(t(exposures_samps[[k]]*prop_exp_u)),1,cumsum))
            backcalc_dt[country==cntry & age_group==age_grp,`:=`(cum_exp_u_q95l=apply(cum_exp_u_samps[[k]],2,function(x) quantile(x,probs = 0.025)),
                                                                 cum_exp_u_q95u=apply(cum_exp_u_samps[[k]],2,function(x) quantile(x,probs = 0.975)))] #cum_exp_med=apply(cum_exp_samps[[k]],2,function(x) quantile(x,probs = 0.5)),
        }
    }
    backcalc_dt[,`:=`(cum_prop_exp_u_q95l=cum_exp_u_q95l/population,cum_prop_exp_u_q95u=cum_exp_u_q95u/population)] #cum_prop_exp_med=cum_exp_med/population,
    # backcalc_dt[,`:=`(cum_prop_exp_q95l=pmin(cum_exp_q95l/population,1),cum_prop_exp_q95u=pmin(cum_exp_q95u/population,1))] #cum_prop_exp_med=cum_exp_med/population,
    # TODO - think about interval (quantile/HPDI) to use for uncertainty as posteriors for cumulative proportion infected for some countries and age groups are very skewed
    backcalc_dt[,exposures_u:=dt[,1-(1-ei)*cum_prop_v/(1-ei*cum_prop_v)]*exposures]
    backcalc_dt[,cum_exp_u:=cumsum(exposures_u),by=.(country,age_group)]
    backcalc_dt[,cum_prop_exp_u:=cum_exp_u/population]
    ggplot() +
        geom_line(aes(x=date,y=cum_prop_exp_u,group=age_group,color=age_group),backcalc_dt) +
        geom_ribbon(aes(x=date,ymin=cum_prop_exp_u_q95l,ymax=cum_prop_exp_u_q95u,fill=age_group),backcalc_dt,linetype=0,alpha=0.5) +
        geom_point(aes(x=date,y=pct(Central.estimate),shape=scope),sero_data[country %in% backcalc_dt[,unique(country)]],size=0.8) +
        # geom_errorbar(aes(x=date,ymin=pct(Lower.bound),ymax=pct(Upper.bound)),sero_data[country %in% backcalc_dt[,unique(country)]]) +
        ylim(c(0,1)) +
        labs(x="Date",y="Cumulative proportion unvaccinated infected",color="Age group",fill="Age group",shape="Scope") +
        facet_wrap(~country)
    ggsave(paste0(dir_out,"cum_prop_infected_vs_seroprev_",source_deaths,".png"),width = 10,height = 8)
        
    # Plot cumulative proportion infected or vaccinated vs seroprevalence
    # Add 18 days to exposure dates for infection-to-seroconversion delay: 
    # 13.3 days from onset to seroconversion [Borremans eLife 2020] + 5 days incubation period
    backcalc_dt = merge(backcalc_dt[,date:=date+18],dt[,.(country,date,age_group,cum_prop_v)],by=c("country","date","age_group"))
    p2 = ggplot() +
        geom_line(aes(x=date,y=cum_prop_exp_u+cum_prop_v,group=age_group,color=age_group),backcalc_dt[cum_prop_v==0]) +
        geom_line(aes(x=date,y=cum_prop_exp_u+cum_prop_v,group=age_group,color=age_group),backcalc_dt[cum_prop_v!=0],linetype="dashed") +
        geom_ribbon(aes(x=date,ymin=cum_prop_exp_u_q95l+cum_prop_v,ymax=cum_prop_exp_u_q95u+cum_prop_v,fill=age_group),backcalc_dt,linetype=0,alpha=0.5) + 
        geom_point(aes(x=date,y=pct(Central.estimate),shape=scope),sero_data[country %in% backcalc_dt[,unique(country)]],size=0.8) +
        # geom_errorbar(aes(x=date,ymin=pct(Lower.bound),ymax=pct(Upper.bound)),sero_data[country %in% backcalc_dt[,unique(country)]]) +
        ylim(c(0,1)) +
        labs(x="Date",y="Cumulative proportion infected or vaccinated",color="Age group",fill="Age group",shape="Scope") +
        facet_wrap(~country)
    # ggsave(paste0(dir_out,"cum_prop_infected_or_vaccinated_vs_seroprev_",source_deaths,".png"),width = 10,height = 8)
    p = plot_grid(p1+theme(legend.position="none",axis.text.x=element_text(angle=45,hjust=1)),
                  p2+theme(legend.position="none",axis.text.x=element_text(angle=45,hjust=1)),
                  labels=c("A","B"))
    l = get_legend(p2 + theme(legend.box.margin = margin(0,0,0,12)))
    ggsave(paste0(dir_out,"infections_and_cum_prop_infected_or_vaccinated_",method,"_",source_deaths,".png"),plot_grid(p,l,rel_widths = c(3,.4)),width = 10,height = 4)
}

# calc_init_condns = function(inc_dt,vax_dt_wide){
#     # Assume initially no infected individuals, so numbers in compartments are 
#     # all 0 apart from S
#     # FOR NOW - assume no waning
#     # TODO - include L
#     # TODO - include waning. N.B. will also need to update model for IFR
#     
#     inc_dt1 = copy(inc_dt)
#     
#     # Calculate differences in numbers entering and leaving compartments at each
#     # time point 
#     inc_dt1[,`:=`(diffE=nS_E-nE_I,
#                  diffIp=nE_Ip-nIp_Is,
#                  diffIs=nIp_Is-nIs_R,
#                  diffIa=nE_Ia-nIa_R,
#                  diffR=nIs_R+nIa_R)]
#     # Calculate numbers in each compartment at each time point
#     cols = c("E","Ip","Is","Ia","R")
#     setorder(inc_dt1,country,vrnt,age_group_model,date)
#     inc_dt1[,(cols):=lapply(.SD,cumsum),.SDcols=paste0("diff",cols),by=.(country,vrnt,age_group_model)]
#     # # Calculate compartment prevalences
#     # inc_dt1[,paste0("prev",cols):=lapply(.SD,function(x) x/population),.SDcols=cols]
#     
#     incS_dt = inc_dt1[,.(nS_E=sum(nS_E)),by=.(country,age_group_model,date)]
#     incS_dt = merge(incS_dt,vax_dt_wide[,!c("prop_va","prop_vb","cum_prop_va","cum_prop_vb")],by=c("country","date","age_group_model"))
#     # TODO - update when L is included
#     incS_dt[,`:=`(diffS=-nS_E-nS_Va1-nS_Vb1,
#                   diffVa1=nS_Va1-nVa1_Va2,
#                   diffVa2=nVa1_Va2,
#                   diffVb1=nS_Vb1-nVb1_Vb2,
#                   diffVb2=nVb1_Vb2)]
#     setorder(incS_dt,country,age_group_model,date)
#     cols1 = c("S","Va1","Va2","Vb1","Vb2")
#     incS_dt[,(cols1):=lapply(.SD,cumsum),.SDcols=paste0("diff",cols1),by=.(country,age_group_model)]
#     # N.B. S in line above is still cumulative number of susceptibles that 
#     # have lost susceptibility so need to add number initially susceptible 
#     # (population)
#     incS_dt[,S:=pmax(0,population+S)] # minumum of 0 once everybody has been infected/vaccinated
#     # incS_dt[,paste0("prev",cols1):=lapply(.SD,function(x) x/population),.SDcols=cols1]
#     
#     inc_dt1[,vrnt:=fcase(vrnt=="Other","",
#                         vrnt=="Alpha","2",
#                         vrnt=="Delta","3")]
#     inc_dt_wide = dcast(inc_dt1,country+age_group_model+date+population ~ vrnt,sep = "",value.var = cols)
#     
#     # Extract current prevalences
#     cols2 = c("country","age_group_model","date",cols1)
#     prev_dt = merge(inc_dt_wide,incS_dt[,..cols2],by=c("country","age_group_model","date"))
#     
#     # Calculate compartment prevalences
#     cols3 = setdiff(names(prev_dt),c("country","age_group_model","date","population"))
#     prev_dt[,paste0("prev",cols3):=lapply(.SD,function(x) x/population),.SDcols=cols3]
#     
#     ggplot(melt(prev_dt[age_group_model=="50-54"],measure.vars=patterns("prevE"),value.name = "prevE"),aes(x=date,y=prevE,color=variable)) + 
#         geom_line() + 
#         facet_wrap(~country)
#     ggplot(melt(prev_dt[age_group_model=="50-54"],measure.vars=patterns("prevR"),value.name = "prevR"),aes(x=date,y=prevR,color=variable)) + 
#         geom_line() + 
#         facet_wrap(~country)
#     
#     return(prev_dt)
# }

calc_init_condns = function(inc_dt,vax_dt_wide,agegroups_model,covy_dt,vrnt_prop,ve_params,dE,dIp,dIs,dIa){
    prev_dt = copy(inc_dt)
    
    countries = prev_dt[,unique(country)]
    ncountries = length(countries)
    
    # Cast variant proportion data table to wide format
    vrnt_prop_wide = dcast(vrnt_prop,country+date~vrnt,value.var = "prop_vrnt")
    cols = paste0("prop_vrnt",c("","2","3"))
    setnames(vrnt_prop_wide,c("Other","Alpha","Delta"),cols)
    
    # Merge with overall data table
    prev_dt = merge(prev_dt,vrnt_prop_wide,by=c("country","date"),all.x=T)
    # Backfill variant proportions with first non-NA observation
    prev_dt[,(cols):=nafill(.SD,type="nocb"),.SDcols=cols,by=.(country)]
    
    # Merge with vaccination data
    prev_dt = merge(prev_dt,vax_dt_wide[,!"population"],by=c("country","age_group_model","date"),all.x=T)
    
    # Calculate relative proportions of each vaccine type
    prev_dt[,`:=`(p_va1=pmax((cum_prop_va_1-cum_prop_va_2)/(cum_prop_va_1+cum_prop_vb_1),0),
             p_va2=pmin(cum_prop_va_2/(cum_prop_va_1+cum_prop_vb_1),1),
             p_vb1=pmax((cum_prop_vb_1-cum_prop_vb_2)/(cum_prop_va_1+cum_prop_vb_1),0),
             p_vb2=pmin(cum_prop_vb_2/(cum_prop_va_1+cum_prop_vb_1),1))]
    
    # Calculate average vaccine efficacy over time against 3 variants, 
    # accounting for changing vaccine type proportions and variant proportions
    prev_dt[,`:=`(ei=fifelse(!(cum_prop_va_1==0 & cum_prop_vb_1==0),
                            p_va1*(prop_vrnt*ve_params$ei_va1+prop_vrnt2*ve_params$ei2_va1+prop_vrnt3*ve_params$ei3_va1) +
                            p_va2*(prop_vrnt*ve_params$ei_va2+prop_vrnt2*ve_params$ei2_va2+prop_vrnt3*ve_params$ei3_va2) +
                            p_vb1*(prop_vrnt*ve_params$ei_vb1+prop_vrnt2*ve_params$ei2_vb1+prop_vrnt3*ve_params$ei3_vb1) +
                            p_vb2*(prop_vrnt*ve_params$ei_vb2+prop_vrnt2*ve_params$ei2_vb2+prop_vrnt3*ve_params$ei3_vb2),0),
             ed=fifelse(!(cum_prop_va_1==0 & cum_prop_vb_1==0),
                            p_va1*(prop_vrnt*ve_params$ed_va1i+prop_vrnt2*ve_params$ed_va1i2+prop_vrnt3*ve_params$ed_va1i3) +
                            p_va2*(prop_vrnt*ve_params$ed_va2i+prop_vrnt2*ve_params$ed_va2i2+prop_vrnt3*ve_params$ed_va2i3) +
                            p_vb1*(prop_vrnt*ve_params$ed_vb1i+prop_vrnt2*ve_params$ed_vb1i2+prop_vrnt3*ve_params$ed_vb1i3) +
                            p_vb2*(prop_vrnt*ve_params$ed_vb2i+prop_vrnt2*ve_params$ed_vb2i2+prop_vrnt3*ve_params$ed_vb2i3),0),
             eh=fifelse(!(cum_prop_va_1==0 & cum_prop_vb_1==0),
                            p_va1*(prop_vrnt*ve_params$eh_va1d+prop_vrnt2*ve_params$eh_va1d2+prop_vrnt3*ve_params$eh_va1d3) +
                            p_va2*(prop_vrnt*ve_params$eh_va2d+prop_vrnt2*ve_params$eh_va2d2+prop_vrnt3*ve_params$eh_va2d3) +
                            p_vb1*(prop_vrnt*ve_params$eh_vb1d+prop_vrnt2*ve_params$eh_vb1d2+prop_vrnt3*ve_params$eh_vb1d3) +
                            p_vb2*(prop_vrnt*ve_params$eh_vb2d+prop_vrnt2*ve_params$eh_vb2d2+prop_vrnt3*ve_params$eh_vb2d3),0),
             em=fifelse(!(cum_prop_va_1==0 & cum_prop_vb_1==0),
                            p_va1*(prop_vrnt*ve_params$em_va1d+prop_vrnt2*ve_params$em_va1d2+prop_vrnt3*ve_params$em_va1d3) +
                            p_va2*(prop_vrnt*ve_params$em_va2d+prop_vrnt2*ve_params$em_va2d2+prop_vrnt3*ve_params$em_va2d3) +
                            p_vb1*(prop_vrnt*ve_params$em_vb1d+prop_vrnt2*ve_params$em_vb1d2+prop_vrnt3*ve_params$em_vb1d3) +
                            p_vb2*(prop_vrnt*ve_params$em_vb2d+prop_vrnt2*ve_params$em_vb2d2+prop_vrnt3*ve_params$em_vb2d3),0),
             et=fifelse(!(cum_prop_va_1==0 & cum_prop_vb_1==0),
                            p_va1*(prop_vrnt*ve_params$et_va1i+prop_vrnt2*ve_params$et_va1i2+prop_vrnt3*ve_params$et_va1i3) +
                            p_va2*(prop_vrnt*ve_params$et_va2i+prop_vrnt2*ve_params$et_va2i2+prop_vrnt3*ve_params$et_va2i3) +
                            p_vb1*(prop_vrnt*ve_params$et_vb1i+prop_vrnt2*ve_params$et_vb1i2+prop_vrnt3*ve_params$et_vb1i3) +
                            p_vb2*(prop_vrnt*ve_params$et_vb2i+prop_vrnt2*ve_params$et_vb2i2+prop_vrnt3*ve_params$et_vb2i3),0))]
    
    # TODO - THIS LINE NEEDS CORRECTING ONCE ISSUE WITH et AND ed VALUES IS RESOLVED
    # prev_dt[,`:=`(q=(1-ei)*(1-ed),r=(1-ei)*ed)]
    # prev_dt[,`:=`(q=(1-ei)*(1-ed),r=(1-ei)*(ed-et),s=(1-ei)*et)]
    
    # Calculate cumulative proportion vaccinated with either type A or type B vaccine
    prev_dt[,cum_prop_v:=cum_prop_va_1+cum_prop_vb_1]
    
    # Calculate numbers vaccinated
    # TODO - update this to include first doses
    prev_dt[,nS_V:=nS_Va1+nS_Vb1]
    
    # Calculate numbers of susceptibles and vaccinated individuals infected
    prev_dt[,`:=`(nS_E=(1-(1-ei)*cum_prop_v/(1-ei*cum_prop_v))*exposures,
                  nV_E=(1-ei)*(1-ed)*cum_prop_v/(1-ei*cum_prop_v)*exposures,
                  nV_L=(1-ei)*ed*cum_prop_v/(1-ei*cum_prop_v)*exposures)]
    # prev_dt[,`:=`(nS_E=(1-(1-ei)*cum_prop_v/(1-ei*cum_prop_v))*exposures,
    #               nV_E=(1-ei)*(1-ed)*cum_prop_v/(1-ei*cum_prop_v)*exposures,
    #               nV_L=(1-ei)*(ed-et)*cum_prop_v/(1-ei*cum_prop_v)*exposures,
    #               nV_R=(1-ei)*et*cum_prop_v/(1-ei*cum_prop_v)*exposures)]
    
    # Merge with age-dependent symptomatic fraction
    prev_dt = merge(prev_dt,covy_dt,by=c("country","age_group_model"))
    
    # Calculate number leaving susceptible compartment at each time point
    prev_dt[,diffS:=-nS_E-nS_V]
    # Calculate number of susceptibles at each time point ensuring it can't go
    # negative
    setorder(prev_dt,country,age_group_model,date)
    prev_dt[,S:=pmax(population+cumsum(diffS),0),by=.(country,age_group_model)]
    
    # Where susceptibles are fully depleted, set exposures and vaccinations among susceptibles to 0
    # print(prev_dt[S==0,.(nS_E=sum(nS_E)),by=.(country)])
    prev_dt[S==0,`:=`(nS_E=0,nS_V=0)]
    
    # Convolve exposures with latent period distribution to get infections
    prev_dt[,nE_I:=disc_conv(nS_E+nV_E,dE),by=.(country,age_group_model)]
    prev_dt[,nL_Ia:=disc_conv(nV_L,dE),by=.(country,age_group_model)]
    
    # Calculate symptomatic and asymptomatic infections
    prev_dt[,`:=`(nE_Ip=y*nE_I,nE_Ia=(1-y)*nE_I)]
    
    # Convolve infections to get numbers entering different infection states
    prev_dt[,nIp_Is:=disc_conv(nE_Ip,dIp),by=.(country,age_group_model)]
    prev_dt[,nIs_R:=disc_conv(nIp_Is,dIs),by=.(country,age_group_model)]
    prev_dt[,nIa_R:=disc_conv(nE_Ia+nL_Ia,dIa),by=.(country,age_group_model)]
    
    # Calculate differences in numbers entering and leaving compartments at each
    # time point
    prev_dt[,`:=`(diffV=nS_V-nV_E-nV_L,#-nV_R,
                  diffE=nS_E+nV_E-nE_I,
                  diffL=nV_L-nL_Ia,
                  diffIp=nE_Ip-nIp_Is,
                  diffIs=nIp_Is-nIs_R,
                  diffIa=nE_Ia+nL_Ia-nIa_R,
                  diffR=nIs_R+nIa_R)]#+nV_R)]
    
    # Calculate numbers in each compartment at each time point
    cols1 = c("V","E","L","Ip","Is","Ia","R")
    
    # Calculate numbers in each compartment at each time point
    prev_dt[,(cols1):=lapply(.SD,cumsum),.SDcols=paste0("diff",cols1),by=.(country,age_group_model)]
    
    # Remove difference columns
    prev_dt[,(paste0("diff",c("S",cols1))):=NULL]
    
    # # Calculate prevalences
    # prev_dt[,(c("prevS",paste0("prev",cols1))):=lapply(.SD,function(x) x/population),.SDcols=c("S",cols1)]
    # 
    # # Plot to check
    # print(ggplot(prev_dt,aes(x=date,y=prevS,color=age_group_model)) +
    #     geom_line() +
    #     facet_wrap(~country))
    
    # Calculate cumulative proportion exposed
    prev_dt[,cum_prop_exp:=cumsum(exposures)/population,by=.(country,age_group_model)]
    
    # Extract current prevalence
    max_dates = prev_dt[,.(date=max(date)),by=.(country)]
    prev_dt = merge(max_dates,prev_dt,by=c("country","date"),all.x=T)
    
    return(prev_dt)
}

calc_rem_burden = function(prev_dt,ihr,ifr_dt){
    # max_dates = prev_dt[,.(date=max(date)),by=.(country)]
    # 
    # rem_burden_dt = merge(max_dates,prev_dt,by=c("country","date"),all.x=T)
    
    # Merge with IHR and IFR data tables
    # rem_burden_dt = merge(rem_burden_dt,ihr,by="age_group_model")
    rem_burden_dt = merge(prev_dt,ihr,by="age_group_model")
    rem_burden_dt = merge(rem_burden_dt,ifr_dt,by=c("country","age_group_model"))
    
    # Calculate maximum number of hospitalisations and deaths among remaining 
    # susceptibles if they were all exposed now
    rem_burden_dt[,`:=`(cum_hosp_u=y*ihr*S,cum_deaths_u=y*ifr*S)]
    # Calculate maximum number of hospitalisations and deaths among vaccinated
    # individuals if they were all exposed now
    # rem_burden_dt[,exposures_v:=fifelse(nS_E!=0,pmin(nV_E/nS_E*S,V),0)]
    # rem_burden_dt[,`:=`(cum_hosp_v=y*(1-eh)*ihr*exposures_v,cum_deaths_v=y*(1-em)*ifr*exposures_v)]
    rem_burden_dt[,`:=`(cum_hosp_v=y*(1-ei)*(1-ed)*(1-eh)*ihr*V,cum_deaths_v=y*(1-ei)*(1-ed)*(1-em)*ifr*V)]
    # Calculate maximum number of hospitalisations and deaths among previously infected
    rem_burden_dt[,`:=`(cum_hosp_i=y*(1-ei)*(1-ed)*(1-eh)*ihr*R,cum_deaths_i=y*(1-ei)*(1-ed)*(1-em)*ifr*R)]
    
    # Calculate total maximum hospitalisations and deaths
    rem_burden_dt[,`:=`(cum_hosp=cum_hosp_u+cum_hosp_v+cum_hosp_i,cum_deaths=cum_deaths_u+cum_deaths_v+cum_deaths_i)]
    # rem_burden_dt[,`:=`(cum_inc_hosp=cum_hosp/population,cum_inc_deaths=cum_deaths/population)]
    cols = c("cum_hosp_u","cum_deaths_u","cum_hosp_v","cum_deaths_v","cum_hosp_i","cum_deaths_i","cum_hosp","cum_deaths")
    rem_burden_dt[,(sub("cum","cum_inc",cols)):=lapply(.SD,function(x) x/population),.SDcols=cols]
    
    return(rem_burden_dt)
}
