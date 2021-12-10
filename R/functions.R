# Functions for getting min and max age from age group string
get_min_age = function(x){
    as.numeric(sub("-.*","",sub("\\+|<","-",x)))  
} 

get_max_age = function(x){
    as.numeric(sub(".*-","",sub("\\+|<","-",x)))    
}

get_country_iso_codes = function(source="who"){
    data = get_national_data(source=source)
    setDT(data)
    country_iso_codes = unique(data[!is.na(country),.(country,iso_code)])
    saveRDS(country_iso_codes,"./data/country_iso_codes.rds")
    return(country_iso_codes)
}

# construct a delay distribution following a gamma distribution with mean mu and shape parameter shape.
cm_delay_gamma = function(mu, shape, t_max, t_step)
{
    scale = mu / shape;
    t_points = seq(0, t_max, by = t_step);
    heights = pgamma(t_points + t_step/2, shape, scale = scale) - 
        pgamma(pmax(0, t_points - t_step/2), shape, scale = scale);
    return (data.table(t = t_points, p = heights / sum(heights)))
}

# construct a delay distribution following a lognormal distribution with true mean mu and coefficient of variation cv.
cm_delay_lnorm = function(mu, cv, t_max, t_step)
{
    meanlog = log(mu / sqrt(1 + cv^2));
    sdlog = sqrt(log(1 + cv^2));
    t_points = seq(0, t_max, by = t_step);
    heights = plnorm(t_points + t_step/2, meanlog, sdlog) - 
        plnorm(pmax(0, t_points - t_step/2), meanlog, sdlog);
    return (data.table(t = t_points, p = heights / sum(heights)))
}

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

download_death_data = function(source_deaths="ined"){
    age_strat_data_path = paste0("./data/",source_deaths,"_data/")
    dir.create(age_strat_data_path,recursive = T)
    datapath = function(x) paste0(age_strat_data_path,x)
    if (source_deaths=="coverage"){
        # Download data from OSF
        tbl = osf_retrieve_file("9dsfk")
        osf_download(tbl, path = age_strat_data_path, conflicts = "overwrite")
        out10_tbl = osf_retrieve_file("43ucn")
        osf_download(out10_tbl, path = age_strat_data_path, conflicts = "overwrite")
    } else if (source_deaths=="ined"){
        fnm = "AgeSex.zip"
        download.file("https://www.ined.fr/fichier/rte/166/Page%20Data/Pooled%20Datasets/AgeSex.zip",datapath(fnm))
        unzip(datapath(fnm),exdir = age_strat_data_path)
    }
}

read_death_data = function(source_deaths="coverage"){
    if (source_deaths=="coverage"){
        deaths_raw = fread(cmd = "unzip -cq ./data/coverage_data/Output_10.zip", skip = 3)
    } else if (source_deaths=="ined"){
        deaths_raw = fread("./data/ined_data/AgeSex/Cum_deaths_by_age_sex.csv")
    } else {
        stop("Source of age-stratified death data currently unsupported. Please choose 'coverage' or 'ined'.")
    }
}

get_data = function(url,ft,dir_out,fnm,sheet=NULL,skip=0){
    if (ft=="csv"){
        x = read.csv(url,na.strings="")
    } else if (ft=="tsv"){
        x = read.delim(url,na.strings="")
    } else if (ft=="xlsx"){
        download.file(url,paste0(dir_out,fnm))
        x = read_xlsx(paste0(dir_out,fnm),sheet=sheet,skip=skip)
    }
    setDT(x)
    dir.create(dir_out,recursive = T)
    if (ft=="csv"){
        write.csv(x,paste0(dir_out,fnm),row.names=F)
    } else if (ft=="tsv"){
        write.table(x,paste0(dir_out,fnm),sep="\t",row.names=F)
    }
    return(x)
}

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

clean_ined_death_data = function(deaths_raw,who_deaths,cols){
    deaths = copy(deaths_raw) 
        
    # # Number of datasets per country
    # deaths[,length(unique(excelsheet)),by=.(country)]
    # # Countries with more than 1 dataset
    # deaths[,length(unique(excelsheet)),by=.(country)][V1>1,country]
    # 
    # # See what the unique age groups are
    # deaths[,unique(age_group)]
    # # by country
    # deaths[,table(country,age_group)]
    
    # Change country name for USA to match WHO deaths data
    deaths[country=="United States of America",country:="United States"]
    
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
            (country=="England & Wales" & grepl("NHS",excelsheet)) |
            (country=="Netherlands" & grepl("RIVM_Death",excelsheet))|
            (country=="Italy" & grepl("ISS_Age&Sex",excelsheet)))]
    
    # Remove what appear to be erroneous entries for Switzerland, USA, France, Italy and Portugal
    deaths = deaths[!((country=="Switzerland" & date %in% as.Date(c("2020-11-08","2020-11-14"))) |
                          (country=="United States" & date %in% as.Date(c("2021-05-02","2021-05-09"))) |
                          (country=="France" & date %in% as.Date(c("2020-12-03","2020-12-07","2020-12-08","2021-06-04"))) |
                          (country=="Italy" & between(date,as.Date("2020-08-11"),as.Date("2020-09-07"))) |
                          (country=="Portugal" & date %in% c(as.Date("2021-01-24"),as.Date("2021-01-25"))))]
    # Remove first day of counts for Italy as the count for individuals 90+ is missing
    deaths = deaths[!(country=="Italy" & date==as.Date("2020-03-09"))]
    
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
    max_dates = deaths[,.(date=max(date)),by=.(country)]
    countries = deaths[,unique(country)]
    total_deaths_list = vector("list",length(countries))
    cols3 = c("country","date","deaths_total")
    for (i in 1:length(countries)){
        cntry = countries[i]
        tmp = who_deaths[country==cntry & (date<min_dates[country==cntry,date] | date>max_dates[country==cntry,date]),..cols3]
        setnames(tmp,"deaths_total","cum_deaths_both")
        tmp[,age_group:="Total"]
        total_deaths_list[[i]] = tmp
    }
    total_deaths = rbindlist(total_deaths_list)
    deaths = rbind(deaths,total_deaths,fill=T)
    
    for (i in 1:length(countries)){
        cntry = countries[i]
        if (cntry %in% who_deaths[,unique(country)]){
            # Fill in early deaths for each country by scaling the cumulative number of 
            # deaths up to the first date with disaggregated age groups for each age 
            # group by the distribution of the overall number of daily deaths in the WHO
            # data over time up to that point
            deaths = split_deaths(deaths,who_deaths,cntry,who_deaths[country==cntry,min(date)],min_dates[country==cntry,date]-1,Inf,cols)
            # Fill in recent deaths for each country by scaling the cumulative number of 
            # deaths up to the last date with disaggregated age groups for each age 
            # group by the distribution of the overall number of daily deaths in the WHO
            # data over time since then
            deaths = split_deaths(deaths,who_deaths,cntry,max_dates[country==cntry,date]+1,who_deaths[country==cntry,max(date)],Inf,cols,direction = "forward")
        }
    }
    
    # Exclude totals
    # TODO - update this to incorporate deaths missing ages
    deaths = deaths[!(age_group %in% c("Total","Total known","Total unknown"))]
    
    # Reorder
    setorder(deaths,country,date,age_group)
}

clean_coverage_death_data = function(deaths_raw,who_deaths,cols){
    deaths = copy(deaths_raw)
    
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
    deaths = deaths[!(country=="Italy" & date %in% as.Date(c("2021-09-19","2021-09-21")))] # values are incorrect (lower than previous date)
    
    # deaths = deaths[!(country=="Italy" & between(date,as.Date("2020-10-06")+1,as.Date("2021-02-18")))] # drop values as they're only weekly
    
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
                          #(country=="Italy" & between(date,as.Date("2020-10-06"),as.Date("2021-02-18"))),..cols3]
                         (country=="Italy" & between(date,as.Date("2021-01-12"),as.Date("2021-02-18"))),..cols3]
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
    # deaths = split_deaths(deaths,who_deaths,"Italy",as.Date("2020-10-06")+1,as.Date("2021-02-18"),Inf,paste0(cols1,"_",c("both","male","female")),direction = "forward")
    deaths = split_deaths(deaths,who_deaths,"Italy",as.Date("2021-01-12"),as.Date("2021-02-18"),Inf,paste0(cols1,"_",c("both","male","female")))
    
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
    vrnt_data1 = vrnt_data1[variant!="UNK"]
    
    # Classify variants into non-Alpha-Delta, Alpha and Delta
    vrnt_data1[,vrnt:="Other"]
    vrnt_data1[variant=="B.1.1.7",vrnt:="Alpha"]
    vrnt_data1[variant %in% c("B.1.617.2","AY.4.2"),vrnt:="Delta"]

    # Remove Delta sequences mis-assigned to 2020 week 42 for Germany
    vrnt_data1[country=="Germany" & year_week %in% c("2020-42","2020-44") & variant=="B.1.617.2",number_detections_variant:=0]
    vrnt_data1[country=="Germany" & year_week %in% c("2020-42","2020-44"),`:=`(number_sequenced=sum(number_detections_variant)+number_sequenced-number_sequenced_known_variant,number_sequenced_known_variant=sum(number_detections_variant)),by=.(year_week)]
    
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
    
    # Drop data for last couple of weeks for Netherlands as number sequenced is very low
    vrnt_data1 = vrnt_data1[!(country=="Netherlands" & number_sequenced_known_variant<20)]
    
    # Drop data for Hungary for recent weeks as high numbers of "Other" sequences seem spurious
    vrnt_data1 = vrnt_data1[!(country=="Hungary" & year_week %in% c("2021-41","2021-43","2021-44","2021-45"))]
    
    return(vrnt_data1)
}

process_cog_variant_data = function(cog_vrnt_data){
    cog_vrnt_data1 = copy(cog_vrnt_data)
    cog_vrnt_data1[,vrnt:="Other"]
    cog_vrnt_data1[Lineage=="B.1.1.7" | grepl("Q.",Lineage),vrnt:="Alpha"]
    cog_vrnt_data1[Lineage=="B.1.617.2" | grepl("AY.",Lineage),vrnt:="Delta"]
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
        theme(axis.text.x=element_text(angle=90,hjust=1)) + 
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
    
    # Store early data for total vaccinations in Romania for use later
    tmp = vax[country=="Romania" & age_group=="ALL" & YearWeekISO<="2021-W21",.(country,country_code2=ReportingCountry,year_week_iso=YearWeekISO,date,age_group,vaccine=Vaccine,first=FirstDose,second=SecondDose)]
    tmp[,type:=fifelse(vaccine %in% c("COM","MOD"),"vb",fifelse(vaccine!="UNK","va","vu"))]
    tmp = tmp[,.(first = sum(first), second = sum(second)), by = setdiff(names(tmp),c("vaccine","first","second"))]
    tmp[,`:=`(p_first=first/sum(first),p_second=second/sum(second)),by=.(type)]
    
    # Drop data disaggregated by age for those <18 for now and data aggregated for 
    # <60 and 60+
    # TODO - revisit this
    # print(vax[age_group=="UNK",sum(FirstDose+SecondDose+UnknownDose),by=.(ReportingCountry)])
    vax = vax[!(age_group %in% c("0-4","5-9","10-14","15-17","<60","60+"))]
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
    
    # Split cumulative vaccinations by age group up to 2021-05-24 in Romania by
    # distribution of total vaccinations over time up to that point
    vaxROM = CJ(year_week_iso=vax[country=="Romania" & year_week_iso<="2021-W21",unique(year_week_iso)],age_group=vax[!(age_group %in% c("ALL","UNK")),unique(age_group)],type=vax[,unique(type)])
    vaxROM = merge(vaxROM,tmp[,!c("age_group","first","second")],by=c("year_week_iso","type"))
    vaxROM = merge(vaxROM,vax[country=="Romania" & year_week_iso=="2021-W21" & age_group!="ALL",.(age_group,type,first,second)],by=c("age_group","type"))
    vaxROM[,`:=`(first=p_first*first,second=p_second*second)]
    vaxROM[,`:=`(p_first=NULL,p_second=NULL)]
    
    # Bind with main data table and overwrite first week that has finer age groups
    # as values are cumulative vaccinations up to that point
    vax = rbind(vax[!(country=="Romania" & year_week_iso=="2021-W21" & age_group %in% vaxROM[,unique(age_group)])],vaxROM[year_week_iso<="2021-W21"])
    
    # Get rollout start dates
    start_dates = vax[sum(first)!=0,.(start_date=min(date)),by=.(country)]
    vax = merge(vax,start_dates,by="country")
    # Create rollout week variable
    vax[,rollout_week:=as.numeric(date-start_date)/7+1]
    
    # Split into all ages and age-stratified data
    vax_all_ages = vax[age_group=="ALL"]
    # vax_all_ages[,`:=`(cum_first=cumsum(first),cum_second=cumsum(second)),by=.(country,type)]
    vax = vax[age_group!="ALL"]
    
    # Calculate mean coverage by age from first week of rollout in each country:
    # Sum doses over vaccine types
    agg_vax = vax[!(country %in% c("Bulgaria","Croatia","Estonia","Latvia","Poland","Romania","Slovakia","Slovenia")),.(first=sum(first),second=sum(second)),by=.(country,rollout_week,age_group)]
    # Get vaccination age groups
    agegroups_vax = vax[,sort(unique(age_group))]
    agegroups_vax = agegroups_vax[agegroups_vax!="UNK"]
    min_ages_vax = get_min_age(agegroups_vax)
    # Add populations to data table
    pop_vax = pop[,.(country,age,population,age_group=cut(age,c(min_ages_vax,Inf),labels=agegroups_vax,right=F))]
    pop_vax = pop_vax[,.(population=sum(population)),by=.(country,age_group)]
    agg_vax = merge(agg_vax,pop_vax,by=c("country","age_group"))
    # # Calculate coverage by country, age and rollout week
    # agg_vax[,`:=`(cum_prop_first=cumsum(first)/population,
    #               cum_prop_second=cumsum(second)/population),by=.(country,age_group)]
    # # Calculate mean coverage by age and rollout week across countries
    # agg_vax = agg_vax[,.(cum_prop_first=mean(cum_prop_first),
    #                      cum_prop_second=mean(cum_prop_second)),by=.(age_group,rollout_week)]
    agg_vax[,`:=`(prop_first=first/population,
                  prop_second=second/population),by=.(country,age_group)]
    # Calculate mean coverage by age and rollout week across countries
    agg_vax = agg_vax[,.(prop_first=mean(prop_first),
                         prop_second=mean(prop_second)),by=.(age_group,rollout_week)]
    
    # Make data table of total vaccinations split by population-weighted mean 
    # coverage by age and rollout week across all countries
    vax1 = merge(agg_vax,vax_all_ages[,!"age_group"],by="rollout_week",allow.cartesian=T)
    vax1 = merge(vax1,pop_vax,by=c("country","age_group"))
    # cols = c("first","second","cum_first","cum_second")
    cols = c("first","second")
    vax1[,(cols):=lapply(.SD,as.numeric),.SDcols=cols]
    # # vax1[,`:=`(cum_first=population*cum_prop_first/sum(population*cum_prop_first)*cum_first,
    # #            cum_second=population*cum_prop_second/sum(population*cum_prop_second)*cum_second),
    # #      by=.(country,rollout_week,type)]
    # vax1[,`:=`(cum_first=cum_prop_first/sum(cum_prop_first)*cum_first,
    #            cum_second=cum_prop_second/sum(cum_prop_second)*cum_second),
    #      by=.(country,rollout_week,type)]
    vax1[,`:=`(first=population*prop_first/sum(population*prop_first)*first,
               second=population*prop_second/sum(population*prop_second)*second),
         by=.(country,rollout_week,type)]
    # vax1[,`:=`(cum_prop_first=NULL,cum_prop_second=NULL)]
    vax1[,`:=`(prop_first=NULL,prop_second=NULL)]
    # # vax1[,`:=`(first=pmax(diff(c(0,cum_first)),0),
    # #            second=pmax(diff(c(0,cum_second)),0)),by=.(country,age_group,type)]
    # vax1[,`:=`(first=diff(c(0,cum_first)),
    #            second=diff(c(0,cum_second))),by=.(country,age_group,type)]
    
    # Exclude Netherlands from vax as data is only for 0-17-year-olds
    vax = vax[country!="Netherlands"]
    # Bind countries without age-stratified data into main data table
    # vax = rbind(vax,vax1[!(country %in% vax[,unique(country)]),!c("population","cum_first","cum_second")])
    vax = rbind(vax,vax1[!(country %in% c(vax[,unique(country)],"Germany")),!"population"])
    vax[,`:=`(start_date=NULL,rollout_week=NULL)]

    # Calculate proportion of vaccine types over all countries over time for each age group and dose
    # vax_type = vax[type!="vu",lapply(.SD,function(x) sum(x,na.rm=T)),.SDcols = c("first","second"),by=.(date,age_group,type)]
    # vax_type[,`:=`(p_first=first/sum(first),p_second=second/sum(second)),by=.(date,age_group)]
    # Calculate proportion of vaccine types in each country over time for each age group and dose
    vax_type = vax[type!="vu",.(type,p_first=first/sum(first,na.rm=T),p_second=second/sum(second,na.rm=T)),by=.(country,date,age_group)]
    # Calculate mean proportion over all countries
    cols1 = c("p_first","p_second")
    vax_type = vax_type[,lapply(.SD,function(x) mean(x,na.rm=T)),.SDcols=cols1,by=.(date,age_group,type)]
    # Normalise so that proportions add up to 1
    vax_type[,(cols1):=lapply(.SD,function(x) x/sum(x,na.rm=T)),.SDcols=cols1,by=.(date,age_group)]
    setnafill(vax_type,fill=0,cols=cols1)
    
    # Split vaccines of unknown type according to average proportion of each type in
    # other countries (should really also be by age group) FOR NOW - can correct with OWID data
    vax_wide = dcast(vax,... ~ type,value.var = c("first","second"))
    cols2 = paste0(rep(c("first","second"),each=3),"_",c("va","vb","vu"))
    vax_wide[,(cols2):=lapply(.SD,as.numeric),.SDcols=cols2]
    
    vax_type_wide = dcast(vax_type,... ~ type,value.var = c("p_first","p_second"))
    vax_wide = merge(vax_wide,vax_type_wide,by=c("date","age_group"),all.x=T)
    
    # vax_wide[!is.na(first_vu),`:=`(first_va=fifelse(is.na(first_va),first_vu*vax_wide[is.na(first_vu),mean(first_va/(first_va+first_vb),na.rm=T)],first_va+first_vu*vax_wide[is.na(first_vu),mean(first_va/(first_va+first_vb),na.rm=T)]),
    #                          first_vb=fifelse(is.na(first_vb),first_vu*vax_wide[is.na(first_vu),mean(first_vb/(first_va+first_vb),na.rm=T)],first_vb+first_vu*vax_wide[is.na(first_vu),mean(first_vb/(first_va+first_vb),na.rm=T)]))]
    # vax_wide[!is.na(second_vu),`:=`(second_va=fifelse(is.na(second_va),second_vu*vax_wide[is.na(second_vu),mean(second_va/(second_va+second_vb),na.rm=T)],second_va+second_vu*vax_wide[is.na(second_vu),mean(second_va/(second_va+second_vb),na.rm=T)]),
    #                          second_vb=fifelse(is.na(second_vb),second_vu*vax_wide[is.na(second_vu),mean(second_vb/(second_va+second_vb),na.rm=T)],second_vb+second_vu*vax_wide[is.na(second_vu),mean(second_vb/(second_va+second_vb),na.rm=T)]))]
    vax_wide[!is.na(first_vu),`:=`(first_va=fifelse(is.na(first_va),p_first_va*first_vu,first_va+p_first_va*first_vu),
                                   first_vb=fifelse(is.na(first_vb),p_first_vb*first_vu,first_vb+p_first_vb*first_vu))]
    vax_wide[!is.na(second_vu),`:=`(second_va=fifelse(is.na(second_va),p_second_va*second_vu,second_va+p_second_va*second_vu),
                                   second_vb=fifelse(is.na(second_vb),p_second_vb*second_vu,second_vb+p_second_vb*second_vu))]
    vax_wide[,`:=`(first_vu=NULL,second_vu=NULL,p_first_va=NULL,p_first_vb=NULL,p_second_va=NULL,p_second_vb=NULL)]
    
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
    
    # return(list(vax=vax,num_type=num_type))
    return(list(vax=vax,vax_type=vax_type,agg_vax=agg_vax))
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
    vaxENG[,country:="England"]
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
    
    # Add ISO week
    vaxENG_long[,year_week_iso:=date2ISOweek(date)]
    
    return(vaxENG_long)
}

process_gov_uk_vaccination_data = function(vaccENG,vax_type){
    # Convert column type
    vaccENG[,date:=as.IDate(date)]
    
    # Extract required columns from raw data and rename
    vaxENG = vaccENG[,.(date,age_group=sub("_","-",age),country=areaName,
                        first=newPeopleVaccinatedFirstDoseByVaccinationDate,
                        second=newPeopleVaccinatedSecondDoseByVaccinationDate)]
    
    # Get vaccination age groups
    agegroups_vax = vax_type[,sort(unique(age_group))]
    agegroups_vax = agegroups_vax[agegroups_vax!="UNK"]
    min_ages_vax = get_min_age(agegroups_vax)
    vaxENG[,age_group_vax:=cut(get_min_age(age_group),c(min_ages_vax,Inf),labels=agegroups_vax,right=F)]
    
    # Cast vaccine type frequency data table to wide format
    vax_type_wide = dcast(vax_type,... ~ type,value.var = c("p_first","p_second"))
    # Merge with vaccination data table
    vaxENG = merge(vaxENG,vax_type_wide,by.x=c("date","age_group_vax"),by.y=c("date","age_group"),all.x=T)
    vaxENG[,age_group_vax:=NULL]
    
    # Remove incorrect early proportions of different vaccine types
    vaxENG[date<as.Date("2020-12-21"),`:=`(p_first_va=NA,p_first_vb=NA,p_second_va=NA,p_second_vb=NA)]
    # Forward fill and back fill missing proportions of vaccine types for each age group with nearest observations
    cols=c("p_first_va","p_first_vb","p_second_va","p_second_vb")
    vaxENG[,(cols):=lapply(.SD,function(x) nafill(x,"locf")),.SDcols=cols,by=.(age_group)]
    vaxENG[,(cols):=lapply(.SD,function(x) nafill(x,"nocb")),.SDcols=cols,by=.(age_group)]
    
    # Calculate number of first and second doses of each vaccine type
    vaxENG[,`:=`(first_va=p_first_va*first,first_vb=p_first_vb*first,
                 second_va=p_second_va*second,second_vb=p_second_vb*second)]
    vaxENG[,(c("first","second",cols)):=NULL]
    
    # Melt to long format
    vaxENG_long = melt(vaxENG,measure.vars = c("first_va","first_vb","second_va","second_vb"),variable.name = "dose",value.name = "count")
    vaxENG_long[,`:=`(type=sub("[a-z]+_","",dose),dose=sub("_[a-z]+","",dose))]
    vaxENG_long[,dose:=fifelse(dose=="first",1,2)]
    
    # Add ISO week
    vaxENG_long[,year_week_iso:=ISOweek(date)]
    
    return(vaxENG_long)
}

process_rki_vaccination_data = function(rki_vax,ecdc_vax,pop,agg_vax){
    vaccDEU = copy(rki_vax)
    
    # Rename variables
    setnames(vaccDEU,c("Impfdatum","LandkreisId_Impfort","Altersgruppe","Impfschutz","Anzahl"),
             c("date","county","age_group","dose","count"))
    
    # Convert date column to date format
    vaccDEU[,date:=as.IDate(date)]
    
    # Sum over counties
    vaccDEU = vaccDEU[,.(count=sum(count)),by=.(date,age_group,dose)]
    # Print out proportion of vaccinations with unknown age group
    print(vaccDEU[,sum(count[age_group=="u"])/sum(count)])
    # Drop vaccinations with unknown age group
    vaccDEU = vaccDEU[age_group!="u"]
    
    # Get age groups in raw vaccination data
    agegroups_vax = vaccDEU[,sort(unique(age_group))]
    min_ages_vax = get_min_age(vaccDEU[,sort(unique(age_group))])

    ###    
    # Calculate relative coverages in finer age groups based on data from RKI COVIMO phone survey
    # https://www.rki.de/DE/Content/InfAZ/N/Neuartiges_Coronavirus/Projekte_RKI/COVIMO_Reports/covimo_studie_bericht_8.pdf?__blob=publicationFile
    agegroups = c("12-17","18-29","30-39","40-49","50-59","60-69","70-79","80+")
    min_ages = get_min_age(agegroups)
    rel_cov = c(100,88.2,81.9,87.8,88.9,93.5,95.5,96.6)
    rel_cov_dt = data.table(age_group=agegroups,rel_cov=rel_cov)
    
    popDEU = pop[country=="Germany"]
    popDEU[,age_group:=cut(age,c(min_ages,Inf),labels=agegroups,right=F)]
    popDEU = popDEU[!is.na(age_group),.(population=sum(population)),by=.(age_group)]
    
    rel_cov_dt = merge(rel_cov_dt,popDEU,by="age_group")
    # Add vaccination age groups in raw data
    rel_cov_dt[,age_group_vax:=cut(get_min_age(age_group),c(min_ages_vax,Inf),labels=agegroups_vax,right=F)]
    # Calculate relative coverages within vaccination age groups in raw data
    rel_cov_dt[,rel_cov:=rel_cov*population/sum(rel_cov*population),by=.(age_group_vax)]
    rel_cov_dt[,`:=`(age_group_vax=NULL,population=NULL)]
    
    # ## 2nd attempt by fitting sigmoid curves to reported coverages in different age groups over time in COVIMO
    # ## Doesn't work as gives too low coverage among those 70-79 and 80+
    # # Read in COVIMO coverage data
    # rel_cov = fread("./data/covimo_rel_cov.csv")
    # # Calculate mid-date of each COVIMO
    # rel_cov[,date:=as.IDate(start_date+(end_date-start_date)/2)]
    # rel_cov_dt = CJ(date=vaccDEU[,unique(date)],age_group=agegroups)
    # rel_cov_dt = merge(rel_cov_dt,rel_cov[,!c("start_date","end_date")],by=c("date","age_group"),all.x=T)
    # # rel_cov_dt[,rel_cov_i:=approx(date,rel_cov,date)$y,by=.(age_group)]
    # 
    # rel_cov_dt[,datenum:=as.numeric(date)]
    # 
    # for (i in 2:length(agegroups)){
    #     age_grp = agegroups[i]
    #     # Fit sigmoidal curve to COVIMO coverages
    #     sig <- nls(rel_cov ~ SSlogis(datenum,Asym,xmid,scal),rel_cov_dt[age_group==age_grp])
    #     rel_cov_dt[age_group==age_grp,rel_cov_i:=predict(sig,newdata = rel_cov_dt[age_group==age_grp,.(datenum)])]
    # }
    # rel_cov_dt[,`:=`(rel_cov=NULL,datenum=NULL)]
    # # ggplot(rel_cov_dt,aes(group=age_group,color=age_group)) + geom_point(aes(x=date,y=rel_cov)) + geom_line(aes(x=date,y=rel_cov_i))
    # 
    # rel_cov_dt[,diff_rel_cov:=diff(c(0,rel_cov_i)),by=.(age_group)]
    # rel_cov_dt[,rel_cov_i:=NULL]
    # # ggplot(rel_cov_dt[age_group!="12-17"],aes(group=age_group,color=age_group)) + geom_line(aes(x=date,y=diff_rel_cov))
    # 
    # rel_cov_dt = merge(rel_cov_dt,popDEU,by="age_group")
    # # Add vaccination age groups in raw data
    # rel_cov_dt[,age_group_vax:=cut(get_min_age(age_group),c(min_ages_vax,Inf),labels=agegroups_vax,right=F)]
    # # Calculate relative coverages within vaccination age groups in raw data
    # rel_cov_dt[,diff_rel_cov:=diff_rel_cov*population/sum(diff_rel_cov*population),by=.(age_group_vax,date)]
    # rel_cov_dt[age_group=="12-17",diff_rel_cov:=1]
    # rel_cov_dt[,`:=`(age_group_vax=NULL,population=NULL)]
    # 
    # # ggplot(rel_cov_dt,aes(date,diff_rel_cov,color=age_group)) + geom_line()
    
    # agegroups = agg_vax[,unique(age_group)]
    # min_ages = get_min_age(agegroups)
    ###
    
    # Make data table for all dates and finer age groups and first and second doses
    vaxDEU = CJ(date=vaccDEU[,unique(date)],age_group=agegroups,dose=c(1,2))
    
    ###
    # Add vaccination age groups in raw data
    vaxDEU[,age_group_vax:=cut(get_min_age(age_group),c(min_ages_vax,Inf),labels=agegroups_vax,right=F)]
    
    # vaxDEU[,age_group_vax:=cut(get_min_age(age_group),c(0,min_ages_vax[2:length(min_ages_vax)],Inf),labels=agegroups_vax,right=F)]
    ###
    
    # Merge with raw data
    vaxDEU = merge(vaxDEU,vaccDEU,by.x=c("date","age_group_vax","dose"),by.y=c("date","age_group","dose"),all.x=T)
    setnafill(vaxDEU,fill=0,cols = "count")
    
    ###
    vaxDEU[,age_group_vax:=NULL]
    # Merge with relative coverage data table
    vaxDEU = merge(vaxDEU,rel_cov_dt,by=c("age_group"))
    # vaxDEU = merge(vaxDEU,rel_cov_dt,by=c("age_group","date"))
    # Split infections across age groups according to relative coverage
    vaxDEU[,count:=count*rel_cov]
    vaxDEU[,rel_cov:=NULL]
    # vaxDEU[,count:=count*diff_rel_cov]
    # vaxDEU[,diff_rel_cov:=NULL]
    
    # agg_vax_long = melt(agg_vax,measure.vars = c("prop_first","prop_second"),variable.name = "dose",value.name = "prop")
    # agg_vax_long[,dose:=ifelse(dose=="prop_first",1,2)]
    # 
    # start_date = vaxDEU[,min(date)]
    # vaxDEU[,rollout_week:=floor((date-start_date)/7)+1]
    # vaxDEU = merge(vaxDEU,agg_vax_long,by=c("age_group","rollout_week","dose"),all.x=T)
    # 
    # popDEU = pop[country=="Germany"]
    # popDEU[,age_group:=cut(age,c(min_ages,Inf),labels=agegroups,right=F)]
    # popDEU = popDEU[,.(population=sum(population)),by=.(age_group)]
    # vaxDEU = merge(vaxDEU,popDEU,by="age_group")
    # vaxDEU[,count:=as.numeric(count)]
    # vaxDEU[,count:=prop*population/sum(prop*population)*count,by=.(date,age_group_vax,dose)]
    # setnafill(vaxDEU,fill=0,cols="count")
    # # Drop unwanted variables
    # vaxDEU[,(c("rollout_week","age_group_vax","prop","population")):=NULL]
    ###
    
    # Get proportions of vaccine types from ECDC vaccination data
    vax = copy(ecdc_vax)
    vax = vax[!(TargetGroup %in% c("HCW","LTCF"))]
    vax_typeDEU = vax[ReportingCountry=="DE" & TargetGroup=="ALL",.(year_week_iso=YearWeekISO,vaccine=Vaccine,first=FirstDose,second=SecondDose)] #,population=Denominator
    # Classify vaccine types
    vax_typeDEU[,type:=fifelse(vaccine %in% c("COM","MOD"),"vb",fifelse(vaccine!="UNK","va","vu"))]
    # Sum by vaccine type
    vax_typeDEU = vax_typeDEU[,.(first=sum(first),second=sum(second)),by=.(year_week_iso,type)]
    # Calculate proportions of vaccine types
    vax_typeDEU = vax_typeDEU[type!="vu",.(type,p_first=first/sum(first,na.rm=T),p_second=second/sum(second,na.rm=T)),by=.(year_week_iso)]
    # Add rows for dates missing both vaccine types
    vax_typeDEU = merge(CJ(year_week_iso=vax_typeDEU[,unique(year_week_iso)],type=vax_typeDEU[,unique(type)]),vax_typeDEU,by=c("year_week_iso","type"),all.x=T)
    setnafill(vax_typeDEU,fill=0,cols=c("p_first","p_second"))
    
    # Melt to long format
    vax_typeDEU_long = melt(vax_typeDEU,measure.vars=c("p_first","p_second"),variable.name="dose",value.name="p")
    vax_typeDEU_long[,dose:=ifelse(dose=="p_first",1,2)]
    
    # Merge with vaccinations data table
    vaxDEU[,year_week_iso:=ISOweek(date)]
    vaxDEU = merge(vaxDEU,vax_typeDEU_long,by=c("year_week_iso","dose"),allow.cartesian=T)
    
    # Split vaccinations by type
    vaxDEU[,count:=p*count]
    vaxDEU[,p:=NULL]
    
    # Add country name
    vaxDEU[,country:="Germany"]
    
    # ggplot(vaxDEU,aes(x=date,y=count,group=age_group,color=age_group)) + geom_line() + facet_wrap(~dose)
    return(vaxDEU)
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
    if (vax[,length(unique(country))]==1 && vax[,unique(country)] %in% c("England","Germany")){
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
    # ifr_dt = copy(base_dt)
    ifr_dt = data.table(age=0:100)
    
    names(ifr) = tolower(names(ifr))
    min_ages_ifr = ifr[,get_min_age(age_group)]
    
    # Add IFR age groups
    ifr_dt[,age_group:=cut(age,c(min_ages_ifr,Inf),labels=ifr[,age_group],right=F)]
    
    # Merge with IFR data table
    ifr_dt = merge(ifr_dt,ifr[,.(age_group,ifr=median_perc_mean/100,ifr_lb=ci_95_lb_mean/100,ifr_ub=ci_95_ub_mean/100)],by="age_group")
    
    # Change age groups
    ifr_dt[,age_group := cut(age,c(min_ages,Inf),labels=agegroups,right=F)]
    
    # Calculate population-weighted average for each age group
    cols = c("ifr","ifr_lb","ifr_ub")
    # ifr_dt = ifr_dt[!is.na(age_group),lapply(.SD,function(x) sum(x*population)/sum(population)),.SDcols=cols,by=.(country,age_group)] 
    ifr_dt = ifr_dt[!is.na(age_group),lapply(.SD,mean),.SDcols=cols,by=.(age_group)]
    
    # Calculate standard deviation of posterior distribution of IFR assuming it is normal
    ifr_dt[,sigma:=(ifr-ifr_lb)/qnorm(0.975)]
    
    return(ifr_dt)
}

construct_ihr_data_table = function(ihr,base_dt,min_ages,agegroups){
    # ihr_dt = copy(base_dt)
    ihr_dt = data.table(age=0:100)
    
    min_ages_ihr = ihr[,get_min_age(age_group)]
    
    # Add IFR age groups
    ihr_dt[,age_group:=cut(age,c(min_ages_ihr,Inf),labels=ihr[,age_group],right=F)]
    
    # Merge with IFR data table
    ihr_dt = merge(ihr_dt,ihr[,.(age_group,ihr=median_perc_mean/100,ihr_lb=ci_95_lb_mean/100,ihr_ub=ci_95_ub_mean/100)],by="age_group")
    
    # Change age groups
    ihr_dt[,age_group := cut(age,c(min_ages,Inf),labels=agegroups,right=F)]
    
    # Calculate population-weighted average for each age group
    cols = c("ihr","ihr_lb","ihr_ub")
    # ihr_dt = ihr_dt[!is.na(age_group),lapply(.SD,function(x) sum(x*population)/sum(population)),.SDcols=cols,by=.(country,age_group)]    
    ihr_dt = ihr_dt[!is.na(age_group),lapply(.SD,mean),.SDcols=cols,by=.(age_group)]    
    
    # Calculate standard deviation of posterior distribution of IFR assuming it is normal
    ihr_dt[,sigma:=(ihr-ihr_lb)/qnorm(0.975)]
    
    return(ihr_dt)
}

# construct_data_table = function(agegroups,deaths,pop,cols,ltc_deaths,vax,num_type,ifr){
construct_data_table = function(agegroups,deaths,pop,cols,ltc_deaths,vax,Ab_delay1,Ab_delay2,vaxENG,vaxDEU,ifr,vrnt_prop,ve_params){
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
    
    # Fill forward missing cumulative deaths in age groups under 40 to avoid losing deaths through rounding of fractional deaths generated by small death numbers
    deaths_dt_long[(age_group %in% agegroups[max_ages<40]) | (country=="Netherlands" & age_group %in% agegroups[max_ages<50]),cum_deaths_i:=nafill(cum_deaths,"locf"),by=.(country,sex,age_group)]
    # Interpolate missing cumulative deaths for age groups over 40
    deaths_dt_long[(age_group %in% agegroups[max_ages>=40]) & !(country=="Netherlands" & age_group %in% agegroups[max_ages<50]),cum_deaths_i := approx(date,cum_deaths,date)$y,by=.(country,sex,age_group)]
    # deaths_dt[,(paste0(cols,"_i")) := lapply(.SD,function(yi) approx(date,yi,dates)$y),.SDcols=cols,by=.(country,age_group)]

    # Calculate new daily deaths
    # FOR NOW - as large negative counts have been removed, treat remaining small
    # negative counts from differencing as 0s
    # TODO - sort out whether diff should be padded with 0 (correct only if 
    # there have been no deaths by date of first record) or cum_deaths[1]
    # (correct if death data is complete) or NA
    deaths_dt_long[,deaths_i := c(0,pmax(diff(cum_deaths_i),0)),by=.(country,sex,age_group)]
    
    # Bernoulli resample fractional deaths between 0 and 1 in COVerAGE data for 40-49-year-olds in Slovenia to avoid instability in deconvolution
    set.seed(123)
    deaths_dt_long[country=="Slovenia" & age_group=="40-49" & deaths_i>0 & deaths_i<1,deaths_i:=as.numeric(deaths_i>runif(.N))]

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
    print(ggplot(deaths_dt[!(country %in% c("Ecuador","United States"))],aes(x=date,y=deaths_i_both,group=age_group,color=age_group)) +
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
    
    # Construct vaccination data table for Germany
    vaxDEU_dt = construct_vax_data_table(vaxDEU,dates1,agegroups,pop)
    
    # Bind together
    vax_dt = rbind(vax_dt,vaxENG_dt,vaxDEU_dt)
    
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
    # ggplot(ifr_dt,aes(x=age_group,y=ifr,group=country,color=country)) + 
    #     geom_line() +
    #     # scale_y_log10() +
    #     theme(legend.position = "none")
    ggplot(ifr_dt,aes(x=age_group,y=ifr)) + 
        geom_line() +
        # scale_y_log10() +
        theme(legend.position = "none")
        
    # Merge IFR with death and vaccination data
    # dt = merge(deaths_vax_dt,ifr_dt,by=c("country","age_group"))
    dt = merge(deaths_vax_dt,ifr_dt,by="age_group")
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
    # Calculate relative proportions of individuals who have been vaccinated with 
    # each vaccine type over time
    # cap at 0 and 1 as proportions can't be outside these bounds and for a few 
    # dates in a few places 2nd dose coverage is higher than 1st dose coverage
    dt[,`:=`(p_va1=pmax((cum_prop_va_1-cum_prop_va_2)/(cum_prop_va_1+cum_prop_vb_1),0),
             p_va2=pmin(cum_prop_va_2/(cum_prop_va_1+cum_prop_vb_1),1),
             p_vb1=pmax((cum_prop_vb_1-cum_prop_vb_2)/(cum_prop_va_1+cum_prop_vb_1),0),
             p_vb2=pmin(cum_prop_vb_2/(cum_prop_va_1+cum_prop_vb_1),1))]
    # Normalise to ensure proportions sum to 1
    cols1 = c("p_va1","p_va2","p_vb1","p_vb2")
    dt[,(cols1):=lapply(.SD,function(x) x/(p_va1+p_va2+p_vb1+p_vb2)),.SDcols=cols1]
    
    # Calculate average vaccine efficacy against different outcomes over time, 
    # accounting for changing vaccine type proportions and variant proportions
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
    
    # Calculate cumulative proportion vaccinated
    dt[,cum_prop_v:=pmin(cum_prop_va_1+cum_prop_vb_1,1)] # cap at 1, so IFR can't go negative
    # Calculate IFR over time
    # dt[,ifr_t := ifr]
    dt[,ifr_t:=(1-(1-ei)*cum_prop_v/(1-ei*cum_prop_v)+(1-ei)*cum_prop_v/(1-ei*cum_prop_v)*(1-ed)*(1-em))*ifr]
    # dt[,ifr_t := (1-cum_prop_va-cum_prop_vb)*ifr + cum_prop_va*ifr_va + cum_prop_vb*ifr_vb]
    # dt[,ifr_t := (1-cum_prop_va-cum_prop_vb)*(prop_vrnt+prop_vrnt2+prop_vrnt3)*ifr + 
    #        cum_prop_va*ifr_va + 
    #        cum_prop_vb*ifr_vb]
    
    # # Calculate uncertainty bounds on time-varying IFR
    # dt[,ifr_t_lb:=(1-(1-ei)*cum_prop_v/(1-ei*cum_prop_v)+(1-ei)*cum_prop_v/(1-ei*cum_prop_v)*(1-ed)*(1-em))*ifr_lb]
    # dt[,ifr_t_ub:=(1-(1-ei)*cum_prop_v/(1-ei*cum_prop_v)+(1-ei)*cum_prop_v/(1-ei*cum_prop_v)*(1-ed)*(1-em))*ifr_ub]
    
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

plot_ovrl_vax_cov = function(dt,line_labels=T){
    dt1 = dt[,.(cum_prop_v=sum(cum_prop_v*population)/sum(population)),by=.(country,date)]
    p = ggplot() +
        geom_line(aes(x=date,y=cum_prop_v,group=country,color=country),dt1)
    
    if (line_labels){
        p = p + geom_text_repel(aes(x=date,y=cum_prop_v,label=country),
                          dt1[date==max(date),.(country,date,cum_prop_v)],size=2,hjust=-2,seed=123) +
            labs(x="Date",y="Proportion who have received at least one dose") +
            xlim(NA_Date_,dt[,max(date)+60]) +
            theme(legend.position = "none")
    } else {
        p = p + labs(x="Date",y="Proportion who have received at least one dose",color="Country")
    }

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

plot_avg_ifr = function(dt,line_labels=T){
    dt[,`:=`(ifr_t_lb=ifr_t-qnorm(0.975)*sigma*ifr_t/ifr,
             ifr_t_ub=ifr_t+qnorm(0.975)*sigma*ifr_t/ifr)]
    cols = c("ifr_t","ifr_t_lb","ifr_t_ub")
    dt1 = dt[,lapply(.SD,function(x) sum(x*population)/sum(population)),.SDcols=cols,by=.(country,date)]
    p = ggplot() + 
        geom_line(aes(x=date,y=ifr_t,group=country,color=country),dt1) +
        geom_ribbon(aes(x=date,ymin=ifr_t_lb,ymax=ifr_t_ub,fill=country),dt1,linetype=0,alpha=0.1) +
        ylim(0,NA)
    
    if (line_labels){
        p = p + geom_text_repel(aes(x=date,y=ifr_t,label=country),
                          dt1[date==max(date),.(country,date,ifr_t)],size=2,min.segment.length=0,hjust=-2,seed=123) +
            labs(x="Date",y="IFR") +
            xlim(NA_Date_,dt[,max(date)+60]) +
            theme(legend.position = "none")
    } else {
        p = p + labs(x="Date",y="IFR",color="Country",fill="Country")
    }
    print(dt1[,.(ifr=max(ifr_t)),by=.(country)])
    return(p)
}

convert_to_nat_mean = function(mu,sigma){
    exp(mu+sigma^2/2)
}

convert_to_nat_sd = function(mu,sigma){
    sqrt((exp(sigma^2)-1)*exp(2*mu+sigma^2))
}

# Deconvolution function
deconv = function(dt,dDeath,method = "ride",ip=NULL,odd=NULL){
    dt1 = copy(dt)
    
    # Make sure that rows of data table are in the correct order so that 
    # deconvolution output can just be added to table as a column
    # N.B. This is essential!
    setorder(dt1,country,age_group,date)
    
    countries = dt1[,unique(country)]
    
    dt_list = split(dt1[,.(country,age_group,date,deaths_i_both)],by=c("country","age_group"))
    
    if (method == "backproj"){ # Backprojection
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
        
        mean_dDeath = sum((0:(length(dDeath)-1))*dDeath)
        
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
        
    } else if (method == "ride"){ # RIDE algorithm
        # exposures_dead_list = vector("list",length(dt_list))
        # for (i in 1:length(dt_list)){
        exposures_dead_list = foreach(i=1:length(dt_list)) %dopar% {
            # TODO - check if delay distribution is only defined from day 1 onwards in incidental
            exposures_model = 
                fit_incidence(dt_list[[i]][,as.integer(round(deaths_i_both))],
                              dDeath[2:length(dDeath)]/sum(dDeath[2:length(dDeath)]),
                              dof_grid = seq(6,30,by=2),linear_tail = 28,
                              extrapolation_prior_precision = 100)
            
            # exposures_dead_list[[i]] = list(exposures_model$Ihat,exposures_model$Isamps)
            list(exposures_model$Ihat,exposures_model$Isamps)
        }
        
    } else if (method == "epinow2"){ # EpiNow2 backcalculation
        generation_time = get_generation_time(disease="SARS-CoV-2",source="ganyani")
        exposures_dead_list = vector("list",length(dt_list))
        for (i in 1:length(dt_list)){
        # exposures_dead_list = foreach(i=1:length(dt_list)) %dopar% {
            reported_deaths = dt_list[[i]][,.(date,confirm=as.integer(round(deaths_i_both)))]
            exposures_model = 
                estimate_infections(reported_deaths,
                                    generation_time=generation_time,
                                    delays=delay_opts(ip,odd),rt=NULL,
                                    backcalc=backcalc_opts(),horizon=0,
                                    filter=F,zero_threshold=Inf)
            tmp = exposures_model$summarised[variable=="infections",median]
            exposures_dead_i = c(rep(0,nrow(dt_list[[i]])-length(tmp)),tmp)
            tmp = as.matrix(dcast(exposures_model$samples[parameter=="infections",.(date,sample,value)],sample~date)[,!"sample"])
            exposures_dead_samps_i = cbind(matrix(0,nrow=nrow(tmp),ncol=nrow(dt_list[[i]])-ncol(tmp)),tmp)
            exposures_dead_list[[i]] = list(exposures_dead_i,exposures_dead_samps_i)
            # list(exposures_dead_i,exposures_dead_samps_i)
        }
    } else {
        stop("method must be 'backproj', 'ride', or 'epinow2'.")
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


backcalc = function(dt,dDeath,dIncub,method = "ride",ip=NULL,odd=NULL){
    dt1 = copy(dt)
    
    # Deconvolve deaths to get IFR-scaled exposures
    out = deconv(dt1,dDeath,method = method,ip=ip,odd=odd)
    backcalc_dt = out[[1]]
    backcalc_samps = out[[2]]
    
    # Calculate exposures and infections
    backcalc_dt = calc_exposures_and_infections(backcalc_dt,dIncub)
    
    return(list(backcalc_dt,backcalc_samps))
}

run_backcalculation = function(dt,dDeath,dIncub,frlty_idx,method = "ride",ip=NULL,odd=NULL){
    # Extract non-LTC and LTC data
    dt_non_ltc = dt[,.(country,age_group,date,deaths_i_both=deaths_i_both_non_ltc,ifr,ifr_t,sigma)]
    dt_ltc = dt[get_min_age(age_group)>=60,.(country,age_group,date,deaths_i_both=deaths_i_both_ltc,ifr=frlty_idx*ifr,ifr_t=frlty_idx*ifr_t,sigma=frlty_idx*sigma)]
    
    tstart = Sys.time()
    # non-LTC
    out = backcalc(dt_non_ltc,dDeath,dIncub,method = method,ip=ip,odd=odd)
    backcalc_dt_non_ltc = out[[1]]
    backcalc_samps_non_ltc = out[[2]]
    # LTC
    out = backcalc(dt_ltc,dDeath,dIncub,method = method,ip=ip,odd=odd)
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

process_seroprevalence_data = function(sero_data){
    # sero_data[,country_code:=dt[match(sero_data[,Country],country),country_code]]
    # Filter to only include national and regional surveys with household and community and residual sera samples
    sero_data = sero_data[`Grade of Estimate Scope` %in% c("National","Regional") &
                              `Sample Frame (groups of interest)` %in% c("Household and community samples","Residual sera"),
                          .(#country_code,
                              scope=`Grade of Estimate Scope`,
                              country=Country,
                              state=`State/Province`,
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
    # Drop pre-pandemic sera for Norway
    sero_data = sero_data[!(country=="Norway" & Start.date<as.Date("2020-01-01"))]
    # Drop data from year-long study in Portugal as data in SeroTracker appears to be incorrect and individuals self-reported to lab
    sero_data = sero_data[!(country=="Portugal" & Start.date==as.Date("2020-04-15") & End.date==as.Date("2021-02-15"))]
        
    # Change country name for surveys in England so that England data is plotted
    sero_data[state=="England",country:="England"]
    
    # Calculate mid-date of each survey
    sero_data[,date := Start.date + (End.date - Start.date)/2]

    
    return(sero_data)
}

process_ons_seroprevalence_data = function(ons_sero_data,popENG,agegroups){
    setDT(ons_sero_data)
    ons_sero_data = ons_sero_data[1:48]
    ons_sero_data[,start_date:=as.Date(sub("to .*","",`Weekly period`),"%d %B %Y")]
    ons_sero_data[,end_date:=as.Date(sub(".* to","",`Weekly period`),"%d %B %Y")]
    ons_sero_data_long = melt(ons_sero_data,measure.vars = setdiff(names(ons_sero_data),c("Weekly period","start_date","end_date")))
    ons_sero_data_long[,age_group_num:=as.numeric(sub(".*\\.\\.\\.([0-9]+)","\\1",variable))]
    ons_sero_data_long[,age_group_num:=((age_group_num-2) %/% 5)+1]
    ons_sero_data_long[,variable:=sub("(.*)\\.\\.\\.[0-9]+","\\1",variable)]
    ons_sero_data_long[,age_group:=fcase(age_group_num==1,"16-24",
                             age_group_num==2,"25-34",
                             age_group_num==3,"35-49",
                             age_group_num==4,"50-59",
                             age_group_num==5,"60-64",
                             age_group_num==6,"65-69",
                             age_group_num==7,"70-74",
                             age_group_num==8,"75-79",
                             age_group_num==9,"80+")]
    ons_sero_data_long[,variable:=fcase(grepl("Modelled % testing positive",variable),"perc_pos",
                            grepl("95% Lower credible interval",variable),"ci_95_lb",
                            grepl("95% Upper credible interval",variable),"ci_95_ub",
                            grepl("Number of people testing",variable),"n_pos",
                            grepl("Number of people in",variable),"n")]
    
    # ggplot(ons_sero_data_long[variable=="perc_pos"],aes(start_date,value,color=age_group)) + geom_line()
    
    ons_sero_data_wide = dcast(ons_sero_data_long,start_date+end_date+age_group_num+age_group~variable)
    
    sero_dataENG = popENG[,.(age,population)]
    
    agegroups_ons = ons_sero_data_long[,unique(age_group)]
    min_ages_ons = get_min_age(agegroups_ons)
    
    sero_dataENG[,age_group:=cut(age,c(min_ages_ons,Inf),labels=agegroups_ons,right=F)]
    
    sero_dataENG = merge(sero_dataENG,ons_sero_data_wide[,.(age_group,start_date,end_date,perc_pos,ci_95_lb,ci_95_ub)],by="age_group",allow.cartesian=T)
    
    # Change age groups
    min_ages = get_min_age(agegroups)
    sero_dataENG[,age_group:=cut(age,c(min_ages,Inf),labels=agegroups,right=F)]
    
    pop_age_groups = unique(sero_dataENG[,.(age_group,age,population)])[,.(population=sum(population)),by=.(age_group)]
    sero_dataENG = sero_dataENG[,lapply(.SD,function(x) sum(x*population)/sum(population)),.SDcols=c("perc_pos","ci_95_lb","ci_95_ub"),by=.(start_date,end_date,age_group)]
    
    sero_dataENG[,date:=start_date+(end_date-start_date)/2]
    
    # ggplot(sero_dataENG,aes(start_date,perc_pos,color=age_group)) + geom_line()
    return(sero_dataENG)
}

pct = function(x){
    as.numeric(str_replace_all(x, "%", "")) / 100  
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
        geom_line(aes(y=exposures/1000)) +
        geom_ribbon(aes(ymin=exposures_q95l/1000,ymax=exposures_q95u/1000,fill=age_group),linetype=0,alpha=0.5) +
        labs(x="Date",y="Infections (thousands)",color="Age group",fill="Age group") 
    
    if (length(countries)>1){
        p = p + facet_wrap(~country)
    }
    
    if (!is.null(yl)){
        p = p + coord_cartesian(ylim=c(0,yl/1000)) #dt[!(country %in% c("Denmark","Norway")),max(exposures_q95u)] 
    }
    
    # # Plot incidence of inferred exposures
    # p = ggplot(dt,aes(x=date,group=age_group,color=age_group)) +
    #     geom_line(aes(y=exposures/population*1e5)) +
    #     geom_ribbon(aes(ymin=exposures_q95l/population*1e5,ymax=exposures_q95u/population*1e5,fill=age_group),linetype=0,alpha=0.5) +
    #     ylim(0,dt[country!="Denmark",max(exposures_q95u/population*1e5)]) +
    #     facet_wrap(~country)
    
    # # Plot inferred infections
    # p = ggplot(dt,aes(x=date,group=age_group,color=age_group)) + 
    #     geom_line(aes(y=infections)) +
    #     geom_ribbon(aes(ymin=infections_q95l,ymax=infections_q95u,fill=age_group),linetype=0,alpha=0.5) +
    #     ylim(0,dt[!(country %in% c("Denmark","Norway")),max(exposures_q95u)]) +
    #     facet_wrap(~country)
    
    return(p)
}

plot_output = function(fnm,pop,ecdc_cases_by_age,cases_by_ageENG,popUK,sero_data,source_deaths,method,dir_fig = "./figs/"){
    # Load backcalculation output
    load(fnm)
    
    # Plot inferred exposures
    max_dates = backcalc_dt[,.(date=max(date)),by=.(country)]
    lessthan2m = merge(max_dates,backcalc_dt,by=c("country","date"),all.x=T)[,.(sum(cum_exp)),by=.(country)][V1<2e6,country]
    p1 = plot_infections(backcalc_dt[country %in% lessthan2m],pop,5e3)
    p1a = plot_infections(backcalc_dt[!(country %in% lessthan2m)],pop,5e4)
    # ggsave(paste0(dir_fig,"infections_by_age_",method,"_",source_deaths,".png"),width = 10,height = 8)
    p1_non_ltc = plot_infections(backcalc_dt_non_ltc,pop,5e4)
    # ggsave(paste0(dir_fig,"infections_by_age_non_ltc_",method,"_",source_deaths,".png"),width = 10,height = 8)
    p1_ltc = plot_infections(backcalc_dt_ltc,pop,backcalc_dt_ltc[date<max(date)-28,max(exposures_q95u)])
    # ggsave(paste0(dir_fig,"infections_by_age_ltc_",method,"_",source_deaths,".png"),width = 10,height = 8)
    p = (p1+theme(legend.position="none",axis.text.x=element_text(angle=90,hjust=1)) + 
             p1a+theme(legend.position="none",axis.text.x=element_text(angle=90,hjust=1)))
    
    # Compare with ECDC age-stratified case data
    # Process ECDC case data
    ecdc_cases_by_age1 = copy(ecdc_cases_by_age)
    ecdc_cases_by_age1[,age_group:=sub("yr","",age_group)]
    ecdc_cases_by_age1[age_group=="<15",age_group:="0-14"]
    ecdc_cases_by_age1[,year_week:=sub("-","-W",year_week)]
    
    # Process data for England
    cases_by_ageENG1 = copy(cases_by_ageENG)
    cases_by_ageENG1 = cases_by_ageENG1[!(age %in% c("00_59","60+","unassigned"))]
    cases_by_ageENG1[,age:=sub("_","-",gsub("0([0-9])", "\\1",age))]
    cases_by_ageENG1[,year_week:=ISOweek(date)]
    cases_by_ageENG1 = cases_by_ageENG1[,.(new_cases=sum(cases)),by=.(areaName,age,year_week)]
    cases_by_ageENG1 = merge(cases_by_ageENG1,popUK[name=="England",.(age,population=1000*(f+m))],by="age")
    setnames(cases_by_ageENG1,c("age","areaName"),c("age_group","country"))
    
    # Bind ECDC and England data
    ecdc_cases_by_age1 = rbind(ecdc_cases_by_age1,cases_by_ageENG1,fill=T)
    
    # Convert ISO week to date
    ecdc_cases_by_age1[,date:=ISOweek2date(paste0(year_week,"-7"))]
    
    # Plot ECDC data
    ggplot(ecdc_cases_by_age1[country %in% dt[,unique(country)]],aes(x=date,y=new_cases,group=age_group,color=age_group)) +
        geom_line() +
        labs(x="Date",y="Cases",color="Age group") +
        theme(axis.text.x=element_text(angle=90,hjust=1)) +
        facet_wrap(~country)
    ggsave(paste0(dir_fig,"ecdc_cases_by_age.png"),width = 10,height = 8)
    ggplot(ecdc_cases_by_age1[country_code %in% c("DK","NO")],aes(x=date,y=new_cases,group=age_group,color=age_group)) +
        geom_line() +
        labs(x="Date",y="Cases",color="Age group") +
        theme(axis.text.x=element_text(angle=90,hjust=1)) +
        facet_wrap(~country)
    ggsave(paste0(dir_fig,"ecdc_cases_by_age_DNK_NOR.png"),width = 6,height = 3)
    
    # Plot estimated infections vs reported cases
    # overall
    ggplot() +
        geom_line(aes(x=date,y=infections),backcalc_dt[,.(infections=sum(infections)),by=.(country,date)]) +
        # geom_ribbon(aes(x=date,ymin=infections_q95l,ymax=infections_q95u),backcalc_dt[,.(infections_q95l=sum(infections_q95l),infections_q95u=sum(infections_q95u)),by=.(country,date)],linetype=0,alpha=0.5) +
        geom_line(aes(x=date,y=new_cases/7),ecdc_cases_by_age1[country %in% backcalc_dt[,unique(country)],.(new_cases=sum(new_cases)),by=.(country,date)],linetype="dashed") +
        # coord_cartesian(ylim=c(0,1e5)) +
        labs(x="Date",y="Infections") +
        theme(axis.text.x=element_text(angle=90,hjust=1)) + 
        facet_wrap(~country)
    ggsave(paste0(dir_fig,"infections_vs_obs_cases_",method,"_",source_deaths,".png"),width = 6,height = 4.8)
    
    # by age
    agegroups_comp = c("0-49","50-79","80+")
    min_ages_comp = get_min_age(agegroups_comp)
    
    backcalc_dt[,age_group_comp:=cut(get_min_age(age_group),c(min_ages_comp,Inf),labels=agegroups_comp,right=F)]
    ecdc_cases_by_age1[,age_group_comp:=cut(get_min_age(age_group),c(min_ages_comp,Inf),labels=agegroups_comp,right=F)]
    ggplot() +
        geom_line(aes(x=date,y=1e5*infections/population,group=age_group_comp),backcalc_dt[,.(population=sum(population),infections=sum(infections)),by=.(country,date,age_group_comp)]) +
        geom_line(aes(x=date,y=1e5*(new_cases/7)/population,group=age_group_comp),ecdc_cases_by_age1[country %in% backcalc_dt[,unique(country)],.(population=sum(population),new_cases=sum(new_cases)),by=.(country,date,age_group_comp)],linetype="dashed") +
        # geom_line(aes(x=date,y=infections/population,group=age_group_comp,color=age_group_comp),backcalc_dt[,.(population=sum(population),infections=sum(infections)),by=.(country,date,age_group_comp)]) +
        # geom_line(aes(x=date,y=(new_cases/7)/population,group=age_group_comp,color=age_group_comp),ecdc_cases_by_age1[country %in% backcalc_dt[,unique(country)],.(population=sum(population),new_cases=sum(new_cases)),by=.(country,date,age_group_comp)],linetype="dashed") +
        labs(x="Date",y="Infections/100,000 population") +
        theme(axis.text.x=element_text(angle=90,hjust=1)) +
        ylim(0,650) +
        theme(axis.text.x = element_text(size=10),
              axis.text.y = element_text(size=12),
              axis.title = element_text(size=16),
              strip.text = element_text(size=11)) + 
        facet_grid(age_group_comp~country)
    ggsave(paste0(dir_fig,"infections_vs_obs_cases_by_age_",method,"_",source_deaths,".png"),width = 18,height = 8)
    
    # Denmark and Norway
    ggplot() +
        geom_line(aes(x=date,y=infections,group=age_group_comp,color=age_group_comp),backcalc_dt[country %in% c("Denmark","Norway"),.(infections=sum(infections)),by=.(country,date,age_group_comp)]) +
        geom_line(aes(x=date,y=new_cases/7,group=age_group_comp,color=age_group_comp),ecdc_cases_by_age1[country %in% c("Denmark","Norway"),.(new_cases=sum(new_cases)),by=.(country,date,age_group_comp)],linetype="dashed") +
        labs(x="Date",y="Infections",color="Age group") +
        theme(axis.text.x=element_text(angle=90,hjust=1)) +
        facet_wrap(~country)
    ggsave(paste0(dir_fig,"infections_vs_obs_cases_by_age_",method,"_",source_deaths,"_DNK_NOR.png"),width = 6,height = 3)
    
    return(p)
    
    # # Compare with SeroTracker seroprevalence data
    # # Plot cumulative proportion exposed vs seroprevalence
    # # print(merge(backcalc_dt[,.(date=max(date)),by=.(country)],backcalc_dt,by=c("country","date"),all.x=T)[,.(cum_prop_exp=sum(cum_exp)/sum(population)),by=.(country)])
    # countries = backcalc_dt[,unique(country)]
    # # cum_exp_samps = vector("list",length(exposures_samps))
    # cum_exp_u_samps = vector("list",length(exposures_samps))
    # for (i in 1:length(countries)){
    #     for (j in 1:length(agegroups)){
    #         k = (i-1)*length(agegroups) + j
    #         cntry = countries[i]
    #         age_grp = agegroups[j]
    #         # cum_exp_samps[[k]] = t(apply(exposures_samps[[k]],1,cumsum))
    #         # Calculate proportions of infections that are among unvaccinated individuals
    #         p_exp_u = dt[country==cntry & age_group==age_grp,1-(1-ei)*cum_prop_v/(1-ei*cum_prop_v)]
    #         cum_exp_u_samps[[k]] = t(apply(t(t(exposures_samps[[k]]*p_exp_u)),1,cumsum))
    #         backcalc_dt[country==cntry & age_group==age_grp,`:=`(cum_exp_u_q95l=apply(cum_exp_u_samps[[k]],2,function(x) quantile(x,probs = 0.025)),
    #                                                              cum_exp_u_med=apply(cum_exp_u_samps[[k]],2,function(x) quantile(x,probs = 0.5)),
    #                                                              cum_exp_u_q95u=apply(cum_exp_u_samps[[k]],2,function(x) quantile(x,probs = 0.975)))]
    #     }
    # }
    # backcalc_dt[,`:=`(cum_prop_exp_u_q95l=cum_exp_u_q95l/population,
    #                   cum_prop_exp_u_med=cum_exp_u_med/population,
    #                   cum_prop_exp_u_q95u=cum_exp_u_q95u/population)]
    # # backcalc_dt[,`:=`(cum_prop_exp_q95l=pmin(cum_exp_q95l/population,1),cum_prop_exp_q95u=pmin(cum_exp_q95u/population,1))] #cum_prop_exp_med=cum_exp_med/population,
    # # TODO - think about interval (quantile/HPDI) to use for uncertainty as posteriors for cumulative proportion infected for some countries and age groups are very skewed
    # backcalc_dt[,exposures_u:=dt[,1-(1-ei)*cum_prop_v/(1-ei*cum_prop_v)]*exposures]
    # backcalc_dt[,cum_exp_u:=cumsum(exposures_u),by=.(country,age_group)]
    # backcalc_dt[,cum_prop_exp_u:=cum_exp_u/population]
    # ggplot() +
    #     geom_line(aes(x=date,y=cum_prop_exp_u_med,group=age_group,color=age_group),backcalc_dt) +
    #     geom_ribbon(aes(x=date,ymin=cum_prop_exp_u_q95l,ymax=cum_prop_exp_u_q95u,fill=age_group),backcalc_dt,linetype=0,alpha=0.5) +
    #     geom_point(aes(x=date,y=pct(Central.estimate),shape=scope),sero_data[country %in% backcalc_dt[,unique(country)]],size=0.8) +
    #     # geom_errorbar(aes(x=date,ymin=pct(Lower.bound),ymax=pct(Upper.bound)),sero_data[country %in% backcalc_dt[,unique(country)]]) +
    #     ylim(0,1) +
    #     labs(x="Date",y="Cumulative proportion unvaccinated infected",color="Age group",fill="Age group",shape="Scope") +
    #     facet_wrap(~country)
    # ggsave(paste0(dir_fig,"cum_prop_infected_vs_seroprev_",source_deaths,".png"),width = 10,height = 8)
    #     
    # # Plot cumulative proportion infected or vaccinated vs seroprevalence
    # # Add 18 days to exposure dates for infection-to-seroconversion delay: 
    # # 13.3 days from onset to seroconversion [Borremans eLife 2020] + 5 days incubation period
    # backcalc_dt = merge(backcalc_dt[,date:=date+18],dt[,.(country,date,age_group,cum_prop_v)],by=c("country","date","age_group"))
    # p2 = ggplot() +
    #     geom_line(aes(x=date,y=cum_prop_exp_u_med+cum_prop_v,group=age_group,color=age_group),backcalc_dt[cum_prop_v==0]) +
    #     geom_line(aes(x=date,y=cum_prop_exp_u_med+cum_prop_v,group=age_group,color=age_group),backcalc_dt[cum_prop_v!=0],linetype="dashed") +
    #     geom_ribbon(aes(x=date,ymin=cum_prop_exp_u_q95l+cum_prop_v,ymax=cum_prop_exp_u_q95u+cum_prop_v,fill=age_group),backcalc_dt,linetype=0,alpha=0.5) + 
    #     geom_point(aes(x=date,y=pct(Central.estimate),shape=scope),sero_data[country %in% backcalc_dt[,unique(country)]],size=1) +
    #     # geom_errorbar(aes(x=date,ymin=pct(Lower.bound),ymax=pct(Upper.bound)),sero_data[country %in% backcalc_dt[,unique(country)]]) +
    #     ylim(0,1) +
    #     labs(x="Date",y="Cumulative proportion infected or vaccinated",color="Age group",fill="Age group",shape="Scope") +
    #     facet_wrap(~country)
    # # # ggsave(paste0(dir_fig,"cum_prop_infected_or_vaccinated_vs_seroprev_",source_deaths,".png"),width = 10,height = 8)
    # # p = plot_grid(p1+theme(legend.position="none",axis.text.x=element_text(angle=90,hjust=1)),
    # #               p1a+theme(legend.position="none",axis.text.x=element_text(angle=90,hjust=1)),
    # #               labels=c("A","B"))
    # # l = get_legend(p2) #+ theme(legend.box.margin = margin(0,0,0,12))
    # # ggsave(paste0(dir_fig,"infections_and_cum_prop_infected_or_vaccinated_",method,"_",source_deaths,".png"),plot_grid(p,l,nrow=2,rel_heights = c(1,.1)),width = 10,height = 5)
    # p = (p1+theme(legend.position="none",axis.text.x=element_text(angle=90,hjust=1)) + 
    #          p1a+theme(legend.position="none",axis.text.x=element_text(angle=90,hjust=1))) / 
    #     p2 + theme(legend.position="bottom",axis.text.x=element_text(angle=90,hjust=1)) + plot_layout(heights = c(3,4)) + plot_annotation(tag_levels = 'A')
    # ggsave(paste0(dir_fig,"infections_and_cum_prop_infected_or_vaccinated_",method,"_",source_deaths,".png"),plot=p,width = 10,height = 10)
    # ggsave(paste0(dir_fig,"infections_and_cum_prop_infected_or_vaccinated_",method,"_",source_deaths,".pdf"),plot=p,width = 10,height = 10)
    # 
    # ggplot() +
    #     geom_line(aes(x=date,y=cum_prop_exp_u_med+cum_prop_v,group=age_group,color=age_group),backcalc_dt[country=="England" & cum_prop_v==0]) +
    #     geom_line(aes(x=date,y=cum_prop_exp_u_med+cum_prop_v,group=age_group,color=age_group),backcalc_dt[country=="England" & cum_prop_v!=0]) +
    #     geom_ribbon(aes(x=date,ymin=cum_prop_exp_u_q95l+cum_prop_v,ymax=cum_prop_exp_u_q95u+cum_prop_v,fill=age_group),backcalc_dt[country=="England"],linetype=0,alpha=0.5) + 
    #     geom_point(aes(x=date,y=pct(Central.estimate),shape=scope),sero_data[country=="England"],size=1) +
    #     geom_line(aes(x=start_date+(end_date-start_date)/2,y=perc_pos/100,color=age_group),sero_dataENG[age_group!="0-39"],linetype="longdash") +
    #     ylim(0,1) +
    #     labs(x="Date",y="Cumulative proportion infected or vaccinated",color="Age group",fill="Age group",shape="Scope")
    # ggsave(paste0(dir_fig,"cum_prop_infected_or_vaccinated_vs_ons_seroprev_",method,"_",source_deaths,"_ENG.png"),width = 5,height = 4)
    
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

# calc_curr_prev = function(inc_dt,vax_dt_wide,agegroups_model,covy_dt,vrnt_prop,ve_params,dE,dIp,dIs,dIa,dir_out,i){
#     prev_dt = copy(inc_dt)
#     
#     countries = prev_dt[,unique(country)]
#     ncountries = length(countries)
#     
#     # Cast variant proportion data table to wide format
#     vrnt_prop_wide = dcast(vrnt_prop,country+date~vrnt,value.var = "prop_vrnt")
#     cols = paste0("prop_vrnt",c("","2","3"))
#     setnames(vrnt_prop_wide,c("Other","Alpha","Delta"),cols)
#     
#     # Merge with overall data table
#     prev_dt = merge(prev_dt,vrnt_prop_wide,by=c("country","date"),all.x=T)
#     # Backfill variant proportions with first non-NA observation
#     prev_dt[,(cols):=nafill(.SD,type="nocb"),.SDcols=cols,by=.(country)]
#     
#     # Merge with vaccination data
#     prev_dt = merge(prev_dt,vax_dt_wide[,!"population"],by=c("country","age_group_model","date"),all.x=T)
#     
#     # Calculate relative proportions of individuals who have been vaccinated with 
#     # each vaccine type over time
#     # cap at 0 and 1 as proportions can't be outside these bounds and for a few 
#     # dates in a few places 2nd dose coverage is higher than 1st dose coverage
#     prev_dt[,`:=`(p_va1=pmax((cum_prop_va_1-cum_prop_va_2)/(cum_prop_va_1+cum_prop_vb_1),0),
#              p_va2=pmin(cum_prop_va_2/(cum_prop_va_1+cum_prop_vb_1),1),
#              p_vb1=pmax((cum_prop_vb_1-cum_prop_vb_2)/(cum_prop_va_1+cum_prop_vb_1),0),
#              p_vb2=pmin(cum_prop_vb_2/(cum_prop_va_1+cum_prop_vb_1),1))]
#     # Normalise to ensure proportions sum to 1
#     cols1 = c("p_va1","p_va2","p_vb1","p_vb2")
#     prev_dt[,(cols1):=lapply(.SD,function(x) x/(p_va1+p_va2+p_vb1+p_vb2)),.SDcols=cols1]
#     
#     # Calculate average vaccine efficacy against different outcomes over time, 
#     # accounting for changing vaccine type proportions and variant proportions
#     prev_dt[,`:=`(ei=fifelse(!(cum_prop_va_1==0 & cum_prop_vb_1==0),
#                             p_va1*(prop_vrnt*ve_params$ei_va1+prop_vrnt2*ve_params$ei2_va1+prop_vrnt3*ve_params$ei3_va1) +
#                             p_va2*(prop_vrnt*ve_params$ei_va2+prop_vrnt2*ve_params$ei2_va2+prop_vrnt3*ve_params$ei3_va2) +
#                             p_vb1*(prop_vrnt*ve_params$ei_vb1+prop_vrnt2*ve_params$ei2_vb1+prop_vrnt3*ve_params$ei3_vb1) +
#                             p_vb2*(prop_vrnt*ve_params$ei_vb2+prop_vrnt2*ve_params$ei2_vb2+prop_vrnt3*ve_params$ei3_vb2),0),
#              ed=fifelse(!(cum_prop_va_1==0 & cum_prop_vb_1==0),
#                             p_va1*(prop_vrnt*ve_params$ed_va1i+prop_vrnt2*ve_params$ed_va1i2+prop_vrnt3*ve_params$ed_va1i3) +
#                             p_va2*(prop_vrnt*ve_params$ed_va2i+prop_vrnt2*ve_params$ed_va2i2+prop_vrnt3*ve_params$ed_va2i3) +
#                             p_vb1*(prop_vrnt*ve_params$ed_vb1i+prop_vrnt2*ve_params$ed_vb1i2+prop_vrnt3*ve_params$ed_vb1i3) +
#                             p_vb2*(prop_vrnt*ve_params$ed_vb2i+prop_vrnt2*ve_params$ed_vb2i2+prop_vrnt3*ve_params$ed_vb2i3),0),
#              eh=fifelse(!(cum_prop_va_1==0 & cum_prop_vb_1==0),
#                             p_va1*(prop_vrnt*ve_params$eh_va1d+prop_vrnt2*ve_params$eh_va1d2+prop_vrnt3*ve_params$eh_va1d3) +
#                             p_va2*(prop_vrnt*ve_params$eh_va2d+prop_vrnt2*ve_params$eh_va2d2+prop_vrnt3*ve_params$eh_va2d3) +
#                             p_vb1*(prop_vrnt*ve_params$eh_vb1d+prop_vrnt2*ve_params$eh_vb1d2+prop_vrnt3*ve_params$eh_vb1d3) +
#                             p_vb2*(prop_vrnt*ve_params$eh_vb2d+prop_vrnt2*ve_params$eh_vb2d2+prop_vrnt3*ve_params$eh_vb2d3),0),
#              em=fifelse(!(cum_prop_va_1==0 & cum_prop_vb_1==0),
#                             p_va1*(prop_vrnt*ve_params$em_va1d+prop_vrnt2*ve_params$em_va1d2+prop_vrnt3*ve_params$em_va1d3) +
#                             p_va2*(prop_vrnt*ve_params$em_va2d+prop_vrnt2*ve_params$em_va2d2+prop_vrnt3*ve_params$em_va2d3) +
#                             p_vb1*(prop_vrnt*ve_params$em_vb1d+prop_vrnt2*ve_params$em_vb1d2+prop_vrnt3*ve_params$em_vb1d3) +
#                             p_vb2*(prop_vrnt*ve_params$em_vb2d+prop_vrnt2*ve_params$em_vb2d2+prop_vrnt3*ve_params$em_vb2d3),0),
#              et=fifelse(!(cum_prop_va_1==0 & cum_prop_vb_1==0),
#                             p_va1*(prop_vrnt*ve_params$et_va1i+prop_vrnt2*ve_params$et_va1i2+prop_vrnt3*ve_params$et_va1i3) +
#                             p_va2*(prop_vrnt*ve_params$et_va2i+prop_vrnt2*ve_params$et_va2i2+prop_vrnt3*ve_params$et_va2i3) +
#                             p_vb1*(prop_vrnt*ve_params$et_vb1i+prop_vrnt2*ve_params$et_vb1i2+prop_vrnt3*ve_params$et_vb1i3) +
#                             p_vb2*(prop_vrnt*ve_params$et_vb2i+prop_vrnt2*ve_params$et_vb2i2+prop_vrnt3*ve_params$et_vb2i3),0))]
#     
#     # TODO - THIS LINE NEEDS CORRECTING ONCE ISSUE WITH et AND ed VALUES IS RESOLVED
#     # prev_dt[,`:=`(q=(1-ei)*(1-ed),r=(1-ei)*ed)]
#     # prev_dt[,`:=`(q=(1-ei)*(1-ed),r=(1-ei)*(ed-et),s=(1-ei)*et)]
#     
#     # Calculate cumulative proportion vaccinated with either type A or type B vaccine
#     prev_dt[,cum_prop_v:=pmin(cum_prop_va_1+cum_prop_vb_1,1)] # cap at 1
#     
#     # Calculate numbers vaccinated
#     # TODO - update this to include first doses
#     prev_dt[,nS_V:=nS_Va1+nS_Vb1]
#     
#     # Calculate numbers of susceptibles and vaccinated individuals infected
#     prev_dt[,`:=`(nS_E=(1-(1-ei)*cum_prop_v/(1-ei*cum_prop_v))*pmin(exposures,population),
#                   nV_E=(1-ei)*(1-ed)*cum_prop_v/(1-ei*cum_prop_v)*pmin(exposures,population),
#                   nV_L=(1-ei)*ed*cum_prop_v/(1-ei*cum_prop_v)*pmin(exposures,population))]
#     # prev_dt[,`:=`(nS_E=(1-(1-ei)*cum_prop_v/(1-ei*cum_prop_v))*exposures,
#     #               nV_E=(1-ei)*(1-ed)*cum_prop_v/(1-ei*cum_prop_v)*exposures,
#     #               nV_L=(1-ei)*(ed-et)*cum_prop_v/(1-ei*cum_prop_v)*exposures,
#     #               nV_R=(1-ei)*et*cum_prop_v/(1-ei*cum_prop_v)*exposures)]
#     
#     # Merge with age-dependent symptomatic fraction
#     prev_dt = merge(prev_dt,covy_dt,by="age_group_model")
#     
#     # # Calculate number leaving susceptible compartment at each time point and 
#     # # change in number in vaccinated state
#     # prev_dt[,`:=`(diffS=-nS_E-nS_V,diffV=nS_V-nV_E-nV_L)]
#     # # Calculate number of susceptibles and number in vaccinated state at each 
#     # # time point ensuring it can't go negative
#     setorder(prev_dt,country,age_group_model,date)
#     # prev_dt[,`:=`(S=pmax(population+cumsum(diffS),0),
#     #               V=pmax(cumsum(diffV),0)),by=.(country,age_group_model)]
#     prev_dt[,diffS:=-nS_E-nS_V]
#     prev_dt[,S:=pmax(population+c(0,cumsum(diffS[-.N])),0),by=.(country,age_group_model)]
#     prev_dt[S==0,`:=`(nS_E=0,nS_V=0)]
#     
#     prev_dt[,diffV:=nS_V-nV_E-nV_L]
#     prev_dt[,V:=pmax(c(0,cumsum(diffV[-.N])),0),by=.(country,age_group_model)]
#     prev_dt[V==0,`:=`(nV_E=0,nV_L=0)]
#     
#     # # Where susceptibles or individuals in vaccinated state are fully depleted, 
#     # # set exposures and vaccinations among susceptibles to 0
#     # # print(prev_dt[S==0,.(nS_E=sum(nS_E)),by=.(country)])
#     # prev_dt[S==0,`:=`(nS_E=0,nS_V=0)]
#     # prev_dt[V==0,`:=`(nV_E=0,nV_L=0)]
#     
#     # Convolve exposures with latent period distribution to get infections
#     prev_dt[,nE_I:=disc_conv(nS_E+nV_E,dE),by=.(country,age_group_model)]
#     prev_dt[,nL_Ia:=disc_conv(nV_L,dE),by=.(country,age_group_model)]
#     
#     # Calculate symptomatic and asymptomatic infections
#     prev_dt[,`:=`(nE_Ip=y*nE_I,nE_Ia=(1-y)*nE_I)]
#     
#     # Convolve infections to get numbers entering different infection states
#     prev_dt[,nIp_Is:=disc_conv(nE_Ip,dIp),by=.(country,age_group_model)]
#     prev_dt[,nIs_R:=disc_conv(nIp_Is,dIs),by=.(country,age_group_model)]
#     prev_dt[,nIa_R:=disc_conv(nE_Ia+nL_Ia,dIa),by=.(country,age_group_model)]
#     
#     # Calculate differences in numbers entering and leaving compartments at each
#     # time point
#     prev_dt[,`:=`(#diffV=nS_V-nV_E-nV_L,#-nV_R,
#                   diffE=nS_E+nV_E-nE_I,
#                   diffL=nV_L-nL_Ia,
#                   diffIp=nE_Ip-nIp_Is,
#                   diffIs=nIp_Is-nIs_R,
#                   diffIa=nE_Ia+nL_Ia-nIa_R,
#                   diffR=nIs_R+nIa_R)]#+nV_R)]
#     
#     # Calculate numbers in each compartment at each time point
#     # cols2 = c("V","E","L","Ip","Is","Ia","R")
#     cols2 = c("E","L","Ip","Is","Ia","R")
#     
#     # Calculate numbers in each compartment at each time point
#     # prev_dt[,(cols2):=lapply(.SD,function(x) pmax(cumsum(x),0)),.SDcols=paste0("diff",cols2),by=.(country,age_group_model)] # pmax to stop numbers going negative
#     prev_dt[,(cols2):=lapply(.SD,function(x) c(0,cumsum(x[-.N]))),.SDcols=paste0("diff",cols2),by=.(country,age_group_model)]
#     
#     # Remove difference columns
#     prev_dt[,(paste0("diff",c("S",cols2))):=NULL]
#     
#     # # Calculate prevalences
#     # prev_dt[,(c("prevS",paste0("prev",cols2))):=lapply(.SD,function(x) x/population),.SDcols=c("S",cols2)]
#     # 
#     # # Plot to check
#     # print(ggplot(prev_dt,aes(x=date,y=prevIs,color=age_group_model)) +
#     #     geom_line() +
#     #     facet_wrap(~country))
#     
#     # Calculate cumulative proportion exposed
#     prev_dt[,cum_prop_exp:=cumsum(exposures)/population,by=.(country,age_group_model)]
#     
#     # Save full prevalence data table
#     # qsave(prev_dt[,.(country,age_group_model,date,population,S,V,E,L,Ip,Is,Ia,R)],paste0(dir_out,"prev/prev_",i,".qs"))
#     saveRDS(prev_dt[,.(country,age_group_model,date,population,S,V,E,L,Ip,Is,Ia,R)],paste0(dir_out,"prev/prev_",i,".RDS"))
#     
#     # Extract current prevalence
#     max_dates = prev_dt[,.(date=max(date)),by=.(country)]
#     prev_dt = merge(max_dates,prev_dt,by=c("country","date"),all.x=T)
#     
#     return(prev_dt)
# }

add = function(x,t,n,d){
    idx = 1:min(length(d),length(x)-t)
    x[t+idx-1] = x[t+idx-1] + n*d[idx]
    return(x)
}

solve_diff_eqns = function(exposures,N,nt,q,r,nV,y,dE,dIp,dIs,dIa){
    
    na = length(N)
    
    S = matrix(nrow = na,ncol = nt+1)
    EV = matrix(nrow = na,ncol = nt+1)
    V = matrix(nrow = na,ncol = nt+1)
    E = matrix(nrow = na,ncol = nt+1)
    L = matrix(nrow = na,ncol = nt+1)
    Ia = matrix(nrow = na,ncol = nt+1)
    Ip = matrix(nrow = na,ncol = nt+1)
    Is = matrix(nrow = na,ncol = nt+1)
    R = matrix(nrow = na,ncol = nt+1)
    
    # Set initial conditions
    S[,1] = N
    EV[,1] = nV[,1]
    V[,1] = rep(0,na)
    E[,1] = rep(0,na)
    L[,1] = rep(0,na)
    Ia[,1] = rep(0,na)
    Ip[,1] = rep(0,na)
    Is[,1] = rep(0,na)
    R[,1] = rep(0,na)
    
    # Convolve exposures to get numbers entering each state
    nS_E = matrix(0,nrow = na,ncol = nt)
    nS_V = matrix(0,nrow = na,ncol = nt)
    nV_E = matrix(0,nrow = na,ncol = nt)
    nV_L = matrix(0,nrow = na,ncol = nt)
    # nV_R = matrix(0,nrow = na,ncol = nt)
    nE_I = matrix(0,nrow = na,ncol = nt)
    nE_Ip = matrix(0,nrow = na,ncol = nt)
    nE_Ia = matrix(0,nrow = na,ncol = nt)
    nL_Ia = matrix(0,nrow = na,ncol = nt)
    nIp_Is = matrix(0,nrow = na,ncol = nt)
    nIs_R = matrix(0,nrow = na,ncol = nt)
    nIa_R = matrix(0,nrow = na,ncol = nt)
    # for (a in 1:na){
    #     nE_I[a,] = disc_conv(exposures[a,],dE)
    #     nE_Ip[a,] = y[a]*nE_I[a,]
    #     nE_Ia[a,] = (1-y[a])*nE_I[a,]
    #     nL_Ia[a,] = disc_conv()
    #     nIp_Is[a,] = disc_conv(nE_Ip[a,],dIp)
    #     nIs_R[a,] = disc_conv(nIp_Is[a,],dIs)
    #     nIa_R[a,] = disc_conv(nE_Ia[a,],dIa)
    # }
    
    # Iterate system to update numbers in each state
    for (t in 1:nt){
        for (a in 1:na){
            nS_V[a,t] = min(S[a,t],S[a,t]/(N[a]-EV[a,t])*nV[a,t])
            if (S[a,t]==0){
                nS_E[a,t] = 0
                # nS_V[a,t] = 0
                if (V[a,t]==0){
                    nV_E[a,t] = 0
                    nV_L[a,t] = 0                     
                }
            } else {
                nS_E[a,t] = S[a,t]/(S[a,t]+(q[a,t]+r[a,t])*V[a,t])*exposures[a,t]
                nV_E[a,t] = q[a,t]*V[a,t]/(S[a,t]+(q[a,t]+r[a,t])*V[a,t])*exposures[a,t]
                nV_L[a,t] = r[a,t]*V[a,t]/(S[a,t]+(q[a,t]+r[a,t])*V[a,t])*exposures[a,t]
            }
            nE_I[a,] = add(nE_I[a,],t,nS_E[a,t]+nV_E[a,t],dE)
            nE_Ip[a,] = y[a]*nE_I[a,]
            nE_Ia[a,] = (1-y[a])*nE_I[a,]
            nL_Ia[a,] = add(nL_Ia[a,],t,nV_L[a,t],dE)
            nIp_Is[a,] = add(nIp_Is[a,],t,nE_Ip[a,t],dIp)
            nIs_R[a,] = add(nIs_R[a,],t,nIp_Is[a,t],dIs)
            nIa_R[a,] = add(nIa_R[a,],t,nE_Ia[a,t],dIa)
            
            # S[a,t+1] = max(S[a,t] - nS_E[a,t] - nS_V[a,t],0) # + nV_S[a,t] + nR_S[a,t]
            # V[a,t+1] = V[a,t] + nS_V[a,t] - nV_E[a,t] - nV_L[a,t] #- nV_R[a,t] - nV_S[a,t]
            # E[a,t+1] = E[a,t] + nS_E[a,t] + nV_E[a,t] - nE_Ip[a,t] - nE_Ia[a,t]
            # L[a,t+1] = L[a,t] + nV_L[a,t] - nL_Ia[a,t]
            # Ia[a,t+1] = Ia[a,t] + nE_Ia[a,t] + nL_Ia[a,t] - nIa_R[a,t]
            # Ip[a,t+1] = Ip[a,t] + nE_Ip[a,t] - nIp_Is[a,t]
            # Is[a,t+1] = Is[a,t] + nIp_Is[a,t] - nIs_R[a,t]
            # R[a,t+1] = R[a,t] + nIa_R[a,t] + nIs_R[a,t] #+ nV_R[a,t] - nR_S[a,t]             
        }
        
        S[,t+1] = pmax(S[,t] - nS_E[,t] - nS_V[,t],0) # + nV_S[,t] + nR_S[,t]
        EV[,t+1] = EV[,t] + nV[,t]
        V[,t+1] = V[,t] + nS_V[,t] - nV_E[,t] - nV_L[,t] #- nV_R[,t] #- nV_S[,t]
        E[,t+1] = E[,t] + nS_E[,t] + nV_E[,t] - nE_Ip[,t] - nE_Ia[,t]
        L[,t+1] = L[,t] + nV_L[,t] - nL_Ia[,t]
        Ia[,t+1] = Ia[,t] + nE_Ia[,t] + nL_Ia[,t] - nIa_R[,t]
        Ip[,t+1] = Ip[,t] + nE_Ip[,t] - nIp_Is[,t]
        Is[,t+1] = Is[,t] + nIp_Is[,t] - nIs_R[,t]
        R[,t+1] = R[,t] + nIa_R[,t] + nIs_R[,t] #+ nV_R[,t] #- nR_S[,t]   
    }
    
    return(list(S=S[,1:nt],EV=EV[,1:nt],V=V[,1:nt],E=E[,1:nt],L=L[,1:nt],Ia=Ia[,1:nt],
                Ip=Ip[,1:nt],Is=Is[,1:nt],R=R[,1:nt],nS_E=nS_E,nS_V=nS_V,
                nV_E=nV_E,nV_L=nV_L,nE_I=nE_I,nE_Ip=nE_Ip,nE_Ia=nE_Ia,
                nL_Ia=nL_Ia,nIp_Is=nIp_Is,nIs_R=nIs_R,nIa_R=nIa_R))
    # return(list(S=S[,1:nt],V=V[,1:nt],E=E[,1:nt],L=L[,1:nt],Ia=Ia[,1:nt],
    #              Ip=Ip[,1:nt],Is=Is[,1:nt],R=R[,1:nt],nS_E=nS_E,nS_V=nS_V,
    #              nV_E=nV_E,nV_L=nV_L,nV_R=nV_R,nE_I=nE_I,nE_Ip=nE_Ip,nE_Ia=nE_Ia,
    #              nL_Ia=nL_Ia,nIp_Is=nIp_Is,nIs_R=nIs_R,nIa_R=nIa_R))
}

calc_curr_prev = function(inc_dt,vax_dt_wide,agegroups_model,covy_dt,vrnt_prop,ve_params,dE,dIp,dIs,dIa,dir_out,samp){
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
    
    # Calculate relative proportions of individuals who have been vaccinated with
    # each vaccine type over time
    # cap at 0 and 1 as proportions can't be outside these bounds and for a few 
    # dates in a few places 2nd dose coverage is higher than 1st dose coverage
    prev_dt[,`:=`(p_va1=pmax((cum_prop_va_1-cum_prop_va_2)/(cum_prop_va_1+cum_prop_vb_1),0),
                  p_va2=pmin(cum_prop_va_2/(cum_prop_va_1+cum_prop_vb_1),1),
                  p_vb1=pmax((cum_prop_vb_1-cum_prop_vb_2)/(cum_prop_va_1+cum_prop_vb_1),0),
                  p_vb2=pmin(cum_prop_vb_2/(cum_prop_va_1+cum_prop_vb_1),1))]
    # Normalise to ensure proportions sum to 1
    cols1 = c("p_va1","p_va2","p_vb1","p_vb2")
    prev_dt[,(cols1):=lapply(.SD,function(x) x/(p_va1+p_va2+p_vb1+p_vb2)),.SDcols=cols1]
    
    # Calculate average vaccine efficacy against different outcomes over time, 
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
    # Calculate scaling factors for rates of symptomatic and asymptomatic infection for vaccinated individuals
    prev_dt[,`:=`(q=(1-ei)*(1-ed),r=(1-ei)*ed)]
    
    # Solve difference equations to calculate current prevalences
    prev_list = vector("list",ncountries)
    for (i in 1:ncountries){
        cntry = countries[i]
        exposures = as.matrix(dcast(prev_dt[country==cntry],age_group_model ~ date,value.var = "exposures")[,-1])
        population = unique(prev_dt[country==cntry,.(age_group_model,population)])[,population]
        nt = ncol(exposures)
        nV = as.matrix(dcast(vax_dt_wide[country==cntry,.(date,age_group_model,nV=nS_Va1+nS_Vb1)],age_group_model ~ date,value.var = "nV")[,-1])[,1:nt]
        q = as.matrix(dcast(prev_dt[country==cntry],age_group_model ~ date,value.var = "q")[,-1])
        r = as.matrix(dcast(prev_dt[country==cntry],age_group_model ~ date,value.var = "r")[,-1])
        
        Y = solve_diff_eqns(exposures,population,nt,q,r,nV,covy_dt[,y],dE,dIp,dIs,dIa)
        
        age_grp = agegroups_model[reshape2::melt(Y[[1]])$Var1]
        dts = prev_dt[country==cntry,unique(date)][reshape2::melt(Y[[1]])$Var2]
        pops = population[reshape2::melt(Y[[1]])$Var1]
        for (j in seq_along(Y)){
            Y[[j]] = reshape2::melt(Y[[j]],value.name = names(Y)[j])[names(Y)[j]]
        }
        Y = do.call(cbind,Y)
        setDT(Y)
        Y[,`:=`(country=cntry,age_group_model=age_grp,date=dts,population=pops)]
        
        prev_list[[i]] = Y
    }
    prev_dt1 = rbindlist(prev_list)
    saveRDS(prev_dt1[,.(country,age_group_model,date,population,S,V,E,L,Ip,Is,Ia,R)],paste0(dir_out,"prev/prev_",samp,".RDS"))
    
    # Merge with age-dependent symptomatic fraction
    prev_dt1 = merge(prev_dt1,covy_dt,by="age_group_model")
    # Merge with original prevalence data table
    prev_dt1 = merge(prev_dt1[,!"population"],prev_dt,by=c("country","age_group_model","date"))
    
    # Calculate cumulative proportion exposed
    prev_dt1[,cum_prop_exp:=cumsum(exposures)/population,by=.(country,age_group_model)]
    
    # Calculate cumulative proportion vaccinated with either type A or type B vaccine
    prev_dt1[,cum_prop_v:=pmin(cum_prop_va_1+cum_prop_vb_1,1)] # cap at 1
    
    # Extract current prevalence
    max_dates = prev_dt1[,.(date=max(date)),by=.(country)]
    prev_dt1 = merge(max_dates,prev_dt1,by=c("country","date"),all.x=T)
    
    return(prev_dt1)
}

plot_prevalence = function(dir_out,countries,sero_data,sero_dataENG){
    # Get list of file names for different posterior samples
    fnms = list.files(paste0(dir_out,"prev/"))
    
    ncountries = length(countries)
    
    # # Create list of data tables with samples for each country
    # prev_cntry_list = vector("list",ncountries)
    # for (i in 1:length(fnms)){
    #     tmp = readRDS(paste0(dir_out,"prev/",fnms[i]))
    #     for (j in 1:ncountries){
    #         prev_cntry_list[[j]] = rbind(prev_cntry_list[[j]],tmp[country==countries[j]])
    #     }
    # }
    # # prev_samp_list = foreach(i=1:length(fnms)) %dopar% {
    # #     readRDS(paste0(dir_out,"prev/",fnms[i]))
    # # }
    # # 
    # # prev_cntry_list = foreach(j=1:ncountries) %dopar% {
    # #     tmp = list()
    # #     for (i in 1:length(fnms)){
    # #         tmp = rbind(tmp,prev_samp_list[[i]][country==countries[j]])
    # #     }
    # # }
    #     
    # qprev_list = vector("list",ncountries)
    # for (j in 1:ncountries){
    #     cols = c("S","V","E","L","Ip","Is","Ia","R")
    #     qprev_list[[j]] = prev_cntry_list[[j]][,unlist(lapply(.SD,function(x) list(q95l=quantile(x,probs=0.025),
    #                                                                           med=quantile(x,probs=0.5),
    #                                                                           q95u=quantile(x,probs=0.975))),
    #                                               recursive = F),.SDcols=cols,by=.(country,age_group_model,date,population)]
    # }
    
    # Calculate median and 95% CIs for prevalences for each country
    qprev_list = vector("list",ncountries)
    for (j in 1:ncountries){
        # Extract estimated prevalences for country j from files with posterior samples
        prev_cntry = list()
        for (i in 1:length(fnms)){
            tmp = readRDS(paste0(dir_out,"prev/",fnms[i]))
            tmp = tmp[country==countries[j]]
            prev_cntry = rbind(prev_cntry,tmp)
        }
        
        # Calcuate median and 95% CI
        cols = c("S","V","E","L","Ip","Is","Ia","R")
        qprev_list[[j]] = prev_cntry[,unlist(lapply(.SD,function(x) list(q95l=quantile(x,probs=0.025),
                                                                         med=quantile(x,probs=0.5),
                                                                         q95u=quantile(x,probs=0.975))),
                                             recursive = F),.SDcols=cols,by=.(country,age_group_model,date,population)]
    }
    
    # Bind into one data data table
    qprev_dt = rbindlist(qprev_list)
    
    # Plot
    p1  = ggplot() + 
        geom_line(aes(x=date,y=(V.med+R.med)/population,color=age_group_model),qprev_dt[V.med==0])+
        geom_line(aes(x=date,y=(V.med+R.med)/population,color=age_group_model),qprev_dt[V.med!=0],linetype="longdash")+ 
        geom_ribbon(aes(x=date,ymin=(V.q95l+R.q95l)/population,ymax=(V.q95u+R.q95u)/population,fill=age_group_model),qprev_dt[!(country=="Denmark" & age_group_model=="60-69")],alpha=0.3,linetype=0) +
        geom_point(aes(x=date,y=pct(Central.estimate),shape=scope),sero_data[country %in% countries],size=0.8) +
        labs(x="Date",y="Cumulative proportion infected or vaccinated",fill="Age group",color="Age group",shape="Scope") +
        coord_cartesian(ylim=c(0,1)) +
        facet_wrap(~country)
    
    # Plot estimates for England vs ONS modelled estimates
    p2 = ggplot() +
        geom_line(aes(x=date,y=(V.med+R.med)/population),qprev_dt[country=="England"]) + 
        geom_ribbon(aes(x=date,ymin=(V.q95l+R.q95l)/population,ymax=(V.q95u+R.q95u)/population),qprev_dt[country=="England"],alpha=0.5,linetype=0) +
        geom_line(aes(x=date,y=pct(perc_pos)),sero_dataENG[age_group_model!="0-39"],linetype="longdash") +
        geom_ribbon(aes(x=date,ymin=pct(ci_95_lb),ymax=pct(ci_95_ub)),sero_dataENG[age_group_model!="0-39"],alpha=0.3,linetype=0) +
        labs(x="Date",y="Cumulative proportion infected or vaccinated") +
        theme(axis.text.x = element_text(angle=45,hjust=1)) +
        coord_cartesian(ylim=c(0,1)) +
        facet_wrap(~age_group_model)
    
    return(list(p1=p1,p2=p2,qprev_dt=qprev_dt))
}

calc_rem_burden = function(prev_dt,ihr_dt,ifr_dt,frlty_idx,ve_params,imm_esc=0,svrty=1){
    # Merge with IHR and IFR data tables
    rem_burden_dt = merge(prev_dt,ihr_dt,by="age_group_model")
    # rem_burden_dt = merge(rem_burden_dt,ifr_dt,by="age_group_model")
    
    # Make vector of random draws from truncated normal distribution for IHR
    ihr_vec = numeric(nrow(rem_burden_dt))
    while (sum(ihr_vec<=0)>0){
        ihr_vec = rem_burden_dt[,ihr+rnorm(.N,0,sigma)]    
    }
    rem_burden_dt[,ihr_samp:=ihr_vec]
    
    # Calculate maximum number of hospitalisations and deaths among remaining 
    # susceptibles if they were all exposed now
    rem_burden_dt[,p_exp_ltc:=exposures_ltc/exposures]
    rem_burden_dt[,f:=(1-p_exp_ltc)+p_exp_ltc*frlty_idx]
    rem_burden_dt[,`:=`(cum_hosp_u=f*svrty*ihr_samp*S,
                        cum_deaths_u=f*svrty*ifr*S)]
    # Calculate maximum number of hospitalisations and deaths among vaccinated
    # individuals if they were all exposed now
    rem_burden_dt[,`:=`(cum_hosp_v=f*svrty*(1-(1-imm_esc)*ei)*(1-ed)*(1-eh)*ihr_samp*V,
                        cum_deaths_v=f*svrty*(1-(1-imm_esc)*ei)*(1-ed)*(1-em)*ifr*V)]
    # Calculate maximum number of hospitalisations and deaths among previously infected
    rem_burden_dt[,`:=`(
        ei_i=(prop_vrnt*ve_params$ei_va2+prop_vrnt2*ve_params$ei2_va2+prop_vrnt3*ve_params$ei3_va2+
                  prop_vrnt*ve_params$ei_vb2+prop_vrnt2*ve_params$ei2_vb2+prop_vrnt3*ve_params$ei3_vb2)/2,
        ed_i=(prop_vrnt*ve_params$ed_va2i+prop_vrnt2*ve_params$ed_va2i2+prop_vrnt3*ve_params$ed_va2i3+
                  prop_vrnt*ve_params$ed_vb2i+prop_vrnt2*ve_params$ed_vb2i2+prop_vrnt3*ve_params$ed_vb2i3)/2,
        eh_i=(prop_vrnt*ve_params$eh_va2d+prop_vrnt2*ve_params$eh_va2d2+prop_vrnt3*ve_params$eh_va2d3+
                  prop_vrnt*ve_params$eh_vb2d+prop_vrnt2*ve_params$eh_vb2d2+prop_vrnt3*ve_params$eh_vb2d3)/2,
        em_i=(prop_vrnt*ve_params$em_va2d+prop_vrnt2*ve_params$em_va2d2+prop_vrnt3*ve_params$em_va2d3+
                  prop_vrnt*ve_params$em_vb2d+prop_vrnt2*ve_params$em_vb2d2+prop_vrnt3*ve_params$em_vb2d3)/2)]
    rem_burden_dt[,`:=`(cum_hosp_i=f*svrty*(1-(1-imm_esc)*ei_i)*(1-ed_i)*(1-eh_i)*ihr_samp*R,
                        cum_deaths_i=f*svrty*(1-(1-imm_esc)*ei_i)*(1-ed_i)*(1-em_i)*ifr*R)]
    
    # Calculate total maximum hospitalisations and deaths
    rem_burden_dt[,`:=`(cum_hosp=cum_hosp_u+cum_hosp_v+cum_hosp_i,cum_deaths=cum_deaths_u+cum_deaths_v+cum_deaths_i)]
    cols = c("cum_hosp_u","cum_deaths_u","cum_hosp_v","cum_deaths_v","cum_hosp_i","cum_deaths_i","cum_hosp","cum_deaths")
    rem_burden_dt[,(sub("cum","cum_inc",cols)):=lapply(.SD,function(x) x/population),.SDcols=cols]
    
    return(rem_burden_dt)
}
