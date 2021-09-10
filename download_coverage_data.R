# library(tidyverse)
# library(readr)
library(data.table)
library(lubridate)
library(ggplot2)
library(osfr)
library(covidAgeData)

age_strat_data_path = "../COVerAGE_data/"
datapath = function(x) paste0(age_strat_data_path,x)

# Download input data from OSF
tbl = osf_retrieve_file("9dsfk")
osf_download(tbl, path = age_strat_data_path, conflicts = "overwrite")

# Read in input data
inputDB = fread(cmd = paste("unzip -cq",datapath("inputDB.zip")), skip = 1)
head(inputDB)

# Drop extra columns
inputDB[,`:=`(y=NULL,`2499`=NULL)]

# See what variables there are
inputDB[,unique(Measure)]

# See what countries are covered
inputDB[,unique(Country)]

# See breakdown by sex
inputDB[,table(Country,Sex)]

# See breakdown by variable and country
inputDB[,table(Country,Measure)]

# See breakdown of vaccination data by age
inputDB[Measure=="Vaccinations",table(Country,Age)]

# Convert dates
inputDB[,Date:=as.Date(Date,"%d.%m.%Y")]

## Plot
# Deaths by age group:
# both sexes
ggplot(inputDB[!(Country %in% c("USA","Brazil","India","United Kingdom","England","England and Wales","Italy")) & Measure=="Deaths" & Region=="All" & !(Age %in% c("TOT","UNK","Unknown")) & Sex=="b"],aes(x=Date,y=Value,group=Age,color=Age)) +
    geom_line() +
    facet_wrap(~Country)
# males
ggplot(inputDB[!(Country %in% c("USA","Brazil","India","United Kingdom","England","England and Wales","Italy")) & Measure=="Deaths" & Region=="All" & !(Age %in% c("TOT","UNK","Unknown")) & Sex=="m"],aes(x=Date,y=Value,group=Age,color=Age)) +
    geom_line() +
    facet_wrap(~Country)
# females
ggplot(inputDB[!(Country %in% c("USA","Brazil","India","United Kingdom","England","England and Wales","Italy")) & Measure=="Deaths" & Region=="All" & !(Age %in% c("TOT","UNK","Unknown")) & Sex=="f"],aes(x=Date,y=Value,group=Age,color=Age)) +
    geom_line() +
    facet_wrap(~Country)

# Vaccinations by age group:
# both sexes
ggplot(inputDB[!(Country %in% c("USA","Brazil","India","United Kingdom","England","England and Wales","Italy")) & Measure=="Vaccinations" & Region=="All" & !(Age %in% c("TOT","UNK","Unknown")) & Sex=="b"],aes(x=Date,y=Value,group=Age,color=Age)) +
    geom_line() +
    facet_wrap(~Country)
# males
ggplot(inputDB[!(Country %in% c("USA","Brazil","India","United Kingdom","England","England and Wales","Italy")) & Measure=="Vaccinations" & Region=="All" & !(Age %in% c("TOT","UNK","Unknown")) & Sex=="m"],aes(x=Date,y=Value,group=Age,color=Age)) +
    geom_line() +
    facet_wrap(~Country)
# females
ggplot(inputDB[!(Country %in% c("USA","Brazil","India","United Kingdom","England","England and Wales","Italy")) & Measure=="Vaccinations" & Region=="All" & !(Age %in% c("TOT","UNK","Unknown")) & Sex=="f"],aes(x=Date,y=Value,group=Age,color=Age)) +
    geom_line() +
    facet_wrap(~Country)

# Download output data tables
out10_tbl = osf_retrieve_file("43ucn")
osf_download(out10_tbl, path = age_strat_data_path, conflicts = "overwrite")

# Read in output data
out10 = fread(cmd = paste("unzip -cq",datapath("Output_10.zip")), skip = 3)
head(out10)

out10[,Date:=as.Date(Date,"%d.%m.%Y")]

out10[Region=="All" & !(Age %in% c("TOT","UNK","Unknown")),table(Country,Sex)]
num_obs = out10[Region=="All" & !(Age %in% c("TOT","UNK","Unknown")) & !is.na(Deaths),.(n=length(unique(Date))),by=.(Country)]
only1obs = num_obs[n==1,Country]
lessthan100obs = num_obs[n<100,Country]
lessthan300obs = num_obs[n<300,Country]

ggplot(out10[!is.na(Deaths) & !(Country %in% lessthan100obs) & Region=="All" & !(Age %in% c("TOT","UNK","Unknown")) & Sex=="b"],aes(x=Date,y=Deaths,group=Age,color=Age)) +
    geom_line() + 
    theme(legend.position = "none") +
    facet_wrap(~Country)

names(out10) = tolower(names(out10))
overlapping_countries = intersect(out10[,unique(country)],vax[,unique(country)])
ggplot() +
    geom_line(aes(x=date,y=deaths,group=factor(age),color=factor(age)),out10[!is.na(deaths) & (country %in% overlapping_countries) & region=="All" & !(age %in% c("TOT","UNK","Unknown")) & sex=="b"]) + 
    geom_point(aes(x=date,y=cum_death_i_both,group=age_group,color=age_group),dt[country %in% overlapping_countries],shape=3,size=0.2) +
    # theme(legend.position = "none") +
    facet_wrap(~country)

