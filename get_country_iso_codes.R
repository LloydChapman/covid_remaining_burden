library(data.table)
library(covidregionaldata)

death_data = get_national_data(source="who")
setDT(death_data)
country_iso_codes = unique(death_data[!is.na(country),.(country,iso_code)])
saveRDS(country_iso_codes,"./fitting_data/country_iso_codes.rds")
