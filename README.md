# Potential COVID-19 burden in European countries

This repository contains code and data for the analyses in 'Unexposed populations and potential COVID-19 burden in European countries as of 21st November 2021' [1]. The aim of this work was to estimate the potential remaining burden of COVID-19 hospitalisations and deaths in European countries in a consistent way, given limited available data in some countries. We do this by:

* deconvolving age-stratified time series of deaths from the COVerAGE [2,3] and INED [4] databases, currently using the `incidental` R package [5,6], to obtain the underlying infection time series
* dividing the deconvolved infection time series by an estimate of the time-varying infection fatality rate (IFR) that accounts for the decrease in the IFR as a result of vaccination (and variation in the average vaccine efficacies against different outcomes over time with vaccine type and changes in circulating variants)
* solving a set of difference equations for the movement of individuals between susceptible, vaccinated and infected states over time to calculate the current numbers of susceptible, vaccinated & uninfected, and previously infected individuals in each age group in each country (see Supplement of [1] for details)
* using these estimates and the age-dependent infection hospitalisation rate [7] and infection fatality rate [8] to calculate the number of hospitalisations and deaths that would occur in each age group in each country if the entire population were to be exposed right now.

The repository contains code for downloading the data required for the analysis, cleaning and processing this data, backcalculating the infection time series, and calculating the maximum remaining burden of hospitalisations and deaths by country and age. Currently 19 countries are included in the analysis as they have the required data (age-stratified death and vaccination data and data on variant proportions).

## Installation

Clone/download this project onto your machine.

The following R packages are required to run the code: 

```
data.table, qs, readxl, ggplot2, ggrepel, osfr, covidregionaldata, ISOweek, lubridate, surveillance, incidental, remotes, EpiNow2, mgcv, nnet, splines, effects, doParallel, cowplot, patchwork, stringr, grDevices
```

and can be installed in R by running:

```R
install.packages(c("data.table","qs","readxl","ggplot2","ggrepel","osfr","covidregionaldata","ISOweek","lubridate","surveillance","incidental","remotes","EpiNow2","mgcv","nnet","splines","effects","doParallel","cowplot","patchwork","stringr","grDevices"))
```

## Data

The analysis uses COVID-19 data from various sources:

* Age-stratified death data:
  * [COVerAGE database](https://osf.io/mpwjq/)
  * [INED database](https://dc-covid.site.ined.fr/en/)
* Vaccination data by age and vaccine product:
  * EU countries (excl. Germany): [ECDC data](https://www.ecdc.europa.eu/en/publications-data/data-covid-19-vaccination-eu-eea)
  * England: [UK government COVID-19 dashboard](https://coronavirus.data.gov.uk/details/download)
  * Germany: [RKI data](https://github.com/robert-koch-institut/COVID-19-Impfungen_in_Deutschland)
* Variant data:
  * EU countries: [ECDC data](https://www.ecdc.europa.eu/en/publications-data/data-virus-variants-covid-19-eueea) (compiled from [GISAID](https://www.gisaid.org/) and [TESSY)](https://www.ecdc.europa.eu/en/publications-data/european-surveillance-system-tessy)
  * England: [COG UK Consortium data for England](https://covid19.sanger.ac.uk/downloads)
* Seroprevalence data:
  * EU countries: [SeroTracker database](https://serotracker.com/)
  * England: [ONS Coronavirus Infection Survey](https://www.ons.gov.uk/peoplepopulationandcommunity/healthandsocialcare/conditionsanddiseases/datasets/coronaviruscovid19antibodydatafortheuk/2021)

The code includes functions to download and save this data locally, but the input data for the analysis (available under [data](data)) has currently been fixed at 21st November 2021 due to the emergence of Omicron, which is not yet included in our analysis.

## Running the code

The full analysis to calculate the remaining burden estimates can then be run in R by changing directory to where the code was downloaded and entering

```R
source("run_remaining_burden_analysis.R")
```

or from the command line on a Mac/Linux machine with:

```
Rscript run_remaining_burden_analysis.R
```

Note that the code takes a long time to run (~6hrs running in parallel in R4.1.0 on an 8-core MacBook Pro with an M1 chip and 16GB RAM) and requires a large amount of RAM (at least 16GB), so you may wish to only run parts of the analysis at a time or need to obtain access to high-performance computing resources.

The estimates of the potential remaining burden by country and age can be plotted with:

```
source("plot_remaining_burden_estimates.R")
```

## Output

The full infection backcalculation output is available on Zenodo ([https://doi.org/10.5281/zenodo.5772163](https://doi.org/10.5281/zenodo.5772163)), and the [output](output) folder contains the potential remaining burden estimates at [country level](output/2021-11-30/ovrl_rem_burden_output.csv) and by [country and age group](output/2021-11-30/rem_burden_output.csv).

## Built With

* [R version 4.1.0 (2021-05-18)](https://www.r-project.org/)

## Authors

* Lloyd Chapman: <lloyd.chapman1@lshtm.ac.uk> 

## License

This project is licensed under the GNU General Public License v3.0 - see the [LICENSE.txt](LICENSE.txt) file for details

## References
1. Chapman LAC, Barnard RC, Russell TW, Abbott S, Van Zandvoort K, Davies NG, Kucharski AJ. Unexposed populations and potential COVID-19 burden in European countries as of 21st November 2021 [https://doi.org/10.1101/2021.11.10.21266166](https://doi.org/10.1101/2021.11.10.21266166).

2. Riffe T, Acosta E, the COVerAGE-DB team, Acosta EJ, Manuel Aburto D, Alburez-Gutierrez A, et al. Data Resource Profile: COVerAGE-DB: a global demographic database of COVID-19 cases and deaths. Int J Epidemiol. 2021 May 15;50(2):390–390f.

3. Riffe T, Acosta E, Schöley J, Donzowa J, Kniffka MS. COVerAGE-DB: A database of COVID-19 cases and deaths by age. 2020 Apr 9 [cited 2021 Oct 21]; Available from: [https://osf.io/mpwjq/](https://osf.io/mpwjq/)

4. Demographics of COVID-19 Deaths [Internet]. [cited 2021 Oct 21]. Available from: [https://dc-covid.site.ined.fr/en/](https://dc-covid.site.ined.fr/en/)

5. Miller AC, Hannah L, Futoma J, Foti NJ, Fox EB, D’Amour A, et al. Statistical deconvolution for inference of infection time series. medRxiv. 2020 Oct 20;2020.10.16.20212753.

6. Miller A, Hannah L, Foti N, Futoma J, Apple Inc. incidental R package version 0.1. 2020 Sep 16 [cited 2021 Oct 20]; Available from: [https://CRAN.R-project.org/package=incidental](https://CRAN.R-project.org/package=incidental)

7. Salje H, Tran Kiem C, Lefrancq N, Courtejoie N, Bosetti P, Paireau J, et al. Estimating the burden of SARS-CoV-2 in France. Science. 2020 Jul 10;369(6500):208–11.

8. O’Driscoll M, Ribeiro Dos Santos G, Wang L, Cummings DAT, Azman AS, Paireau J, et al. Age-specific mortality and immunity patterns of SARS-CoV-2. Nature. 2021 Feb;590(7844):140–5.
