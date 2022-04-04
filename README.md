# Mapping social distancing measures to the reproduction number for COVID-19
Ellen Brooks-Pollock, Jonathan M. Read, Angela R. McLean, Matt J. Keeling and Leon Danon (2021) Philosophical Transactions of the Royal Society B, https://doi.org/10.1098/rstb.2020.0276

This is the baseline code for recreating ready reckoners. 

## What are ready reckoners?
Ready reckoners is a nickname given to a method used for relating interventions aimed at limiting infectious disease transmission to the effective reproduction number, R. The ready reckoners were developed during the COVID-19 pandemic for understanding the relationship between interventions and SARS-CoV-2 transmission. 

The basic principle behind the ready reckoners is to use social contact data to calculate individual reproduction numbers, which are combined to calculate a population-level reproduction number. 

The social contact data we used here is from the Social Contact Survey described by Danon et al. in Danon et al. 2012 doi:10.1098/rsif.2012.0357
 and Danon et al. 2013  doi:10.1098/rspb.2013.1037. 
 
## What does the code do?

1. Load in the social contact data
2. Set up the functions for calculating an individual's reproduction number
3. Calibrate to a desired baseline basic reproduction number
4. Re-calculate the reproduction numbers for a range of scenarios (different levels of social distancing, COVID security, contact tracing etc)
5. Plot the results on a single figure 

## How to run the code
All the code is contained in the file "limiting contacts_3by3_comments.R". You should be able to source the file and it will produce the 3-by-3 figure in the output/ directory. 
