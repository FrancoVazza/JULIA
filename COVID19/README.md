

Simplistic monte-carlo model for the spread of the COVID19 virus in Italy (starting from Feb.21st).

The Julia code requires Julia v.1.0+ and a bunch of packages (see beginning of the macro). 

The model "predicts" the evolution of the number of infected, diagnosed and casualities based on the interplay of 3 parameters:

- **pinf0**= probability that an infected (but undiagnosed) person can infect someone else on each day;
- **pdead**= probability that an infected person dies;
- **t_incub** = mean incubiation period (after which one gets **always** diagnosed); assuming a gaussian distribution with FWHM t_incub/2. 

Ref. numbers that seems to work against real data are pinf0=0.4, pdead=0.01 and t_incub=5.5 days 

Main assumptions:

- we assume a single large population (i.e. no clusterting in cities etc) with a single age. 
- an initial population of **n_inf0** infected (larger than the initial number of reported cases - but to be calibrated against them). n_inf0=400 works well in matching the initial normalisation of italian data. 
- a person is only infective while its personal counter has a time < its personal incubation time (randomly drawn); after this the person is diagnosed/isolated in 100% of cases (but can still die). 
- all infected people can infect others and will eventually develop a desease (questionable). 

The model can also account for the long-term effect of an effective reduced contagion rate, by reducing pinf0 at a given day; various histories can be compared. 

The model runs in a few tens of seconds to generate ~10 realization of the evolution over ~20 days, processing ~1e5 infected. 

FV, 6 March 2020 
