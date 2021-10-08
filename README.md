# monte_carlo_simulations

This repository includes a monte carlo experiment (an R file with the simulations and its related report). The experiment focuses on the usage of IV-analysis with panel data. In a working paper, I use instrumental variable analysis to re-evaluate the effect of gender inequality indicators on armed conflict on- set. The instruments, geo-climatic suitability for plough-positive cereals or plough-negative cereals are on the country-level. The gender inequality, conflict onset, and control variables are all available on the country-year level (panel data). The data are thus clustered by design. Clustered data can be seen as a special form of omitted variable bias where the cluster(s) account for within-group variation. This can lead to bias and inefficiency. The question posed here is whether one must account for clustered data in a instrumental variable framework.

The monte carlo experiment suggests that, if the instrument is clustered and thus assigns the treatment in a clustered manner, the standard errors need to be adjusted.




