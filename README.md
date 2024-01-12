# New_Haven_Sewage
Additional analyses of coronavirus in Sewage in New Haven for Peccia et al, published in Nature Biotechnology: https://www.nature.com/articles/s41587-020-0684-z

These are additional analyses of data examining increases in SARS-CoV-2 RNA concentrations in primary sewage sludge in 
relation to COVID-19 hospitalizations and cases. The original pre-print ( 
https://www.medrxiv.org/content/10.1101/2020.05.19.20105999v1)  included a correlation of smoothed versions of the viral RNA data 
and the epidemiological indicators. In the analyses presented here, we use an error-in-variables 
time series model to get an estimate of the underlying viral dynamics in sewage and the relationship of this with hospitalizations. 
This analysis is carried out in the Bayesian framework, allowing us to correctly quantify uncertainty in the estimated associations. 
These additional analyses were performed by Dan Weinberger (Epidemiology of Microbial Diseases, Yale School of Public Health),
with input from Josh Warren (Biostatistics, Yale School of Public Health) and the rest of the original study team.

The results of the analyses can be found here: https://weinbergerlab.github.io/New_Haven_Sewage/

On the data file:
'new.cases' is based on positive tests as reported on the department of public health website based on date of report on the website

'hospt.admit' is number of admissions to the hospital in the catchment area

N_positive_reported and N_test_reported is based on positive tests and total tests performed based on the data when the lab reported to DPH

N_positive_received and N_test_received is based on positive samples and tests performed based on date of sample collection.
