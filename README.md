# HNSCC-Cetuximab-Resistance
Data and code associated with "Using Mathematical Modeling to Distinguish Intrinsic and Acquired Targeted Therapeutic Resistance in Head and Neck Cancer"
Raw experimental data can be found in Raw_Data folder, with control mice found in Control_data.xlsx and cetuximab-treated mice found in CTX_data.xlsx
Censored experimental data, used for fitting, can be found in the Censored_Data folder. 
Control_Fitting folder contains all code to individually fit the 25 control mice: 1) Fit_Exponential_Analytic.m fits to an exponential curve, 2) Fit_Logistic.m fits to a logistic curve, and 3) Fit_Allee fits to an allee curve
Treatment_Fitting folder contains all code to individually fit the 29 treatment mice:
