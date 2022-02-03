# HNSCC-Cetuximab-Resistance
Data and code associated with "Using Mathematical Modeling to Distinguish Intrinsic and Acquired Targeted Therapeutic Resistance in Head and Neck Cancer"

- Raw experimental data can be found in Raw_Data folder, with control mice found in Control_data.xlsx and cetuximab-treated mice found in CTX_data.xlsx
- Censored experimental data, used for fitting, can be found in the Censored_Data folder
- Control_Fitting folder contains all code to individually fit the 25 control mice: 1) Fit_Exponential_Analytic.m fits to an exponential curve, 2) Fit_Logistic.m fits to a logistic curve, and 3) Fit_Allee fits to an allee curve
- Treatment_Fitting folder contains all code to individually fit the 29 treatment mice (assuming exponential growth): 1) model1_1.m fits to drug model with no resistance, 2) model1_2.m fits to drug model with only pre-existing resistance, 3) model2_1.m fits to drug model with only randomly acquired resistance, 4) model2_2. fits to drug model with both pre-existing and randomly acquired resistance, 5) model3_1.m fits to drug model wtih only drug-induced resistance, and 6) model3_2.m fits to drug model with both pre-existing and drug-induced resistance
