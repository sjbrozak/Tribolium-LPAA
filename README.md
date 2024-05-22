# Tribolium-LPAA
Code used in manuscript titled, "Dynamics of an LPAA model: Insights into population chaos".

## Data
The data is located in ``tribolium_data_pilot.xlsx`` and sorted into sheets by experimental group. The spreadsheet contains a README with additional information about how the spreadsheet is organized.

## MATLAB code

- ``onestep_fit_expgroups.m``: fits experimental subgroups by minimizing sum of squares of one-step residuals. Generates plots comparing the LPA and LPAA models, as well as a QQ plot of the LPAA residuals. Figures are not automatically saved due to formatting bugs.
- ``X_bif_lyap.m``: generates bifurcation diagrams and Lyapunov exponent diagrams.

Final parameter estimation results are located in ``one-step_parameters.xlsx``.
