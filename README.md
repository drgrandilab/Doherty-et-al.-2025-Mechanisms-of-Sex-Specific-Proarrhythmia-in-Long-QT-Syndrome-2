Matlab code recreating the figures from Doherty et al. (Am J Physiol Heart Circle Physical. 2025).

Updated model versions including male and female parameterized versions of the Tomek-Rodriguez (ToR-ORd), O’Hara-Rudy (ORd), and Grandi-Bers (GB). Additionally, INaL current formulation from ORd model was integrated into GB model.

_____________________________________________________________________________________________________________
Contents:

Folder “ToR-ORd” : 

Torord_singlerunner	        loads initial conditions and runs baseline simulations
Torord_script_ICaL	        loads initial conditions and runs ICaL Gain-of-Function simulations
Torord_script_IKr	        loads initial conditions and runs IKr Loss-of-Function simulations
model_Torord	                excitation-contraction coupling model
modelRunner	                computes currents 
getCurrentStructure	        assembles the current structure of ToR-ORd model
getStartingState	        initializes initial conditions
ratescript_torord	        plots rate dependent biomarkers as shown in Fig 1B, Fig 2B, and Fig 4B

Folder “ORd” : 

main_ORd	                loads initial conditions and runs baseline simulations
model_ORd	                excitation-contraction coupling model
ratescript_ORd	                plots rate dependent figures of biomarkers

Folder “Grandi” : 

Grandi_main_LQT	                loads initial conditions and runs baseline simulations
Grandi_model_LQT	        excitation-contraction coupling model
yfin_endo_0p5Hz	                initial condition male (2000 BCL where LQT simulation is ran)
yfin_endo_female_0p5Hz	        initial condition female (2000 BCL where LQT simulation is ran)
ratescript_grandi	        plots rate dependent figures of biomarkers

Folder “Population-level Analysis”:

logistic_regression_analysis_LQT	performs logistic regression analysis
populationplots	                    	plots Fig 3A,B and Fig 5A,B
regressionplots	                   	plots results of regression coefficients as shown in Fig 3C and 5C

.mat files in sub-folder “Population-level analysis/Conditions” contains data extracted from populations for each model
________________________________________________________________________________________________________________________________

Reference: 

I. Doherty, R. Shetty, H. Ni, S. Morotti, E. Grandi. Exploring the mechanisms of sex-specific proarrhythmia in long QT syndromethrough computational modeling. Am J Physiol Heart Circle Physiol. 2025.

