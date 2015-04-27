# Artificial-Pancreas-2015-Mathematical-Model

Description of how the model is solved for: https://www.dropbox.com/home/2015_Artificial_Pancreas_Project/01_Mathematical_Modeling/Documentation?preview=Hovorka+Mathematical+Modelling.pdf

Mathematical Model of Artificial Pancreas for a development of a "virtual copy" of a type 1 diabetes patient for decision-making, used for in-silico testing in combination with a nonlinear model predictive (glucose) controller (to be created in the future) in MATLAB with the MATLAB Statistics Toolbox and the MATLAB Simulink extension. Based on Stochastic Virtual Population of Subjects With Type 1 Diabetes for the Assessment of Closed-Loop Glucose Controllers. Requires MATLAB, MATLAB Statistics Toolbox, and MATLAB Simulink extension to run. This is work in progress. 

# Main.m File
This is where the program executes from, and where the mathematical modeling occurs. The ultimate desired result with this particular file is to generate a joint posterior distribution of 22 paremeters, which represents a single diabetes patient. Future developments include creation of a virtual library of type 1 diabetes patients, for evaluation of differences with different subsets within the population and further analysis.


