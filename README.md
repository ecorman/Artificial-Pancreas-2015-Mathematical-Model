# Artificial-Pancreas-2015-Mathematical-Model

Description of how the model is solved for: https://www.dropbox.com/home/2015_Artificial_Pancreas_Project/01_Mathematical_Modeling/Documentation?preview=Hovorka+Mathematical+Modelling.pdf

Mathematical model for an artificial pancreas system for development of a "virtual copy" of a type 1 diabetes patient for decision-making, or a more complex full-scale #OpenAPS which requires interfacing with the MATLAB Database Toolbox plus NightScout, used for in-silico testing in combination with a nonlinear model predictive (glucose) controller (to be created in the future), which builds upon the mathematical models, in MATLAB with the MATLAB Statistics Toolbox and the MATLAB Simulink extension. Based on Stochastic Virtual Population of Subjects With Type 1 Diabetes for the Assessment of Closed-Loop Glucose Controllers. Requires MATLAB, MATLAB Statistics Toolbox, and MATLAB Simulink extension to run this repository. This is work in progress. 

# Main.m File
This is where the program executes from. The ultimate desired result with this particular file is to generate a joint posterior distribution of 22 parameters unique to a particular type 1 diabetes patient, so that the patient can be represented virtually for in-silico simulation. Future developments include creation of a virtual library of type 1 diabetes patients, for evaluation of differences with different subsets within the population and further analysis for estabilishing safety and efficacy of a nonlinear model predictive (glucose) controller.


