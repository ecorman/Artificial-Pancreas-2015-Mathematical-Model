% The result of running this file will be obtaining (blood) plasma insulin
% concentrations which will be used as inputs for solving for the 22
% parameters in Stochastic Virtual Population of Subjects
% NOTE: You cannot have insulin-on-board (IOB) from residual insulin
% boluses either during the experiment or when solving for a solution

% clear, clc, close
clear, clc, close

% format long
format long

% Initial value for i_sub_1_of_0, [Units]
% NOTE: You cannot have insulin-on-board (IOB) from residual insulin
% boluses either during the experiment or when solving for a solution
i_sub_1_of_t = 0; % [Units]

% Initial value for i_sub_2_of_0, [Units]
% NOTE: You cannot have insulin-on-board (IOB) from residual insulin
% boluses either during the experiment or when solving for a solution
i_sub_2_of_t = 0; % [Units]

% t_sub_max_prior_distribution
t_sub_max_prior_distribution = makedist('Normal','mu',log(60),'sigma',100);

% MCR_prior_distribution
MCR_prior_distribution = makedist('Normal','mu',log(0.01),'sigma',100);

% ins_sub_c_prior_distribution
ins_sub_c_prior_distribution = makedist('Normal','mu',log(36),'sigma',100);

% Initial guesses for t_sub_max
t_sub_max_initial_guess_1 = 60;
t_sub_max_initial_guess_2 = 60;
t_sub_max_initial_guess_3 = 60;

% Initial guesses for MCR
MCR_initial_guess_1 = 0.01;
MCR_initial_guess_2 = 0.01;
MCR_initial_guess_3 = 0.01;

% Initial guesses for ins_sub_c
ins_sub_c_initial_guess_1 = 36;
ins_sub_c_initial_guess_2 = 36;
ins_sub_c_initial_guess_3 = 36;

% Initialize iteration at t = 1
t = 1;
