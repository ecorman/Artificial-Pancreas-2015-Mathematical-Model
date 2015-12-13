# Importing relevant packages
import pymc
import numpy as pymc

# c_sub_i, [U/min], background insulin appearance
# Metropolis-Hastings sampling as a Markov chain Monte Carlo
# From: Stochastic Virtual Population of Subjects with Type 
# 1 Diabetes for the Assessment of Closed-Loop Glucose Controllers

c_sub_i = pymc.distributions.truncated_normal_like('c_sub_i',mu=1,tau=1.491824697641270,a=0,b=float("inf"))

# log_k_sub_e, [1/min], fractional clearance rate 
# represented as a normal distribution for Metropolis-Hastings
# sampling as a Markov chain Monte Carlo method
# From: Stochastic Virtual Population of Subjects with Type
# 1 Diabetes for the Assessment of Closed-Loop Glucose Controllers

log_k_sub_e = pymc.Normal('log_k_sub_e',mu=-2.63,tau=4.49)

k_sub_e = numpy.exp(log_k_sub_e)

# log_Q_sub_b, [U], insulin-on-board due to a preceding insulin 
# delivery represented as a normal distribution for 
# Metropolis-Hastings sampling as a Markov chain Monte Carlo   
# method 
# From: Stochastic Virtual Population of Subjects with Type 
# 1 Diabetes for the Assessment of Closed-Loop Glucose Controllers

log_Q_sub_b = numpy.Normal('log_Q_sub_b',mu=-0.7,tau=1.5)

Q_sub_b = numpy.exp(log_Q_sub_b)

# p_sub_i, [unitless], portion of insulin absorbed in the slow
# channel, in the subcutaneous insulin absorption subchannel,
# represented as a uniform prior distribution for 
# Metropolis-Hastings sampling as a Markov chain Monte Carlo 
# method
# From: Stochastic Virtual Population of Subjects with Type
# 1 Diabetes for the Assessment of Closed-Loop Glucose Controllers

log_p_sub_i =pymc.distributions.uniform_like('log_p_sub_i',1,2.718281828459046)

p_sub_i = numpy.exp(log_p_sub_i)

# p_sub_m, [unitless], uniform prior distribution

log_p_sub_m = pymc.distributions.uniform_like('log_p_sub_m',1,2.718281828459046)

p_sub_m = numpy.exp(log_p_sub_m)

# log_k_sub_i_mean, [log(1/min)] [log(1/min)] [log(1/min)],      
# contains the population mean values for log_k_sub_is1 
# (fractional transfer rate), log_k_sub_is2 (fractional 
# transfer rate), and log_k_sub_if 
# (shared fractional transfer rate)

log_k_sub_is1_mean = -3.912

log_k_sub_is2_mean = -3.912

log_k_sub_if_mean = -2.708

log_k_sub_is1_variance = 1/6.25

log_k_sub_is2_variance = 1/6.25

log_k_sub_if_variance = 1/6.24

log_k_sub_is1 = pymc.Normal('log_k_sub_is1', mu=log_k_sub_is1_mean, tau=log_k_sub_is1_variance)

k_sub_is1 = numpy.exp(log_k_sub_is1)

log_k_sub_is2 = pymc.Normal('log_k_sub_is2', mu=log_k_sub_is2_mean, tau=log_k_sub_is2_variance)

k_sub_is2 = numpy.exp(log_k_sub_is2)

log_k_sub_if = pymc.Normal('log_k_sub_if', mu=log_k_sub_if_mean, tau=log_k_sub_if_variance)

k_sub_if = numpy.exp(log_k_sub_if)

# Dealing with dinner distributions

log_k_sub_m_mean_dinner = -3.9

log_d_mean_dinner = 2.3

log_k_sub_m_variance_dinner = 1/12.84

log_d_variance_dinner = (CHO_dinner/100)^2

log_k_sub_m_dinner = pymc.Normal('log_k_sub_m_dinner', mu=log_k_sub_m_mean_dinner, tau=log_k_sub_m_variance_dinner)

k_sub_m_dinner = numpy.exp(log_k_sub_m_dinner)

log_d_dinner = pymc.Normal('log_d_dinner', mu=log_d_mean_dinner, tau=log_d_variance_dinner)

d_dinner = numpy.exp(log_d_dinner)

# dealing with breakfast distributions

log_k_sub_m_mean_breakfast = -3.9

log_d_mean_breakfast = 2.3

log_k_sub_m_variance_breakfast = 1/12.84

log_d_variance_breakfast = (CHO_breakfast/100)^2

log_k_sub_m_breakfast = pymc.Normal('log_k_sub_m_breakfast', mu=log_k_sub_m_mean_breakfast, tau=log_k_sub_m_variance_breakfast)

k_sub_m_breakfast = numpy.exp(log_k_sub_m_breakfast)

log_d_breakfast = pymc.Normal('log_d_breakfast', mu=log_d_mean_breakfast, tau=log_d_variance_breakfast)

d_breakfast = numpy.exp(log_d_breakfast)

# k_sub_12

log_k_sub_12_mean = -2.8

k_sub_12_mean = numpy.exp(log_k_sub_12_mean)

k_sub_12_variance = 1.564600225344889

k_sub_12 = pymc.Normal('k_sub_12', mu=k_sub_12_mean, tau=k_sub_12_variance)

# k_sub_a1

log_k_sub_a1_mean = -5.7

k_sub_a1_mean = numpy.exp(log_k_sub_a1_mean)

k_sub_a1_variance = 1.564082426075518

k_sub_a1 = pymc.Normal('k_sub_a1', mu=k_sub_a1_mean, tau=k_sub_a1_variance)

# k_sub_a2

log_k_sub_a2_mean = -2.9

k_sub_a2_mean = numpy.exp(log_k_sub_a2_mean)

k_sub_a2_variance = 1.564733448492373

k_sub_a2 = pymc.Normal('k_sub_a2', mu=k_sub_a2_mean, tau=k_sub_a2_variance)


# k_sub_a3

log_k_sub_a3_mean = -3.7

k_sub_a3_mean = numpy.exp(log_k_sub_a3_mean)

k_sub_a3_variance = 1.564082370774399

k_sub_a3 = pymc.Normal('k_sub_a3', mu=k_sub_a3_mean, tau=k_sub_a3_variance)

# S_sub_t

log_S_sub_t_mean = 3.7

S_sub_t = numpy.exp(log_S_sub_t_mean)

S_sub_t_variance = 58.238249021516120

S_sub_t = pymc.Normal('S_sub_t', mu=S_sub_t_mean, tau=S_sub_t_variance)

# S_sub_d

log_S_sub_d_mean = 1.6

S_sub_d_mean = numpy.exp(log_S_sub_d_mean)

S_sub_d_variance = 62.391934897146620

S_sub_d = pymc.Normal('S_sub_d', mu=S_sub_d_mean, tau=S_sub_d_variance)

# S_sub_e 

log_S_sub_e_mean = 6

S_sub_e_mean = numpy.exp(log_S_sub_e_mean)

S_sub_e_variance = 30.166534590087520

S_sub_e = pymc.Normal('S_sub_e', mu=S_sub_e_mean, tau=S_sub_e_variance)

# EGP_0

EGP_0_mean = 9.7

EGP_0_variance = 4.561666666666667

EGP_0 = pymc.Normal('EGP_0', mu=EGP_0_mean, tau=EGP_0_variance)

# F_01

F_01_mean = 16.1

F_01_variance = 11.918835953735060

F_01 = pymc.Normal('F_01', mu=F_01_mean, tau=F_01_variance)



