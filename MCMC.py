# Importing relevant packages
import pymc3 as pymc
import numpy
import scipy

# c_sub_i, [U/min], background insulin appearance
# Metropolis-Hastings sampling as a Markov chain Monte Carlo
# From: Stochastic Virtual Population of Subjects with Type 
# 1 Diabetes for the Assessment of Closed-Loop Glucose Controllers

# c_sub_i = pymc.distributions.truncated_normal_like('c_sub_i',mu=1,tau=1.491824697641270,a=0,b=0)


BoundedNormal = pymc.Bound(pymc.Normal, lower=0, upper=numpy.inf)
c_sub_i = BoundedNormal('c_sub_i', mu=0, tau=0.4)
  
	
# log_k_sub_e, [1/min], fractional clearance rate 
# represented as a normal distribution for Metropolis-Hastings
# sampling as a Markov chain Monte Carlo method
# From: Stochastic Virtual Population of Subjects with Type
# 1 Diabetes for the Assessment of Closed-Loop Glucose Controllers

k_sub_e = pymc.Normal('k_sub_e',mu=0.072078,tau=89.12144)

# log_Q_sub_b, [U], insulin-on-board due to a preceding insulin 
# delivery represented as a normal distribution for 
# Metropolis-Hastings sampling as a Markov chain Monte Carlo   
# method 
# From: Stochastic Virtual Population of Subjects with Type 
# 1 Diabetes for the Assessment of Closed-Loop Glucose Controllers

Q_sub_b = pymc.Normal('Q_sub_b',mu=0.4965853,tau=4.481689)

# p_sub_i, [unitless], portion of insulin absorbed in the slow
# channel, in the subcutaneous insulin absorption subchannel,
# represented as a uniform prior distribution for 
# Metropolis-Hastings sampling as a Markov chain Monte Carlo 
# method
# From: Stochastic Virtual Population of Subjects with Type
# 1 Diabetes for the Assessment of Closed-Loop Glucose Controllers

p_sub_i = pymc.Uniform('p_sub_i',0,1)

# p_sub_m, [unitless], uniform prior distribution

p_sub_m = pymc.Uniform('p_sub_m',0,1)

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

k_sub_is1 = pymc.Normal('k_sub_is1', mu=numpy.exp(log_k_sub_is1_mean), tau=numpy.exp(log_k_sub_is1_variance))

k_sub_is2 = pymc.Normal('k_sub_is2', mu=numpy.exp(log_k_sub_is2_mean), tau=numpy.exp(log_k_sub_is2_variance))

log_k_sub_if = pymc.Normal('k_sub_if', mu=numpy.exp(log_k_sub_if_mean), tau=numpy.exp(log_k_sub_if_variance))

# Dealing with dinner distributions

log_k_sub_m_mean_dinner = 0.020241911445804

log_d_mean_dinner = 9.974182454814718

log_k_sub_m_variance_dinner = 1.080994682877431

log_d_variance_dinner = (100/CHO_dinner)^(1/2)

k_sub_m_dinner = pymc.Normal('k_sub_m_dinner', mu=numpy.exp(log_k_sub_m_mean_dinner), tau=numpy.exp(log_k_sub_m_variance_dinner))

d_dinner = pymc.Normal('d_dinner', mu=numpy.exp(log_d_mean_dinner), tau=numpy.exp(log_d_variance_dinner))


# dealing with breakfast distributions

log_k_sub_m_mean_breakfast = -3.9

log_d_mean_breakfast = 2.3

log_k_sub_m_variance_breakfast = 1/12.84

log_d_variance_breakfast = (CHO_breakfast/100)^2

k_sub_m_breakfast = pymc.Normal('k_sub_m_breakfast', mu=numpy.exp(log_k_sub_m_mean_breakfast), tau=numpy.exp(log_k_sub_m_variance_breakfast))

d_breakfast = pymc.Normal('d_breakfast', mu=numpy.exp(log_d_mean_breakfast), tau=numpy.exp(log_d_variance_breakfast))



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

# insulin infusion rate--hourly, iir_h

iir_h = 1.4

# insulin infusion rate--minutes, u_sub_i

u_sub_i = iir_h/60

################################################################
##### Start of CAR Stuff #######################################
################################################################

adj = numpy.array([[2,4,5],
	[1,4,5,6,3],
	[2,5,6],
	[1,2,5,8,7],
	[1,2,3,4,6,7,8,9],
	[3,2,5,8,9],
	[4,5,8,11,10],
	[4,5,6,7,9,10,11,12],
	[6,5,8,11,12],
	[7,8,11,14,13],
	[7,8,9,10,12,13,14,15],
	[9,8,11,14,15],
	[10,11,14],
	[13,10,11,12,15],
	[12,11,14]])


# y0 

Q_sub_is1_0 = (u_sub_i*p_sub_i.item())/(k_sub_is1.item())

Q_sub_is2_0 = (Q_sub_is1_0*k_sub_is1.item())/(k_sub_is2.item())

Q_sub_if1_0 = 

# Gut Glucose Absorption Subsystem


# deterministic compartmental model
# @pymc.deterministic
# def Sto(c_sub_i_final = c_sub_i, k_sub_e_final = k_sub_e, Q_sub_b_final = Q_sub_b, p_sub_i_final = p_sub_i, p_sub_m_final = p_sub_m, k_sub_is1_final = k_sub_is1, k_sub_is2_final = k_sub_is2, k_sub_if_final = k_sub_if, k_sub_m_dinner_final = k_sub_m_dinner, d_dinner_final = d_dinner, k_sub_m_breakfast_final = k_sub_m_breakfast, d_breakfast_final = d_breakfast, k_sub_12_final = k_sub_12, k_sub_a1_final = k_sub_a1, k_sub_a2_final = k_sub_a2, k_sub_a3_final = k_sub_a3, S_sub_t_final = S_sub_t, S_sub_d_final = S_sub_d, S_sub_e_final = S_sub_e, EGP_0_final = EGP_0, F_01_final = F_01, u_sub_i_final = u_sub_i, I_sub_m_final = I_sub_m, I_sub_p_final = I_sub_p, f_sub_g_final = f_sub_g):
#    def glu_model(t):
#        dQ_sub_is1 = ((u_sub_i_final.item())*p_sub_i_final.item()) - ((Q_sub_i)*k_sub_is1_final.item())
#        dQ_sub_is2 = ((Q_sub_is1)*(k_sub_is1_final.item()))-((Q_sub_is2)*k_sub_is2_final.item())
#        dQ_sub_if1 = ((u_sub_i_final.item())*(1-p_sub_i_final.item())) - ((Q_sub_if1)*k_sub_if_final.item())
#        dQ_sub_if2 = ((Q_sub_if1)*(k_sub_if_final.item())) - ((Q_sub_if2)*(k_sub_if_final.item()))
#        dQ_sub_i = (((Q_sub_is2)*(k_sub_is2_final.item()) + ((Q_sub_if2)*(k_sub_if_final.item())))*I_sub_m_final.item()) - ((Q_sub_i)*(k_sub_e_final.item())) + c_sub_i_final.item()
#        dx_sub_1 = ((-k_sub_a1_final.item())*x_sub_1) + ((k_sub_a1_final.item())*(S_sub_t_final.item())*(I_sub_p_final.item()))
#        dx_sub_2 = ((-k_sub_a2_final.item())*x_sub_2) + ((k_sub_a2_final.item())*(S_sub_d_final.item())*(I_sub_p_final.item()))
#        dx_sub_3 = ((-k_sub_a3_final.item())*x_sub_3) + ((k_sub_a3_final.item())*(S_sub_e_final.item())*(I_sub_p_final.item()))
#        dQ_sub_1 = (-F_01)*((Q_sub_1/V)/(1+(Q_sub_1/V))) - (x_sub_1*Q_sub_1) + ((k_sub_12_final.item())*Q_sub_2) + ((EGP_0_final.item())*(1-x_sub_3)) + f_sub_g_final.item()
#        dydt = [dQ_sub_is1, dQ_sub_is2, dQ_sub_if1, dQ_sub_if2, dQ_sub_i, dx_sub_1, dx_sub_2, dx_sub_3, d_Q_sub_1, ]
#        return dydt
#    soln = odeint(glu_model, y0, tspan)
#    cc, cp = soln[:,0], soln[:,1]
#    return [cc, cp]
# cc = pc.Lambda('cc', lambda PK=PK: PK[0])
# cp = pc.Lambda('cp', lambda PK=PK: PK[1])



