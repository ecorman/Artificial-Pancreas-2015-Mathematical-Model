from __future__ import division
import pymc
import numpy
from scipy.integrate import odeint


W = 110 # kg

log_t_max = pymc.Normal('log_t_max', mu=numpy.log(60), tau=100)

t_max = pymc.exp(log_t_max)

log_MCR_sub_I = pymc.Normal('log_MCR_sub_I', mu=numpy.log(0.01),tau=100)

MCR_sub_I = pymc.exp(log_MCR_sub_I)

log_ins_sub_c = pymc.Normal('log_ins_sub_c', mu=numpy.log(36),tau=100)

ins_sub_c = pymc.exp(log_ins_sub_c)

u_of_t_min_0_15 = 1.4/60

bolus_delivery_at_t_equals_0 = 0

IOB = 1.3

i_sub_1_of_t_0 = (u_of_t_min_0_15*t_max) + bolus_delivery_at_t_equals_0

i_sub_2_of_t_0 = i_sub_1_of_t_0 + IOB

tspan_0_15 = 0,15

y = i_sub_1_of_t_0, i_sub_2_of_t_0

# deterministic compartmental model
@pymc.deterministic
def Sto_0_15(t_max = t_max, MCR_sub_I = MCR_sub_I, ins_sub_c = ins_sub_c):
	def ins_model_0_15(y,t):
		di_sub_1_of_t_0_15 = u_of_t_min_0_15-(i_sub_1_of_t_0/t_max.item())
		di_sub_2_of_t_0_15 = (-1/t_max.item())*(i_sub_2_of_t_0-i_sub_1_of_t_0)
		dydt_0_15=[di_sub_1_of_t_0_15,di_sub_2_of_t_0_15]
		return dydt_0_15
	soln_0_15 = odeint(ins_model_0_15, y, tspan_0_15)
	i_sub_1_of_t_15, i_sub_2_of_t_15 = soln_0_15[:,0], soln_0_15[:,1]
	I_plasma_of_t = ((1000/(MCR_sub_I.item()*W))*i_sub_2_of_t_0) + ins_sub_c.item()
	return [i_sub_1_of_t_15, i_sub_2_of_t_15]







