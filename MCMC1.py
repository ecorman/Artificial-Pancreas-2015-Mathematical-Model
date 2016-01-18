from __future__ import division
import pymc
import PlasmaInsulinMCMC as model

mc = pymc.MCMC(model)

mc.use_step_method(pymc.Metropolis, model.log_MCR_sub_I, model.log_ins_sub_c)

mc.sample(iter=100, burn=20, thin=1, verbose=1)
