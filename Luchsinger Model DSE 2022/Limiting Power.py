#import modules and input data
import numpy as np
from Input import *
from Optimal_gamma import *

#assume Power and Tether Limit are reached at the same time
#get gamma data from Optimal gamme for v_w=<v_n

gamma_out_n,gamma_in_n = calculate_opt_gamma_nominal() #nominal reel out from optimal_gamma.py


v_out_n = gamma_out_n*v_w_n
T_out_n = 0.5 * rho * v_w_n**2*A_proj*(1-gamma_out_n)**2*F_out
P_out_n = T_out_n*v_out_n

#when v_w>v_n

v_w = np.linspace(v_w_n,2.5*v_w_n,100)
gamma_in = np.linspace(gamma_in_n, 0.25, 100)
gamma_out = np.linspace(gamma_out_n, 0.05, 100)
gamma_out_max_f_c = []
f_c_mu = np.zeros(100)

ci = 0
cj = 0

for i in v_w:
    mu = i/v_w_n
    for j in gamma_in:
        gamma_in = j
        f_c_mu[cj] = ((1/(mu**2))*(1-gamma_out_n)**2-(F_in/F_out)*(1+gamma_in)**2)*((gamma_out_n*gamma_in)/(gamma_out_n+mu*gamma_in))
        cj +=1
    max_f_c = np.amax(f_c_mu)
    a = np.where(f_c_mu == max_f_c)
    gamma_in_max_f_c = np.append(gamma_in[a])
    ci +=1
    cj = 0

