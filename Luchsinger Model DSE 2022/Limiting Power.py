#import modules and input data
import numpy as np
from Input import *
from Optimal_gamma import *
import matplotlib.pyplot as plt

#assume Power and Tether Limit are reached at the same time
#get gamma data from Optimal gamme for v_w=<v_n
##test

gamma_out_n,gamma_in_n = calculate_opt_gamma_nominal() #nominal reel out from optimal_gamma.py


v_out_n = gamma_out_n*v_w_n
T_out_n = 0.5 * rho * v_w_n**2*A_proj*(1-gamma_out_n)**2*F_out
P_out_n = T_out_n*v_out_n

#when v_w>v_n

v_w = np.linspace(v_w_n,2.5*v_w_n,100)
gamma_in = np.linspace(gamma_in_n, 0.25, 100)
gamma_out = np.linspace(gamma_out_n, 0.05, 100)
gamma_in_max_f_c = np.zeros(100)
gamma_out_max_f_c = np.zeros(100)
f_c_mu = np.zeros(100)

ci = 0
cj = 0

for i in v_w:
    mu = i/v_w_n
    for j in gamma_in:
        f_c_mu[cj] = ((1/(mu**2))*(1-gamma_out_n)**2-(F_in/F_out)*(1+j)**2)*((gamma_out_n*j)/(gamma_out_n+mu*j))
        cj +=1
    max_f_c = np.amax(f_c_mu)
    a = np.where(f_c_mu == max_f_c)
    gamma_in_max_f_c[ci] = gamma_in[a]
    gamma_out_max_f_c[ci] = gamma_out[ci]/mu
    ci +=1
    cj = 0
print (gamma_in_max_f_c,gamma_out_max_f_c,v_w)
    
plt.plot(v_w, gamma_in_max_f_c, label = 'Reel-in speed')
plt.plot(v_w,gamma_out_max_f_c,'--', label = 'Reel-out speed')
plt.xlabel('Wind speed',fontsize = 16)
plt.ylabel('Reel speeds',fontsize = 16)
plt.legend(fontsize = 16)
plt.grid()
plt.show()


