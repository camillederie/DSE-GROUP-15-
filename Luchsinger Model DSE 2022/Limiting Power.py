#import modules and input data
import numpy as np
from Input import *
from Optimal_gamma import *
import matplotlib.pyplot as plt

#assume Power and Tether Limit are reached at the same time
    
def calculate_power_limit_reel_speeds():

    #call optimal nominal gamma
    gamma_out_n,gamma_in_n = calculate_opt_gamma_nominal() #nominal reel out from optimal_gamma.py
    
    #create lists to have running variables
    v_w = np.linspace(v_w_n,2.5*v_w_n,100)
    gamma_in = np.linspace(gamma_in_n, 0.25, 100)
    gamma_in_max_f_c = np.zeros(100)
    gamma_out_max_f_c = np.zeros(100)
    f_c_mu = np.zeros(100)

    #create counters and set to zero
    ci = 0
    cj = 0

    #loop to go through all the different gamma_in for all v_w
    #also calculates gamma_out
    for i in v_w:
        mu = i/v_w_n
        for j in gamma_in:
            f_c_mu[cj] = ((1/(mu**2))*(1-gamma_out_n)**2-(F_in/F_out)*(1+j)**2)*((gamma_out_n*j)/(gamma_out_n+mu*j))
            cj +=1
        max_f_c = np.amax(f_c_mu)
        a = np.where(f_c_mu == max_f_c)
        gamma_in_max_f_c[ci] = gamma_in[a]
        gamma_out_max_f_c[ci] = gamma_out_n/mu
        ci +=1
        cj = 0
    print ("gamma_in=",gamma_in_max_f_c,"gamma_out",gamma_out_max_f_c,"wind_speed=",v_w)
    return gamma_in_max_f_c,gamma_out_max_f_c,v_w

def plot_power_limit_reel_speeds():
    gamma_in_max_f_c,gamma_out_max_f_c,v_w = calculate_power_limit_reel_speeds()   
    plt.plot(v_w, gamma_in_max_f_c, label = 'Reel-in speed')
    plt.plot(v_w,gamma_out_max_f_c,'--', label = 'Reel-out speed')
    plt.xlabel('Wind speed',fontsize = 16)
    plt.ylabel('Reel speeds',fontsize = 16)
    plt.legend(fontsize = 16)
    plt.grid()
    plt.show()

plot_power_limit_reel_speeds()

