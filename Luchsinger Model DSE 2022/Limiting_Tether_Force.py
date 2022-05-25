import numpy as np
from Input import *
from Optimal_gamma import *
import matplotlib.pyplot as plt
from Limiting_Power import *



def calculate_tether_limit_reel_speeds(step):
    
    gamma_out_n,gamma_in_n = calculate_opt_gamma_nominal() #nominal reel out from optimal_gamma.py

    F_in = 0.07   
    F_out = 5.4

    v_w = np.linspace(v_w_n,2.5*v_w_n,step)
    gamma_in = np.linspace(gamma_in_n, 0.25, step)
    gamma_out = np.linspace(gamma_out_n, 0.05, step)
    gamma_in_max_f_c = np.zeros(step)
    gamma_out_max_f_c = np.zeros(step)
    f_c_mu = np.zeros(step)
    P_c = np.zeros(step)

    ci = 0
    cj = 0

    for i in v_w:
        mu = i/v_w_n
        P_w = 0.5*rho*i**3
        for j in gamma_in:
            f_c_mu[cj] = ((1/(mu**2))*(1-gamma_out_n)**2-(F_in/F_out)*(1+j)**2)*((j*(mu-1+gamma_out_n))/(mu*j+mu-1+gamma_out_n))
            cj +=1
        max_f_c = np.amax(f_c_mu)
        a = np.where(f_c_mu == max_f_c)
        gamma_in_max_f_c[ci] = gamma_in[a]
        gamma_out_max_f_c[ci] = 1-((1-gamma_out_n)/mu)
        P_c [ci] = P_w*max_f_c
        ci +=1
        cj = 0
    #print (gamma_in_max_f_c,gamma_out_max_f_c,v_w)
    return gamma_in_max_f_c,gamma_out_max_f_c,v_w,P_c

def plot_tether_limit_reel_speeds(step):
    gamma_in_max_f_c,gamma_out_max_f_c,v_w,P_c = calculate_tether_limit_reel_speeds(step)    
    plt.plot(v_w, gamma_in_max_f_c, label = 'Reel-in speed')
    plt.plot(v_w,gamma_out_max_f_c,'--', label = 'Reel-out speed')
    plt.xlabel('Wind speed',fontsize = 16)
    plt.ylabel('Reel speeds',fontsize = 16)
    plt.legend(fontsize = 16)
    plt.grid()
    plt.show()

def plot_tether_limit_cycle_power(step):
    gamma_in_max_f_c,gamma_out_max_f_c,v_w,P_c = calculate_tether_limit_reel_speeds(step)
    power_before_v_w_n,v_w_before_v_n = calculate_power_before_v_w_n(step)    
    plt.plot(v_w, P_c, label = 'Cycle Power')
    plt.plot(v_w_before_v_n,power_before_v_w_n)
    plt.xlabel('Wind speed',fontsize = 16)
    plt.ylabel('P_c',fontsize = 16)
    plt.legend(fontsize = 16)
    plt.grid()
    plt.show()

#plot_tether_limit_reel_speeds(100)
plot_tether_limit_cycle_power(100)