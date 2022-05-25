import numpy as np
from Optimal_gamma import *

F_in = 0.07   
F_out = 5.4
v_t_n = 8
v_p_n = 12

def calculate_three_phase(step):
    gamma_out_n,gamma_in_n = calculate_opt_gamma_nominal() #nominal reel out from optimal_gamma.py

    v_w = np.linspace(0,2.5*v_w_n,step)
    gamma_in = np.linspace(gamma_in_n, 0.25, step)
    gamma_in_max_f_c = np.zeros(step)
    gamma_out_max_f_c = np.zeros(step)
    f_c_mu = np.zeros(step)
    P_c = np.zeros(step)

    ci = 0
    cj = 0

    for i in v_w:

        if i <= v_t_n:

            P_w = 0.5*rho*i**3
            P_c [ci] = (1/F_out)*P_w*(F_out*(1-gamma_out_n)**2-(F_in*(1+gamma_in_n)**2))*((gamma_out_n*gamma_in_n)/(gamma_out_n+gamma_in_n))
        elif i <= v_p_n:

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
            cj = 0    
        else:

            #loop to go through all the different gamma_in for all v_w
            #also calculates gamma_out
            mu = i/v_w_n
            P_w = 0.5*rho*i**3
            #F_out_mu = (F_out/(mu**2))*(((1-gamma_out_n)**2)/((1-(gamma_out_n/mu))**2))
            for j in gamma_in:
                f_c_mu[cj] = ((1/(mu**2))*(1-gamma_out_n)**2-(F_in/F_out)*(1+j)**2)*((gamma_out_n*j)/(gamma_out_n+mu*j))
                cj +=1
            max_f_c = np.amax(f_c_mu)
            a = np.where(f_c_mu == max_f_c)
            gamma_in_max_f_c[ci] = gamma_in[a]
            gamma_out_max_f_c[ci] = gamma_out_n/mu
            P_c [ci] = P_w*max_f_c
            cj = 0
        ci +=1
    return gamma_in_max_f_c,gamma_out_max_f_c,v_w,P_c

def plot_three_phase_reel_speeds(step):
    gamma_in_max_f_c,gamma_out_max_f_c,v_w,P_c = calculate_three_phase(step)  
    plt.plot(v_w, gamma_in_max_f_c, label = 'Reel-in speed')
    plt.plot(v_w,gamma_out_max_f_c,'--', label = 'Reel-out speed')
    plt.xlabel('Wind speed',fontsize = 16)
    plt.ylabel('Reel speeds_three_phase',fontsize = 16)
    plt.legend(fontsize = 16)
    plt.grid()
    plt.show()

def plot_three_phase_cycle_power(step):
    gamma_in_max_f_c,gamma_out_max_f_c,v_w,P_c = calculate_three_phase(step)
    plt.plot(v_w, P_c)
    plt.xlabel('Wind speed',fontsize = 16)
    plt.ylabel('P_c_three_phase',fontsize = 16)
    plt.legend(fontsize = 16)
    plt.grid()
    plt.show()

plot_three_phase_reel_speeds(100)
plot_three_phase_cycle_power(100)    
