from cmath import cos
#from shlex import join

import numpy as np
import matplotlib.pyplot as plt
from InputV2 import *

data = get_initial_data()
#data = run_nominal_analysis(data)

F_in = 0.103
F_out = 55.72
A_proj = 15.19
#v_t_n = 7.5
P_out_max = 50000
T_out_n = 11100


def calculate_three_phase(step):
    # gamma_out_n,gamma_in_n = data['gamma_out_n'],data['gamma_in_n'] #nominal reel out from optimal_gamma.py

    v_w = np.linspace(0,2.5*data['v_w_n'],step)
    gamma_in = np.linspace(0.1, 4, step)
    gamma_out = np.linspace(0.1,1,step)
    gamma_in_max_f_c = np.zeros(step)
    gamma_out_max_f_c = np.zeros(step)
    f_c = np.zeros((step,step))
    f_c_mu = np.zeros(step)
    P_c = np.zeros(step)
    reel_out_speed = np.zeros(step)
    reel_in_speed = np.zeros(step)
    
    ci = 0
    cj = 0
    ck = 0
    d = 0

    # gamma_out_n = 0.43
    # gamma_in_n = 1.7706060606060607

    for j in gamma_in:
        for k in gamma_out:
            f_c[cj][ck] = (data['eff_out']*((np.cos(data['a_elev_out'])-k)**2)-((F_in/F_out)*(j**2+2*np.cos(data['a_elev_in'])*j+1))/data['eff_in'])*((k*j)/(k+j))
            ck +=1
        ck = 0
        cj +=1
    max_f_c = np.amax(f_c)
    (a,b) = np.where(f_c == max_f_c)
    gamma_in_n = gamma_in[a]
    gamma_out_n = gamma_out[b]
    print(gamma_out_n,gamma_in_n)
    cj = 0
    v_t_n = np.sqrt((T_out_n/(0.5*data['rho']*A_proj*((np.cos(data['a_elev_out'])-gamma_out_n)**2)*F_out)))
    print(v_t_n)
    v_p_n = 10
    for i in v_w:

        if i <= v_t_n:

            P_w = 0.5*data['rho']*i**3
            #P_w = 0.5*rho*i**3
            f_c_mu[ci] = (data['eff_out']*((np.cos(data['a_elev_out'])- gamma_out_n)**2)-((F_in/F_out)*(gamma_in_n**2+2*np.cos(data['a_elev_in'])*gamma_in_n+1))/data['eff_in'])*(( gamma_out_n*gamma_in_n)/( gamma_out_n+gamma_in_n))
            P_c [ci] = P_w*f_c_mu[ci]*F_out*A_proj    
            gamma_in_max_f_c [ci] = gamma_in_n
            gamma_out_max_f_c [ci] =  gamma_out_n
            # reel_in_speed [ci] = gamma_in_n*data['v_w_n']
            # reel_out_speed [ci] = gamma_out_n*data['v_w_n']
            reel_in_speed [ci] = gamma_in_n*i
            reel_out_speed [ci] =  gamma_out_n*i
        # if i <= v_t_n:

        #     P_w = 0.5*data['rho']*i**3
        #     #P_w = 0.5*rho*i**3
        #     for j in gamma_in:
        #         for k in gamma_out:
        #             f_c[cj][ck] = ((np.cos(data['a_elev_out'])-k)**2-((F_in/F_out)*(j**2+2*np.cos(data['a_elev_in'])*j+1)))*((k*j)/(k+j))
        #             ck +=1
        #         cj +=1
        #     max_f_c = np.amax(f_c_mu)
        #     (a,b) = np.where(f_c_mu == max_f_c)
        #     gamma_in_max_f_c [ci] = 
        #     gamma_out_max_f_c [ci] = 
        #     # reel_in_speed [ci] = gamma_in_n*data['v_w_n']
        #     # reel_out_speed [ci] = gamma_out_n*data['v_w_n']
        #     reel_in_speed [ci] = data['gamma_in_n']*i
        #     reel_out_speed [ci] = data['gamma_out_n']*i
        #     cj = 0
        elif i <= v_p_n:

            mu = i/v_t_n
            P_w = 0.5*data['rho']*i**3
            #P_w = 0.5*rho*i**3
            for j in gamma_in:
                f_c_mu[cj] = (data['eff_out']*((1/(mu**2))*(np.cos(data['a_elev_out'])- gamma_out_n)**2)-(((F_in/F_out)*(1+j*np.cos(data['a_elev_in'])*2+j**2))/data['eff_in']))*((j*(mu-1+ gamma_out_n))/(mu*j+mu-1+ gamma_out_n))
                cj +=1
            max_f_c = np.amax(f_c_mu)
            a = np.where(f_c_mu == max_f_c)
            gamma_in_max_f_c[ci] = gamma_in[a]
            gamma_out_max_f_c[ci] = 1-((1- gamma_out_n)/mu)
            P_c [ci] = P_w*max_f_c*F_out*A_proj
            reel_in_speed [ci] = gamma_in_max_f_c[ci]*i
            reel_out_speed [ci] = gamma_out_max_f_c[ci]*i
            cj = 0
            c = ci    
        else:

            #loop to go through all the different gamma_in for all v_w
            #also calculates gamma_out

            mu = i/v_p_n
            P_w = 0.5*data['rho']*i**3
            #F_out_mu = (F_out/(mu**2))*(((1-gamma_out_n)**2)/((1-(gamma_out_n/mu))**2))
            gamma_out_max_f_c[ci] = np.amax(gamma_out_max_f_c)/mu
            for j in gamma_in:
                #f_c_mu[cj] = ((data['eff_out']*((1/((mu**2)))*(np.cos(data['a_elev_out'])- gamma_out_n)**2)-(F_in/F_out)*(1+j*np.cos(data['a_elev_in'])*2+j**2)/data['eff_in'])*(((gamma_out_n)*j)/(gamma_out_n+mu*j)))
                f_c_mu[cj] = ((data['eff_out']*((1/((mu**2)*((v_p_n/v_t_n)**2)))*(np.cos(data['a_elev_out'])- gamma_out_n)**2)-(F_in/F_out)*(1+j*np.cos(data['a_elev_in'])*2+j**2)/data['eff_in'])*((((v_p_n/v_t_n)-1+gamma_out_n)*j)/(gamma_out_n+mu*(v_p_n/v_t_n)*j+(v_p_n/v_t_n)-1)))
                # f_c_mu[cj] = ((1/(mu**2*(v_p_n/v_t_n)))*(np.cos(data['a_elev_out'])-data['gamma_out_n'])**2-(F_in/F_out)*(1+j*np.cos(data['a_elev_in'])*2+j**2))*((((v_p_n/v_t_n)-1+data['gamma_out_n'])*j)/(data['gamma_out_n']+mu*(v_p_n/v_t_n)*j+(v_p_n/v_t_n)-1))
                #f_c_mu[cj] = ((1/(mu**2*(v_p_n/v_t_n)))*(1-gamma_out_p)**2-(F_in/F_out)*(1+j)**2)*(((gamma_out_p)*j)/(gamma_out_p+mu*j)
                cj +=1
            max_f_c = np.amax(f_c_mu)
            a = np.where(f_c_mu == max_f_c)
            gamma_in_max_f_c[ci] = gamma_in[a]
            reel_in_speed [ci] = gamma_in_max_f_c[ci]*i
            reel_out_speed [ci] = gamma_out_max_f_c[ci]*i
            P_c [ci] = P_w*max_f_c*F_out*A_proj
            cj = 0
            d +=1
        ci +=1
    return gamma_in_max_f_c,gamma_out_max_f_c,v_w,P_c,reel_in_speed,reel_out_speed



def plot_three_phase_gamma(step):
    gamma_in_max_f_c,gamma_out_max_f_c,v_w,P_c,reel_in_speed,reel_out_speed = calculate_three_phase(step)
    plt.plot(v_w,gamma_in_max_f_c, label = 'Reel-in gamma')
    plt.plot(v_w,gamma_out_max_f_c,'--', label = 'Reel-out gamma')
    plt.xlabel('Wind speed [m/s]',fontsize = 16)
    plt.ylabel('gamma []',fontsize = 16)
    plt.legend(fontsize = 12)
    plt.grid()
    plt.show()

def plot_three_phase_reel_speeds(step):
    gamma_in_max_f_c,gamma_out_max_f_c,v_w,P_c,reel_in_speed,reel_out_speed = calculate_three_phase(step)
    plt.plot(v_w,reel_in_speed, label = 'Reel-in speed')
    plt.plot(v_w,reel_out_speed,'--', label = 'Reel-out speed')
    plt.xlabel('Wind speed [m/s]',fontsize = 16)
    plt.ylabel('Reeling speed [m/s]',fontsize = 16)
    plt.legend(fontsize = 12)
    plt.grid()
    plt.show()

def plot_three_phase_cycle_power(step):
    gamma_in_max_f_c,gamma_out_max_f_c,v_w,P_c,reel_in_speed,reel_out_speed = calculate_three_phase(step)
    plt.plot(v_w, P_c)
    plt.xlabel('Wind speed [m/s]',fontsize = 16)
    plt.ylabel('Power Output [W]',fontsize = 16)
    plt.legend(fontsize = 12)
    plt.grid()
    plt.show()

plot_three_phase_gamma(1000)
plot_three_phase_reel_speeds(1000)
plot_three_phase_cycle_power(1000)    
