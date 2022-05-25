
import numpy as np 
from Input import *

output = {}
### This function finds the optimal gamma for nominal flight conditions defined in Input.py ###
def calculate_opt_gamma_nominal():
    ## Define gamma range ##
    gamma_in  = np.linspace(1,2,100)
    gamma_out = np.linspace(0,1,100)
    ## Set empty arrays ##
    power_array_m = np.zeros((100,100))
    power_array_e = np.zeros((100,100))
    ## Initiate counters ##
    ci = 0
    cj = 0
    for j in gamma_out:
        for i in gamma_in: 
            
            power_array_m[cj][ci] = P_w*A_proj*(F_out*(1-j)**2-(F_in*(1+i)**2))*((j*i)/(j+i))
            power_array_e[cj][ci] = P_w*A_proj*(eff_out*F_out*(1-j)**2-(F_in*(1+i)**2)/eff_in)*((j*i)/(j+i))
            ci  += 1
    
        cj +=1
        ci = 0 


    ## Find maximal mechanical power  ##     
    max_power_m = np.amax(power_array_m)
    max_power_e = np.amax(power_array_e)
    #(a,b) = np.where(power_array_m == max_power_m)
    (a,b) = np.where(power_array_e == max_power_e)
    print(gamma_out[a],gamma_in[b])
    gamma_out = gamma_out[a][0]
    gamma_in = gamma_in[b][0]
    #print(gamma_out, gamma_in)
    print(max_power_m,max_power_e)
    return gamma_out, gamma_in



### This function calculates the traction forces for nominal flight conditions ###

def calculate_nominal_tractionF():
    gamma_out, gamma_in = calculate_opt_gamma_nominal()

    T_out_n = 0.5*rho*v_w_n**2*A_proj*(1-gamma_out)**2*F_out
    T_in_n = 0.5*rho*v_w_n**2*A_proj*(1+gamma_in)**2*F_in

    return T_out_n, T_in_n

def calculate_nominal_powers():

    T_out_n, T_in_n  = calculate_nominal_tractionF()
    gamma_out, gamma_in = calculate_opt_gamma_nominal()

    P_out = T_out_n*gamma_out*v_w_n
    P_out_e = P_out * eff_out

    P_in = T_in_n*gamma_in*v_w_n
    P_in_e = P_in / eff_in

    ## Sanity check ##
    P_avg_mech = P_out*(gamma_in)/(gamma_in + gamma_out) - P_in*gamma_out/(gamma_in + gamma_out)
    P_avg_elec = P_out_e*(gamma_in)/(gamma_in + gamma_out) - P_in_e*gamma_out/(gamma_in + gamma_out) 
    return T_in_n,T_out_n

def calculated_updated_projected_area():
    
    gamma_out, gamma_in = calculate_opt_gamma_nominal()
    A_proj = P_avg_e_n/P_w/((eff_out*F_out*(1-gamma_out)**2-(F_in*(1+gamma_in)**2)/eff_in)*((gamma_out*gamma_in)/(gamma_out+gamma_in)))
    print(A_proj)


calculate_nominal_powers()
print(calculate_nominal_tractionF())
A_proj = calculated_updated_projected_area()
#calculate_opt_gamma_nominal()

       
        