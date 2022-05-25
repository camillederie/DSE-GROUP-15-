
import numpy as np 
from Input import *
from matplotlib import cm
from matplotlib.colors import ListedColormap,LinearSegmentedColormap
import matplotlib.pyplot as plt

output = {}
### This function finds the optimal gamma for nominal flight conditions defined in Input.py ###
def calculate_opt_gamma_nominal():
    ## Define gamma range ##
    # Prohibits reel-in speed from exceeding max reeling speed # 
    if max_reel_speed <= 2*v_w_n: 
        lim = max_reel_speed/v_w_n
    else:
        lim = 2
    
    gamma_in  = np.linspace(0.01,lim,100)
    gamma_out = np.linspace(0.01,1,100)

    # gamma_in = np.linspace(1,3,3)
    # gamma_out = np.linspace(1,3,3)

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
    gamma_out_n = gamma_out[a]
    gamma_in_n = gamma_in[b]
    #print(gamma_out, gamma_in)
    print(max_power_m,max_power_e)
    np.savetxt('Power.txt',power_array_e)
    ## Plot ##
    test =[[1,2,3],[4,5,6],[7,8,9]]
    #new_inferno = cm.get_cmap('inferno', 5)# visualize with the new_inferno colormaps
    hsv_modified = cm.get_cmap('hsv', 256)# create new hsv colormaps in range of 0.3 (green) to 0.7 (blue)
    newcmp = ListedColormap(hsv_modified(np.linspace(0.1, 1.0, 256)))# show figure
    plt.pcolormesh(gamma_out, gamma_in,np.transpose(power_array_e), cmap = newcmp)
    #plt.pcolormesh(gamma_in,gamma_out,test, cmap = newcmp)
    plt.colorbar()
    plt.show()

    return gamma_out_n, gamma_in_n

calculate_opt_gamma_nominal()

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
    A_proj_u = P_avg_e_n/P_w/((eff_out*F_out*(1-gamma_out)**2-(F_in*(1+gamma_in)**2)/eff_in)*((gamma_out*gamma_in)/(gamma_out+gamma_in)))
    print(A_proj_u)
    return A_proj_u

# calculate_nominal_powers()
# print(calculate_nominal_tractionF())
# A_proj_u = calculated_updated_projected_area()


       
        