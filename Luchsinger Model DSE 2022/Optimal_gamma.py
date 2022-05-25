
from InputV2 import *
import numpy as np 
from Input import *
from matplotlib import cm
from matplotlib.colors import ListedColormap,LinearSegmentedColormap
import matplotlib.pyplot as plt


### This function finds the optimal gamma for nominal flight conditions defined in Input.py and plots this on a heat map ###
def calculate_opt_gamma_nominal(data):

    plot_gamma_data = {}
    ## Define gamma range ##
    # Prohibits reel-in speed from exceeding max reeling speed # 
    if max_reel_speed <= 2*data['v_w_n']: 
        lim = data['max_reel_speed']/data['v_w_n']
    else:
        lim = 2
    
    plot_gamma_data['gamma_in']  = np.linspace(0.01,lim,100)
    plot_gamma_data['gamma_out'] = np.linspace(0.01,1,100)

    # gamma_in = np.linspace(1,3,3)
    # gamma_out = np.linspace(1,3,3)

    ## Set empty arrays ##
    plot_gamma_data['power_array_m'] = np.zeros((100,100))
    plot_gamma_data['power_array_e'] = np.zeros((100,100))
    ## Initiate counters ##
    ci = 0
    cj = 0
    for j in plot_gamma_data['gamma_out']:
        for i in plot_gamma_data['gamma_in']: 
            
            plot_gamma_data['power_array_m'][cj][ci] = data['P_w']*data['A_proj'] *(data['F_out']*(1-j)**2-(data['F_in']*(1+i)**2))*((j*i)/(j+i))
            plot_gamma_data['power_array_e'][cj][ci] = data['P_w']*data['A_proj']*(data['eff_out']*data['F_out']*(1-j)**2-(data['F_in']*(1+i)**2)/data['eff_in'])*((j*i)/(j+i))
            ci  += 1
    
        cj +=1
        ci = 0 


    ## Find maximal mechanical power  ##    
     
    data['max_power_m'] = np.amax(plot_gamma_data['power_array_m'])
    data['max_power_e'] = np.amax(plot_gamma_data['power_array_e'])
    #(a,b) = np.where(power_array_m == max_power_m)
    (a,b) = np.where(plot_gamma_data['power_array_e'] == data['max_power_e'])
    
    data['gamma_out_n'] = plot_gamma_data['gamma_out'][a][0]
    data['gamma_in_n'] = plot_gamma_data['gamma_in'][b][0]
    
    #print(gamma_out[a],gamma_in[b])
    #print(data['max_power_m'],data['max_power_e'])
    #np.savetxt('Power.txt',power_array_e)


    
    return data, plot_gamma_data

def plot_gamma_power():
    plot_data = calculate_opt_gamma_nominal(data)[1]
    hsv_modified = cm.get_cmap('hsv', 256)# create new hsv colormaps in range of 0.3 (green) to 0.7 (blue)
    newcmp = ListedColormap(hsv_modified(np.linspace(0.2, 1.0, 256)))# show figure
    plt.pcolormesh(plot_data['gamma_out'], plot_data['gamma_in'],np.transpose(plot_data['power_array_e']), cmap = newcmp)
    plt.colorbar()
    plt.show()

### This function calculates the traction forces for nominal flight conditions ###

def calculate_nominal_tractionF(data):
    
    data['T_out_n'] = 0.5*data['rho']*data['v_w_n']**2*data['A_proj']*(1-data['gamma_out_n'])**2*data['F_out']
    data['T_in_n'] = 0.5*data['rho']*data['v_w_n']**2*data['A_proj']*(1+data['gamma_in_n'])**2*data['F_in']

    return data

def calculate_nominal_powers(data):

    data['P_out'] = data['T_out_n']*data['gamma_out_n']*data['v_w_n']
    data['P_out_e'] = data['P_out'] * data['eff_out']

    data['P_in'] = data['T_in_n']*data['gamma_in_n']*data['v_w_n']
    data['P_in_e'] = data['P_in'] / data['eff_in']

    ## Sanity check, code verification ##
    # P_avg_mech = P_out*(gamma_in)/(gamma_in + gamma_out) - P_in*gamma_out/(gamma_in + gamma_out)
    # P_avg_elec = P_out_e*(gamma_in)/(gamma_in + gamma_out) - P_in_e*gamma_out/(gamma_in + gamma_out) 
    return data

def calculated_updated_projected_area(data):
    
    data['A_proj_u'] =  data['P_avg_e_n']/ data['P_w']/(( data['eff_out']* data['F_out']*(1- data['gamma_out_n'])**2-( data['F_in']*(1+ data['gamma_in_n'])**2)/ data['eff_in'])*(( data['gamma_out_n']* data['gamma_in_n'])/( data['gamma_out_n']+ data['gamma_in_n'])))
    print(data['A_proj_u'])
    return data

def run_nominal_analysis(data):

     
    data = calculate_opt_gamma_nominal(data)[0]
    data = calculate_nominal_tractionF(data)
    data = calculate_nominal_powers(data)
    data = calculated_updated_projected_area(data)
    plot_gamma_power()
    # Write to file #
    file = open("data.txt","w") 
    for key, value in data.items(): 
        file.write('%s:%s\n' % (key, value))
    file.close()

data = get_initial_data()
run_nominal_analysis(data)

    



       
        