
#from Luchsinger.luchsingermodel.InputV2 import *
from InputV2 import *
import numpy as np 
from matplotlib import cm
from matplotlib.colors import ListedColormap,LinearSegmentedColormap
import matplotlib.pyplot as plt

#Luchsinger.luchsingermodel.
### This function finds the optimal gamma for nominal flight conditions defined in Input.py and plots this on a heat map ###
def calculate_opt_gamma_nominal(data):

    plot_gamma_data = {}
    # Define gamma range ##
    #Prohibits reel-in speed from exceeding max reeling speed # 
    if data['max_reel_speed'] <= 2*data['v_w_n']: 
        lim = data['max_reel_speed']/data['v_w_n']
    else:
        lim = 2
    lim = 2.0
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
     
    #data['max_power_m'] = np.amax(plot_gamma_data['power_array_m'])
    data['P_elec_opt_gamma'] = np.amax(plot_gamma_data['power_array_e'])
    (a,b) = np.where(plot_gamma_data['power_array_e'] == data['P_elec_opt_gamma'])
    
    data['gamma_out_n'] = plot_gamma_data['gamma_out'][a][0]
    data['gamma_in_n'] = plot_gamma_data['gamma_in'][b][0]
    
    return data, plot_gamma_data


def calculate_opt_gamma_nominal_elev(data):

        plot_gamma_data = {}
        # Define gamma range ##
        #Prohibits reel-in speed from exceeding max reeling speed # 
        if data['max_reel_speed'] <= 2*data['v_w_n']:
            lim = data['max_reel_speed']/data['v_w_n']
        else:
            lim = 2.5
        
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
                
                #plot_gamma_data['power_array_m'][cj][ci] = data['P_w']*data['A_proj']*(data['F_out']*(np.cos(data['a_elev_out'])-j)**2-(data['F_in']*(i**2+2*np.cos(data['a_elev_in'])*i+1)))*((j*i)/(j+i))
                plot_gamma_data['power_array_e'][cj][ci] = data['A_proj']*(data['eff_out']*data['F_out']*(np.cos(data['a_elev_out'])-j)**2-(data['F_in']*(i**2+2*np.cos(data['a_elev_in'])*i+1))/data['eff_in'])*((j*i)/(j+i))
                ci  += 1
        #data['P_w']*
            cj +=1
            ci = 0 


        ## Find maximal mechanical power  ##    
        
        #data['max_power_m'] = np.amax(plot_gamma_data['power_array_m'])
        data['P_elec_opt_gamma'] = np.amax(plot_gamma_data['power_array_e'])
        #(a,b) = np.where(power_array_m == max_power_m)
        (a,b) = np.where(plot_gamma_data['power_array_e'] == data['P_elec_opt_gamma'])
        
        data['gamma_out_n'] = plot_gamma_data['gamma_out'][a][0]
        data['gamma_in_n'] = plot_gamma_data['gamma_in'][b][0]

    
        #data['max_cycle_power'] = data['P_w']*data['A_proj']*data['F_out']*(np.cos(data['a_elev_out']))**3*4/27
        #print(gamma_out[a],gamma_in[b])

        #print(data['max_power_m'],data['max_power_e'])
        #np.savetxt('Power.txt',power_array_e)
    
        return data, plot_gamma_data

### This function finds the optimal gamma in when a fixed gamma out is set ###

def calculate_opt_gamma_in(data):

    plot_gamma_data = {}
    ## Define gamma range ##
    #Prohibits reel-in speed from exceeding max reeling speed # 
    if data['max_reel_speed'] <= 2*data['v_w_n']:
        lim = data['max_reel_speed']/data['v_w_n']
    else:
        lim = 2.5
    
    plot_gamma_data['gamma_in']  = np.linspace(0.01,lim,100)
    #plot_gamma_data['gamma_out'] = np.linspace(0.01,1,100)

    # gamma_in = np.linspace(1,3,3)
    # gamma_out = np.linspace(1,3,3)

    ## Set empty arrays ##
    plot_gamma_data['power_array_m'] = np.zeros(100)
    plot_gamma_data['power_array_e'] = np.zeros(100)
    ## Initiate counters ##
    ci = 0
    cj = 0
    
    for i in plot_gamma_data['gamma_in']: 
            
        #plot_gamma_data['power_array_m'][cj][ci] = data['P_w']*data['A_proj']*(data['F_out']*(np.cos(data['a_elev_out'])-j)**2-(data['F_in']*(i**2+2*np.cos(data['a_elev_in'])*i+1)))*((j*i)/(j+i))
        plot_gamma_data['power_array_e'][ci] = data['P_w']*data['A_proj']*(data['eff_out']*data['F_out']*(np.cos(data['a_elev_out'])-data['gamma_out_n'])**2-(data['F_in']*(i**2+2*np.cos(data['a_elev_in'])*i+1))/data['eff_in'])*((data['gamma_out_n']*i)/(data['gamma_out_n']+i))
        ci  += 1

    ## Find maximal mechanical power  ##    
     
    #data['max_power_m'] = np.amax(plot_gamma_data['power_array_m'])
    data['P_elec_opt_gamma_in'] = np.amax(plot_gamma_data['power_array_e'])
    #(a,b) = np.where(power_array_m == max_power_m)
    b = np.where(plot_gamma_data['power_array_e'] == data['P_elec_opt_gamma_in'])
    
    #data['gamma_out_n'] = plot_gamma_data['gamma_out'][a][0]
    data['gamma_in_n'] = plot_gamma_data['gamma_in'][b][0]
    
    return data

### This function plots the Power for gamma in and out combinations ###

def plot_gamma_power(data_plot):
    plot_data = data_plot
    hsv_modified = cm.get_cmap('hsv', 256)# create new hsv colormaps in range of 0.3 (green) to 0.7 (blue)
    newcmp = ListedColormap(hsv_modified(np.linspace(0.2, 1.0, 256)))# show figure
    plt.pcolormesh(plot_data['gamma_out'], plot_data['gamma_in'],np.transpose(plot_data['power_array_e']), cmap = newcmp, shading='auto')
    cbar = plt.colorbar()
    plt.xlabel(r'$\gamma_{out}$')
    plt.ylabel(r'$\gamma_{in}$')
    #plt.colorbar()
    cbar.set_label('Average Output Power')
    
    # plt.show()

### This function calculates the traction forces for nominal flight conditions ###

def calculate_nominal_tractionF(data):
    
    data['T_out_n'] = 0.5*data['rho']*data['v_w_n']**2*data['A_proj']*(1-data['gamma_out_n'])**2*data['F_out']
    data['T_out_elev_n'] = 0.5*data['rho']*data['v_w_n']**2*data['A_proj']*(np.cos(data['a_elev_out'])-data['gamma_out_n'])**2*data['F_out']

    data['T_in_n'] = 0.5*data['rho']*data['v_w_n']**2*data['A_proj']*(1+data['gamma_in_n'])**2*data['F_in']
    data['T_in_elev_n'] =  0.5*data['rho']*data['v_w_n']**2*data['A_proj']*(1+2*data['gamma_in_n']*np.cos(data['a_elev_in'])+data['gamma_in_n']**2)*data['F_in']

    return data

### This function calculates the powers during all phases of the power cycle ###

def calculate_nominal_powers(data):

    data['P_out'] = data['T_out_n']*data['gamma_out_n']*data['v_w_n']
    data['P_out_e'] = data['P_out'] * data['eff_out']

    data['P_out_elev'] = data['T_out_elev_n']*data['gamma_out_n']*data['v_w_n']
    data['P_out_e_elev'] = data['P_out_elev'] * data['eff_out']

    data['P_in'] = data['T_in_n']*data['gamma_in_n']*data['v_w_n']
    data['P_in_e'] = data['P_in'] / data['eff_in']

    data['P_in_elev'] = data['T_in_elev_n']*data['gamma_in_n']*data['v_w_n']
    data['P_in_e_elev'] = data['P_in_elev'] / data['eff_in']

    data['P_avg_mech']= data['P_w']*data['A_proj']*(data['F_out']*(1-data['gamma_out_n'])**2-(data['F_in']*(data['gamma_in_n']+1)**2))*(( data['gamma_out_n']* data['gamma_in_n'])/( data['gamma_out_n']+ data['gamma_in_n']))
    data['P_avg_elec']= data['P_w']*data['A_proj']*(data['eff_out']*data['F_out']*(1-data['gamma_out_n'])**2-(data['F_in']*(data['gamma_in_n']+1)**2)/data['eff_in'])*(( data['gamma_out_n']* data['gamma_in_n'])/( data['gamma_out_n']+ data['gamma_in_n']))

    data['P_avg_mech_elev']= data['P_w']*data['A_proj']*(data['F_out']*(np.cos(data['a_elev_out'])-data['gamma_out_n'])**2-(data['F_in']*(data['gamma_in_n']**2+2*np.cos(data['a_elev_in'])*data['gamma_in_n']+1)))*(( data['gamma_out_n']* data['gamma_in_n'])/( data['gamma_out_n']+ data['gamma_in_n']))
    data['P_avg_elec_elev']= data['P_w']*data['A_proj']*(data['eff_out']*data['F_out']*(np.cos(data['a_elev_out'])-data['gamma_out_n'])**2-(data['F_in']*(data['gamma_in_n']**2+2*np.cos(data['a_elev_in'])*data['gamma_in_n']+1))/data['eff_in'])*(( data['gamma_out_n']* data['gamma_in_n'])/( data['gamma_out_n']+ data['gamma_in_n']))

    data['max_cycle_power_elev'] = data['P_w']*data['A_proj']*data['F_out']*(np.cos(data['a_elev_out']))**3*4/27

    ## Sanity check, code verification ##
    data['P_avg_mech_verif'] = data['P_out']*(data['gamma_in_n'])/(data['gamma_in_n'] + data['gamma_out_n']) - data['P_in']*data['gamma_out_n']/(data['gamma_in_n'] + data['gamma_out_n'])
    data['P_avg_elec_verif'] = data['P_out_e']*(data['gamma_in_n'])/(data['gamma_in_n'] + data['gamma_out_n']) - data['P_in_e']*data['gamma_out_n']/(data['gamma_in_n'] + data['gamma_out_n'])
    
    data['P_avg_mech_elev_verif'] = data['P_out_elev']*(data['gamma_in_n'])/(data['gamma_in_n'] + data['gamma_out_n']) - data['P_in_elev']*data['gamma_out_n']/(data['gamma_in_n'] + data['gamma_out_n'])
    data['P_avg_elec_elev_verif'] = data['P_out_e_elev']*(data['gamma_in_n'])/(data['gamma_in_n'] + data['gamma_out_n']) - data['P_in_e_elev']*data['gamma_out_n']/(data['gamma_in_n'] + data['gamma_out_n'])
    
    return data

### This function updates the projected area based on the average power requirement ###

def calculate_updated_projected_area(data):
    
    data['A_proj_u'] = data['P_avg_e_req']/ data['P_w']/((data['eff_out']*data['F_out']*(np.cos(data['a_elev_out'])-data['gamma_out_n'])**2-(data['F_in']*(data['gamma_in_n']**2+2*np.cos(data['a_elev_in'])*data['gamma_in_n']+1))/data['eff_in'])*(( data['gamma_out_n']* data['gamma_in_n'])/( data['gamma_out_n']+ data['gamma_in_n'])))
    #data['A_proj']= data['A_proj_u']
    return data

### This function calculates the cycle times ###

def calculate_cycle_param(data):

    data['cycle_time'] = (data['lc']/data['v_w_n'])*((data['gamma_out_n']+data['gamma_in_n'])/(data['gamma_out_n']*data['gamma_in_n']))
    data['t_in'] = data['cycle_time']*data['gamma_out_n']/(data['gamma_in_n'] + data['gamma_out_n'])
    data['t_out'] = data['cycle_time']*data['gamma_in_n']/(data['gamma_in_n'] + data['gamma_out_n'])

    return data

### This function looks for gamma outs to lower the traction force, and directly updates the area accordingly
def evaluate_tether_force(data): 
     
    data['gamma_out_n_init'] = data['gamma_out_n']
    data['gamma_in_n_init'] = data['gamma_in_n']
    data['A_proj_init'] = data['A_proj']

    TF_an = {'gamma_out':[],'force': [],'area':[]}
    while data['gamma_out_n'] < .6:
        data['gamma_out_n'] += 0.01
        data = calculate_opt_gamma_in(data)
        data = calculate_updated_projected_area(data)
        data['A_proj'] = data['A_proj_u']

        data = calculate_nominal_tractionF(data)
        
        TF_an['gamma_out'].append(data['gamma_out_n'])
        TF_an['force'].append(data['T_out_elev_n'])
        TF_an['area'].append(data['A_proj'])

    #plot_TF_an(TF_an)
    # data['gamma_out_n'] 
    # data['gamma_out_n'] = 0.43#float(input('Enter the chosen gamma reel-out to find the correspinding optimal gamma reel-in: '))
    # data = calculate_opt_gamma_in(data)
    # data = calculate_updated_projected_area(data)
    #data['A_proj'] = data['A_proj_u']
    data['A_proj'] =  data['A_proj_init']
    data['gamma_out_n'] = np.cos(data['a_elev_out']) - np.sqrt(data['T_out_target']*2/(data['rho']*data['v_w_n']**2*data['A_proj']*data['F_out']))
    data = calculate_opt_gamma_in(data)
    data = calculate_updated_projected_area(data)
    if data['A_proj']< data['A_proj_u']:
        print('ok')

    else:
        print('not ok')
    print(data['A_proj']) 
    print(data['A_proj_u'])
    data['A_proj'] = data['A_proj_u']
    data = calculate_nominal_tractionF(data)
    data = calculate_nominal_powers(data)
    
    return data

def evaluate_adj_wind_areafix(data,v_w_adj):
    data['gamma_out_v_w_adj'] = np.cos(data['a_elev_out']) - np.sqrt(data['T_out_elev_n']/(0.5*data['rho']*v_w_adj**2*data['A_proj']*data['F_out']))
    return data
         
def plot_TF_an(TF_an):
    fig, ax1 = plt.subplots()
    ax2 = ax1.twiny()
    ax1.plot(TF_an['gamma_out'],TF_an['force'], color = 'r')
    ax2.plot(TF_an['area'],TF_an['force'], color = 'b')
    ax1.set_ylabel('Tether Force (N)')
    ax1.set_xlabel('Gamma Reel-out', color = 'r')
    ax2.set_xlabel('Projected Area (m2)', color = 'b')
    plt.grid()
    # plt.show()
def three_phase_an(data):
    datasens= {}
    datasens['F_out_list'] = data['F_out_list']
    datasens['A_proj_list'] = data['A_proj_list']
    datasens['v_w_list'] = data['v_w_adj']

    datasens['T_out_list_VW'] =[]
    datasens['P_avg_e_list_VW'] = []
    datasens['gamma_out_list_VW'] =[]
    datasens['gamma_in_list_VW'] = []
    datasens['cycle_time_list_VW'] = []
    datasens['supercap_list_VW'] = []
    datasens['elevation_angle_VW'] = []
    datasens['Power_reel_out_list_VW'] = []
    data['a_elev_init'] = data['a_elev_out']
    for w in data['v_w_adj']:

        data['v_w_n'] = float(w)
        data['P_w'] = 0.5*data['v_w_n'] **3*data['rho'] 
        data['a_elev_out'] = data['a_elev_init']
        data = calculate_opt_gamma_nominal_elev(data)[0]
        data = calculate_nominal_tractionF(data)
        data = calculate_nominal_powers(data)
        if data['T_out_elev_n'] > data['T_out_max']:
            
            data['gamma_out_n'] = np.cos(data['a_elev_out']) - np.sqrt(data['T_out_target']*2/(data['rho']*data['v_w_n']**2*data['A_proj']*data['F_out']))
            data = calculate_opt_gamma_in(data)
            data = calculate_nominal_tractionF(data)
            data = calculate_nominal_powers(data)
            if data['P_out_e_elev']>data['Generator_lim']:
                
                print('Generator limit reached.')
                data['rpm_out'] = data['v_w_n']*data['gamma_out_n']/(data['drum_circum'])*60/0.105
                print(data['rpm_out'])
                data['gamma_out_n'] = 3000*(data['drum_circum'])/60/data['v_w_n']*0.105
                data['a_elev_out_new'] = np.arccos( data['gamma_out_n'] + np.sqrt(data['T_out_target']*2/(data['rho']*data['v_w_n']**2*data['A_proj']*data['F_out'])))
                print(data['a_elev_out_new']*180/np.pi)
                
                data['a_elev_out'] = data['a_elev_out_new']
                
                data = calculate_opt_gamma_in(data)
                data = calculate_nominal_tractionF(data)
                data = calculate_nominal_powers(data)
                if data['a_elev_out']<data['a_elev_init']:
                    data['a_elev_out'] = data['a_elev_init']
        
        data = calculate_cycle_param(data)
        data['a_elev_out'] =  data['a_elev_out']*180/np.pi
        datasens['T_out_list_VW'].append(data['T_out_elev_n'])
        datasens['P_avg_e_list_VW'].append(data['P_avg_elec_elev'])
        datasens['gamma_out_list_VW'].append(data['gamma_out_n'])
        datasens['gamma_in_list_VW'].append(data['gamma_in_n'])
        datasens['cycle_time_list_VW'].append(data['cycle_time'])
        
        datasens['elevation_angle_VW'].append(data['a_elev_out'])
        datasens['Power_reel_out_list_VW'].append(data['P_out_e_elev'])
                
    plt.plot(datasens["v_w_list"],datasens['P_avg_e_list_VW'])
    plt.grid()
    plt.show()
    plt.plot(datasens["v_w_list"],datasens['T_out_list_VW'])
    plt.grid()
    plt.show()
    plt.plot(datasens["v_w_list"],datasens['elevation_angle_VW'])
    plt.grid()
    plt.show()
        
    #data = size_supercap(data)
    return data
### This function calculates the apparent and kite cross wind speed required ###

def calculate_apparent_speed(data):
    v_out = data['v_w_n']*data['gamma_out_n']
    v_in  = data['v_w_n']*data['gamma_in_n']
    v_w   = data['v_w_n']
    
    data['v_a_in'] = np.sqrt(v_w**2+2*v_w*v_in*np.cos(data['a_elev_in'])+v_in**2)
    data['v_a_out']= np.sqrt(data['T_out_elev_n']/(0.5*data['rho']*np.sqrt(data['CL_out']**2+data['CD_out']**2)*data['A_proj']))
    data['v_kc'] = np.sqrt(-(v_w**2-2*v_w*v_out*np.cos(data['a_elev_out'])+v_out**2-data['v_a_out']**2))
    return data

### Size Power Components ###
def size_supercap(data):

    data['E_out'] = data['P_out_e_elev']*data['t_out'] *0.000277777778
    data['E_in'] = data['P_in_e_elev']*data['t_in'] *0.000277777778
    data['SC_cap'] = data['E_in']*data['SF_supercap']  #wh
    return data

def size_generator(data):
    data['rpm_n_out'] = data['v_w_n']*data['gamma_out_n']/(data['drum_circum'])*60
    data['rpm_n_in'] = data['v_w_n']*data['gamma_in_n']/(data['drum_circum'])*60
    # print(data['rpm_n_out'],data['rpm_n_in'])
    data['GR_min'] = data['rpm_n_out']/data['rpm_min']
    data['GR_n'] = data['rpm_n_out']/data['rpm_n']
    data['GR_max'] = data['rpm_n_out']/data['rpm_max']
    data['GR_motor'] = data['rpm_n_in']/data['rpm_motor']


    data = evaluate_adj_wind_areafix(data,data['v_w_adj'])
    v_r_out = data['gamma_out_v_w_adj']*data['v_w_adj']
    rpm_max = np.ones(len(v_r_out))*data['rpm_max']
    
    rpm = (v_r_out/(data['drum_circum'])*60)/data['GR_n']
    data['P_out_e_adj_wind'] = (((rpm+436.25)/69.124))

    fig, ax1 = plt.subplots()
    ax2 = ax1.twinx()

    ax1.plot(data['v_w_adj'],data['P_out_e_adj_wind'], color = 'g')
    ax2.plot(data['v_w_adj'],rpm_max, color = 'r', label = 'Maximal rotational speed')
    ax2.plot(data['v_w_adj'],rpm, color = 'b', label = 'Rotational speed and output power vs wind speed')
    ax2.set_ylabel('Rotational speed of generator (rpm)')
    ax1.set_xlabel('Wind speed (m/s)')
    ax1.set_ylabel('Electrical output power (kW)')
    plt.legend()
    plt.grid()
    
    plt.show()

    return data
    
### Do Analysis ###
def run_nominal_analysis(data):
    ip = 1
    #ip = int(input('Enter 1 if you want to take the tether elevation into account for finding the optimal reeling speeds, 0 for ignoring it: '))
    if ip == 0:
        data = calculate_opt_gamma_nominal(data)[0]
        data_plot = calculate_opt_gamma_nominal(data)[1]
    elif ip == 1:
        data = calculate_opt_gamma_nominal_elev(data)[0]
        data_plot = calculate_opt_gamma_nominal_elev(data)[1]
    else: 
        print('Enter a valid number!')
        quit()
    data = calculate_nominal_tractionF(data)
    data = calculate_nominal_powers(data)
    data = calculate_updated_projected_area(data)
    if abs(data['A_proj_u']-data['A_proj'])>0.01:
        print('The projected area of the kite should be iterated on!')
    else: 
        print('The area of the kite is optimal for the required power output.')
    #data['A_proj'] = data['A_proj_u']
    data = calculate_nominal_tractionF(data)
    data = calculate_nominal_powers(data)
    data = calculate_cycle_param(data)
    data = calculate_apparent_speed(data)
    data = size_supercap(data)
    
    #plot_gamma_power(data_plot)


    # Write to file #
    file = open("Luchsinger\data.txt","w") 
    for key, value in data.items(): 
        file.write('%s:%s\n' % (key, value))
    file.close()
    print('The extended results of the analysis can be found in the data file added to the directory.')
    #ip2 = input('The tetherforce is now', data['T_out_elev_n'], 'enter 1 if you want to analyse how to lower it, else enter 0: ')
    ip2 = 1
    if ip2 == 0:
        quit()
    elif ip2 == 1:
        run_TF_anal(data)
    
    else: 
        print('Enter a valid number!')
        quit()
    
    return data

def run_TF_anal(data):
    data = evaluate_tether_force(data)
    data = calculate_cycle_param(data)
    data = calculate_apparent_speed(data)
    data = size_supercap(data)
    data = size_generator(data)
   
    # Write to file #
    file = open("Luchsinger\data_TF_an.txt","w") 
    for key, value in data.items(): 
        file.write('%s:%s\n' % (key, value))
    file.close()
    print('The extended results of the analysis can be found in the data file added to the directory.')

def sensitivity_analysis(data):
    datasens ={}
    datasens['F_out_list'] = data['F_out_list']
    datasens['A_proj_list'] = data['A_proj_list']
    datasens['v_w_list'] = data['v_w_adj']

    datasens['T_out_list_VW'] =[]
    datasens['P_avg_e_list_VW'] = []
    datasens['gamma_out_list_VW'] =[]
    datasens['gamma_in_list_VW'] = []
    datasens['cycle_time_list_VW'] = []
    datasens['supercap_list_VW'] = []
    datasens['Power_reel_out_list_VW'] = []

    datasens['T_out_list_A'] = []
    datasens['P_avg_e_list_A'] = []
    datasens['supercap_list_A'] = []
    datasens['Power_reel_out_list_A'] = []

    datasens['T_out_list_FO'] =[]
    datasens['P_avg_e_list_FO'] = []
    datasens['gamma_out_list_FO'] =[]
    datasens['gamma_in_list_FO'] = []
    datasens['cycle_time_list_FO'] = []
    datasens['supercap_list_FO'] = []
    datasens['Power_reel_out_list_FO'] = []

   
    for w in data['v_w_adj']:
        data['v_w_n'] = float(w)
        data['P_w'] = 0.5*data['v_w_n'] **3*data['rho'] 
       
        data = calculate_opt_gamma_nominal_elev(data)[0]
        data = calculate_nominal_tractionF(data)
        if data['T_out_elev_n'] > data['T_out_max']:
            data['gamma_out_n'] = np.cos(data['a_elev_out']) - np.sqrt(data['T_out_target']*2/(data['rho']*data['v_w_n']**2*data['A_proj']*data['F_out']))
            data = calculate_opt_gamma_in(data)
            data = calculate_nominal_tractionF(data)
        data = calculate_nominal_powers(data)
        data = calculate_cycle_param(data)
        data = size_supercap(data)

        datasens['T_out_list_VW'].append(data['T_out_elev_n'])
        datasens['P_avg_e_list_VW'].append(data['P_avg_elec_elev'])
        datasens['gamma_out_list_VW'].append(data['gamma_out_n'])
        datasens['gamma_in_list_VW'].append(data['gamma_in_n'])
        datasens['cycle_time_list_VW'].append(data['cycle_time'])
        datasens['supercap_list_VW'].append(data['SC_cap'])
        datasens['Power_reel_out_list_VW'].append(data['P_out_e_elev'])

    data = get_initial_data()
    
    for a in data['A_proj_list']: 
        data['A_proj'] = a
        data = calculate_nominal_tractionF(data)
        data = calculate_nominal_powers(data)
        data = size_supercap(data)
        datasens['T_out_list_A'].append(data['T_out_elev_n'])
        datasens['P_avg_e_list_A'].append(data['P_avg_elec_elev'])
        datasens['supercap_list_A'].append(data['SC_cap'])
        datasens['Power_reel_out_list_VW'].append(data['P_out_e_elev'])
        

    data = get_initial_data()
    
    for c in data['F_out_list']:
        data['F_out'] = c
        data = calculate_opt_gamma_nominal_elev(data)[0]
        data = calculate_nominal_tractionF(data)
        if data['T_out_elev_n'] > data['T_out_max']:
            data['gamma_out_n'] = np.cos(data['a_elev_out']) - np.sqrt(data['T_out_target']*2/(data['rho']*data['v_w_n']**2*data['A_proj']*data['F_out']))
            data = calculate_opt_gamma_in(data)
            data['gamma_out_n'] = np.cos(data['a_elev_out']) - np.sqrt(data['T_out_target']*2/(data['rho']*data['v_w_n']**2*data['A_proj']*data['F_out']))
    
        # data = calculate_updated_projected_area(data)
        # data['A_proj'] = data['A_proj_u']
        data = calculate_nominal_tractionF(data)
        data = calculate_nominal_powers(data)
        data = calculate_cycle_param(data)
        data = size_supercap(data)
        datasens['T_out_list_FO'].append(data['T_out_elev_n'])
        datasens['P_avg_e_list_FO'].append(data['P_avg_elec_elev'])
        datasens['gamma_out_list_FO'].append(data['gamma_out_n'])
        datasens['gamma_in_list_FO'].append(data['gamma_in_n'])
        datasens['cycle_time_list_FO'].append(data['cycle_time'])
        datasens['supercap_list_FO'].append(data['SC_cap'])
        datasens['Power_reel_out_list_VW'].append(data['P_out_e_elev'])
    print(len(datasens['T_out_list_FO']),len(datasens['F_out_list']))
    # file = open("Luchsinger\datasens.txt","w")
    # for key, value in datasens.items():
    #     file.write('%s:%s\n' % (key, value))
    # file.close()
    # print('The extended results of the analysis can be found in the data file added to the directory.')
    plt.plot(datasens["v_w_list"],datasens['gamma_out_list_VW'],color = 'r')
    plt.plot(datasens['v_w_list'],datasens['gamma_in_list_VW'],color = 'b')
    # plt.show()
    plt.plot(datasens["v_w_list"],datasens['P_avg_e_list_VW'])
    # plt.show()
    plt.plot(datasens['F_out_list'],datasens['T_out_list_FO'])
    # plt.show()

    datasens['T_out_list_VW'] =    np.array( datasens['T_out_list_VW'] )
    datasens['P_avg_e_list_VW'] =   np.array( datasens['P_avg_e_list_VW'] )
    datasens['gamma_out_list_VW'] = np.array( datasens['gamma_out_list_VW'] )
    datasens['gamma_in_list_VW'] =  np.array( datasens['gamma_in_list_VW'] )
    datasens['cycle_time_list_VW'] =  np.array( datasens['cycle_time_list_VW'] )
    datasens['supercap_list_VW'] = np.array(datasens['supercap_list_VW'] )
    datasens['Power_reel_out_list_VW'] = np.array( datasens['Power_reel_out_list_VW'] )

    datasens['T_out_list_A'] =  np.array( datasens['T_out_list_A'] )
    datasens['P_avg_e_list_A'] =  np.array( datasens['P_avg_e_list_A'] )
    datasens['supercap_list_A'] = np.array(datasens['supercap_list_A'] )
    datasens['Power_reel_out_list_A'] = np.array( datasens['Power_reel_out_list_A'] )

    datasens['T_out_list_FO'] =  np.array( datasens['T_out_list_FO'] )
    datasens['P_avg_e_list_FO'] =  np.array( datasens['P_avg_e_list_FO'] )
    datasens['gamma_out_list_FO'] =  np.array( datasens['gamma_out_list_FO'] )
    datasens['gamma_in_list_FO'] =  np.array( datasens['gamma_in_list_FO'] )
    datasens['cycle_time_list_FO'] =  np.array( datasens['cycle_time_list_FO'] )
    datasens['supercap_list_FO'] = np.array(datasens['supercap_list_FO'] )
    datasens['Power_reel_out_list_FO'] = np.array( datasens['Power_reel_out_list_FO'] )

    return datasens

def sanety_check(data):
    data = calculate_nominal_tractionF(data)
    data = calculate_nominal_powers(data)
    return data
#data = get_initial_data()

#data = sensitivity_analysis(get_initial_data())
#data = run_nominal_analysis(get_initial_data()) 
#print(sanety_check(get_initial_data()))
data = three_phase_an(get_initial_data())




    



       
        