import numpy as np

def get_initial_data():
    data = {}

    data['v_w_n'] = 10 #np.linspace(5,20,50)
    data['CL_out'] = 1.0 #np.linspace(0.6,1.5,50)
    data['A_proj'] = 21.952#24.19 #8.18 #24.19 #20.79 #21.57 #np.linspace(15,35,50) #19.8 #

    data['rho'] = 1.18
    data['lc'] = 250
    data['CD_out'] = 0.2
    data['CL_in'] = 0.14
    data['CD_in'] = 0.07
    data['eff_in'] = 0.652
    data['eff_out'] = 0.639

    data['P_avg_e_req'] = 20000 #Nominal electrical power (W)
    data['max_reel_speed'] = 25 #m/s

    data['a_elev_out'] = 20*np.pi/180 
    data['a_elev_in'] = 70*np.pi/180
    data['a_elev_in'] = 70*np.pi/180
    ## Intermediate calculations

    data['F_out'] = data['CL_out']**3/data['CD_out']**2
    data['F_in'] = data['CD_in'] 
    data['P_w'] = 0.5*data['v_w_n'] **3*data['rho']  # Wind Power 
    return data 