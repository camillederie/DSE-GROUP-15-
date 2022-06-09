import numpy as np

def get_initial_data():
    data = {}
    data['v_w_adj'] = np.linspace(4, 20, 30)
    data['A_proj_list'] = np.linspace(7, 25, 25)
    data['F_out_list'] = np.linspace(10, 150, 20)  # np.linspace(10, 100, 90)

    data['v_w_n'] = 10 #np.linspace(5,15,20) #10 #10.44#np.linspace(5,20,50)
    data['CL_out'] = 1.0720464426912995 #np.linspace(0.6,1.5,50)
    data['A_proj'] = 15.284289430824494 #12.302 #10.04#15.681 #9.723#9.34#21.54#24.19 #8.18 #24.19 #20.79 #21.57 #np.linspace(15,35,50) #19.8 #
    data['T_out_target'] = 10405
    data['T_out_max'] = 10500
    data['rho'] = 1.18
    data['lc'] = 250
    data['CD_out'] = 0.14927169826336975
    #data['CL_in'] = 0.1
    data['CD_in'] = 0.10293127475722369
    data['eff_in'] = 0.639#0.652 #were switched
    data['eff_out'] = 0.652#0.639

    data['P_avg_e_req'] = 20000 #Nominal electrical power (W)
    data['max_reel_speed'] = 25 #m/s

    data['a_elev_out'] = 30*np.pi/180 
    data['a_elev_in'] = 70*np.pi/180

    data['SF_supercap'] =1.2 #Safety factor supercap
    data['diameter_drum'] = 0.46586
    data['drum_circum'] = data['diameter_drum']*np.pi
    data['SF_force'] = 1.35
    data['rpm_min'] = 1500
    data['rpm_n'] = 1525
    data['rpm_max'] = 3000
    data['rpm_motor'] = 750 
    

    ## Intermediate calculations ##

    data['F_out'] = data['CL_out']**3/data['CD_out']**2
    data['F_in'] = data['CD_in'] 
    data['P_w'] = 0.5*data['v_w_n'] **3*data['rho']  # Wind Power

    '''Data for sensitivity analysis: fix some coefficients'''
    data['gamma_out_n'] = 0.3926
    data['gamma_in_n'] = 1.8973

    return data 