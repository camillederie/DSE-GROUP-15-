#Optimisation of the Main System Parameters
import numpy as np

from Drum_and_tether_design import structures_calculation
from Luchsinger.luchsingermodel.Nominal_power_cycle import  run_nominal_analysis
from Luchsinger.luchsingermodel.InputV2 import get_initial_data
from Aero.aero_main_function import main_aero_function

A_proj = 16.65
Points = 1000000
Kite_segments = 12
N_split = 5
AoA_range_out = np.arange(8, 15.5, 0.5)
AoA_range_in = np.arange(-1, 3.5, 0.5)
TAS = 32.5



A_proj_last = 0
c=0
area_diff = 5
while abs(A_proj_last-A_proj) > area_diff:
    A_proj_last = A_proj
    c+= 1
    print('Hey Bradda, the optimisation has started. Leggo! This is cycle: ',c)
    ## AERO ##
    print('Aero started')
    CL_average_out, CD_average_out, CL3_CD2_average_out, CD_average_in, A_proj, Strut_area_av = main_aero_function(A_proj_last, Points, Kite_segments, N_split, AoA_range_out, AoA_range_in, TAS, Print=False)
    # CL_average_in, CD_average_in, CL3_CD2_average_in, A_proj, Strut_area_av = main_aero_function(A_proj, Points, Kite_segments, N_split, AoA_range_in, TAS, Print=False)
    print('Aero part finished')
    ## POWER ##
    data = get_initial_data()
    data['F_out'] = CL3_CD2_average_out
    data['F_in'] = CD_average_in
    data['CL_out'] = CL_average_out
    data['CD_out'] = CD_average_out
    # data['CL_in'] = CL_average_in
    data['CD_in'] = CD_average_in
    data =  run_nominal_analysis(data) 
    A_proj = data['A_proj']
    T_F_out = data['T_out_elev_n']
    TAS = data['v_a_out']
    gamma_out = data['gamma_out_n']
    gamma_in = data['gamma_in_n']
    print('Power part finished, the old area was: ', A_proj_last, 'The new one is:', A_proj)
    ## STRUCTURES ##

print('Optimisation completed!')
# Write to file #
file = open("data_optim.txt","w") 
for key, value in data.items(): 
    file.write('%s:%s\n' % (key, value))
file.close()
print('The extended results of the analysis can be found in the data file added to the directory.')
    
#DEFINITIONS
len_drum_in = 1.10
extra_len_in = 75
nom_load_in = data['T_out_elev_n']
saf_fac_in = 2
kite_mass_margin = 1.05
t_out = data['t_out']
t_in = data['t_in']
structures_calculation(A_proj, Strut_area_av, len_drum_in, data['a_elev_out'], extra_len_in, nom_load_in, saf_fac_in, kite_mass_margin, data['t_out'], data['t_in'])
