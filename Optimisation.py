#Optimisation of the Main System Parameters
import numpy as np

from Drum_and_tether_design import *
from Luchsinger.luchsingermodel.Nominal_power_cycle import  run_nominal_analysis
from Luchsinger.luchsingermodel.InputV2 import get_initial_data
from Aero.aero_main_function import main_aero_function

A_proj = 16.65
Points = 10000
Kite_segments = 12
N_split = 5
AoA_range_out = np.arange(8, 15.5, 0.5)
AoA_range_in = np.arange(-1, 3.5, 0.5)
TAS = 32.5

CL_average_out, CD_average_out, CL3_CD2_average_out, A_proj, Strut_area_av = main_aero_function(A_proj, Points, Kite_segments, N_split, AoA_range_out, TAS, Print=False)
CL_average_in, CD_average_in, CL3_CD2_average_in, A_proj, Strut_area_av = main_aero_function(A_proj, Points, Kite_segments, N_split, AoA_range_in, TAS, Print=False)

data = get_initial_data()
data['F_out'] = CL3_CD2_average_out
data['F_in'] = CD_average_in
data['CL_out'] = CL_average_out
data['CD_out'] = CD_average_out
data['CL_in'] = CL_average_in
data['CD_in'] = CD_average_in
data =  run_nominal_analysis(data) 
A_proj = data['A_proj']
T_F_out = data['T_out_elev_n']


#DEFINITIONS
kite_area_in = A_proj
avg_strut_in = 0.6867
len_drum_in = 1.10
angle_in = 30
extra_len_in = 75
nom_load_in = 10329.316
saf_fac_in = 2
kite_mass_margin = 1.05
t_out = 103.5
t_in = 11.5

structures_calculation(kite_area_in, avg_strut_in, len_drum_in, angle_in, extra_len_in, nom_load_in, saf_fac_in, kite_mass_margin, t_out, t_in)