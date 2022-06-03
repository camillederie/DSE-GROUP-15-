#Optimisation of the Main System Parameters
import numpy as np

from Drum_and_tether_design import
from Nominal_power_cycle import  run_nominal_analysis, get_initial_data 
from aero_main_function import main_aero_function



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
data =  run_nominal_analysis(data) 
A_proj = data['A_proj']
T_F_out = data['T_out_elev_n']