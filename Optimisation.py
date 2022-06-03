#Optimisation of the Main System Parameters

from Drum_and_tether_design import
from Nominal_power_cycle
from aero_main_function import main_aero_function



A_proj = 16.65
Points = 10000
Kite_segments = 12
N_split = 5
AoA_range = np.arange(8, 10.5, 0.5)
TAS = 32.5

CL_average, CD_average, CL3_CD2_average, A_proj, Strut_area_av = main_aero_function(A_proj, Points, Kite_segments, N_split, AoA_range, TAS, Print=False)