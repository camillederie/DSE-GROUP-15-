# Optimisation of the Main System Parameters
import numpy as np
import matplotlib.pyplot as plt
from Drum_and_tether_design import structures_calculation
from Luchsinger.luchsingermodel.Nominal_power_cycle import sensitivity_analysis
from Luchsinger.luchsingermodel.InputV2 import get_initial_data
from Aero.aero_main_function import main_aero_function
from Anchoring_design import anchoring_info
from LaunchingEquipment import launching_equipment_info


'''Aero and Power inputs'''
A_proj = 15.19
A_proj = np.linspace(10,25,20)
Points = 1000000
Kite_segments = 12
N_split = 5
AoA_range_out = np.arange(8, 15.5, 0.5)
AoA_range_in = np.arange(-1, 3.5, 0.5)
TAS = 32.5

'''Structures Input'''
len_drum = 1.10
extra_len = 75
saf_fac = 2
kite_mass_margin = 1.05

'''Inputs for Anchoring'''
Riv_w = 3692
Safety_f = 1.2
ref_kin_fric_coeff = 0.35
force_duration = 0.1
min_kite_angle = 25
max_kite_angle = 35
road_angle = 2 #degrees

'''Inputs Launching System'''
vw_ground = 10.  # m/s

# data = {}

# Open Data
def import_data(file_name):
    file = open(file_name)
    data = {}
    entries = ''
    for line in file:
        entries = entries + line.replace("\n", "")
    entries = entries[:-1]
    entries = entries.split(',')

    for entry in entries:
        key, value = entry.split(':')
        if '[' in value:
            value = value.replace('[', '')
            value = value.replace(']', '')
            value = value.split()
            data[key] = [float(x) for x in value]
        else:
            data[key] = float(value)
    return data

def sens_iteration_aero_power(run_aero, A_proj, TAS):
    Points = 1000000
    Kite_segments = 12
    N_split = 5
    AoA_range_out = np.arange(8, 15.5, 0.5)
    AoA_range_in = np.arange(-1, 3.5, 0.5)

    # print('A_proj=', A_proj)
    if run_aero == True:
        print('Aero started')
        CL_average_out, CD_average_out, CL3_CD2_average_out, CD_average_in, A_proj, Strut_area_av, flat_area, flat_area_span, chords = main_aero_function(
            A_proj, Points, Kite_segments, N_split, AoA_range_out, AoA_range_in, TAS, Print=False)
            # CL_average_in, CD_average_in, CL3_CD2_average_in, A_proj, Strut_area_av, flat_area, flat_area_span, chords = main_aero_function(A_proj, Points, Kite_segments, N_split, AoA_range_in, TAS, Print=False)
        print('Aero part finished')
        print('A_proj=', A_proj)
        ## POWER ##
        data = get_initial_data()
        data['F_out'] = CL3_CD2_average_out
        data['F_in'] = CD_average_in
        data['CL_out'] = CL_average_out
        data['CD_out'] = CD_average_out
        # data['CL_in'] = CL_average_in
        data['CD_in'] = CD_average_in
        data['Strut_area_av'] = Strut_area_av
        data['flat_area'] = flat_area
        data['flat_area_span'] = flat_area_span
        data['chords'] = chords

    data = get_initial_data()
    datasens = sensitivity_analysis(data)
    # print("The average power =", data['P_avg_elec_elev_verif'])
    # print("The average tension force =", data['T_out_elev_n'])


    # Write to file #
    datasens_save = True
    if datasens_save == True:
        file = open("data_sens.txt", "w")
        for key, value in datasens.items():
            file.write('%s:%s,\n' % (key, value))
        file.close()
        print('The extended results of the analysis can be found in the data file added to the directory.')

    return datasens

# datasens = sens_iteration_aero_power(False, A_proj, TAS)
datasens = import_data("data_sens.txt")

def sensitivity_plots(datasens):
    # datasens['F_out_list'] = data['F_out_list']
    # datasens['A_proj_list'] = data['A_proj_list']
    # datasens['v_w_list'] = data['v_w_adj']

    'FOR WINDSPEED'
    # datasens['T_out_list_VW'] = []
    # datasens['P_avg_e_list_VW'] = []
    # datasens['gamma_out_list_VW'] = []
    # datasens['gamma_in_list_VW'] = []
    # datasens['cycle_time_list_VW'] = []
    plot_wind = False
    if plot_wind == True:
        colors = datasens['v_w_list']
        fig, ax1 = plt.subplots()
        ax1.set_xlabel('Cycle Time (s)')
        ax1.set_ylabel('Nominal Reel-Out Tension Force (N)')
        plt.scatter(datasens['cycle_time_list_VW'], datasens['T_out_list_VW'],
                    c = colors, sizes = 30*datasens['v_w_list'],
                    alpha=0.5, cmap = 'viridis')
        clb = plt.colorbar()
        plt.tight_layout()
        clb.set_label('Nominal Wind Speed (m/s)', fontsize = 10)
        plt.show()

        colors = datasens['v_w_list']
        fig, ax1 = plt.subplots()
        ax1.set_xlabel('Gamma out (-)')
        ax1.set_ylabel('Average Electrical Power (W)')
        plt.scatter(datasens['gamma_out_list_VW'], datasens['P_avg_e_list_VW'],
                    c = colors, sizes = 30*datasens['v_w_list'],
                    alpha=0.5, cmap = 'viridis')
        clb = plt.colorbar()
        plt.tight_layout()
        clb.set_label('Nominal Wind Speed (m/s)', fontsize = 10)
        plt.show()

        # colors = datasens['F_out_list']
        # fig, ax1 = plt.subplots()
        # ax1.set_xlabel('Gamma out (-)')
        # ax1.set_ylabel('Gamma in (-)')
        # plt.scatter(datasens['gamma_out_list_FO'], datasens['gamma_in_list_FO'],
        #             c=colors, sizes=datasens['F_out_list'],
        #             alpha=0.5, cmap='viridis')
        # clb = plt.colorbar()  # show color scaleax1.plot(datasens['F_out_list'], datasens['gamma_out_list_FO'], color='g')
        # plt.tight_layout()
        # clb.set_label('F_out', fontsize = 10)
        # plt.show()

    'CL^3/CD^2'
    # datasens['T_out_list_FO'] = []
    # datasens['P_avg_e_list_FO'] = []
    # datasens['gamma_out_list_FO'] = []
    # datasens['gamma_in_list_FO'] = []
    # datasens['cycle_time_list_FO'] = []
    colors = datasens['F_out_list']
    fig, ax1 = plt.subplots()
    ax1.set_xlabel('Gamma out (-)')
    ax1.set_ylabel('Nominal Average Electric Power (W)')
    plt.scatter(datasens['gamma_out_list_FO'], datasens['P_avg_e_list_FO'],
                c=colors, sizes=50 * datasens['F_out_list'],
                alpha=0.5, cmap='viridis')
    clb = plt.colorbar()
    plt.tight_layout()
    clb.set_label('$CL^3/CD^2$ (m)', fontsize=12)
    plt.show()

    colors = datasens['F_out_list']
    fig, ax1 = plt.subplots()
    ax1.set_xlabel('Cycle Time (s)')
    ax1.set_ylabel('Nominal Reel-Out Tension Force (N)')
    plt.scatter(datasens['cycle_time_list_FO'], datasens['T_out_list_FO'],
                c=colors, sizes=30 * datasens['F_out_list'],
                alpha=0.5, cmap='viridis')
    clb = plt.colorbar()
    plt.tight_layout()
    clb.set_label('$CL^3/CD^2$ (m)', fontsize = 10)
    plt.show()


    fig, ax1 = plt.subplots()
    ax1.set_xlabel('CL^3/CD^2 out (-)')
    ax1.plot(datasens['F_out_list'], datasens['gamma_out_list_FO'], color='g')
    ax1.plot(datasens['F_out_list'], datasens['gamma_in_list_FO'], color='r')
    ax1.set_ylabel('y_out (-)', color='g')
    ax1.set_ylabel('y_in (-)', color='r')
    ax2 = ax1.twinx()
    ax2.plot(datasens['F_out_list'], datasens['gamma_in_list_FO'], color='b')
    ax1.legend(['Gamma Out', 'Gamma In'])
    plt.grid()
    plt.tight_layout()
    plt.show()

    fig, ax1 = plt.subplots()
    ax1.plot(datasens['F_out_list'],datasens['P_avg_e_list_FO'], color = 'g')
    ax1.set_ylabel('Average Electric Power', color = 'g')
    ax2 = ax1.twinx()
    ax2.plot(datasens['F_out_list'],datasens['T_out_list_FO'], color = 'b')
    ax1.set_xlabel('CL^3/CD^2 out (-)')
    ax2.set_ylabel('Tether Force (N)', color = 'b')
    plt.grid()
    plt.tight_layout()
    plt.show()

    'AREA'
    # datasens['T_out_list_A'] = []
    # datasens['P_avg_e_list_A'] = []

    fig, ax1 = plt.subplots()
    ax1.plot(datasens['A_proj_list'],datasens['P_avg_e_list_A'], color = 'g')
    ax1.set_ylabel('Average Electric Power', color = 'g')
    ax2 = ax1.twinx()
    ax2.plot(datasens['A_proj_list'],datasens['T_out_list_A'], color = 'b')
    ax1.set_xlabel('Projected Area (m^2)')
    ax2.set_ylabel('Tether Force (N)', color = 'b')
    plt.grid()
    plt.tight_layout()
    plt.show()

# sensitivity_plots(datasens)


data = import_data("data_optim.txt")

# T_F_out = data['T_out_elev_n']

'Inputs for structures, outputs from aero_power'
kite_mass_list = []
tether_d_list = []
drum_d_list = []
req_fric_coef_list = []
F_clamp_list = []

i = 0
for area in datasens['A_proj_list']:
    nom_load = datasens['T_out_list_A'][i]

    if area == data['A_proj']:
        factor = 1
    else:
        factor = area/data['A_proj']
    flat_area = factor*data['flat_area']
    Strut_area_av = factor*data['Strut_area_av']
    'there is no cycle time for changing area'
    # t_out = datasens['cycle_time'][i] * (1 / datasens['gamma_out_n'][i]) / ((1 / datasens['gamma_out_n'][i]) + 1 / datasens['gamma_in_n'][i])
    # t_in = datasens['cycle_time'][i] * (1 / datasens['gamma_in_n'][i]) / ((1 / datasens['gamma_in_n'][i]) + 1 / datasens['gamma_out_n'][i])

    t_out = data['t_out']
    t_in = data['t_in']

    #functions
    kite_mass_ALUULA, tether_diameter, tether_mass, load, drum_D = structures_calculation(area, Strut_area_av, len_drum,
                                                                                          data['a_elev_out'], extra_len,
                                                                                          nom_load, saf_fac,
                                                                                          kite_mass_margin, t_out,
                                                                                          t_in)

    req_fric_coeff, sliding_d_mm = anchoring_info(nom_load * Safety_f, Riv_w, road_angle, min_kite_angle,
                                                  max_kite_angle, ref_kin_fric_coeff, force_duration)

    chord = factor * max(data['chords'])
    F_clamp = launching_equipment_info(vw_ground, req_fric_coeff, flat_area, chord, kite_mass_ALUULA)

    kite_mass_list.append(kite_mass_ALUULA)
    tether_d_list.append(tether_diameter)
    drum_d_list.append(drum_D)
    req_fric_coef_list.append(req_fric_coeff)
    F_clamp_list.append(F_clamp)
    #
    i += 1

fig, ax1 = plt.subplots()
ax1.plot(datasens['T_out_list_A'], req_fric_coef_list, color = 'g')
ax1.set_xlabel('Tension Force (N)')
ax1.set_ylabel('Require Friction Coefficient (-)', color='g')
plt.grid()
plt.tight_layout()

fig, ax1 = plt.subplots()
ax1.plot(datasens['A_proj_list'],kite_mass_list, color = 'g')
ax1.set_ylabel('Kite Mass (kg)', color = 'g')
ax2 = ax1.twinx()
ax2.plot(datasens['A_proj_list'],tether_d_list, color = 'b')
ax1.set_xlabel('Projected Area (m^2)')
ax2.set_ylabel('Tether Diameter (m)', color = 'b')
plt.grid()
plt.tight_layout()
plt.show()

fig, ax1 = plt.subplots()
ax1.plot(datasens['A_proj_list'],datasens['T_out_list_A'], color = 'g')
ax1.set_xlabel('Kite Area ($m^2$)')
ax2 = ax1.twinx()
ax2.plot(datasens['A_proj_list'],np.array(tether_d_list)*1000, color = 'b')
ax1.set_ylabel('Tension force (N)', color='g')
ax2.set_ylabel('Tether Diameter (mm)', color = 'b')
plt.grid()
plt.tight_layout()
plt.show()












'Inputs for structures, outputs from aero_power'
nom_load = data['T_out_elev_n']
t_out = data['t_out']
t_in = data['t_in']
Strut_area_av = data['Strut_area_av']
A_proj = data['A_proj']

''''functions for optimised design (one data array)'''
# ite_mass_ALUULA, Tether_diameter, Tether_mass, Load, D_drum = structures_calculation(A_proj, Strut_area_av, len_drum,
#                                                                                       data['a_elev_out'], extra_len,
#                                                                                       nom_load, saf_fac,
#                                                                                       kite_mass_margin, data['t_out'],
#                                                                                       data['t_in'])
#
# req_fric_coeff, sliding_d_mm = anchoring_info(nom_load * Safety_f, Riv_w, road_angle, min_kite_angle, max_kite_angle,
#                                               ref_kin_fric_coeff, force_duration)
#
# launching_equipment_info(vw_ground, req_fric_coeff, data['flat_area'], max(data['chords']), Kite_mass_ALUULA)