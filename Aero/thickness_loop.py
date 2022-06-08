# -*- coding: utf-8 -*-
"""
Created on Mon Jan 10 18:09:16 2022

@author: oriol2
"""

import os
import numpy as np
import matplotlib.pyplot as plt
import time
import functions_VSM_LLT as VSM
import Coordinates_generation as CGEN

script_dir = os.path.dirname(__file__)
results_dir = os.path.join(script_dir, 'Results/')

''' Code iteration part for polar graphs generation '''

#%%  Input DATA

Polars = True  # Set parameter for loop running
Kite15 = True   # Set parameter for own kite or v3 kite

# Definition of the geometry
if Kite15 == True:
    Name1 = 'Kite15'
    Atot = 21.54  # Projected area
    Segments = 10  # Number of kite segments
    Points = 1000000  # Number of points for coord function
    Plotting = False
    # CAD = VSM.get_Kite15_coords()  # Geometry nodes locations
    coords, chords, MAC, arclength = CGEN.Generate_Kite15_coords(Atot, Segments, Points, Plotting)
    # coords = VSM.struct2aero_geometry(CAD)  # Change geometry to definition
    N = int(len(coords) / 2)  # Number of sections defined

    # LE thickness at each section [m]
    t = [0.104, 0.16, 0.18, 0.19, 0.195, 0.195, 0.19, 0.18, 0.16, 0.104]
    t = [0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2]
    t = 0.1*chords

    # Camber for each section (ct in my case)
    # If you want it to vary do the same as the thickness
    k = 0.08

else:
    Name1 = 'v3'
    Atot = 19.753  # Projected area
    CAD = VSM.get_CAD_matching_uri()  # Geometry nodes locations
    coords = VSM.struct2aero_geometry(CAD) / 1000  # Change geometry to definition
    N = int(len(coords) / 2)  # Number of sections defined

    # LE thickness at each section [m]
    t = [0.118753, 0.151561, 0.178254, 0.19406, 0.202418, 0.202418, 0.19406, 0.178254, 0.151561, 0.118753]

    # Camber for each section (ct in my case)
    # If you want it to vary do the same as the thickness
    k = 0.095

#   Model and program specifics
ring_geo = '5fil' # System of vorticity defined with 5 or 3 filaments per wing section
model = 'VSM'     # Choose between Vortex Step Method (VSM) or Lifting Line Method (LLT)
# Convergence criteria
conv_crit = {
    'Niterations': 1500,
    'error' : 1e-5,
    'Conv_weight': 0.03
    }

thicknesses = np.arange(0.15, 0.25, 0.05)
thicknesses = [0.08*chords, 0.1*chords, 0.12*chords, 0.14*chords,
               0.04*MAC*np.ones(len(chords)), 0.06*MAC*np.ones(len(chords)), 0.08*MAC*np.ones(len(chords)), 0.1*MAC*np.ones(len(chords))]
CD_loop = []
CL_loop = []
LD_loop = []
Power_loop = []
aoa_0 = []

if Polars == True:
    for i in range(len(thicknesses)):
        # t = [thicknesses[i], thicknesses[i], thicknesses[i], thicknesses[i], thicknesses[i], thicknesses[i], thicknesses[i], thicknesses[i],thicknesses[i], thicknesses[i]]
        t = thicknesses[i]
        fig = plt.figure(figsize=(6, 5))

        aoa_0 = np.arange(6, 19.5, 1.5)
        CD_x = []
        CL_x = []

        for angle in aoa_0:

            # Wind speed vector definition
            Umag = 22.5               # Wind speed magnitude
            aoa = angle*np.pi/180      # Angle of attack
            sideslip = 0/180*np.pi  # Sideslip angle
            Uinf = np.array([np.cos(aoa)*np.cos(sideslip),np.sin(sideslip),np.sin(aoa)])*Umag # Wind speed vector
            N = int(len(coords) / 2)  # Number of sections defined
            # Plot thickness distribution
            indexes = np.empty(N)
            for i in range(N):
                indexes[i]=coords[2*i,1]
            # fig = plt.figure(figsize = (6,5))
            # plt.plot(indexes,t,'ro',label='CAD Drawing')


            # Number of splits per panel
            N_split = 4
            # Refine structrural mesh into more panels
            coord = VSM.refine_LEI_mesh(coords, N-1, N_split)
            # Define system of vorticity
            controlpoints,rings,bladepanels,ringvec,coord_L = VSM.create_geometry_LEI(
                coord, Uinf, int((N-1)*N_split+1),ring_geo,model)
            N = int(len(coord)/2)           # Number of section after refining the mesh


            #%% Airfoil Coefficients

            # Definition of the thickness distribution for the refined mesh
            thicc = np.array([])
            for i in range(Segments-1):
                temp = np.linspace(t[i],t[i+1],N_split+1)
                temp1 = []
                for a in range(len(temp)-1):
                    temp1.append((temp[a] +temp[a+1])/2)
                thicc = np.append(thicc,temp1)

            # Definition of airfoil coefficients
            # Based on Breukels (2011) correlation model
            aoas = np.arange(-20,21,1)
            data_airf = np.empty((len(aoas),4,N-1))
            t_c = np.empty(N-1)
            for i in range(N-1):
                for j in range(len(aoas)):
                    t_c[i] = thicc[i]/controlpoints[i]['chord']
                    alpha = aoas[j]
                    Cl,Cd,Cm = VSM.LEI_airf_coeff(t_c[i], k, alpha)
                    data_airf[j,0,i] = alpha
                    data_airf[j,1,i] = Cl
                    data_airf[j,2,i] = Cd
                    data_airf[j,3,i] = Cm

            #%% SOLVER
            # Define system of vorticity
            controlpoints,rings,bladepanels,ringvec,coord_L = VSM.create_geometry_LEI(
                coord, Uinf, N,ring_geo,model)
            # Solve for Gamma
            Fmag, Gamma,aero_coeffs = VSM.solve_lifting_line_system_matrix_approach_semiinfinite(
                ringvec, controlpoints, rings,Uinf,data_airf,conv_crit,model)
            #%OUTPUT Results
            F_rel,F_gl,Ltot,Dtot,CL,CD,Lift = VSM.output_results(Fmag,aero_coeffs,ringvec,Uinf,controlpoints,Atot)

            print(angle, CL,CD)

            CD_x.append(CD)
            CL_x.append(CL)


        CD_x = np.array(CD_x)
        CL_x = np.array(CL_x)
        LD = CL_x/CD_x
        Power = CL_x**3 / CD_x**2

        CD_loop.append(CD_x)
        CL_loop.append(CL_x)
        LD_loop.append(LD)
        Power_loop.append(Power)

    # Plotting
    linestyles = ["-",":","--","-.","-",":","--","-."]
    fig, ax1 = plt.subplots()
    fig, ax2 = plt.subplots()
    fig, ax3 = plt.subplots()
    fig, ax4 = plt.subplots()

    labels = ["8% c", '10% c', '12% c', '14% c', '4% MAC', '6% MAC', '8% MAC', '10% MAC']

    for i in range(len(CL_loop)):
        ax1.plot(aoa_0, LD_loop[i], linestyle=linestyles[i], label=f'{labels[i]}')
        ax2.plot(aoa_0, Power_loop[i], linestyle=linestyles[i], label =f'{labels[i]}')
        ax3.plot(aoa_0, CL_loop[i], linestyle=linestyles[i], label =f'{labels[i]}')
        ax4.plot(aoa_0, CD_loop[i], linestyle=linestyles[i], label =f'{labels[i]}')

    ax1.set_xlabel(r'$\alpha$ [deg]')
    ax1.set_ylabel(r'$C_L/C_D$ [-]')
    ax1.legend()
    fig.tight_layout()

    ax2.set_xlabel(r'$\alpha$ [deg]')
    ax2.set_ylabel(r'$C_L^3/C_D^2$ [-]')
    ax2.legend()
    fig.tight_layout()
    # plt.savefig(results_dir + Name1 + '_Power_coeff_a.png')

    ax3.set_xlabel(r'$\alpha$ [deg]')
    ax3.set_ylabel(r'$C_{L}$ [-]')
    ax3.legend()
    plt.tight_layout()
    # plt.savefig(results_dir + Name1 + '_CL_a.png')
    # plt.legend()

    ax4.set_xlabel(r'$\alpha$ [deg]')
    ax4.set_ylabel(r'$C_{D}$ [-]')
    ax4.legend()
    plt.tight_layout()
    # plt.savefig(results_dir + Name1 + '_CD_a.png')
    # plt.legend()
    for el in thicknesses:
        print(el[0])
    plt.show()

else:
    print('No implementation')
