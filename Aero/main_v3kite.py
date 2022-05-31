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
    Atot = 13.74  # Projected area
    Segments = 10  # Number of kite segments
    Points = 1000000  # Number of points for coord function
    Plotting = False
    # CAD = VSM.get_Kite15_coords()  # Geometry nodes locations
    coords, chords = CGEN.Generate_Kite15_coords(Atot, Segments, Points, Plotting)
    # coords = VSM.struct2aero_geometry(CAD)  # Change geometry to definition
    N = int(len(coords) / 2)  # Number of sections defined

    # LE thickness at each section [m]
    t = [0.104, 0.16, 0.18, 0.19, 0.195, 0.195, 0.19, 0.18, 0.16, 0.104]
    t = [0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2]
    # t = 0.1*chords

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

if Polars == True:

    aoa_0 = np.arange(0, 17.5, 0.5)
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
        for i in range(9):
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
        F_rel,F_gl,Ltot,Dtot,CL,CD = VSM.output_results(Fmag,aero_coeffs,ringvec,Uinf,controlpoints,Atot)

        print(angle, CL,CD)

        CD_x.append(CD)
        CL_x.append(CL)


    print(CD_x)
    print(CL_x)
    CD_x = np.array(CD_x)
    CL_x = np.array(CL_x)
    LD = CL_x/CD_x
    Power = CL_x**3 / CD_x**2
    print(np.max(LD))
    print(np.max(Power))

    fig, ax1 = plt.subplots()
    ax1.plot(aoa_0, LD, color='red', label=r'$C_{L}/C_{D}$')
    ax2 = ax1.twinx()
    ax2.plot(aoa_0, Power, linestyle='--', color='blue', label=r'$C_{L}^3/C_{D}^2$')
    ax1.set_xlabel(r'$\alpha$ [deg]')
    ax1.set_ylabel(r'$C_L/C_D$ [-]', color='red')
    ax2.set_ylabel(r'$C_L^3/C_D^2$ [-]', color='blue')
    # fig.legend()
    ax1.legend()
    ax2.legend(loc=4)
    fig.tight_layout()
    #plt.savefig(results_dir + Name1 + '_Power_coeff_a.png')

    fig = plt.figure(figsize=(6, 5))
    plt.plot(aoa_0, CL_x)
    plt.xlabel(r'$\alpha$ [deg]')
    plt.ylabel(r'$C_{L}$ [-]')
    plt.tight_layout()
    # plt.savefig(results_dir + Name1 + '_CL_a.png')
    # plt.legend()

    fig = plt.figure(figsize=(6, 5))
    plt.plot(aoa_0, CD_x)
    plt.xlabel(r'$\alpha$ [deg]')
    plt.ylabel(r'$C_{D}$ [-]')
    plt.tight_layout()
    # plt.savefig(results_dir + Name1 + '_CD_a.png')
    # plt.legend()
    print('CD = ', np.average(CD_x))
    print('CL = ', np.average(CL_x))

    plt.show()

else:
    ''' Code for individual case analysis + geometry implementation '''
    # Wind speed vector definition
    Umag = 22.5               # Wind speed magnitude
    aoa = 15*np.pi/180      # Angle of attack
    sideslip = 0/180*np.pi  # Sideslip angle
    Uinf = np.array([np.cos(aoa)*np.cos(sideslip),np.sin(sideslip),np.sin(aoa)])*Umag # Wind speed vector

    # Plot thickness distribution
    indexes = np.empty(N)
    for i in range(N):
        indexes[i]=coords[2*i,1]
    # fig = plt.figure(figsize = (6,5))
    # plt.plot(indexes,t,'ro',label='CAD Drawing')


    # Number of splits per panel
    N_split = 8
    # Refine structural mesh into more panels
    coord = VSM.refine_LEI_mesh(coords, N-1, N_split)
    # Define system of vorticity
    controlpoints,rings,bladepanels,ringvec,coord_L = VSM.create_geometry_LEI(
        coord, Uinf, int((N-1)*N_split+1),ring_geo,model)
    N = int(len(coord)/2)           # Number of section after refining the mesh


    #%% Airfoil Coefficients

    # Definition of the thickness distribution for the refined mesh
    thicc = np.array([])
    for i in range(9):
        temp = np.linspace(t[i],t[i+1],N_split+1)
        temp1 = []
        for a in range(len(temp)-1):
            temp1.append((temp[a] +temp[a+1])/2)
        thicc = np.append(thicc,temp1)
    #Plot thickness distribution
    plt.plot(coord_L[:,1],thicc,color='blue',label='Fit')
    plt.xlabel(r'$y$ [m]')
    plt.ylabel(r'$t$ [m]')
    plt.legend()

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

    # Plot thickness ratio along the span
    fig = plt.figure(figsize = (6,5))
    plt.plot(coord_L[:,1],t_c)
    plt.ylabel(r'Thickness ratio ($t/c$)')
    plt.xlabel(r'$y$ [m]')

    #%% SOLVER
    # Define system of vorticity
    controlpoints,rings,bladepanels,ringvec,coord_L = VSM.create_geometry_LEI(
        coord, Uinf, N,ring_geo,model)
    # Solve for Gamma
    Fmag, Gamma,aero_coeffs = VSM.solve_lifting_line_system_matrix_approach_semiinfinite(
        ringvec, controlpoints, rings,Uinf,data_airf,conv_crit,model)
    #%OUTPUT Results
    F_rel,F_gl,Ltot,Dtot,CL,CD = VSM.output_results(Fmag,aero_coeffs,ringvec,Uinf,controlpoints,Atot)

    print(f'alpha = {aoa*180/np.pi}, CL = {CL}, CD = {CD}')

    # #%% PLOTS

    fig = plt.figure(figsize = (6,5))
    plt.plot(coord_L[:,1]/max(coord_L[:,1]),aero_coeffs[:,1], label = model)
    plt.xlabel(r'$y/s$')
    plt.ylabel(r'$C_l$')
    plt.legend()


    fig = plt.figure(figsize = (6,5))
    plt.plot(coord_L[:,1]/max(coord_L[:,1]),aero_coeffs[:,0]*180/np.pi, label = model)
    plt.xlabel(r'$y/s$')
    plt.ylabel(r'$\alpha$')
    plt.legend()


    fig = plt.figure(figsize = (6,5))
    plt.plot(coord_L[:,1]/max(coord_L[:,1]),aero_coeffs[:,2], label = model)
    plt.xlabel(r'$y/s$')
    plt.ylabel(r'Viscous $C_d$')
    plt.legend()


    fig = plt.figure(figsize = (6,5))
    plt.plot(coord_L[:,1]/max(coord_L[:,1]),Gamma, label = model)
    plt.xlabel(r'$y/s$')
    plt.ylabel(r'$\Gamma$')
    plt.legend()



    #%% PLOT GEOMETRY

    # xy plane
    fig, ax = plt.subplots(figsize = (10,6))
    for panel in bladepanels:
            coord = np.array([panel['p1'],panel['p2'],panel['p3'],panel['p4']])
            plt.plot(coord[:,1],coord[:,0],  'k',linestyle = '--')
            plt.plot(panel['p1'][1],panel['p1'][0],'o',mfc = 'white',c = 'k')
            plt.plot(panel['p4'][1],panel['p4'][0],'o',mfc = 'white',c = 'k')
    for cp in controlpoints:
        plt.plot(cp['coordinates'][1],cp['coordinates'][0], c = 'tab:green', marker = '*',markersize = 12,label = 'cp VSM')
        plt.plot(cp['coordinates_aoa'][1],cp['coordinates_aoa'][0],c = 'tab:red', marker = 'x', markersize = 12,label = 'cp LLT')

    # ax.set_ylim(0.8, -0.3)  # decreasing time
    plt.grid(linestyle = ':', color = 'gray')

    #zy plane
    fig, ax = plt.subplots(figsize = (10,6))
    for panel in bladepanels:
            coord = np.array([panel['p1'],panel['p2'],panel['p3'],panel['p4']])
            plt.plot(coord[:,1],coord[:,2],  'k',linestyle = '--')
            # plt.plot(panel['p1'][1],panel['p1'][2],'o',mfc = 'white',c = 'k')
            # plt.plot(panel['p4'][1],panel['p4'][2],'o',mfc = 'white',c = 'k')
    for cp in controlpoints:
        plt.plot(cp['coordinates'][1],cp['coordinates'][2], c = 'tab:green', marker = '*',markersize = 12,label = 'cp VSM')
        plt.plot(cp['coordinates_aoa'][1],cp['coordinates_aoa'][2],c = 'tab:red', marker = 'x', markersize = 12,label = 'cp LLT')

    # ax.set_ylim(0.8, -0.3)  # decreasing time
    plt.grid(linestyle = ':', color = 'gray')
    VSM.plot_geometry(bladepanels,controlpoints,rings,F_gl,coord_L,'True')
    plt.show()
