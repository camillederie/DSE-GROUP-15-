from cProfile import label
import numpy as np
import matplotlib.pyplot as plt

## Parameters
v_w = 10 #np.linspace(5,20,50)
CL_out = 1.0 #np.linspace(0.6,1.5,50)
A_proj = 21.87 #np.linspace(15,35,50) #19.8 #

rho = 1.18
g_out = 0.35
g_in = 1.6
lc = 250
CD_out = 0.2
CL_in = 0.14
CD_in = 0.07
eff_in = 0.579
eff_out = 0.591

P_avg_e_n = 20000 #Nominal electrical power (W)

## Intermediate calculations

F_out = CL_out**3/CD_out**2
F_in = CD_in

tc = (lc/v_w)*((g_out+g_in)/(g_out*g_in))
t_in = 1/6*tc
t_out = 5/6*tc
print(tc)

E_control = 150*tc #Energy usage from control system (W)

T_out = 0.5*rho*v_w**2*A_proj*(1-g_out)**2*F_out
T_in = 0.5*rho*v_w**2*A_proj*(1+g_in)**2*F_in

P_out = T_out*g_out*v_w
P_out_e = P_out * eff_out

P_in = T_in*g_in*v_w
P_in_e = P_in / eff_in

P_avg_mech = P_out*(g_in)/(g_in + g_out) - P_in*g_out/(g_in + g_out)
P_avg_elec = P_out_e*(g_in)/(g_in + g_out) - P_in_e*g_out/(g_in + g_out) 
P_M_G = P_out*0.77*0.85*0.95

## Energy Calculations ##

Ec_e = (T_out*eff_out - T_in*eff_in)*lc
Ec = (T_out - T_in)*lc 
E_in_e = P_in_e * t_in #in Joules
E_out_e = P_out_e *t_out #in Joules
E_supcap_pc = (E_control + E_in_e)*0.000277777778 #in Wh
R_E = E_in_e/ E_out_e

## Power Calculations ##

P_w = 0.5*v_w**3*rho # Wind Power 
P_avg_m = P_w*A_proj*(F_out*(1-g_out)**2-(F_in*(1+g_in)**2))*((g_out*g_in)/(g_out+g_in)) #Ec/tc # Average Mechanical Power 
P_avg_e = P_w*A_proj*(eff_out*F_out*(1-g_out)**2-(F_in*(1+g_in)**2)/eff_in)*((g_out*g_in)/(g_out+g_in)) # Average Electrical Power

## OUTPUT ##

print('The traction force (N) reel out and in are: ',T_out,T_in)
print('The mechanical power (W) reel out and in are: ',P_out,P_in)
print('The electrical power (W) reel out and in are: ',P_out_e, P_in_e)
print('The total electrical energy (J) reel out and in and the in/out ratio is: ', E_out_e, E_in_e, R_E)
print('The energy stored in supercap per cycle (Wh): ', E_supcap_pc)
print('The average mechanical power (W) and electrical power is: ',P_avg_m, P_avg_e, P_avg_elec)
print('Motor generator power (W): ', P_M_G)

## Plots ##

def plotcl():
    
    CL_out = np.linspace(0.6,1.4,50)
    F_out = CL_out**3/CD_out**2

    # POWER
    P_w = 0.5*v_w**3*rho # Wind Power 
    P_avg_m = P_w*A_proj*(F_out*(1-g_out)**2-(F_in*(1+g_in)**2))*((g_out*g_in)/(g_out+g_in)) #Ec/tc # Average Mechanical Power 
    P_avg_e = P_w*A_proj*(eff_out*F_out*(1-g_out)**2-(F_in*(1+g_in)**2)/eff_in)*((g_out*g_in)/(g_out+g_in)) # Average Electrical Power
    
    plt.plot(CL_out, P_avg_m, label = 'Potential Mechanical Power')
    plt.plot(CL_out,P_avg_e,'--', label = 'Potential Electrical Output Power')
    plt.xlabel('Reel-out $C_{L}$',fontsize = 16)
    plt.ylabel('Power (W)',fontsize = 16)
    plt.legend(fontsize = 16)
    plt.grid()
    plt.show()

def plotvw():

    v_w = np.linspace(5,16,50)

    # POWER
    P_w = 0.5*v_w**3*rho # Wind Power 
    P_avg_m = P_w*A_proj*(F_out*(1-g_out)**2-(F_in*(1+g_in)**2))*((g_out*g_in)/(g_out+g_in)) #Ec/tc # Average Mechanical Power 
    P_avg_e = P_w*A_proj*(eff_out*F_out*(1-g_out)**2-(F_in*(1+g_in)**2)/eff_in)*((g_out*g_in)/(g_out+g_in)) # Average Electrical Power

    plt.plot(v_w, P_avg_m, label = 'Potential Mechanical Power')
    plt.plot(v_w,P_avg_e,'--', label = 'Potential Electrical Output Power')
    plt.xlabel('Wind Speed (m/s)',fontsize = 16)
    plt.ylabel('Power (W)',fontsize = 16)
    plt.legend(fontsize = 16)
    plt.grid()
    plt.show()

def plotA():
    
    A_proj = np.linspace(10,28,50) #19.8 #

    # POWER
    P_w = 0.5*v_w**3*rho # Wind Power 
    P_avg_m = P_w*A_proj*(F_out*(1-g_out)**2-(F_in*(1+g_in)**2))*((g_out*g_in)/(g_out+g_in)) #Ec/tc # Average Mechanical Power 
    P_avg_e = P_w*A_proj*(eff_out*F_out*(1-g_out)**2-(F_in*(1+g_in)**2)/eff_in)*((g_out*g_in)/(g_out+g_in)) # Average Electrical Power

    plt.plot(A_proj, P_avg_m, label = 'Potential Mechanical Power')
    plt.plot(A_proj,P_avg_e,'--', label = 'Potential Electrical Output Power')
    plt.xlabel('Projected Area of the Kite ($m^2$)',fontsize = 16)
    plt.ylabel('Power (W)',fontsize = 16)
    plt.legend(fontsize = 16)
    plt.grid()
    plt.show()

def plotA_wind():

    v_w = np.linspace(6.5,18,50)
    P_w = 0.5*v_w**3*rho
    A_proj = P_avg_e_n/P_w/((eff_out*F_out*(1-g_out)**2-(F_in*(1+g_in)**2)/eff_in)*((g_out*g_in)/(g_out+g_in)))
    plt.plot(v_w,A_proj)
    plt.ylabel('Projected Area of the Kite ($m^2$)', fontsize = 16)
    plt.xlabel('Wind Speed (m/s)',fontsize = 16)
    plt.grid()
    plt.show()
 
def plotA_CL():

    CL_out = np.linspace(0.73,1.5,50)
    F_out = CL_out**3/CD_out**2
    P_w = 0.5*v_w**3*rho
    A_proj = P_avg_e_n/P_w/((eff_out*F_out*(1-g_out)**2-(F_in*(1+g_in)**2)/eff_in)*((g_out*g_in)/(g_out+g_in)))
    plt.plot(CL_out,A_proj)
    
    plt.ylabel('Projected Area of the Kite ($m^2$)',fontsize = 16)
    plt.xlabel('Reel-out $C_{L}$',fontsize = 16)
    plt.grid()
    plt.show()



# plotcl()
# plotvw()
# plotA()
# plotA_wind()
# plotA_CL()
