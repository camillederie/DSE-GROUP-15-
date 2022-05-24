

v_w_n = 10 #np.linspace(5,20,50)
CL_out = 0.6 #np.linspace(0.6,1.5,50)
A_proj =  20.79 #21.57 #np.linspace(15,35,50) #19.8 #

rho = 1.18
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

P_w = 0.5*v_w_n**3*rho # Wind Power 