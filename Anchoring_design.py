from math import cos,sin,radians
import numpy as np

#Functions 
def req_fric_coeff(alfa,beta,Riv_w,force): #ranges of beta
    alfa_list = np.arange(-alfa,alfa+0.5,0.5)
    beta_list = np.arange(beta[0],beta[1]+0.5,0.5)
    req_fric_coeff = 0
    alfa_extr = 0
    beta_extr = 0
    for i in alfa_list:
        for x in beta_list:
            fric_coeff = (force*cos(radians(i+x))+Riv_w*9.81*sin(radians(i)))/(Riv_w*9.81*cos(radians(i))-force*sin(radians(i+x)))
            if fric_coeff > req_fric_coeff:
                req_fric_coeff = fric_coeff
                alfa_extr = i
                beta_extr = x
    return req_fric_coeff, alfa_extr, beta_extr

def sliding_dist(ref_fric_coeff,force_duration,Riv_w,force,alfa,beta):
    fric_force = (force*cos(radians(alfa+beta))+Riv_w*9.81*sin(radians(alfa)))
    Norm_force = (Riv_w*9.81*cos(radians(alfa))-force*sin(radians(alfa+beta)))
    f_res = fric_force - Norm_force*ref_fric_coeff
    d = 0.5*f_res/Riv_w*force_duration**2
    return d

def vert_axis_tilting(force,kite_angle,road_angle,weight,wheels_dist_len,loc_cg,dist_cg_ground): #kite angle in range
    alfa_list = np.arange(-road_angle,road_angle+0.5,0.5)
    beta_list = np.arange(kite_angle[0],kite_angle[1]+0.5,0.5)
    rear_tire_force_lst = []
    for i in alfa_list:
        for x in beta_list:
            ver_t_force = force*sin(radians(x-i))
            rear_tire_force = (cos(radians(i))*weight*9.81*(wheels_dist_len-loc_cg)-ver_t_force*wheels_dist_len-sin(radians(i))*weight*9.81*(dist_cg_ground))/wheels_dist_len
            rear_tire_force_lst.append(rear_tire_force)
    if min(rear_tire_force_lst)>= 0:
        return True, min(rear_tire_force_lst)
    else:
        return False, min(rear_tire_force_lst)

def long_axis_tilting(force,kite_angle,road_angle,weight,wheels_dist_width,dist_drum_cg,dist_cg_ground): #kite angle in range
    alfa_list = np.arange(-road_angle,road_angle+0.5,0.5)
    beta_list = np.arange(kite_angle[0],kite_angle[1]+0.5,0.5)
    left_tire_force_lst = []
    for i in alfa_list:
        for x in beta_list:
            left_tire_force = weight*9.81*0.5*cos(radians(i))-0.5*force*sin(radians(x+i))-(dist_drum_cg+dist_cg_ground)*force*cos(radians(x+i))/wheels_dist_width-weight*9.81*sin(radians(i))*dist_cg_ground
            left_tire_force_lst.append(left_tire_force)
    if min(left_tire_force_lst)>= 0:
        return True, min(left_tire_force_lst)
    else:
        return False, min(left_tire_force_lst)

if __name__ == "__main__":

    #Inputs 
    Riv_w = 3692 
    Kite_f_without_sf = 10329.32
    Safety_f = 1.2
    Kite_f = Kite_f_without_sf * Safety_f

    ref_kin_fric_coeff = 0.35
    force_duration = 0.1

    min_kite_angle = 25
    max_kite_angle = 35
    road_angle = 2

    loc_cg = 1.536
    wheels_dist_len = 3.43
    height_riv = 1.83 
    dist_drum_cg = 0.7
    dist_cg_ground = 0.4
    wheels_dist_width = 1.9 

    #Sliding
    req_fr_coeff, alfa_extr, beta_extr = req_fric_coeff(road_angle,[min_kite_angle,max_kite_angle],Riv_w,Kite_f)
    sliding_d = sliding_dist(ref_kin_fric_coeff,force_duration,Riv_w,Kite_f,alfa_extr,beta_extr)

    #Tilting vertical axis
    stationary_bol_vert, min_tireforce_vert = vert_axis_tilting(Kite_f,[min_kite_angle,max_kite_angle],road_angle,Riv_w,wheels_dist_len,loc_cg,dist_cg_ground)

    #Tilting longitudinal axis
    stationary_bol_long, min_tireforce_long = long_axis_tilting(Kite_f,[min_kite_angle,max_kite_angle],road_angle,Riv_w,wheels_dist_width,dist_drum_cg,dist_cg_ground)

    print(f"Required Friction Coefficient = {req_fr_coeff}, currently it slides {sliding_d*1000} mm")
    print(f"Tilting over vertical axis: Stationary: {stationary_bol_vert}, with minimal normal force = {min_tireforce_vert} N")
    print(f"Tilting over longitudinal axis: Stationary: {stationary_bol_long}, with minimal normal force = {min_tireforce_long} N")

