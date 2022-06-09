from math import *
# This code calculates the downforce the clamp needs to provide to hold the kite still during launching operations
vw_ground = 8  # m/s, ground speed
angle_kite_ground = pi / 4
Cd = 0.9  # -, drag coefficient of an angled plate
rho = 1.225  # kg/m^3, air density
friction_coefficient = 0.5  # minimum ground coefficient coming from anchoring
S = 20.  # m^2, surface area of the kite
c = 2.  # m, chord length in middle = largest
kite_mass = 6.  # kg, mass of kite only
g = 9.81  # gravitational constant

def launching_equipment_info(vw_ground, friction_coefficient, S, c, kite_mass):
    l_1 = c * cos(angle_kite_ground)  # m, projected horizontal distance/width of the kite
    l_2 = c * sin(angle_kite_ground)  # m, projected vertical distance/height of the kite
    #Since an angle of 45 degrees is chosen, l1 = l2, and they can be dropped for the moment equations.

    F_wind = 0.5*rho*sin(angle_kite_ground)*S*Cd*vw_ground**2 #Wind Force on kite

    """"
    For friction: F_friction = friction_coefficient * normal force
    Where normal force = F_clamp + weight_kite * g.
    
    Moment Equation gives:
    F_tension = Fclamp + 0.5 *Fwind
    """

    F_clamp = (F_wind/2 + friction_coefficient * kite_mass * g)/ (1+friction_coefficient)

    print(f'The downforce that needs to be provided by the launching equipment is {F_clamp} N \n'
          f'This means a weight of at least {F_clamp/g} kg of the launching equipment.')
    return F_clamp

if __name__ == "__main__":
    launching_equipment_info(vw_ground, friction_coefficient, S, c, kite_mass)