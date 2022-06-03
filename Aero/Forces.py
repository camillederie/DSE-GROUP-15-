from numpy import *
F_total = 10329
F_total_power = 0.7*F_total #[N]
alpha_0 = radians(13)
beta_0 = radians(13)
alpha_1 = radians(33)
beta_1 = radians(21)
alpha_2 = radians(26)
beta_2 = radians(20)
alpha_3 = radians(34)
beta_3 = radians(25)

print(F_total_power)
#cross-section_0
F_symmetry = (cos(alpha_0)*sin(beta_0)/sin(alpha_0)+cos(beta_0))**(-1)*F_total_power
F_symmetry_1 = sin(beta_0)/sin(alpha_0)*F_symmetry
print(F_symmetry)
print(F_symmetry_1)

#cross-section_1
F_2 = (cos(alpha_1)*sin(beta_1)/sin(alpha_1)+cos(beta_1))**(-1)*F_symmetry
F_1 = sin(beta_1)/sin(alpha_1)*F_2
print(F_2)
print(F_1)

#cross_section_2
F_3 = (cos(alpha_2)*sin(beta_2)/sin(alpha_2)+cos(beta_2))**(-1)*F_2
F_4 = sin(beta_2)/sin(alpha_2)*F_3
print(F_3)
print(F_4)

#cross_section_3
F_5 = (cos(alpha_3)*sin(beta_3)/sin(alpha_3)+cos(beta_3))**(-1)*F_1
F_6 = sin(beta_3)/sin(alpha_3)*F_5
print(F_5)
print(F_6)


list_bridle_diameters = [0.003, 0.003, 0.003, 0.002, 0.002, 0.002, 0.002]
list_bridle_lengths = [17.841458317180013, 3.2319258311826022, 2.8447202154208635, 3.9423003053088035, 3.0836634428987466, 3.0658060311942874, 2.54608450485428]
print(sum(list_bridle_lengths))
list_tether_mass = []
for i in range(len(list_bridle_diameters)):
    if list_bridle_diameters[i] == 0.003:
        tether_mass = list_bridle_lengths[i]/100*0.5
        list_tether_mass.append(tether_mass)
    elif list_bridle_diameters[i] == 0.002:
        tether_mass = list_bridle_lengths[i]/100*0.2
        list_tether_mass.append(tether_mass)

#print(list_tether_mass)
print(sum(list_tether_mass))