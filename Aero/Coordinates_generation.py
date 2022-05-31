import numpy as np
import matplotlib.pyplot as plt

def Generate_Kite15_coords(Projected_area,segments, Points_no, Plotting):
    ''' a is an ellipse's semi-major axis
        b is an ellipse's semi-minor axis '''

    Top_ellipse_ratio = 3
    Front_ellipse_ratio = 2
    AR_Projected = 1/((1/(2*Top_ellipse_ratio) * np.pi + 4/(3*Top_ellipse_ratio))/4)    # Based on 0.4 Taper ratio
    print("# ----- Kite Data Printing ----- #\nProjected AR = ", AR_Projected)

    # Inputs
    Points_N = Points_no
    segments = segments
    # AR = 6                  # -
    Projected_area = Projected_area     # m^2
    Taper = 0.4             # -

    # First Calculations
    a_proj = np.sqrt(AR_Projected*Projected_area)/2
    b_proj = a_proj/Top_ellipse_ratio
    # print(b_proj*2/3)
    c_tip = 2/3 * b_proj
    c_root = c_tip/Taper

    y = []
    y2 = []


    # Front view Kite
    x = np.linspace(-1,1,Points_N)
    x = a_proj*x
    a = a_proj*1
    b = a_proj*1/Front_ellipse_ratio
    projected_span = np.max(x)-np.min(x)

    # Generate semi ellipse
    for i in x:
        if np.abs(i) <= x[-1]:
            y.append(1/Front_ellipse_ratio * np.sqrt(a**2-i**2))
            y2.append(-1/Front_ellipse_ratio * np.sqrt(a**2-i**2))
        else:
            y.append(0)

    # Introduce empty lists
    LE_x_segmented = []
    z_segmented = []
    y_segmented = []
    TE_x_segmented = []

    angles = []
    arclength = 0

    for i in np.arange(1,Points_N, 1):
        arclength += np.sqrt(((x[i]-x[i-1])**2 + (y[i]-y[i-1])**2))
    print('arclength', arclength)

    # Subdivide based on arclength
    arclist = []
    LE_x_segmented.append(a)
    z_segmented.append(0)
    for i in range(0,segments,1):
        angles.append(np.pi - i*np.pi/(segments-1))
        arclist.append(arclength - i*arclength/(segments-1))

    for length in arclist:
        arc = 0
        for i in range(1,len(x), 1):
            if arc < length:
                arc += np.sqrt(((x[i] - x[i - 1]) ** 2 + (y[i] - y[i - 1]) ** 2))
            else:
                LE_x_segmented.append(x[i])
                z_segmented.append(y[i])
                break

    # Subdivide based on angles
    # angles = np.array(angles[::-1])
    # for angle in angles:
    #     for i in range(len(x)):
    #         if np.arctan2(y[i], x[i]) > angle:
    #             continue
    #         else:
    #             LE_x_segmented.append(x[i])
    #             z_segmented.append(y[i])
    #             break

    # print('z coords = ', z_segmented)
    # plt.figure()
    # plt.plot(x,y,label='Front View')
    # plt.plot(LE_x_segmented, z_segmented,label='Front View Segmented')
    # plt.legend()

    span_scaling = projected_span/arclength
    print("Span Scaling = ", span_scaling)
    flat_area_span = 1/span_scaling * projected_span
    flat_area = Projected_area/span_scaling

    # Projected top view Kite
    y = []
    y2 = []

    '''# TOP view Kite'''
    x = np.linspace(-1,1,Points_N)
    x = a_proj * x
    a = a_proj
    b = a_proj *  1/Top_ellipse_ratio

    LE_x_segmented = []
    y_segmented1 = []

    for i in x:
        if np.abs(i) <= x[-1]:
            y.append(b/a * np.sqrt(a**2-i**2))
            y2.append(-b/a * np.sqrt(a**2-i**2))
        else:
            y.append(0)

    arclength = 0

    for i in np.arange(1,Points_N, 1):
        arclength += np.sqrt(((x[i]-x[i-1])**2 + (y[i]-y[i-1])**2))
    # print('arclength', arclength)

    # Subdivide based on arclength
    arclist = []
    LE_x_segmented.append(a)
    y_segmented1.append(0)
    for i in range(0,segments,1):
        arclist.append(arclength - i*arclength/(segments-1))

    for length in arclist:
        arc = 0
        for i in range(1,len(x), 1):
            if arc < length:
                arc += np.sqrt(((x[i] - x[i - 1]) ** 2 + (y[i] - y[i - 1]) ** 2))
            else:
                LE_x_segmented.append(x[i])
                y_segmented1.append(y[i])
                break


    # for angle in angles:
    #     for i in range(len(x)):
    #         if np.arctan2(y[i], x[i]) > angle:
    #             continue
    #         else:
    #             LE_x_segmented.append(x[i])
    #             y_segmented1.append(y[i])
    #             break
    # chords = []
    # for i in range(0,5):
    #     chords.append(c_tip + (c_root-c_tip)/a_proj * (a_proj-LE_x_segmented[i]))
    #
    # chords.extend(chords[::-1])
    # chords = np.array(chords)
    # print(chords)
    # print(c_tip)
    # print(c_root)

    y_segmented1 = np.array(y_segmented1)

    LE_y_segmented = -y_segmented1 - c_tip
    # print('chords = ', -np.array(LE_y_segmented))
    # print(LE_y_segmented[0]/LE_y_segmented[4])
    y_segmented = np.zeros(segments)
    y_segmented = np.array(y_segmented)
    z_segmented = np.array(z_segmented)


    # plt.figure()
    # plt.plot(LE_x_segmented, y_segmented, label='Planform View Segmented', color='r')
    # plt.plot(LE_x_segmented, -LE_y_segmented, color='r')
    # plt.plot([LE_x_segmented[0], LE_x_segmented[0]], [y_segmented[0], -LE_y_segmented[0]], color='r')
    # plt.plot([LE_x_segmented[-1], LE_x_segmented[-1]], [y_segmented[-1], -LE_y_segmented[-1]] , color='r')
    # plt.legend()
    # plt.show()

    # Convert to Uri compatible coordinates

    # coord = np.empty((21, 3))
    # coord[0, :] = [0, 0, 0]
    # coord[1, :] = [y_segmented[0], LE_x_segmented[0], z_segmented[0]]
    # for i in np.arange(2,10,1):
    #     coord[i, :] = [LE_y_segmented[i-1], LE_x_segmented[i-1], z_segmented[i-1]]
    #     coord[8 + i, :] = [y_segmented[i-2], LE_x_segmented[-i+1], z_segmented[i-2]]
    # coord[18, :] = [y_segmented[8], LE_x_segmented[-9], z_segmented[8]]
    # coord[19, :] = [LE_y_segmented[0], LE_x_segmented[0], z_segmented[0]]
    # coord[20, :] = [LE_y_segmented[0], LE_x_segmented[9], z_segmented[9]]

    coord = np.empty((2*len(LE_y_segmented), 3))
    for i in range(len(LE_y_segmented)):
        coord[2*i, 0] = LE_y_segmented[i]
        coord[2*i, 1] = LE_x_segmented[-i-1]
        coord[2*i, 2] = z_segmented[i]

        coord[2*i + 1, 0] = y_segmented[i]
        coord[2*i + 1, 1] = LE_x_segmented[-i-1]
        coord[2*i + 1, 2] = z_segmented[i]
    chords = -np.array(LE_y_segmented)

    print('Flat area = ', flat_area)
    print('Flat span = ', flat_area_span)
    print('AR_flat = ', flat_area_span ** 2 / flat_area,
          '\nCoordinate Generation Finished\n# ------------------------------ #\n')

    # print(coord)
    if Plotting:
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        ax.plot(coord[::2, 0], coord[::2, 1], coord[::2, 2], label='Leading Edge')
        ax.plot(coord[1::2, 0], coord[1::2, 1], coord[1::2, 2], label='Trailing Edge')
        plt.legend()
        plt.show()
    return coord, chords


if __name__ == "__main__":
    A_proj = 16.65
    segments = 10
    plotting = True
    Points = 500000
    Generate_Kite15_coords(A_proj, segments, Points, plotting)
