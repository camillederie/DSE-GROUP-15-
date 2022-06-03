import numpy as np
import matplotlib.pyplot as plt


# %% MATHEMATICAL FUNCTIONS

def interp1d(xnew, xarray, yarray):
    """
  Interpolate a 1D array

  """
    i1 = None
    for i in range(len(xarray) - 1):
        if (xarray[i] <= xnew):
            if (xarray[i + 1] >= xnew):
                i1 = i;
                i2 = i + 1
                # if i1 is None:        # In case it is not in the array, take last value
    #     i1 = i
    #     i2 = i-2
    result = (xnew - xarray[i1]) / (xarray[i2] - xarray[i1]) * (yarray[i2] - yarray[i1]) + yarray[i1];
    return result


def cosspace(min, max, n_points):
    """
    Create an array with cosine spacing, from min to max values, with n points

    """
    mean = (max + min) / 2
    amp = (max - min) / 2

    return mean + amp * np.cos(np.linspace(np.pi, 0, n_points))


def cross_product(r1, r2):
    """
    Cross product between r1 and r2

    """

    return np.array([r1[1] * r2[2] - r1[2] * r2[1], r1[2] * r2[0] - r1[0] * r2[2], r1[0] * r2[1] - r1[1] * r2[0]])


def vec_norm(v):
    """
    Norm of a vector

    """
    return np.sqrt(v[0] ** 2 + v[1] ** 2 + v[2] ** 2)


def dot_product(r1, r2):
    """
    Dot product between r1 and r2

    """
    return r1[0] * r2[0] + r1[1] * r2[1] + r1[2] * r2[2]


def vector_projection(v, u):
    """
    Find the projection of a vector into a direction

    Parameters
    ----------
    v : vector to be projected
    u : direction

    Returns
    -------
    proj : projection of the vector v onto u

    """
    # Inputs:
    #     u = direction vector
    #     v = vector to be projected

    unit_u = u / np.linalg.norm(u)
    proj = np.dot(v, unit_u) * unit_u

    return proj


# %% INDUCED VELOCITIES FUNCTIONS

def velocity_induced_single_ring_semiinfinite(ring, controlpoint, model, Uinf):
    """
    Calculates the velocity induced by a ring at a certain controlpoint

    Parameters
    ----------
    ring : List of dictionaries defining the filaments of a vortex ring
    controlpoint : Dictionary defining a controlpoint
    model : VSM: Vortex Step method/ LLT: Lifting Line Theory
    Uinf : Wind speed vector

    Returns
    -------
    velind : Induced velocity

    """

    velind = [0, 0, 0]
    for filament in ring:
        GAMMA = filament['Gamma']
        XV1 = filament['x1']

        XVP = controlpoint
        if filament['id'] == 'trailing_inf1':
            Vf = filament['dir']
            tempvel = velocity_3D_trailing_vortex_semiinfinite(XV1, Vf, XVP, GAMMA, Uinf)
        elif filament['id'] == 'trailing_inf2':
            Vf = filament['dir']
            tempvel = velocity_3D_trailing_vortex_semiinfinite(XV1, Vf, XVP, -GAMMA, Uinf)
        elif filament['id'] == 'bound':
            if model == 'VSM':
                XV2 = filament['x2']
                tempvel = velocity_3D_bound_vortex(XV1, XV2, XVP, GAMMA)
            else:
                tempvel = [0, 0, 0]
        else:
            XV2 = filament['x2']
            tempvel = velocity_3D_trailing_vortex(XV1, XV2, XVP, GAMMA, Uinf)

        velind[0] += tempvel[0]
        velind[1] += tempvel[1]
        velind[2] += tempvel[2]

    return velind


def update_Gamma_single_ring(ring, GammaNew, WeightNew):
    """
    Update Gamma of all the filaments in a horshoe ring

    """
    # Runs through all filaments
    for filament in ring:
        filament['Gamma'] = filament['Gamma'] * (1 - WeightNew) + WeightNew * GammaNew  # Update each Gamma

    return ring


def velocity_3D_bound_vortex(XV1, XV2, XVP, gamma):
    """
    Calculate the velocity induced by a bound vortex filament in a point in space ()

    Vortex core correction from:
        Rick Damiani et al. “A vortex step method for nonlinear airfoil polar data as implemented in
KiteAeroDyn”.

    ----------
    XV1 : Point A of Bound vortex (array)
    XV2 : Point B of Bound vortex (array)
    XVP : Control point (array)
    gamma : Strength of the vortex (scalar)

    Returns
    -------
    vel_ind : Induced velocity (array)

    """
    r0 = XV2 - XV1  # Vortex filament
    r1 = XVP - XV1  # Controlpoint to one end of the vortex filament
    r2 = XVP - XV2  # Controlpoint to one end of the vortex filament

    # Cross products used for later computations
    r1Xr0 = cross_product(r1, r0)
    r2Xr0 = cross_product(r2, r0)

    epsilon = 0.05 * vec_norm(r0)  # Cut-off radius

    if vec_norm(r1Xr0) / vec_norm(r0) > epsilon:  # Perpendicular distance from XVP to vortex filament (r0)
        r1Xr2 = cross_product(r1, r2)
        vel_ind = gamma / (4 * np.pi) * r1Xr2 / (vec_norm(r1Xr2) ** 2) * dot_product(r0,
                                                                                     r1 / vec_norm(r1) - r2 / vec_norm(
                                                                                         r2))
    else:
        # The control point is placed on the edge of the radius core
        # proj stands for the vectors respect to the new controlpoint
        r1_proj = dot_product(r1, r0) * r0 / (vec_norm(r0) ** 2) + epsilon * r1Xr0 / vec_norm(r1Xr0)
        r2_proj = dot_product(r2, r0) * r0 / (vec_norm(r0) ** 2) + epsilon * r2Xr0 / vec_norm(r2Xr0)
        r1Xr2_proj = cross_product(r1_proj, r2_proj)
        vel_ind_proj = gamma / (4 * np.pi) * r1Xr2_proj / (vec_norm(r1Xr2_proj) ** 2) * dot_product(r0,
                                                                                                    r1_proj / vec_norm(
                                                                                                        r1_proj) - r2_proj / vec_norm(
                                                                                                        r2_proj))
        vel_ind = vec_norm(r1Xr0) / (vec_norm(r0) * epsilon) * vel_ind_proj
    return vel_ind


def velocity_3D_trailing_vortex(XV1, XV2, XVP, gamma, Uinf):
    """
    Calculate the velocity induced by a trailing vortex filament in a point in space

    Vortex core correction from:
        Rick Damiani et al. “A vortex step method for nonlinear airfoil polar data as implemented in
KiteAeroDyn”.
    ----------
    XV1 : Point A of the vortex filament (array)
    XV2 : Point B of the vortex filament (array)
    XVP : Controlpoint (array)
    gamma : Strength of the vortex (scalar)
    Uinf : Inflow velocity modulus (scalar)

    Returns
    -------
    vel_ind : induced velocity by the trailing fil. (array)

    """
    r0 = XV2 - XV1  # Vortex filament
    r1 = XVP - XV1  # Controlpoint to one end of the vortex filament
    r2 = XVP - XV2  # Controlpoint to one end of the vortex filament

    alpha0 = 1.25643  # Oseen parameter
    nu = 1.48e-5  # Kinematic viscosity of air
    r_perp = dot_product(r1, r0) * r0 / (vec_norm(r0) ** 2)  # Vector from XV1 to XVP perpendicular to the core radius
    epsilon = np.sqrt(4 * alpha0 * nu * vec_norm(r_perp) / Uinf)  # Cut-off radius

    # Cross products used for later computations
    r1Xr0 = cross_product(r1, r0)
    r2Xr0 = cross_product(r2, r0)

    if vec_norm(r1Xr0) / vec_norm(r0) > epsilon:  # Perpendicular distance from XVP to vortex filament (r0)
        r1Xr2 = cross_product(r1, r2)
        vel_ind = gamma / (4 * np.pi) * r1Xr2 / (vec_norm(r1Xr2) ** 2) * dot_product(r0,
                                                                                     r1 / vec_norm(r1) - r2 / vec_norm(
                                                                                         r2))
    else:
        # The control point is placed on the edge of the radius core
        # proj stands for the vectors respect to the new controlpoint
        r1_proj = dot_product(r1, r0) * r0 / (vec_norm(r0) ** 2) + epsilon * r1Xr0 / vec_norm(r1Xr0)
        r2_proj = dot_product(r2, r0) * r0 / (vec_norm(r0) ** 2) + epsilon * r2Xr0 / vec_norm(r2Xr0)
        r1Xr2_proj = cross_product(r1_proj, r2_proj)
        vel_ind_proj = gamma / (4 * np.pi) * r1Xr2_proj / (vec_norm(r1Xr2_proj) ** 2) * dot_product(r0,
                                                                                                    r1_proj / vec_norm(
                                                                                                        r1_proj) - r2_proj / vec_norm(
                                                                                                        r2_proj))
        vel_ind = vec_norm(r1Xr0) / (vec_norm(r0) * epsilon) * vel_ind_proj
    return vel_ind


def velocity_3D_trailing_vortex_semiinfinite(XV1, Vf, XVP, GAMMA, Uinf):
    """
    Calculate the velocity induced by a semiinfinite trailing vortex filament in a point in space

    Vortex core correction from:
        Rick Damiani et al. “A vortex step method for nonlinear airfoil polar data as implemented in
KiteAeroDyn”.
    ----------
    XV1 : Point A of the vortex filament (array)
    XV2 : Point B of the vortex filament (array)
    XVP : Controlpoint (array)
    gamma : Strength of the vortex (scalar)
    Uinf : Inflow velocity modulus (scalar)

    Returns
    -------
    vel_ind : induced velocity by the trailing fil. (array)

    """
    r1 = XVP - XV1  # Vector from XV1 to XVP
    r1XVf = cross_product(r1, Vf)

    alpha0 = 1.25643  # Oseen parameter
    nu = 1.48e-5  # Kinematic viscosity of air
    r_perp = dot_product(r1, Vf) * Vf  # Vector from XV1 to XVP perpendicular to the core radius
    epsilon = np.sqrt(4 * alpha0 * nu * vec_norm(r_perp) / Uinf)  # Cut-off radius

    if vec_norm(r1XVf) / vec_norm(Vf) > epsilon:
        # determine scalar
        K = GAMMA / 4 / np.pi / vec_norm(r1XVf) ** 2 * (1 + dot_product(r1, Vf) / vec_norm(r1))
        # determine the three velocity components
        vel_ind = K * r1XVf
    else:
        r1_proj = dot_product(r1, Vf) * Vf + epsilon * (r1 / vec_norm(r1) - Vf) / vec_norm(r1 / vec_norm(r1) - Vf)
        r1XVf_proj = cross_product(r1_proj, Vf)
        K = GAMMA / 4 / np.pi / vec_norm(r1XVf_proj) ** 2 * (1 + dot_product(r1_proj, Vf) / vec_norm(r1_proj))
        # determine the three velocity components
        vel_ind = K * r1XVf_proj
    # output results, vector with the three velocity components
    return vel_ind


def velocity_induced_bound_2D(ringvec):
    r0 = ringvec['r0']
    r3 = ringvec['r3']

    cross = [r0[1] * r3[2] - r0[2] * r3[1], r0[2] * r3[0] - r0[0] * r3[2], r0[0] * r3[1] - r0[1] * r3[0]]

    ind_vel = cross / (cross[0] ** 2 + cross[1] ** 2 + cross[2] ** 2) / 2 / np.pi * np.linalg.norm(r0)

    return ind_vel


# %% SOLVER and ITERATION FUNCTIONS

def solve_lifting_line_system_matrix_approach_semiinfinite(ringvec, controlpoints,
                                                           rings, Uinf, data_airf,
                                                           conv_crit, model):
    """
    Solve the VSM or LLM by finding the distribution of Gamma

    Parameters
    ----------
    ringvec : List of dictionaries containing the vectors that define each ring
    controlpoints :List of dictionaries with the variables needed to define each wing section
    rings : List of list with the definition of each vortex filament
    Uinf : Wind speed velocity vector
    data_airf : 2D airfoil data with alpha, Cl, Cd, Cm
    recalc_alpha : True if you want to recalculate the induced angle of attack at 1/4 of the chord (VSM)
    Gamma0 : Initial Guess of Gamma
    model : VSM: Vortex Step method/ LLT: Lifting Line Theory

    Returns
    -------
    F: Lift, Drag and Moment at each section
    Gamma: Gamma at each section
    aero_coeffs: alpha, cl, cd, cm at each wing section

    """

    nocore = False  # To shut down core corrections input True
    # Initialization of the parameters
    velocity_induced = []
    u = 0;
    v = 0;
    w = 0
    Gamma = np.zeros(len(rings))
    GammaNew = np.zeros(len(rings))
    Lift = np.zeros(len(rings))
    Drag = np.zeros(len(rings))
    Ma = np.zeros(len(rings))
    alpha = np.zeros(len(rings))
    cl = np.zeros(len(rings))
    cd = np.zeros(len(rings))
    cm = np.zeros(len(rings))
    MatrixU = np.empty((len(controlpoints), len(rings)))
    MatrixV = np.empty((len(controlpoints), len(rings)))
    MatrixW = np.empty((len(controlpoints), len(rings)))
    U_2D = np.zeros((len(controlpoints), 3))

    # Number of iterations and convergence criteria
    Niterations = conv_crit['Niterations']
    errorlimit = conv_crit['error']
    ConvWeight = conv_crit['Conv_weight']

    for icp in range(len(controlpoints)):

        if model == 'VSM':
            # Velocity induced by a infinte bound vortex with Gamma = 1
            U_2D[icp] = velocity_induced_bound_2D(ringvec[icp])

        for jring in range(len(rings)):
            rings[jring] = update_Gamma_single_ring(rings[jring], 1, 1)
            # Calculate velocity induced by a ring to a control point
            velocity_induced = velocity_induced_single_ring_semiinfinite(rings[jring],
                                                                         controlpoints[icp]['coordinates'],
                                                                         model, vec_norm(Uinf))
            # If CORE corrections are deactivated
            if nocore == True:
                # Calculate velocity induced by a ring to a control point
                velocity_induced = velocity_induced_single_ring_semiinfinite_nocore(rings[jring],
                                                                                    controlpoints[icp]['coordinates'],
                                                                                    model)

            # AIC Matrix
            MatrixU[icp, jring] = velocity_induced[0]
            MatrixV[icp, jring] = velocity_induced[1]
            MatrixW[icp, jring] = velocity_induced[2]

    # Start solving iteratively
    for kiter in range(Niterations):
        for ig in range(len(Gamma)):
            Gamma[ig] = GammaNew[ig]

        for icp in range(len(controlpoints)):
            # Initialize induced velocity to 0
            u = 0;
            v = 0;
            w = 0;
            # Compute induced velocities with previous Gamma distribution
            for jring in range(len(rings)):
                u = u + MatrixU[icp][jring] * Gamma[jring];  # x-component of velocity
                v = v + MatrixV[icp][jring] * Gamma[jring];  # y-component of velocity
                w = w + MatrixW[icp][jring] * Gamma[jring];  # z-component of velocity

            u = u - U_2D[icp, 0] * Gamma[icp]
            v = v - U_2D[icp, 1] * Gamma[icp]
            w = w - U_2D[icp, 2] * Gamma[icp]

            # Calculate relative velocity and angle of attack
            Urel = [Uinf[0] + u, Uinf[1] + v, Uinf[2] + w]
            vn = dot_product(controlpoints[icp]['normal'], Urel)
            vtan = dot_product(controlpoints[icp]['tangential'], Urel)
            alpha[icp] = np.arctan(vn / vtan)

            # Calculate terms of induced corresponding to the z-airfoil direction
            z_airf = controlpoints[icp]['airf_coord'][:, 2]
            Urelcrossz = np.cross(Urel, z_airf)
            Umag = np.linalg.norm(Urelcrossz)
            Uinfcrossz = np.cross(Uinf, ringvec[icp]['r0'] / np.linalg.norm(ringvec[icp]['r0']))
            Umagw = np.linalg.norm(Uinfcrossz)

            # Look-up airfoil 2D coefficients
            if data_airf.ndim == 3:
                cl[icp], cd[icp], cm[icp] = airfoil_coeffs(alpha[icp], data_airf[:, :, icp])
            else:
                cl[icp], cd[icp], cm[icp] = airfoil_coeffs(alpha[icp], data_airf)

            # Retrieve forces and moments
            Lift[icp] = 0.5 * Umag ** 2 * cl[icp] * controlpoints[icp]['chord']
            Drag[icp] = 0.5 * Umag ** 2 * cd[icp] * controlpoints[icp]['chord']
            Ma[icp] = 0.5 * Umag ** 2 * cm[icp] * controlpoints[icp]['chord'] ** 2

            # Find the new gamma using Kutta-Joukouski law
            GammaNew[icp] = 0.5 * Umag ** 2 / Umagw * cl[icp] * controlpoints[icp]['chord']

        # check convergence of solution
        refererror = np.amax(np.abs(GammaNew));
        refererror = np.amax([refererror, 0.001]);  # define scale of bound circulation
        error = np.amax(np.abs(GammaNew - Gamma));  # difference betweeen iterations
        error = error / refererror;  # relative error
        if (error < errorlimit):
            # if error smaller than limit, stop iteration cycle
            print('Converged N = ' + str(kiter))
            break

        # If not converged, apply convergence weighting and continue the iteration
        for ig in range(len(Gamma)):
            GammaNew[ig] = (1 - ConvWeight) * Gamma[ig] + ConvWeight * GammaNew[ig]

    # In case VSM, calculate the effective angle of attack at a 1/4 chord
    if model == 'VSM':

        for ig in range(len(Gamma)):
            Gamma[ig] = GammaNew[ig]
        for icp in range(len(controlpoints)):
            # Compute induced velocities at 1/4 chord
            for jring in range(len(rings)):
                rings[jring] = update_Gamma_single_ring(rings[jring], 1, 1)
                velocity_induced = velocity_induced_single_ring_semiinfinite(rings[jring],
                                                                             controlpoints[icp]['coordinates_aoa'],
                                                                             'LLT', vec_norm(Uinf))
                if nocore == True:
                    velocity_induced = velocity_induced_single_ring_semiinfinite_nocore(rings[jring],
                                                                                        controlpoints[icp][
                                                                                            'coordinates_aoa'],
                                                                                        model)
                MatrixU[icp, jring] = velocity_induced[0]
                MatrixV[icp, jring] = velocity_induced[1]
                MatrixW[icp, jring] = velocity_induced[2]

        for icp in range(len(controlpoints)):
            u = 0;
            v = 0;
            w = 0;
            for jring in range(len(rings)):
                u = u + MatrixU[icp][jring] * Gamma[jring];  # x-component of velocity
                v = v + MatrixV[icp][jring] * Gamma[jring];  # y-component of velocity
                w = w + MatrixW[icp][jring] * Gamma[jring];  # z-component of velocity

            Urel = [Uinf[0] + u, Uinf[1] + v, Uinf[2] + w]
            vn = np.dot(controlpoints[icp]['normal'], Urel)
            vtan = np.dot(controlpoints[icp]['tangential'], Urel)
            # New relative angle of attack
            alpha[icp] = np.arctan(vn / vtan)

    aero_coeffs = np.column_stack([alpha, cl, cd, cm])
    F = np.column_stack([Lift, Drag, Ma])

    return F, Gamma, aero_coeffs


def output_results(Fmag, aero_coeffs, ringvec, Uinf, controlpoints, Atot):
    """
    Post-process results to get global forces and aerodynamic coefficients

    Parameters
    ----------
    Fmag : Lift, Drag and Moment magnitudes
    aero_coeffs : alpha, cl, cd, cm
    ringvec : List of dictionaries containing the vectors that define each ring
    Uinf : Wind speed velocity vector
    controlpoints : List of dictionaries with the variables needed to define each wing section
    Atot : Planform area

    Returns
    -------
    F_rel : Lift and drag forces relative to the local angle of attack
    F_gl : Lift and drag forces relative to the wind direction
    Ltot : Total lift
    Dtot : Total drag
    CL : Global CL
    CD : Global CD

    """

    alpha = aero_coeffs[:, 0]
    F_rel = []
    F_gl = []
    Fmag_gl = []
    Ltot = 0
    Dtot = 0
    for i in range(len(alpha)):
        r0 = ringvec[i]['r0']
        # Relative wind speed direction
        dir_urel = np.cos(alpha[i]) * controlpoints[i]['tangential'] + np.sin(alpha[i]) * controlpoints[i]['normal']
        dir_urel = dir_urel / np.linalg.norm(dir_urel)
        # Lift direction relative to Urel
        dir_L = np.cross(dir_urel, r0)
        dir_L = dir_L / np.linalg.norm(dir_L)
        # Drag direction relative to Urel
        dir_D = np.cross([0, 1, 0], dir_L)
        dir_D = dir_D / np.linalg.norm(dir_D)
        # Lift and drag relative to Urel
        L_rel = dir_L * Fmag[i, 0]
        D_rel = dir_D * Fmag[i, 1]
        F_rel.append([L_rel, D_rel])
        # Lift direction relative to the wind speed
        dir_L_gl = np.cross(Uinf, [0, 1, 0])
        dir_L_gl = dir_L_gl / vec_norm(dir_L_gl)
        # Lift and drag relative to the windspeed
        L_gl = vector_projection(L_rel, dir_L_gl) + vector_projection(D_rel, dir_L_gl)
        D_gl = vector_projection(L_rel, Uinf) + vector_projection(D_rel, Uinf)
        F_gl.append([L_gl, D_gl])
        Fmag_gl.append([dot_product(L_rel, dir_L_gl) + dot_product(D_rel, dir_L_gl) \
                           , dot_product(L_rel, Uinf / vec_norm(Uinf)) + dot_product(D_rel, Uinf / vec_norm(Uinf))])

    # Calculate total aerodynamic forces
    Lift = []
    for i in range(len(Fmag_gl)):
        Ltot += Fmag_gl[i][0] * np.linalg.norm(ringvec[i]['r0'])
        Lift.append(Fmag_gl[i][0] * np.linalg.norm(ringvec[i]['r0']))
        Dtot += Fmag_gl[i][1] * np.linalg.norm(ringvec[i]['r0'])

    Umag = np.linalg.norm(Uinf)
    CL = Ltot / (0.5 * Umag ** 2 * Atot)
    CD = Dtot / (0.5 * Umag ** 2 * Atot)

    return F_rel, F_gl, Ltot, Dtot, CL, CD, Lift


def create_geometry_general(coordinates, Uinf, N, ring_geo, model):
    """
    Create geometry structures necessary for solving the system of circualtion

    Parameters
    ----------
    coordinates : coordinates the nodes (each section is defined by two nodes,
                                         the first is the LE, so each section
                                         defined by a pair of coordinates)
    Uinf : Wind speed vector
    N : Number of sections
    ring_geo :  - '3fil': Each horsehoe is defined by 3 filaments
                - '5fil': Each horseshoe is defined by 5 filaments
    model : VSM: Vortex Step method/ LLT: Lifting Line Theory

    Returns
    -------
    controlpoints :  List of dictionaries with the variables needed to define each wing section
    rings : List of list with the definition of each vortex filament
    wingpanels : List with the points defining each wing pannel
    ringvec : List of dictionaries containing the vectors that define each ring
    coord_L : coordinates of the aerodynamic centers of each wing panel

   """

    filaments = []
    controlpoints = []
    rings = []
    wingpanels = []
    ringvec = []
    coord_L = []

    # Go through all wing panels
    for i in range(N - 1):

        # Identify points defining the panel
        section = {
            'p1': coordinates[2 * i, :],
            'p2': coordinates[2 * i + 2, :],
            'p3': coordinates[2 * i + 3, :],
            'p4': coordinates[2 * i + 1, :]
        }
        wingpanels.append(section)

        di = vec_norm(coordinates[2 * i, :] * 0.75 + coordinates[2 * i + 1, :] * 0.25 - (
                    coordinates[2 * i + 2, :] * 0.75 + coordinates[2 * i + 3, :] * 0.25))
        if i == 0:
            diplus = vec_norm(coordinates[2 * (i + 1), :] * 0.75 + coordinates[2 * (i + 1) + 1, :] * 0.25 - (
                        coordinates[2 * (i + 1) + 2, :] * 0.75 + coordinates[2 * (i + 1) + 3, :] * 0.25))
            ncp = di / (di + diplus)
        elif i == N - 2:
            dimin = vec_norm(coordinates[2 * (i - 1), :] * 0.75 + coordinates[2 * (i - 1) + 1, :] * 0.25 - (
                        coordinates[2 * (i - 1) + 2, :] * 0.75 + coordinates[2 * (i - 1) + 3, :] * 0.25))
            ncp = dimin / (dimin + di)
        else:
            dimin = vec_norm(coordinates[2 * (i - 1), :] * 0.75 + coordinates[2 * (i - 1) + 1, :] * 0.25 - (
                        coordinates[2 * (i - 1) + 2, :] * 0.75 + coordinates[2 * (i - 1) + 3, :] * 0.25))
            diplus = vec_norm(coordinates[2 * (i + 1), :] * 0.75 + coordinates[2 * (i + 1) + 1, :] * 0.25 - (
                        coordinates[2 * (i + 1) + 2, :] * 0.75 + coordinates[2 * (i + 1) + 3, :] * 0.25))
            ncp = 0.25 * (dimin / (dimin + di) + di / (di + diplus) + 1)

        ncp = 1 - ncp
        chord = np.linalg.norm((section['p2'] + section['p1']) / 2 - (section['p3'] + section['p4']) / 2)
        LLpoint = (section['p2'] * (1 - ncp) + section['p1'] * ncp) * 3 / 4 + (
                    section['p3'] * (1 - ncp) + section['p4'] * ncp) * 1 / 4
        VSMpoint = (section['p2'] * (1 - ncp) + section['p1'] * ncp) * 1 / 4 + (
                    section['p3'] * (1 - ncp) + section['p4'] * ncp) * 3 / 4
        coord_L.append(LLpoint)

        # Define bound vortex filament
        bound = {'id': 'bound',
                 'x1': section['p1'] * 3 / 4 + section['p4'] * 1 / 4,
                 'x2': section['p2'] * 3 / 4 + section['p3'] * 1 / 4,
                 'Gamma': 0}
        filaments.append(bound)

        x_airf = np.cross(VSMpoint - LLpoint, section['p2'] - section['p1'])
        x_airf = x_airf / np.linalg.norm(x_airf)
        y_airf = VSMpoint - LLpoint
        y_airf = y_airf / np.linalg.norm(y_airf)
        z_airf = bound['x2'] - bound['x1']
        # z_airf[0] = 0
        z_airf = z_airf / np.linalg.norm(z_airf)
        airf_coord = np.column_stack([x_airf, y_airf, z_airf])

        normal = x_airf
        tangential = y_airf
        if model == 'VSM':
            cp = {'coordinates': VSMpoint, 'chord': chord,
                  'normal': normal,
                  'tangential': tangential,
                  'airf_coord': airf_coord,
                  'coordinates_aoa': LLpoint}
            controlpoints.append(cp)
        elif model == 'LLT':

            cp = {'coordinates': LLpoint, 'chord': chord,
                  'normal': normal,
                  'tangential': tangential,
                  'airf_coord': airf_coord}
            controlpoints.append(cp)

        temp = {'r0': bound['x2'] - bound['x1'],
                'r1': cp['coordinates'] - bound['x1'],
                'r2': cp['coordinates'] - bound['x2'],
                'r3': cp['coordinates'] - (bound['x2'] + bound['x1']) / 2}
        ringvec.append(temp)

        temp = Uinf / np.linalg.norm(Uinf)
        if ring_geo == '3fil':
            # create trailing filaments, at x1 of bound filament
            temp1 = {'dir': temp, 'id': 'trailing_inf1',
                     'x1': bound['x1'], 'Gamma': 0}
            filaments.append(temp1)

            # create trailing filaments, at x2 of bound filament
            temp1 = {'x1': bound['x2'],
                     'dir': temp, 'id': 'trailing_inf2', 'Gamma': 0}
            filaments.append(temp1)
        elif ring_geo == '5fil':
            temp1 = {'x1': section['p4'],
                     'x2': bound['x1'], 'Gamma': 0,
                     'id': 'trailing1'}
            filaments.append(temp1)

            temp1 = {'dir': temp, 'id': 'trailing_inf1',
                     'x1': section['p4'], 'Gamma': 0}
            filaments.append(temp1)

            # create trailing filaments, at x2 of bound filament
            temp1 = {'x2': section['p3'],
                     'x1': bound['x2'], 'Gamma': 0,
                     'id': 'trailing1'}
            filaments.append(temp1)

            temp1 = {'x1': section['p3'],
                     'dir': temp, 'id': 'trailing_inf2', 'Gamma': 0}
            filaments.append(temp1)

        #

        rings.append(filaments)
        filaments = []

    coord_L = np.array(coord_L)
    return controlpoints, rings, wingpanels, ringvec, coord_L


def create_geometry_LEI(coordinates, Uinf, N, ring_geo, model):
    filaments = []
    controlpoints = []
    rings = []
    wingpanels = []
    ringvec = []
    coord_L = []
    for i in range(N - 1):

        section = {
            'p1': coordinates[2 * i, :],
            'p2': coordinates[2 * i + 2, :],
            'p3': coordinates[2 * i + 3, :],
            'p4': coordinates[2 * i + 1, :]
        }
        wingpanels.append(section)
        di = vec_norm(coordinates[2 * i, :] * 0.75 + coordinates[2 * i + 1, :] * 0.25 - (
                    coordinates[2 * i + 2, :] * 0.75 + coordinates[2 * i + 3, :] * 0.25))
        if i == 0:
            diplus = vec_norm(coordinates[2 * (i + 1), :] * 0.75 + coordinates[2 * (i + 1) + 1, :] * 0.25 - (
                        coordinates[2 * (i + 1) + 2, :] * 0.75 + coordinates[2 * (i + 1) + 3, :] * 0.25))
            ncp = di / (di + diplus)
        elif i == N - 2:
            dimin = vec_norm(coordinates[2 * (i - 1), :] * 0.75 + coordinates[2 * (i - 1) + 1, :] * 0.25 - (
                        coordinates[2 * (i - 1) + 2, :] * 0.75 + coordinates[2 * (i - 1) + 3, :] * 0.25))
            ncp = dimin / (dimin + di)
        else:
            dimin = vec_norm(coordinates[2 * (i - 1), :] * 0.75 + coordinates[2 * (i - 1) + 1, :] * 0.25 - (
                        coordinates[2 * (i - 1) + 2, :] * 0.75 + coordinates[2 * (i - 1) + 3, :] * 0.25))
            diplus = vec_norm(coordinates[2 * (i + 1), :] * 0.75 + coordinates[2 * (i + 1) + 1, :] * 0.25 - (
                        coordinates[2 * (i + 1) + 2, :] * 0.75 + coordinates[2 * (i + 1) + 3, :] * 0.25))
            ncp = 0.25 * (dimin / (dimin + di) + di / (di + diplus) + 1)

        ncp = 1 - ncp
        chord = np.linalg.norm((section['p2'] + section['p1']) / 2 - (section['p3'] + section['p4']) / 2)
        LLpoint = (section['p2'] * (1 - ncp) + section['p1'] * ncp) * 3 / 4 + (
                    section['p3'] * (1 - ncp) + section['p4'] * ncp) * 1 / 4
        VSMpoint = (section['p2'] * (1 - ncp) + section['p1'] * ncp) * 1 / 4 + (
                    section['p3'] * (1 - ncp) + section['p4'] * ncp) * 3 / 4
        coord_L.append(LLpoint)

        # Define bound vortex filament
        bound = {'id': 'bound',
                 'x1': section['p1'] * 3 / 4 + section['p4'] * 1 / 4,
                 'x2': section['p2'] * 3 / 4 + section['p3'] * 1 / 4,
                 'Gamma': 0}
        filaments.append(bound)

        x_airf = np.cross(VSMpoint - LLpoint, section['p2'] - section['p1'])
        x_airf = x_airf / np.linalg.norm(x_airf)
        y_airf = VSMpoint - LLpoint
        y_airf = y_airf / np.linalg.norm(y_airf)
        z_airf = bound['x2'] - bound['x1']
        # z_airf[0] = 0
        z_airf = z_airf / np.linalg.norm(z_airf)
        airf_coord = np.column_stack([x_airf, y_airf, z_airf])

        normal = x_airf
        tangential = y_airf
        if model == 'VSM':
            cp = {'coordinates': VSMpoint, 'chord': chord,
                  'normal': normal,
                  'tangential': tangential,
                  'airf_coord': airf_coord,
                  'coordinates_aoa': LLpoint}
            controlpoints.append(cp)
        elif model == 'LLT':

            cp = {'coordinates': LLpoint, 'chord': chord,
                  'normal': normal,
                  'tangential': tangential,
                  'airf_coord': airf_coord}
            controlpoints.append(cp)

        temp = {'r0': bound['x2'] - bound['x1'],
                'r1': cp['coordinates'] - bound['x1'],
                'r2': cp['coordinates'] - bound['x2'],
                'r3': cp['coordinates'] - (bound['x2'] + bound['x1']) / 2}
        ringvec.append(temp)

        temp = Uinf / np.linalg.norm(Uinf)
        if ring_geo == '3fil':
            # create trailing filaments, at x1 of bound filament
            temp1 = {'dir': temp, 'id': 'trailing_inf1',
                     'x1': bound['x1'], 'Gamma': 0}
            filaments.append(temp1)

            # create trailing filaments, at x2 of bound filament
            temp1 = {'x1': bound['x2'],
                     'dir': temp, 'id': 'trailing_inf2', 'Gamma': 0}
            filaments.append(temp1)
        elif ring_geo == '5fil':
            temp1 = {'x1': section['p4'],
                     'x2': bound['x1'], 'Gamma': 0,
                     'id': 'trailing1'}
            filaments.append(temp1)

            temp1 = {'dir': temp, 'id': 'trailing_inf1',
                     'x1': section['p4'], 'Gamma': 0}
            filaments.append(temp1)

            # create trailing filaments, at x2 of bound filament
            temp1 = {'x2': section['p3'],
                     'x1': bound['x2'], 'Gamma': 0,
                     'id': 'trailing1'}
            filaments.append(temp1)

            temp1 = {'x1': section['p3'],
                     'dir': temp, 'id': 'trailing_inf2', 'Gamma': 0}
            filaments.append(temp1)

        #

        rings.append(filaments)
        filaments = []

    coord_L = np.array(coord_L)
    return controlpoints, rings, wingpanels, ringvec, coord_L


# %% 2D airfoil functions

def airfoil_coeffs(alpha, coeffs):
    cl = np.interp(alpha * 180 / np.pi, coeffs[:, 0], coeffs[:, 1])
    cd = np.interp(alpha * 180 / np.pi, coeffs[:, 0], coeffs[:, 2])
    cm = np.interp(alpha * 180 / np.pi, coeffs[:, 0], coeffs[:, 3])

    return cl, cd, cm


def LEI_airf_coeff(t, k, alpha):
    """
    ----------
    t : TYPE
        DESCRIPTION.
    k : TYPE
        DESCRIPTION.
    alpha : TYPE
        DESCRIPTION.

    Returns
    -------
    Cl : TYPE
        DESCRIPTION.
    Cd : TYPE
        DESCRIPTION.
    Cm : TYPE
        DESCRIPTION.

    """
    C20 = -0.008011
    C21 = -0.000336
    C22 = 0.000992
    C23 = 0.013936
    C24 = -0.003838
    C25 = -0.000161
    C26 = 0.001243
    C27 = -0.009288
    C28 = -0.002124
    C29 = 0.012267
    C30 = -0.002398
    C31 = -0.000274
    C32 = 0
    C33 = 0
    C34 = 0
    C35 = -3.371000
    C36 = 0.858039
    C37 = 0.141600
    C38 = 7.201140
    C39 = -0.676007
    C40 = 0.806629
    C41 = 0.170454
    C42 = -0.390563
    C43 = 0.101966
    C44 = 0.546094
    C45 = 0.022247
    C46 = -0.071462
    C47 = -0.006527
    C48 = 0.002733
    C49 = 0.000686
    C50 = 0.123685
    C51 = 0.143755
    C52 = 0.495159
    C53 = -0.105362
    C54 = 0.033468
    C55 = -0.284793
    C56 = -0.026199
    C57 = -0.024060
    C58 = 0.000559
    C59 = -1.787703
    C60 = 0.352443
    C61 = -0.839323
    C62 = 0.137932

    S9 = C20 * t ** 2 + C21 * t + C22
    S10 = C23 * t ** 2 + C24 * t + C25
    S11 = C26 * t ** 2 + C27 * t + C28
    S12 = C29 * t ** 2 + C30 * t + C31
    S13 = C32 * t ** 2 + C33 * t + C34
    S14 = C35 * t ** 2 + C36 * t + C37
    S15 = C38 * t ** 2 + C39 * t + C40
    S16 = C41 * t ** 2 + C42 * t + C43

    lambda5 = S9 * k + S10
    lambda6 = S11 * k + S12
    lambda7 = S13 * k + S14
    lambda8 = S15 * k + S16

    Cl = lambda5 * alpha ** 3 + lambda6 * alpha ** 2 + lambda7 * alpha + lambda8
    Cd = ((C44 * t + C45) * k ** 2 + (C46 * t + C47) * k + (C48 * t + C49)) * alpha ** 2 + \
         (C50 * t + C51) * k + (C52 * t ** 2 + C53 * t + C54)
    Cm = ((C55 * t + C56) * k + (C57 * t + C58)) * alpha ** 2 + (C59 * t + C60) * k + (C61 * t + C62)

    if alpha > 20 or alpha < -20:
        Cl = 2 * np.cos(alpha * np.pi / 180) * np.sin(alpha * np.pi / 180) ** 2
        Cd = 2 * np.sin(alpha * np.pi / 180) ** 3

    return Cl, Cd, Cm


# %% PLOT FUNCTIONS

def plot_geometry(wingpanels, controlpoints, rings, F, coord_L, plot):
    if plot == 'True':
        plt.figure()
        ax = plt.axes(projection='3d')

        for panel in wingpanels:
            coord = np.array([panel['p1'], panel['p2'], panel['p3'], panel['p4']])
            ax.plot3D(coord[:, 0], coord[:, 1], coord[:, 2], 'blue')
        # for cp in controlpoints:
        #     ax.plot3D(cp['coordinates'][0], cp['coordinates'][1], cp['coordinates'][2], '*')

        for ring in rings:
            for filament in ring:

                if filament['id'] == 'trailing_inf1' or filament['id'] == 'trailing_inf2':
                    coord = np.array([filament['x1'], filament['x1'] + filament['dir'] * 2])
                    # ax.plot3D(coord[:, 0], coord[:, 1], coord[:, 2], 'gray')
                else:
                    coord = np.array([filament['x1'], filament['x2']])
                    # ax.plot3D(coord[:, 0], coord[:, 1], coord[:, 2])
        # ax.plot3D([-1.048, 0], [-4.716, -4.716], [0, 0], 'blue')
        ax.set_xlabel("x [m]")
        ax.set_ylabel("y [m]")
        ax.set_zlabel("z [m]")
        ax.set_box_aspect([ub - lb for lb, ub in (getattr(ax, f'get_{a}lim')() for a in 'xyz')])

        # for i in range(len(F)):
        # coord = np.array([coord_L[i],coord_L[i]+F[i][0]])
        # ax.plot3D(coord[:,0], coord[:,1], coord[:,2], 'green')
        # coord = np.array([coord_L[i],coord_L[i]++F[i][1]])
        # ax.plot3D(coord[:,0], coord[:,1], coord[:,2], 'red')

        # plt.xlim((-0.25,1.5))


# %% GEOMETRY DEFINITION FUNCTIONS

def generate_coordinates_rect_wing(chord, span, twist, beta, N, dist):
    coord = np.empty((2 * N, 3))
    if dist == 'cos':
        span = cosspace(-span / 2, span / 2, N)
    elif dist == 'lin':
        span = np.linspace(-span / 2, span / 2, N)

    for i in range(N):
        coord[2 * i, :] = np.array([-0 * chord[i] * np.cos(twist[i]), span[i],
                                    0 * chord[i] * np.sin(twist[i]) - abs(span[i] * np.sin(beta[i]))])
        coord[2 * i + 1, :] = np.array([1 * chord[i] * np.cos(twist[i]), span[i],
                                        -1 * chord[i] * np.sin(twist[i]) - abs(span[i] * np.sin(beta[i]))])

    return coord


def generate_coordinates_curved_wing(chord, span, theta, R, N, dist):
    coord = np.empty((2 * N, 3))
    if dist == 'cos':
        theta = cosspace(-theta, theta, N)
    elif dist == 'lin':
        theta = np.linspace(-theta, theta, N)
    elif dist == 'cos2':
        theta1 = cosspace(-theta, -theta / N / 10, int(N / 2))
        theta2 = cosspace(theta / N / 10, theta, int(N / 2))
        theta = np.concatenate((theta1, theta2))

    for i in range(N):
        coord[2 * i, :] = np.array([0, R * np.sin(theta[i]), R * np.cos(theta[i])])
        coord[2 * i + 1, :] = np.array([chord, R * np.sin(theta[i]), R * np.cos(theta[i])])

    return coord


def generate_coordinates_el_wing(max_chord, span, N, dist):
    coord = np.empty((2 * N, 3))
    start = span * 1e-5
    if dist == 'cos':
        y_arr = cosspace(-span / 2 + start, span / 2 - start, N)
    elif dist == 'lin':
        y_arr = np.linspace(-span / 2 + start, span / 2 - start, N)

    c_arr = 2 * np.sqrt(1 - (y_arr / (span / 2)) ** 2) * max_chord / 2

    for i in range(N):
        coord[2 * i, :] = [-0.25 * c_arr[i], y_arr[i], 0]
        coord[2 * i + 1, :] = [0.75 * c_arr[i], y_arr[i], 0]

    return coord


def create_geometry_multiplanewing(span, h, chord, Uinf, N, ring_geo, model):
    coord = np.empty((2 * N, 3))

    span = np.linspace(-span / 2, span / 2, N)

    for i in range(N):
        coord[2 * i, :] = np.array([-0.25 * chord, span[i], 0])
        coord[2 * i + 1, :] = np.array([0.75 * chord, span[i], 0])

    controlpoints1, rings1, wingpanels1, ringvec1, coord_L = create_geometry_general(
        coord, Uinf, N, ring_geo, model)

    coord = np.empty((2 * N, 3))
    for i in range(N):
        coord[2 * i, :] = np.array([-0.25 * chord, span[i], h])
        coord[2 * i + 1, :] = np.array([0.75 * chord, span[i], h])

    controlpoints2, rings2, wingpanels2, ringvec2, coord_L = create_geometry_general(
        coord, Uinf, N, ring_geo, model)

    controlpoints = controlpoints1 + controlpoints2
    rings = rings1 + rings2
    wingpanels = wingpanels1 + wingpanels2
    ringvec = ringvec1 + ringvec2

    return controlpoints, rings, wingpanels, ringvec, coord_L


def create_geometry_boxwing(b, h, chord, Uinf, N, ring_geo, model):
    coord = np.empty((2 * N, 3))

    span = np.linspace(-b / 2, b / 2, N)

    for i in range(N):
        coord[2 * i, :] = np.array([-0.25 * chord, span[i], 0])
        coord[2 * i + 1, :] = np.array([0.75 * chord, span[i], 0])

    controlpoints1, rings1, wingpanels1, ringvec1, coord_L = create_geometry_general(
        coord, Uinf, N, ring_geo, model)

    height = np.linspace(h, 0, int(N * 0.2))
    coord = np.empty((int(2 * N / 5), 3))
    for i in range(int(N * 0.2)):
        coord[2 * i, :] = np.array([-0.25 * chord, b / 2, height[i]])
        coord[2 * i + 1, :] = np.array([0.75 * chord, b / 2, height[i]])

    controlpoints2, rings2, wingpanels2, ringvec2, coord_L = create_geometry_general(
        coord, Uinf, int(N * 0.2), ring_geo, model)

    span = np.linspace(-b / 2, b / 2, N)
    coord = np.empty((2 * N, 3))
    for i in range(N):
        coord[2 * i, :] = np.array([-0.25 * chord, span[i], h])
        coord[2 * i + 1, :] = np.array([0.75 * chord, span[i], h])

    controlpoints3, rings3, wingpanels3, ringvec3, coord_L = create_geometry_general(
        coord, Uinf, N, ring_geo, model)

    height = np.linspace(h, 0, int(N * 0.2))
    coord = np.empty((int(2 * N / 5), 3))
    for i in range(int(N * 0.2)):
        coord[2 * i, :] = np.array([-0.25 * chord, -b / 2, height[i]])
        coord[2 * i + 1, :] = np.array([0.75 * chord, -b / 2, height[i]])

    controlpoints4, rings4, wingpanels4, ringvec4, coord_L = create_geometry_general(
        coord, Uinf, int(N * 0.2), ring_geo, model)

    controlpoints = controlpoints1 + controlpoints2 + controlpoints3 + controlpoints4
    rings = rings1 + rings2 + rings3 + rings4
    wingpanels = wingpanels1 + wingpanels2 + wingpanels3 + wingpanels4
    ringvec = ringvec1 + ringvec2 + ringvec3 + ringvec4

    return controlpoints, rings, wingpanels, ringvec, coord_L


def generate_coordinates_swept_wing(chord, offset, span, twist, beta, N, dist):
    gamma = np.arctan(offset / (span / 2))
    coord = np.empty((2 * N, 3))
    if dist == 'cos':
        span = cosspace(-span / 2, span / 2, N)
    elif dist == 'lin':
        span = np.linspace(-span / 2, span / 2, N)

    for i in range(N):
        coord[2 * i, :] = np.array([-0.25 * chord[i] * np.cos(twist[i]) + abs(span[i]) * np.tan(gamma), span[i],
                                    0.25 * chord[i] * np.sin(twist[i]) - abs(span[i] * np.sin(beta[i]))])
        coord[2 * i + 1, :] = np.array([0.75 * chord[i] * np.cos(twist[i]) + abs(span[i]) * np.tan(gamma), span[i],
                                        -0.75 * chord[i] * np.sin(twist[i]) - abs(span[i] * np.sin(beta[i]))])

    return coord


def get_CAD_matching_jelle():
    CAD_matching_pos_values = np.array([[-1.20214652e-16, -6.38069517e-17, 0.00000000e+00],
                                        [1.36553385e+03, 4.14236984e+03, 7.09946855e+03],
                                        [4.76961681e+02, 3.95206256e+03, 8.00551817e+03],
                                        [2.30478291e+02, 3.06880325e+03, 8.90485914e+03],
                                        [6.88394771e+01, 1.94187448e+03, 9.54755529e+03],
                                        [-1.08570395e-15, 6.61184742e+02, 9.80037146e+03],
                                        [1.97694803e-15, -6.61184742e+02, 9.80036480e+03],
                                        [6.88430793e+01, -1.94187280e+03, 9.54754337e+03],
                                        [2.30488064e+02, -3.06880550e+03, 8.90485558e+03],
                                        [4.76981710e+02, -3.95206989e+03, 8.00552253e+03],
                                        [1.36556870e+03, -4.14235323e+03, 7.09948160e+03],
                                        [1.80144187e+03, -3.94674601e+03, 8.02389506e+03],
                                        [2.05717785e+03, -3.13539849e+03, 9.05882988e+03],
                                        [2.14511732e+03, -1.97920884e+03, 9.69846947e+03],
                                        [2.19982664e+03, -6.65928192e+02, 9.84278303e+03],
                                        [2.20022427e+03, 6.65926335e+02, 9.84273675e+03],
                                        [2.14511591e+03, 1.97920888e+03, 9.69845624e+03],
                                        [2.05717029e+03, 3.13539331e+03, 9.05880827e+03],
                                        [1.80142229e+03, 3.94673112e+03, 8.02386871e+03],
                                        [8.18537746e+02, 4.12869016e+03, 7.09106489e+03],
                                        [8.18570922e+02, -4.12875938e+03, 7.09108662e+03],
                                        [2.38494496e+02, -1.07895023e-02, 1.07401457e+03],
                                        [6.58990189e+02, 1.20894080e+03, 2.98123837e+03],
                                        [1.20275287e+03, 1.64585189e+03, 5.41623819e+03],
                                        [1.55498756e+03, 2.90548261e+03, 6.93216120e+03],
                                        [1.59398956e+03, 1.38635885e+03, 7.16856870e+03],
                                        [6.58961647e+02, -1.20896114e+03, 2.98123627e+03],
                                        [1.20269368e+03, -1.64592536e+03, 5.41623253e+03],
                                        [1.55496820e+03, -2.90550736e+03, 6.93218692e+03],
                                        [1.59385914e+03, -1.38646987e+03, 7.16858435e+03],
                                        [8.81276486e+01, 1.76230705e+03, 5.22297239e+03],
                                        [1.81420157e+02, 2.71347427e+03, 6.73271250e+03],
                                        [6.34674771e+01, 1.49914947e+03, 7.39093555e+03],
                                        [8.81345703e+01, -1.76232753e+03, 5.22296592e+03],
                                        [1.81431606e+02, -2.71349494e+03, 6.73270580e+03],
                                        [6.34713584e+01, -1.49913942e+03, 7.39092539e+03]])

    return CAD_matching_pos_values


def get_CAD_matching_uri():
    CAD_matching_pos_values = np.array([[0.000000E+00, 0.000000E+00, 0.000000E+00],
                                        [1.538773E+03, 4.113307E+03, 5.530496E+03],
                                        [-1.762300E+01, 3.967978E+03, 6.471622E+03],
                                        [-2.376830E+02, 3.134335E+03, 7.476759E+03],
                                        [-3.837330E+02, 1.959729E+03, 8.078914E+03],
                                        [-4.565620E+02, 6.642520E+02, 8.339101E+03],
                                        [-4.565620E+02, -6.642520E+02, 8.339101E+03],
                                        [-3.837330E+02, -1.959729E+03, 8.078914E+03],
                                        [-2.376830E+02, -3.134335E+03, 7.476759E+03],
                                        [-1.762300E+01, -3.967978E+03, 6.471622E+03],
                                        [1.538773E+03, -4.113307E+03, 5.530496E+03],
                                        [1.703467E+03, -3.955506E+03, 6.467819E+03],
                                        [2.002516E+03, -3.116753E+03, 7.454254E+03],
                                        [2.110145E+03, -1.946587E+03, 8.041739E+03],
                                        [2.158559E+03, -6.607210E+02, 8.294064E+03],
                                        [2.158559E+03, 6.607210E+02, 8.294064E+03],
                                        [2.110145E+03, 1.946587E+03, 8.041739E+03],
                                        [2.002516E+03, 3.116753E+03, 7.454254E+03],
                                        [1.703467E+03, 3.955506E+03, 6.467819E+03],
                                        [8.595800E+02, 4.139660E+03, 5.654227E+03],
                                        [8.595800E+02, -4.139660E+03, 5.654227E+03],
                                        [3.138480E+02, 0.000000E+00, 1.252129E+03],
                                        [3.875460E+02, 6.909170E+02, 1.479174E+03],
                                        [1.053321E+03, 1.772499E+03, 4.343864E+03],
                                        [1.441205E+03, 2.708913E+03, 5.601785E+03],
                                        [1.528031E+03, 1.338349E+03, 5.966178E+03],
                                        [3.875460E+02, -6.909170E+02, 1.479174E+03],
                                        [1.053321E+03, -1.772499E+03, 4.343864E+03],
                                        [1.441205E+03, -2.708913E+03, 5.601785E+03],
                                        [1.528031E+03, -1.338349E+03, 5.966178E+03],
                                        [-6.552600E+01, 1.321471E+03, 4.213046E+03],
                                        [-4.248900E+01, 2.046976E+03, 5.525725E+03],
                                        [-9.173800E+01, 1.262274E+03, 5.961848E+03],
                                        [-6.552600E+01, -1.321471E+03, 4.213046E+03],
                                        [-4.248900E+01, -2.046976E+03, 5.525725E+03],
                                        [-9.173800E+01, -1.262274E+03, 5.961848E+03]])

    return CAD_matching_pos_values


def get_Kite15_coords():

    # coords = np.array([[ 0.,          0.,          0.        ],
    # [ 0.,          4.71710426,  0.        ],
    # [ 2.20739359,  3.18723261,  1.38374481],
    # [ 2.5094198,   1.74235383,  2.02382997],
    # [ 2.59221754,  0.89242513,  2.26470515],
    # [ 2.6178689,   0.27858774,  2.34917373],
    # [ 2.61805196, -0.26914409,  2.34959129],
    # [ 2.59282071, -0.88298148,  2.26607316],
    # [ 2.51066709, -1.73291017,  2.02664824],
    # [ 2.21027092, -3.17778896,  1.39023904],
    # [ 0.,         -4.71710426,  0.        ],
    # [ 0.,         -3.17778896,  1.38374481],
    # [ 0.,         -1.73291017,  2.02382997],
    # [ 0.,         -0.88298148,  2.26470515],
    # [ 0.,         -0.26914409,  2.34917373],
    # [ 0.,          0.27858774,  2.34959129],
    # [ 0.,          0.89242513,  2.26607316],
    # [ 0.,          1.74235383,  2.02664824],
    # [ 0.,          3.18723261,  1.39023904],
    # [ 1.04824539,  4.71710426,  0.        ],
    # [ 1.04824539, -4.71710426,  0.        ]])

    # Equal angle division
    coords = np.array([[ 0.,          0.,          0.,        ],
 [ 0.,          4.71710426,  0.        ],
 [-2.20739359,  3.18723261,  1.38374481],
 [-2.5094198,   1.74235383,  2.02382997],
 [-2.59221754,  0.89242513,  2.26470515],
 [-2.6178689,   0.27858774,  2.34917373],
 [-2.61805196, -0.26914409,  2.34959129],
 [-2.59282071, -0.88298148,  2.26607316],
 [-2.51066709, -1.73291017,  2.02664824],
 [-2.21027092, -3.17778896,  1.39023904],
 [ 0.,         -4.71710426,  0.        ],
 [ 0.,         -3.17778896,  1.38374481],
 [ 0.,         -1.73291017,  2.02382997],
 [ 0.,         -0.88298148,  2.26470515],
 [ 0.,         -0.26914409,  2.34917373],
 [ 0.,          0.27858774,  2.34959129],
 [ 0.,          0.89242513,  2.26607316],
 [ 0.,          1.74235383,  2.02664824],
 [ 0.,          3.18723261,  1.39023904],
 [-1.04824539,  4.71710426,  0.        ],
 [-1.04824539, -4.71710426,  0.        ]])

    # Equal arc division
    coords = np.array([[ 0.,          0.,          0.        ],
 [ 0.,          4.71710426,  0.        ],
 [-1.88898095,  3.98616936,  1.10483762],
 [-2.2903054,   2.89251708,  1.77829233],
 [-2.50888284,  1.74640111,  2.16089092],
 [-2.6085338,   0.58358657,  2.33713317],
 [-2.60853498, -0.58355827,  2.33713509],
 [-2.5088866,  -1.7463728,   2.16089711],
 [-2.29031272, -2.89248877,  1.77830466],
 [-1.88899586, -3.98614106,  1.10486431],
 [ 0.,         -4.71709483,  0.        ],
 [ 0.,         -3.98614106,  1.10483762],
 [ 0.,         -2.89248877,  1.77829233],
 [ 0.,         -1.7463728,   2.16089092],
 [ 0.,         -0.58355827,  2.33713317],
 [ 0.,          0.58358657,  2.33713509],
 [ 0.,          1.74640111,  2.16089711],
 [ 0.,          2.89251708,  1.77830466],
 [ 0.,          3.98616936,  1.10486431],
 [-1.04824539,  4.71710426,  0.        ],
 [-1.04824539, -4.71709483,  0.0047171 ]])
    return coords

def struct2aero_geometry(coord_struc):
    coord = np.empty((20, 3))

    coord[0, :] = coord_struc[20, :]
    coord[1, :] = coord_struc[10, :]

    coord[2, :] = coord_struc[9, :]
    coord[3, :] = coord_struc[11, :]

    coord[4, :] = coord_struc[8, :]
    coord[5, :] = coord_struc[12, :]

    coord[6, :] = coord_struc[7, :]
    coord[7, :] = coord_struc[13, :]

    coord[8, :] = coord_struc[6, :]
    coord[9, :] = coord_struc[14, :]

    coord[10, :] = coord_struc[5, :]
    coord[11, :] = coord_struc[15, :]

    coord[12, :] = coord_struc[4, :]
    coord[13, :] = coord_struc[16, :]

    coord[14, :] = coord_struc[3, :]
    coord[15, :] = coord_struc[17, :]

    coord[16, :] = coord_struc[2, :]
    coord[17, :] = coord_struc[18, :]

    coord[18, :] = coord_struc[19, :]
    coord[19, :] = coord_struc[1, :]

    return coord


def refine_LEI_mesh(coord, N_sect, N_split):
    refined_coord = []

    for i_sec in range(N_sect):
        temp_coord = np.empty((int(N_split * 2), 3))
        for i_spl in range(N_split):
            temp_coord[2 * i_spl] = coord[2 * i_sec, :] * (N_split - i_spl) / N_split + \
                                    coord[2 * (i_sec + 1), :] * (i_spl) / N_split
            temp_coord[2 * i_spl + 1] = coord[2 * i_sec + 1, :] * (N_split - i_spl) / N_split + \
                                        coord[2 * (i_sec + 1) + 1, :] * (i_spl) / N_split
        if i_sec == 0:
            refined_coord = temp_coord
        else:
            refined_coord = np.append(refined_coord, temp_coord, axis=0)

    refined_coord = np.append(refined_coord, [coord[2 * N_sect, :], coord[2 * N_sect + 1, :]], axis=0)

    return refined_coord


# %% FUNCTIONS NOT CURRENTLY USED

def velocity_3D_from_vortex_filament(XV1, XV2, XVP, GAMMA, CORE):
    # function to calculate the velocity induced by a straight 3D vortex filament
    # with circulation GAMMA at a point VP1. The geometry of the vortex filament
    # is defined by its edges: the filaments start at XV1 and ends at XV2.
    # the input CORE defines a vortex core radius, inside which the velocity
    # is defined  as a solid body rotation.
    # The function is adapted from the algorithm presented in:
    #                Katz, Joseph, and Allen Plotkin. Low-speed aerodynamics.
    #                Vol. 13. Cambridge university press, 2001.

    # read coordinates that define the vortex filament
    X1 = XV1[0];
    Y1 = XV1[1];
    Z1 = XV1[2]  # start point of vortex filament
    X2 = XV2[0];
    Y2 = XV2[1];
    Z2 = XV2[2]  # end point of vortex filament
    # read coordinates of target point where the velocity is calculated
    XP = XVP[0];
    YP = XVP[1];
    ZP = XVP[2]
    # calculate geometric relations for integral of the velocity induced by filament
    R1 = np.sqrt((XP - X1) ** 2 + (YP - Y1) ** 2 + (ZP - Z1) ** 2)
    R2 = np.sqrt((XP - X2) ** 2 + (YP - Y2) ** 2 + (ZP - Z2) ** 2)
    R1XR2_X = (YP - Y1) * (ZP - Z2) - (ZP - Z1) * (YP - Y2)
    R1XR2_Y = -(XP - X1) * (ZP - Z2) + (ZP - Z1) * (XP - X2)
    R1XR2_Z = (XP - X1) * (YP - Y2) - (YP - Y1) * (XP - X2)
    R1XR_SQR = R1XR2_X ** 2 + R1XR2_Y ** 2 + R1XR2_Z ** 2
    R0R1 = (X2 - X1) * (XP - X1) + (Y2 - Y1) * (YP - Y1) + (Z2 - Z1) * (ZP - Z1)
    R0R2 = (X2 - X1) * (XP - X2) + (Y2 - Y1) * (YP - Y2) + (Z2 - Z1) * (ZP - Z2)
    # check if target point is in the vortex filament core,
    # and modify to solid body rotation
    if (R1XR_SQR < CORE ** 2):
        R1XR_SQR = CORE ** 2
    # GAMMA = 0;
    if (R1 < CORE):
        R1 = CORE
    # GAMMA = 0;
    if (R2 < CORE):
        R2 = CORE
    # GAMMA = 0
    # determine scalar
    K = GAMMA / 4 / np.pi / R1XR_SQR * (R0R1 / R1 - R0R2 / R2)
    # determine the three velocity components
    U = K * R1XR2_X
    V = K * R1XR2_Y
    W = K * R1XR2_Z
    # output results, vector with the three velocity components
    results = [U, V, W]
    return results


def velocity_induced_single_ring_semiinfinite_nocore(ring, controlpoint, model):
    CORE = 1e-20
    velind = [0, 0, 0]
    for filament in ring:
        GAMMA = filament['Gamma']
        XV1 = filament['x1']

        XVP = controlpoint
        if filament['id'] == 'trailing_inf1':
            Vf = filament['dir']
            tempvel = velocity_3D_from_vortex_filament_semiinfinite(XV1, Vf, XVP, GAMMA, CORE)
        elif filament['id'] == 'trailing_inf2':
            Vf = filament['dir']
            tempvel = velocity_3D_from_vortex_filament_semiinfinite(XV1, Vf, XVP, -GAMMA, CORE)
        elif filament['id'] == 'bound':
            if model == 'VSM':
                XV2 = filament['x2']
                tempvel = velocity_3D_from_vortex_filament(XV1, XV2, XVP, GAMMA, CORE)
            else:
                tempvel = [0, 0, 0]
        else:
            XV2 = filament['x2']
            tempvel = velocity_3D_from_vortex_filament(XV1, XV2, XVP, GAMMA, CORE)

        velind[0] += tempvel[0]
        velind[1] += tempvel[1]
        velind[2] += tempvel[2]

    return velind


def velocity_3D_from_vortex_filament_semiinfinite(XV1, Vf, XVP, GAMMA, CORE):
    # function to calculate the velocity induced by a straight 3D vortex filament
    # with circulation GAMMA at a point VP1. The geometry of the vortex filament
    # is defined by its edges: the filaments start at XV1 and ends at XV2.
    # the input CORE defines a vortex core radius, inside which the velocity
    # is defined  as a solid body rotation.
    # The function is adapted from the algorithm presented in:
    #                Katz, Joseph, and Allen Plotkin. Low-speed aerodynamics.
    #                Vol. 13. Cambridge university press, 2001.

    # read coordinates that define the vortex filament
    X1 = XV1[0];
    Y1 = XV1[1];
    Z1 = XV1[2]  # start point of vortex filament
    Vfx = Vf[0];
    Vfy = Vf[1];
    Vfz = Vf[2]  # end point of vortex filament
    # read coordinates of target point where the velocity is calculated
    XP = XVP[0];
    YP = XVP[1];
    ZP = XVP[2]
    # calculate geometric relations for integral of the velocity induced by filament
    R1 = np.sqrt((XP - X1) ** 2 + (YP - Y1) ** 2 + (ZP - Z1) ** 2)

    R1XV_X = (YP - Y1) * Vfz - (ZP - Z1) * Vfy
    R1XV_Y = -(XP - X1) * Vfz + (ZP - Z1) * Vfx
    R1XV_Z = (XP - X1) * Vfy - (YP - Y1) * Vfx
    R1XR_SQR = R1XV_X ** 2 + R1XV_Y ** 2 + R1XV_Z ** 2
    VfR1 = Vfx * (XP - X1) + Vfy * (YP - Y1) + Vfz * (ZP - Z1)
    # check if target point is in the vortex filament core,
    # and modify to solid body rotation
    if (R1XR_SQR < CORE ** 2):
        R1XR_SQR = CORE ** 2
    # GAMMA = 0;
    if (R1 < CORE):
        R1 = CORE

    # GAMMA = 0
    # determine scalar
    K = GAMMA / 4 / np.pi / R1XR_SQR * (1 + VfR1 / R1)
    # determine the three velocity components
    U = K * R1XV_X
    V = K * R1XV_Y
    W = K * R1XV_Z
    # output results, vector with the three velocity components
    results = [U, V, W]
    return results


def output_results_LEI(Fmag, aero_coeffs, ringvec, Uinf, controlpoints, Atot):
    N_struct = 9  # Number of panels in the deformation model
    alpha = aero_coeffs[:, 0]
    F_rel = []
    F_gl = []
    Fmag_gl = []
    F_panel = np.empty((N_struct, 3))
    M_panel = np.empty((N_struct))
    Ltot = 0
    Dtot = 0
    N_split = int(len(alpha) / N_struct)
    for i in range(len(alpha)):
        # check directions
        # twist = np.arcsin(controlpoints[i]['normal'][0])
        # twist = 0
        # dir_urel = np.array([np.cos(alpha[i]-twist),0,np.sin(alpha[i]-twist)])
        r0 = ringvec[i]['r0']

        dir_urel = np.cos(alpha[i]) * controlpoints[i]['tangential'] + np.sin(alpha[i]) * controlpoints[i]['normal']
        dir_urel = dir_urel / np.linalg.norm(dir_urel)

        dir_L = np.cross(dir_urel, r0)
        dir_L = dir_L / np.linalg.norm(dir_L)
        dir_D = np.cross([0, 1, 0], dir_L)
        dir_D = dir_D / np.linalg.norm(dir_D)
        L_rel = dir_L * Fmag[i, 0]
        D_rel = dir_D * Fmag[i, 1]
        F_rel.append([L_rel, D_rel])

        dir_L_gl = np.cross(Uinf, [0, 1, 0])
        dir_L_gl = dir_L_gl / vec_norm(dir_L_gl)
        L_gl = vector_projection(L_rel, dir_L_gl) + vector_projection(D_rel, dir_L_gl)
        D_gl = vector_projection(L_rel, Uinf) + vector_projection(D_rel, Uinf)
        F_gl.append([L_gl, D_gl])
        Fmag_gl.append([dot_product(L_rel, dir_L_gl) + dot_product(D_rel, dir_L_gl) \
                           , dot_product(L_rel, Uinf / vec_norm(Uinf)) + dot_product(D_rel, Uinf / vec_norm(Uinf))])

    for i in range(len(Fmag_gl)):
        Ltot += Fmag_gl[i][0] * np.linalg.norm(ringvec[i]['r0'])
        Dtot += Fmag_gl[i][1] * np.linalg.norm(ringvec[i]['r0'])
        sec = (N_struct - 1) - int((i + 1) / N_split - 0.01)
        F_panel[sec] += (F_rel[i][0] + F_rel[i][1]) * np.linalg.norm(ringvec[i]['r0'])
        M_panel[sec] += Fmag[i, 2] * np.linalg.norm(ringvec[i]['r0'])

    for i in range(len(F_panel)):
        if abs(F_panel[i, 1]) < 1e-10:
            F_panel[i, 1] = 0

    Umag = np.linalg.norm(Uinf)
    CL = Ltot / (0.5 * Umag ** 2 * Atot)
    CD = Dtot / (0.5 * Umag ** 2 * Atot)

    return F_rel, F_gl, Ltot, Dtot, CL, CD


def create_geometry_general_frozen_wake(coordinates, Uinf, N, ring_geo, model, vortex_N):
    filaments = []
    controlpoints = []
    rings = []
    wingpanels = []
    ringvec = []
    coord_L = []
    for i in range(N - 1):

        section = {

           'p1': coordinates[2 * i, :],
            'p2': coordinates[2 * i + 2, :],
            'p3': coordinates[2 * i + 3, :],
            'p4': coordinates[2 * i + 1, :]
        }
        wingpanels.append(section)

        chord = np.linalg.norm((section['p2'] + section['p1']) / 2 - (section['p3'] + section['p4']) / 2)
        LLpoint = (section['p2'] + section['p1']) / 2 * 3 / 4 + (section['p3'] + section['p4']) / 2 * 1 / 4
        VSMpoint = (section['p2'] + section['p1']) / 2 * 1 / 4 + (section['p3'] + section['p4']) / 2 * 3 / 4
        coord_L.append(LLpoint)

        # Define bound vortex filament
        bound = {'id': 'bound',
                 'x1': section['p1'] * 3 / 4 + section['p4'] * 1 / 4,
                 'x2': section['p2'] * 3 / 4 + section['p3'] * 1 / 4,
                 'Gamma': 0}
        filaments.append(bound)

        x_airf = np.cross(VSMpoint - LLpoint, section['p2'] - section['p1'])
        x_airf = x_airf / np.linalg.norm(x_airf)
        y_airf = VSMpoint - LLpoint
        y_airf = y_airf / np.linalg.norm(y_airf)
        z_airf = bound['x2'] - bound['x1']
        # z_airf[0] = 0
        z_airf = z_airf / np.linalg.norm(z_airf)
        airf_coord = np.column_stack([x_airf, y_airf, z_airf])

        normal = x_airf
        tangential = y_airf
        if model == 'VSM':
            cp = {'coordinates': VSMpoint, 'chord': chord,
                  'normal': normal,
                  'tangential': tangential,
                  'airf_coord': airf_coord,
                  'coordinates_aoa': LLpoint}
            controlpoints.append(cp)
        elif model == 'LLT':

            cp = {'coordinates': LLpoint, 'chord': chord,
                  'normal': normal,
                  'tangential': tangential,
                  'airf_coord': airf_coord}
            controlpoints.append(cp)

        temp = {'r0': bound['x2'] - bound['x1'],
                'r1': cp['coordinates'] - bound['x1'],
                'r2': cp['coordinates'] - bound['x2'],
                'r3': cp['coordinates'] - (bound['x2'] + bound['x1']) / 2}
        ringvec.append(temp)

        temp = Uinf / np.linalg.norm(Uinf)

        dt = 1
        dx = Uinf[0] * dt
        dy = Uinf[1] * dt
        dz = Uinf[2] * dt
        if ring_geo == '3fil':

            temp1 = {'x1': np.array([bound['x1'][0] + dx, bound['x1'][1] + dy, bound['x1'][2] + dz]),
                     'x2': bound['x1'], 'Gamma': 0,
                     'id': 'trailing1'}
            filaments.append(temp1)

            for j in range(vortex_N):
                xt = filaments[len(filaments) - 1]['x1'][0];
                yt = filaments[len(filaments) - 1]['x1'][1];
                zt = filaments[len(filaments) - 1]['x1'][2];

                temp1 = {'x1': np.array([xt + dx, yt + dy, zt + dz]), 'x2': np.array([xt, yt, zt]), 'Gamma': 0,
                         'id': 'trailing2'}
                filaments.append(temp1)

            temp1 = {'x2': np.array([bound['x2'][0] + dx, bound['x2'][1] + dy, bound['x2'][2] + dz]),
                     'x1': bound['x2'], 'Gamma': 0,
                     'id': 'trailing2'}
            filaments.append(temp1)

            for j in range(vortex_N):
                xt = filaments[len(filaments) - 1]['x2'][0];
                yt = filaments[len(filaments) - 1]['x2'][1];
                zt = filaments[len(filaments) - 1]['x2'][2];

                temp1 = {'x1': np.array([xt, yt, zt]), 'x2': np.array([xt + dx, yt + dy, zt + dz]), 'Gamma': 0,
                         'id': 'trailing2'}
                filaments.append(temp1)


        elif ring_geo == '5fil':
            temp1 = {'x1': section['p4'],
                     'x2': bound['x1'], 'Gamma': 0,
                     'id': 'trailing1'}
            filaments.append(temp1)

            for j in range(vortex_N):
                xt = filaments[len(filaments) - 1]['x1'][0];
                yt = filaments[len(filaments) - 1]['x1'][1];
                zt = filaments[len(filaments) - 1]['x1'][2];

                temp1 = {'x1': np.array([xt + dx, yt + dy, zt + dz]), 'x2': np.array([xt, yt, zt]), 'Gamma': 0,
                         'id': 'trailing2'}
                filaments.append(temp1)

                # create trailing filaments, at x2 of bound filament
            temp1 = {'x2': section['p3'],
                     'x1': bound['x2'], 'Gamma': 0,
                     'id': 'trailing2'}
            filaments.append(temp1)

            for j in range(vortex_N):
                xt = filaments[len(filaments) - 1]['x2'][0];
                yt = filaments[len(filaments) - 1]['x2'][1];
                zt = filaments[len(filaments) - 1]['x2'][2];

                temp1 = {'x1': np.array([xt, yt, zt]), 'x2': np.array([xt + dx, yt + dy, zt + dz]), 'Gamma': 0,
                         'id': 'trailing2'}
                filaments.append(temp1)

        #

        #

        rings.append(filaments)
        filaments = []

    coord_L = np.array(coord_L)
    return controlpoints, rings, wingpanels, ringvec, coord_L