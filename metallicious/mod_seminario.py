#   Program to implement the Modified Seminario Method
#   Written by Alice E. A. Allen, TCM, University of Cambridge
#   Adapted for metallicious by TK Piskorz, 2023/06/28
#   Reference using AEA Allen, MC Payne, DJ Cole, J. Chem. Theory Comput. (2018), doi:10.1021/acs.jctc.7b00785

import numpy as np
from operator import itemgetter

def modified_seminario_method(hessian, coords, atom_names, bond_list, angle_list, vibrational_scaling = 1):
    N = len(coords)
    vibrational_scaling_squared = vibrational_scaling**2

    bond_lengths = np.zeros((N, N))

    for i in range (0,N):
        for j in range(0,N):
            diff_i_j = np.array(coords[i,:]) - np.array(coords[j,:])
            bond_lengths[i][j] =  np.linalg.norm(diff_i_j)

    eigenvectors = np.empty((3, 3, N, N), dtype=complex)
    eigenvalues = np.empty((N, N, 3), dtype=complex)

    for i in range(0,N):
        for j in range(0,N):
            partial_hessian = hessian[(i * 3):((i + 1)*3),(j * 3):((j + 1)*3)]
            [a, b] = np.linalg.eig(partial_hessian)
            eigenvalues[i,j,:] = (a)
            eigenvectors[:,:,i,j] = (b)

    # The bond values are calculated and written to file
    bonds = bonds_calculated_printed(vibrational_scaling_squared, bond_list, bond_lengths, atom_names, eigenvalues, eigenvectors, coords )

    # The angle values are calculated and written to file
    angles = angles_calculated_printed(vibrational_scaling_squared, angle_list, bond_lengths, atom_names, eigenvalues, eigenvectors, coords )

    return bonds, angles


def bonds_calculated_printed(vibrational_scaling_squared, bond_list, bond_lengths, atom_names, eigenvalues, eigenvectors, coords):
    bonds = {}

    #This function uses the Seminario method to find the bond
    #parameters and print them to file

    #Open output file bond parameters are written to
    fid = open('Modified_Seminario_Bonds', 'w')

    k_b = np.zeros(len(bond_list))
    bond_length_list = np.zeros(len(bond_list))
    unique_values_bonds = [] # Used to find average values

    for i in range(0, len(bond_list)):
        AB = force_constant_bond(bond_list[i][0], bond_list[i][1],eigenvalues, eigenvectors, coords)
        BA = force_constant_bond(bond_list[i][1], bond_list[i][0],eigenvalues, eigenvectors, coords)

        # Order of bonds sometimes causes slight differences, find the mean
        k_b[i] = np.real(( AB + BA ) /2);

        # Vibrational_scaling takes into account DFT deficities/ anharmocity
        k_b[i] = k_b[i] * vibrational_scaling_squared

        bond_length_list[i] =  bond_lengths[bond_list[i][0]][bond_list[i][1]]
        fid.write(atom_names[bond_list[i][0]] + '-' + atom_names[bond_list[i][1]] + '  ')
        fid.write(str("%#.5g" % k_b[i])+ '   ' + str("%#.4g" % bond_length_list[i]) +  '   ' +
                  str(bond_list[i][0] + 1) +  '   ' + str(bond_list[i][1] + 1))
        fid.write('\n')

        #unique_values_bonds.append([atom_names[bond_list[i][0]], atom_names[bond_list[i][1]], , , 1 ])

        order = np.argsort(bond_list[i])
        sorted_bond_list = (bond_list[i][order[0]], bond_list[i][order[1]])
        sorted_atom_name = (atom_names[bond_list[i][order[0]]], atom_names[bond_list[i][order[1]]])

        bonds[(sorted_bond_list, sorted_atom_name)]=(bond_length_list[i], k_b[i])

    fid.close()

    return bonds

def force_constant_bond(atom_A, atom_B, eigenvalues, eigenvectors, coords):
    #Force Constant - Equation 10 of Seminario paper - gives force
    #constant for bond


    #Eigenvalues and eigenvectors calculated
    eigenvalues_AB = eigenvalues[atom_A,atom_B,:]

    eigenvectors_AB = eigenvectors[:,:,atom_A,atom_B]

    #Vector along bond
    diff_AB = np.array(coords[atom_B,:]) - np.array(coords[atom_A,:])
    norm_diff_AB = np.linalg.norm(diff_AB)

    unit_vectors_AB = diff_AB / norm_diff_AB

    k_AB = 0

    #Projections of eigenvalues
    for i in range(0,3):
        dot_product = abs(np.dot(unit_vectors_AB, eigenvectors_AB[:,i]))
        k_AB = k_AB + ( eigenvalues_AB[i] * dot_product )

    k_AB = -k_AB * 0.5 # Convert to OPLS form

    return k_AB



def angles_calculated_printed(vibrational_scaling_squared, angle_list, bond_lengths, atom_names, eigenvalues,
                              eigenvectors, coords):
    # This function uses the modified Seminario method to find the angle
    # parameters and print them to file

    # Open output file angle parameters are written to
    fid = open('Modified_Seminario_Angle', 'w')

    k_theta = np.zeros(len(angle_list))
    theta_0 = np.zeros(len(angle_list))
    unique_values_angles = []  # Used to find average values

    # Modified Seminario part goes here ...
    # Connectivity information for Modified Seminario Method
    central_atoms_angles = []

    angles = {}

    # A structure is created with the index giving the central atom of the
    # angle, an array then lists the angles with that central atom.
    # ie. central_atoms_angles{3} contains an array of angles with central atom
    # 3
    for i in range(0, len(coords)):
        central_atoms_angles.append([])
        for j in range(0, len(angle_list)):
            if i == angle_list[j][1]:
                # For angle ABC, atoms A C are written to array
                AC_array = [angle_list[j][0], angle_list[j][2], j]
                central_atoms_angles[i].append(AC_array)

                # For angle ABC, atoms C A are written to array
                CA_array = [angle_list[j][2], angle_list[j][0], j]
                central_atoms_angles[i].append(CA_array)

    # Sort rows by atom number
    for i in range(0, len(coords)):
        central_atoms_angles[i] = sorted(central_atoms_angles[i], key=itemgetter(0))

    # Find normals u_PA for each angle
    unit_PA_all_angles = []

    for i in range(0, len(central_atoms_angles)):
        unit_PA_all_angles.append([])
        for j in range(0, len(central_atoms_angles[i])):
            # For the angle at central_atoms_angles[i][j,:] the corresponding
            # u_PA value is found for the plane ABC and bond AB, where ABC
            # corresponds to the order of the arguements
            # This is why the reverse order was also added
            unit_PA_all_angles[i].append(
                u_PA_from_angles(central_atoms_angles[i][j][0], i, central_atoms_angles[i][j][1], coords))

    # Finds the contributing factors from the other angle terms
    # scaling_factor_all_angles = cell(max(max(angle_list))); %This array will contain scaling factor and angle list position
    scaling_factor_all_angles = []

    for i in range(0, len(central_atoms_angles)):
        scaling_factor_all_angles.append([])
        for j in range(0, len(central_atoms_angles[i])):
            n = 1
            m = 1
            angles_around = 0
            additional_contributions = 0
            scaling_factor_all_angles[i].append([0, 0])

            # Position in angle list
            scaling_factor_all_angles[i][j][1] = central_atoms_angles[i][j][2]

            # Goes through the list of angles with the same central atom
            # And computes the term need for the modified Seminario method

            # Forwards directions, finds the same bonds with the central atom i
            while (((j + n) < len(central_atoms_angles[i])) and central_atoms_angles[i][j][0] ==
                   central_atoms_angles[i][j + n][0]):
                additional_contributions = additional_contributions + (
                    abs(np.dot(unit_PA_all_angles[i][j][:], unit_PA_all_angles[i][j + n][:]))) ** 2
                n = n + 1
                angles_around = angles_around + 1

            # Backwards direction, finds the same bonds with the central atom i
            while (((j - m) >= 0) and central_atoms_angles[i][j][0] == central_atoms_angles[i][j - m][0]):
                additional_contributions = additional_contributions + (
                    abs(np.dot(unit_PA_all_angles[i][j][:], unit_PA_all_angles[i][j - m][:]))) ** 2
                m = m + 1
                angles_around = angles_around + 1

            if (n != 1 or m != 1):
                # Finds the mean value of the additional contribution
                # To change to normal Seminario method comment out + part
                scaling_factor_all_angles[i][j][0] = 1 + (additional_contributions / (m + n - 2))
            else:
                scaling_factor_all_angles[i][j][0] = 1

    scaling_factors_angles_list = []
    for i in range(0, len(angle_list)):
        scaling_factors_angles_list.append([])

        # Orders the scaling factors according to the angle list
    for i in range(0, len(central_atoms_angles)):
        for j in range(0, len(central_atoms_angles[i])):
            scaling_factors_angles_list[scaling_factor_all_angles[i][j][1]].append(scaling_factor_all_angles[i][j][0])

            # Finds the angle force constants with the scaling factors included for each angle
    for i in range(0, len(angle_list)):
        # Ensures that there is no difference when the ordering is changed
        [AB_k_theta, AB_theta_0] = force_angle_constant(angle_list[i][0], angle_list[i][1], angle_list[i][2],
                                                        bond_lengths, eigenvalues, eigenvectors, coords,
                                                        scaling_factors_angles_list[i][0],
                                                        scaling_factors_angles_list[i][1])
        [BA_k_theta, BA_theta_0] = force_angle_constant(angle_list[i][2], angle_list[i][1], angle_list[i][0],
                                                        bond_lengths, eigenvalues, eigenvectors, coords,
                                                        scaling_factors_angles_list[i][1],
                                                        scaling_factors_angles_list[i][0])
        k_theta[i] = (AB_k_theta + BA_k_theta) / 2
        theta_0[i] = (AB_theta_0 + BA_theta_0) / 2

        # Vibrational_scaling takes into account DFT deficities/ anharmocity
        k_theta[i] = k_theta[i] * vibrational_scaling_squared

        fid.write(atom_names[angle_list[i][0]] + '-' + atom_names[angle_list[i][1]] + '-' + atom_names[
            angle_list[i][2]] + '  ')
        fid.write(str("%#.4g" % k_theta[i]) + '   ' + str("%#.4g" % theta_0[i]) + '   ' + str(
            angle_list[i][0] + 1) + '   ' + str(angle_list[i][1] + 1) + '   ' + str(angle_list[i][2] + 1))
        fid.write('\n')

        indecies = angle_list[i]

        order = [0, 1, 2]
        if indecies[0] > indecies[2]:  # we have indecies from smaler to larger in the angles
            order = [2, 1, 0]

        sorted_angle_list = (angle_list[i][order[0]], angle_list[i][order[1]], angle_list[i][order[2]])
        sorted_atom_name = (
        atom_names[angle_list[i][order[0]]], atom_names[angle_list[i][order[1]]], atom_names[angle_list[i][order[2]]])

        angles[(sorted_angle_list, sorted_atom_name)] = (theta_0[i], k_theta[i])

    fid.close()

    return angles



def u_PA_from_angles(atom_A, atom_B, atom_C, coords):
    # This gives the vector in the plane A,B,C and perpendicular to A to B

    diff_AB = coords[atom_B, :] - coords[atom_A, :]
    norm_diff_AB = np.linalg.norm(diff_AB)
    u_AB = diff_AB / norm_diff_AB

    diff_CB = coords[atom_B, :] - coords[atom_C, :]
    norm_diff_CB = np.linalg.norm(diff_CB)
    u_CB = diff_CB / norm_diff_CB

    u_N = unit_vector_N(u_CB, u_AB)

    u_PA = np.cross(u_N, u_AB)
    norm_PA = np.linalg.norm(u_PA)
    u_PA = u_PA / norm_PA;

    return u_PA


def force_angle_constant(atom_A, atom_B, atom_C, bond_lengths, eigenvalues, eigenvectors, coords, scaling_1, scaling_2):
    # Force Constant- Equation 14 of seminario calculation paper - gives force
    # constant for angle (in kcal/mol/rad^2) and equilibrium angle in degrees

    # Vectors along bonds calculated
    diff_AB = coords[atom_B, :] - coords[atom_A, :]
    norm_diff_AB = np.linalg.norm(diff_AB)
    u_AB = diff_AB / norm_diff_AB

    diff_CB = coords[atom_B, :] - coords[atom_C, :]
    norm_diff_CB = np.linalg.norm(diff_CB)
    u_CB = diff_CB / norm_diff_CB

    # Bond lengths and eigenvalues found
    bond_length_AB = bond_lengths[atom_A, atom_B]
    eigenvalues_AB = eigenvalues[atom_A, atom_B, :]
    eigenvectors_AB = eigenvectors[0:3, 0:3, atom_A, atom_B]

    bond_length_BC = bond_lengths[atom_B, atom_C]
    eigenvalues_CB = eigenvalues[atom_C, atom_B, :]
    eigenvectors_CB = eigenvectors[0:3, 0:3, atom_C, atom_B]

    # Normal vector to angle plane found
    u_N = unit_vector_N(u_CB, u_AB)

    u_PA = np.cross(u_N, u_AB)
    norm_u_PA = np.linalg.norm(u_PA)
    u_PA = u_PA / norm_u_PA

    u_PC = np.cross(u_CB, u_N)
    norm_u_PC = np.linalg.norm(u_PC)
    u_PC = u_PC / norm_u_PC

    sum_first = 0
    sum_second = 0

    # Projections of eigenvalues
    for i in range(0, 3):
        eig_AB_i = eigenvectors_AB[:, i]
        eig_BC_i = eigenvectors_CB[:, i]
        sum_first = sum_first + (eigenvalues_AB[i] * abs(dot_product(u_PA, eig_AB_i)))
        sum_second = sum_second + (eigenvalues_CB[i] * abs(dot_product(u_PC, eig_BC_i)))

        # Scaling due to additional angles - Modified Seminario Part
    sum_first = sum_first / scaling_1
    sum_second = sum_second / scaling_2

    # Added as two springs in series
    k_theta = (1 / ((bond_length_AB ** 2) * sum_first)) + (1 / ((bond_length_BC ** 2) * sum_second))
    k_theta = 1 / k_theta

    k_theta = - k_theta  # Change to OPLS form
    k_theta = abs(k_theta * 0.5)  # Change to OPLS form

    # Equilibrium Angle
    theta_0 = np.degrees(np.arccos(np.dot(u_AB, u_CB)))

    # If the vectors u_CB and u_AB are linearly dependent u_N cannot be defined
    # This case is dealt with here
    if abs(sum((u_CB) - (u_AB))) < 0.01 or (abs(sum((u_CB) - (u_AB))) > 1.99 and abs(sum((u_CB) - (u_AB))) < 2.01):
        scaling_1 = 1
        scaling_2 = 1
        [k_theta, theta_0] = force_angle_constant_special_case(atom_A, atom_B, atom_C, bond_lengths, eigenvalues,
                                                               eigenvectors, coords, scaling_1, scaling_2)

    return k_theta, theta_0


def dot_product(u_PA, eig_AB):
    x = 0
    for i in range(0, 3):
        x = x + u_PA[i] * eig_AB[i].conjugate()
    return x


def force_angle_constant_special_case(atom_A, atom_B, atom_C, bond_lengths, eigenvalues, eigenvectors, coords,
                                      scaling_1, scaling_2):
    # Force Constant- Equation 14 of seminario calculation paper - gives force
    # constant for angle (in kcal/mol/rad^2) and equilibrium angle in degrees
    # This function deals with cases when u_N cannot be defined and instead
    # takes samples of u_N across a unit sphere.

    # Vectors along bonds calculated
    diff_AB = coords[atom_B, :] - coords[atom_A, :]
    norm_diff_AB = np.linalg.norm(diff_AB)
    u_AB = diff_AB / norm_diff_AB

    diff_CB = coords[atom_B, :] - coords[atom_C, :]
    norm_diff_CB = np.linalg.norm(diff_CB)
    u_CB = diff_CB / norm_diff_CB

    # Bond lengths and eigenvalues found
    bond_length_AB = bond_lengths[atom_A, atom_B]
    eigenvalues_AB = eigenvalues[atom_A, atom_B, :]
    eigenvectors_AB = eigenvectors[0:3, 0:3, atom_A, atom_B]

    bond_length_BC = bond_lengths[atom_B, atom_C]
    eigenvalues_CB = eigenvalues[atom_C, atom_B, :]
    eigenvectors_CB = eigenvectors[0:3, 0:3, atom_C, atom_B]

    k_theta_array = np.zeros((180, 360))

    # Find force constant with varying u_N (with vector uniformly sampled across a sphere)
    for theta in range(0, 180):
        for phi in range(0, 360):
            r = 1
            u_N = [r * np.sin(np.deg2rad(theta)) * np.cos(np.deg2rad(theta)),
                   r * np.sin(np.deg2rad(theta)) * np.sin(np.deg2rad(theta)), r * np.cos(np.deg2rad(theta))]

            u_PA = np.cross(u_N, u_AB)
            u_PA = u_PA / np.linalg.norm(u_PA)

            u_PC = np.cross(u_CB, u_N)
            u_PC = u_PC / np.linalg.norm(u_PC)

            sum_first = 0
            sum_second = 0

            # Projections of eigenvalues
            for i in range(0, 3):
                eig_AB_i = eigenvectors_AB[:, i]
                eig_BC_i = eigenvectors_CB[:, i]
                sum_first = sum_first + (eigenvalues_AB[i] * abs(dot_product(u_PA, eig_AB_i)))
                sum_second = sum_second + (eigenvalues_CB[i] * abs(dot_product(u_PC, eig_BC_i)))

                # Added as two springs in series
            k_theta_ij = (1 / ((bond_length_AB ** 2) * sum_first)) + (1 / ((bond_length_BC ** 2) * sum_second))
            k_theta_ij = 1 / k_theta_ij

            k_theta_ij = - k_theta_ij  # Change to OPLS form

            k_theta_ij = abs(k_theta_ij * 0.5)  # Change to OPLS form
            k_theta_array[theta, phi] = k_theta_ij

    # Removes cases where u_N was linearly dependent of u_CB or u_AB

    # Force constant used is taken as the mean
    k_theta = np.mean(np.mean(k_theta_array))
    # Equilibrium Angle independent of u_N
    theta_0 = np.degrees(np.arccos(np.dot(u_AB, u_CB)))

    return k_theta, theta_0


def unit_vector_N(u_BC, u_AB):
    # Calculates unit normal vector which is perpendicular to plane ABC

    cross_product = np.cross(u_BC, u_AB)
    norm_u_N = np.linalg.norm(cross_product)
    u_N = cross_product / norm_u_N
    return u_N
