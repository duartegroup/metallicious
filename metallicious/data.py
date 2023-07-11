# Taken form: https://gist.github.com/lukasrichters14/
# Dictionary of all elements matched with their atomic masses.

name2mass = {'H' : 1.008,'He' : 4.003, 'Li' : 6.941, 'BE' : 9.012,\
                 'B' : 10.811, 'C' : 12.011, 'N' : 14.007, 'O' : 15.999,\
                 'F' : 18.998, 'Ne' : 20.180, 'Na' : 22.990, 'Mg' : 24.305,\
                 'Al' : 26.982, 'Si' : 28.086, 'P' : 30.974, 'S' : 32.066,\
                 'Cl' : 35.453, 'Ar' : 39.948, 'K' : 39.098, 'Ca' : 40.078,\
                 'Sc' : 44.956, 'Ti' : 47.867, 'V' : 50.942, 'Cr' : 51.996,\
                 'Mn' : 54.938, 'Fe' : 55.845, 'Co' : 58.933, 'Ni' : 58.693,\
                 'Cu' : 63.546, 'Zn' : 65.38, 'Ga' : 69.723, 'Ge' : 72.631,\
                 'As' : 74.922, 'Se' : 78.971, 'Br' : 79.904, 'Kr' : 84.798,\
                 'Rb' : 84.468, 'Sr' : 87.62, 'Y' : 88.906, 'Zr' : 91.224,\
                 'Nb' : 92.906, 'Mo' : 95.95, 'Tc' : 98.907, 'Ru' : 101.07,\
                 'Rh' : 102.906, 'Pd' : 106.42, 'Ag' : 107.868, 'Cd' : 112.414,\
                 'In' : 114.818, 'Sn' : 118.711, 'Sb' : 121.760, 'Te' : 126.7,\
                 'I' : 126.904, 'Xe' : 131.294, 'Cs' : 132.905, 'Ba' : 137.328,\
                 'La' : 138.905, 'Ce' : 140.116, 'Pr' : 140.908, 'Nd' : 144.243,\
                 'Pm' : 144.913, 'Sm' : 150.36, 'Eu' : 151.964, 'Gd' : 157.25,\
                 'Tb' : 158.925, 'Dy': 162.500, 'Ho' : 164.930, 'Er' : 167.259,\
                 'Tm' : 168.934, 'Yb' : 173.055, 'Lu' : 174.967, 'Hf' : 178.49,\
                 'Ta' : 180.948, 'W' : 183.84, 'Re' : 186.207, 'Os' : 190.23,\
                 'Ir' : 192.217, 'Pt' : 195.085, 'Au' : 196.967, 'Hg' : 200.592,\
                 'Tl' : 204.383, 'Pb' : 207.2, 'Bi' : 208.980, 'Po' : 208.982,\
                 'At' : 209.987, 'Rn' : 222.081, 'Fr' : 223.020, 'Ra' : 226.025,\
                 'Ac' : 227.028, 'Th' : 232.038, 'Pa' : 231.036, 'U' : 238.029,\
                 'Np' : 237, 'Pu' : 244, 'Am' : 243, 'Cm' : 247, 'Bk' : 247,\
                 'Ct' : 251, 'Es' : 252, 'Fm' : 257, 'Md' : 258, 'No' : 259,\
                 'Lr' : 262, 'Rf' : 261, 'Db' : 262, 'Sg' : 266, 'Bh' : 264,\
                 'Hs' : 269, 'Mt' : 268, 'Ds' : 271, 'Rg' : 272, 'Cn' : 285,\
                 'Nh' : 284, 'Fl' : 289, 'Mc' : 288, 'Lv' : 292, 'Ts' : 294,\
                 'Og' : 294}


name_to_atomic_number = {'H': 1, 'He': 2, 'Li': 3, 'Be': 4, 'B': 5, 'C': 6, 'N': 7, 'O': 8, 'F': 9, 'Ne': 10, 'Na': 11,
                         'Mg': 12, 'Al': 13, 'Si': 14, 'P': 15, 'S': 16, 'Cl': 17, 'Ar': 18, 'K': 19, 'Ca': 20,
                         'Sc': 21, 'Ti': 22, 'V': 23, 'Cr': 24, 'Mn': 25, 'Fe': 26, 'Co': 27, 'Ni': 28, 'Cu': 29,
                         'Zn': 30, 'Ga': 31, 'Ge': 32, 'As': 33, 'Se': 34, 'Br': 35, 'Kr': 36, 'Rb': 37, 'Sr': 38,
                         'Y': 39, 'Zr': 40, 'Nb': 41, 'Mo': 42, 'Tc': 43, 'Ru': 44, 'Rh': 45, 'Pd': 46, 'Ag': 47,
                         'Cd': 48, 'In': 49, 'Sn': 50, 'Sb': 51, 'Te': 52, 'I': 53, 'Xe': 54, 'Cs': 55, 'Ba': 56,
                         'La': 57, 'Ce': 58, 'Pr': 59, 'Nd': 60, 'Pm': 61, 'Sm': 62, 'Eu': 63, 'Gd': 64, 'Tb': 65,
                         'Dy': 66, 'Ho': 67, 'Er': 68, 'Tm': 69, 'Yb': 70, 'Lu': 71, 'Hf': 72, 'Ta': 73, 'W': 74,
                         'Re': 75, 'Os': 76, 'Ir': 77, 'Pt': 78, 'Au': 79, 'Hg': 80, 'Tl': 81, 'Pb': 82, 'Bi': 83,
                         'Po': 84, 'At': 85, 'Rn': 86, 'Fr': 87, 'Ra': 88, 'Ac': 89, 'Th': 90, 'Pa': 91, 'U': 92,
                         'Np': 93, 'Pu': 94, 'Am': 95, 'Cm': 96, 'Bk': 97, 'Cf': 98, 'Es': 99, 'Fm': 100, 'Md': 101,
                         'No': 102, 'Lr': 103, 'Rf': 104, 'Db': 105, 'Sg': 106, 'Bh': 107, 'Hs': 108, 'Mt': 109,
                         'Ds': 110, 'Rg': 111, 'Uub': 112}

vdw_data = {}

vdw_data['merz-tip3p'] = {'Li1': (0.00312065, 1.267), 'Na1': (0.03171494, 1.475), 'K1': (0.15131351, 1.719), 'Rb1': (0.24140216, 1.834), 'Cs1': (0.37853483, 1.988), 'Tl1': (0.14021803, 1.703), 'Cu1': (0.001123, 1.201), 'Ag1': (0.00818431, 1.341), 'Be2': (3.95e-06, 0.956), 'Cu2': (0.00148497, 1.218), 'Ni2': (0.0026232, 1.255), 'Pt2': (0.00307642, 1.266), 'Zn2': (0.00330286, 1.271), 'Co2': (0.00483892, 1.299), 'Pd2': (0.00509941, 1.303), 'Ag2': (0.00770969, 1.336), 'Cr2': (0.00868178, 1.346), 'Fe2': (0.00941798, 1.353), 'Mg2': (0.01020237, 1.36), 'V2': (0.01067299, 1.364), 'Mn2': (0.0168671, 1.407), 'Hg2': (0.0168671, 1.407), 'Cd2': (0.01773416, 1.412), 'Yb2': (0.10185975, 1.642), 'Ca2': (0.1059287, 1.649), 'Sn2': (0.11617738, 1.666), 'Pb2': (0.17018074, 1.745), 'Eu2': (0.21475916, 1.802), 'Sr2': (0.22132374, 1.81), 'Sm2': (0.22878796, 1.819), 'Ba2': (0.40664608, 2.019), 'Ra2': (0.40664608, 2.019), 'Al3': (8.32e-06, 0.981), 'Fe3': (0.00011017, 1.082), 'Cr3': (0.00089969, 1.188), 'In3': (0.00114198, 1.202), 'Tl3': (0.00122067, 1.206), 'Y3': (0.02639002, 1.454), 'La3': (0.09399072, 1.628), 'Ce3': (0.07688443, 1.595), 'Pr3': (0.06441235, 1.568), 'Nd3': (0.05605698, 1.548), 'Sm3': (0.04630154, 1.522), 'Eu3': (0.03994409, 1.503), 'Gd3': (0.03745682, 1.495), 'Tb3': (0.03336723, 1.481), 'Dy3': (0.02986171, 1.468), 'Er3': (0.02133669, 1.431), 'Tm3': (0.01937874, 1.421), 'Lu3': (0.01937874, 1.421), 'Hf4': (0.00012321, 1.087), 'Zr4': (0.00036479, 1.139), 'Ce4': (0.00941798, 1.353), 'U4': (0.00128267, 1.209), 'Pu4': (0.0033972, 1.273), 'Th4': (0.0285863, 1.463)}

vdw_data['merz-opc3'] = {'Li1': (0.0032565, 1.27), 'Na1': (0.0301216, 1.469), 'K1': (0.14021803, 1.703), 'Rb1': (0.21312875, 1.8), 'Cs1': (0.35762995, 1.965), 'Tl1': (0.12628793, 1.682), 'Cu1': (0.001123, 1.201), 'Ag1': (0.00761745, 1.335), 'Be2': (6.21e-06, 0.971), 'Cu2': (0.0017408, 1.228), 'Ni2': (0.00247282, 1.251), 'Pt2': (0.00334975, 1.272), 'Zn2': (0.00374505, 1.28), 'Co2': (0.00530214, 1.306), 'Pd2': (0.00551135, 1.309), 'Ag2': (0.00799176, 1.339), 'Cr2': (0.00899152, 1.349), 'Fe2': (0.00974813, 1.356), 'Mg2': (0.01055378, 1.363), 'V2': (0.01091456, 1.366), 'Mn2': (0.0173834, 1.41), 'Hg2': (0.0173834, 1.41), 'Cd2': (0.01791152, 1.413), 'Yb2': (0.09072908, 1.622), 'Ca2': (0.09399072, 1.628), 'Sn2': (0.10128575, 1.641), 'Pb2': (0.15415012, 1.723), 'Eu2': (0.19078645, 1.772), 'Sr2': (0.1954949, 1.778), 'Sm2': (0.2002477, 1.784), 'Ba2': (0.35853865, 1.966), 'Ra2': (0.35853865, 1.966), 'Al3': (1.631e-05, 1.005), 'Fe3': (0.00025449, 1.121), 'Cr3': (0.00146124, 1.217), 'In3': (0.00212802, 1.241), 'Tl3': (0.00254709, 1.253), 'Y3': (0.02735051, 1.458), 'La3': (0.10417397, 1.646), 'Ce3': (0.08440707, 1.61), 'Pr3': (0.07117158, 1.583), 'Nd3': (0.06182717, 1.562), 'Sm3': (0.04953859, 1.531), 'Eu3': (0.04354662, 1.514), 'Gd3': (0.04058327, 1.505), 'Tb3': (0.03280986, 1.479), 'Dy3': (0.02960343, 1.467), 'Er3': (0.02174524, 1.433), 'Tm3': (0.01975917, 1.423), 'Lu3': (0.01975917, 1.423), 'Hf4': (0.0003306, 1.134), 'Zr4': (0.00064221, 1.169), 'Ce4': (0.01141046, 1.37), 'U4': (0.0023297, 1.247), 'Pu4': (0.00544088, 1.308), 'Th4': (0.03171494, 1.475)}

vdw_data['merz-spce'] = {'Li1': (0.00274091, 1.258), 'Na1': (0.02639002, 1.454), 'K1': (0.12693448, 1.683), 'Rb1': (0.20665151, 1.792), 'Cs1': (0.34673208, 1.953), 'Tl1': (0.11741683, 1.668), 'Cu1': (0.00096394, 1.192), 'Ag1': (0.0071693, 1.33), 'Be2': (4.6e-06, 0.961), 'Cu2': (0.0016086, 1.223), 'Ni2': (0.00254709, 1.253), 'Pt2': (0.00334975, 1.272), 'Zn2': (0.00354287, 1.276), 'Co2': (0.00523385, 1.305), 'Pd2': (0.00523385, 1.305), 'Ag2': (0.00780282, 1.337), 'Cr2': (0.00888732, 1.348), 'Fe2': (0.00952704, 1.354), 'Mg2': (0.01020237, 1.36), 'V2': (0.01079325, 1.365), 'Mn2': (0.0166976, 1.406), 'Hg2': (0.0166976, 1.406), 'Cd2': (0.01773416, 1.412), 'Yb2': (0.09731901, 1.634), 'Ca2': (0.09788018, 1.635), 'Sn2': (0.10710756, 1.651), 'Pb2': (0.1598965, 1.731), 'Eu2': (0.2018416, 1.786), 'Sr2': (0.20826406, 1.794), 'Sm2': (0.21312875, 1.8), 'Ba2': (0.37126402, 1.98), 'Ra2': (0.37126402, 1.98), 'Al3': (1.107e-05, 0.991), 'Fe3': (0.00013462, 1.091), 'Cr3': (0.00103208, 1.196), 'In3': (0.00128267, 1.209), 'Tl3': (0.00136949, 1.213), 'Y3': (0.02759452, 1.459), 'La3': (0.09454081, 1.629), 'Ce3': (0.07786298, 1.597), 'Pr3': (0.0657303, 1.571), 'Nd3': (0.0572627, 1.551), 'Sm3': (0.04772212, 1.526), 'Eu3': (0.04122946, 1.507), 'Gd3': (0.03868661, 1.499), 'Tb3': (0.03450196, 1.485), 'Dy3': (0.03091095, 1.472), 'Er3': (0.02236885, 1.436), 'Tm3': (0.02034021, 1.426), 'Lu3': (0.02034021, 1.426), 'Hf4': (0.00015685, 1.098), 'Zr4': (0.00044254, 1.149), 'Ce4': (0.01020237, 1.36), 'U4': (0.00148497, 1.218), 'Pu4': (0.00379705, 1.281), 'Th4': (0.02986171, 1.468)}

vdw_data['merz-tip3p-fb'] = {'Li1': (0.0028218, 1.26), 'Na1': (0.02759452, 1.459), 'K1': (0.14158262, 1.705), 'Rb1': (0.21475916, 1.802), 'Cs1': (0.36217399, 1.97), 'Tl1': (0.12244452, 1.676), 'Cu1': (0.00096394, 1.192), 'Ag1': (0.00699604, 1.328), 'Be2': (4.89e-06, 0.963), 'Cu2': (0.0016086, 1.223), 'Ni2': (0.0027012, 1.257), 'Pt2': (0.0034452, 1.274), 'Zn2': (0.00354287, 1.276), 'Co2': (0.00516628, 1.304), 'Pd2': (0.00544088, 1.308), 'Ag2': (0.00799176, 1.339), 'Cr2': (0.00899152, 1.349), 'Fe2': (0.00974813, 1.356), 'Mg2': (0.01055378, 1.363), 'V2': (0.01091456, 1.366), 'Mn2': (0.01755812, 1.411), 'Hg2': (0.01755812, 1.411), 'Cd2': (0.01827024, 1.415), 'Yb2': (0.09454081, 1.629), 'Ca2': (0.09788018, 1.635), 'Sn2': (0.1076997, 1.652), 'Pb2': (0.15917293, 1.73), 'Eu2': (0.19865859, 1.782), 'Sr2': (0.20424131, 1.789), 'Sm2': (0.20907204, 1.795), 'Ba2': (0.37399087, 1.983), 'Ra2': (0.37399087, 1.983), 'Al3': (6.39e-06, 0.972), 'Fe3': (0.00016028, 1.099), 'Cr3': (0.00104974, 1.197), 'In3': (0.00132548, 1.211), 'Tl3': (0.00132548, 1.211), 'Y3': (0.02133669, 1.431), 'La3': (0.09180886, 1.624), 'Ce3': (0.07399405, 1.589), 'Pr3': (0.06140287, 1.561), 'Nd3': (0.05292325, 1.54), 'Sm3': (0.04287573, 1.512), 'Eu3': (0.03655251, 1.492), 'Gd3': (0.03450196, 1.485), 'Tb3': (0.02759452, 1.459), 'Dy3': (0.02454281, 1.446), 'Er3': (0.01827024, 1.415), 'Tm3': (0.01773416, 1.412), 'Lu3': (0.01773416, 1.412), 'Hf4': (0.00019838, 1.109), 'Zr4': (0.00045105, 1.15), 'Ce4': (0.00930991, 1.352), 'U4': (0.0016086, 1.223), 'Pu4': (0.00401101, 1.285), 'Th4': (0.02639002, 1.454)}

vdw_data['merz-opc'] = {'Li1': (0.00216058, 1.242), 'Na1': (0.02960343, 1.467), 'K1': (0.13953816, 1.702), 'Rb1': (0.2279546, 1.818), 'Cs1': (0.35308749, 1.96), 'Tl1': (0.11250137, 1.66), 'Cu1': (0.00078213, 1.18), 'Ag1': (0.00602547, 1.316), 'Be2': (1.28e-06, 0.921), 'Cu2': (0.0007549, 1.178), 'Ni2': (0.00104974, 1.197), 'Pt2': (0.00150903, 1.219), 'Zn2': (0.00150903, 1.219), 'Co2': (0.00294683, 1.263), 'Pd2': (0.00321068, 1.269), 'Ag2': (0.00523385, 1.305), 'Cr2': (0.00602547, 1.316), 'Fe2': (0.00657749, 1.323), 'Mg2': (0.0071693, 1.33), 'V2': (0.00752608, 1.334), 'Mn2': (0.0128646, 1.381), 'Hg2': (0.0128646, 1.381), 'Cd2': (0.01400886, 1.389), 'Yb2': (0.08034231, 1.602), 'Ca2': (0.08337961, 1.608), 'Sn2': (0.09235154, 1.625), 'Pb2': (0.14295367, 1.707), 'Eu2': (0.17618319, 1.753), 'Sr2': (0.18612361, 1.766), 'Sm2': (0.18612361, 1.766), 'Ba2': (0.35308749, 1.96), 'Ra2': (0.35308749, 1.96), 'Al3': (1.42e-06, 0.924), 'Fe3': (2.635e-05, 1.023), 'Cr3': (0.00040986, 1.145), 'In3': (0.00051476, 1.157), 'Tl3': (0.00051476, 1.157), 'Y3': (0.0155469, 1.399), 'La3': (0.06841702, 1.577), 'Ce3': (0.05292325, 1.54), 'Pr3': (0.04188268, 1.509), 'Nd3': (0.03566355, 1.489), 'Sm3': (0.02735051, 1.458), 'Eu3': (0.02279185, 1.438), 'Gd3': (0.02133669, 1.431), 'Tb3': (0.01919059, 1.42), 'Dy3': (0.01721, 1.409), 'Er3': (0.01272679, 1.38), 'Tm3': (0.01166488, 1.372), 'Lu3': (0.01166488, 1.372), 'Hf4': (2.635e-05, 1.023), 'Zr4': (0.00013761, 1.092), 'Ce4': (0.00530214, 1.306), 'U4': (0.00049581, 1.155), 'Pu4': (0.00155814, 1.221), 'Th4': (0.01460944, 1.393)}

vdw_data['merz-tip4p-fb'] = {'Li1': (0.00209587, 1.24), 'Na1': (0.02499549, 1.448), 'K1': (0.1282327, 1.685), 'Rb1': (0.20907204, 1.795), 'Cs1': (0.34401021, 1.95), 'Tl1': (0.11617738, 1.666), 'Cu1': (0.00067804, 1.172), 'Ag1': (0.0058006, 1.313), 'Be2': (1.42e-06, 0.924), 'Cu2': (0.00089969, 1.181), 'Ni2': (0.00136949, 1.213), 'Pt2': (0.00176831, 1.229), 'Zn2': (0.00191142, 1.234), 'Co2': (0.00359255, 1.277), 'Pd2': (0.00384964, 1.282), 'Ag2': (0.00602547, 1.316), 'Cr2': (0.00691068, 1.327), 'Fe2': (0.00752608, 1.334), 'Mg2': (0.00828195, 1.342), 'V2': (0.00858042, 1.345), 'Mn2': (0.01476261, 1.394), 'Hg2': (0.01476261, 1.394), 'Cd2': (0.01538757, 1.398), 'Yb2': (0.09019198, 1.621), 'Ca2': (0.09344247, 1.627), 'Sn2': (0.10359269, 1.645), 'Pb2': (0.15773029, 1.728), 'Eu2': (0.19865859, 1.782), 'Sr2': (0.20424131, 1.789), 'Sm2': (0.20988115, 1.796), 'Ba2': (0.38216886, 1.992), 'Ra2': (0.38216886, 1.992), 'Al3': (1.08e-06, 0.916), 'Fe3': (2.024e-05, 1.013), 'Cr3': (0.00037198, 1.14), 'In3': (0.00043416, 1.148), 'Tl3': (0.00046851, 1.152), 'Y3': (0.01430674, 1.391), 'La3': (0.06706518, 1.574), 'Ce3': (0.05177853, 1.537), 'Pr3': (0.04090549, 1.506), 'Nd3': (0.03450196, 1.485), 'Sm3': (0.02759452, 1.459), 'Eu3': (0.02113456, 1.43), 'Gd3': (0.01937874, 1.421), 'Tb3': (0.01721, 1.409), 'Dy3': (0.01586934, 1.401), 'Er3': (0.01179373, 1.373), 'Tm3': (0.01067299, 1.364), 'Lu3': (0.01067299, 1.364), 'Hf4': (2.079e-05, 1.014), 'Zr4': (0.00012321, 1.087), 'Ce4': (0.00496778, 1.301), 'U4': (0.00043416, 1.148), 'Pu4': (0.00141473, 1.215), 'Th4': (0.01386171, 1.388)}

vdw_data['merz-tip4-ew'] = {'Li1': (0.00168686, 1.226), 'Na1': (0.02154025, 1.432), 'K1': (0.11803919, 1.669), 'Rb1': (0.18689752, 1.767), 'Cs1': (0.33132862, 1.936), 'Tl1': (0.10417397, 1.646), 'Cu1': (0.0005052, 1.156), 'Ag1': (0.00452863, 1.294), 'Be2': (1.16e-06, 0.918), 'Cu2': (0.00101467, 1.195), 'Ni2': (0.00155814, 1.221), 'Pt2': (0.00247282, 1.251), 'Zn2': (0.00250973, 1.252), 'Co2': (0.00417787, 1.288), 'Pd2': (0.00417787, 1.288), 'Ag2': (0.00657749, 1.323), 'Cr2': (0.00743559, 1.333), 'Fe2': (0.00838052, 1.343), 'Mg2': (0.00941798, 1.353), 'V2': (0.00941798, 1.353), 'Mn2': (0.01586934, 1.401), 'Hg2': (0.01586934, 1.401), 'Cd2': (0.0166976, 1.406), 'Yb2': (0.10888937, 1.654), 'Ca2': (0.11068733, 1.657), 'Sn2': (0.1186633, 1.67), 'Pb2': (0.1799796, 1.758), 'Eu2': (0.2321311, 1.823), 'Sr2': (0.2354895, 1.827), 'Sm2': (0.24480038, 1.838), 'Ba2': (0.43454345, 2.05), 'Ra2': (0.43454345, 2.05), 'Al3': (2.6e-07, 0.876), 'Fe3': (9.07e-06, 0.984), 'Cr3': (0.00015019, 1.096), 'In3': (0.0002026, 1.11), 'Tl3': (0.00022027, 1.114), 'Y3': (0.01205473, 1.375), 'La3': (0.05807581, 1.553), 'Ce3': (0.04525501, 1.519), 'Pr3': (0.03655251, 1.492), 'Nd3': (0.03064622, 1.471), 'Sm3': (0.02431873, 1.445), 'Eu3': (0.02014513, 1.425), 'Gd3': (0.01863432, 1.417), 'Tb3': (0.01619682, 1.403), 'Dy3': (0.01400886, 1.389), 'Er3': (0.00909668, 1.35), 'Tm3': (0.00808758, 1.34), 'Lu3': (0.00808758, 1.34), 'Hf4': (7.41e-06, 0.977), 'Zr4': (3.24e-05, 1.031), 'Ce4': (0.0027012, 1.257), 'U4': (0.00018227, 1.105), 'Pu4': (0.00067804, 1.172), 'Th4': (0.01141046, 1.37)}

vdw_data['zhang-tip3p'] = {'Li1': (1.5076, 0.7333), 'Na1': (0.0127, 1.6), 'K1': (0.2249, 1.7348), 'Rb1': (0.1769, 1.8774), 'Cs1': (0.284, 2.0537), 'Tl1': (2.3692, 1.6319), 'Cu1': (0.5857, 0.6199), 'Ag1': (14.2922, 0.9151), 'Ba2': (1.6919, 1.7195), 'Ca2': (4.4453, 1.185), 'Cd2': (76.5949, 0.7363), 'Co2': (198.9299, 0.5382), 'Cr2': (71.9284, 0.5855), 'Cu2': (446.5808, 0.4908), 'Fe2': (90.9076, 0.602), 'Hg2': (168.6165, 0.8205), 'Mg2': (51.3807, 0.6105), 'Mn2': (46.5908, 0.6944), 'Ni2': (298.4696, 0.486), 'Sn2': (26.5766, 1.1771), 'Sr2': (2.2044, 1.4518), 'V2': (120.0052, 0.6414), 'Zn2': (295.5289, 0.5152), 'Al3': (2666.2447, 0.212), 'Cr3': (415.1452, 0.3928), 'Fe3': (1710.803, 0.3503), 'In3': (605.8987, 0.4894), 'Lu3': (82.1297, 0.8412), 'Tl3': (834.0651, 0.5756), 'Tm3': (100.4847, 0.8466), 'Er3': (88.736, 0.8639), 'Y3': (60.7435, 0.9015), 'Dy3': (46.5479, 0.9256), 'Tb3': (54.5507, 0.9512), 'Gd3': (32.2478, 0.9842), 'Eu3': (65.8415, 0.9903), 'Sm3': (45.4255, 1.0266), 'Nd3': (41.7061, 1.0674), 'Pr3': (57.2532, 1.1017), 'Ce3': (39.1021, 1.1483), 'La3': (13.2221, 1.2345), 'Hf4': (3926.4494, 0.4671), 'Zr4': (2844.4611, 0.5073), 'Ce4': (999.3095, 0.8001), 'Pu4': (2644.8436, 0.7077), 'U4': (1907.2165, 0.7667), 'Th4': (426.5795, 0.9035)}

vdw_data['zhang-opc3'] = {'Li1': (0.6703, 0.8077), 'Na1': (0.0133, 1.6106), 'K1': (0.1771, 1.7663), 'Rb1': (0.0552, 2.0085), 'Cs1': (0.2087, 2.0848), 'Tl1': (1.7422, 1.6426), 'Cu1': (0.118, 0.7778), 'Ag1': (11.5027, 0.9236), 'Ba2': (0.7829, 1.805), 'Ca2': (1.3095, 1.295), 'Cd2': (52.2878, 0.7709), 'Co2': (108.1185, 0.5648), 'Cr2': (22.9192, 0.6567), 'Cu2': (301.3006, 0.5112), 'Fe2': (41.4381, 0.6565), 'Hg2': (113.6842, 0.8218), 'Mg2': (22.1361, 0.6804), 'Mn2': (23.19, 0.7607), 'Ni2': (114.8418, 0.4954), 'Sn2': (19.8335, 1.2056), 'Sr2': (1.5918, 1.5022), 'V2': (78.2348, 0.6715), 'Zn2': (150.2449, 0.5312), 'Al3': (1678.804, 0.2243), 'Cr3': (161.7335, 0.4436), 'Fe3': (1192.6143, 0.3648), 'In3': (352.9393, 0.4981), 'Lu3': (41.0866, 0.9002), 'Tl3': (550.9345, 0.5915), 'Tm3': (53.2599, 0.899), 'Er3': (43.2813, 0.9195), 'Y3': (23.9828, 0.9757), 'Dy3': (19.3108, 1.0048), 'Tb3': (23.3776, 1.0238), 'Gd3': (8.5409, 1.0948), 'Eu3': (31.8273, 1.0525), 'Sm3': (13.6458, 1.1178), 'Nd3': (12.471, 1.1646), 'Pr3': (23.3131, 1.1752), 'Ce3': (15.4739, 1.2307), 'La3': (4.5113, 1.3516), 'Hf4': (2876.0736, 0.4774), 'Zr4': (2061.1045, 0.5243), 'Ce4': (636.7955, 0.8143), 'Pu4': (1984.2665, 0.7149), 'U4': (1412.8628, 0.785), 'Th4': (366.7874, 0.9678)}

vdw_data['zhang-spce'] = {'Li1': (2.1717, 0.71), 'Na1': (0.0114, 1.6033), 'K1': (0.2769, 1.7032), 'Rb1': (0.2084, 1.8578), 'Cs1': (0.3099, 2.0312), 'Tl1': (2.1449, 1.6194), 'Cu1': (0.3061, 0.6745), 'Ag1': (15.7109, 0.9128), 'Ba2': (1.6244, 1.7202), 'Ca2': (3.4348, 1.2015), 'Cd2': (68.1868, 0.7385), 'Co2': (157.362, 0.545), 'Cr2': (35.0026, 0.6096), 'Cu2': (411.2444, 0.493), 'Fe2': (66.6807, 0.6125), 'Hg2': (169.1219, 0.8158), 'Mg2': (33.5274, 0.6326), 'Mn2': (44.1469, 0.7099), 'Ni2': (255.2114, 0.489), 'Sn2': (30.3319, 1.1766), 'Sr2': (2.7271, 1.4326), 'V2': (118.2225, 0.6463), 'Zn2': (247.5141, 0.5198), 'Al3': (2335.6082, 0.2179), 'Cr3': (315.6458, 0.4059), 'Fe3': (1580.5201, 0.3505), 'In3': (531.0067, 0.4869), 'Lu3': (72.2271, 0.8519), 'Tl3': (738.7543, 0.5767), 'Tm3': (85.6249, 0.8563), 'Er3': (74.7997, 0.8734), 'Y3': (46.7197, 0.9156), 'Dy3': (33.6589, 0.946), 'Tb3': (39.7283, 0.971), 'Gd3': (22.3409, 1.0068), 'Eu3': (55.5265, 1.0091), 'Sm3': (35.9501, 1.0476), 'Nd3': (31.7687, 1.0901), 'Pr3': (49.8884, 1.1145), 'Ce3': (31.5355, 1.1641), 'La3': (8.1527, 1.271), 'Hf4': (3480.1663, 0.4618), 'Zr4': (2521.1576, 0.5049), 'Ce4': (889.2011, 0.8033), 'Pu4': (2384.5147, 0.7035), 'U4': (1731.0116, 0.768), 'Th4': (443.4044, 0.9297)}

vdw_data['zhang-spceb'] = {'Li1': (1.8382, 0.7203), 'Na1': (0.0125, 1.5896), 'K1': (0.3077, 1.6896), 'Rb1': (0.2127, 1.8428), 'Cs1': (0.293, 2.0377), 'Tl1': (2.4508, 1.5983), 'Cu1': (0.2065, 0.6903), 'Ag1': (15.2405, 0.9), 'Ba2': (2.2121, 1.6888), 'Ca2': (3.2486, 1.1965), 'Cd2': (77.6426, 0.7482), 'Co2': (159.1842, 0.5503), 'Cr2': (31.681, 0.6154), 'Cu2': (399.7606, 0.4964), 'Fe2': (64.4615, 0.6179), 'Hg2': (172.4249, 0.823), 'Mg2': (28.0156, 0.6398), 'Mn2': (39.6096, 0.7157), 'Ni2': (239.3867, 0.493), 'Sn2': (27.5296, 1.1676), 'Sr2': (3.3551, 1.4155), 'V2': (110.9685, 0.6542), 'Zn2': (233.9915, 0.5228), 'Al3': (2238.2057, 0.22), 'Cr3': (297.5089, 0.4152), 'Fe3': (1533.9104, 0.3545), 'In3': (510.7402, 0.4826), 'Lu3': (63.183, 0.8594), 'Tl3': (751.9691, 0.5894), 'Tm3': (77.4462, 0.8624), 'Er3': (65.7204, 0.8806), 'Y3': (42.1114, 0.9234), 'Dy3': (34.5223, 0.95), 'Tb3': (41.505, 0.9739), 'Gd3': (27.6376, 1.0021), 'Eu3': (56.0661, 1.0121), 'Sm3': (38.7258, 1.0467), 'Nd3': (32.4788, 1.0901), 'Pr3': (47.1954, 1.1209), 'Ce3': (27.919, 1.1723), 'La3': (11.1506, 1.257), 'Hf4': (3346.5707, 0.4611), 'Zr4': (2510.7299, 0.5113), 'Ce4': (851.5301, 0.8044), 'Pu4': (2314.1953, 0.7044), 'U4': (1713.1682, 0.7766), 'Th4': (452.3759, 0.956)}

vdw_data['zhang-tip3p-fb'] = {'Li1': (2.1503, 0.7046), 'Na1': (0.0154, 1.5659), 'K1': (0.3673, 1.6659), 'Rb1': (0.2871, 1.8016), 'Cs1': (0.431, 1.9952), 'Tl1': (2.5811, 1.5803), 'Cu1': (0.2563, 0.6758), 'Ag1': (14.5345, 0.8941), 'Ba2': (1.8134, 1.6944), 'Ca2': (3.6191, 1.1859), 'Cd2': (82.8324, 0.7333), 'Co2': (161.9198, 0.5401), 'Cr2': (56.7545, 0.5986), 'Cu2': (412.3823, 0.4903), 'Fe2': (75.0758, 0.6116), 'Hg2': (177.787, 0.8174), 'Mg2': (44.0859, 0.6229), 'Mn2': (46.7627, 0.7084), 'Ni2': (232.0598, 0.4801), 'Sn2': (28.2944, 1.1647), 'Sr2': (4.254, 1.3975), 'V2': (119.9223, 0.644), 'Zn2': (242.8287, 0.512), 'Al3': (2162.2206, 0.2043), 'Cr3': (311.9608, 0.4068), 'Fe3': (1589.6445, 0.3437), 'In3': (533.2121, 0.4839), 'Lu3': (72.8115, 0.8485), 'Tl3': (793.5969, 0.5865), 'Tm3': (88.9816, 0.8529), 'Er3': (75.0931, 0.8697), 'Y3': (43.3112, 0.9102), 'Dy3': (30.8674, 0.9412), 'Tb3': (34.435, 0.967), 'Gd3': (23.142, 0.9987), 'Eu3': (45.6562, 1.006), 'Sm3': (34.914, 1.0417), 'Nd3': (35.2371, 1.082), 'Pr3': (51.2271, 1.111), 'Ce3': (27.9319, 1.1629), 'La3': (9.822, 1.2581), 'Hf4': (3624.0986, 0.4711), 'Zr4': (2511.3081, 0.4969), 'Ce4': (896.1898, 0.8029), 'Pu4': (2380.6743, 0.6951), 'U4': (1709.228, 0.7645), 'Th4': (460.6807, 0.9341)}

vdw_data['zhang-opc'] = {'Li1': (1.3624, 0.7035), 'Na1': (0.0155, 1.5203), 'K1': (0.4372, 1.6202), 'Rb1': (0.2929, 1.7728), 'Cs1': (0.5326, 1.9621), 'Tl1': (2.6164, 1.5967), 'Cu1': (0.5213, 0.5928), 'Ag1': (10.9169, 0.9009), 'Ba2': (1.5747, 1.6708), 'Ca2': (2.9567, 1.161), 'Cd2': (86.7561, 0.7194), 'Co2': (158.2705, 0.5308), 'Cr2': (41.3618, 0.5765), 'Cu2': (362.5768, 0.4883), 'Fe2': (65.7961, 0.5924), 'Hg2': (145.68, 0.7967), 'Mg2': (37.0937, 0.6011), 'Mn2': (48.3504, 0.6836), 'Ni2': (201.1407, 0.4657), 'Sn2': (22.8192, 1.1514), 'Sr2': (2.2331, 1.4013), 'V2': (121.5906, 0.6326), 'Zn2': (225.8396, 0.5089), 'Al3': (2032.357, 0.1756), 'Cr3': (363.58, 0.3811), 'Fe3': (1655.6733, 0.3565), 'In3': (646.2492, 0.47), 'Lu3': (97.4092, 0.8127), 'Tl3': (759.4514, 0.5578), 'Tm3': (112.3311, 0.8221), 'Er3': (97.0063, 0.8358), 'Y3': (65.4787, 0.8687), 'Dy3': (56.3378, 0.8898), 'Tb3': (60.062, 0.9186), 'Gd3': (39.9852, 0.9317), 'Eu3': (63.3578, 0.9627), 'Sm3': (38.1505, 0.9931), 'Nd3': (32.9837, 1.041), 'Pr3': (55.335, 1.0821), 'Ce3': (38.9762, 1.1209), 'La3': (12.9688, 1.2008), 'Hf4': (3388.4415, 0.4367), 'Zr4': (2609.1557, 0.4988), 'Ce4': (907.6115, 0.7823), 'Pu4': (2329.1636, 0.701), 'U4': (1604.3533, 0.7429), 'Th4': (483.9494, 0.9182)}

vdw_data['zhang-tip4p2005'] = {'Li1': (5.9279, 0.5808), 'Na1': (0.0943, 1.3029), 'K1': (0.9701, 1.5297), 'Rb1': (0.6275, 1.6748), 'Cs1': (1.1184, 1.886), 'Tl1': (3.9264, 1.5499), 'Cu1': (9.7163, 0.3825), 'Ag1': (20.54, 0.8595), 'Ba2': (4.5405, 1.5377), 'Ca2': (14.7028, 1.0174), 'Cd2': (160.8051, 0.6664), 'Co2': (351.4795, 0.4887), 'Cr2': (159.9926, 0.5046), 'Cu2': (653.2809, 0.4598), 'Fe2': (190.9414, 0.5299), 'Hg2': (245.5274, 0.7654), 'Mg2': (131.22, 0.522), 'Mn2': (120.9205, 0.6095), 'Ni2': (427.0709, 0.4234), 'Sn2': (47.5773, 1.0937), 'Sr2': (10.4159, 1.2608), 'V2': (232.8091, 0.5781), 'Zn2': (435.6121, 0.46), 'Al3': (4048.5546, 0.1798), 'Cr3': (966.0509, 0.3233), 'Fe3': (2855.6172, 0.3142), 'In3': (1245.3746, 0.4389), 'Lu3': (246.4336, 0.7435), 'Tl3': (1355.8137, 0.5405), 'Tm3': (270.0847, 0.7581), 'Er3': (245.1884, 0.7694), 'Y3': (174.5823, 0.7751), 'Dy3': (157.2896, 0.8035), 'Tb3': (156.9278, 0.8345), 'Gd3': (118.7682, 0.8394), 'Eu3': (158.4528, 0.8828), 'Sm3': (123.8226, 0.9023), 'Nd3': (111.9438, 0.9446), 'Pr3': (143.1528, 0.9986), 'Ce3': (109.5973, 1.0344), 'La3': (61.4328, 1.0647), 'Hf4': (5747.7782, 0.436), 'Zr4': (4388.3373, 0.4683), 'Ce4': (1623.6783, 0.7527), 'Pu4': (3659.3167, 0.692), 'U4': (2618.183, 0.7153), 'Th4': (910.1228, 0.8618)}

vdw_data['zhang-tip4p-d'] = {'Li1': (3.4229, 0.6141), 'Na1': (0.0172, 1.4781), 'K1': (0.5961, 1.5629), 'Rb1': (0.4431, 1.7104), 'Cs1': (0.7912, 1.9023), 'Tl1': (2.8987, 1.5536), 'Cu1': (3.0395, 0.4394), 'Ag1': (14.1254, 0.8734), 'Ba2': (4.4504, 1.5596), 'Ca2': (9.1517, 1.0465), 'Cd2': (123.1119, 0.6868), 'Co2': (231.8996, 0.4973), 'Cr2': (99.106, 0.5212), 'Cu2': (458.4585, 0.4547), 'Fe2': (125.5741, 0.5469), 'Hg2': (182.1378, 0.7693), 'Mg2': (91.7487, 0.5443), 'Mn2': (91.243, 0.6374), 'Ni2': (280.0271, 0.4251), 'Sn2': (34.9382, 1.1112), 'Sr2': (6.4491, 1.2967), 'V2': (175.4284, 0.5987), 'Zn2': (302.1343, 0.4726), 'Al3': (2894.6762, 0.1605), 'Cr3': (602.4208, 0.3304), 'Fe3': (2260.4765, 0.3338), 'In3': (986.7338, 0.4377), 'Lu3': (174.7029, 0.7533), 'Tl3': (1004.3845, 0.5158), 'Tm3': (196.1553, 0.7668), 'Er3': (175.0653, 0.7767), 'Y3': (130.1368, 0.799), 'Dy3': (110.8664, 0.8156), 'Tb3': (111.7892, 0.8512), 'Gd3': (79.5426, 0.8586), 'Eu3': (111.5578, 0.9053), 'Sm3': (82.054, 0.9242), 'Nd3': (71.5319, 0.9654), 'Pr3': (99.2659, 1.0217), 'Ce3': (70.0487, 1.0555), 'La3': (33.266, 1.0936), 'Hf4': (4467.8646, 0.4208), 'Zr4': (3488.1889, 0.478), 'Ce4': (1249.9712, 0.7483), 'Pu4': (2934.2695, 0.6821), 'U4': (2070.1413, 0.7163), 'Th4': (734.8522, 0.8976)}

vdw_data['zhang-tip4p-fb'] = {'Li1': (2.9909, 0.6524), 'Na1': (0.024, 1.4713), 'K1': (0.5925, 1.5857), 'Rb1': (0.3025, 1.7573), 'Cs1': (0.6637, 1.9342), 'Tl1': (3.1739, 1.5623), 'Cu1': (2.045, 0.483), 'Ag1': (18.5353, 0.8894), 'Ba2': (3.2606, 1.5981), 'Ca2': (6.6145, 1.0986), 'Cd2': (110.23, 0.7002), 'Co2': (246.3769, 0.5197), 'Cr2': (79.5243, 0.5475), 'Cu2': (473.9145, 0.4716), 'Fe2': (111.9438, 0.5694), 'Hg2': (197.0153, 0.7914), 'Mg2': (68.1868, 0.5699), 'Mn2': (66.4966, 0.66), 'Ni2': (311.6736, 0.4578), 'Sn2': (36.1826, 1.127), 'Sr2': (6.1475, 1.327), 'V2': (160.1769, 0.6032), 'Zn2': (319.0803, 0.4959), 'Al3': (2953.9286, 0.1988), 'Cr3': (524.8074, 0.3596), 'Fe3': (2065.8558, 0.3344), 'In3': (809.655, 0.4636), 'Lu3': (123.3105, 0.7864), 'Tl3': (959.8425, 0.5585), 'Tm3': (142.7579, 0.7977), 'Er3': (122.4616, 0.8103), 'Y3': (86.0994, 0.8427), 'Dy3': (71.796, 0.8642), 'Tb3': (80.2417, 0.8962), 'Gd3': (57.2137, 0.9095), 'Eu3': (104.472, 0.9465), 'Sm3': (76.4892, 0.9672), 'Nd3': (59.6348, 1.0089), 'Pr3': (79.7995, 1.0544), 'Ce3': (57.0821, 1.0932), 'La3': (25.1594, 1.1473), 'Hf4': (4257.9445, 0.4412), 'Zr4': (3207.7462, 0.4869), 'Ce4': (1163.858, 0.783), 'Pu4': (2866.1571, 0.694), 'U4': (2011.8704, 0.7375), 'Th4': (642.3918, 0.903)}

vdw_data['zhang-tip4p-ew'] = {'Li1': (4.9774, 0.5963), 'Na1': (0.0807, 1.3369), 'K1': (0.8204, 1.5657), 'Rb1': (0.5415, 1.7111), 'Cs1': (0.9262, 1.928), 'Tl1': (4.1448, 1.5581), 'Cu1': (5.6624, 0.4182), 'Ag1': (21.4388, 0.8691), 'Ba2': (5.6208, 1.5526), 'Ca2': (15.9845, 1.0419), 'Cd2': (145.3784, 0.6793), 'Co2': (338.6882, 0.4938), 'Cr2': (140.1522, 0.514), 'Cu2': (648.0373, 0.457), 'Fe2': (169.4728, 0.5378), 'Hg2': (240.4916, 0.7687), 'Mg2': (129.5686, 0.5343), 'Mn2': (109.4964, 0.6158), 'Ni2': (423.4479, 0.4349), 'Sn2': (47.1954, 1.1031), 'Sr2': (7.9763, 1.2804), 'V2': (211.7386, 0.5849), 'Zn2': (424.3263, 0.4663), 'Al3': (2668.0871, 0.1565), 'Cr3': (534.1952, 0.3334), 'Fe3': (2144.8651, 0.3446), 'In3': (936.9146, 0.4408), 'Lu3': (162.6672, 0.7556), 'Tl3': (944.4957, 0.5165), 'Tm3': (183.4003, 0.7704), 'Er3': (162.8921, 0.7802), 'Y3': (122.377, 0.8023), 'Dy3': (99.2888, 0.8196), 'Tb3': (104.1118, 0.8554), 'Gd3': (75.9277, 0.8674), 'Eu3': (102.5416, 0.9106), 'Sm3': (73.3499, 0.9203), 'Nd3': (63.9882, 0.9709), 'Pr3': (86.4569, 1.0264), 'Ce3': (60.562, 1.0579), 'La3': (31.9669, 1.0971), 'Hf4': (4140.9501, 0.4234), 'Zr4': (3181.2661, 0.4714), 'Ce4': (1177.606, 0.7574), 'Pu4': (2747.8942, 0.6848), 'U4': (1927.0811, 0.7218), 'Th4': (673.2867, 0.8963)}

vdw_data['uff'] = {'Ac': (0.033, 1.739), 'Ag': (0.036, 1.574), 'Al': (0.505, 2.2495), 'Am': (0.014, 1.6905), 'Ar': (0.185, 1.934),
                        'As': (0.309, 2.115), 'At': (0.284, 2.375), 'Au': (0.039, 1.6465), 'B': (0.18, 2.0415), 'Ba': (0.364, 1.8515),
                        'Be': (0.085, 1.3725), 'Bi': (0.518, 2.185), 'Bk': (0.013, 1.6695), 'Br': (0.251, 2.0945), 'C': (0.105, 1.9255),
                        'Ca': (0.238, 1.6995), 'Cd': (0.228, 1.424), 'Ce': (0.013, 1.778), 'Cf': (0.013, 1.6565), 'Cl': (0.227, 1.9735),
                        'Cm': (0.013, 1.663), 'Co': (0.014, 1.436), 'Cr': (0.015, 1.5115), 'Cs': (0.045, 2.2585),'Cu': (0.005, 1.7475),
                        'Dy': (0.007, 1.714), 'Er': (0.007, 1.6955), 'Es': (0.012, 1.6495), 'Eu': (0.008, 1.7465),'F': (0.05, 1.682),
                        'Fe': (0.013, 1.456), 'Fm': (0.012, 1.643), 'Fr': (0.05, 2.45), 'Ga': (0.415, 2.1915),'Gd': (0.009, 1.684),
                        'Ge': (0.379, 2.14), 'Hf': (0.072, 1.5705), 'Hg': (0.385, 1.3525), 'Ho': (0.007, 1.7045),'I': (0.339, 2.25),
                        'In': (0.599, 2.2315), 'Ir': (0.073, 1.42), 'Kr': (0.22, 2.0705), 'La': (0.017, 1.761),'Li': (0.025, 1.2255),
                        'Lu': (0.041, 1.82), 'Lw': (0.011, 1.618), 'Md': (0.011, 1.637), 'Mg': (0.111, 1.5105),'Mn': (0.013, 1.4805),
                        'Mo': (0.056, 1.526), 'N': (0.069, 1.83), 'Na': (0.03, 1.4915), 'Nb': (0.059, 1.5825),'Nd': (0.01, 1.7875),
                        'Ni': (0.015, 1.417), 'No': (0.011, 1.624), 'Np': (0.019, 1.712), 'O': (0.06, 1.75),'Os': (0.037, 1.56),
                        'Pa': (0.022, 1.712), 'Pb': (0.663, 2.1485), 'Pd': (0.048, 1.4495), 'Pm': (0.009, 1.7735),'Po': (0.325, 2.3545),
                        'Pr': (0.01, 1.803), 'Pt': (0.08, 1.377), 'Pu': (0.016, 1.712), 'Ra': (0.404, 1.8385),'Rb': (0.04, 2.057),
                        'Re': (0.066, 1.477), 'Rh': (0.053, 1.4645), 'Rn': (0.248, 2.3825), 'Ru': (0.056, 1.4815),'Sb': (0.449, 2.21),
                        'Sc': (0.019, 1.6475), 'Se': (0.291, 2.1025), 'Si': (0.402, 2.1475), 'Sm': (0.008, 1.76),'Sn': (0.567, 2.196),
                        'Sr': (0.291, 1.8205), 'Ta': (0.081, 1.585), 'Tb': (0.007, 1.7255), 'Tc': (0.048, 1.499),'Te': (0.398, 2.235),
                        'Th': (0.026, 1.698), 'Ti': (0.017, 1.5875), 'Tl': (0.68, 2.1735), 'Tm': (0.006, 1.687),'U': (0.022, 1.6975),
                        'V': (0.016, 1.572), 'W': (0.067, 1.5345), 'Xe': (0.332, 2.202), 'Y': (0.072, 1.6725),'Yb': (0.228, 1.6775),
                        'Zn': (0.124, 1.3815), 'Zr': (0.069, 1.562)}