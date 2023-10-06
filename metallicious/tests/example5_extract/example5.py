from metallicious.extract_metal_site import extract_metal_structure

f='pd2l4_tall.xyz'
metal='Pd'
extract_metal_structure(f, metal, "pd2l4_tall_fp1")

f='pd2l4_tall.pdb'
metal='Pd'
extract_metal_structure(f, metal, "pd2l4_tall_fp2")

f='pd2l4_tall.gro'
metal='Pd'
extract_metal_structure(f, metal, "pd2l4_tall_fp3")

f='Fe4L6.pdb'
metal='Fe'
extract_metal_structure(f, metal, "Fe4L6_fp")

f='Fe4L6_2.pdb'
metal='Fe'
extract_metal_structure(f, metal, "Fe4L6_2_fp")

f='Ga4L6.xyz'
metal='Ga'
extract_metal_structure(f, metal, "Ga4L6_fp")

f='Pd4L6.xyz'
metal='Pd'
extract_metal_structure(f, metal, "Pd4L6_fp")

f='knot_4n.pdb'
metal='Fe'
extract_metal_structure(f, metal, "knot_4n_fp")

f='knot_3n.pdb'
metal='Zn'
extract_metal_structure(f, metal, "knot_3n_fp")

f='Co8L16.xyz'
metal='Co'
extract_metal_structure(f, metal, "Co8L16_fp")

f='ru_pd.xyz'
metal='Ru'
extract_metal_structure(f, metal, "ru_fp")
