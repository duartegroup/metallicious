
try:
    from metallicious.main import patcher
except:
    from cgbind2pmd.main import cgbind2pmd



try:
    from cgbind2pmd.load_fingerprint import load_fingerprint_from_file, guess_fingerprint, reduce_site_to_fingerprint
    from cgbind2pmd.prepare_initial_topology import prepare_initial_topology
except:
    from metallicious.load_fingerprint import load_fp_from_file, guess_fingerprint, reduce_site_to_fingerprint
#f='start.pdb'

f='cage2.xyz'
linker_topol="linker0.top"

metal='Pd'
metal_charge=2
fingerprint='Pd2d'
fingerprint_style='dih'

cgbind2gmx = patcher()




fp1 = guess_fingerprint('cage.xyz', 0, metal_name = metal, fingerprint_guess_list=['Pd2d', 'PdB1', 'PdB2'])
fp2 = guess_fingerprint('cage2.xyz', 184, metal_name = metal, fingerprint_guess_list=['Pd2d', 'PdB1', 'PdB2'])
fp3 = guess_fingerprint('cage3.xyz', 160, metal_name = metal, fingerprint_guess_list=['Pd2d', 'PdB1', 'PdB2'])


if fp1 == 'Pd2d':
    print("[+] FP1 correct!")
else:
    print("[-] FP1 incorrect!")


if fp2 == 'PdB1':
    print("[+] FP2 correct!")
else:
    print("[-] FP2 incorrect!")



if fp3 == 'Pd2d':
    print("[+] FP3 correct!")
else:
    print("[-] FP3 incorrect!")

