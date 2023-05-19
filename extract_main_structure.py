import argparse
import MDAnalysis
import networkx as nx
import re
from mapping import map_two_structures


def extract_structure(filename, output=None):  # TODO this will be delated
    metals = ['pd', 'co', 'fe']
    syst = MDAnalysis.Universe(filename)

    # Remove numbers from name
    for atom in syst.atoms:
        print(atom.type, atom.name, atom.type.lower())
        atom.name = re.match('([a-zA-Z]+)', atom.name).group(1)
        if '.' in atom.type:
            atom.type = re.match('([a-zA-Z]+)\.', atom.type).group(1)

        if atom.name.lower() in metals:
            atom.type = "M"
            print("YES")

    G_fingerprint = nx.Graph(
        MDAnalysis.topology.guessers.guess_bonds(syst.atoms, syst.atoms.positions, vdwradii={"M": 3.0, "Cl": 1.735}))
    connected = nx.connected_components(G_fingerprint)
    largest_cc = max(nx.connected_components(G_fingerprint), key=len)

    if output is not None:
        syst.atoms[list(largest_cc)].write(output)

    return syst.atoms[list(largest_cc)]






def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", help="Metaloorganic structre (*gro, *pdb, etc. all supported by MDAnalysis)")
    parser.add_argument("-o", help="Clean structure")
    return parser.parse_args()



if __name__ == '__main__':
    args = get_args()
    if args.o is not None:
        extract_structure(args.f, args.o)


