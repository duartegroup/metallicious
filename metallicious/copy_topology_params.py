import parmed as pmd
import re

# try:
#    from cgbind2pmd.log import logger
# except:
from metallicious.log import logger
from metallicious.utils import strip_numbers_from_atom_name
import networkx as nx


'''
def adjust_charge(topol_new, fp_topol, mapping_fp_to_new):
    logger.info("   [ ] Changing charges and atomtypes")
    #sum_of_charge_diffrences = 0

    for a in mapping_fp_to_new:
        logger.info(f"          {topol_new.atoms[mapping_fp_to_new[a]].type:s} "
                    f"{topol_new.atoms[mapping_fp_to_new[a]].name:s} "
                    f"{topol_new.atoms[mapping_fp_to_new[a]].charge:} --> "
                    f"{fp_topol.atoms[a].type, fp_topol.atoms[a].name:} {topol_new.atoms[mapping_fp_to_new[a]].charge + fp_topol.atoms[a].charge:}")

        atom = fp_topol.atoms[a]

        if atom.type not in topol_new.parameterset.atom_types.keys():
            logger.info(f"      [^] Adding new atomtype: {atom.type:s}")
            atomtype = pmd.topologyobjects.AtomType(name=atom.type, number=atom.number, mass=atom.mass)
            topol_new.parameterset.atom_types[atom.type] = atomtype

        topol_new.atoms[mapping_fp_to_new[a]].type = fp_topol.atoms[a].type
        topol_new.atoms[mapping_fp_to_new[a]].epsilon = fp_topol.atoms[a].epsilon
        topol_new.atoms[mapping_fp_to_new[a]].sigma = fp_topol.atoms[a].sigma
        topol_new.atoms[mapping_fp_to_new[a]].rmin = fp_topol.atoms[a].rmin
        topol_new.atoms[mapping_fp_to_new[a]].charge += fp_topol.atoms[a].charge  # topol_fp.atoms[a].charge

        #sum_of_charge_diffrences += fp_topol.atoms[a].charge

    return topol_new
'''
def adjust_charge(topol_new, fp_topol, mapping_fp_to_new):
    logger.info("   [ ] Changing charges ")
    sum_of_charge_diffrences = 0

    for a in mapping_fp_to_new:
        logger.info(f"          {topol_new.atoms[mapping_fp_to_new[a]].type:s} "
                    f"{topol_new.atoms[mapping_fp_to_new[a]].name:s} "
                    f"{topol_new.atoms[mapping_fp_to_new[a]].charge:} --> "
                    f"{fp_topol.atoms[a].type, fp_topol.atoms[a].name:} {topol_new.atoms[mapping_fp_to_new[a]].charge + fp_topol.atoms[a].charge:}")

        atom = fp_topol.atoms[a]

        topol_new.atoms[mapping_fp_to_new[a]].charge += fp_topol.atoms[a].charge  # topol_fp.atoms[a].charge


    return topol_new


def adjust_bonds(topol_new, topol_fp, mapping_fp_to_new):
    logger.info("   [ ] Adding new bonds to topology")
    for bond_fp in topol_fp.bonds:
        found = False
        for bond_new in topol_new.bonds:

            # Check if atoms in mapping (that they are not the additional atoms)
            if bond_fp.atom1.idx in mapping_fp_to_new and bond_fp.atom2.idx in mapping_fp_to_new:
                if ((bond_new.atom1.idx == mapping_fp_to_new[bond_fp.atom1.idx] and bond_new.atom2.idx ==
                     mapping_fp_to_new[bond_fp.atom2.idx]) or
                        ((bond_new.atom1.idx == mapping_fp_to_new[bond_fp.atom2.idx] and bond_new.atom2.idx ==
                          mapping_fp_to_new[bond_fp.atom1.idx]))):
                    if (bond_fp.type != bond_new.type):
                        logger.info(f"      [o] Diffrent bond type  {bond_new.type:} -> {bond_fp.type:}")
                        logger.info(
                            f"          Mapping ({bond_fp.atom1.idx},{bond_fp.atom2.idx}) to ({bond_new.atom1.idx:},{bond_new.atom2.idx:})")
                        logger.info(
                            f"          Mapping ({bond_fp.atom1.name},{bond_fp.atom2.name}) to ({bond_new.atom1.name:},{bond_new.atom2.name:})")
                        type_to_assign = pmd.topologyobjects.BondType(bond_fp.type.k, bond_fp.type.req,
                                                                      list=topol_new.bond_types)
                        # if type_to_assign not in topol_new.bond_types:
                        topol_new.bond_types.append(type_to_assign)

                        # bond_new.type = deepcopy(bond_fp.type)
                        bond_new.type = type_to_assign  # bond_fp.type

    logger.info("   [ ] Adding new bonds to topology")

    for bond_fp in topol_fp.bonds:
        found = False

        for bond_new in topol_new.bonds:
            # if bond_fp.atom1.name=="ZN":
            #    print(bond_fp.atom1.idx , bond_new)

            # print(bond_new, bond_new.atom1.name)
            # Check if atoms in mapping (that they are not the additional atoms)
            if bond_fp.atom1.idx in mapping_fp_to_new and bond_fp.atom2.idx in mapping_fp_to_new:
                # print("It is in the mapping")

                if ((bond_new.atom1.idx == mapping_fp_to_new[bond_fp.atom1.idx] and bond_new.atom2.idx ==
                     mapping_fp_to_new[bond_fp.atom2.idx]) or
                        ((bond_new.atom1.idx == mapping_fp_to_new[bond_fp.atom2.idx] and bond_new.atom2.idx ==
                          mapping_fp_to_new[bond_fp.atom1.idx]))):
                    if (bond_fp.type != bond_new.type):
                        logger.info(f"      [o] Diffrent bond type {bond_new.type:} {bond_fp.type:}")
                    found = True
            else:
                found = True

        if not found:
            # type_to_assign= bond_fp.type
            # type_to_assign = pmd.topologyobjects.AngleType(angle_fp.type.k, angle_fp.type.theteq)
            # print("AAA", bond_fp.type.k, bond_fp.type.req, len(topol_new.bond_types))

            type_to_assign = pmd.topologyobjects.BondType(bond_fp.type.k, bond_fp.type.req,
                                                          list=topol_new.bond_types)

            # if type_to_assign not in topol_new.bond_types:
            topol_new.bond_types.append(type_to_assign)

            atom1 = topol_new.atoms[mapping_fp_to_new[bond_fp.atom1.idx]]
            atom2 = topol_new.atoms[mapping_fp_to_new[bond_fp.atom2.idx]]
            topol_new.bonds.append(pmd.topologyobjects.Bond(atom1, atom2, type=type_to_assign))

            logger.info(
                f"      [o] New bond: {mapping_fp_to_new[bond_fp.atom1.idx]:}  {type_to_assign:} {atom1.name:}-{atom2.name:} ({mapping_fp_to_new[bond_fp.atom2.idx]:})")

    return topol_new


def adjust_angles(topol_new, topol_fp, mapping_fp_to_new):
    logger.info("   [ ] Adding new angles to topology")
    for angle_fp in topol_fp.angles:
        found = False
        if angle_fp.atom1.idx in mapping_fp_to_new and angle_fp.atom2.idx in mapping_fp_to_new and angle_fp.atom3.idx in mapping_fp_to_new:
            for angle_new in topol_new.angles:
                # Check if atoms in mapping (that they are not the additional atoms)
                if ((angle_new.atom1.idx == mapping_fp_to_new[angle_fp.atom1.idx] and angle_new.atom2.idx ==
                     mapping_fp_to_new[angle_fp.atom2.idx] and angle_new.atom3.idx == mapping_fp_to_new[
                         angle_fp.atom3.idx]) or
                        ((angle_new.atom1.idx == mapping_fp_to_new[angle_fp.atom3.idx] and angle_new.atom2.idx ==
                          mapping_fp_to_new[angle_fp.atom2.idx] and angle_new.atom3.idx == mapping_fp_to_new[
                              angle_fp.atom1.idx]))):
                    if (angle_fp.type != angle_new.type):
                        logger.info(f"      [o] Diffrent angle type {angle_new.type:} -> {angle_fp.type:}")
                        # angle_new.funct = angle_fp.funct
                        # angle_new.type.k = angle_fp.type.k
                        # angle_new.type.theteq = angle_fp.type.theteq

                        type_to_assign = pmd.topologyobjects.AngleType(angle_fp.type.k, angle_fp.type.theteq,
                                                                       list=topol_new.angle_types)

                        # if type_to_assign not in topol_new.angle_types:
                        topol_new.angle_types.append(type_to_assign)
                        angle_new.type = type_to_assign

                    found = True
        else:
            logger.info(
                f"[-] Not found in fingerprint: {angle_fp.atom1.idx + 1:d} {angle_fp.atom2.idx + 1:d} {angle_fp.atom3.idx + 1:d}")
            logger.info(f"              {topol_fp.atoms[angle_fp.atom1.idx]:}")
            logger.info(f"              {topol_fp.atoms[angle_fp.atom2.idx]:}")
            logger.info(f"              {topol_fp.atoms[angle_fp.atom3.idx]:}")
            logger.info(f"But it is standard bond so it should be there")
            logger.info(f"And you should be worry if it is not the end of fingerprint")
            found = True

        if not found:
            # type_to_assign= angle_fp.type
            type_to_assign = pmd.topologyobjects.AngleType(angle_fp.type.k, angle_fp.type.theteq,
                                                           list=topol_new.angle_types)
            # if type_to_assign not in topol_new.angle_types:
            topol_new.angle_types.append(type_to_assign)

            logger.info(f"      [x] New angle:{angle_fp.atom1.name:}-{angle_fp.atom2.name:}-{angle_fp.atom3.name:} "
                        f"{mapping_fp_to_new[angle_fp.atom1.idx]:} {mapping_fp_to_new[angle_fp.atom2.idx]:} "
                        f"{mapping_fp_to_new[angle_fp.atom3.idx]:} {type_to_assign:}")
            # print("              ", topol_fp.atoms[angle_fp.atom1.idx])
            # print("              ", topol_fp.atoms[angle_fp.atom2.idx])
            # print("              ", topol_fp.atoms[angle_fp.atom3.idx])
            atom1 = topol_new.atoms[mapping_fp_to_new[angle_fp.atom1.idx]]
            atom2 = topol_new.atoms[mapping_fp_to_new[angle_fp.atom2.idx]]
            atom3 = topol_new.atoms[mapping_fp_to_new[angle_fp.atom3.idx]]

            topol_new.angles.append(pmd.topologyobjects.Angle(atom1, atom2, atom3, type=type_to_assign))

    return topol_new


def adjust_dihedrals(topol_new, topol_fp, mapping_fp_to_new):
    logger.info("   [ ] Adding new dihedrals to topology")
    for dihedral_fp in topol_fp.dihedrals:
        found = False

        # One combinations of atoms can have several dihedrals, therfore we need to iterate them (template)
        if isinstance(dihedral_fp.type, list):
            dihedral_fp_types = dihedral_fp.type
        else:
            dihedral_fp_types = [dihedral_fp.type]

        for dihedral_fp_type in dihedral_fp_types:

            if dihedral_fp.atom1.idx in mapping_fp_to_new and dihedral_fp.atom2.idx in mapping_fp_to_new and dihedral_fp.atom3.idx in mapping_fp_to_new and dihedral_fp.atom4.idx in mapping_fp_to_new:
                for dihedral_new in topol_new.dihedrals:
                    # Check if atoms in mapping (that they are not the additional atoms)

                    # One combinations of atoms can have several dihedrals, therfore we need to iterate them (binding site)
                    if isinstance(dihedral_new.type, list):
                        dihedral_new_types = dihedral_new.type
                    else:
                        dihedral_new_types = [dihedral_new.type]

                    for dihedral_new_type in dihedral_new_types:

                        if (((dihedral_new.atom1.idx == mapping_fp_to_new[
                            dihedral_fp.atom1.idx] and dihedral_new.atom2.idx ==
                              mapping_fp_to_new[dihedral_fp.atom2.idx] and dihedral_new.atom3.idx == mapping_fp_to_new[
                                  dihedral_fp.atom3.idx] and dihedral_new.atom4.idx == mapping_fp_to_new[
                                  dihedral_fp.atom4.idx]) or
                             ((dihedral_new.atom1.idx == mapping_fp_to_new[
                                 dihedral_fp.atom4.idx] and dihedral_new.atom2.idx ==
                               mapping_fp_to_new[dihedral_fp.atom3.idx] and dihedral_new.atom3.idx == mapping_fp_to_new[
                                   dihedral_fp.atom2.idx] and dihedral_new.atom4.idx == mapping_fp_to_new[
                                   dihedral_fp.atom1.idx]))) and
                                (dihedral_new_type.per == dihedral_fp_type.per)):
                            if (dihedral_fp_type != dihedral_new_type):
                                logger.info(
                                    f"      [a] Diffrent Dihedral type: {dihedral_fp.atom1.name:}-{dihedral_fp.atom2.name:}-{dihedral_fp.atom3.name:}-{dihedral_fp.atom4.name:}")
                                logger.info(f"           {dihedral_new_type:}")
                                logger.info(f"        -> {dihedral_fp_type:}")

                                type_to_assign = pmd.topologyobjects.DihedralType(phi_k=dihedral_fp_type.phi_k,
                                                                                  per=dihedral_fp_type.per,
                                                                                  phase=dihedral_fp_type.phase,
                                                                                  scee=dihedral_fp_type.scee,
                                                                                  scnb=dihedral_fp_type.scnb,
                                                                                  list=topol_new.dihedral_types)
                                # if type_to_assign not in topol_new.dihedral_types:
                                topol_new.dihedral_types.append(type_to_assign)
                                dihedral_new_type.type = type_to_assign

                            found = True
            else:
                logger.info(
                    f"[-] Not found in finger print: {dihedral_fp.atom1.idx + 1:} {dihedral_fp.atom2.idx + 1:} {dihedral_fp.atom3.idx + 1:} {dihedral_fp.atom4.idx + 1:}")
                logger.info(f"               {topol_fp.atoms[dihedral_fp.atom1.idx]:}")
                logger.info(f"               {topol_fp.atoms[dihedral_fp.atom2.idx]:}")
                logger.info(f"               {topol_fp.atoms[dihedral_fp.atom3.idx]:}")
                logger.info(f"               {topol_fp.atoms[dihedral_fp.atom4.idx]:}")
                logger.info("But it is standard so it should be there")
                logger.info("And you should be worry if it is not the end of fingerprint")
                found = True

                '''
                for dihedral in topol_fp.dihedrals:
                if dihedral.atom1.idx==0:
                print(dihedral.atom1.idx+1, dihedral.atom2.idx+1, dihedral.atom3.idx+1, dihedral.atom4.idx+1)
                '''

            if not found:
                type_to_assign = pmd.topologyobjects.DihedralType(phi_k=dihedral_fp.type.phi_k,
                                                                  per=dihedral_fp.type.per,
                                                                  phase=dihedral_fp.type.phase,
                                                                  scee=dihedral_fp.type.scee,
                                                                  scnb=dihedral_fp.type.scnb,
                                                                  list=topol_new.dihedral_types)
                # if type_to_assign not in topol_new.dihedral_types:
                topol_new.dihedral_types.append(type_to_assign)
                logger.info(
                    f"      [b] New dihedral: {dihedral_fp.atom1.name:}-{dihedral_fp.atom2.name:}-{dihedral_fp.atom3.name:}-{dihedral_fp.atom4.name:} "
                    f"{mapping_fp_to_new[dihedral_fp.atom1.idx]:} {mapping_fp_to_new[dihedral_fp.atom2.idx]:} "
                    f"{mapping_fp_to_new[dihedral_fp.atom3.idx]:} {type_to_assign:}")

                atom1 = topol_new.atoms[mapping_fp_to_new[dihedral_fp.atom1.idx]]
                atom2 = topol_new.atoms[mapping_fp_to_new[dihedral_fp.atom2.idx]]
                atom3 = topol_new.atoms[mapping_fp_to_new[dihedral_fp.atom3.idx]]
                atom4 = topol_new.atoms[mapping_fp_to_new[dihedral_fp.atom4.idx]]

                topol_new.dihedrals.append(
                    pmd.topologyobjects.Dihedral(atom1, atom2, atom3, atom4, type=type_to_assign))

    return topol_new


def adjust_impropers(topol_new, topol_fp, mapping_fp_to_new):
    # This is almost the same as above, with the diffrence on being improper dihedrals

    logger.info("   [ ] Adding new improper dihedrals to topology")
    for dihedral_fp in topol_fp.impropers:
        found = False
        if dihedral_fp.atom1.idx in mapping_fp_to_new and dihedral_fp.atom2.idx in mapping_fp_to_new and dihedral_fp.atom3.idx in mapping_fp_to_new and dihedral_fp.atom4.idx in mapping_fp_to_new:
            for dihedral_new in topol_new.impropers:
                # Check if atoms in mapping (that they are not the additional atoms)
                if (((dihedral_new.atom1.idx == mapping_fp_to_new[dihedral_fp.atom1.idx] and dihedral_new.atom2.idx ==
                      mapping_fp_to_new[dihedral_fp.atom2.idx] and dihedral_new.atom3.idx == mapping_fp_to_new[
                          dihedral_fp.atom3.idx] and dihedral_new.atom4.idx == mapping_fp_to_new[
                          dihedral_fp.atom4.idx]) or
                     ((dihedral_new.atom1.idx == mapping_fp_to_new[dihedral_fp.atom4.idx] and dihedral_new.atom2.idx ==
                       mapping_fp_to_new[dihedral_fp.atom3.idx] and dihedral_new.atom3.idx == mapping_fp_to_new[
                           dihedral_fp.atom2.idx] and dihedral_new.atom4.idx == mapping_fp_to_new[
                           dihedral_fp.atom1.idx]))) and
                        (dihedral_new.type.per == dihedral_fp.type.per)):
                    if (dihedral_fp.type != dihedral_new.type):
                        logger.info(
                            f"      [a] Diffrent Dihedral type: {dihedral_fp.atom1.name:}-{dihedral_fp.atom2.name:}-{dihedral_fp.atom3.name:}-{dihedral_fp.atom4.name:}")
                        logger.info(f"           {dihedral_new.type:}")
                        logger.info(f"        -> {dihedral_fp.type:}")

                        type_to_assign = pmd.topologyobjects.ImproperType(psi_k=dihedral_fp.type.psi_k,
                                                                          psi_eq=dihedral_fp.type.psi_eq,
                                                                          list=topol_new.improper_types)
                        # if type_to_assign not in topol_new.dihedral_types:
                        topol_new.improper_types.append(type_to_assign)
                        dihedral_new.type = type_to_assign

                    found = True
        else:
            logger.info(
                f"[-] Not found in finger print: {dihedral_fp.atom1.idx + 1:} {dihedral_fp.atom2.idx + 1:} {dihedral_fp.atom3.idx + 1:} {dihedral_fp.atom4.idx + 1:}")
            logger.info(f"               {topol_fp.atoms[dihedral_fp.atom1.idx]:}")
            logger.info(f"               {topol_fp.atoms[dihedral_fp.atom2.idx]:}")
            logger.info(f"               {topol_fp.atoms[dihedral_fp.atom3.idx]:}")
            logger.info(f"               {topol_fp.atoms[dihedral_fp.atom4.idx]:}")
            logger.info("But it is standard so it should be there")
            logger.info("And you should be worry if it is not the end of fingerprint")
            found = True

            '''
            for dihedral in topol_fp.dihedrals:
            if dihedral.atom1.idx==0:
            print(dihedral.atom1.idx+1, dihedral.atom2.idx+1, dihedral.atom3.idx+1, dihedral.atom4.idx+1)
            '''

        if not found:
            type_to_assign = pmd.topologyobjects.ImproperType(psi_k=dihedral_fp.type.psi_k,
                                                              psi_eq=dihedral_fp.type.psi_eq,
                                                              list=topol_new.improper_types)

            # if type_to_assign not in topol_new.dihedral_types:
            topol_new.improper_types.append(type_to_assign)
            logger.info(
                f"      [b] New dihedral: {dihedral_fp.atom1.name:}-{dihedral_fp.atom2.name:}-{dihedral_fp.atom3.name:}-{dihedral_fp.atom4.name:} "
                f"{mapping_fp_to_new[dihedral_fp.atom1.idx]:} {mapping_fp_to_new[dihedral_fp.atom2.idx]:} "
                f"{mapping_fp_to_new[dihedral_fp.atom3.idx]:} {type_to_assign:}")

            atom1 = topol_new.atoms[mapping_fp_to_new[dihedral_fp.atom1.idx]]
            atom2 = topol_new.atoms[mapping_fp_to_new[dihedral_fp.atom2.idx]]
            atom3 = topol_new.atoms[mapping_fp_to_new[dihedral_fp.atom3.idx]]
            atom4 = topol_new.atoms[mapping_fp_to_new[dihedral_fp.atom4.idx]]

            topol_new.impropers.append(
                pmd.topologyobjects.Improper(atom1, atom2, atom3, atom4, type=type_to_assign))

    return topol_new

'''
def adjust_bonds(topol_new, topol_fp, mapping_fp_to_new):
    logger.info("   [ ] Adding new bonds to topology")

    for bond_fp in topol_fp.bonds:
        found = False

        for bond_new in topol_new.bonds:
            # if bond_fp.atom1.name=="ZN":
            #    print(bond_fp.atom1.idx , bond_new)

            # print(bond_new, bond_new.atom1.name)
            # Check if atoms in mapping (that they are not the additional atoms)
            if bond_fp.atom1.idx in mapping_fp_to_new and bond_fp.atom2.idx in mapping_fp_to_new:
                # print("It is in the mapping")

                if ((bond_new.atom1.idx == mapping_fp_to_new[bond_fp.atom1.idx] and bond_new.atom2.idx ==
                     mapping_fp_to_new[bond_fp.atom2.idx]) or
                        ((bond_new.atom1.idx == mapping_fp_to_new[bond_fp.atom2.idx] and bond_new.atom2.idx ==
                          mapping_fp_to_new[bond_fp.atom1.idx]))):
                    if (bond_fp.type != bond_new.type):
                        logger.info(f"      [o] Diffrent bond type {bond_new.type:} {bond_fp.type:}")
                    found = True
            else:
                found = True

        if not found:
            # type_to_assign= bond_fp.type
            # type_to_assign = pmd.topologyobjects.AngleType(angle_fp.type.k, angle_fp.type.theteq)
            # print("AAA", bond_fp.type.k, bond_fp.type.req, len(topol_new.bond_types))

            type_to_assign = pmd.topologyobjects.BondType(bond_fp.type.k, bond_fp.type.req,
                                                          list=topol_new.bond_types)

            # if type_to_assign not in topol_new.bond_types:
            topol_new.bond_types.append(type_to_assign)

            atom1 = topol_new.atoms[mapping_fp_to_new[bond_fp.atom1.idx]]
            atom2 = topol_new.atoms[mapping_fp_to_new[bond_fp.atom2.idx]]
            topol_new.bonds.append(pmd.topologyobjects.Bond(atom1, atom2, type=type_to_assign))

            logger.info(
                f"      [o] New bond: {mapping_fp_to_new[bond_fp.atom1.idx]:}  {type_to_assign:} {atom1.name:}-{atom2.name:} ({mapping_fp_to_new[bond_fp.atom2.idx]:})")

    return topol_new

'''
def adjust_pair_exclusions(topol_new, topol_fp, mapping_fp_to_new):
    logger.info("   [ ] Adding new pair exclusions to topology")
    for pair_fp in topol_fp.adjusts:
        found = False
        for pair_new in topol_new.adjusts:
            # Check if atoms in mapping (that they are not the additional atoms)
            if pair_fp.atom1.idx in mapping_fp_to_new and pair_fp.atom2.idx in mapping_fp_to_new:
                if ((pair_new.atom1.idx == mapping_fp_to_new[pair_fp.atom1.idx] and pair_new.atom2.idx ==
                     mapping_fp_to_new[pair_fp.atom2.idx]) or
                        ((pair_new.atom1.idx == mapping_fp_to_new[pair_fp.atom2.idx] and pair_new.atom2.idx ==
                          mapping_fp_to_new[pair_fp.atom1.idx]))):
                    found = True
            else:
                found = True

        if not found:
            atom1 = topol_new.atoms[mapping_fp_to_new[pair_fp.atom1.idx]]
            atom2 = topol_new.atoms[mapping_fp_to_new[pair_fp.atom2.idx]]

            type_to_assign = pmd.topologyobjects.NonbondedExceptionType(atom1.rmin + atom2.rmin,
                                                                        0.5 * (atom1.epsilon + atom2.epsilon))

            topol_new.adjust_types.append(type_to_assign)
            topol_new.adjusts.append(pmd.topologyobjects.NonbondedException(atom1, atom2))

            logger.info(
                f"      [o] New pair exclusion: {mapping_fp_to_new[pair_fp.atom1.idx]:}  {type_to_assign:} {atom1.name:}-{atom2.name:} ({mapping_fp_to_new[pair_fp.atom2.idx]:})")
    return topol_new


def add_1_4_metal_pairs(topol_new, metal_index):
    '''
    We add all the 1-4 interaction which are around the metal
    :param topol_new: (parmed topology) contains topology of the site
    :return: (parmed topology) modified topology
    '''
    orginal_pairs = []
    for pair in topol_new.adjusts:
        orginal_pairs.append((pair.atom1.idx, pair.atom2.idx))

    # We include also all the pairs
    G = nx.Graph([(bond.atom1.idx, bond.atom2.idx) for bond in topol_new.bonds])
    pairs_1_4 = []
    all_pairs_1_4 = list(nx.generators.ego_graph(G, metal_index, radius=3).nodes)
    G_cut = G.subgraph(all_pairs_1_4)

    for idx, node1 in enumerate(all_pairs_1_4):
        for node2 in all_pairs_1_4[idx+1:]:
            paths = list(nx.all_simple_paths(G_cut, node1, node2))
            if sum([1 for path in paths if len(set(path)) == 4 and metal_index in path ])>0:
                pairs_1_4.append([node1, node2])

    for pair in pairs_1_4:
        atom1 = topol_new.atoms[pair[0]]
        atom4 = topol_new.atoms[pair[1]]
        if (atom1.idx, atom4.idx) not in orginal_pairs or (atom4.idx, atom1.idx) not in orginal_pairs:
            type_to_assign = pmd.topologyobjects.NonbondedExceptionType(atom1.rmin + atom4.rmin,
                                                                        0.5 * (atom1.epsilon + atom4.epsilon))
            topol_new.adjust_types.append(type_to_assign)
            topol_new.adjusts.append(pmd.topologyobjects.NonbondedException(atom1, atom4))
    return topol_new




def copy_bonds(topol_new, bonds, metal_name):
    '''
    Copies bonds into topology

    :param topol_new:
    :param bonds:
    :return:
    '''
    orginal_bonds = []

    for bond in topol_new.bonds:
        orginal_bonds.append((bond.atom1.idx, bond.atom2.idx))

    for bond in bonds:
        find = [a for a, bond_o in enumerate(orginal_bonds) if (bond_o == bond or bond_o == bond[::-1])]
        if len(find) == 1:
            bond_new = topol_new.bonds[find[0]]

            type_to_assign = pmd.topologyobjects.BondType(bonds[bond][1], bonds[bond][0], list=topol_new.bond_types)
            # if type_to_assign not in self.topol_new.bond_types:
            topol_new.bond_types.append(type_to_assign)

            # bond_new.type = deepcopy(bond_fp.type)
            bond_new.type = type_to_assign  # bond_fp.type
        else:
            if strip_numbers_from_atom_name(
                    topol_new.atoms[bond[0]].name).title() == metal_name.title() or strip_numbers_from_atom_name(
                topol_new.atoms[bond[1]].name).title() == metal_name.title():
                type_to_assign = pmd.topologyobjects.BondType(bonds[bond][1], bonds[bond][0], list=topol_new.bond_types)
                topol_new.bond_types.append(type_to_assign)
                atom1 = topol_new.atoms[bond[0]]
                atom2 = topol_new.atoms[bond[1]]
                topol_new.bonds.append(pmd.topologyobjects.Bond(atom1, atom2, type=type_to_assign))
            else:
                print("this is wierd bond", topol_new.atoms[bond[0]].name, topol_new.atoms[bond[1]].name, bond,
                      bonds[bond])
    return topol_new


def copy_angles(topol_new, angles, metal_name):
    '''
    Copy angles into topology

    :param topol_new:
    :param angles:
    :return:
    '''
    orginal_angles = []

    for angle in topol_new.angles:
        orginal_angles.append((angle.atom1.idx, angle.atom2.idx, angle.atom3.idx))

    for angle in angles:
        find = [a for a, angle_o in enumerate(orginal_angles) if (angle_o == angle or angle_o == angle[::-1])]
        if len(find) == 1:
            angle_new = topol_new.angles[find[0]]
            type_to_assign = pmd.topologyobjects.AngleType(angles[angle][1], angles[angle][0],
                                                           list=topol_new.angle_types)

            topol_new.angle_types.append(type_to_assign)
            angle_new.type = type_to_assign
        else:
            if strip_numbers_from_atom_name(
                    topol_new.atoms[angle[0]].name).title() == metal_name.title() or strip_numbers_from_atom_name(
                topol_new.atoms[angle[1]].name).title() == metal_name.title() or strip_numbers_from_atom_name(
                topol_new.atoms[angle[2]].name).title() == metal_name.title():

                type_to_assign = pmd.topologyobjects.AngleType(angles[angle][1], angles[angle][0],
                                                               list=topol_new.angle_types)

                topol_new.angle_types.append(type_to_assign)

                atom1 = topol_new.atoms[angle[0]]
                atom2 = topol_new.atoms[angle[1]]
                atom3 = topol_new.atoms[angle[2]]

                topol_new.angles.append(pmd.topologyobjects.Angle(atom1, atom2, atom3, type=type_to_assign))
            else:
                print('nope', angle)

    return topol_new


def copy_dihedrals(topol_new, dihedrals, metal_name):
    '''
    Copy angles into topology

    :param topol_new:
    :param angles:
    :return:
    '''
    orginal_dihedrals = []

    for dihedral in topol_new.dihedrals:
        orginal_dihedrals.append((dihedral.atom1.idx, dihedral.atom2.idx, dihedral.atom3.idx, dihedral.atom4.idx))

    for dihedral in dihedrals:
        find = [a for a, dihedral_o in enumerate(orginal_dihedrals) if
                (dihedral_o == dihedral or dihedral_o == dihedral[::-1])]
        if len(find) == 1:
            raise Exception(
                "We suppose to add dihedral containing metal. The nonbondate template has them... it should not!")
            # This procedures are not done in reality
            # dihedral_new = topol_new.dihedrals[find[0]]
            # type_to_assign = pmd.topologyobjects.DihedralType(dihedrals[dihedral][1], dihedrals[dihedral][0],list=topol_new.angle_types)
            # topol_new.angle_types.append(type_to_assign)
            # dihedral_new.type = type_to_assign
        else:
            if strip_numbers_from_atom_name(
                    topol_new.atoms[dihedral[0]].name).title() == metal_name.title() or strip_numbers_from_atom_name(
                topol_new.atoms[dihedral[1]].name).title() == metal_name.title() or strip_numbers_from_atom_name(
                topol_new.atoms[dihedral[2]].name).title() == metal_name.title() or strip_numbers_from_atom_name(
                topol_new.atoms[dihedral[3]].name).title() == metal_name.title():

                # we do use periodicy=2 and use defualt values for scee=1.2 and scnb=2 (screening of electrostatics and noboned)
                type_to_assign = pmd.topologyobjects.DihedralType(phi_k=dihedrals[dihedral][1], per=1,
                                                                  phase=dihedrals[dihedral][0],
                                                                  list=topol_new.dihedral_types)

                topol_new.dihedral_types.append(type_to_assign)

                atom1 = topol_new.atoms[dihedral[0]]
                atom2 = topol_new.atoms[dihedral[1]]
                atom3 = topol_new.atoms[dihedral[2]]
                atom4 = topol_new.atoms[dihedral[3]]
                topol_new.dihedrals.append(
                    pmd.topologyobjects.Dihedral(atom1, atom2, atom3, atom4, type=type_to_assign))
            else:
                print('nope', dihedral)

    return topol_new


def copy_impropers(topol_new, dihedrals, metal_name):
    '''
    Copy angles into topology

    :param topol_new:
    :param angles:
    :return:
    '''
    orginal_dihedrals = []

    for dihedral in topol_new.impropers:
        orginal_dihedrals.append((dihedral.atom1.idx, dihedral.atom2.idx, dihedral.atom3.idx, dihedral.atom4.idx))

    for dihedral in dihedrals:
        find = [a for a, dihedral_o in enumerate(orginal_dihedrals) if
                (dihedral_o == dihedral or dihedral_o == dihedral[::-1])]
        if len(find) == 1:
            raise Exception(
                "We suppose to add dihedral containing metal. The nonbondate template has them... it should not!")
        else:
            if strip_numbers_from_atom_name(
                    topol_new.atoms[dihedral[0]].name).title() == metal_name.title() or strip_numbers_from_atom_name(
                topol_new.atoms[dihedral[1]].name).title() == metal_name.title() or strip_numbers_from_atom_name(
                topol_new.atoms[dihedral[2]].name).title() == metal_name.title() or strip_numbers_from_atom_name(
                topol_new.atoms[dihedral[3]].name).title() == metal_name.title():

                type_to_assign = pmd.topologyobjects.ImproperType(psi_k=dihedrals[dihedral][1],
                                                                  psi_eq=dihedrals[dihedral][0],
                                                                  list=topol_new.improper_types)

                topol_new.improper_types.append(type_to_assign)

                atom1 = topol_new.atoms[dihedral[0]]
                atom2 = topol_new.atoms[dihedral[1]]
                atom3 = topol_new.atoms[dihedral[2]]
                atom4 = topol_new.atoms[dihedral[3]]
                topol_new.impropers.append(
                    pmd.topologyobjects.Improper(atom1, atom2, atom3, atom4, type=type_to_assign))
            else:
                print('nope', dihedral)

    return topol_new


def copy_pair_exclusions(topol_new, pairs):
    '''
    Copy pairs exclusions to the topology
    :param topol_new: (parmed topology) contains topology of the site
    :param pairs: (list) list of pairs of atoms for pair exclusion
    :return: (parmed topology) modified topology
    '''
    orginal_pairs = []
    for pair in topol_new.adjusts:
        orginal_pairs.append((pair.atom1.idx, pair.atom2.idx))

    for pair in pairs:
        atom1 = topol_new.atoms[pair[0]]
        atom4 = topol_new.atoms[pair[1]]
        if (atom1.idx, atom4.idx) not in orginal_pairs or (atom4.idx, atom1.idx) not in orginal_pairs:
            type_to_assign = pmd.topologyobjects.NonbondedExceptionType(atom1.rmin + atom4.rmin,
                                                                        0.5 * (atom1.epsilon + atom4.epsilon))
            topol_new.adjust_types.append(type_to_assign)
            topol_new.adjusts.append(pmd.topologyobjects.NonbondedException(atom1, atom4))
    return topol_new



def update_pairs(topol):
    '''
    Removes the exclusion pairs which are 2 bonds away
    :param topol: (parmed topology)
    :return: (parmed topology) modified topology
    '''
    new_pairs = []
    angles = [[angle.atom1.idx, angle.atom2.idx, angle.atom3.idx] for angle in topol.angles]
    for pair in topol.adjusts:
        found = False
        for angle in angles:
            if pair.atom1.idx in angle and pair.atom2.idx in angle:
                found = True
                break

        # If the atoms belonging to pair are not part of angles then save this
        if found == False:
            new_pairs.append(pair)

    topol.adjusts = new_pairs
    return topol
