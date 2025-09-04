#!/usr/bin/env python

# Brandon L

# (C) 2022 Cadence Design Systems, Inc. (Cadence) 
# All rights reserved.
# TERMS FOR USE OF SAMPLE CODE The software below ("Sample Code") is
# provided to current licensees or subscribers of Cadence products or
# SaaS offerings (each a "Customer").
# Customer is hereby permitted to use, copy, and modify the Sample Code,
# subject to these terms. Cadence claims no rights to Customer's
# modifications. Modification of Sample Code is at Customer's sole and
# exclusive risk. Sample Code may require Customer to have a then
# current license or subscription to the applicable Cadence offering.
# THE SAMPLE CODE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
# EXPRESS OR IMPLIED.  OPENEYE DISCLAIMS ALL WARRANTIES, INCLUDING, BUT
# NOT LIMITED TO, WARRANTIES OF MERCHANTABILITY, FITNESS FOR A
# PARTICULAR PURPOSE AND NONINFRINGEMENT. In no event shall Cadence be
# liable for any damages or liability in connection with the Sample Code
# or its use.

# Requires openeye license to run
# To run: 
# conda env create -f openeye.yml
# conda activate openeye
# source /path/to/openeye/license

#############################################################################
# Align two compounds based on the maximum common substructure
#############################################################################
import sys
from openeye import oechem
import math

'''
def MCSAlign(refmol, fitmol, ofs):
    atomexpr = oechem.OEExprOpts_AtomicNumber | oechem.OEExprOpts_Aromaticity
    bondexpr = 0
    mcss = oechem.OEMCSSearch(oechem.OEMCSType_Exhaustive)
    mcss.Init(refmol, atomexpr, bondexpr)
    mcss.SetMCSFunc(oechem.OEMCSMaxBondsCompleteCycles())

    rmat = oechem.OEDoubleArray(9)
    trans = oechem.OEDoubleArray(3)
    unique = True
    overlay = True
    for match in mcss.Match(fitmol, unique):
        rms = oechem.OERMSD(mcss.GetPattern(), fitmol, match, overlay, rmat, trans)
        if rms < 0.0:
            oechem.OEThrow.Warning("RMS overlay failure")
            continue
        oechem.OERotate(fitmol, rmat)
        oechem.OETranslate(fitmol, trans)
        oechem.OEWriteMolecule(ofs, fitmol)
'''

import math
from openeye import oechem

_AMIDE_SMARTS = "[#16:6]-[#6:5]-[#6:4]-[#7:3]-[#6:2]=[#8:1]"

def _get_amide_atoms_from_match(match):
    """
    From an MCS or SMARTS match, return OEAtomBases for O, C, N, C (map 1,2,3,4).
    """
    atom_map = {}
    for mp in match.GetAtoms():
        mi = mp.pattern.GetMapIdx()
        atom_map[mi] = mp.target
    # Confirm all indices are present
    if not all(i in atom_map for i in (1,2,3,4)):
        raise RuntimeError("Match does not cover all four mapped amide atoms")
    return atom_map[1], atom_map[2], atom_map[3], atom_map[4]

def _get_amide_dihedral_from_match(mol, match):
    """
    Get the O–C–N–C torsion from a match (using mapped atoms).
    """
    atom_o, atom_c, atom_n, atom_r = _get_amide_atoms_from_match(match)
    return oechem.OEGetTorsion(mol, atom_o, atom_c, atom_n, atom_r), (atom_o, atom_c, atom_n, atom_r)

def _get_amide_match(mol):
    ss1 = oechem.OESubSearch("[#8:1]=[#6:2]-[#7:3]-[#6:4]-[#6:5]-[#16:6]")
    ss2 = oechem.OESubSearch("[#16:6]-[#6:5]-[#6:4]-[#7:3]-[#6:2]=[#8:1]")
    n1 = sum(1 for _ in ss1.Match(mol))
    n2 = sum(1 for _ in ss2.Match(mol))
    print(f"Matches in refmol: direction1={n1}, direction2={n2}")
    if n1:
        for match in ss1.Match(mol):
            return match
    if n2:
        for match in ss2.Match(mol):
            return match
    raise RuntimeError("No amide-next-to-sulfur fragment found")

def MCSAlign(refmol, fitmol, ofs):
    """
    Aligns fitmol onto refmol by MCS+RMSD, flipping amide to match refmol's cis geometry.
    """
    to_deg = lambda x: x * 180.0 / math.pi

    # 1) get reference torsion (using SMARTS match)
    ref_match = _get_amide_match(refmol)
    phi_ref, (atom_o_ref, atom_c_ref, atom_n_ref, atom_r_ref) = _get_amide_dihedral_from_match(refmol, ref_match)
    print('Reference dihedral (deg):', to_deg(phi_ref))

    # 2) set up MCS search
    atomexpr = oechem.OEExprOpts_AtomicNumber | oechem.OEExprOpts_Aromaticity
    mcss = oechem.OEMCSSearch(oechem.OEMCSType_Exhaustive)
    mcss.Init(refmol, atomexpr, 0)
    mcss.SetMCSFunc(oechem.OEMCSMaxBondsCompleteCycles())

    rmat = oechem.OEDoubleArray(9)
    trans = oechem.OEDoubleArray(3)
    unique = True
    overlay = True

    # 3) loop over matches
    for match in mcss.Match(fitmol, unique):
        try:
            # Get dihedral for these exact matched atoms
            phi_fit, (atom_o, atom_c, atom_n, atom_r) = _get_amide_dihedral_from_match(fitmol, match)
        except RuntimeError:
            continue  # no motif to measure

        print('fitmol dihedral (deg):', to_deg(phi_fit))

        # 4) Flip amide to match reference (cis) if needed
        delta = phi_ref - phi_fit
        if abs(to_deg(delta)) > 1.0:  # skip if already basically matching, or always rotate
            # Rotate the N–C bond (atom_n–atom_r) by delta
            oechem.OERotateBond(fitmol, atom_n, atom_r, delta, oechem.OERotateBond_All)

            # Optional: remeasure dihedral to confirm
            phi_fit_new = oechem.OEGetTorsion(fitmol, atom_o, atom_c, atom_n, atom_r)
            print('After rotation, fitmol dihedral (deg):', to_deg(phi_fit_new))

        # 5) RMSD overlay
        rms = oechem.OERMSD(mcss.GetPattern(),
                            fitmol, match,
                            overlay, rmat, trans)
        if rms < 0.0:
            oechem.OEThrow.Warning("RMS overlay failure")
            continue

        oechem.OERotate(fitmol, rmat)
        oechem.OETranslate(fitmol, trans)
        oechem.OEWriteMolecule(ofs, fitmol)

def main(argv=[__name__]):
    if len(argv) != 4:
        oechem.OEThrow.Usage("%s <refmol> <fitmol> <outfile>" % argv[0])

    reffs = oechem.oemolistream()
    if not reffs.open(argv[1]):
        oechem.OEThrow.Fatal("Unable to open %s for reading" % argv[1])
    if not oechem.OEIs3DFormat(reffs.GetFormat()):
        oechem.OEThrow.Fatal("Invalid input format: need 3D coordinates")
    refmol = oechem.OEGraphMol()
    if not oechem.OEReadMolecule(reffs, refmol):
        oechem.OEThrow.Fatal("Unable to read molecule in %s" % argv[1])
    if not refmol.GetDimension() == 3:
        oechem.OEThrow.Fatal("%s doesn't have 3D coordinates" % refmol.GetTitle())

    fitfs = oechem.oemolistream()
    if not fitfs.open(argv[2]):
        oechem.OEThrow.Fatal("Unable to open %s for reading" % argv[2])
    if not oechem.OEIs3DFormat(fitfs.GetFormat()):
        oechem.OEThrow.Fatal("Invalid input format: need 3D coordinates")

    ofs = oechem.oemolostream()
    if not ofs.open(argv[3]):
        oechem.OEThrow.Fatal("Unable to open %s for writing" % argv[3])
    if not oechem.OEIs3DFormat(ofs.GetFormat()):
        oechem.OEThrow.Fatal("Invalid output format: need 3D coordinates")

    #oechem.OEWriteConstMolecule(ofs, refmol)
    oechem.OESuppressHydrogens(refmol)

    for fitmol in fitfs.GetOEGraphMols():
        if not fitmol.GetDimension() == 3:
            oechem.OEThrow.Warning("%s doesn't have 3D coordinates" % fitmol.GetTitle())
            continue
        MCSAlign(refmol, fitmol, ofs)


if __name__ == "__main__":
    sys.exit(main(sys.argv))
