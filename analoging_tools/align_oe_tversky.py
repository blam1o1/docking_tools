#!/usr/bin/env python

# Brandon L
#
# Script to calculate the Tverksy fit score between reference(reffile) and queries(fitfile).
# Current parameters for tversky calcuation = alpha:0.95 beta:0.05
# Read more at: https://docs.eyesopen.com/applications/rocs/theory/measure_similarity.html 
#
# Requires openeye license to run
# To run: 
# conda env create -f openeye.yml
# conda activate openeye
# source /path/to/openeye/license


import sys
import csv
import argparse
from openeye import oechem, oeshape

def main(argv=None):
    parser = argparse.ArgumentParser(
        description='Calculate shape and color Tversky scores and save passing entries to a CSV file.'
    )
    parser.add_argument('reffile', help='Path to reference molecule file')
    parser.add_argument('fitfile', help='Path to fit molecules file')
    parser.add_argument('--cutoff', type=float, default=0.0,
                        help='Tversky score cutoff; only entries with score ≥ cutoff are saved')
    parser.add_argument('--output', default='results.csv',
                        help='Output CSV filename')
    args = parser.parse_args(argv[1:] if argv is not None else None)

    reffs = oechem.oemolistream()
    if not reffs.open(args.reffile):
        oechem.OEThrow.Fatal(f"Can't open reference file {args.reffile}")
    fitfs = oechem.oemolistream()
    if not fitfs.open(args.fitfile):
        oechem.OEThrow.Fatal(f"Can't open fit file {args.fitfile}")

    refmol = oechem.OEGraphMol()
    if not oechem.OEReadMolecule(reffs, refmol):
        oechem.OEThrow.Fatal('Failed to read reference molecule')
    prep = oeshape.OEOverlapPrep()
    prep.Prep(refmol)
    shapeFunc = oeshape.OEExactShapeFunc()
    shapeFunc.SetupRef(refmol)
    results = oeshape.OEOverlapResults()

    entries = []
    for fitmol in fitfs.GetOEGraphMols():
        title = fitmol.GetTitle() or '<untitled>'
        prep.Prep(fitmol)
        shapeFunc.Overlap(fitmol, results)
        score = results.GetTversky(alpha=0.95, beta=0.05) + results.GetColorTversky(alpha=0.95, beta=0.05)
        num_atoms = fitmol.NumAtoms()
        if score >= args.cutoff:
            entries.append((title, score, num_atoms))

    entries.sort(key=lambda x: x[1], reverse=True)

    with open(args.output, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(['Title', 'Score', 'NumAtoms'])
        for title, score, num_atoms in entries:
            writer.writerow([title, f"{score:.6f}", num_atoms])

    print(f"Results saved to '{args.output}' (cutoff={args.cutoff}), sorted by score descending")

if __name__ == '__main__':
    sys.exit(main(sys.argv))
