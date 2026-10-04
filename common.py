"""Shared I/O, validation and chemistry helpers. Units are explicit at boundaries."""
import argparse
import csv
import json
import math
from pathlib import Path

R_KJ = 0.00831446261815324
EV_KJ = 96.4853321233
HARTREE_KJ = 2625.49963948

def positive(value):
    x = float(value)
    if not math.isfinite(x) or x <= 0:
        raise argparse.ArgumentTypeError('Expected a finite positive number')
    return x

def finite(value):
    x = float(value)
    if not math.isfinite(x):
        raise argparse.ArgumentTypeError('Expected a finite number')
    return x

def parser(description):
    p = argparse.ArgumentParser(description=description)
    p.add_argument('--out', default='results', help='Output directory (created if needed)')
    return p

def output_dir(args):
    p = Path(args.out)
    p.mkdir(parents=True, exist_ok=True)
    return p

def dump(path, obj):
    Path(path).write_text(json.dumps(obj, indent=2, allow_nan=False) + '\n')

def read_json(path):
    return json.loads(Path(path).read_text())

def rows(path):
    with open(path, newline='') as f:
        return list(csv.DictReader(f))

def table(path, data, fields=None):
    data = list(data)
    if fields is None:
        if not data:
            raise ValueError('No rows to write')
        fields = list(data[0])
    with open(path, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(data)

def molecule(smiles, neutral=False):
    from rdkit import Chem
    m = Chem.MolFromSmiles(smiles)
    if m is None or len(Chem.GetMolFrags(m)) != 1:
        raise ValueError('Supply a valid, single-component SMILES; split salts into components')
    if neutral and Chem.GetFormalCharge(m) != 0:
        raise ValueError('This model requires the neutral molecule')
    return m

def descriptors(smiles):
    from rdkit.Chem import Descriptors, Crippen, Lipinski, rdMolDescriptors
    m = molecule(smiles)
    return dict(smiles=smiles, mw_g_mol=Descriptors.MolWt(m),
                logP=Crippen.MolLogP(m), tpsa_A2=rdMolDescriptors.CalcTPSA(m),
                hbd=Lipinski.NumHDonors(m), hba=Lipinski.NumHAcceptors(m),
                rotatable_bonds=Lipinski.NumRotatableBonds(m),
                aromatic_rings=rdMolDescriptors.CalcNumAromaticRings(m))

def read_xyz(path):
    from rdkit import Chem
    lines = Path(path).read_text().splitlines()
    try:
        n = int(lines[0])
        data = lines[2:]
        if n < 1 or len(data) != n:
            raise ValueError('XYZ atom count does not match file')
        pt = Chem.GetPeriodicTable()
        for line in data:
            fields = line.split()
            if len(fields) != 4 or pt.GetAtomicNumber(fields[0]) < 1:
                raise ValueError('Invalid XYZ atom record')
            if not all(math.isfinite(float(x)) for x in fields[1:]):
                raise ValueError('Nonfinite XYZ coordinate')
    except (IndexError, TypeError) as e:
        raise ValueError('Invalid XYZ file') from e
    return '\n'.join(data)

def quantum_args(p):
    p.add_argument('xyz', help='Explicit-H geometry in Angstrom')
    p.add_argument('--charge', type=int, default=0)
    p.add_argument('--multiplicity', type=int, default=1)
    p.add_argument('--method', default='b3lyp-d3bj/def2-svp')
    p.add_argument('--threads', type=int, default=4)
    p.add_argument('--memory', default='4 GB')

def psi4_setup(args):
    import psi4
    from rdkit import Chem
    body = read_xyz(args.xyz)
    electrons = sum(Chem.GetPeriodicTable().GetAtomicNumber(l.split()[0])
                    for l in body.splitlines()) - args.charge
    if args.multiplicity < 1 or args.threads < 1 or electrons < args.multiplicity - 1 or (electrons - (args.multiplicity-1)) % 2:
        raise ValueError('Invalid electron count, spin multiplicity or thread count')
    out = output_dir(args)
    psi4.set_num_threads(args.threads)
    psi4.set_memory(args.memory)
    psi4.core.set_output_file(str(out / 'psi4.log'), False)
    mol = psi4.geometry(f'{args.charge} {args.multiplicity}\n{body}\nunits angstrom\nsymmetry c1\nno_reorient\nno_com')
    psi4.set_options({'reference': 'rhf' if args.multiplicity == 1 else 'uhf',
                      'scf_type': 'df', 'e_convergence': 1e-8, 'd_convergence': 1e-8})
    return psi4, mol, out

def broaden(peaks, intensities, lower, upper, width, points=4001):
    import numpy as np
    if width <= 0 or upper <= lower or points < 2:
        raise ValueError('Invalid spectrum bounds or width')
    x = np.linspace(lower, upper, points)
    y = np.zeros_like(x)
    for peak, amplitude in zip(peaks, intensities, strict=True):
        y += amplitude * (width / 2 / np.pi) / ((x - peak)**2 + (width/2)**2)
    return x, y

def spectrum_csv(path, x, y, axis):
    table(path, ({axis: float(a), 'relative_intensity': float(b)} for a,b in zip(x,y,strict=True)))

def activity(value):
    x = positive(value)
    if x > 1:
        raise argparse.ArgumentTypeError('Activity must be in (0,1] for this workflow')
    return x
