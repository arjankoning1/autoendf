#!/usr/bin/env python3
"""Compare new neutron ENDF evaluation with installed AUTOTALYS world libraries."""
import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import sys
import traceback
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import openmc
import openmc.data

WORLD_LIBRARIES = ('endfb8.1', 'jeff4.0', 'jendl5.0', 'cendl3.2', 'tendl.2025')


def checksum(path):
    h = hashlib.sha256()
    with path.open('rb') as inp:
        for chunk in iter(lambda: inp.read(1024 * 1024), b''):
            h.update(chunk)
    return h.hexdigest()


def find_libraries(root, element, mass, selection):
    isotope = f'{element}{mass:03d}'
    projnuc = f'n-{isotope}'
    labels = WORLD_LIBRARIES if selection == 'all' else tuple(x.strip() for x in selection.split(',') if x.strip())
    for label in labels:
        if not re.fullmatch(r'[A-Za-z0-9_.-]+', label):
            raise ValueError(f'Invalid library name: {label}')
        path = root / 'n' / isotope / label / 'files' / f'{projnuc}.{label}'
        yield label, path


def process_endf(label, src, cache, nuclide, njoy, temperature):
    digest = checksum(src)
    # Settings are part of the key; preserve each processing provenance manifest.
    key = hashlib.sha256(f'{digest}:{temperature}:openmc-{openmc.__version__}:njoy-{njoy}'.encode()).hexdigest()[:16]
    directory = cache / label / key
    directory.mkdir(parents=True, exist_ok=True)
    h5 = directory / f'{nuclide}.h5'
    xs = directory / 'cross_sections.xml'
    manifest = directory / 'processing.json'
    if not (h5.is_file() and xs.is_file() and manifest.is_file()):
        print(f'[{label}] NJOY processing: {src}', flush=True)
        evaluated = openmc.data.IncidentNeutron.from_njoy(
            str(src), temperatures=[temperature], njoy_exec=njoy)
        if evaluated.name != nuclide:
            raise ValueError(f'ENDF target {evaluated.name} is not requested {nuclide}')
        evaluated.export_to_hdf5(str(h5))
        lib = openmc.data.DataLibrary()
        lib.register_file(str(h5))
        lib.export_to_xml(str(xs))
        manifest.write_text(json.dumps({
            'label': label, 'source_endf': str(src), 'sha256': digest,
            'njoy_executable': njoy, 'temperature_K': temperature,
            'openmc_version': openmc.__version__, 'hdf5': str(h5),
        }, indent=2) + '\n')
    else:
        print(f'[{label}] Reusing processed HDF5: {h5}', flush=True)
    return xs, digest


def run_case(label, xs, output, nuclide, args, edges):
    output.mkdir(parents=True, exist_ok=True)
    material = openmc.Material(name=nuclide)
    material.add_nuclide(nuclide, 1.)
    material.set_density('g/cm3', args.density)
    material.temperature = args.temperature
    materials = openmc.Materials([material])
    materials.cross_sections = str(xs)
    surface = openmc.Sphere(r=args.radius, boundary_type='vacuum')
    cell = openmc.Cell(fill=material, region=-surface)
    geometry = openmc.Geometry([cell])
    source = openmc.IndependentSource(
        space=openmc.stats.Point((0., 0., 0.)),
        angle=openmc.stats.Isotropic(),
        energy=openmc.stats.Discrete([args.energy * 1e6], [1.]),
        particle='neutron')
    settings = openmc.Settings()
    settings.run_mode = 'fixed source'
    settings.source = source
    settings.batches = args.batches
    settings.particles = args.particles
    settings.seed = args.seed
    tally = openmc.Tally(name='neutron_flux')
    tally.filters = [openmc.CellFilter(cell), openmc.ParticleFilter(['neutron']), openmc.EnergyFilter(edges)]
    tally.scores = ['flux']
    model = openmc.Model(geometry=geometry, materials=materials,
                         settings=settings, tallies=openmc.Tallies([tally]))
    print(f'[{label}] OpenMC transport ({args.particles * args.batches} histories)', flush=True)
    statepoint = model.run(cwd=str(output), threads=args.threads)
    with openmc.StatePoint(statepoint) as sp:
        result = sp.get_tally(name='neutron_flux')
        v = result.mean.ravel().copy()
        e = result.std_dev.ravel().copy()
    volume = 4 * np.pi * args.radius ** 3 / 3
    widths = np.diff(edges) / 1e6
    energies = (edges[:-1] + edges[1:]) / 2e6
    spectrum = v / volume / widths
    error = e / volume / widths
    np.savetxt(output / 'spectrum.dat', np.column_stack([energies, spectrum, error]),
               header='Energy_MeV Flux_cm-2_MeV-1_per_source Std_dev')
    return spectrum, error


def make_plots(outdir, energies, results, new_label):
    fig, (ax, ratio_ax) = plt.subplots(2, 1, figsize=(10, 8), sharex=True,
                                       gridspec_kw={'height_ratios': [3, 1]})
    ref, _ = results[new_label]
    for label, (values, errs) in results.items():
        ax.step(energies, np.where(values > 0, values, np.nan), where='mid', label=label,
                linewidth=2 if label == new_label else 1.1)
        if label != new_label:
            ratio = np.divide(values, ref, out=np.full_like(values, np.nan), where=ref > 0)
            ratio_ax.step(energies, ratio, where='mid', label=label)
    ax.set(yscale='log', ylabel='Volume-averaged flux / source\n(cm$^{-2}$ MeV$^{-1}$)',
           title='OpenMC neutron spectra: same geometry and source, different evaluations')
    ax.legend(fontsize=9)
    ax.grid(alpha=.3)
    ratio_ax.axhline(1., color='black', linestyle='--', linewidth=.8)
    ratio_ax.set(xlabel='Neutron energy (MeV)', ylabel='Library / new TENDL')
    ratio_ax.grid(alpha=.3)
    fig.tight_layout()
    fig.savefig(outdir / 'comparison.png', dpi=180)
    plt.close(fig)
    with (outdir / 'comparison.csv').open('w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['Energy_MeV'] + [part for label in results for part in (f'{label}_flux', f'{label}_stderr')])
        for i, e in enumerate(energies):
            writer.writerow([e] + [part for spec, err in results.values() for part in (spec[i], err[i])])


def main():
    parser = argparse.ArgumentParser(description='Compare own ENDF neutron evaluation against installed AUTOTALYS world libraries using OpenMC')
    parser.add_argument('-file', required=True)
    parser.add_argument('-element', required=True)
    parser.add_argument('-mass', type=int, required=True)
    parser.add_argument('-outdir', default=None)
    parser.add_argument('-libraries', default=os.environ.get('LIBRARIES_DIR', str(Path.home() / 'libraries')),
                        help='LIBRARIES_DIR root; defaults to environment or ~/libraries')
    parser.add_argument('-worldlibs', default='all', help='all, none, or comma-separated names')
    parser.add_argument('-njoybin', default='njoy')
    parser.add_argument('-energy', type=float, default=14.)
    parser.add_argument('-radius', type=float, default=10.)
    parser.add_argument('-density', type=float, default=7.874)
    parser.add_argument('-temperature', type=float, default=293.6)
    parser.add_argument('-particles', type=int, default=5000)
    parser.add_argument('-batches', type=int, default=20)
    parser.add_argument('-threads', type=int, default=4)
    parser.add_argument('-seed', type=int, default=12345)
    args = parser.parse_args()
    if args.mass < 1 or min(args.energy, args.radius, args.density, args.temperature) <= 0 or min(args.particles, args.batches, args.threads) < 1:
        parser.error('Mass, energy, radius, density, temperature, particles, batches and threads must be positive')
    element = args.element.capitalize()
    nuclide = f'{element}{args.mass}'
    src = Path(args.file).expanduser().resolve(strict=True)
    outdir = (Path(args.outdir).expanduser() if args.outdir else Path.cwd() / ('openmc_' + nuclide)).resolve()
    outdir.mkdir(parents=True, exist_ok=True)
    root = Path(args.libraries).expanduser().resolve()
    njoy = shutil.which(args.njoybin)
    if njoy is None:
        parser.error(f'NJOY executable not found: {args.njoybin}. Supply -njoybin /full/path/to/njoy')
    paths = [('new_tendl', src)]
    missing = []
    if args.worldlibs != 'none':
        for label, path in find_libraries(root, element, args.mass, args.worldlibs):
            if path.is_file():
                paths.append((label, path.resolve()))
            else:
                missing.append(f'{label}: {path}')
    for x in missing:
        print(f'[SKIP] Missing {x}', flush=True)
    edges = np.linspace(0, max(15., args.energy + 1.) * 1e6, 151)
    centers = (edges[:-1] + edges[1:]) / 2e6
    results = {}
    manifest = {'nuclide': nuclide, 'source_energy_MeV': args.energy,
                'density_g_cm3': args.density, 'radius_cm': args.radius,
                'temperature_K': args.temperature, 'libraries_dir': str(root),
                'runs': {}, 'skipped_missing': missing, 'failed': {}}
    for label, path in paths:
        try:
            xs, digest = process_endf(label, path, outdir / 'cache', nuclide, njoy, args.temperature)
            run_dir = outdir / 'runs' / label
            values, errors = run_case(label, xs, run_dir, nuclide, args, edges)
            results[label] = (values, errors)
            manifest['runs'][label] = {'source_endf': str(path), 'sha256': digest,
                                      'cross_sections_xml': str(xs), 'spectrum': str(run_dir / 'spectrum.dat')}
        except Exception as exc:
            manifest['failed'][label] = str(exc)
            print(f'[ERROR] {label}: {exc}', file=sys.stderr, flush=True)
            if label == 'new_tendl':
                (outdir / 'manifest.json').write_text(json.dumps(manifest, indent=2) + '\n')
                raise
    make_plots(outdir, centers, results, 'new_tendl')
    (outdir / 'manifest.json').write_text(json.dumps(manifest, indent=2) + '\n')
    print(f'Finished: {len(results)} successful evaluations, {len(missing)} missing, '
          f'{len(manifest["failed"])} failed. Plot: {outdir / "comparison.png"}', flush=True)
    return 0 if not manifest['failed'] else 2


if __name__ == '__main__':
    raise SystemExit(main())
