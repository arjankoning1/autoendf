#!/usr/bin/env python3
"""Compare neutron and neutron-induced photon leakage from evaluated ENDF libraries."""
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
import xml.etree.ElementTree as ET
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import openmc
import openmc.data

WORLD_LIBRARIES = ('endfb8.1', 'jeff4.0', 'jendl5.0', 'tendl.2025')
PLOT_STYLES = {
    'new_tendl': {'label': 'NEW', 'color': 'black'},
    'endfb8.1': {'label': 'ENDF/B-VIII.1', 'color': 'blue'},
    'jeff4.0': {'label': 'JEFF-4.0', 'color': 'red'},
    'jendl5.0': {'label': 'JENDL-5.0', 'color': 'green'},
    'tendl.2025': {'label': 'TENDL-2025', 'color': 'orange', 'linestyle': '--'},
}


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


def find_photon_hdf5(index_path, element):
    """Return elemental photoatomic HDF5 from an existing OpenMC XML index."""
    index_path = Path(index_path).expanduser().resolve()
    if not index_path.is_file():
        raise FileNotFoundError(f'Photon cross-section index not found: {index_path}')
    root = ET.parse(index_path).getroot()
    directory = root.findtext('directory')
    base = Path(directory).expanduser() if directory else index_path.parent
    if not base.is_absolute():
        base = (index_path.parent / base).resolve()
    for entry in root.findall('library'):
        if entry.get('type') != 'photon':
            continue
        if element not in entry.get('materials', '').split():
            continue
        path = Path(entry.attrib['path'])
        path = path if path.is_absolute() else base / path
        path = path.resolve()
        if not path.is_file():
            raise FileNotFoundError(f'Photon HDF5 file referenced by {index_path} is missing: {path}')
        return path
    raise ValueError(f'No photon HDF5 entry for element {element} in {index_path}')


def transport_index(xs, photon_hdf5, output):
    """Keep each evaluation's own neutron data, with common photoatomic data."""
    if photon_hdf5 is None:
        return xs
    target = output / 'cross_sections_transport.xml'
    # Generate index using HDF5 metadata rather than editing paths by hand.
    neutron_lib = openmc.data.DataLibrary.from_xml(str(xs))
    neutron_files = [Path(item['path']) for item in neutron_lib if item['type'] == 'neutron']
    if len(neutron_files) != 1:
        raise ValueError(f'Expected one neutron entry in {xs}, found {len(neutron_files)}')
    neutron_h5 = neutron_files[0]
    if not neutron_h5.is_absolute():
        neutron_h5 = (xs.parent / neutron_h5).resolve()
    library = openmc.data.DataLibrary()
    library.register_file(str(neutron_h5))
    library.register_file(str(photon_hdf5))
    library.export_to_xml(str(target))
    return target


def run_case(label, xs, output, nuclide, args, neutron_edges, photon_edges, photon_hdf5):
    output.mkdir(parents=True, exist_ok=True)
    material = openmc.Material(name=nuclide)
    material.add_nuclide(nuclide, 1.)
    material.set_density('g/cm3', args.density)
    material.temperature = args.temperature
    materials = openmc.Materials([material])
    materials.cross_sections = str(transport_index(xs, photon_hdf5, output))
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
    settings.photon_transport = photon_hdf5 is not None

    tallies = openmc.Tallies()
    for particle, edges in [('neutron', neutron_edges), ('photon', photon_edges)]:
        if particle == 'photon' and photon_hdf5 is None:
            continue
        tally = openmc.Tally(name=f'{particle}_leakage')
        tally.filters = [openmc.SurfaceFilter(surface),
                         openmc.ParticleFilter([particle]),
                         openmc.EnergyFilter(edges)]
        tally.scores = ['current']
        tallies.append(tally)
    model = openmc.Model(geometry=geometry, materials=materials,
                         settings=settings, tallies=tallies)
    print(f'[{label}] OpenMC transport ({args.particles * args.batches} histories, '
          f'photons={settings.photon_transport})', flush=True)
    statepoint = model.run(cwd=str(output), threads=args.threads)
    spectra = {}
    with openmc.StatePoint(statepoint) as sp:
        for particle, edges in [('neutron', neutron_edges), ('photon', photon_edges)]:
            if particle == 'photon' and photon_hdf5 is None:
                continue
            result = sp.get_tally(name=f'{particle}_leakage')
            values = result.mean.ravel().copy()
            errors = result.std_dev.ravel().copy()
            widths = np.diff(edges) / 1e6
            energies = (edges[:-1] + edges[1:]) / 2e6
            # Vacuum outer surface: net outward current = escaping particles.
            spectrum = values / widths
            error = errors / widths
            np.savetxt(output / f'{particle}_spectrum.dat',
                       np.column_stack([energies, spectrum, error]),
                       header='Energy_MeV Leakage_MeV-1_per_source Std_dev')
            spectra[particle] = (spectrum, error)
    return spectra


def make_plot(outdir, particle, nuclide, energies, results, new_label):
    fig, (ax, ratio_ax) = plt.subplots(2, 1, figsize=(10, 8), sharex=True,
                                       gridspec_kw={'height_ratios': [3, 1]})
    ref, _ = results[new_label]
    for label, (values, errs) in results.items():
        style = PLOT_STYLES.get(label, {'label': label})
        ax.step(energies, np.where(values > 0, values, np.nan), where='mid',
                linewidth=2 if label == new_label else 1.1, **style)
        if label != new_label:
            ratio = np.divide(values, ref, out=np.full_like(values, np.nan),
                              where=(ref > 0) & (values > 0))
            ratio_ax.step(energies, ratio, where='mid', **style)
    ax.set(yscale='log', ylabel='Leakage / source neutron / MeV',
           title=f'OpenMC {particle} leakage: 14 MeV neutron source in {nuclide} sphere of 10 cm radius')
    ax.legend(fontsize=9)
    ax.grid(alpha=.3)
    ratio_ax.axhline(1., color='black', linestyle='--', linewidth=.8)
    ratio_ax.set(xlabel=f'{particle.capitalize()} energy (MeV)', ylabel='Library / NEW')
    ratio_ax.grid(alpha=.3)
    fig.tight_layout()
    stem = f'{particle[0]}-{nuclide}_comparison'
    plot_path = outdir / f'{stem}.png'
    fig.savefig(plot_path, dpi=180)
    plt.close(fig)
    with (outdir / f'{stem}.csv').open('w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['Energy_MeV'] + [part for label in results
                         for part in (f'{label}_leakage', f'{label}_stderr')])
        for i, energy in enumerate(energies):
            writer.writerow([energy] + [part for spec, err in results.values()
                           for part in (spec[i], err[i])])
    return plot_path


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
    parser.add_argument('-photon-xs', default=os.environ.get('OPENMC_CROSS_SECTIONS'),
                        help='Cross sections XML with photoatomic HDF5 for element (required unless -no-photons)')
    parser.add_argument('-no-photons', action='store_true',
                        help='Run neutron leakage only, without coupled neutron-photon transport')
    parser.add_argument('-energy', type=float, default=14.)
    parser.add_argument('-radius', type=float, default=10.)
    parser.add_argument('-density', type=float, default=7.874)
    parser.add_argument('-temperature', type=float, default=293.6)
    parser.add_argument('-particles', type=int, default=50000)
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
    photon_hdf5 = None
    if not args.no_photons:
        if not args.photon_xs:
            parser.error('Photon transport requires -photon-xs /path/to/cross_sections.xml '
                         '(or OPENMC_CROSS_SECTIONS). Alternatively use -no-photons.')
        photon_hdf5 = find_photon_hdf5(args.photon_xs, element)
        print(f'Using common {element} photoatomic data: {photon_hdf5}', flush=True)

    neutron_edges = np.linspace(0, max(15., args.energy + 1.) * 1e6, 151)
    photon_edges = np.linspace(0, max(15., args.energy + 1.) * 1e6, 151)
    centers = {
        'neutron': (neutron_edges[:-1] + neutron_edges[1:]) / 2e6,
        'photon': (photon_edges[:-1] + photon_edges[1:]) / 2e6,
    }
    results = {'neutron': {}, 'photon': {}}
    manifest = {'nuclide': nuclide, 'source_energy_MeV': args.energy,
                'density_g_cm3': args.density, 'radius_cm': args.radius,
                'temperature_K': args.temperature, 'libraries_dir': str(root),
                'photoatomic_hdf5': str(photon_hdf5) if photon_hdf5 else None,
                'observable': 'outward surface current per source neutron per MeV',
                'runs': {}, 'skipped_missing': missing, 'failed': {}}
    for label, path in paths:
        try:
            xs, digest = process_endf(label, path, outdir / 'cache', nuclide,
                                      njoy, args.temperature)
            run_dir = outdir / 'runs' / label
            spectra = run_case(label, xs, run_dir, nuclide, args,
                               neutron_edges, photon_edges, photon_hdf5)
            for particle, (values, errors) in spectra.items():
                results[particle][label] = (values, errors)
            manifest['runs'][label] = {
                'source_endf': str(path), 'sha256': digest,
                'cross_sections_xml': str(xs),
                'spectra': {particle: str(run_dir / f'{particle}_spectrum.dat')
                            for particle in spectra}}
        except Exception as exc:
            manifest['failed'][label] = str(exc)
            print(f'[ERROR] {label}: {exc}', file=sys.stderr, flush=True)
            if label == 'new_tendl':
                (outdir / 'manifest.json').write_text(json.dumps(manifest, indent=2) + '\n')
                raise
    for particle in ('neutron', 'photon'):
        if results[particle]:
            plot_path = make_plot(outdir, particle, nuclide, centers[particle],
                                  results[particle], 'new_tendl')
            print(f'Plot: {plot_path}', flush=True)
    (outdir / 'manifest.json').write_text(json.dumps(manifest, indent=2) + '\n')
    print(f'Finished: {len(results["neutron"])} successful evaluations, '
          f'{len(missing)} missing, {len(manifest["failed"])} failed.', flush=True)
    return 0 if not manifest['failed'] else 2


if __name__ == '__main__':
    raise SystemExit(main())
