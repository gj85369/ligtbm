#!/usr/bin/env python3

import click
import json
import logging
import numpy as np
import multiprocessing
import pybel
import os
import signal
import subprocess
import logging.config


ATLAS_TIMEOUT = 2400  # seconds
OUT_DIR = 'parameterized'

FILE_DIR = os.path.dirname(os.path.abspath(__file__))
logging.config.fileConfig(os.path.join(FILE_DIR, 'logging.conf'), disable_existing_loggers=False)
LOG_USR = logging.getLogger('user')
LOG_DEV = logging.getLogger('dev')


def run_atlas(input_mol2, output_json, options):
    # The `timeout` parameter of subprocess.Popen (and its siblings) only kills the called process, not its children.
    # To kill children (they're the ones that hang) we have to do the process group magick
    args = [
        options['ATLAS_BIN'],
        '--dont-minimize',
        '--pH', '7.4',
        input_mol2, output_json]
    LOG_DEV.info('Starting {}'.format(' '.join(args)))
    proc = subprocess.Popen(args, preexec_fn=os.setpgrp)
    try:
        proc.wait(ATLAS_TIMEOUT)
    except subprocess.TimeoutExpired:
        os.killpg(os.getpgid(proc.pid), signal.SIGTERM)
        LOG_DEV.warning('ATLAS for input file {} timed out'.format(input_mol2))
        raise RuntimeError('Unable to run ATLAS, timed out (max {} seconds)'.format(ATLAS_TIMEOUT))
    if proc.returncode != 0:
        LOG_DEV.warning('ATLAS for input file {} failed'.format(input_mol2))
        raise RuntimeError('Unable to run ATLAS, return code {}'.format(proc.returncode))
    if not os.path.exists(output_json):
        LOG_DEV.warning('ATLAS for input file {} did not produce output'.format(input_mol2))
        raise RuntimeError('Unable to run ATLAS, output file does not exist')


def run_pdbprep(input_pdb, output_prefix, options):
    args = [
        options['SBLU_BIN'], 'pdb', 'prep',
        input_pdb,
        '--no-minimize',
        '--prm', options['CHARMM_PRM'],
        '--rtf', options['CHARMM_RTF'],
        '--out-prefix', output_prefix]
    LOG_DEV.info('Starting {}'.format(' '.join(args)))
    subprocess.check_call(args)


def read_mol2(mol2_filename):
    mol = next(pybel.readfile('mol', mol2_filename))
    assert mol.dim == 3  # We managed to get coordinates from mol2
    coords = [at.coords for at in mol.atoms]
    return np.array(coords)


def read_json(json_filename):
    with open(json_filename) as f:
        json_data = json.load(f)['atoms']
    json_names = [i['name'] for i in json_data]
    json_crds = [(i['x'], i['y'], i['z']) for i in json_data]
    return json_names, np.array(json_crds)


def map_mol2_to_json_by_coords(input_json_filename, input_mol2_filename):
    def crd_eq(a, b):
        DELTA_MAX = 0.01
        a = np.array(a)
        b = np.array(b)
        d = np.sum((a-b)**2) ** 0.5
        return d < DELTA_MAX

    json_names, json_crds = read_json(input_json_filename)
    assert len(json_crds) == len(json_names)  # A little sanity check
    mol2_crds = read_mol2(input_mol2_filename)
    assert sum([n[0] != 'H' for n in json_names]) == len(mol2_crds)

    mapping = {}  # mol2 atom index to json name
    # Atom names are total mess. We just compare atoms by their coordinates
    for json_crd, json_name in zip(json_crds, json_names):
        flags = [crd_eq(json_crd, i) for i in mol2_crds]
        if sum(flags) == 1:
            matching_index = [idx for idx, f in enumerate(flags) if f][0]
            mapping[matching_index] = json_name
        elif sum(flags) == 0 and json_name[0] == 'H':
            pass  # It's okay, we might have no hydrogens in mol2
        else:
            raise RuntimeError('Unable to find matching atom in MOL2 ({}} for atom {} from JSON ({})'.format(
                input_mol2_filename, json_name, input_json_filename
            ))
    return mapping


def make_restraints(input_json, input_mol2, ref_parts, out_pdb):
    mapping_json = map_mol2_to_json_by_coords(input_json, input_mol2)
    ref_coords = {}

    for ref_part in ref_parts:
        ref_crds = read_mol2(ref_part['ref_lig'])
        for atompair in ref_part['mapping']:
            ref_atom_id = atompair['ref_atom_id']
            target_atom_id = atompair['target_atom_id']
            x, y, z = ref_crds[ref_atom_id]
            name = mapping_json[target_atom_id]
            ref_coords[name] = (x, y, z)

    with open(out_pdb, 'w') as fo:
        num = 0
        for atom_name in ref_coords.keys():
            atom_coords = ref_coords[atom_name]
            num += 1
            line = 'HETATM  {: 3d} {:4s} LIG B   1    {:8.3f}{:8.3f}{:8.3f}  1.00  0.00\n'.format(
                num, atom_name, atom_coords[0], atom_coords[1], atom_coords[2])
            fo.write(line)


def do_single_case(record, name, options):
    receptor_prefix = os.path.join(OUT_DIR, '{}.rec'.format(name))
    ligand_json = os.path.join(OUT_DIR, '{}.lig.json'.format(name))
    restraints_pdb = os.path.join(OUT_DIR, '{}.restraints.pdb'.format(name))
    try:
        run_pdbprep(record['rec'], receptor_prefix, options)
        run_atlas(record['lig'], ligand_json, options)
        make_restraints(ligand_json, record['lig'], record['ref_parts'], restraints_pdb)
    except Exception as exc:
        LOG_DEV.exception('Error processing {}'.format(name), exc_info=exc)
        return None
    else:
        return dict(
            name=name,
            receptor_pdb=receptor_prefix+'.pdb',
            receptor_psf=receptor_prefix+'.psf',
            ligand_json=ligand_json,
            restraints_pdb=restraints_pdb,
            template_properties=record
        )


def main(data, options):
    output = []
    if not os.path.exists(OUT_DIR):
        os.mkdir(OUT_DIR)
    with multiprocessing.Pool(options.get('num_threads', 1)) as pool:
        for record_name, record in data.items():
            output.append(pool.apply_async(
                do_single_case, (record, record_name, options)
            ))
        output = [r.get() for r in output]  # Waiting for all results to finish
    output_good = {r['name']: r for r in output if r is not None}
    return output_good


@click.command()
@click.argument('input', type=click.File('r'))
@click.argument('output', type=click.File('w'))
@click.argument('options', type=click.File('r'))
@click.option('--num_threads', type=int, default=1, help='Number of threads for parallel processing.')
def cli(input, output, options, num_threads):
    input_data = json.load(input)
    input_len = len(input_data.keys())
    if input_len == 0:
        LOG_USR.error('No templates found')  # This string is used to check the job status in runner.py, please change it carefully
        LOG_DEV.info('Got 0 targets to parameterize, nothing to do')
        return

    options_data = json.load(options)
    options_data['num_threads'] = num_threads
    output_data = main(input_data, options_data)
    output_len = len(output_data.keys())

    if input_len != output_len:
        LOG_DEV.warning('Got {} targets to parameterize, {} succeeded'.format(input_len, output_len))
    else:
        LOG_DEV.info('Got {} targets to parameterize, all succeeded'.format(input_len))
    LOG_USR.info(f'We found {input_len} templates, and managed to parameterize {output_len} of them')
    if not output_data:  # If it's empty
        LOG_USR.error('Unfortunately, system preparation failed for all structures')
        raise RuntimeError('Unable to parameterize anything')
    json.dump(output_data, output, indent=4)


if __name__ == '__main__':
    cli()
