#!/usr/bin/env python3

import click
import json
import logging.config
import os
import subprocess


RMIN_TIMEOUT = 120  # seconds
RMIN_RESTRAINT_SCALE = 10.0  # Value used in D3R GC4 Stage 1a. Now rmin does not scale internally, so it's x10 of what was passed as a parameter before
RMIN_DENSITY_SCALE = 10000.0  # TODO: this needs some testing
RMIN_CUTOFF = 10.0  # Angstroms. Value used in D3R GC4 Stage 1a
RMIN_PROTOCOL = 1  # 1 - Simple, 2 - Multistage
OUT_DIR = 'minimized'
MAX_ENERGY_PER_BOND = 5  # kcal/mol. If higher, reject the structure

FILE_DIR = os.path.dirname(os.path.abspath(__file__))
logging.config.fileConfig(os.path.join(FILE_DIR, 'logging.conf'), disable_existing_loggers=False)
LOG_USR = logging.getLogger('user')
LOG_DEV = logging.getLogger('dev')


def run_rmin(receptor_pdb, receptor_psf, ligand_json, restraints_pdb, output_pdb, output_dat, protocol, options):
    # The `timeout` parameter of subprocess.Popen (and its siblings) only kills the called process, not its children.
    # To kill children (they're the ones that hang) we have to do the process group magick
    args = [
        options['RMIN_BIN'],
        receptor_pdb, receptor_psf, ligand_json, restraints_pdb,
        options['CHARMM_PRM'], options['CHARMM_RTF'], options['LIBMOL_PRM'],
        output_pdb, output_dat,
        '--spring', str(RMIN_RESTRAINT_SCALE),
        '--density', str(RMIN_DENSITY_SCALE),
        '--cutoff', str(RMIN_CUTOFF),
        '--protocol', str(protocol)]
    LOG_DEV.info('Starting {}'.format(' '.join(args)))
    subprocess.check_call(args, timeout=RMIN_TIMEOUT)


def do_single_case(record, name, options):
    out_pdb = os.path.join(OUT_DIR, '{}.minimized.pdb'.format(name))
    out_dat = os.path.join(OUT_DIR, '{}.minimized.energy.json'.format(name))
    with open(record['ligand_json']) as fp:
        ligand_data = json.load(fp)

    try:
        run_rmin(
            record['receptor_pdb'], record['receptor_psf'], record['ligand_json'], record['restraints_pdb'],
            out_pdb, out_dat, RMIN_PROTOCOL, options)
    except Exception as exc:
        LOG_DEV.exception('Error processing {}'.format(name), exc_info=exc)
        return None
    if not os.path.exists(out_pdb):
        LOG_DEV.error('RMin did not produce output PDB file for {}'.format(name))
        return None

    with open(out_dat) as fp:
        energies = json.load(fp)

    output = dict(minimized_pdb=out_pdb)
    output.update(energies)
    total_energy = energies['total']

    LOG_DEV.info('Total energy for {} is {}'.format(out_pdb, total_energy))
    if total_energy > 0:
        LOG_DEV.warning('Positive total energy: {}'.format(json.dumps(output)))
        return None

    max_bond_energy = len(ligand_data['bonds']) * MAX_ENERGY_PER_BOND
    if energies['bonds'] > max_bond_energy:
        LOG_DEV.warning('Covalent bond energy > {}: {}'.format(max_bond_energy, json.dumps(output)))
        return None

    return output


def main(data, options):
    output = dict()
    if not os.path.exists(OUT_DIR):
        os.mkdir(OUT_DIR)
    for record_name, record in data.items():
        record = data[record_name]
        try:
            retval = do_single_case(record, record_name, options)
        except Exception:
            LOG_DEV.exception('Unable to process record {}'.format(record_name))
        else:
            if retval is not None:
                output[record_name] = retval
                output[record_name]['template_properties'] = record['template_properties']
            else:
                LOG_DEV.warning('Had issues processing record {}'.format(record_name))
    output_good = {k: v for k, v in output.items() if v is not None}
    return output_good


@click.command()
@click.argument('input', type=click.File('r'))
@click.argument('output', type=click.File('w'))
@click.argument('options', type=click.File('r'))
def cli(input, output, options):
    input_data = json.load(input)
    options_data = json.load(options)
    output_data = main(input_data, options_data)
    if not output_data:  # If it's empty
        LOG_USR.error('Unfortunately, minimization failed for all structures')
        raise RuntimeError('Unable to minimize any structure')
    json.dump(output_data, output, indent=4)


if __name__ == '__main__':
    cli()
