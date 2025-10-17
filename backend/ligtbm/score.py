#!/usr/bin/env python3

import click
import json
import logging.config
import os

CONFIDENCE_THRESHOLD_HIGH = 0.9
CONFIDENCE_THRESHOLD_MEDIUM = 0.65

FILE_DIR = os.path.dirname(os.path.abspath(__file__))
logging.config.fileConfig(os.path.join(FILE_DIR, 'logging.conf'), disable_existing_loggers=False)
LOG_USR = logging.getLogger('user')
LOG_DEV = logging.getLogger('dev')


def _make_model_score(record):
    # get negative average mcs coverage
    coverages_w = [ref_part['refinfo']['mcs']['coverage']
                   for ref_part in record['template_properties']['ref_parts']]
    avg_coverage_w = sum(coverages_w) / len(coverages_w)
    coverages_s = [ref_part['refinfo']['mcs']['strict_coverage']
                   for ref_part in record['template_properties']['ref_parts']]
    avg_coverage_s = sum(coverages_s) / len(coverages_s)
    identities = [ref_part['blast_alignment']['identity']
                  for ref_part in record['template_properties']['ref_parts']]
    avg_identity = sum(identities) / len(identities)
    return -1 * (avg_coverage_w * 100 + avg_coverage_s * 10 + avg_identity)


def _make_model_confidence(record):
    tanimotos = [ref_part['refinfo']['fp_tanimoto']
                 for ref_part in record['template_properties']['ref_parts']]
    avg_tanimoto = sum(tanimotos) / len(tanimotos)
    if avg_tanimoto >= CONFIDENCE_THRESHOLD_HIGH:
        return 'high'
    elif avg_tanimoto >= CONFIDENCE_THRESHOLD_MEDIUM:
        return 'medium'
    else:
        return 'low'


def _make_model_template_id(record):
    ids = [ref_part['blast_alignment']['sbjct_pdbid'] + '_' + ref_part['refinfo']['ref_chemid']
           for ref_part in record['template_properties']['ref_parts']]
    return ';'.join(ids)


def do_single_case(record, name, options):
    if 'template_properties' not in record:
        raise KeyError('Need template_properties for scoring')
    output = dict(
        score=_make_model_score(record),
        confidence=_make_model_confidence(record),
        template_id=_make_model_template_id(record)
    )
    output.update(record)
    return output


def rank_dict(results: dict, key: str):
    results_list = [(k, v) for k, v in results.items()]
    results_sorted = sorted(results_list, key=lambda x: x[1][key])
    results_ranked = {
        k: dict(rank=i, **v)
        for i, (k, v) in enumerate(results_sorted, 1)
    }
    return results_ranked


def main(data, options):
    output = dict()
    for record_name, record in data.items():
        record = data[record_name]
        try:
            retval = do_single_case(record, record_name, options)
        except Exception:
            LOG_DEV.exception('Unable to process record {}'.format(record_name))
        else:
            if retval is not None:
                output[record_name] = retval
                try:
                    output[record_name]['template_properties'] = record['template_properties']
                except KeyError:
                    pass
            else:
                LOG_DEV.warning('Had issues scoring record {}'.format(record_name))
    output_sorted = rank_dict(output, key='score')
    return output_sorted


@click.command()
@click.argument('input', type=click.File('r'))
@click.argument('output', type=click.File('w'))
@click.argument('options', type=click.File('r'))
def cli(input, output, options):
    input_data = json.load(input)
    options_data = json.load(options)
    output_data = main(input_data, options_data)
    if not output_data:  # If it's empty
        LOG_USR.error('Unfortunately, scoring failed for all structures')
        raise RuntimeError('Unable to score any structure')
    json.dump(output_data, output, indent=4)


if __name__ == '__main__':
    cli()
