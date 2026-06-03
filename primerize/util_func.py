import math
import numpy
import os
import re
import xlwt

from . import util


def _primer_suffix(num):
    if num % 2:
        return '\033[95m R\033[0m'
    else:
        return '\033[94m F\033[0m'


def _format_coord(coord):
    coord = re.findall('^([A-H]{1}[0-9]{1,2})$', coord.upper().strip())
    if not coord: return None
    coord = coord[0]
    return coord[0] + coord[1:].zfill(2)


def _get_primer_index(primer_set, sequence):
    N_primers = len(primer_set)
    coverage = numpy.zeros((1, len(sequence)))
    primers = numpy.zeros((3, N_primers))

    for n in range(N_primers):
        primer = util.RNA2DNA(primer_set[n])
        if n % 2:
            i = sequence.find(util.reverse_complement(primer))
        else:
            i = sequence.find(primer)
        if i == -1:
            return ([], False)
        else:
            start_pos = i
            end_pos = i + len(primer_set[n]) - 1
            seq_dir = math.copysign(1, 0.5 - n % 2)
            primers[:, n] = [start_pos, end_pos, seq_dir]
            coverage[0, start_pos:(end_pos + 1)] = 1

    return (primers.astype(int), coverage.all())


def _print_primer_plate(plate, ref_primer):
    if not plate: return '(empty)\n'
    string = ''
    for key in sorted(plate._data):
        string += '\033[94m%s\033[0m' % util.num_to_coord(key).ljust(5)
        mut = plate._data[key][0]
        if not isinstance(mut, str):
            offset = 20 if mut else 30
            lbl = mut.echo()
            string += plate.tag + lbl.ljust(max(offset + 27 * len(mut), len(lbl)))
        else:
            if mut[-2:] == 'WT':
                string += ('%s\033[100m%s\033[0m' % (mut[:-2], mut[-2:])).ljust(28)
            else:
                string += ('%s\033[96m%s\033[0m\033[93m%s\033[0m\033[91m%s\033[0m' % (mut[:5], mut[5], mut[6:-1], mut[-1])).ljust(50)

        if ref_primer:
            for i, ref in enumerate(ref_primer):
                if ref != plate._data[key][1][i]:
                    string += '\033[41m%s\033[0m' % plate._data[key][1][i]
                else:
                    string += plate._data[key][1][i]
        else:
            string += plate._data[key][1]
        string += '\n'

    return string


def _print_pair_mismatch_warning(sequence, warnings, offset):
    string = ''
    for pair in warnings:
        string += '\033[93mWARNING\033[0m: Mismatch in base-pair between \033[96m%s\033[0m\033[93m%s\033[0m and \033[96m%s\033[0m\033[93m%s\033[0m.\n' % (sequence[pair[0] - 1], pair[0] - offset, sequence[pair[1] - 1], pair[1] - offset)
    return string


def _save_plate_layout(plates, ref_primer=[], prefix='', path='./'):
    for k in range(len(plates[0])):
        for p in range(len(plates)):
            primer_sequences = plates[p][k]
            num_primers_on_plate = len(primer_sequences)

            if num_primers_on_plate:
                if num_primers_on_plate == 1 and 'A01' in primer_sequences:
                    tag = primer_sequences.get('A01')[0]
                    if (not tag) or (isinstance(tag, str) and 'WT' in tag): continue

                file_name = os.path.join(path, '%s_plate_%d_primer_%d.svg' % (primer_sequences.tag[:-1], k + 1, p + 1))
                print('Creating plate image: \033[94m%s\033[0m.' % file_name)
                title = '%s_plate_%d_primer_%d' % (prefix, k + 1, p + 1)
                primer_sequences.save(ref_primer[p], file_name, title)


def _save_construct_key(keys, name, path='./', prefix=''):
    prefix = 'Lib%s-' % prefix if prefix else ''
    print('Creating keys file ...')
    lines = keys.echo(prefix)
    lines = lines.replace('\033[100m', '').replace('\033[96m', '').replace('\033[93m', '').replace('\033[91m', '').replace('\033[0m', '')
    open(os.path.join(path, '%s_keys.txt' % name), 'w').write(lines)


def _save_structures(structures, warnings, sequence, offset, name, path='./'):
    print('Creating structures file ...')
    lines = structures
    if warnings:
        lines.extend(['', 'WARNINGS:', ''])
        for pair in warnings:
            lines.append('Mismatch in base-pair between %s%d and %s%d.' % (sequence[pair[0] - 1], pair[0] - offset, sequence[pair[1] - 1], pair[1] - offset))
    open(os.path.join(path, '%s_structures.txt' % name), 'w').write('\n'.join(lines))


def _save_plates_excel(plates, ref_primer=[], prefix='', path='./'):
    for k in range(len(plates[0])):
        file_name = os.path.join(path, '%s_plate_%d.xls' % (prefix, k + 1))
        print('Creating plate file: \033[94m%s\033[0m.' % file_name)
        workbook = xlwt.Workbook()

        for p in range(len(plates)):
            primer_sequences = plates[p][k]
            num_primers_on_plate = len(primer_sequences)

            if num_primers_on_plate:
                if num_primers_on_plate == 1 and 'A01' in primer_sequences:
                    tag = primer_sequences.get('A01')[0]
                    if (not tag) or (isinstance(tag, str) and 'WT' in tag): continue

                sheet = workbook.add_sheet('primer_%d' % (p + 1))
                sheet.col(1).width = 256 * 15
                sheet.col(2).width = 256 * 75

                sheet.write(0, 0, 'WellPosition', xlwt.easyxf('font: bold 1'))
                sheet.write(0, 1, 'Name', xlwt.easyxf('font: bold 1'))
                sheet.write(0, 2, 'Sequence', xlwt.easyxf('font: bold 1'))
                sheet.write(0, 3, 'Notes', xlwt.easyxf('font: bold 1'))

                for i, row in enumerate(sorted(primer_sequences._data)):
                    tag = primer_sequences._data[row][0]
                    primer = primer_sequences._data[row][1]
                    if not isinstance(tag, str):
                        format = 'font: color blue,' if (not tag or primer == ref_primer[p]) else 'font: color black,'
                        tag = ';'.join(tag.list()) if tag else 'WT'
                        tag = primer_sequences.tag + tag
                    else:
                        format = 'font: color blue,' if 'WT' in tag else 'font: color black,'

                    sheet.write(i + 1, 0, util.num_to_coord(row), xlwt.easyxf(format + '  italic 1'))
                    sheet.write(i + 1, 1, tag, xlwt.easyxf(format))
                    sheet.write(i + 1, 2, primer, xlwt.easyxf(format))

        if len(workbook._Workbook__worksheets): workbook.save(file_name)


def _mutate_primers(plates, primers, primer_set, offset, constructs, which_lib=1, is_fillWT=False):
    for i, mut in enumerate(constructs):
        plate_num = int(math.floor(i / 96.0))
        plate_pos = i % 96 + 1
        well_tag = util.num_to_coord(plate_pos)
        # keep track of unmatched mutations
        is_valid = mut.list()

        for p in range(len(primer_set)):
            wt_primer = primer_set[p]
            if mut == 'WT':
                well_name = 'Lib%d-%s' % (which_lib, 'WT')
                plates[p][plate_num].set(well_tag, mut, wt_primer)
                continue

            mut_primer = util.reverse_complement(wt_primer) if primers[2, p] == -1 else wt_primer
            for k, seq in mut:
                k = k + offset - 1
                if (k >= primers[0, p] and k <= primers[1, p]):
                    m_shift = int(k - primers[0, p])
                    mut_primer = list(mut_primer)
                    # check for mismatch in mutations
                    if seq[0] != mut_primer[m_shift]:
                        raise ValueError('\033[41mERROR\033[0m: Mismatch of \033[94mMutation\033[0m %s with input sequence "\033[95m%s\033[0m" at positon \033[92m%d\033[0m.\n' % (mut, mut_primer[m_shift], k))
                    mut_primer[m_shift] = seq[1]
                    mut_primer = ''.join(mut_primer)

                    seq = '%s%d%s' % (seq[0], k - offset + 1, seq[1])
                    if seq in is_valid: is_valid.remove(seq)

            mut_primer = util.reverse_complement(mut_primer) if primers[2, p] == -1 else mut_primer
            if mut_primer != wt_primer or is_fillWT:
                well_name = 'Lib%d-%s' % (which_lib, mut) if mut_primer != wt_primer else 'Lib%d-%s' % (which_lib, 'WT')
                plates[p][plate_num].set(well_tag, mut, mut_primer)

        for m in is_valid:
            print('\033[93mWARNING\033[0m: Unmatched (or out of bound) \033[94mMutation\033[0m position \033[92m%s\033[0m in construct %s.\n' % (m, mut))

    return plates
