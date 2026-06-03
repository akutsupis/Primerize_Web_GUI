import math

from . import util
from .thermo import calc_Tm


def _draw_assembly(sequence, primers, COL_SIZE):
    N_primers = primers.shape[1]
    seq_line_prev = list(' ' * max(len(sequence), COL_SIZE))
    bp_lines = []
    seq_lines = []
    Tms = []

    for i in range(N_primers):
        primer = primers[:, i]
        seq_line = list(' ' * max(len(sequence), COL_SIZE))
        (seg_start, seg_end, seg_dir) = primer

        if (seg_dir == 1):
            seq_line[seg_start: seg_end + 1] = sequence[seg_start: seg_end + 1]

            if (seg_end + 1 < len(sequence)):
                seq_line[seg_end + 1] = '-'
            if (seg_end + 2 < len(sequence)):
                seq_line[seg_end + 2] = '>'

            num_txt = '%d' % (i + 1)
            if (seg_end + 4 + len(num_txt) < len(sequence)):
                offset = seg_end + 4
                seq_line[offset:(offset + len(num_txt))] = num_txt
        else:
            seq_line[seg_start: seg_end + 1] = list(map(util.reverse_complement, sequence[seg_start: seg_end + 1]))

            if (seg_start - 1 >= 0):
                seq_line[seg_start - 1] = '-'
            if (seg_start - 2 >= 0):
                seq_line[seg_start - 2] = '<'

            num_txt = '%d' % (i + 1)
            if (seg_start - 3 - len(num_txt) >= 0):
                offset = seg_start - 3 - len(num_txt)
                seq_line[offset:(offset + len(num_txt))] = num_txt

        bp_line = list(' ' * max(len(sequence), COL_SIZE))
        overlap_seq = ''
        last_bp_pos = 1
        for j in range(len(sequence)):
            if (seq_line_prev[j] in 'ACGT' and seq_line[j] in 'ACGT'):
                bp_line[j] = '|'
                last_bp_pos = j
                overlap_seq += sequence[j]

        if (last_bp_pos > 1):
            Tm = calc_Tm(overlap_seq, 0.2e-6, 0.1, 0.0015)
            Tms.append(Tm)
            Tm_txt = '%2.1f' % Tm
            offset = last_bp_pos + 2
            bp_line[offset:(offset + len(Tm_txt))] = 'x' * len(Tm_txt)

        bp_lines.append(''.join(bp_line))
        seq_lines.append(''.join(seq_line))
        seq_line_prev = seq_line

    print_lines = []
    for i in range(int(math.floor((len(sequence) - 1) / COL_SIZE)) + 1):
        start_pos = COL_SIZE * i
        end_pos = min(COL_SIZE * (i + 1), len(sequence))
        out_line = sequence[start_pos:end_pos]
        print_lines.append(('~', out_line))

        for j in range(len(seq_lines)):
            if (len(bp_lines[j][end_pos:].replace(' ', '')) and ('|' not in bp_lines[j][end_pos:].replace(' ', '')) and (not len(bp_lines[j][:start_pos].replace(' ', '')))):
                bp_line = bp_lines[j][start_pos:].rstrip()
            elif ('|' not in bp_lines[j][start_pos:end_pos]):
                bp_line = ' ' * (end_pos - start_pos + 1)
            else:
                bp_line = bp_lines[j][start_pos:end_pos]
            seq_line = seq_lines[j][start_pos:end_pos]

            if len(bp_line.replace(' ', '')) or len(seq_line.replace(' ', '')):
                print_lines.append(('$', bp_line))
                print_lines.append(('^!'[j % 2], seq_line))
        print_lines.append(('$', ' ' * (end_pos - start_pos + 1)))
        print_lines.append(('=', util.complement(out_line)))
        print_lines.append(('', '\n'))

    return (bp_lines, seq_lines, print_lines, Tms)


def _draw_common(fragments, labels):
    (illustration_1, illustration_2, illustration_3) = ('', '', '')

    if len(fragments[0]) >= len(labels[0]):
        illustration_1 += '\033[91m' + fragments[0][0] + '\033[0m\033[40m' + fragments[0][1:] + '\033[0m'
        illustration_2 += '\033[91m|%s\033[0m' % (' ' * (len(fragments[0]) - 1))
        illustration_3 += '\033[91m%s%s\033[0m' % (labels[0], ' ' * (len(fragments[0]) - len(labels[0])))
    elif fragments[0]:
        illustration_1 += '\033[91m' + fragments[0][0] + '\033[0m\033[40m' + fragments[0][1:] + '\033[0m'
        illustration_2 += '\033[91m%s\033[0m' % (' ' * len(fragments[0]))
        illustration_3 += '\033[91m%s\033[0m' % (' ' * len(fragments[0]))

    if len(fragments[1]) >= len(labels[1]) + len(labels[2]):
        illustration_1 += '\033[44m' + fragments[1][0] + '\033[0m\033[46m' + fragments[1][1:-1] + '\033[0m\033[44m' + fragments[1][-1] + '\033[0m'
        illustration_2 += '\033[92m|%s|\033[0m' % (' ' * (len(fragments[1]) - 2))
        illustration_3 += '\033[92m%s%s%s\033[0m' % (labels[1], ' ' * (len(fragments[1]) - len(labels[1]) - len(labels[2])), labels[2])
    elif fragments[1]:
        if len(fragments[1]) >= len(labels[1]):
            illustration_1 += '\033[44m' + fragments[1][0] + '\033[0m\033[46m' + fragments[1][1:] + '\033[0m'
            illustration_2 += '\033[92m|%s\033[0m' % (' ' * (len(fragments[1]) - 1))
            illustration_3 += '\033[92m%s%s\033[0m' % (labels[1], ' ' * (len(fragments[1]) - len(labels[1])))
        else:
            illustration_1 += '\033[46m' + fragments[1] + '\033[0m'
            illustration_2 += '\033[92m%s\033[0m' % (' ' * len(fragments[1]))
            illustration_3 += '\033[92m%s\033[0m' % (' ' * len(fragments[1]))

    if len(fragments[2]) >= len(labels[3]):
        illustration_1 += '\033[40m' + fragments[2][:-1] + '\033[0m\033[91m' + fragments[2][-1] + '\033[0m'
        illustration_2 += '\033[91m%s|\033[0m' % (' ' * (len(fragments[2]) - 1))
        illustration_3 += '\033[91m%s%s\033[0m' % (' ' * (len(fragments[2]) - len(labels[3])), labels[3])
    elif fragments[2]:
        illustration_1 += '\033[40m' + fragments[2][:-1] + '\033[0m\033[91m' + fragments[2][-1] + '\033[0m'
        illustration_2 += '\033[91m%s\033[0m' % (' ' * len(fragments[2]))
        illustration_3 += '\033[91m%s\033[0m' % (' ' * len(fragments[2]))

    return (illustration_1, illustration_2, illustration_3)


def _draw_region(sequence, params):
    offset = params['offset']
    start = params['which_muts'][0] + offset - 1
    end = params['which_muts'][-1] + offset - 1
    fragments = []

    if start <= 20:
        fragments.append(sequence[:start])
    else:
        fragments.append(sequence[:10] + '......' + sequence[start - 10:start])
    if end - start <= 40:
        fragments.append(sequence[start:end + 1])
    else:
        fragments.append(sequence[start:start + 20] + '......' + sequence[end - 19:end + 1])
    if len(sequence) - end <= 20:
        fragments.append(sequence[end + 1:])
    else:
        fragments.append(sequence[end + 1:end + 11] + '......' + sequence[-10:])

    labels = ['%d' % (1 - offset), '%d' % params['which_muts'][0], '%d' % params['which_muts'][-1], '%d' % (len(sequence) - offset)]
    (illustration_1, illustration_2, illustration_3) = _draw_common(fragments, labels)
    return {'labels': labels, 'fragments': fragments, 'lines': (illustration_1, illustration_2, illustration_3)}


def _draw_str_region(sequence, structures, bps, warnings, params):
    offset = params['offset']
    start = params['which_muts'][0] + offset - 1
    end = params['which_muts'][-1] + offset - 1
    fragments = []

    fragments.append(sequence[:start])
    fragments.append(sequence[start:end + 1])
    fragments.append(sequence[end + 1:])

    labels = ['%d' % (1 - offset), '%d' % params['which_muts'][0], '%d' % params['which_muts'][-1], '%d' % (len(sequence) - offset)]
    (illustration_1, illustration_2, illustration_3) = _draw_common(fragments, labels)

    illustration_str = ''
    bps = [bp for helix in bps for bp in helix]
    for structure in structures:
        this_bps = [bp for helix in util.str_to_bps(structure) for bp in helix]
        this_bps = list(filter(lambda x: (x in bps), this_bps))
        bps = list(filter(lambda x: (x not in this_bps), bps))
        this_nts = [nt for bp in this_bps for nt in bp]
        mismatch = list(filter(lambda x: (x in this_bps), warnings))
        this_mis = [nt for bp in mismatch for nt in bp]

        for i, nt in enumerate(structure):
            if i + 1 in this_nts:
                if i + 1 in this_mis:
                    illustration_str += '\033[41m%s\033[0m' % nt
                else:
                    illustration_str += '\033[43m%s\033[0m' % nt
            else:
                illustration_str += nt
        illustration_str += '\n'

    return {'labels': labels, 'fragments': fragments, 'lines': (illustration_3, illustration_2, illustration_1, illustration_str)}
