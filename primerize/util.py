import re


def DNA2RNA(sequence):
    """Convert a DNA sequence input to RNA.

    Args:
        sequence: ``str``: Input DNA sequence.

    Returns:
        ``str``: String of RNA
    """

    return sequence.upper().replace('T', 'U')


def RNA2DNA(sequence):
    """Convert a RNA sequence input to DNA.

    Args:
        sequence: ``str``: Input RNA sequence.

    Returns:
        ``str``: String of DNA
    """

    return sequence.upper().replace('U', 'T')


def complement(sequence):
    """Convert a DNA sequence input to its complement strand.

    Args:
        sequence: ``str``: Input DNA sequence.

    Returns:
        ``str``: String of complement DNA strand.

    Raises:
        ValueError: For illegal **sequence**.
    """

    rc_dict = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'U': 'A'}
    try:
        sequence = list(map(lambda x: rc_dict[x], list(sequence)))
    except KeyError:
        raise ValueError('\033[41mERROR\033[0m: Illegal sequence value \033[95m%s\033[0m for \033[94mcomplement()\033[0m.\n' % sequence)

    return ''.join(sequence)


def reverse(sequence):
    """Convert a DNA sequence input to its reverse order.

    Args:
        sequence: ``str``: Input DNA sequence.

    returns:
        ``str``: String of reverse DNA strand.
    """

    return sequence[::-1]


def reverse_complement(sequence):
    """Convert a DNA sequence input to its reverse complement strand.

    Args:
        sequence: ``str``: Input DNA sequence.

    Returns:
        ``str``: String of reverse complement DNA strand.
    """

    return complement(reverse(sequence))


def coord_to_num(coord):
    """Convert a 96-Well Coordinate string to number (1-based).

    Args:
        coord: ``str``: Input WellPosition coordinate string, e.g. ``'A01'``.

    Returns:
        ``int`` or ``None`` if illegal input.
    """

    if not isinstance(coord, str): return None
    coord = re.findall('^([A-H]){1}(0[1-9]|1[0-2]){1}$', coord.upper().strip())
    if not coord: return None
    coord = ''.join(coord[0])
    row = 'ABCDEFGH'.find(coord[0])
    col = int(coord[1:])
    return (col - 1) * 8 + row + 1


def num_to_coord(num):
    """Convert a 96-Well Coordinate number (1-based) to string.

    Args:
        num: ``int``: Input WellPosition coordinate number, e.g. ``96``.

    Returns:
        ``str`` or ``None`` if illegal input.
    """

    if not isinstance(num, int) or num < 1 or num > 96: return None
    row = 'ABCDEFGH'[(num - 1) % 8]
    col = (num - 1) / 8 + 1
    return '%s%0*d' % (row, 2, col)


def get_mut_range(mut_start, mut_end, offset, sequence):
    """Validate and calculate mutation range based on input sequence and offset. If mutation range exceeds possible range, the maximum possible range is returned.

    Args:
        mut_start: ``int``: Lower limit of mutation range, should be based on **offset**.
        mut_end: ``int``: Upper limit of mutation range, should be based on **offset**.
        offset: ``int``: Index numbering offset.
        sequence: ``str``: The sequence (length used).

    Returns:
        ``(which_muts, mut_start, mut_end)``

        - **which_muts** - ``list(int)``: The final range of mutations.
        - **mut_start** - ``int``: The valid **mut_start**.
        - **mut_end** - ``int``: The valid **mut_end**.
    """

    if (not mut_start) or (mut_start is None): mut_start = 1 - offset
    mut_start = min(max(mut_start, 1 - offset), len(sequence) - offset)
    if (not mut_end) or (mut_end is None): mut_end = len(sequence) - offset
    mut_end = max(min(mut_end, len(sequence) - offset), 1 - offset)
    which_muts = list(range(mut_start, mut_end + 1))
    return (which_muts, mut_start, mut_end)


def get_mutation(nt, lib):
    """Mutate a single nucleotide.

    Args:
        nt: ``str``: The nucleotide of interest.
        lib: ``int``: The mutation library choice; choose from (``1``, ``2``, ``3``, ``4``)::

            * 1 represents "A->U, U->A, C->G, G->C",
            * 2 represents "A->C, U->C, C->A, G->A",
            * 3 represents "A->G, U->G, C->U, G->U",
            * 4 represents "A->C, U->G, C->A, G->U",
            * 5 represents "A->C, U->G, C->G, G->C".

    Returns:
        ``str``

    Raises:
        ValueError: For illegal **lib** input.
    """

    libs = {1: 'TAGC', 2: 'CCAA', 3: 'GGTT', 4: 'CGAT', 5: 'CGGC'}
    if lib not in libs:
        raise ValueError('\033[41mERROR\033[0m: Illegal value \033[95m%s\033[0m for params \033[92mwhich_lib\033[0m.\n' % lib)
    else:
        idx = 'ATCG'.find(nt)
        return libs[lib][idx]


def str_to_bps(structure, offset=0):
    """Convert a dot-bracket secondary structure into base-pair tuples.

    Args:
        structure: ``str``: Input secondary struture.
        offset: ``int``: `(Optional)` Index numbering offset for output numbers.

    Returns:
        ``list(list(tuple(int, int)))``: List of helices, and each helix is a list of tuple of base-pairs with their ``seqpos``.
    """

    (lbs, lbs_pk, bps, bps_pk, helices) = ([], [], [], [], [])

    for i, char in enumerate(structure):
        if char == '(':
            lbs.append(i + 1 - offset)
        elif char == ')':
            bps.append((lbs[-1], i + 1 - offset))
            lbs.pop(-1)
        elif char == '[':
            lbs_pk.append(i + 1 - offset)
        elif char == ']':
            bps_pk.append((lbs_pk[-1], i + 1 - offset))
            lbs_pk.pop(-1)
        else:
            if len(bps):
                helices.append(sorted(bps, key=lambda x: x[0]))
                bps = []
            elif len(bps_pk):
                helices.append(sorted(bps_pk, key=lambda x: x[0]))
                bps_pk = []

    if lbs or lbs_pk:
        raise ValueError('\033[41mERROR\033[0m: Unbalanced \033[92mstructure\033[0m "\033[95m%s\033[0m".\n' % structure)
    return helices


def diff_bps(structures, offset=0, flag=True):
    """Find base-pairs that are not present in all secondary structure inputs. Each input secondary structure is compared to all the others.

    Args:
        structures: ``list(str)``: Input secondary structures.
        offst: ``int``: `(Optional)` Index numbering offset for output numbers.
        flag: ``bool``: `(Optional)` Overriding flag for excluding shared helices.

    Returns:
        ``list(list(tuple(int, int)))``: List of helices, and each helix is a list of tuple of base-pairs with their ``seqpos``.
    """

    if isinstance(structures, str): structures = [structures]

    if len(structures) == 1:
        return str_to_bps(structures[0], offset)
    else:
        # collapse helices for all base-pairs in one layer
        helix_all = [helix for structure in structures for helix in str_to_bps(structure, offset)]
        # convert to stings of "int@int" for easy filtering of repeats
        bps_all = ['%d@%d' % (bp[0], bp[1]) for helix in helix_all for bp in helix]
        if flag:
            # remove pairs that present in all structures
            bps = list(filter(lambda x: (bps_all.count(x) < len(structures)), set(bps_all)))
        else:
            # remove repeats
            bps = set(bps_all)
        # convert back to tuple(int, int)
        bps = [(int(x[0]), int(x[1])) for x in map(lambda x: x.split('@'), bps)]

        for i in range(len(helix_all)):
            helix = helix_all[i]
            # remove base-pairs that not made through
            helix_all[i] = list(filter(lambda x: x in bps, helix))
            # remove base-pairs in bps that are taken in helix_all
            bps = list(filter(lambda x: x not in helix, bps))
        # remove empty structure []
        return list(filter(len, helix_all))


def valid_WC_pair(nt_1, nt_2):
    """Check if two nucleotides form a valid Watson-Crick base-pair.

    Args:
        nt_1: ``str``: Nucleotide.
        nt_2: ``str``: Nucleotide.

    Returns:
        ``bool``
    """

    pair = ''.join(sorted([DNA2RNA(nt_1), DNA2RNA(nt_2)]))
    return pair in ['AU', 'CG', 'GU']


