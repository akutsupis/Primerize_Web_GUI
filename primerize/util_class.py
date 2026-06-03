import matplotlib
matplotlib.use('SVG')
import matplotlib.pyplot as pyplot
import os
import re

from . import util
from . import util_func
from . import util_server

mpl_version = matplotlib.__version__.split('.')[0]
svg_well_area = 31**2 if mpl_version == '1' else 25**2
svg_well_stroke = 5 if mpl_version == '1' else 4

class Assembly(object):
    """Collection of result data essential for drawing an assembly scheme.

    Args:
        sequence: ``str``: Sequence of assembly design.
        primers: ``list(list(int)``: Numeric representation (end numbering and direction) of primers.
        name: ``str``: Construct prefix/name.
        COL_SIZE: ``int``: `(Optional)` Column width for assembly output. Positive number only.

    Attributes:
        sequence: ``str``: Sequence of assembly design.
        primers: ``list(list(int))``: Numeric representation (end numbering and direction) of primers.
        name: ``str``: Construct prefix/name.
        bp_lines: ``list(str)``: Strings for base-pairing lines (``'|'``).
        seq_lines: ``list(str)``: Strings for primer sequence lines.
        print_lines: ``list(tuple(str, str))``: Strings for all lines assembled, i.e. ``list(tuple('marker', 'print_line'))``.
        Tm_overlaps: ``list(float)``: List of melting temperature for all overlapping regions.

    """

    def __init__(self, sequence, primers, name, COL_SIZE=142):
        self.sequence = sequence
        self.primers = primers
        self.name = name
        (self.bp_lines, self.seq_lines, self.print_lines, self.Tm_overlaps) = util_server._draw_assembly(self.sequence, self.primers, COL_SIZE)

        self.primer_set = []
        for i in range(self.primers.shape[1]):
            primer_seq = self.sequence[self.primers[0, i]:self.primers[1, i] + 1]
            if self.primers[2, i] == -1:
                self.primer_set.append(util.reverse_complement(primer_seq))
            else:
                self.primer_set.append(str(primer_seq))


    def __repr__(self):
        """Representation of the ``Assembly`` class.
        """

        return '\033[94m%s\033[0m {\n    \033[93m\'primers\'\033[0m: %s, \n    \033[93m\'seq_lines\'\033[0m: \033[91mlist\033[0m(\033[91mstring\033[0m * %d), \n    \033[93m\'bp_lines\'\033[0m: \033[91mlist\033[0m(\033[91mstring\033[0m * %d), \n    \033[93m\'print_lines\'\033[0m: \033[91mlist\033[0m(\033[91mtuple\033[0m * %d), \n    \033[93m\'Tm_overlaps\'\033[0m: %s\n}' % (self.__class__, repr(self.primers), len(self.seq_lines), len(self.bp_lines), len(self.print_lines), repr(self.Tm_overlaps))

    def __str__(self):
        """Results of the ``Assembly`` class. Calls ``echo()``.
        """

        return self.echo()


    def echo(self):
        """Print result in rich-text.

        Returns:
            ``str``
        """

        output = ''
        x = 0
        for i, (flag, string) in enumerate(self.print_lines):
            if (flag == '$' and 'xx' in string):
                Tm = '%2.1f' % self.Tm_overlaps[x]
                output += string.replace('x' * len(Tm), '\033[41m%s\033[0m' % Tm) + '\n'
                x += 1
            elif (flag == '^' or flag == '!'):
                num = ''.join(re.findall("[0-9]+", string))
                string = string.replace(num, '\033[100m%s\033[0m' % num) + '\n'
                if flag == '^':
                    string = string.replace('A', '\033[94mA\033[0m').replace('G', '\033[94mG\033[0m').replace('C', '\033[94mC\033[0m').replace('T', '\033[94mT\033[0m')
                else:
                    string = string.replace('A', '\033[95mA\033[0m').replace('G', '\033[95mG\033[0m').replace('C', '\033[95mC\033[0m').replace('T', '\033[95mT\033[0m')
                output += string
            elif (flag == '~'):
                output += '\033[92m%s\033[0m' % string + '\n'
            elif (flag == '='):
                output += '\033[96m%s\033[0m' % string + '\n'
            else:
                output += string + '\n'

        return output[:-1]


    def save(self, path='./', name=None):
        """Save result to text file.

        Args:
            path: ``str``: `(Optional)` Path for file saving. Use either relative or absolute path.
            name: ``str``: `(Optional)` Prefix/name for file name. When nonspecified, current object's name is used.
        """

        name = self.name if name is None else name
        f = open(os.path.join(path, '%s_assembly.txt' % name), 'w')
        lines = self.echo()
        lines += '\n%s%s\tSEQUENCE\n' % ('PRIMERS'.ljust(20), 'LENGTH'.ljust(10))
        for i, primer in enumerate(self.primer_set):
            name = '%s-\033[100m%s\033[0m%s' % (self.name, i + 1, util_func._primer_suffix(i))
            lines += '%s\033[93m%s\033[0m\t%s\n' % (name.ljust(39), str(len(primer)).ljust(10), util_func._primer_suffix(i).replace(' R', primer).replace(' F', primer))

        lines = lines.replace('\033[0m', '').replace('\033[100m', '').replace('\033[92m', '').replace('\033[93m', '').replace('\033[94m', '').replace('\033[95m', '').replace('\033[96m', '').replace('\033[41m', '')
        f.write(lines)
        f.close()



class Plate_96Well(object):
    """Abstraction of 96-well plates.

    Args:
        tag: ``int``: `(Optional)` Mutation library tag. Use **which_lib** number.

    Attributes:
        coords: ``set(str)``: Filled 96-Well Coordinates.
        _data: Data of primers and names, in format of ``dict: { str: tuple(primerize.Mutation, str) }``, i.e. ``dict: {'coord': ('tag', 'primer') }``.
    """

    def __init__(self, tag=1):
        self.coords = set()
        self._data = {}
        self.tag = 'Lib%d-' % tag

    def __repr__(self):
        """Representation of the ``Plate_96Well`` class.
        """

        if len(self):
            return '\033[94m%s\033[0m {\033[93m\'coords\'\033[0m: %s, \033[93m\'data\'\033[0m: \033[91mdict\033[0m(\033[91mtuple\033[0m * %d)}' % (self.__class__, ' '.join(sorted(self.coords)), len(self._data))
        else:
            return '\033[94m%s\033[0m (empty)' % self.__class__

    def __str__(self):
        """Results of the ``Plate_96Well`` class. Calls ``echo()``.
        """

        return self.echo()

    def __len__(self):
        """Number of filled wells.

        Returns:
            ``int``
        """

        return len(self.coords)

    def __contains__(self, coord):
        """Test if data of a given WellPosition is present.

        Args:
            coord: ``str``: WellPosition (e.g. ``'A01'``) for data.

        Returns:
            ``bool``
        """

        return coord in self.coords


    def get(self, coord):
        """Get data of a particular well or number of wells filled.

        Args:
            coord: ``str``: Keyword of parameter. Use WellPosition for well data.

        Returns:
            value of specified **coord**.

        Raises:
            AttributeError: For illegal WellPosition (out of ``range(0, 96) + 1``).
            KeyError: For nonexisted **coord**.
        """

        coord = util_func._format_coord(coord)
        if util.coord_to_num(coord) is None:
            raise AttributeError('\033[41mERROR\033[0m: Illegal coordinate value \033[95m%s\033[0m for \033[94m%s.get()\033[0m.\n' % (coord, self.__class__))
        elif coord in self:
            return self._data[util.coord_to_num(coord)]
        else:
            raise KeyError('\033[41mERROR\033[0m: Non-Existent coordinate value \033[95m%s\033[0m for \033[94m%s.get()\033[0m.\n' % (coord, self.__class__))


    def set(self, coord, tag, primer):
        """Record data of a particular well.

        Args:
            coord: ``str``: WellPosition for data. Use same range as ``get()``. Existing data for the same well is overwritten.
            tag: ``primerize.Mutation``: Mutant representd by ``primerize.Mutation``. ``str`` is only supported for backward compatibility.
            primer: ``str``: Primer seuqence of well. Use sense-strand.

        Raises:
            AttributeError: For illegal WellPosition.
        """

        coord = util_func._format_coord(coord)
        if util.coord_to_num(coord) is None:
            raise AttributeError('\033[41mERROR\033[0m: Illegal coordinate value \033[95m%s\033[0m for \033[94m%s.set()\033[0m.\n' % (coord, self.__class__))
        else:
            self.coords.add(coord)
            self._data[util.coord_to_num(coord)] = (tag, primer)


    def reset(self):
        """Clear current plate data.
        """

        self.coords = set()
        self._data = {}


    def echo(self, ref_primer=''):
        """Print result in rich-text.

        Args:
            ref_primer: ``list(str)``: `(Optional)` List of Wild-type **primer_set** for highlighting. If nonspecified, highlighting is disabled.

        Returns:
            ``str``
        """

        return util_func._print_primer_plate(self, ref_primer)


    def save(self, ref_primer='', file_name='./plate.svg', title=''):
        """Save plate layout to image file (`SVG`).

        Args:
            ref_primer: ``list(str)``: `(Optional)` List of Wild-type primer_set for highlighting. If nonspecified, highlighting is disabled.
            file_name: ``str``: `(Optional)` File name. Include path into **file_name** when specifying. Use either relative or absolute path.
            title: ``str``: `(Optional)` Title to display on image. LaTex NOT supported.
        """

        fig = pyplot.figure()
        ax = pyplot.subplot(111)
        ax.set_aspect('equal')
        pyplot.axis([0, 13.875, 0, 9.375])
        pyplot.xticks([x * 1.125 + 0.75 for x in range(12)], [str(x + 1) for x in range(12)], fontsize=14)
        pyplot.yticks([y * 1.125 + 0.75 for y in range(8)], list('ABCDEFGH'), fontsize=14)
        fig.suptitle(title, fontsize=16, fontweight='bold')

        for edge in ('bottom', 'top', 'left', 'right'):
            ax.spines[edge].set_color('w')
        ax.invert_yaxis()
        ax.xaxis.set_ticks_position('top')
        for tic in ax.xaxis.get_major_ticks():
            tic.tick1On = tic.tick2On = False
        for tic in ax.yaxis.get_major_ticks():
            tic.tick1On = tic.tick2On = False

        (x_green, x_violet, x_gray, y_green, y_violet, y_gray) = ([], [], [], [], [], [])
        for i in range(8):
            for j in range(12):
                num = i + j * 8 + 1
                if util.num_to_coord(num) in self.coords:
                    tag = self._data[num][0]
                    if (isinstance(tag, Mutation) and not tag) or (isinstance(tag, str) and 'WT' in tag):
                        x_green.append(j * 1.125 + 0.75)
                        y_green.append(i * 1.125 + 0.75)
                    else:
                        if ref_primer and ref_primer == self._data[num][1]:
                            x_green.append(j * 1.125 + 0.75)
                            y_green.append(i * 1.125 + 0.75)
                        else:
                            x_violet.append(j * 1.125 + 0.75)
                            y_violet.append(i * 1.125 + 0.75)
                else:
                    x_gray.append(j * 1.125 + 0.75)
                    y_gray.append(i * 1.125 + 0.75)
        ax.scatter(x_gray, y_gray, svg_well_area, c='#ffffff', edgecolor='#333333', linewidth=svg_well_stroke)
        ax.scatter(x_violet, y_violet, svg_well_area, c='#ecddf4', edgecolor='#c28fdd', linewidth=svg_well_stroke)
        ax.scatter(x_green, y_green, svg_well_area, c='#beebde', edgecolor='#29be92', linewidth=svg_well_stroke)

        matplotlib.rcParams['svg.fonttype'] = 'none'
        matplotlib.rcParams['xtick.labelsize'] = 14
        matplotlib.rcParams['ytick.labelsize'] = 14
        pyplot.savefig(file_name, orientation='landscape', format='svg')
        pyplot.close(fig)



class Mutation(object):
    """Collection of mutations for a construct.

    Args:
        mut_list: ``list(str)``: `(Optional)` List of mutations. When nonspecified, an empty instance is created; when specified, it calls ``push()``. An empty instance means no mutations, i.e. Wild-type.

    Attributes:
        _data: Data of mutations, in format of ``dict: { int: tuple(str, str) }``, i.e. ``dict: {'seqpos': ('wt_char', 'mut_char') }``.
    """

    def __init__(self, mut_str=[]):
        self._data = {}
        if mut_str: self.push(mut_str)

    def __repr__(self):
        """Representation of the ``Mutation`` class.
        """

        return '\033[94m%s\033[0m' % self.__class__

    def __str__(self):
        """Results of the ``Mutation`` class. Calls ``echo()``.
        """

        return self.echo()

    def __len__(self):
        """Number of filled wells.

        Returns:
            ``int``
        """

        return len(self._data)

    def __eq__(self, other):
        """Comparison method for whether two ``Mutation`` objects contain the same set of mutations.
        """

        if isinstance(other, str) and other == 'WT': return len(self) == 0
        if isinstance(other, Mutation): other = other.list()
        return other in self and len(self) == len(other)

    def __iter__(self):
        """Iterator through all mutations.
        """

        for k in self._data.keys():
            yield k, self._data[k]

    def __contains__(self, mut_str):
        """Test if a list of given mutation is present.

        Args:
            mut_list: ``list(str)``: Mutations in format of ``'wt_char'``, ``'seq_pos'``, ``'mut_char'``, (e.g. ``['G13C', 'A15T']``).

        Returns:
            ``bool``
        """

        if isinstance(mut_str, str): mut_str = [mut_str]
        if not (mut_str or self._data): return True
        flag = False

        for mut in mut_str:
            seq_org = util.RNA2DNA(mut[0])
            seq_mut = util.RNA2DNA(mut[-1])
            seq_pos = int(mut[1:-1])

            flag = seq_pos in self._data and self._data[seq_pos] == (seq_org, seq_mut)
            if not flag: return flag

        return flag


    def push(self, mut_str):
        """Add a list of mutations.

        Args:
            mut_list: ``list(str)``: Mutations. Valid keywords are the same as ``has()``. Each ``'seq_pos'`` can only be mutated once. Conflicting mutations are overwritten and the most recent one is saved. ``'WT'`` is ignored.

        Raises:
            ValueError: For illegal **mut_str**.
        """

        if isinstance(mut_str, str): mut_str = [mut_str]
        for mut in mut_str:
            if mut == 'WT': continue

            seq_org = util.complement(util.complement(util.RNA2DNA(mut[0])))
            seq_mut = util.complement(util.complement(util.RNA2DNA(mut[-1])))
            seq_pos = int(mut[1:-1])
            if seq_org == seq_mut:
                raise ValueError('\033[41mERROR\033[0m: Unchanged sequence identity by \033[94mMutation\033[0m \033[92m%s\033[0m.\n' % mut)
            self._data[seq_pos] = (seq_org, seq_mut)


    def pop(self, mut_str):
        """Remove a list of mutations.

        Args:
            mut_list: ``list(str)``: Mutations. Valid keywords are the same as ``has()``. Mutations that are not present will result in a premature return with ``False``.

        Returns:
            ``bool``: Whether all mutations in **mut_list** are successfully removed.
        """

        if isinstance(mut_str, str): mut_str = [mut_str]
        for mut in mut_str:
            if mut in self:
                seq_pos = int(mut[1:-1])
                self._data.pop(seq_pos, None)
            else:
                return False
        return True


    def merge(self, other):
        """Merge 2 lists of mutations.

        Args:
            other: ``primerize.Mutation``: Another list of mutations.

        Raises:
            TypeError: For illegal **other**.
        """

        if not isinstance(other, Mutation):
            raise TypeError('\033[41mERROR\033[0m: Illegal input type for \033[94m%s.merge()\033[0m.\n' % self.__class__)
        for mut in other.list():
            self.push(mut)


    def list(self):
        """Return a list of all mutations.

        Returns:
            ``list(str)``
        """

        return list(map(lambda x: '%s%s%s' % (self._data[x][0], x, self._data[x][1]), sorted(self._data.keys())))


    def echo(self):
        """Print result in rich-text, delimited by ``';'``.

        Returns:
            ``str``
        """

        output = []
        for mut in self.list():
            if mut[-2:] == 'WT':
                output.append('\033[100mWT\033[0m')
            else:
                output.append('\033[96m%s\033[0m\033[93m%s\033[0m\033[91m%s\033[0m' % (mut[0], mut[1:-1], mut[-1]))

        output = ';'.join(output) if output else '\033[100mWT\033[0m'
        return output



class Construct_List(object):
    """Collection of mutant constructs. An empty isntance of ``primerize.Mutation`` is always initiated as the first element, i.e. a Wild-type well as the first of list.

    Attributes:
        _data: Data of constructs, in format of ``list(primerize.Mutation)``.
    """

    def __init__(self):
        self._data = [Mutation()]

    def __repr__(self):
        """Representation of the ``Construct_List`` class.
        """

        if len(self):
            return '\033[94m%s\033[0m {\033[91mlist\033[0m(%s * %d)}' % (self.__class__, repr(Mutation()), len(self))
        else:
            return '\033[94m%s\033[0m (empty)' % self.__class__

    def __str__(self):
        """Results of the ``Construct_List`` class. Calls ``echo()``.
        """

        return self.echo()

    def __len__(self):
        """Number of filled wells.

        Returns:
            ``int``
        """

        return len(self._data)

    def __iter__(self):
        """Iterator through all constructs.
        """

        for i in range(len(self._data)):
            yield self._data[i]

    def __contains__(self, mut_list):
        """Test if a list of given mutant construct is present.

        Args:
            mut_list: ``primerize.Mutation``: A mutant represented by ``primerize.Mutation``.

        Returns:
            ``bool``
        """

        if not isinstance(mut_list, Mutation): mut_list = Mutation(mut_list)
        for construct in self._data:
            if construct == mut_list: return True
        return False


    def push(self, mut_list):
        """Add a list of mutations.

        Args:
            mut_list: ``primerize.Mutation``: Mutations. A mutant represented by ``primerize.Mutation``. If the mutant is already present, it will return ``False``.

        Returns:
            ``bool``: Whether **mut_list** is successfully added.
        """

        if not isinstance(mut_list, Mutation): mut_list = Mutation(mut_list)
        if mut_list in self: return False
        self._data.append(mut_list)
        return True


    def pop(self, mut_list):
        """Remove a list of mutations.

        Args:
            mut_list: ``primerize.Mutation``: A mutant represented by ``primerize.Mutation``. Mutant that is not present will result in a premature return with ``False``.

        Returns:
            ``bool``: Whether **mut_list** is successfully removed.
        """

        if not isinstance(mut_list, Mutation): mut_list = Mutation(mut_list)
        for i, construct in enumerate(self._data):
            if construct == mut_list:
                self._data.pop(i)
                return True
        return False


    def merge(self, other):
        """Merge 2 lists of constructs.

        Args:
            other: ``primerize.Construct_List``: Another list of constructs.

        Returns:
            ``primerize.Construct_List``: A list of duplicated constructs between inputs.

        Raises:
            TypeError: For illegal **other**.
        """

        if not isinstance(other, Construct_List):
            raise TypeError('\033[41mERROR\033[0m: Illegal input type for \033[94m%s.merge()\033[0m.\n' % self.__class__)
        repeated = Construct_List()
        for mut in other:
            flag = self.push(mut)
            if not flag: repeated.push(mut)
        repeated.pop('WT')
        return repeated


    def list(self):
        """Return a list of all constructs.

        Returns:
            ``list(list(str))``
        """

        return map(lambda x: x.list(), self._data)


    def echo(self, prefix=''):
        """Print result in rich-text.

        Returns:
            ``str``
        """

        output = ''
        for construct in self._data:
            output += prefix + construct.echo() + '\n'
        return output
