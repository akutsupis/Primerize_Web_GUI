import argparse
import math
import time
import traceback

from . import util
from . import util_class
from . import util_func
from .primerize_1d import Primerize_1D
from .thermo import Singleton
from .wrapper import Design_Single, Design_Plate


class Primerize_2D(Singleton):
    """Construct a worker for 2D Primer Design (Mutate-and-Map Plates).

    Args:
        offset: ``int``: `(Optional)` Sequence numbering offset, which is one minus the final number of the first nucleotide.
        which_muts: ``list(int)``: `(Optional)` Array of mutation positions. Use numbering based on ``offset``. When nonspecified, the entire sequence is included for mutagenesis.
        which_lib: ``int``: `(Optional)` Mutation library choice. Valid choices are (``1``, ``2``, ``3``)::

            * 1 represents "A->U, U->A, C->G, G->C" library;
            * 2 represents "A->C, U->C, C->A, G->A" library;
            * 3 represents "A->G, U->G, C->U, G->U" library.

        COL_SIZE: ``int``: `(Optional)` Column width for assembly output. Positive number only.
        prefix: ``str``: `(Optional)` Construct prefix/name.

    Returns:
        ``primerize.Primerize_2D``

    Note:
        This ``class`` follows the singleton pattern so that only one instance is created. An instance is already initialized as ``primerize.Primerize_2D``.
    """

    def __init__(self, offset=0, which_muts=[], which_lib=1, COL_SIZE=142, prefix='lib'):
        self.prefix = prefix
        self.offset = offset
        self.which_muts = which_muts
        self.which_lib = max(min(which_lib, 3), 1)
        self.COL_SIZE = max(COL_SIZE, 0)

    def __repr__(self):
        return repr(self.__dict__)

    def __str__(self):
        return repr(self.__dict__)


    def get(self, key):
        """Get current worker parameters.

        Args:
            key: ``str``: Keyword of parameter. Valid keywords are ``'offset'``, ``'which_muts'``, ``'which_lib'``, ``'COL_SIZE'``, ``'prefix'``; case insensitive.

        Returns:
            value of specified **key**.

        Raises:
            AttributeError: For illegal **key**.
        """

        key = key.lower()
        if hasattr(self, key):
            return getattr(self, key)
        elif key == 'col_size':
            return self.COL_SIZE
        else:
            raise AttributeError('\033[41mERROR\033[0m: Unrecognized key \033[92m%s\033[0m for \033[94m%s.get()\033[0m.\n' % (key, self.__class__))


    def set(self, key, value):
        """Set current worker parameters.

        Args:
            key: ``str``: Keyword of parameter. Valid keywords are the same as ``get()``.
            value: ``(auto)``: New value for specified keyword. Type of value must match **key**.

        Raises:
            AttributeError: For illegal **key**.
            ValueError: For illegal **value**.
        """

        key = key.lower()
        if hasattr(self, key):
            if key == 'prefix':
                self.prefix = str(value)
            elif key == 'offset' and isinstance(value, (int, float)):
                self.offset = int(value)
            elif key == 'which_lib' and isinstance(value, (float, int)) and value in (1, 2, 3):
                self.which_lib = int(value)
            elif key == 'which_muts' and (isinstance(value, list) and all(isinstance(x, (float, int)) for x in value)):
                self.which_muts = sorted(set(value))
            elif key == 'col_size' and isinstance(value, int) and value > 0:
                self.COL_SIZE = int(value)
            else:
                raise ValueError('\033[41mERROR\033[0m: Illegal value \033[95m%s\033[0m for key \033[92m%s\033[0m for \033[94m%s.set()\033[0m.\n' % (value, key, self.__class__))
        else:
            raise AttributeError('\033[41mERROR\033[0m: Unrecognized key \033[92m%s\033[0m for \033[94m%s.set()\033[0m.\n' % (key, self.__class__))


    def reset(self):
        """Reset current worker parameters to default.
        """

        self.prefix = 'lib'
        self.offset = 0
        self.which_muts = []
        self.which_lib = 1
        self.COL_SIZE = 142


    def design(self, sequence, primer_set=[], offset=None, which_muts=None, which_lib=None, prefix=None, is_force=False):
        """Run design code to get library plates for input sequence according to specified library options. Current worker parameters are used for nonspecified optional arguments.

        Args:
            job_1d: ``primerize.Design_Single``: Result of ``primerize.Primerize_1D.design()``. Its ``sequence``, ``primer_set``, and ``prefix`` are used.
            offset: ``int``: `(Optional)` Sequence numbering offset.
            which_muts: ``list(int)``: `(Optional)` Array of mutation positions.
            which_lib: ``int``: `(Optional)` Mutation library choice.

        Returns:
            ``primerize.Design_Plate``
        """

        if isinstance(sequence, Design_Single):
            design_1d = sequence
            sequence = design_1d.sequence
            primer_set = design_1d.primer_set
            prefix = design_1d.name

        offset = self.offset if offset is None else offset
        which_muts = self.which_muts if which_muts is None else which_muts
        which_lib = self.which_lib if which_lib is None else which_lib
        prefix = self.prefix if prefix is None else prefix

        if len(primer_set) % 2:
            raise ValueError('\033[41mERROR\033[0m: Illegal length \033[95m%s\033[0m of value for params \033[92mprimer_set\033[0m for \033[94m%s.design()\033[0m.\n' % (len(primer_set), self.__class__))

        name = prefix
        sequence = util.RNA2DNA(sequence)
        N_BP = len(sequence)
        params = {'offset': offset, 'which_muts': which_muts, 'which_lib': which_lib, 'N_BP': N_BP, 'type': 'Mutate-and-Map'}
        data = {'plates': [], 'assembly': [], 'constructs': []}

        is_success = True
        primer_set = list(map(util.RNA2DNA, primer_set))
        if not primer_set:
            if is_force:
                prm = Primerize_1D()
                res = prm.design(sequence)
                if res.is_success:
                    primer_set = res.primer_set
                else:
                    is_success = False
                    print('\033[41mFAIL\033[0m: \033[91mNO Solution\033[0m found under given contraints.\n')
            else:
                print('\033[93mWARNING\033[0m: Please run \033[34mPrimerize_1D.design()\033[0m first to get a solution for \033[92mprimer_set\033[0m.\n')
                is_success = False

        if not is_success:
            return Design_Plate({'sequence': sequence, 'name': name, 'is_success': is_success, 'primer_set': primer_set, 'params': params, 'data': data})

        if not which_muts:
            which_muts = list(range(1 - offset, N_BP + 1 - offset))
        else:
            which_muts = list(filter(lambda x: (x >= 1 - offset and x < N_BP + 1 - offset), which_muts))
        which_lib = which_lib[0] if isinstance(which_lib, list) else which_lib

        N_primers = len(primer_set)
        N_constructs = 1 + len(which_muts)
        N_plates = int(math.floor((N_constructs - 1) / 96.0) + 1)
        params.update({'which_muts': which_muts, 'which_lib': which_lib, 'N_PRIMER': N_primers, 'N_PLATE': N_plates, 'N_CONSTRUCT': N_constructs})

        (primers, is_success) = util_func._get_primer_index(primer_set, sequence)
        if not is_success:
            print('\033[41mFAIL\033[0m: \033[91mMismatch\033[0m of given \033[92mprimer_set\033[0m for given \033[92msequence\033[0m.\n')
            return Design_Plate({'sequence': sequence, 'name': name, 'is_success': is_success, 'primer_set': primer_set, 'params': params, 'data': data})

        assembly = util_class.Assembly(sequence, primers, name, self.COL_SIZE)
        constructs = util_class.Construct_List()
        plates = [[util_class.Plate_96Well(which_lib) for i in range(N_plates)] for j in range(N_primers)]
        print('Filling out sequences ...')

        try:
            for m_pos in range(-1, len(which_muts)):
                # m is actual position along sequence
                m = -1 if m_pos == -1 else offset + which_muts[m_pos] - 1
                if m != -1:
                    constructs.push('%s%d%s' % (sequence[m], which_muts[m_pos], util.get_mutation(sequence[m], which_lib)))

            plates = util_func._mutate_primers(plates, primers, primer_set, offset, constructs, which_lib)
            print('\033[92mSUCCESS\033[0m: Primerize 2D design() finished.\n')
        except:
            is_success = False
            print(traceback.format_exc())
            print('\033[41mERROR\033[0m: Primerize 2D design() encountered error.\n')

        data.update({'plates': plates, 'assembly': assembly, 'constructs': constructs})
        return Design_Plate({'sequence': sequence, 'name': name, 'is_success': is_success, 'primer_set': primer_set, 'params': params, 'data': data})


def main():
    parser = argparse.ArgumentParser(description='\033[92mPrimerize 2D Mutate-and-Map Plate Design\033[0m', epilog='\033[94mby Siqi Tian, 2016\033[0m', add_help=False)
    parser.add_argument('sequence', type=str, help='DNA Template Sequence')
    parser.add_argument('-p', metavar='prefix', type=str, help='Display Name of Construct', dest='prefix', default='lib')
    group1 = parser.add_argument_group('advanced options')
    group1.add_argument('-s', metavar='PRIMER_SET', type=str, nargs='+', help='Set of Primers for Assembly (Default runs Primerize 1D)', dest='primer_set', action='append')
    group1.add_argument('-o', metavar='OFFSET', type=int, help='Sequence Numbering Offset', dest='offset', default=0)
    group1.add_argument('-l', metavar='MUT_START', type=int, help='First Position of Mutagenesis (Inclusive), numbering with OFFSET applied', dest='mut_start', default=None)
    group1.add_argument('-u', metavar='MUT_END', type=int, help='Last Position of Mutagenesis (Inclusive), numbering with OFFSET applied', dest='mut_end', default=None)
    group1.add_argument('-w', metavar='LIB', type=int, choices=(1, 2, 3), help='Mutation Library Choices {1, 2, 3}', dest='which_lib', default=1)
    group2 = parser.add_argument_group('commandline options')
    group2.add_argument('-q', '--quiet', action='store_true', dest='is_quiet', help='Suppress Results Printing to stdout')
    group2.add_argument('-e', '--excel', action='store_true', dest='is_excel', help='Write Order Table to Excel File(s)')
    group2.add_argument('-i', '--image', action='store_true', dest='is_image', help='Save Layout to Image File(s)')
    group2.add_argument('-t', '--text', action='store_true', dest='is_text', help='Save Construct and Assembly to Text File(s)')
    group2.add_argument('-h', '--help', action='help', help='Show this Help Message')
    args = parser.parse_args()

    t0 = time.time()
    args.primer_set = [] if args.primer_set is None else args.primer_set[0]
    (which_muts, _, _) = util.get_mut_range(args.mut_start, args.mut_end, args.offset, args.sequence)

    prm = Primerize_2D()
    res = prm.design(args.sequence, args.primer_set, args.offset, which_muts, args.which_lib, args.prefix, True)
    if res.is_success:
        if not args.is_quiet:
            print(res)
        if args.is_excel:
            res.save('table')
        if args.is_image:
            res.save('image')
        if args.is_text:
            res.save('constructs')
            res.save('assembly')

    print('Time elapsed: %.1f s.' % (time.time() - t0))


if __name__ == "__main__":
    main()

