import argparse
import math
import time
import traceback

from . import util
from . import util_class
from . import util_func
from .primerize_1d import Primerize_1D
from .wrapper import Design_Single, Design_Plate
from .thermo import Singleton


class Primerize_Custom(Singleton):
    """Construct a worker for Custom Primer Design (Custom Plates).

    Args:
        offset: ``int``: `(Optional)` Sequence numbering offset, which is one minus the final number of the first nucleotide.
        COL_SIZE: ``int``: `(Optional)` Column width for assembly output. Positive number only.
        prefix: ``str``: `(Optional)` Construct prefix/name.

    Returns:
        ``primerize.Primerize_Custom``

    Note:
        This ``class`` follows the singleton pattern so that only one instance is created. An instance is already initialized as ``primerize.Primerize_Custom``.
    """

    def __init__(self, offset=0, COL_SIZE=142, prefix='lib'):
        self.prefix = prefix
        self.offset = offset
        self.COL_SIZE = max(COL_SIZE, 0)

    def __repr__(self):
        return repr(self.__dict__)

    def __str__(self):
        return repr(self.__dict__)


    def get(self, key):
        """Get current worker parameters.

        Args:
            key: ``str``: Keyword of parameter. Valid keywords are ``'offset'``, ``'COL_SIZE'``, ``'prefix'``; case insensitive.

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
            elif key == 'col_size' and isinstance(value, int) and value > 0:
                self.COL_SIZE = int(value)
            else:
                raise ValueError('\033[41mERROR\033[0m: Illegal value \033[95m%s\033[0m for key \033[92m%s\033[0m for \033[94m%s.set()\033[0m.\n' % (value, key, self.__class__))
        else:
            raise AttributeError('\033[41mERROR\033[0m: Unrecognized key \033[92m%s\033[0m for \033[94m%s.get()\033[0m.\n' % (key, self.__class__))


    def reset(self):
        """Reset current worker parameters to default.
        """

        self.prefix = 'lib'
        self.offset = 0
        self.COL_SIZE = 142


    def design(self, sequence, primer_set=[], mut_list=[], offset=None, prefix=None, is_force=False):
        """Run design code to get customized plates for input sequence and specified list of constructs. Current worker parameters are used for nonspecified optional arguments.

        Args:
            job_1d: ``primerize.Design_Single``: Result of ``primerize.Primerize_1D.design()``. Its ``sequence``, ``primer_set``, and ``prefix`` are used.
            mut_list: ``primerize.Construct_List``: List of constructs for design.
            offset: ``int``: `(Optional)` Sequence numbering offset.

        Returns:
            ``primerize.Design_Plate``
        """

        if isinstance(sequence, Design_Single):
            design_1d = sequence
            sequence = design_1d.sequence
            primer_set = design_1d.primer_set
            prefix = design_1d.name

        offset = self.offset if offset is None else offset
        prefix = self.prefix if prefix is None else prefix

        if len(primer_set) % 2:
            raise ValueError('\033[41mERROR\033[0m: Illegal length \033[95m%s\033[0m of value for params \033[92mprimer_set\033[0m for \033[94m%s.design()\033[0m.\n' % (len(primer_set), self.__class__))

        name = prefix
        sequence = util.RNA2DNA(sequence)
        N_BP = len(sequence)
        params = {'offset': offset, 'N_BP': N_BP, 'type': 'Custom'}
        data = {'plates': [], 'assembly': [], 'constructs': []}

        is_success = True
        if not isinstance(mut_list, util_class.Construct_List):
            is_success = False
            print('\033[41mFAIL\033[0m: \033[91mIllegal\033[0m type of given \033[92mmut_list\033[0m. Should be \033[94mConstruct_List\033[0m.\n')
            return Design_Plate({'sequence': sequence, 'name': name, 'is_success': is_success, 'primer_set': primer_set, 'params': params, 'data': data})

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

        N_primers = len(primer_set)
        params.update({'N_PRIMER': N_primers})
        (primers, is_success) = util_func._get_primer_index(primer_set, sequence)
        if not is_success:
            print('\033[41mFAIL\033[0m: \033[91mMismatch\033[0m of given \033[92mprimer_set\033[0m for given \033[92msequence\033[0m.\n')
            return Design_Plate({'sequence': sequence, 'name': name, 'is_success': is_success, 'primer_set': primer_set, 'params': params, 'data': data})

        assembly = util_class.Assembly(sequence, primers, name, self.COL_SIZE)
        constructs = mut_list
        data.update({'assembly': assembly, 'constructs': constructs})

        which_lib = 0
        N_constructs = len(constructs)
        N_plates = int(math.floor((N_constructs - 1) / 96.0) + 1)
        plates = [[util_class.Plate_96Well(which_lib) for i in range(N_plates)] for j in range(N_primers)]

        try:
            plates = util_func._mutate_primers(plates, primers, primer_set, offset, constructs, which_lib)
            print('\033[92mSUCCESS\033[0m: Primerize Custom design() finished.\n')
        except:
            is_success = False
            print(traceback.format_exc())
            print('\033[41mERROR\033[0m: Primerize Custom design() encountered error.\n')

        params.update({'N_PLATE': N_plates, 'N_CONSTRUCT': N_constructs})
        data.update({'plates': plates, 'constructs': constructs})
        return Design_Plate({'sequence': sequence, 'name': name, 'is_success': is_success, 'primer_set': primer_set, 'params': params, 'data': data})

