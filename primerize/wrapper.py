import math
import os

from . import util_func
from . import util_server


class Design_Single(object):
    """Result of a ``primerize.Primerize_1D.design()`` run.

    Args:
        init_dict: A ``dict`` with the following keys:

        sequence: ``str``: Sequence of assembly design.
        name: ``str``: Construct prefix/name.
        is_success: ``bool``: Flag for whether ``primerize.Primerize_1D.design()`` run successfully found a solution.
        primer_set: ``list(str)``: List of primers for assembly.
        params: ``dict``: Dictionary of parameters used for this result.
        data: ``dict``: Dictionary of result data.

    Attributes:
        sequence: ``str``: Sequence of assembly design.
        name: ``str``: Construct prefix/name.
        is_success: ``bool``: Flag for whether a solution is found.
        primer_set: ``list(str)``: Strings of solution primers.
        _params: Input parameters, in format of ``dict: { 'MIN_TM': float, 'NUM_PRIMERS': int, 'MIN_LENGTH': int, 'MAX_LENGTH': int, 'N_BP': int, 'COL_SIZE': int, 'WARN_CUTOFF': int }``.
        _data: Data of assembly solution, in format of ::

                dict: {
                    'misprime_score': [str, str],
                    'warnings': list(tuple(int, int, int, int)),
                    'assembly': primerize.Assembly
                }
    """

    def __init__(self, init_dict):
        for key in init_dict:
            if key not in ['sequence', 'name', 'is_success', 'primer_set', 'params', 'data']:
                raise ValueError('\033[41mERROR\033[0m: Unrecognized key \033[92m%s\033[0m for \033[94m%s\033[0m.\n' % (key, self.__class__))
            key_rename = '_' + key if key in ['params', 'data'] else key
            setattr(self, key_rename, init_dict[key])

    def __repr__(self):
        """Representation of the ``Design_Single`` class.
        """

        return '\033[94m%s\033[0m {\n\033[95msequence\033[0m = \'%s\', \n\033[95mname\033[0m = \'%s\', \n\033[95mis_success\033[0m = \033[41m%s\033[0m, \n\033[95mprimer_set\033[0m = %s, \n\033[95mparams\033[0m = %s, \n\033[95mdata\033[0m = {\n    \033[92m\'misprime_score\'\033[0m: %s, \n    \033[92m\'assembly\'\033[0m: %s, \n    \033[92m\'warnings\'\033[0m: %s\n}' % (self.__class__, self.sequence, self.name, self.is_success, repr(self.primer_set), repr(self._params), repr(self._data['misprime_score']), repr(self._data['assembly']), repr(self._data['warnings']))

    def __str__(self):
        """Results of the ``Design_Single`` class. Calls ``echo()``.
        """

        return self.echo()


    def get(self, key):
        """Get result parameters.

        Args:
            key: ``str``: Keyword of parameter. Valid keywords are ``'MIN_TM'``, ``'NUM_PRIMERS'``, ``'MIN_LENGTH'``, ``'MAX_LENGTH'``, ``'COL_SIZE'``, ``'WARN_CUTOFF'``, ``'WARNING'``, ``'PRIMER'``, ``'MISPRIME'``; case insensitive.

        Returns:
            value of specified **key**.

        Raises:
            AttributeError: For illegal **key**.
        """

        key = key.upper()
        if key in self._params:
            return self._params[key]
        elif key == 'WARNING':
            return self._data['warnings']
        elif key == 'PRIMER':
            return self._data['asssembly'].primers
        elif key == 'MISPRIME':
            return self._data['misprime_score']
        else:
            raise AttributeError('\033[41mERROR\033[0m: Unrecognized key \033[92m%s\033[0m for \033[94m%s.get()\033[0m.\n' % (key, self.__class__))


    def save(self, path='./', name=None):
        """Save result to text file.

        Args:
            path: ``str``: `(Optional)` Path for file saving. Use either relative or absolute path.
            name: ``str``: `(Optional)` Prefix/name for file name. When nonspecified, current object's name is used.
        """

        if self.is_success:
            name = self.name if name is None else name
            f = open(os.path.join(path, '%s.txt' % name), 'w')

            f.write('Primerize Result\n\nINPUT\n=====\n%s\n' % self.sequence)
            f.write('#\nMIN_TM: %.1f\n' % self._params['MIN_TM'])
            if not self._params['NUM_PRIMERS']:
                f.write('NUM_PRIMERS: auto (unspecified)')
            else:
                f.write('NUM_PRIMERS: %d' % self._params['NUM_PRIMERS'])
            f.write('\nMAX_LENGTH: %d\nMIN_LENGTH: %d\n' % (self._params['MAX_LENGTH'], self._params['MIN_LENGTH']))

            f.write('\n\nOUTPUT\n======\n')
            lines = str(self).replace('\033[0m', '').replace('\033[100m', '').replace('\033[92m', '').replace('\033[93m', '').replace('\033[94m', '').replace('\033[95m', '').replace('\033[96m', '').replace('\033[41m', '')
            f.write(lines)
            f.write('#\n\n------/* IDT USER: for primer ordering, copy and paste to Bulk Input */------\n------/* START */------\n')
            for i in range(len(self.primer_set)):
                suffix = 'FR'[i % 2]
                f.write('%s-%d%s\t%s\t\t25nm\tSTD\n' % (self.name, i + 1, suffix, self.primer_set[i]))
            f.write('------/* END */------\n------/* NOTE: use "Lab Ready" for "Normalization" */------\n')
            f.close()
        else:
            raise UnboundLocalError('\033[41mFAIL\033[0m: Result unavailable for \033[94m%s\033[0m where \033[94mis_success\033[0m = \033[41mFalse\033[0m.\n' % self.__class__)


    def echo(self, key=''):
        """Print part(s) of result in rich-text.

        Args:
            key: ``str``: `(Optional)` Keyword of printing. Valid keywords are ``'misprime'``, ``'warning'``, ``'primer'``, ``'assembly'``; case insensitive. When nonspecified, result of all keywords is returned.

        Returns:
            ``str``

        Raises:
            AttributeError: For illegal **key**.
            UnboundLocalError: When ``is_success = False``.
        """

        if self.is_success:
            key = key.lower()
            if key == 'misprime':
                output = ''
                for i in range(int(math.floor(self._params['N_BP'] / self._params['COL_SIZE'])) + 1):
                    output += '%s\n\033[92m%s\033[0m\n%s\n\n' % (self._data['misprime_score'][0][i * self._params['COL_SIZE']:(i + 1) * self._params['COL_SIZE']], self.sequence[i * self._params['COL_SIZE']:(i + 1) * self._params['COL_SIZE']], self._data['misprime_score'][1][i * self._params['COL_SIZE']:(i + 1) * self._params['COL_SIZE']])
                return output[:-1]

            elif key == 'warning':
                output = ''
                for warning in self._data['warnings']:
                    p_1 = '\033[100m%d\033[0m%s' % (warning[0], util_func._primer_suffix(warning[0] - 1))
                    p_2 = ', '.join('\033[100m%d\033[0m%s' % (x, util_func._primer_suffix(x - 1)) for x in warning[3])
                    output += '\033[93mWARNING\033[0m: Primer %s can misprime with %d-residue overlap to position %s, which is covered by primers: %s\n' % (p_1.rjust(4), warning[1], str(int(warning[2])).rjust(3), p_2)
                return output[:-1]

            elif key == 'primer':
                output = '%s%s\tSEQUENCE\n' % ('PRIMERS'.ljust(20), 'LENGTH'.ljust(10))
                for i, primer in enumerate(self.primer_set):
                    name = '%s-\033[100m%s\033[0m%s' % (self.name, i + 1, util_func._primer_suffix(i))
                    output += '%s\033[93m%s\033[0m\t%s\n' % (name.ljust(39), str(len(primer)).ljust(10), util_func._primer_suffix(i).replace(' R', primer).replace(' F', primer))
                return output[:-1]

            elif key == 'assembly':
                return self._data['assembly'].echo()
            elif not key:
                return self.echo('misprime') + '\n' + self.echo('assembly') + '\n' + self.echo('primer') + '\n\n' + self.echo('warning') + '\n'

            else:
                raise AttributeError('\033[41mERROR\033[0m: Unrecognized key \033[92m%s\033[0m for \033[94m%s.echo()\033[0m.\n' % (key, self.__class__))
        else:
            raise UnboundLocalError('\033[41mFAIL\033[0m: Result of key \033[92m%s\033[0m unavailable for \033[94m%s\033[0m where \033[94mis_cucess\033[0m = \033[41mFalse\033[0m.\n' % (key, self.__class__))



class Design_Plate(object):
    """Result of a ``primerize.Primerize_2D.design()`` or ``primerize.Primerize_3D.design()`` run.

    Args:
        init_dict: A ``dict`` with the following keys:
        sequence: ``str``: Sequence of assembly design.
        name: ``str``: Construct prefix/name.
        is_success: ``bool``: Flag for whether ``primerize.Primerize_2D.design()`` or ``primerize.Primerize_3D.design()`` run successfully found a solution.
        primer_set: ``list(str)``: List of primers for assembly.
        structures: ``list(str)``: `(Optional)` List of secondary structures, only `Required` for ``primerize.Primerize_3D.design()`` results.
        params: ``dict``: Dictionary of parameters used for this result.
        data: ``dict``: Dictionary of result data.

    Attributes:
        sequence: ``str``: Sequence of assembly design.
        name: ``str``: Construct prefix/name.
        is_success: ``bool``: Flag for whether a solution is found.
        primer_set: ``list(str)``: Strings of solution primers.
        structures: ``list(str)``: Strings of input secondary structures.
        _params: Input parameters, in format of ``dict: { 'offset': int, 'which_muts': list(int), 'which_lib': list(int), 'N_PRIMER': int, 'N_PLATE': int, 'N_CONSTRUCT': int, 'N_BP': int, 'type': str }``.

            For ``primerize.Primerize_3D.design()`` results, it also has ``'N_MUTATION': int, 'is_exclude': bool, 'is_single': bool, 'is_fillWT': bool``.

        _data: Data of assembly solution, in format of ::

                dict: {
                    'constructs': primerize.Construct_List,
                    'plates': list(list(primerize.Plate_96Well)),
                    'assembly': primerize.Assembly,
                    'illustration': { 'labels': list(str), 'fragments': list(str), 'lines': tuple(str) },
                }

            For ``primerize.Primerize_3D.design()`` results, it also has ``'bps': list(tuple(int, int)), 'warnings': list(tuple(int, int))``.
    """

    def __init__(self, init_dict):
        for key in init_dict:
            if key not in ['sequence', 'name', 'is_success', 'primer_set', 'structures', 'params', 'data']:
                raise ValueError('\033[41mERROR\033[0m: Unrecognized key \033[92m%s\033[0m for \033[94m%s\033[0m.\n' % (key, self.__class__))
            key_rename = '_' + key if key in ['params', 'data'] else key
            setattr(self, key_rename, init_dict[key])

        if self.get('TYPE') == 'Mutate-and-Map':
            self._data['illustration'] = util_server._draw_region(self.sequence, self._params)
        elif self.get('TYPE') == 'Mutation/Rescue':
            self._data['illustration'] = util_server._draw_str_region(self.sequence, self.structures, self._data['bps'], self._data['warnings'], self._params)
            print(util_func._print_pair_mismatch_warning(self.sequence, self._data['warnings'], self._params['offset']))
        else:
            self._data['illustration'] = {'lines': ''}

    def __repr__(self):
        """Representation of the ``Design_Plate`` class.
        """

        structures = '\033[95mstructures\033[0m = %s, \n' % repr(self.structures) if self.get('TYPE') == 'Mutation/Rescue' else ''
        return '\033[94m%s\033[0m {\n\033[95msequence\033[0m = \'%s\', \n\033[95mname\033[0m = \'%s\', \n\033[95mis_success\033[0m = \033[41m%s\033[0m, \n\033[95mprimer_set\033[0m = %s, \n%s\033[95mparams\033[0m = %s, \n\033[95mdata\033[0m = {\n    \033[92m\'constructs\'\033[0m: %s, \n    \033[92m\'assembly\'\033[0m: %s, \n    \033[92m\'plates\'\033[0m: %s\n}' % (self.__class__, self.sequence, self.name, self.is_success, repr(self.primer_set), structures, repr(self._params), repr(self._data['constructs']), repr(self._data['assembly']), repr(self._data['plates']))

    def __str__(self):
        """Results of the ``Design_Plate`` class. Calls ``echo()``.
        """

        return self.echo()


    def get(self, key):
        """Get result parameters.

        Args:
            key: ``str``: Keyword of parameter. Valid keywords are ``'offset'``, ``'which_muts'``, ``'which_lib'``, ``'N_PRIMER'``, ``'N_PLATE'``, ``'N_CONSTRUCT'``, ``'N_BP'``, ``'PRIMER'``, ``'CONSTRUCT'``, (``'is_exclude'``, ``'is_signle'``, ``'is_fillWT'``, ``'STRUCTURE'`` and ``'WARNING'`` only for ``primerize.Primerize_3D.design()`` results); case insensitive.

        Returns:
            value of specified **key**.

        Raises:
            AttributeError: For illegal **key**.
        """

        key = key.upper()
        if key in self._params:
            return self._params[key]
        elif key.lower() in self._params:
            return self._params[key.lower()]
        elif key == 'PRIMER':
            return self._data['assembly'].primers
        elif key == 'CONSTRUCT':
            return self._data['constructs']
        elif key == 'STRUCTURE' and self.get('TYPE') == 'Mutation/Rescue':
            return self.structures
        elif key == 'WARNING' and self.get('TYPE') == 'Mutation/Rescue':
            return self._data['warnings']
        else:
            raise AttributeError('\033[41mERROR\033[0m: Unrecognized key \033[92m%s\033[0m for \033[94m%s.get()\033[0m.\n' % (key, self.__class__))


    def save(self, key='', path='./', name=None):
        """Save result to text file.

        Args:
            key: ``str``: `(Optional)` Keyword of saving. Valid keywords are ``'table'``, ``'image'``, ``'constructs'``, ``'assembly'``, (``'structures'`` only for ``primerize.Primerize_3D.design()`` results); case insensitive. When nonspecified, files of all keywords are saved.
            path: ``str``: `(Optional)` Path for file saving. Use either relative or absolute path.
            name: ``str``: `(Optional)` Prefix/name for file name. When nonspecified, current object's name is used.

        Raises:
            AttributeError: For illegal **key**.
            UnboundLocalError: When ``is_success = False``.
        """

        if self.is_success:
            name = self.name if name is None else name
            key = key.lower()
            if key == 'table':
                util_func._save_plates_excel(self._data['plates'], self.primer_set, name, path)
            elif key == 'image':
                util_func._save_plate_layout(self._data['plates'], self.primer_set, name, path)
            elif key == 'constructs':
                util_func._save_construct_key(self._data['constructs'], name, path, self._params['which_lib'])
            elif key == 'assembly':
                self._data['assembly'].save(path, name)
            elif key == 'structures' and self.get('TYPE') == 'Mutation/Rescue':
                util_func._save_structures(self.structures, self._data['warnings'], self.sequence, self._params['offset'], name, path)

            elif not key:
                keys = ['table', 'image', 'constructs', 'assembly']
                if self.get('TYPE') == 'Mutation/Rescue': keys.append('structures')
                for key in keys:
                    self.save(key, path, name)
            else:
                raise AttributeError('\033[41mERROR\033[0m: Unrecognized key \033[92m%s\033[0m for \033[94m%s.save()\033[0m.\n' % (key, self.__class__))
        else:
            raise UnboundLocalError('\033[41mFAIL\033[0m: Result of key \033[92m%s\033[0m unavailable for \033[94m%s\033[0m where \033[94mis_cucess\033[0m = \033[41mFalse\033[0m.\n' % (key, self.__class__))


    def echo(self, key=''):
        """Print part(s) of result in rich-text.

        Args:
            key: ``str``: `(Optional)` Keyword of printing. Valid keywords are ``'plate'``, ``'assembly'``, ``'region'``; case insensitive. When nonspecified, result of all keywords is returned.

        Return:
            ``str``

        Raises:
            AttributeError: For illegal **key**.
            UnboundLocalError: When ``is_success = False``.
        """

        if self.is_success:
            key = key.lower()
            if key == 'plate':
                output = ''
                for i in range(len(self._data['plates'][0])):
                    for j in range(len(self._data['plates'])):
                        output += 'Plate \033[95m%d\033[0m; Primer \033[92m%d\033[0m\n' % (i + 1, j + 1)
                        output += self._data['plates'][j][i].echo(self.primer_set[j])
                return output[:-1]
            elif key == 'assembly':
                return self._data['assembly'].echo()
            elif key == 'region':
                structures = '\n'.join(self._data['illustration']['lines'])
                if self.get('TYPE') == 'Mutation/Rescue':
                    warnings = util_func._print_pair_mismatch_warning(self.sequence, self._data['warnings'], self._params['offset'])
                    return structures if not warnings else '%s\n\n%s\n' % (structures, warnings)
                else:
                    return structures

            elif not key:
                return self.echo('assembly') + '\n\n' + self.echo('plate') + '\n\n' + self.echo('region')
            else:
                raise AttributeError('\033[41mERROR\033[0m: Unrecognized key \033[92m%s\033[0m for \033[94m%s.echo()\033[0m.\n' % (key, self.__class__))
        else:
            raise UnboundLocalError('\033[41mFAIL\033[0m: Result of key \033[92m%s\033[0m unavailable for \033[94m%s\033[0m where \033[94mis_cucess\033[0m = \033[41mFalse\033[0m.\n' % (key, self.__class__))


