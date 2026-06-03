from .thermo import Nearest_Neighbor, Singleton
from .util_class import Assembly, Plate_96Well, Mutation, Construct_List
from .wrapper import Design_Single, Design_Plate

from .primerize_1d import Primerize_1D
from .primerize_2d import Primerize_2D
from .primerize_3d import Primerize_3D
from .primerize_custom import Primerize_Custom

__version__ = '2.0.0'


Primerize_1D = Primerize_1D()
Primerize_2D = Primerize_2D()
Primerize_3D = Primerize_3D()
Primerize_Custom = Primerize_Custom()

__all__ = ['Primerize_1D', 'Primerize_2D', 'Primerize_3D', 'Primerize_Custom', 'Design_Single', 'Design_Plate', 'Nearest_Neighbor', 'Singleton', 'Assembly', 'Plate_96Well', 'Mutation', 'Construct_List', 'misprime', 'util']
