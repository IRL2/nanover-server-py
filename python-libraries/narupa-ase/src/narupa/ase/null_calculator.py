import numpy as np
from ase.calculators.calculator import Calculator, all_changes


class NullCalculator(Calculator):
    """
    An empty ASE calculator.
    """

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.implemented_properties = ('forces', 'energy')

    def calculate(self, atoms=None, properties=('forces', 'energy'), system_changes=all_changes):
        self.results = {'forces': np.zeros((len(atoms), 3)), 'energies': 0.0}
