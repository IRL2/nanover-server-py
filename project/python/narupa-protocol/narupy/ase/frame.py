from narupa.protocol.trajectory.frame_pb2 import FrameData
from ase import Atoms

AngToNm = 0.1

def ase_atoms_to_frame_data(ase_atoms: Atoms) -> FrameData:
    data = FrameData()
    array = data.arrays['atom.position'].float_values.values
    floats = ase_atoms.get_positions().flatten() * AngToNm
    array[:] = floats
    data.values['energy.potential'].number_value = ase_atoms.get_potential_energy()

    return data