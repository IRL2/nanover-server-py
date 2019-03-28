from narupa.protocol.trajectory.frame_pb2 import FrameData

def openmm_positions_to_frame_data(positions) -> FrameData:
    data = FrameData()

    array = data.arrays['atom.position'].float_values.values

    floats = [value for position in positions for value in position._value]
    array.extend(floats)

    return data