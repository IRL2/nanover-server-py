from narupa.protocol import trajectory


class FrameData:
    def __init__(self, raw_frame=None):
        if raw_frame is None:
            self.raw = trajectory.FrameData()
        else:
            self.raw = raw_frame

    def __getitem__(self, key):
        if key in self.raw.arrays:
            field = self.raw.arrays[key]
            return _convert_array_to_python(field)
        if key in self.raw.values:
            field = self.raw.values[key]
            return _convert_value_to_python(field)
        raise KeyError(f'No "{key}" key in the f{self.__class__.__name__}.')


def _convert_value_to_python(field):
    return field.ListFields()[0][1]


def _convert_array_to_python(field):
    return field.ListFields()[0][1].values
