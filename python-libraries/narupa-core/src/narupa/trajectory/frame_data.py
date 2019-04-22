from narupa.protocol import trajectory


class FrameData:
    def __init__(self, raw_frame=None):
        if raw_frame is None:
            self.raw = trajectory.FrameData()
        else:
            self.raw = raw_frame
        self.values = ValuesView(self.raw)
        self.arrays = ArraysView(self.raw)

    def __contains__(self, key):
        return key in self.arrays or key in self.values


class RecordView:
    record_name = None  # MUST be overwritten as "arrays" or "values"
    singular = None  # MUST be overwritten as "array" or "value"

    def __init__(self, raw):
        if self.record_name is None or self.singular is None:
            raise NotImplementedError(
                'FieldView must be subclassed; record_name, singular, and'
                '_converter must be overwritten.'
            )
        self._raw_field = getattr(raw, self.record_name)

    def __getitem__(self, key):
        if key in self:
            field = self._raw_field[key]
            return self._converter(field)
        raise KeyError(f'No {self.singular} with the key "{key}".')

    def __contains__(self, key):
        return key in self._raw_field

    def get(self, key, default=None):
        if key in self:
            return self[key]
        return default

    @staticmethod
    def _converter(field):
        raise NotImplementedError('Subclasses must overwrite the _converter method.')


class ValuesView(RecordView):
    record_name = 'values'
    singular = 'value'

    @staticmethod
    def _converter(field):
        return field.ListFields()[0][1]


class ArraysView(RecordView):
    record_name = 'arrays'
    singular = 'array'

    @staticmethod
    def _converter(field):
        return field.ListFields()[0][1].values
