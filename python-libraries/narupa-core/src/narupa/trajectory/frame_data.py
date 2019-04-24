from collections import namedtuple
from narupa.protocol import trajectory

POSITIONS = 'particle.position'
ELEMENTS = 'particle.element'
TYPES = 'particle.type'
BONDS = 'bond'


_Shortcut = namedtuple('_Shortcut', ['name', 'record_type', 'key', 'converter'])


class MissingDataError(KeyError):
    pass


def _as_it(value):
    return value


def _n_by_2(value):
    return list(value[i:i + 2] for i in range(0, len(value), 2))


def _n_by_3(value):
    return list(value[i:i + 3] for i in range(0, len(value), 3))


def _make_getter(shortcut):
    def wrapped(self):
        try:
            value = getattr(self, shortcut.record_type)[shortcut.key]
        except KeyError as error:
            raise MissingDataError(str(error))
        return shortcut.converter(value)

    return wrapped


def _make_shortcut(shortcut):
    return property(fget=_make_getter(shortcut))


class _FrameDataMeta(type):
    _shortcuts = ()

    def __init__(cls, name, bases, nmspc):
        super().__init__(name, bases, nmspc)
        for shortcut in cls._shortcuts:
            setattr(cls, shortcut.name, _make_shortcut(shortcut))


class FrameData(metaclass=_FrameDataMeta):
    _shortcuts = (
        _Shortcut(name='positions', key=POSITIONS,
                  record_type='arrays', converter=_n_by_3),
        _Shortcut(name='elements', key=ELEMENTS,
                  record_type='arrays', converter=_as_it),
        _Shortcut(name='types', key=TYPES,
                  record_type='arrays', converter=_as_it),
        _Shortcut(name='bonds', key=BONDS,
                  record_type='arrays', converter=_n_by_2),
    )

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



