from collections import namedtuple
import itertools
import numbers
from narupa.protocol import trajectory

POSITIONS = 'particle.position'
ELEMENTS = 'particle.element'
TYPES = 'particle.type'
BONDS = 'bond'

PYTHON_TYPES_TO_GRPC_VALUE_ATTRIBUTE = {
    int: 'number_value', float: 'number_value', str: 'string_value',
    bool: 'bool_value',
}

_Shortcut = namedtuple(
    '_Shortcut', ['name', 'record_type', 'key', 'to_python', 'to_raw']
)


class MissingDataError(KeyError):
    pass


def _as_it(value):
    return value


def _n_by_2(value):
    return list(value[i:i + 2] for i in range(0, len(value), 2))


def _n_by_3(value):
    return list(value[i:i + 3] for i in range(0, len(value), 3))


def _flatten_2d(value):
    return list(itertools.chain(*value))


def _make_getter(shortcut):
    def wrapped(self):
        try:
            value = getattr(self, shortcut.record_type)[shortcut.key]
        except KeyError as error:
            raise MissingDataError(str(error))
        return shortcut.to_python(value)

    return wrapped


def _make_setter(shortcut):
    def wrapped(self, value):
        converted_value = shortcut.to_raw(value)
        getattr(self, shortcut.record_type)[shortcut.key] = converted_value

    return wrapped


def _make_shortcut(shortcut):
    return property(fget=_make_getter(shortcut), fset=_make_setter(shortcut))


class _FrameDataMeta(type):
    _shortcuts = ()

    def __init__(cls, name, bases, nmspc):
        super().__init__(name, bases, nmspc)
        for shortcut in cls._shortcuts:
            setattr(cls, shortcut.name, _make_shortcut(shortcut))


class FrameData(metaclass=_FrameDataMeta):
    _shortcuts = (
        _Shortcut(name='positions', key=POSITIONS,
                  record_type='arrays', to_python=_n_by_3, to_raw=_flatten_2d),
        _Shortcut(name='elements', key=ELEMENTS,
                  record_type='arrays', to_python=_as_it, to_raw=_as_it),
        _Shortcut(name='types', key=TYPES,
                  record_type='arrays', to_python=_as_it, to_raw=_as_it),
        _Shortcut(name='bonds', key=BONDS,
                  record_type='arrays', to_python=_n_by_2, to_raw=_flatten_2d),
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

    def __eq__(self, other):
        return self.raw == other.raw


class RecordView:
    record_name = None  # MUST be overwritten as "arrays" or "values"
    singular = None  # MUST be overwritten as "array" or "value"

    def __init__(self, raw):
        if self.record_name is None or self.singular is None:
            raise NotImplementedError(
                'FieldView must be subclassed; record_name, singular, and'
                '_convert_to_python must be overwritten.'
            )
        self._raw_record = getattr(raw, self.record_name)

    def __getitem__(self, key):
        if key in self:
            field = self._raw_record[key]
            return self._convert_to_python(field)
        raise KeyError(f'No {self.singular} with the key "{key}".')

    def __setitem__(self, key, value):
        self.set(key, value)

    def __contains__(self, key):
        return key in self._raw_record

    def get(self, key, default=None):
        if key in self:
            return self[key]
        return default

    def set(self, key, value):
        raise NotImplementedError('Subclasses must overwrite the set method.')

    @staticmethod
    def _convert_to_python(field):
        raise NotImplementedError('Subclasses must overwrite the _convert_to_python method.')


class ValuesView(RecordView):
    record_name = 'values'
    singular = 'value'

    @staticmethod
    def _convert_to_python(field):
        return field.ListFields()[0][1]

    def set(self, key, value):
        type_attribute = PYTHON_TYPES_TO_GRPC_VALUE_ATTRIBUTE[type(value)]
        setattr(self._raw_record[key], type_attribute, value)


class ArraysView(RecordView):
    record_name = 'arrays'
    singular = 'array'

    @staticmethod
    def _convert_to_python(field):
        return field.ListFields()[0][1].values

    def set(self, key, value):
        try:
            reference_value = value[0]
        except IndexError:
            raise ValueError('Cannot decide what type to use for an empty object.')
        except TypeError:
            raise ValueError('Value must be indexable.')

        if isinstance(reference_value, numbers.Integral) and reference_value >= 0:
            type_attribute = 'index_values'
        elif isinstance(reference_value, numbers.Real):
            type_attribute = 'float_values'
        elif isinstance(reference_value, str):
            type_attribute = 'string_values'
        else:
            raise ValueError('Cannot decide what type to use.')

        if key in self:
            while self[key]:
                self[key].pop()

        getattr(self._raw_record[key], type_attribute).values.extend(value)
