from collections import namedtuple
import itertools
import numbers
import numpy as np
from narupa.protocol import trajectory

POSITIONS = 'particle.position'
ELEMENTS = 'particle.element'
TYPES = 'particle.type'
BONDS = 'bond'

# This dictionary matches the python types to the attributes of the GRPC
# values. This is not to do type conversion (which is handled by protobuf),
# but to figure out where to store the data.
PYTHON_TYPES_TO_GRPC_VALUE_ATTRIBUTE = {
    int: 'number_value', float: 'number_value', str: 'string_value',
    bool: 'bool_value', np.float32: 'number_value', np.float64: 'number_value',
}

_Shortcut = namedtuple(
    '_Shortcut', ['name', 'record_type', 'key', 'to_python', 'to_raw']
)


class MissingDataError(KeyError):
    """
    A shortcut does not contain data to return.
    """
    pass


def _as_is(value):
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
    """
    Metaclass that adds shortcuts to the :class:`FrameData` class.

    The shortcuts are defined as a tuple of :class:`_Shortcut` named tuples
    under the :attr:`_shortcuts` class attribute.
    """
    _shortcuts = ()

    def __init__(cls, name, bases, nmspc):
        super().__init__(name, bases, nmspc)
        for shortcut in cls._shortcuts:
            setattr(cls, shortcut.name, _make_shortcut(shortcut))


class FrameData(metaclass=_FrameDataMeta):
    """
    Wrapper around the GRPC FrameData.

    A ``FrameData`` contains two kinds of records: single values of any type,
    or homogeneous arrays. The former kind can be accessed through the
    :attr:`values` attribute, while the later is accessible through the
    :attr:`arrays` one. Both attribute behave like a dictionary. Trying to
    access a key that does not exist raises a :exc:`KeyError`.

    The most common frame properties are accessible as attribute in a
    normalized format. Shortcuts are not guaranteed to contain data. Trying to
    access a shortcut that does not contain data raises a
    :exc:`MissingDataError` that can also be caught as a :exc:`KeyError`.
    """
    _shortcuts = (
        _Shortcut(name='positions', key=POSITIONS,
                  record_type='arrays', to_python=_n_by_3, to_raw=_flatten_2d),
        _Shortcut(name='elements', key=ELEMENTS,
                  record_type='arrays', to_python=_as_is, to_raw=_as_is),
        _Shortcut(name='types', key=TYPES,
                  record_type='arrays', to_python=_as_is, to_raw=_as_is),
        _Shortcut(name='bonds', key=BONDS,
                  record_type='arrays', to_python=_n_by_2, to_raw=_flatten_2d),
    )

    def __init__(self, raw_frame=None):
        if raw_frame is None:
            self._raw = trajectory.FrameData()
        else:
            self._raw = raw_frame
        self.values = ValuesView(self.raw)
        self.arrays = ArraysView(self.raw)

    def __contains__(self, key):
        return key in self.arrays or key in self.values

    def __eq__(self, other):
        return self.raw == other.raw

    @property
    def raw(self):
        """
        Underlying GRPC/protobuf object.
        """
        # Use a property to make self.raw read-only.
        return self._raw


class RecordView:
    """
    Base class that wraps the access to a kind of record.

    This class needs to be subclassed.
    """
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
        """
        Extract the value from a protobuf field so it is usable by python.

        The method needs to be adapted to the type of field that is manipulated.
        """
        raise NotImplementedError('Subclasses must overwrite the _convert_to_python method.')


class ValuesView(RecordView):
    """
    Give access to singular values from a :class:`FrameData`.
    """
    record_name = 'values'
    singular = 'value'

    @staticmethod
    def _convert_to_python(field):
        return field.ListFields()[0][1]

    def set(self, key, value):
        type_attribute = PYTHON_TYPES_TO_GRPC_VALUE_ATTRIBUTE[type(value)]
        setattr(self._raw_record[key], type_attribute, value)


class ArraysView(RecordView):
    """
    Give access to homogeneous arrays from a :class:`FrameData`.
    """
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

        getattr(self._raw_record[key], type_attribute).values[:] = value
