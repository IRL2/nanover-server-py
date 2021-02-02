# Copyright (c) Intangible Realities Lab, University Of Bristol. All rights reserved.
# Licensed under the GPL. See License.txt in the project root for license information.
from collections import namedtuple
from collections.abc import Set
import numbers
from typing import Dict, Optional, List

import numpy as np
from narupa.protocol import trajectory
from narupa.utilities.protobuf_utilities import value_to_object, object_to_value

BOX_VECTORS = 'system.box.vectors'

BOND_PAIRS = 'bond.pairs'
BOND_ORDERS = 'bond.orders'

PARTICLE_POSITIONS = 'particle.positions'
PARTICLE_ELEMENTS = 'particle.elements'
PARTICLE_TYPES = 'particle.types'
PARTICLE_NAMES = 'particle.names'
PARTICLE_RESIDUES = 'particle.residues'  # Index of the residue each particle belongs to.
PARTICLE_COUNT = 'particle.count'

RESIDUE_NAMES = 'residue.names'
RESIDUE_IDS = 'residue.ids'  # Index of the chain each residue belongs to.
RESIDUE_CHAINS = 'residue.chains'
RESIDUE_COUNT = 'residue.count'

CHAIN_NAMES = 'chain.names'
CHAIN_COUNT = 'chain.count'

KINETIC_ENERGY = 'energy.kinetic'
POTENTIAL_ENERGY = 'energy.potential'

SERVER_TIMESTAMP = 'server.timestamp'


_Shortcut = namedtuple(
    '_Shortcut', ['record_type', 'key', 'field_type', 'to_python', 'to_raw']
)


class MissingDataError(KeyError):
    """
    A shortcut does not contain data to return.
    """
    pass


def _as_is(value):
    return value


def _as_int(value):
    return int(value)


def _n_by_2(value):
    return list(value[i:i + 2] for i in range(0, len(value), 2))


def _n_by_3(value):
    return list(value[i:i + 3] for i in range(0, len(value), 3))


def _flatten_array(value):
    return np.asarray(value).ravel()


def _make_getter(shortcut):
    def wrapped(self):
        try:
            value = getattr(self, shortcut.record_type)[shortcut.key]
        except KeyError as error:
            raise MissingDataError(str(error))
        return shortcut.to_python(value)

    return wrapped


def _make_setter(shortcut):
    if shortcut.record_type == 'arrays' and shortcut.field_type is not None:
        method_name = f'set_{shortcut.field_type}_array'

        def wrapped(self, value):
            converted_value = shortcut.to_raw(value)
            getattr(self, method_name)(shortcut.key, converted_value)

    else:

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
    _shortcuts: Dict[str, _Shortcut] = {}

    def __init__(cls, name, bases, nmspc):
        shortcuts = {}
        super().__init__(name, bases, nmspc)
        for attribute_name, attribute in nmspc.items():
            if isinstance(attribute, _Shortcut):
                shortcuts[attribute_name] = attribute
                setattr(cls, attribute_name, _make_shortcut(attribute))
        cls._shortcuts = shortcuts


class FrameData(metaclass=_FrameDataMeta):
    """
    Wrapper around the GRPC FrameData.

    A ``FrameData`` contains two kinds of records: single values of any type,
    or homogeneous arrays. The former kind can be accessed through the
    :attr:`values` attribute, while the later is accessible through the
    :attr:`arrays` one. Both attribute behave like a dictionary. Trying to
    access a key that does not exist raises a :exc:`KeyError`.

    The set of keys with data in the frame is listed by :meth:`value_keys`
    and :meth:`array_keys`.

    The most common frame properties are accessible as attribute in a
    normalized format. Shortcuts are not guaranteed to contain data. Trying to
    access a shortcut that does not contain data raises a
    :exc:`MissingDataError` that can also be caught as a :exc:`KeyError`.

    The available shortcuts can be listed using the :attr:`shortcuts` property.
    The set of shortcuts that contain data is available from the
    :attr:`used_shortcuts`.
    """
    bond_pairs: List[List[int]] = _Shortcut(  # type: ignore[assignment]
        key=BOND_PAIRS, record_type='arrays',
        field_type='index', to_python=_n_by_2, to_raw=_flatten_array)
    bond_orders: List[float] = _Shortcut(  # type: ignore[assignment]
        key=BOND_ORDERS, record_type='arrays',
        field_type='float', to_python=_as_is, to_raw=_as_is)

    particle_positions: List[List[float]] = _Shortcut(  # type: ignore[assignment]
        key=PARTICLE_POSITIONS, record_type='arrays',
        field_type='float', to_python=_n_by_3, to_raw=_flatten_array)
    particle_elements: List[int] = _Shortcut(  # type: ignore[assignment]
        key=PARTICLE_ELEMENTS, record_type='arrays',
        field_type='index', to_python=_as_is, to_raw=_as_is)
    particle_types: List[str] = _Shortcut(  # type: ignore[assignment]
        key=PARTICLE_TYPES, record_type='arrays',
        field_type='string', to_python=_as_is, to_raw=_as_is)
    particle_names: List[str] = _Shortcut(  # type: ignore[assignment]
        key=PARTICLE_NAMES, record_type='arrays',
        field_type='string', to_python=_as_is, to_raw=_as_is)
    particle_residues: List[int] = _Shortcut(  # type: ignore[assignment]
        key=PARTICLE_RESIDUES, record_type='arrays',
        field_type='index', to_python=_as_is, to_raw=_as_is)
    particle_count: int = _Shortcut(  # type: ignore[assignment]
        key=PARTICLE_COUNT, record_type='values',
        field_type='number_value', to_python=_as_int, to_raw=_as_is)

    residue_names: List[str] = _Shortcut(  # type: ignore[assignment]
        key=RESIDUE_NAMES, record_type='arrays',
        field_type='string', to_python=_as_is, to_raw=_as_is)
    residue_ids: List[int] = _Shortcut(  # type: ignore[assignment]
        key=RESIDUE_IDS, record_type='arrays',
        field_type='string', to_python=_as_is, to_raw=_as_is)
    residue_chains: List[int] = _Shortcut(  # type: ignore[assignment]
        key=RESIDUE_CHAINS, record_type='arrays',
        field_type='index', to_python=_as_is, to_raw=_as_is)
    residue_count: int = _Shortcut(  # type: ignore[assignment]
        key=RESIDUE_COUNT, record_type='values',
        field_type='number_value', to_python=_as_int, to_raw=_as_is)

    chain_names: List[str] = _Shortcut(  # type: ignore[assignment]
        key=CHAIN_NAMES, record_type='arrays',
        field_type='string', to_python=_as_is, to_raw=_as_is)
    chain_count: int = _Shortcut(  # type: ignore[assignment]
        key=CHAIN_COUNT, record_type='values',
        field_type='number_value', to_python=_as_int, to_raw=_as_is)

    kinetic_energy: float = _Shortcut(  # type: ignore[assignment]
        key=KINETIC_ENERGY, record_type='values',
        field_type='number_value', to_python=_as_is, to_raw=_as_is)
    potential_energy: float = _Shortcut(  # type: ignore[assignment]
        key=POTENTIAL_ENERGY, record_type='values',
        field_type='number_value', to_python=_as_is, to_raw=_as_is)
    box_vectors: List[List[float]] = _Shortcut(  # type: ignore[assignment]
        key=BOX_VECTORS, record_type='arrays',
        field_type='float', to_python=_n_by_3, to_raw=_flatten_array)

    server_timestamp: float = _Shortcut(  # type: ignore[assignment]
        key=SERVER_TIMESTAMP, record_type='values',
        field_type='number_value', to_python=_as_is, to_raw=_as_is)

    _shortcuts: Dict[str, _Shortcut]
    _raw: trajectory.FrameData

    def __init__(self, raw_frame: trajectory.FrameData = None):
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

    def __repr__(self):
        return repr(self.raw)

    def __delattr__(self, attr):
        if attr in self._shortcuts:
            shortcut = self._shortcuts[attr]
            getattr(self, shortcut.record_type).delete(shortcut.key)
        else:
            super().__delattr__(attr)

    def __delitem__(self, item):
        if item in self.value_keys:
            del self.values[item]
        if item in self.array_keys:
            del self.arrays[item]

    def copy(self):
        copy = FrameData()
        for key in self.value_keys:
            copy.values.set(key, self.values[key])
        for key in self.array_keys:
            copy.arrays.set(key, self.arrays[key])
        return copy

    @property
    def raw(self) -> trajectory.FrameData:
        """
        Underlying GRPC/protobuf object.
        """
        # Use a property to make self.raw read-only.
        return self._raw

    # Methods to match the C# API
    def set_float_array(self, key, value):
        """
        Set an homogeneous array of floats in an existing or a new key.

        :param key: The key under which to store the array.
        :param value: The array to store.
        """
        self.raw.arrays[key].float_values.values[:] = value

    def set_index_array(self, key, value):
        """
        Set an homogeneous array of indices in an existing or a new key.

        :param key: The key under which to store the array.
        :param value: The array to store.
        """
        self.raw.arrays[key].index_values.values[:] = value

    def set_string_array(self, key, value):
        """
        Set an homogeneous array of strings in an existing or a new key.

        :param key: The key under which to store the array.
        :param value: The array to store.
        """
        self.raw.arrays[key].string_values.values[:] = value

    @property
    def value_keys(self) -> Set:
        return self.values.keys()

    @property
    def array_keys(self) -> Set:
        return self.arrays.keys()

    @property
    def shortcuts(self) -> Set:
        return set(self._shortcuts.keys())

    @property
    def used_shortcuts(self) -> Set:
        return set(
            name for name, shortcut in self._shortcuts.items()
            if shortcut.key in self
        )


class RecordView:
    """
    Base class that wraps the access to a kind of record.

    This class needs to be subclassed.
    """
    record_name: Optional[str] = None  # MUST be overwritten as "arrays" or "values"
    singular: Optional[str] = None  # MUST be overwritten as "array" or "value"

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

    def __delitem__(self, key):
        del self._raw_record[key]

    def __contains__(self, key):
        return key in self._raw_record

    def get(self, key, default=None):
        if key in self:
            return self[key]
        return default

    def set(self, key, value):
        raise NotImplementedError('Subclasses must overwrite the set method.')

    def delete(self, key):
        del self[key]

    @staticmethod
    def _convert_to_python(field):
        """
        Extract the value from a protobuf field so it is usable by python.

        The method needs to be adapted to the type of field that is manipulated.
        """
        raise NotImplementedError('Subclasses must overwrite the _convert_to_python method.')

    def keys(self) -> Set:
        return set(self._raw_record.keys())


class ValuesView(RecordView):
    """
    Give access to singular values from a :class:`FrameData`.
    """
    record_name = 'values'
    singular = 'value'

    @staticmethod
    def _convert_to_python(field):
        return value_to_object(field)

    def set(self, key, value):
        self._raw_record[key].CopyFrom(object_to_value(value))


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
