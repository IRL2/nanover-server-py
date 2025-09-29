from collections.abc import Set
import numbers
from .protocol import trajectory
from nanover.utilities.protobuf_utilities import value_to_object, object_to_value


class FrameData:
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

    def __delitem__(self, item):
        if item in self.value_keys:
            del self.values[item]
        if item in self.array_keys:
            del self.arrays[item]

    def copy(self):
        copy = FrameData()
        for key in self.value_keys:
            copy.raw.values[key].CopyFrom(self.raw.values[key])
        for key in self.array_keys:
            copy.raw.arrays[key].CopyFrom(self.raw.arrays[key])
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


class RecordView:
    """
    Base class that wraps the access to a kind of record.

    This class needs to be subclassed.
    """

    record_name: str | None = None  # MUST be overwritten as "arrays" or "values"
    singular: str | None = None  # MUST be overwritten as "array" or "value"

    def __init__(self, raw):
        if self.record_name is None or self.singular is None:
            raise NotImplementedError(
                "FieldView must be subclassed; record_name, singular, and"
                "_convert_to_python must be overwritten."
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
        raise NotImplementedError("Subclasses must overwrite the set method.")

    def delete(self, key):
        del self[key]

    @staticmethod
    def _convert_to_python(field):
        """
        Extract the value from a protobuf field so it is usable by python.

        The method needs to be adapted to the type of field that is manipulated.
        """
        raise NotImplementedError(
            "Subclasses must overwrite the _convert_to_python method."
        )

    def keys(self) -> Set:
        return set(self._raw_record.keys())


class ValuesView(RecordView):
    """
    Give access to singular values from a :class:`FrameData`.
    """

    record_name = "values"
    singular = "value"

    @staticmethod
    def _convert_to_python(field):
        return value_to_object(field)

    def set(self, key: str, value):
        self._raw_record[key].CopyFrom(object_to_value(value))


_EMPTY = object()


class ArraysView(RecordView):
    """
    Give access to homogeneous arrays from a :class:`FrameData`.
    """

    record_name = "arrays"
    singular = "array"

    @staticmethod
    def _convert_to_python(field):
        return field.ListFields()[0][1].values

    def set(self, key: str, value):
        try:
            reference_value = value[0]
        except IndexError:
            reference_value = _EMPTY
        except TypeError:
            raise ValueError(f"Value must be indexable for array {key}.") from None

        if reference_value is _EMPTY:
            type_attribute = "index_values"
        elif (
            isinstance(reference_value, numbers.Integral) and int(reference_value) >= 0
        ):
            type_attribute = "index_values"
        elif isinstance(reference_value, numbers.Real):
            type_attribute = "float_values"
        elif isinstance(reference_value, str):
            type_attribute = "string_values"
        else:
            raise ValueError(f"Cannot decide what type to use for [{key}]={value}")

        getattr(self._raw_record[key], type_attribute).values[:] = value
