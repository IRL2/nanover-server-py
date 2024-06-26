class Unpacker:
    """
    Unpack data representations from a buffer of bytes.

    The unpacking methods return the requested value from the next bytes of the
    buffer and move the cursor forward. They raise an `IndexError` if the
    buffer is too short to fulfill the request.
    """

    _buffer: bytes
    _cursor: int

    @classmethod
    def from_path(cls, path: str):
        with open(path, "rb") as infile:
            data = infile.read()
        return cls(data)

    def __init__(self, data: bytes):
        self._buffer = data
        self._cursor = 0

    def unpack_bytes(self, n_bytes: int) -> bytes:
        """
        Get the next `n_bytes` bytes from the buffer.

        The method raises a `ValueError` if the requested number of bytes is
        negative.
        """
        if n_bytes < 0:
            raise ValueError("Cannot unpack a negative number of bytes.")
        end = self._cursor + n_bytes
        if end >= len(self._buffer):
            raise IndexError("Not enough bytes left in the buffer.")
        bytes_to_return = self._buffer[self._cursor : end]
        self._cursor = end
        return bytes_to_return

    def unpack_u64(self) -> int:
        """
        Get an unsigned 64 bits integer from the next bytes of the buffer.
        """
        buffer = self.unpack_bytes(8)
        return int.from_bytes(buffer, "little", signed=False)

    def unpack_u128(self) -> int:
        """
        Get an unsigned 128 bits integer from the next bytes of the buffer.
        """
        buffer = self.unpack_bytes(16)
        return int.from_bytes(buffer, "little", signed=False)
