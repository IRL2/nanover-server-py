"""
Module providing a server for running a :class: ImdService.
"""

from typing import Optional

from nanover.core import (
    NanoverServer,
    get_requested_port_or_default,
    DEFAULT_SERVE_ADDRESS,
)
from nanover.imd.imd_state import ImdStateWrapper

DEFAULT_PORT = 54322


class ImdServer(NanoverServer):
    """
    Class providing a NanoverServer with an ImdStateWrapper for accessing
    IMD-specific state.

    :param: address: URL or IP address at which to run the server.
    :param: port: Port at which to run the server.
    """

    def __init__(
        self,
        *,
        address: Optional[str] = None,
        port: Optional[int] = None,
    ):
        if address is None:
            address = DEFAULT_SERVE_ADDRESS
        port = get_requested_port_or_default(port, DEFAULT_PORT)

        super().__init__(address=address, port=port)
        self._imd_state = ImdStateWrapper(self._state_service.state_dictionary)

    @property
    def imd_state(self) -> ImdStateWrapper:
        """
        An ImdStateWrapper for accessing the interaction-relevant state of this
        server.
        """
        return self._imd_state
