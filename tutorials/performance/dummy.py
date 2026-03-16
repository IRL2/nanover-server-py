from uuid import uuid4
from random import random

from nanover.utilities.change_buffers import DictionaryChange
from nanover.websocket import NanoverImdClient
from nanover.utilities.cli import suppress_keyboard_interrupt_as_cancellation
from nanover.utilities.timing import yield_interval

userid = str(uuid4())
headset = {
    "name": "headset",
    "position": [0.0, 0.0, 0.0],
    "rotation": [0.0, 0.0, 0.0, 1.0],
}
update = {
    f"avatar.{userid}": {
        "components": [headset],
        "playerid": userid,
        "name": "DUMMY",
    }
}

with suppress_keyboard_interrupt_as_cancellation() as cancellation:
    with NanoverImdClient.from_discovery(server_name="") as client:
        print(client.wait_until_minimum_usable_frame(timeout=30))

        for _ in yield_interval(1 / 30):
            if cancellation.is_cancelled:
                break

            headset["position"] = [random() * 0.5 + 1 for _ in range(3)]
            client.update_state(DictionaryChange(updates=update))
