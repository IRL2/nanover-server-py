from dataclasses import dataclass
from typing import Callable


CommandHandler = Callable[..., dict | None]


@dataclass(kw_only=True)
class CommandRegistration:
    name: str
    label: str
    icon: str = "❓"
    arguments: dict
    handler: CommandHandler

    def run(self, arguments: dict):
        args = {}
        args.update(self.arguments)
        args.update(arguments)
        return self.handler(**args)

    def to_dict(self):
        return {
            "name": self.name,
            "label": self.label,
            "icon": self.icon,
            "arguments": self.arguments,
        }
