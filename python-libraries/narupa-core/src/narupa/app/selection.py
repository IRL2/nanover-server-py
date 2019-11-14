# Copyright (c) Intangible Realities Lab, University Of Bristol. All rights reserved.
# Licensed under the GPL. See License.txt in the project root for license information.
from typing import Dict, Iterable

INTERACTION_SINGLE = 'single'
INTERACTION_GROUP = 'group'
INTERACTION_RESTRAINT = 'restraint'

RENDERER_LIQUORICE = 'liquorice'


class NarupaImdSelection:
    @classmethod
    def from_selection_message(cls, message):
        selection = cls(message['id'], message['name'])
        selection.select_particles(
            get_nested_or_default(
                message,
                set(),
                'selected',
                'particle_ids',
            )
        )
        selection.interaction_method = get_nested_or_default(
            message,
            INTERACTION_SINGLE,
            'properties',
            'narupa.interaction',
            'method',
        )
        selection.velocity_reset = get_nested_or_default(
            False,
            'properties',
            'narupa.interaction',
            'velocity.reset',
        )
        selection.rendering_renderer = get_nested_or_default(
            RENDERER_LIQUORICE,
            'properties',
            'narupa.rendering',
            'renderer',
        )

        return selection

    def __init__(self, id: str, name: str = 'Unnamed Selection'):
        self.selection_id = id
        self.selection_name = name
        self.selected_particle_ids = set()

        self.interaction_method = INTERACTION_SINGLE
        self.velocity_reset = False
        self.rendering_renderer = RENDERER_LIQUORICE

    def clear_particles(self):
        self.selected_particle_ids.clear()

    def select_particles(self, particle_ids: Iterable[int]):
        self.selected_particle_ids.update(particle_ids)

    def to_selection_message(self):
        return {
            'id': self.selection_id,
            'name': self.selection_name,
            'selected': {
                'particle_ids': list(self.selected_particle_ids),
            },
            'properties': {
                'narupa.interaction': {
                    'method': self.interaction_method,
                    'velocity.reset': self.velocity_reset,
                },
                'narupa.rendering': {
                    'renderer': self.rendering_renderer,
                }
            }
        }


def get_nested_or_default(dict: Dict, default, *layers: Iterable[str]):
    for layer in layers:
        try:
            dict = dict[layer]
        except (ValueError, TypeError):  # GRPC raises these
            return default
    return dict
