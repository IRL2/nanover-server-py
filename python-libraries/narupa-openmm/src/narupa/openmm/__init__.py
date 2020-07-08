# Copyright (c) Intangible Realities Lab, University Of Bristol. All rights reserved.
# Licensed under the GPL. See License.txt in the project root for license information.
"""
Interface between Narupa and OpenMM.
"""
from .converter import openmm_to_frame_data
from .narupareporter import NarupaReporter
from .runner import (
    OpenMMRunner,
    GET_FRAME_INTERVAL_COMMAND_KEY, SET_FRAME_INTERVAL_COMMAND_KEY,
    GET_FORCE_INTERVAL_COMMAND_KEY, SET_FORCE_INTERVAL_COMMAND_KEY,
)
