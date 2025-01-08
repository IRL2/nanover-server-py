"""
Interface between NanoVer and OpenMM.
"""

from .converter import openmm_to_frame_data

GET_FRAME_INTERVAL_COMMAND_KEY = "trajectory/get-frame-interval"
SET_FRAME_INTERVAL_COMMAND_KEY = "trajectory/set-frame-interval"
GET_FORCE_INTERVAL_COMMAND_KEY = "imd/get-force-interval"
SET_FORCE_INTERVAL_COMMAND_KEY = "imd/set-force-interval"
GET_INCLUDE_VELOCITIES_COMMAND_KEY = "trajectory/get-include-velocities"
SET_INCLUDE_VELOCITIES_COMMAND_KEY = "trajectory/set-include-velocities"
GET_INCLUDE_FORCES_COMMAND_KEY = "trajectory/get-include-forces"
SET_INCLUDE_FORCES_COMMAND_KEY = "trajectory/set-include-forces"
