# Copyright (c) Intangible Realities Lab, University Of Bristol. All rights reserved.
# Licensed under the GPL. See License.txt in the project root for license information.
import MDAnalysis
from MDAnalysis.tests.datafiles import PSF, DCD  # test trajectory

#U contains both topology and  positions
u = MDAnalysis.Universe(PSF, DCD)  # always start with a Universe

import time

from narupa.trajectory import FrameServer
from narupa.mdanalysis import mdanalysis_to_frame_data


frameServer = FrameServer(address='localhost', port=54321)
###Send topolgy once at the begining of the server
# Get topolgy in the right grpc format
#takes u because it likes things in mdalaysis format
topology_data = mdanalysis_to_frame_data(u, topology=True, positions=False)

frame_index = 0
#Now actually send the frame, only contains the topology at this stage
frameServer.send_frame(frame_index, topology_data)

print("Starting Trajectory Server")

while True:

    for frame in u.trajectory:
        frame_data = mdanalysis_to_frame_data(u, topology=False, positions=True)

        frameServer.send_frame(frame_index, frame_data)
        time.sleep(1.0 / 30.0)
        frame_index = frame_index + 1


