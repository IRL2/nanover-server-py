import sys
from pathlib import Path

from nanover.openmm.serializer import deserialize_simulation
from nanover.openmm.bundles import bundle_openmm_simulation

XML_PATH = Path(" ".join(sys.argv[1:]))
BUNDLE_PATH = f"{XML_PATH.parent}/{XML_PATH.name.rstrip(".xml")}.openmm.zip"

print(XML_PATH, BUNDLE_PATH)

with open(XML_PATH) as xml_file:
    simulation = deserialize_simulation(xml_file.read())
    bundle_openmm_simulation(simulation, outfile=BUNDLE_PATH)