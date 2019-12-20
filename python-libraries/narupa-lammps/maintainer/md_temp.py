# This is a valueing routine for various conversion factiors in mdanalysis,
# this is just used as a useful sanity checking tool for varous conversions.

import MDAnalysis
#print(q)
value = MDAnalysis.units.convert(1,"Newton","kJ/(mol*nm)")
print(value)


value = MDAnalysis.units.convert(1,"Newton","kJ/(mol*nm)")
print("Newton", value)

value = MDAnalysis.units.convert(0.0000001,"Newton","kJ/(mol*nm)")
print("picogram-micrometer/microsecond", value)

value = MDAnalysis.units.convert(0.001,"Newton","kJ/(mol*nm)")
print("atoogram-nanometer/nanosecond", value)

value = MDAnalysis.units.convert(1/100000,"Newton","kJ/(mol*nm)")
print("Dyne", value)

value = MDAnalysis.units.convert(1,"kcal/(mol*Angstrom)","kJ/(mol*nm)")
print("Kcal", value)

value = MDAnalysis.units.convert(1,"eV","kJ/mol")
print("eV", value)
print("eV/nm", value/10)

value = MDAnalysis.units.convert(1,"bohr","nm")
print("bhor- > nm", value)

value = MDAnalysis.units.convert(1,"eV","kJ/mol")
print("eV", value)
