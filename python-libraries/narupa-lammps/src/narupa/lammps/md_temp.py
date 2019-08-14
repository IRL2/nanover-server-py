import MDAnalysis
#import quantities as pq
#q = 1 * pq.newton
#q.units = pq.m
#q9 = q.simplified
#pq.
#print(q)
hold = MDAnalysis.units.convert(1,"Newton","kJ/(mol*nm)")
print(hold)


hold = MDAnalysis.units.convert(1,"Newton","kJ/(mol*nm)")
print("Newton", hold)

hold = MDAnalysis.units.convert(0.0000001,"Newton","kJ/(mol*nm)")
print("picogram-micrometer/microsecond", hold)

hold = MDAnalysis.units.convert(0.001,"Newton","kJ/(mol*nm)")
print("atoogram-nanometer/nanosecond", hold)

hold = MDAnalysis.units.convert(1/100000,"Newton","kJ/(mol*nm)")
print("Dyne", hold)


hold = MDAnalysis.units.convert(1,"kcal/(mol*Angstrom)","kJ/(mol*nm)")
print("Kcal", hold)


hold = MDAnalysis.units.convert(1,"eV","kJ/mol")
print("eV", hold)
print("eV/nm", hold/10)

#hold = MDAnalysis.units.convert(1,"bohr","nm")
#print("bhor- > nm", hold)

hold = MDAnalysis.units.convert(1,"eV","kJ/mol")
print("eV", hold)
