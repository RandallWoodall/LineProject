import json
import numpy
import scipy.stats.mstats as st
import scipy.constants as constants
from Conductor import conductor

# Load in the pre-built conductor table and make a dictionary of Conductor objects.
conductors = {}
conductorsRaw = json.loads(open("conductors.json", "r").read())
for conductorKey in conductorsRaw.keys():
    conductors[conductorKey] = conductor(conductorsRaw[conductorKey])

# Get inputs
nomVolt = float(input("Input nominal voltage (kV): "))
lineLen = float(input("Input line length (km): "))
bundleNum = int(input("Input conductors per bundle: "))
bundleDist = float(input("Input distance between conductors in a bundle (feet): "))
phaseDist = float(input("Input distance between phases (feet): "))
lineName = input("Input the conductor name: ")
# conductorChoice gives us gmr(feet), diam & radius (inches), ac resistance per mile,
# current carrying capacity
conductorChoice = conductors[lineName]
seriesComp = float(input("Input % series compensation: ")) / 100
shuntComp = float(input("Input % shunt compensation: ")) / 100
load = float(input("Input load power (MVA): "))
powerFactor = float(input("Input power factor: "))
leadLag = input("Input lead/lag: ")

# Calculate resistance/km
resistance_per_km = conductorChoice.r_60 * 1.609 / bundleNum
# Calculate Deq and DS, and DSC
# !!ASUMPTION flat (H-frame)
Deq = st.gmean([phaseDist, phaseDist, phaseDist*2])
DS = st.gmean([conductorChoice.gmr] + [bundleDist] * (bundleNum - 1))
DSC = st.gmean([conductorChoice.out_diam / 24] + [bundleDist] * (bundleNum - 1))
# Calculate reactance/km
reactance_per_km = constants.pi * 240000 * 10**(-7) * numpy.log(Deq/DS)
# Calculate suseptance/km
suseptance_per_km = 2000 * constants.pi * constants.epsilon_0 / numpy.log(Deq/DSC)

print("Resistance/km: " + str(resistance_per_km))
print("Reactance/km: " + str(reactance_per_km))
print("Suseptance/km: " + str(suseptance_per_km))
