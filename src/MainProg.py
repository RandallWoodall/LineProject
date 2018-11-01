import json
import numpy
import scipy.stats.mstats as st
import scipy.constants as constants
from Conductor import conductor

# Load in the pre-built conductor table and make a dictionary of Conductor objects.
conductors = {}
conductorsRaw = json.loads(open("conductors.json", "r").read())
for conductorKey in conductorsRaw.keys():
    conductors[conductorKey] = conductor(conductorsRaw[conductorKey], conductorKey)

# Get inputs
nomVolt = float(input("Input nominal voltage (kV): "))
nomVolt = nomVolt * 10**3
lineLen = float(input("Input line length (km): "))
bundleNum = int(input("Input conductors per bundle: "))
bundleDist = float(input("Input distance between conductors in a bundle (feet): "))
phaseDist = float(input("Input distance between phases (feet): "))
lineName = input("Input the conductor name or circular mills: ")
# conductorChoice gives us gmr(feet), diam & radius (inches), ac resistance per mile,
# current carrying capacity
try:
    conductorChoice = conductors[lineName]
except KeyError:
    for key in conductors.keys():
        if conductors[key].circ_mils == int(lineName):
            conductorChoice = conductors[key]
            lineName = conductorChoice.name
            break;
if conductorChoice == None:
    print("Error choosing conductor.")
seriesComp = float(input("Input % series compensation: ")) / 100
shuntComp = float(input("Input % shunt compensation: ")) / 100
load = float(input("Input load power (MW): "))
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
reactance_per_km = numpy.log(Deq/DS) * 754e-4
# Calculate suseptance/km
suseptance_per_km = 240000 * constants.pi**2 * constants.epsilon_0 / numpy.log(Deq/DSC)

# Get ABCD
z = resistance_per_km + 1j*reactance_per_km
y = 1j*suseptance_per_km

gamma = numpy.sqrt(z * y)
Zc = numpy.sqrt(z/y)
Beta = numpy.imag(gamma)
A = numpy.cosh(gamma * lineLen)
D = A
B = Zc * numpy.sinh(gamma*lineLen)
C = numpy.sinh(gamma*lineLen)/Zc
ABCDline = numpy.matrix([[A,B],[C,D]])

# Use compensation to produce capacitor and inductor matrices
compCap = seriesComp / 2 * numpy.imag(z * lineLen)
ABCDcap = numpy.matrix([[1, -1j*compCap], [0, 1]])
compInd = shuntComp / 2 * numpy.imag(y * lineLen)
ABCDind = numpy.matrix([[1, 0], [1j * compInd, 1]])

# Apply compensation
step1 = numpy.matmul(ABCDind, ABCDcap)
step2 = numpy.matmul(step1, ABCDline)
ABCDnoLoad = numpy.matmul(step2, ABCDcap)
ABCDfullLoad = numpy.matmul(ABCDnoLoad, ABCDind)
Afl = ABCDfullLoad[0,0]
Bfl = ABCDfullLoad[0,1]
Cfl = ABCDfullLoad[1,0]
Dfl = ABCDfullLoad[1,1]
Anl = ABCDnoLoad[0,0]
Bnl = ABCDnoLoad[0,1]
Cnl = ABCDnoLoad[1,0]
Dnl = ABCDnoLoad[1,1]

# Get complex power per phase
if leadLag == 'lag':
    S = (load * 10**6) * (powerFactor + 1j * numpy.sin(numpy.arccos(powerFactor)))
else:
    S = (load * 10**6) * (powerFactor - 1j * numpy.sin(numpy.arccos(powerFactor)))
Sphase = S/3
    
# Get line to nuetral voltage
Vrfl = nomVolt / numpy.sqrt(3)

# Calculate line current
Iline = numpy.conj(Sphase / Vrfl)

# Calculate source voltage, current, power, and losses
Vs = Afl * Vrfl + Bfl * Iline
Is = Cfl * Vrfl + Dfl * Iline
Ssfl = 3 * Vs * numpy.conj(Is)
lossfl = numpy.real(Ssfl) - numpy.real(S)

# Get no load voltage + current
Vrnl = Vs / Anl
Isnl = Cnl * Vrnl
Ssnl = 3 * Vs * numpy.conj(Isnl)
lossnl = numpy.real(Ssnl)

# Calculate voltage regulation
Regulation = (numpy.abs(Vrnl) - numpy.abs(Vrfl))/numpy.abs(Vrfl) * 100

# Calculate loadability @ 35 degrees
vMag = numpy.abs(Vrfl)
Pmax = (vMag**2 * numpy.cos(numpy.angle(Bnl) - 35 * constants.pi/180)/numpy.abs(Bnl) 
        - vMag**2 * numpy.real(Anl/Bnl)) * 3

# Output all the desired values
print("Receiving end voltage (kV): " + str(nomVolt*10**-3))
print("Line Length (km): " + str(lineLen))
print("Conductors per bundle: " + str(bundleNum))
print("Distance between lines in bundle: " + str(bundleDist))
print("Conductor name: " + lineName)
print("Conductor gmr (ft): " + str(conductorChoice.gmr))
print("Conductor radius (ft): " + str(conductorChoice.out_diam / 24))
print("AC resistance @ 50C: " + str(conductorChoice.r_60))
print("Ampacity: " + str(conductorChoice.amp_cap))
print("Series compensation: " + str(seriesComp * 100) + "%")
print("Shunt compensation: " + str(shuntComp * 100) + "%")
print("Load (MW): " + str(load))
print("Power factor: " + str(powerFactor) + " " + leadLag)
print("Voltage regulation: " + str(Regulation) + "%")
print("Loadability at 35 degrees (MW): " + str(Pmax * 10**-6))

if Regulation < 10 and Pmax > 1.5 * load and numpy.abs(load)/numpy.imag(Ssnl) \
    and numpy.abs(load)/numpy.imag(Ssfl):
        print("Performance is acceptable!")
else:
    print("Performance is unacceptable!")
    
print("Full load ABCD: ")
print(ABCDfullLoad)
print("Sending end voltage (kV) (fl): " + str(Vs * 10**-3))
print("Sending end current (A) (fl): " + str(Is))
print("Sending end power (MVA) (fl): " + str(Ssfl * 10**-6))
print("Receiving end voltage (kV) (fl): " + str(Vrfl * 10**-6))
print("Receiving end current (A) (fl): " + str(Iline))
print("Receiving end power (MVA) (fl): " + str(S * 10**-6))
print("Full load losses (MW) (fl): " + str(lossfl * 10**-6))

print("No load ABCD: ")
print(ABCDnoLoad)
print("Sending end voltage (kV) (nl): " + str(Vs * 10**-6))
print("Sending end current (A) (nl): " + str(Isnl))
print("Sending end power (MVA) (nl): " + str(Ssnl * 10**-6))
print("Receiving end voltage (kV) (nl): " + str(Vrnl * 10**-3))
print("Receiving end current (A) (nl): 0 (disconnected)")
print("Receiving end power (VA) (nl): 0 (disconnected)")
print("Full load losses (MW) (nl): " + str(lossnl * 10**-6))
