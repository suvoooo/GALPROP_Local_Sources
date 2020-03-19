import math
import numpy as np
import matplotlib.pyplot as plt


filepath = '/media/abbl/d17db44d-1f0c-4659-891d-b7ff26b02041/suvo/Jupyter/'

def hnr1fileread():
    hnrfile = open(filepath+'HNR_Np_esc_n1_check.txt', 'r')
    hnrvals = {}
    epochline = hnrfile.readline()
    timeline = hnrfile.readline()
    mombinline = hnrfile.readline()
    pline = hnrfile.readline()
    npdpline = hnrfile.readline()
    Escenline = hnrfile.readline()
    note1line = hnrfile.readline()
    note2line = hnrfile.readline()
    print ("noteline2 :", note2line)
    while True:
        stringline = hnrfile.readline()
        # print ("stringline :", stringline)
        if stringline=='':
            break
        stringlist = stringline.split()
        epoch = float(stringlist.pop(0))
        time = float(stringlist.pop(0))
        pbin = float(stringlist.pop(0))
        momentum = float(stringlist.pop(0)) # given as log 10 (p) [GeV]
        npdp = float(stringlist.pop(0)) # given as log(n)
        Eesc = float(stringlist.pop(0)) # given in erg
        hnrvals[momentum] = [epoch, npdp, time, Eesc]
        #if epoch =2:
            # break
    return hnrvals


hnrvals = hnr1fileread()
hnrkeys = sorted(hnrvals.keys())
print (len(hnrkeys))
# print (hnrvals)    

for e in sorted(hnrvals.keys()):
    print (10**e, e)






def plotcheckhnr1(fig, hnrvals = hnr1fileread()):
    Epoch = []
    Time = []
    NpDp = []
    Momentum = []
    for point in sorted(hnrvals.keys()):
    	logMomentum = point
        # print (logMomentum)
    	Momentum.append(10**(logMomentum))
    	lognpdp = hnrvals[point][1]
    	NpDp.append(10**(lognpdp) * (10**logMomentum))
        # print (lognpdp, logMomentum, hnrvals[point][0])
    plt.plot(Momentum, NpDp, marker='*', markersize=3, linestyle=':', linewidth=2, color='red', label='Original') 








Egalp = []
starten = 1.7e7 # start galprop energy [MeV]
enden = 4.9e6 # end galprop energy [GeV]
bin_number = 81
en_factor = 1.1

for i in range(bin_number):
    En = (starten * (en_factor ** i))/(10**3)
    if En < enden :
        Egalp.append(En)

print ("check the galprop energy bins in GeV: ")
print ("\n")
print (Egalp)
print ("number of bins: ", len(Egalp))



def hnrinterpol(Egal):
    intpolflux = 0. 
    for e in range(len(hnrkeys)):
        logE2 = hnrkeys[e]
        # print ("logE2: ", logE2)
        E2 = 10**(logE2)
        if E2 > Egal:
            logE1 = hnrkeys[e-1]
            # print ("logE1: ", logE1)
            E1 = 10**(logE1)
            lognpdp1 = hnrvals[logE1][1]
            npdp1 = 10**(lognpdp1)
            lognpdp2 = hnrvals[logE2][1]
            npdp2 = 10**(lognpdp2)
            intpolflux = npdp1 + ((npdp2 - npdp1) * (Egal - E1))/(E2 - E1)
            return intpolflux
        # else: 
            # break    


interpolnpdp = []
for i in range(len(Egalp)):
    intNp = hnrinterpol(Egalp[i])
    interpolnpdp.append(intNp*(Egalp[i]))


fig = plt.figure(figsize=(10, 7))
plt.xlabel('Energy [GeV]', fontsize=12)
plt.ylabel('Npdp', fontsize=12)
plt.xscale('log')
plt.yscale('log')
plt.plot(Egalp, interpolnpdp, linestyle='--', color='navy', label='Interpolated Flux')
plotcheckhnr1(fig)
plt.legend(fontsize=12)
plt.show()


#finished