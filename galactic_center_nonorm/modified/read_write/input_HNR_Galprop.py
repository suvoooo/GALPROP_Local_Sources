import math
import matplotlib.pyplot as plt

import numpy as np

print ("version of numpy: ", np.__version__)




def gaussdens(x, y, z):
    density =1. 
    R_sph = np.sqrt(x**2 + y**2 + z**2)
    if (-0.1<x<=0.1 and -0.1<y<=0.1 and -0.1<z<=0.1):
        return 1.
    else:    
        density = np.exp(-R_sph**2)
    return density


def gaussdens1(d1):
    density1 =1. 
    R_sph1 = np.sqrt(d1**2)
    if (-0.2<d<=0.2):
        return density1
    else:    
        density1 = np.exp(-R_sph1)
    return density1


def const(x, y, z):
    density = 1
    return density





######
filepath = '/media/abbl/d17db44d-1f0c-4659-891d-b7ff26b02041/suvo/Jupyter/'






####################################################
#+ Interpolate flux to match Galprop Energy bins
####################################################
def readexpnorm1():
    expnorm1file = open('./54_HypernovaGCNearbySource50Age1d1e4yrs_xyz2Rsun0d25Ecut200TeVEmax10PeVExpDecaySx0kpcEcutaddedTau1e5.d', 'r')
    expnorm1vals = {}
    while True:
        stringline = expnorm1file.readline()
        if stringline=='':
            break
        stringlist = stringline.split()
        E = float(stringlist.pop(0))/(10**3)
        flux = float(stringlist.pop(0)) * (10**7)
        expnorm1vals[E] = [flux]
    return expnorm1vals


def plotexpnorm1(fig, expnorm1vals=readexpnorm1()):
    Eexpnorm1 = []
    fluxexpnorm1 = []
    for point in sorted(expnorm1vals.keys()):
        Eexpnorm1.append(point)
        fluxexpnorm1.append(expnorm1vals[point][0] * (point**1))
    plt.plot(Eexpnorm1, fluxexpnorm1, linestyle='--', color='blue', label = r'$\tau = 1e4,$')








expnorm1vals = readexpnorm1()
expkeys = sorted(expnorm1vals.keys())
print ("len(expkeys): ", len(expkeys))
# print (hnrvals)    
# print (expkeys)
# for e in sorted(hnrvals.keys()):
    # print (10**e, e)






Egalp = []
starten = 1.0e4 # start galprop energy [MeV]
enden = 0.8e7 # end galprop energy [GeV]
bin_number = 81
en_factor = 1.4

for i in range(bin_number):
    En = (starten * (en_factor ** i))/(10**3)
    if En <= enden :
        Egalp.append(En)

print ("check the galprop energy bins in GeV: ")
print ("\n")
print (Egalp)
print ("number of bins: ", len(Egalp))



def hnrinterpol(Egal):
    intpolflux = 0. 
    for e in range(len(expkeys)):
        E2 = expkeys[e]
        # print ("logE2: ", logE2)
        if E2 > Egal:
            E1 = expkeys[e-1]
            # print ("logE1: ", logE1)
            npdp1 = expnorm1vals[E1][0]
            npdp2 = expnorm1vals[E2][0]
            intpolflux = npdp1 + ((npdp2 - npdp1) * (Egal - E1))/(E2 - E1)
            # print (npdp2, npdp1, e)
            return intpolflux
        # else: 
            # break

interpolnpdp = []
for i1 in range(len(Egalp)):
    intNp = hnrinterpol(Egalp[i1])
    # print (intNp, i1, Egalp[i1])
    interpolnpdp.append(intNp)
print ("length of interpolnpdp: ", len(interpolnpdp))

################################################


#+++++++++++++++++++++++++++++++++++++++++++++++++
#+ Create the input file for GALPROP
#+++++++++++++++++++++++++++++++++++++++++++++++++


d_file = open('./HNR_source.dat', 'w')
# 
xcord = []
ycord = []
zcord = []
Conv_HNRflux = []
for i in range(-40, 41, 1):
    xc = i/20
    xcord.append(xc)
    for j in range(-40, 41, 1):
        yc = j/20
        ycord.append(yc)
        for k in range(-40, 41, 1):
            zc = k/20
            zcord.append(zc)
            for l in range(len(Egalp)):
                intHNRflux = hnrinterpol(Egalp[l])
                intHNRfluxsp = ( intHNRflux * gaussdens(xc, yc, zc) )
                Conv_HNRflux.append(intHNRfluxsp)
                d_file.write('{0:2}' .format(intHNRfluxsp) + "\n")
# 
# 
# 
# print ("spatial bins: x, y, z: ",  len(xcord), len(ycord), len(zcord))

dist = [np.sqrt(i**2 + j**2 + k**2) for i, j, k in zip(xcord, xcord, zcord)]
print (len(dist))

# for d in range(len(dist)):
#     conv_flux = [fl*gaussdens1(dist[i]) for fl in interpolnpdp] 
# print ("cnv_flux: ",  len(conv_flux))

print ("len Conv_HNRflux: ", len(Conv_HNRflux))

# fig = plt.figure(figsize=(10, 7))
# plt.xlabel('Radius [Kpc]', fontsize=12)
# plt.ylabel('Npdp', fontsize=12)
# plt.plot(dist, , linestyle='-.', color='pink', alpha=0.8)
# plt.show()
