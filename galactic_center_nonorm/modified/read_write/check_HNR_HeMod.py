import math 
import matplotlib.pyplot as plt

filepath = '/media/abbl/d17db44d-1f0c-4659-891d-b7ff26b02041/suvo/Jupyter/'

def hnr1fileread():
    hnrfile = open(filepath+'HNR_Np_esc_n1_copy.txt', 'r')
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
        logmomentum = float(stringlist.pop(0)) # given as log 10 (p) [GeV]
        lognpdp = float(stringlist.pop(0)) # given as log(n)
        Eesc = float(stringlist.pop(0)) # given in erg
        hnrvals[epoch] = [time, pbin, logmomentum, lognpdp, Eesc]
        #if epoch =2:
            # break
    return hnrvals


hnrvals = hnr1fileread()
epochs = sorted(hnrvals.keys())
print ("length of epochs :", len(epochs))
print ("\n")
print ("epochs :", epochs) 

for i in range(len(sorted(hnrvals.keys()))):
    while epochs[i]==i+1:
        logpdp = hnrvals[epochs[i]][2]
        lognp = hnrvals[epochs[i]][3]
        


# def plotcheckhnr1(fig, hnrvals = hnr1fileread()):
#     Epoch = []
#     Time = []
#     NpDp = []
#     Momentum = []
#     for point in sorted(hnrvals.keys()):
#     	logMomentum = point
#         # print (logMomentum)
#     	Momentum.append(10**(logMomentum))
#     	lognpdp = hnrvals[point][1]
#     	NpDp.append(10**(lognpdp) * (10**logMomentum))
#         # print (lognpdp, logMomentum, hnrvals[point][0])
#     plt.plot(Momentum, NpDp, marker='*', markersize=3, linestyle='-', linewidth=2, color='red')


# fig = plt.figure(figsize=(12, 9))
# plt.xlabel('p [GeV]')
# plt.ylabel('Npdp')
# plt.yscale('log')
# plt.xscale('log')
# plotcheckhnr1(fig)
# plt.show()
