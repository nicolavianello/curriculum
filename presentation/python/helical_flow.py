# programma per il calcolo del profilo in u della velocita in diversi impulsi in QSH con diversi valori di n/ng

import sys
import os
outputDir = "../"
if sys.platform == 'darwin':
    sys.path.append("/Users/nicola/.pythonlib")
    inputDir = "/Users/nicola/LN/linuxHome/idl/rfx/isis/topologyflow/dati/"
else:
    sys.path.append("/home/vianello/.pythonlib")
    inputDir = "/home/vianello/idl/rfx/isis/topologyflow/dati/"

# first personal library
import bin_by
from bin_by import *
import smooth
from smooth import * 
import congrid
from congrid import * 
# then generic library of matplotlib
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import rc
from matplotlib import *
from matplotlib.patches import Polygon
from matplotlib import pylab
rc('font',**{'family':'sans-serif','sans-serif':['Tahoma']})
rc("font", size=18)
## params = {'legend.fontsize': 10}
## pylab.rcParams.update(params)

# numpy libraries
import numpy as np
# scipy libraries
from scipy import stats
from scipy import ndimage
from scipy import io
# others
outputDir = "../"

# restore number of shots and corresponding values of n/ng
nngData = np.loadtxt("/home/vianello/idl/rfx/isis/topologyflow/dati/valori_nng.txt")
nng = nngData[:,1]
gShot = nngData[:, 0]

# now choosing the three range of data

shotGlobal = [gShot[(nng <= 0.1)], gShot[((nng > 0.1) &  (nng < 0.2))], gShot[(nng >= 0.2)]]
nngAvg = [np.average(nng[(nng <= 0.1)]), np.average(nng[((nng > 0.1) &  (nng < 0.2))]),np.average(nng[(nng >= 0.2)])]

fig = plt.figure(figsize=(6.38,5.38))
fig.subplots_adjust(bottom=0.12,left=0.16)
ax = fig.add_subplot(111)
colori = ['b', 'r', 'm', 'g']
colori2 = ['cyan', 'orange', 'orchid', 'lime']
costante = np.array([3, 2., 2., 1])

i=0
shotlist = np.array(shotGlobal[i])
for j in range(len(shotlist)):
    if os.path.isfile(inputDir+'angles/uAngleGianThetaShot' + shotlist[j].astype('|S5')+'.sav'):
        tmp = io.readsav(inputDir + 'flow/toroidal/'+shotlist[j].astype('|S5')+'.sav')
        v = tmp.out.velocity[0]
        a = tmp.out.angle[0]
        ind = np.argsort(a)
        t= tmp.out.time[0]
        v= v[np.argsort(a),:]
        a= np.sort(a)
        # now gaussian laplacian filter
        vfilt = ndimage.gaussian_laplace(v, [1.76,2.82], mode = 'constant')
        ## vfilt = ndimage.median_filter(v, size=3, mode='constant')
        # now rebin on the same number of point of angle
        tmp= io.readsav(inputDir+'angles/uAngleGianThetaShot' + shotlist[j].astype('|S5')+'.sav')
        uangle = tmp.uangle[5,: ,:]
        utime = tmp.utime
        # now verify that velocity time and uangle time have the same minimum and maximum
        tmn = np.choose(np.less(np.min(t), np.min(utime)), (np.min(t), np.min(utime)))
        tmx = np.choose(np.greater(np.max(t), np.max(utime)), (np.max(t), np.max(utime)))
        i1 = ((t >= tmn) & (t <= tmx))
        i2 = ((utime >= tmn) &  (utime <= tmx))
        vcon = np.reshape(congrid(vfilt[:,i1],
                                  [len(a), len(utime[i2])],
                                  method = 'spline'),- 1) # transformed in 1D array
        # now concatenate
        if j == 0:
            vout = vcon
            uout = np.reshape(uangle[:, i2],- 1)
        else:
            vout = np.append(vout, vcon)
            uout = np.append(uout, np.reshape(uangle[:, i2],- 1))

            
            # now transform in [0,2pi] range
res = np.fmod(uout, 2.* np.pi)
uout[(res <= 0)] += 2 * np.pi
vbin, bins, bm, bw = bin_by(uout, vout, nbins = 8)
vavg = np.array([0.0] * len(bm))
vstd = np.array([0.0] * len(bm))
for k in range(len(bm)):
    vavg[k] = np.average(np.asarray(vbin[k]))
    vstd[k] = stats.sem(np.asarray(vbin[k])) 
    
plt.plot(bm, costante[i] * vavg, colori[i] + '*--',  markeredgecolor='k',
         markersize=25,alpha=0.7)
ax.fill_between(bm, costante[i] * (vavg - vstd), costante[i] * (vavg + vstd),
                facecolor = colori2[i], alpha = 0.5,
                interpolate = True)
        
# ax.set_xlim([ - np.pi, np.pi])
ax.set_xlim([ 0, 2 * np.pi])
ax.set_ylim([ - 300, 300])
ax.set_xlabel(r'u', fontsize = 18)
ax.set_ylabel(r'$\delta$v$_{\perp}$ [m/s]', fontsize = 18)
lab = [r'0', r'$\pi/2 $', '$\pi$', r'$3\pi/2$', r'$2\pi$']
plt.xticks([ 0, np.pi / 2., np.pi, 3 * np.pi / 2., 2 * np.pi], lab)
plt.savefig(outputDir+"pdf_box/helical-flow.pdf",transparent=True)



