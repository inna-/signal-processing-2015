import scipy
import numpy as np
import scipy.io.wavfile
import pylab
import matplotlib
import operator
from dtw import *
from erdp import *

import matplotlib.pyplot as plt

# Computes the Short-Time Fourier Transform (STFT) of a signal, with a given
# window length, and shift between adjacent windows
def stft(x, window_len=4096, window_shift=2048):
    w = scipy.hamming(window_len)
    X = scipy.array([scipy.fft(w*x[i:i+window_len])
        for i in range(0, len(x)-window_len, window_shift)])
    return scipy.absolute(X[:,0:window_len/2])

# Plot a transformed signal X, i.e., if X = stft(x), then
# plot_transform(X) plots the spectrogram of x
def plot_transform(X):
    pylab.ion()
    pylab.figure()
    pylab.imshow(scipy.log(X.T), origin='lower', aspect='auto', interpolation='nearest', norm=matplotlib.colors.Normalize())
    pylab.xlabel('Window index')
    pylab.ylabel('Transform coefficient')
    pylab.ioff()

# Plot a list of peaks in the form [(s1, f1), (s2, f2), ...]
def plot_peaks(peak_list,xsample,ysample):
    fig = matplotlib.pyplot.figure()
    ax = fig.add_subplot(1,1,1)    
    s_list, f_list = zip(*peak_list)    
    matplotlib.pyplot.plot(s_list, f_list, 'bo')    
    if ysample == (0,0):
        ymin, ymax = ax.get_ylim()    
        xmin, xmax = ax.get_xlim()
    else:
        ymin, ymax = ysample
        xmin, xmax = xsample
    ax.vlines((xsample),ymin,ymax,'red')
    ax.hlines((ysample),xmin,xmax,'red')
    matplotlib.pyplot.xlabel('Window index')
    matplotlib.pyplot.ylabel('Transform coefficient')

def graph_accel(signal):
    pylab.ion()
    pylab.figure()
    pylab.plot(signal)
    pylab.ioff()

def grid_peaks(X, K):
    XP = np.array(X)
    N = len(X) 
    M = len(X[0])
    peaks = []
    for n in range(K,N-K):
        for m in range(K,M-K):
            XX = XP[n-K:n+K, m-K:m+K]
            pn, pm = np.unravel_index(XX.argmax(), XX.shape)
            if XX[pn][pm] == XP[n][m]:
                peaks.append((n,m))
    return peaks


def hash_pmaps(peaks):
    F = 200
    T = 10
    N = 10
    hashtable = {}
    for p in peaks:
        n, m = p   # anchor
        for q in peaks:
            qn, qm = q
            if qn <= n     or qn >= n + T: continue      # target zone
            if qm <= m - F or qm >= m + F: continue
            hm = (m+2)/4                                 # allow collisions
            hpn = (qn-n+2)/4                      
            hpm = (qm-m+F+2)/4
            h = hm + hpn * 1000 + hpm * 10000000         # approximation hash

            #h = m + (qn-n) * 1000 + (qm-m+F) * 10000000  
            #h = m * 1000000 + (qn-n) * 1000 + (qm-m+F)   
            h = (qn-n) * 1000000 + m * 1000 + (qm-m+F)   # hash for exact match
            if n in hashtable:
                hashtable[n] += [h]
            else:
                hashtable[n]  = [h]
    return hashtable

def conv_to_1D(y, sample):   # B.1
    z = []
    ysort = sorted(y.items(), key=operator.itemgetter(0))
    for ys in ysort:
        if ys[0] >= sample[0] and ys[0] <= sample[1]:   # sub sequence
            y1 = ys[1]
            y1.sort()
            z += y1
    return z    
                
if __name__ == '__main__':
    # start and end point of waves
    N = 200
    samples = [(x/2048, y/2048) for (x,y) in
               [[135387,181251],[0,54684],[40572,86877],[185661,234612],
               [336483,388521],[302967,356328],[83349,140238],[239022,292824]]
              ]

    yy = {}
    tracks = [0,3,5]                              # waves 1,4,6 (audio 1,2,3) 

    for k in tracks: 
        rate, data = scipy.io.wavfile.read(str(k+1)+'.wav')
        if (len(data.shape) > 1):
            data = data[:,0]
        xsize = min(10*rate, 407500)
        x = data[0:xsize]

        X = stft(x)
        peaks = grid_peaks(X, 10)
        yy[k] = hash_pmaps(peaks)

    print 'DTW & EDRP by entire sequences'         # B.4 (dtw by entire sequence)
    for n in tracks: 
        for m in tracks:
            if n == m: continue
            z1 = conv_to_1D(yy[n], (0,N)) # get 1-D of entire sequence
            z2 = conv_to_1D(yy[m], (0,N))   
            dist1, cost, path = dtw(z1, z2)  
            # dist2 = edrp1(z1, z2, len(z1), len(z2))  #causes maximum recursion depth :(              
            dist3 = edrp2(z1, z2)                
            print n+1, m+1, dist1, dist3

    print '\nDTW by window approach'        # B.5 (dtw by sliding window)
    n = 0                                   # wave 1 (audio 1)
    m = 3                                   # wave 4 (audio 2)
    distance = []
    z1 = conv_to_1D(yy[n], samples[n])      # mutual sentence
    L = samples[n][1] - samples[n][0]       # start & end point
    for l in range(N-L):                    # window approach
        z2 = conv_to_1D(yy[m], (l, l+L))
        dist, cost, path = dtw(z1, z2)  
        distance += [dist]
    graph_accel(distance)
    pylab.savefig('dtw_'+str(n+1)+'_'+str(m+1)+'.png')            
    idx = distance.index(min(distance))
    print n+1, m+1, idx, distance[idx]      # minimum distance
    
    print '\nDTW by window approach'        # B.6 (dtw by sliding window)
    n = 3                                   # wave 4 (audio 2)
    m = 0                                   # wave 1 (audio 1)
    distance = []
    z1 = conv_to_1D(yy[n], samples[n])      # mutual sentence
    L = samples[n][1] - samples[n][0]       # start & end point
    for l in range(N-L):                    # window approach
        z2 = conv_to_1D(yy[m], (l, l+L))
        dist, cost, path = dtw(z1, z2)  
        distance += [dist]
    graph_accel(distance)
    pylab.savefig('dtw_'+str(n+1)+'_'+str(m+1)+'.png')            
    idx = distance.index(min(distance))
    print n+1, m+1, idx, distance[idx]      # minimum distance

