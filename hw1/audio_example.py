import scipy
import numpy as np
import scipy.io.wavfile
import pylab
import matplotlib
from haar import Haar

# Computes the Short-Time Fourier Transform (STFT) of a signal, with a given
# window length, and shift between adjacent windows
def stft(x, window_len=4096, window_shift=2048):
    w = scipy.hamming(window_len)
    X = scipy.array([scipy.fft(w*x[i:i+window_len])
        for i in range(0, len(x)-window_len, window_shift)])
    return scipy.absolute(X[:,0:window_len/2])

def short_time_haar(x, window_len=4096, window_shift=2048):
    X = scipy.array([Haar(x[i:i+window_len])
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
def plot_peaks(peak_list,start,end):
    fig = matplotlib.pyplot.figure()
    ax = fig.add_subplot(1,1,1)    
    s_list, f_list = zip(*peak_list)    
    matplotlib.pyplot.plot(s_list, f_list, 'bo')    
    ymin, ymax = ax.get_ylim()    
    ax.vlines((start,end),ymin,ymax,'red')
    matplotlib.pyplot.xlabel('Window index')
    matplotlib.pyplot.ylabel('Transform coefficient')

def graph_accel(signal):
    pylab.ion()
    pylab.figure()
    pylab.plot(signal)
    pylab.ioff()

# Finds peak in every 20x20 grid
def grid_peaks(X, G):
    peaks = []
    N = len(X) / G
    M = len(X[0]) / G
    XP = np.array(X)
    # array for peaks
    y = [0] * len(X)  
    for n in range(N):
        for m in range(M):
            # G x G grid 
            XX = XP[n*G:(n+1)*G, m*G:(m+1)*G]
            # index in grid
            pn, pm = np.unravel_index(XX.argmax(), XX.shape)  
            peak = XX[pn][pm]             # peak value
            pn += n*G                     # index in t
            pm += m*G
            peaks += [(pn, pm)]
            y[pn] = peak
    return y, peaks 

if __name__ == '__main__':
    # start and end point of 4 waves
    samples = [[135387,181251],[0,54684],[40572,86877],[185661,234612],
               [336483,388521],[302967,356328],[83349,140238],[239022,292824]]

    yy = []  # array for all y(t)

    # k+1 : 1,2,3,4
    for k in range(5): 
        rate, data = scipy.io.wavfile.read('../Data/tracks/%s' % str(k+1)+'.wav')
        if (len(data.shape) > 1):
            data = data[:,0]

        x = data[0:10*rate]

        X = stft(x)
        # X = short_time_haar(x)

        plot_transform(X)
        pylab.savefig('spectrogram'+str(k+1)+'.png')


        sample = samples[k]

        # start time of 2048 samples
        pstart = sample[0] / 2048   
        pend   = sample[1] / 2048

        # array of y(t)
        y, peaks = grid_peaks(X, 20)
        yy += [y]   
          
        plot_peaks(peaks,pstart,pend)
        pylab.savefig('peaks_'+str(k+1)+'.png')

    for n in range(0,5):         # n : 0,1,2,3
        for m in range(n+1, 5):  # m : n+1, n+2, ..
            yn=np.array(yy[n])
            ym=np.array(yy[m])
            corr = np.correlate(yn, ym, 'full')
            graph_accel(corr)
            pylab.savefig('corr_'+str(n+1)+str(m+1)+'.png')            
