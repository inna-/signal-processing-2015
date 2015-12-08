import pylab

def trim(x):
    N = len(x)
    M = 1
    while M < N:
        M *= 2
    x = x[:M/2]
    return x

def trunc(x, k):
    N = len(x)
    M = N / (2**k)
    # x = x[:M] + [0.0]*(N-M)
    return x[:M]
    
def Haar(x):
    N = len(x)
    har = [0.0]*N

    M = N / 2
    while True:
        for i in range(M):
            har[i]   = (x[i*2] + x[i*2+1]) / 2
            har[M+i] = (x[i*2] - x[i*2+1]) / 2

        if M == 1:
            return har

        x = har[:M*2]
        M = M / 2

def invHaar(x):
    N = len(x)
    har = [0.0]*N

    M = 1
    while M < N:
        for i in range(M):
            har[i*2]   = x[i] + x[M+i]
            har[i*2+1] = x[i] - x[M+i]

        M = M * 2
        x[:M] = har[:M]
        
    return har

def graph_accel(signal):
    pylab.ion()
    pylab.figure()
    pylab.plot(signal)
    pylab.ioff()

if __name__ == '__main__':

    # A. Haar transformed data
    # Get just the displacement in the x coordinate
    f = open('accel.csv', 'r')
    x = [float(line.split(',')[6]) for line in f]
    x = trim(x)  # trim to length M 
    M = len(x)

    graph_accel(x)
    pylab.savefig("2A_HaarTransform.png")
    
    # B. Haar transforms and “edges”
    hx = Haar(x)
    hf2 = hx[(M/2):]  # second half (X1)
    
    hpp = [(k+abs(k))/2 for k in hf2]  # positive part
    pent = [hpp.index(k) for k in sorted(hpp, reverse=True)]
    
    print 'Data set'
    print 'positive entries'
    for n in range(5):
        # pent[n] is the index for positive entry 
        # hpp[pent[n]] is the Haar at every index
        print n, x[pent[n]*2], x[pent[n]*2+1], hpp[pent[n]]

    hnn = [(k-abs(k))/2 for k in hf2]  # negative part
    nent = [hnn.index(k) for k in sorted(hnn, reverse=False)]
    print 'negative entries'
    for n in range(5):
        print n, x[nent[n]*2], x[nent[n]*2+1], hnn[nent[n]]

    # C. Smoothing
    for g in range(5):
        hh = trunc(hx, g)   # truncate, k=0 means no truncation
        x2 = invHaar(hh)
        graph_accel(x2)
        pylab.savefig('accel_'+'_'+str(g)+'.png')
        #pylab.savefig("2C_Smoothing.png")
        pylab.show()
