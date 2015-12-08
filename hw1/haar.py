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

dat = [9., 7., 3., 5., 11., 7., 13., 17.]
dat = [100., 200., 44., 50., 20., 20., 4., 2.]
dat = [9., 7., 3., 5.]
dat = [9, 7, 3, 5, 6, 10, 2, 6]
# print dat

res = Haar(dat)
# print res

dat2 = invHaar(res)
# print dat2

