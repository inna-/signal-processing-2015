from numpy import array, zeros, argmin, inf
from numpy.linalg import norm

##########################
### PYTHON DTW LIBRARY ###
##########################
def dtw(x, y):                                      # B.2
    dist=lambda x, y: norm(x - y, ord=1)

    x = array(x)
    if len(x.shape) == 1:
        x = x.reshape(-1, 1)
    y = array(y)
    if len(y.shape) == 1:
        y = y.reshape(-1, 1)

    r, c = len(x), len(y)

    D = zeros((r + 1, c + 1))
    D[0, 1:] = inf
    D[1:, 0] = inf

    for i in range(r):
        for j in range(c):
            D[i+1, j+1] = dist(x[i], y[j])

    for i in range(r):
        for j in range(c):
            D[i+1, j+1] += min(D[i, j], D[i, j+1], D[i+1, j])

    D = D[1:, 1:]
    dist = D[-1, -1] / sum(D.shape)

    return dist, D, _trackeback(D)

def _trackeback(D):
    i, j = array(D.shape) - 1
    p, q = [i], [j]
    while (i > 0 and j > 0):
        tb = argmin((D[i-1, j-1], D[i-1, j], D[i, j-1]))

        if (tb == 0):
            i = i - 1
            j = j - 1
        elif (tb == 1):
            i = i - 1
        elif (tb == 2):
            j = j - 1

        p.insert(0, i)
        q.insert(0, j)

    p.insert(0, 0)
    q.insert(0, 0)
    return (array(p), array(q))

# Attempt 1
def dtw2(z1, z2):
    n, m = len(z1), len(z2)
    costs = [[0 for x in range(m)] for y in range(n)]
 
    costs[0][0] = abs(z1[0]-z2[0])
    for i in range(1, n):
        costs[i][0] = costs[i-1][0] + abs(z1[i]-z2[0])
 
    for j in range(1, m):
        costs[0][j] = costs[0][j-1] + abs(z1[0]-z2[j])
 
    for i in range(1, n):
        for j in range(1, m):
            temp = costs[i-1][j], costs[i][j-1], costs[i-1][j-1]
            costs[i][j] = min(temp) + abs(z1[i]-z2[j])
    return costs[-1][-1]

# Attempt 2
def dtw3(x, y):
    dtw = []
    for y_i in range(len(y)):
        dtw.append([0]*len(x))
    for i in range(len(x)):
        dtw[i][0] = inf
    for j in range(len(y)):
        dtw[0][j] = inf
    dtw[0][0] = 0

    for i in range(1, len(x)):
        for j in range(1, len(y)):
            cost = norm(x[i] - y[j])
            dtw[i][j] = cost + min(dtw[i-1][j], dtw[i][j-1], dtw[i-1][j-1])
    return dtw[-1][-1]

