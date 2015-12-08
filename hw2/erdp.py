def edrp1(z1, z2, m , n):                         # B.3   recursive
    if m==0:        
        return n
    if n==0:        
        return m
 
    if z1[m-1]==z2[n-1]:
        return edrp1(z1,z2,m-1,n-1)
 
    return 1 + min(edrp1(z1, z2, m, n-1),          # Insert
                   edrp1(z1, z2, m-1, n),          # Remove
                   edrp1(z1, z2, m-1, n-1)         # Replace
                   )    

def edrp2(z1, z2):                                # B.3   dynamic programming
    n = len(z1)
    m = len(z2)
    temp = [[0 for x in range(m+1)] for x in range(n+1)]
    for i in range(n):
        temp[i][0] = i
    for i in range(m):
        temp[0][i] = i

    for i in range(1,n+1):        
        for j in range(1,m+1):
            if z1[i-1] == z2[j-1]:
                temp[i][j] = temp[i-1][j-1]
            else:
                temp[i][j] = min(temp[i-1][j-1], temp[i-1][j], temp[i][j-1]) + 1
    return temp[n][m]
