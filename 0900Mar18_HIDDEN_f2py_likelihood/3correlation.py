##=============================================================================
## 2018.02.23: Correlation
##=============================================================================
import numpy as np
import sys
import timeit

np.random.seed(1)

#=============================================================================

N=100
fN2=2
nens=10
#L0 = int(fN2*N**2)+1
L0=943*nens
L=L0-1

H0=np.loadtxt('H0.dat',unpack=True)
H0 = H0.T
H0 = H0[:,1]
#print(H0)

W=np.loadtxt('W.dat',unpack=True)
W = W.T
W = W[:,2].reshape((N,N))
#print(W)

S0=np.loadtxt('Gasper_sort-after-clean.txt',unpack=True)
S0=S0.T
S0 = S0[0:L0,0:N]
S10 = np.copy(S0[1:])

def cross_cov(a,b):
   da = a - np.mean(a, axis=0)
   db = b - np.mean(b, axis=0)
   return np.matmul(da.T,db)/a.shape[0]

C0 = cross_cov(S10,S0[:L])
#print(C0)

start_time = timeit.default_timer()

S = np.ones((L0, N))
def generate_config(H0,W,L0):
   for t in range(L0-1):
      H = H0+ np.sum(W[:,:]*S[t,:],axis=1)
      P = 1/(1+np.exp(-2*H))
      S[t+1,:]= np.sign(P-np.random.rand(N))
   return S
   
generate_config(H0,W,L0)

S1 = np.copy(S[1:])
C = cross_cov(S1,S[:L])
#print(C)

stop_time = timeit.default_timer()
run_time = stop_time - start_time
print('run_time:',run_time)

C_out = open('C.dat','w')
for i in range(N):
    for j in range(N):
        C_out.write("%i %i %f %f \n" % (i+1, j+1, C0[i,j], C[i,j]))
C_out.close()

MSE = np.mean((C-C0)**2)
print('MSE:',MSE)



