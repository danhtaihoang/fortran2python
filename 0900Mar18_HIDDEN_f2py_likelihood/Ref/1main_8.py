##========================================================================================
## 2018.01.23: Network reconstruction with latent variables
##========================================================================================
import numpy as np
import sys
import timeit
from scipy import linalg

np.random.seed(1)

generate_data = 'YES'

nrepeat = 10   # update hidden spin
nloop = 1000     # predict W

# parameters:
g=2.0 ; N0=300 ; fN2=2.0

#g = sys.argv[1] ; fN2 = sys.argv[2]
print(g,N0,fN2)
g = float(g) ; fN2=float(fN2)
L0 = int(fN2*N0**2)+1
L = L0-1

Nh = 0 ; Na = Nh
N = N0-Nh ; N2 =N + Na
print('N2:',N2)

if Na ==0:
    nrepeat = 1

##========================================================================================

#file_out = 'W_FEM_%.01f_%0.2f.dat'%(g,fL)
#W_out = open(file_out,'w')
#file_out = 'MSE_FEM_%.01f_%0.2f.dat'%(g,fL)
#MSE_out = open(file_out,'w')
W_out = open('W.dat','w')
H0_out = open('H0.dat','w')
MSE_out = open('MSE.dat','w')
cost_out = open('cost_iter.dat','w')
Wtest_out = open('Wtest.dat','w')
##========================================================================================
# generate sequences:
##========================================================================================
S0 = np.ones((L0, N0))
def generate_config(W0_all,L0):
   for t in range(L0-1):
      H = np.sum(W0_all[:,:]*S0[t,:],axis=1)
      P = 1/(1+np.exp(-2*H))
      S0[t+1,:]= np.sign(P-np.random.rand(N0))
   return S0
##========================================================================================
## predict W:
##========================================================================================
W_all=np.empty((N2,N2))
#def coupling_inference(N,S1,C_inv):
def coupling_inference():
    S1 = np.copy(S[1:]) ; H = S1.copy() # initial value
    iloop = 1 ; cost_tmp = 100.
    W=np.zeros((N2,N2)) ; stop_iter=np.zeros(N2)

    while iloop < nloop and np.mean(stop_iter)<1:
        H_av = np.mean(H,axis=0)                   # (n,)
        HS = np.matmul(dS[0:L,:].T, H-H_av.T)/L     # (N,L) x (L,n) = (N,n)
        w_tmp = np.matmul(HS.T, C_inv)               # (n,N) x (N,N) = (n,N)

        H0 = H_av - np.matmul(M,w_tmp.T)
        H = np.matmul(w_tmp,S[0:L,:].T)
        H = H0+H.T

        cost=np.mean((S1-np.tanh(H))**2,axis=0)
        stop_iter=np.ceil((np.sign(np.around(cost-cost_tmp,decimals=8))+1)/2)

        #for i0 in range(N2):
        #    W[i0,:]=stop_iter[i0]*W[i0,:]+(1-stop_iter[i0])*w_tmp[i0,:]

        stop_iter = stop_iter.reshape((N2,1))
        W=stop_iter*(W-w_tmp)+w_tmp

        H*=S1/np.tanh(H)

        cost_tmp = cost

        iloop += 1

    print('iloop:',iloop-1)
    #niter=iloop-2
    #print('i0:',i0,'niter:',niter)

    MSE = np.mean((W0_all-W)**2)
    slope =np.sum(W0_all*W)/np.sum(W0_all**2)
    print('MSE, slope:',np.mean(cost),MSE,slope)

    return W,np.mean(cost)
##========================================================================================
## update hidden:
##========================================================================================
def update_hidden():
    H1 = np.empty(N2) ; P1 = np.empty(N2) ; H2 = np.empty(N2) ; P2 = np.empty(N2)
    for t in range(2,L):
        # P(S_hidden(t)):
        H12 = H0_all+np.sum(W_all[:,:]*S[t-1,:],axis=1)
        P12 = 1/(1+np.exp(-2*S[t,:]*H12))
        for i in range(N,N2):
            H1[:] = H0_all + np.sum(W_all[:,:]*S[t,:],axis=1)
            P1[:] = 1/(1+np.exp(-2*S[t+1,:]*H1[:]))

            H2[:] = H1[:] - 2*W_all[:,i]*S[t,i]
            P2[:] = 1/(1+np.exp(-2*S[t+1,:]*H2[:]))

            P11 = np.prod(P1)*P12[i] ; P21 = np.prod(P2)*(1-P12[i])

            S[t,i] *= -np.sign(P21/(P21+P11)-np.random.rand())

##========================================================================================
## hidden_cordinate
##========================================================================================
MSE_final=np.empty(5) ; slope_final=np.empty(5) ; Wfinal = np.empty((N0,N0)) ; costS = np.empty(N0)
def hidden_cordinate():
    Wfinal_out = open('Wfinal.dat', 'w')
    cost1=np.empty((N0,N0)) ; cost2=np.empty((N0,N0)) ; i_tab=np.empty(N0) ; i_sign=np.ones(N0)

    for i in range(N0):
        i_tab[i]=i

    for i in range(N,N0):
        for j in range(N,N0):
            cost1[i,j]=np.sum(np.abs(S0[:L,i]-S[:L,j]))
            cost2[i,j]=np.sum(np.abs(S0[:L,i]+S[:L,j]))
        j1 = N+np.argmin(cost1[i,N:N0])
        j2 = N+np.argmin(cost2[i,N:N0])

        if cost1[i,j1]<cost2[i,j2]:
            i_tab[i]=j1   ; i_sign[i]=1 ; costS[i] = cost1[i,j1]
        else:
            i_tab[i] = j2 ; i_sign[i] = -1 ; costS[i] = cost2[i,j2]

        print(i,i_tab[i],i_sign[i])

    # correct W:
    for i in range(N0):
        for j in range(N0):
            Wfinal[i,j]=W_all[int(i_tab[i]),int(i_tab[j])]*i_sign[i]*i_sign[j]

    #print(W_all)
    #print(Wfinal)
    # all:
    MSE_final[0] = np.mean((W0_all[0:N0,0:N0]-Wfinal[0:N0,0:N0])**2)
    slope_final[0]=np.sum(W0_all[0:N0,0:N0]*Wfinal[0:N0,0:N0])/np.sum(W0_all[0:N0,0:N0]**2)
    # obs --> obs:
    MSE_final[1] = np.mean((W0_all[0:N,0:N]-Wfinal[0:N,0:N])**2)
    slope_final[1]=np.sum(W0_all[0:N,0:N]*Wfinal[0:N,0:N])/np.sum(W0_all[0:N,0:N]**2)
    # hidden --> obs:
    MSE_final[2] = np.mean((W0_all[0:N,N:N0]-Wfinal[0:N,N:N0])**2)
    slope_final[2]=np.sum(W0_all[0:N,N:N0]*Wfinal[0:N,N:N0])/np.sum(W0_all[0:N,N:N0]**2)
    # obs --> hidden:
    MSE_final[3] = np.mean((W0_all[N:N0,0:N]-Wfinal[N:N0,0:N])**2)
    slope_final[3]=np.sum(W0_all[N:N0,0:N]*Wfinal[N:N0,0:N])/np.sum(W0_all[N:N0,0:N]**2)
    # hidden --> hidden:
    MSE_final[4] = np.mean((W0_all[N:N0,N:N0]-Wfinal[N:N0,N:N0])**2)
    slope_final[4]=np.sum(W0_all[N:N0,N:N0]*Wfinal[N:N0,N:N0])/np.sum(W0_all[N:N0,N:N0]**2)
    #print(cost[0,0,1])

    fS_correct = 1-np.mean(costS[N:N0])/(2*L)
    print(MSE_final,slope_final,fS_correct)

    for i in range(N2):
        for j in range(N2):
            Wfinal_out.write("%i %i %f %f \n" % (i + 1, j + 1, W0_all[i, j], Wfinal[i, j]))
    Wfinal_out.close()

##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## main program:
##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## read data:
if generate_data=='YES':
   print('generate NEW data')
   W0_all = np.random.normal(0.0, g/np.sqrt(N0), size=(N0,N0))
   generate_config(W0_all,L0)
else:
   print('read experiment data')
   W0_all=np.ones((N0,N0))
   #S0=np.loadtxt('S0.txt',unpack=True)
   S0=np.loadtxt('Gasper_sort-after-clean.txt',unpack=True)
   S0=S0.T
#------------------------------------------------------------------------------------------
start_time = timeit.default_timer()

S=np.ones((L0,N2))
S[:,0:N] = S0[:,0:N] # observed sequences
#initial hidden:
if Na>0:
    S[0:L0,N:N2] = np.sign(np.random.rand(L0,Na)-0.5)
#print(S)

## observed part and initial hidden part:
M = np.mean(S,axis=0)
dS = S - M
C = np.cov(dS,rowvar=False,bias=True)
C_inv = linalg.inv(C)

cost_repeat = np.empty(nrepeat)
for irepeat in range(nrepeat):
    if Na>0:
        # hidden part:
        dS[:,N:N2] = S[:,N:N2]-np.mean(S[:,N:N2],axis=0)
        C = np.cov(dS, rowvar=False, bias=True)
        C_inv = linalg.inv(C)

    W_all,cost=coupling_inference()
    cost_repeat[irepeat] = cost

    if irepeat % 2 == 0:
        MSE_oo = np.mean((W0_all[0:N,0:N] - W_all[0:N,0:N])**2)
        slope_oo = np.sum(W0_all[0:N,0:N]*W_all[0:N,0:N])/np.sum(W0_all[0:N,0:N]**2)
        MSE = np.mean((W0_all - W_all)**2)
        slope = np.sum(W0_all*W_all)/np.sum(W0_all**2)
        print(irepeat,cost_repeat[irepeat],MSE_oo,slope_oo)

    cost_out.write("%i %f %f %f \n" % (irepeat,cost_repeat[irepeat],MSE_oo,slope_oo))

    if Nh>0:
        update_hidden()

if Nh == Na and Nh > 0:
    hidden_cordinate()
##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##========================================================================================
stop_time = timeit.default_timer()
run_time = stop_time - start_time
print('run_time:',run_time)

MSE_out.write("%f %f %f %f %f %f \n" % (g,float(L)/(N**2),MSE_oo,slope_oo,MSE,slope))
print (MSE_oo,slope_oo,MSE,slope)

for i in range(N2):
    for j in range(N2):
        W_out.write("%i %i %f %f \n" % (i+1, j+1, W0_all[i,j], W_all[i,j]))
W_out.close()

#for i in range(N2):
#    H0_out.write("%i %f \n" % (i+1, H0_all[i]))
#H0_out.close()

##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
