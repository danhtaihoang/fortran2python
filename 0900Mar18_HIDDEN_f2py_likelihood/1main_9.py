##========================================================================================
## 2018.01.23: Network reconstruction with latent variables
## 2018.02.28: speed up by using multiprocessing
## 2018.03.01: check: update hidden spin in parallel (?)
##========================================================================================
import numpy as np
import sys
import timeit
from scipy import linalg
import multiprocessing
import fortran_lib as ftlib

np.random.seed(1)
generate_data = 'YES'
save_Shidden = 'NO'
nPC = 8

print("There are %d CPUs available on this machine" % multiprocessing.cpu_count())
print("You are using %d CPU(s)" % nPC)

nrepeat = 5   # update hidden spin
nloop = 500     # predict W

# parameters:
#g=2.0 ; N0=20 ; Nh = 6 ; Nb = Nh ; fN2=4.0
g = sys.argv[1] ; N0 = sys.argv[2] ; Nh = sys.argv[3] ; Nb = sys.argv[4] ; ln2 = sys.argv[5]
print(g, N0, Nh, Nb, ln2)
g = float(g) ; N0 = int(N0) ; Nh = int(Nh) ;  Nb = int(Nb) ; ln2=float(ln2)
L0 =int(ln2*N0**2)+1
L = L0-1

#Nh = 10 ; Nb = Nh
N = N0-Nh
N2 = N+Nb

if Nb == 0:
    nrepeat = 1

##========================================================================================
ext_name = '%02d_%02d_%02d_%.2f.dat'%(N0,Nh,Nb,ln2)

W_out = open('W/W_%s'%ext_name,'w')
H0_out = open('W/H0_%s'%ext_name,'w')
cost_out = open('cost/cost_%s'%ext_name,'w')
MSE_out = open('cost/MSE_%s'%ext_name,'w')
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
W=np.empty((N2, N2)) ; cost = np.empty(N2) ; H0 = np.zeros(N2)
def coupling_inference(i0):
    s = S[0:L,:].copy() ; s1 = S[1:,i0].copy()
    c_inv = C_inv.copy() ; ds = dS[0:L,:].T
    h = s1.copy() # initial value
    w1=np.empty((nloop,N2)) ; h0=np.empty(nloop)
    iloop=1 ; stop_iloop=0 ; cost1 = np.full(nloop+1,100.)
    while iloop < nloop and stop_iloop == 0:
        h_av = np.mean(h)
        hs_av = np.matmul(ds,h-h_av)/L
        w = np.matmul(hs_av,c_inv)
        h0[iloop] = h_av-np.sum(w*M)
        h[0:L] = h0[iloop]+np.matmul(s,w)

        cost1[iloop] = np.mean((s1 - np.tanh(h))**2)
        stop_iloop = (np.sign(cost1[iloop]-cost1[iloop-1])+1)/2

        h *= s1/np.tanh(h)
        w1[iloop,:] = w[:]

        iloop +=1

    niter=iloop-2

    W[i0,:] = w1[niter,:]
    cost[i0] = cost1[niter]
    H0[i0]=h0[niter]

    return W[i0,:],cost[i0],H0[i0]
##========================================================================================
## update hidden:
# 2018.03.01: parallel
##========================================================================================
##========================================================================================
##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
def update_hidden():  ## NOT parallel
    n = N ; n2 = N2 ; h0 = H0.copy() ; w=W.copy() ; s = S.copy()
    h11=np.empty(n2); p11=np.empty(n2); p12=np.empty(n2)
    h12=np.empty((n2, n2))
    for t in range(2, L):
        # P(S_hidden(t)):
        h11[n:n2]=h0[n:n2]+np.sum(w[n:n2, :]*s[t-1, :], axis=1)
        p11[n:n2]=1/(1+np.exp(-2*h11[n:n2]))
        p12[n:n2]=1-p11[n:n2]

        h120=h0+np.sum(w[:,0:n]*s[t,0:n],axis=1)

        for i in range(n,n2):
            s[t,i]=1
            h12[:,i]=h120+np.sum(w[:,n:n2]*s[t,n:n2],axis=1)
            p1=p11[i]/np.prod(1+np.exp(-2*s[t+1,:]*h12[:,i]))
            p2=p12[i]/np.prod(1+np.exp(-2*s[t+1,:]*(h12[:,i]-2*w[:,i])))
            s[t,i]=np.sign(p1/(p1+p2)-np.random.rand())
    S[:,N:N2] = s[:,N:N2]
##========================================================================================
def update_hidden_no_parallel():  ## NOT parallel
    h11=np.empty(N2);
    p11=np.empty(N2);
    p12=np.empty(N2)
    h12=np.empty((N2, N2))
    for t in range(2, L):
        # P(S_hidden(t)):
        h11[N:N2]=H0[N:N2]+np.sum(W[N:N2, :]*S[t-1, :], axis=1)
        p11[N:N2]=1/(1+np.exp(-2*h11[N:N2]))
        p12[N:N2]=1-p11[N:N2]

        h120=H0[0:N2]+np.sum(W[0:N2, 0:N]*S[t, 0:N], axis=1)

        for i in range(N, N2):
            S[t, i]=1
            h12[:, i]=h120+np.sum(W[:, N:N2]*S[t, N:N2], axis=1)
            p1=p11[i]/np.prod(1+np.exp(-2*S[t+1, :]*h12[:, i]))
            p2=p12[i]/np.prod(1+np.exp(-2*S[t+1, :]*(h12[:, i]-2*W[:, i])))
            S[t, i]=np.sign(p1/(p1+p2)-np.random.rand())
#---------------------------------------------------------------
def update_hidden_parallel():  # parallel
    h11=np.empty(N2) ; h12=np.empty((N2,N2)) ; p1 = np.empty(N2) ; p2 = np.empty(N2)
    for t in range(2,L):
        # p here is inverse of probability (1/p)
        h11[N:N2]=H0[N:N2]+np.sum(W[N:N2,:]*S[t-1,:],axis=1)
        h12[:,N:N2]= (H0[:]+np.transpose(np.sum(W[:,:]*S[t,:],axis=1))).reshape((-1,1))\
                     +W[:,N:N2]*(1-S[t,N:N2])
        p1[N:N2]=(1+np.exp(-2*h11[N:N2]))*np.prod(1+np.exp(-2*S[t+1,:]*(h12[:,N:N2].T)),axis=1)
        p2[N:N2]=(1+np.exp(2*h11[N:N2]))*np.prod(1+np.exp(-2*S[t+1,:]*(h12[:,N:N2]-2*W[:,N:N2]).T),axis=1)

        S[t,N:N2]=np.sign(p2[N:N2]/(p1[N:N2]+p2[N:N2])-np.random.rand(Nb))

##========================================================================================
## hidden_cordinate
##========================================================================================
MSE_final=np.empty(5) ; slope_final=np.empty(5) ; Wfinal = np.empty((N0,N0)) ; costS = np.empty(N0)
def hidden_cordinate():
    Wfinal_out = open('Wfinal/Wfinal_%s'%ext_name,'w')
    Woo_out=open('Wfinal/Woo_%s'%ext_name, 'w')
    Woh_out=open('Wfinal/Woh_%s'%ext_name, 'w')
    Who_out=open('Wfinal/Who_%s'%ext_name, 'w')
    Whh_out=open('Wfinal/Whh_%s'%ext_name, 'w')

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
            i_tab[i] = j1   ; i_sign[i] = 1 ; costS[i] = cost1[i,j1]
        else:
            i_tab[i] = j2 ; i_sign[i] = -1 ; costS[i] = cost2[i,j2]

        print(i,i_tab[i],i_sign[i])

    # correct W:
    for i in range(N0):
        for j in range(N0):
            Wfinal[i,j]= W[int(i_tab[i]), int(i_tab[j])]*i_sign[i]*i_sign[j]

    #print(W_all)
    #print(Wfinal)
    # all:
    MSE_final[0] = np.mean((W0[0:N0,0:N0]-Wfinal[0:N0,0:N0])**2)
    slope_final[0]= np.sum(W0[0:N0,0:N0]*Wfinal[0:N0,0:N0])/np.sum(W0[0:N0,0:N0]**2)
    # obs --> obs:
    MSE_final[1] = np.mean((W0[0:N,0:N]-Wfinal[0:N,0:N])**2)
    slope_final[1]= np.sum(W0[0:N,0:N]*Wfinal[0:N,0:N])/np.sum(W0[0:N,0:N]**2)
    # hidden --> obs:
    MSE_final[2] = np.mean((W0[0:N,N:N0]-Wfinal[0:N,N:N0])**2)
    slope_final[2]= np.sum(W0[0:N,N:N0]*Wfinal[0:N,N:N0])/np.sum(W0[0:N,N:N0]**2)
    # obs --> hidden:
    MSE_final[3] = np.mean((W0[N:N0,0:N]-Wfinal[N:N0,0:N])**2)
    slope_final[3]= np.sum(W0[N:N0,0:N]*Wfinal[N:N0,0:N])/np.sum(W0[N:N0,0:N]**2)
    # hidden --> hidden:
    MSE_final[4] = np.mean((W0[N:N0,N:N0]-Wfinal[N:N0,N:N0])**2)
    slope_final[4]= np.sum(W0[N:N0,N:N0]*Wfinal[N:N0,N:N0])/np.sum(W0[N:N0,N:N0]**2)

    fS_correct = 1-np.mean(costS[N:N0])/(2*L)
    print('MSE_final:',MSE_final)
    print('slope_final:',slope_final)
    print('fS_correct:',fS_correct)

    res_final = [[MSE_final[i],slope_final[i]] for i in range(5)]
    res_final=np.reshape(res_final,(10,1)).T

    np.savetxt('Wfinal/MSEfinal_%s'%ext_name,res_final,fmt= '% f')
    np.savetxt('Wfinal/fScorrect_%s'%ext_name,(N0,Nh,Nb,ln2,fS_correct),fmt='% f',newline='')

    # all:
    for i in range(N0):
        for j in range(N0):
            Wfinal_out.write("% i % i % f % f \n" % (i+1, j+1, W0[i,j], Wfinal[i,j]))
    Wfinal_out.close()

    # obs --> obs
    for i in range(N):
        for j in range(N):
            Woo_out.write("% i % i % f % f \n" % (i+1, j+1, W0[i,j], Wfinal[i,j]))
    Woo_out.close()

    # hidden --> obs
    for i in range(N):
        for j in range(N,N0):
            Who_out.write("% i % i % f % f \n" % (i+1, j+1, W0[i,j], Wfinal[i,j]))
    Who_out.close()

    # obs --> hidden
    for i in range(N,N0):
        for j in range(N):
            Woh_out.write("% i % i % f % f \n" % (i+1, j+1, W0[i,j], Wfinal[i,j]))
    Woh_out.close()

    # hidden --> hidden
    for i in range(N,N0):
        for j in range(N,N0):
            Whh_out.write("% i % i % f % f \n" % (i+1, j+1, W0[i,j], Wfinal[i,j]))
    Whh_out.close()

##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## main program:
##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## read data:
if generate_data=='YES':
   print('generate NEW data')
   W0 = np.random.normal(0.0,g/np.sqrt(N0),size=(N0, N0))
   generate_config(W0,L0)
   #np.savetxt('S0.dat',S0,fmt='% 2d')
   #np.savetxt('W0.dat',W0,fmt='% 2.5f')

else:
   print('read experiment data')
   W0=np.ones((N0,N0))
   #S0=np.loadtxt('S0.txt',unpack=True)
   S0=np.loadtxt('Gasper_sort-after-clean.txt',unpack=True)
   S0=S0.T
#------------------------------------------------------------------------------------------

S=np.ones((L0,N2))
S[0:L0,0:N] = S0[0:L0,0:N] # observed sequences
S1=np.copy(S[1:])

#initial hidden:
if Nb>0:
    S[0:L0,N:N2] = np.sign(np.random.rand(L0, Nb)-0.5)
#print(S)

## observed part and initial hidden part:
M = np.mean(S,axis=0)
dS = S - M
C = np.cov(dS,rowvar=False,bias=True)
C_inv = linalg.inv(C)

start_time = timeit.default_timer()

cost_repeat = np.empty(nrepeat) ; cost_repeat_oo = np.empty(nrepeat)
for irepeat in range(nrepeat):
    if Nb>0:
        # hidden part:
        dS[:,N:N2] = S[:,N:N2]-np.mean(S[:,N:N2],axis=0)
        C = np.cov(dS, rowvar=False, bias=True)
        C_inv = linalg.inv(C)

    #coupling_inference:
    if (nPC>1):
        pool=multiprocessing.Pool(processes=nPC)
        res=pool.map(coupling_inference,list(range(N2)))
        pool.close()

        W=np.array([r[0] for r in res])
        cost=np.array([r[1] for r in res])
        H0=np.array([r[2] for r in res])
    else:
        [coupling_inference(i0) for i0 in range(N2)]

    cost_repeat[irepeat] = np.mean(cost[0:N2])
    cost_repeat_oo[irepeat] = np.mean(cost[0:N])

    MSE_oo = np.mean((W0[0:N,0:N]-W[0:N,0:N])**2)
    slope_oo = np.sum(W0[0:N,0:N]*W[0:N,0:N])/np.sum(W0[0:N,0:N]**2)

    if irepeat%1 == 0:
        print(irepeat,cost_repeat[irepeat],MSE_oo,slope_oo)

    cost_out.write("%i %f %f %f %f \n" % (irepeat,cost_repeat[irepeat],cost_repeat_oo[irepeat],\
            MSE_oo,slope_oo))

##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    #start_time=timeit.default_timer()

    if Nb>0:
        S_out=ftlib.update_hidden(n1=N+1, w=W.T, h0=H0, s=S) # use fortran lib
        S[:,N:N2]=S_out[:,0:Nb]
        #update_hidden()  # use python def

##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##========================================================================================
stop_time=timeit.default_timer()
run_time=stop_time-start_time
print('run_time:', run_time)
np.savetxt('W/time_%s'%ext_name,(run_time,run_time/86400.,run_time/3600.,run_time/60.),\
            fmt='% 8.3f', newline='')
##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

if Nh == Nb and Nh > 0:
    hidden_cordinate()

if save_Shidden == 'YES':
   Shidden = S[:,N:N2].copy()
   np.savetxt('W/Shidden_%s'%ext_name,Shidden,fmt='% 2d')

##========================================================================================
MSE_out.write("%i %i %i %f %f %f %f %f \n" % (N0,Nh,Nb,ln2,cost_repeat[irepeat],\
cost_repeat_oo[irepeat],MSE_oo,slope_oo))
for i in range(N2):
    for j in range(N2):
        if generate_data == 'YES':
            W_out.write("% i % i % f % f \n" % (i+1,j+1,W[i,j],W0[i,j]))
        else:
            W_out.write("% i % i % f \n"%(i+1,j+1,W[i,j]))
W_out.close()

for i in range(N2):
    H0_out.write("% i % f \n"%(i+1,H0[i]))
H0_out.close()

##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
