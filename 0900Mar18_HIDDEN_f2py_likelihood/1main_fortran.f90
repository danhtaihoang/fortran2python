!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!! Danh-Tai Hoang - LBM/NIDDK/NIH 
!-----------------------------------------------------------------------------------------
!! Network Reconstruction: Kinectic Ising Model with Hidden Spin
!!======================================================================================== 
!! 2017.06.26: flip hidden spins to get maximum probability
!! 2017.07.06: hidden cordinate
!! 2018.01.25: speed up the code
!! 2018.01.25: With H0
!! 2018.02.02: simulation
!! 2018.02.15: circle + block
!! 2018.02.16: parallel_i0 (for N=500, 1000)
!! 2018.02.22: speed up update
!! 2018.03.05: Gasper with ensembles
!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   PROGRAM main  
   IMPLICIT NONE

   INTEGER (KIND=8),PARAMETER :: seed=1,nloop=300,nloop2=1,ns0=943
   INTEGER (KIND=8) :: nrepeat=4000
   CHARACTER (LEN=50),PARAMETER :: NEW_CONFIG='NO'
   REAL (KIND=8),PARAMETER :: Wmax=0.2

   INTEGER (KIND=8) :: i,j,k,na0,na,ns,i0,na1,irepeat,na2,n_iter,stop_iteration,iloop
   INTEGER (KIND=8) :: nb,nh,na3,iloop2,iens,ns_nens,nens
   REAL (KIND=8) :: g,rdn_S,P11,P21,fs_correct,ln2,H_av
   REAL (KIND=8) :: time_2,time_3,time_5,time_6,time_7

   CHARACTER (LEN=50) :: arg_g,arg_na0,arg_nh,arg_ln2,arg_nb,arg_nens
   CHARACTER (LEN=50) :: name_g,name_na0,name_nh,name_ln2,name_nb,name_nens
   CHARACTER (LEN=200) :: file_name,sequence,sequence_S0
   
   INTEGER (KIND=8),DIMENSION(50) :: seed1
   REAL (KIND=8),DIMENSION(:),ALLOCATABLE :: W0,W2,cost,MSE,slope,cost_repeat1,H12,P12,M,HS
   REAL (KIND=8),DIMENSION(:,:,:),ALLOCATABLE :: S,S0,dS
   REAL (KIND=8),DIMENSION(:,:),ALLOCATABLE :: W0_all,W_all,W1,H0_repeat,H,St1
   REAL (KIND=8),DIMENSION(:),ALLOCATABLE :: H1,H2,cost_repeat,cost1,W,H0,H0_all
   REAL (KIND=8),DIMENSION(:,:),ALLOCATABLE :: Wpred,Wfinal
   REAL (KIND=8),DIMENSION(:,:,:),ALLOCATABLE :: W_repeat
   REAL (KIND=8),DIMENSION(:,:),ALLOCATABLE :: C,C_inv
   REAL (KIND=8),DIMENSION(5) :: MSE_final,slope_final
  
!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
!!========================================================================================
!!========================================================================================
   !CALL system('rm *.dat*')     
   !!======================================
   !! get parameters from argument:
   CALL getarg(1,arg_na0)
   READ(arg_na0, *)na0

   CALL getarg(2,arg_g)
   READ(arg_g, *)g

   CALL getarg(3,arg_nh)
   READ(arg_nh, *)nh

   CALL getarg(4,arg_nb)
   READ(arg_nb, *)nb

   CALL getarg(5,arg_ln2)
   READ(arg_ln2, *)ln2

   CALL getarg(6,arg_nens)
   READ(arg_nens, *)nens
  
   seed1(1) = seed
   CALL RANDOM_SEED(PUT=seed1)         

   !!!------------------------
   !! set parameter:  
   na=na0-nh ; na1=na+1  ; na2=na+nb  ; na3=min(na0,na2)
   !ns=int(ln2*(na0**2)) ; ns0=ns+1
   ns=ns0-1
   ns_nens=ns*nens

   IF (nb==0) THEN
      nrepeat=1
   END IF     
 
   !!!------------------------
   WRITE(*,*)'na0,na,nb,ln2,nens:'
   WRITE(*,'(3I6,2X,F5.2,I6)')na0,na,nb,ln2,nens

   ALLOCATE(W0(na0),S(nens,ns0,na2),S0(nens,ns0,na0),W0_all(na0,na0),dS(nens,ns,na2))
   ALLOCATE(St1(nens,ns),H0(nloop),H0_all(na2),H0_repeat(na2,nrepeat))
   ALLOCATE(W2(na2),W_all(na2,na2),cost(na2))
   ALLOCATE(MSE(nrepeat),slope(nrepeat))
   ALLOCATE(W_repeat(na2,na2,nrepeat),H1(na2),H2(na2))
   ALLOCATE(cost_repeat(nrepeat),Wpred(na2,na2),cost_repeat1(nrepeat))
   ALLOCATE(Wfinal(na2,na2),M(na2),C(na2,na2),C_inv(na2,na2))
   ALLOCATE(H(nens,ns),cost1(0:nloop),HS(na2),W(na2),W1(nloop,na2),H12(na1:na2),P12(na1:na2))

   WRITE(name_na0,'(I3.3)')na0
   WRITE(name_g,'(F3.1)')g
   WRITE(name_nh,'(I2.2)')nh
   WRITE(name_nb,'(I2.2)')nb
   WRITE(name_ln2,'(F5.2)')ln2
   WRITE(name_nens,'(I3.3)')nens

   !file_name='_'//trim(name_nh)//'_'//trim(name_nb)//'_'//trim(name_ln2)
   !file_name='_'//trim(name_na0)//'_'//trim(name_nb)//'_'//trim(name_ln2)
   file_name='_'//trim(name_na0)//'_'//trim(name_nb)//'_'//trim(name_nens)
   
   OPEN(21,file='cost/cost_iter'//trim(file_name)//'.dat')
   OPEN(22,file='cost/cost'//trim(file_name)//'.dat')   
   OPEN(23,file='cost/nit'//trim(file_name)//'.dat')
   OPEN(24,file='cost/time_run'//trim(file_name)//'.dat')
   OPEN(25,file='cost/number_hidden'//trim(file_name)//'.dat')
   !OPEN(26,file='cost/i0'//trim(file_name)//'.dat')

   OPEN(32,file='W/W'//trim(file_name)//'.dat')
   OPEN(33,file='W/H0'//trim(file_name)//'.dat')
   OPEN(34,file='W/S_hidden'//trim(file_name)//'.dat')

!!!=======================================================================================
!!! Generate a new configuration or Load a given configuration
!!!=======================================================================================

   IF (NEW_CONFIG=='YES') THEN
      WRITE(*,*)'Generate NEW coupling and config'
      CALL generate_coupling(na0,g,W0_all)
     
      !CALL generate_coupling_circle(na0,Wmax,W0_all)
      !CALL generate_coupling_checkerboard(na0,Wmax,W0_all)

      CALL generate_config(na0,ns,W0_all,S0)
      !WRITE(*,*)'S0:',S0(1,ns0,1:na0)

      IF (nb==nh) THEN
         OPEN(14,file='W/S0'//trim(file_name)//'.dat')
         WRITE(sequence_S0,'("("I4,"F4.0)")')na0
         DO iens=1,nens
         DO k=1,ns0   
            WRITE(14,sequence_S0)S0(iens,k,1:na0)
         END DO
         END DO
      END IF
      CLOSE(14)

      !IF (nb==nh) THEN  
      !OPEN(15,file='W/W0'//trim(file_name)//'.dat')
      !DO i=1,na0
      !   DO j=1,na0
      !   WRITE(15,'(2I6,F25.12)')j,i,W0_all(j,i)
      !   END DO
      !END DO  
      !END IF
      !CLOSE(15)

   ELSE
      WRITE(*,*)'Load config'      
      CALL read_Gasper_data()  
      W0_all(:,:)=0.
   END IF           
                   
   CALL initial_CPU_time()

!!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!!!=======================================================================================
!!! Predict couplings:
   CALL flip_hidden_based()

   CALL computation_time()            

   IF ((nb.eq.nh).and.(nh>0).and.(NEW_CONFIG=='YES')) THEN            
      CALL hidden_cordinate_main()
   END IF 

!!! write hidden sequences:
   IF (nb>0) THEN           
      WRITE(sequence,'("("I4,"F4.0)")')nb

      DO iens=1,nens
      DO k=1,ns0   
         WRITE(34,sequence)S(iens,k,na1:na2)
      END DO
      END DO

   END IF

!!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
      
   CLOSE(12) ; CLOSE(21) ; CLOSE(22) ; CLOSE(23) ; CLOSE(26)
   CLOSE(31) ; CLOSE(32) ; CLOSE(33) ; CLOSE(34) 

   WRITE(*,*)'Finished'                

   CONTAINS
!!!======================================================================================= 
!!!======================================================================================= 
   SUBROUTINE flip_hidden_based()
   IMPLICIT NONE   
   
   !! observed config:     
   S(1:nens,1:ns0,1:na)=S0(1:nens,1:ns0,1:na)
   
   !! initial hidden:
   S(:,:,na1:na2)=1.   
   DO iens=1,nens
   DO k=1,ns0 
      DO i=na1,na2        
         CALL random_number(rdn_S)         
         IF (rdn_S>0.5) THEN
            S(iens,k,i)=-1.
         END IF
      END DO         
   END DO
   END DO
   
   !!-------------------------------------------
   !! observed part first 
   DO i=1,na
      M(i)=sum(S(:,1:ns,i))/ns_nens
      dS(:,1:ns,i) = S(:,1:ns,i)-M(i)      
   END DO
       
   DO i=1,na
   DO j=i,na
      C(i,j)=sum(dS(:,1:ns,i)*dS(:,1:ns,j))/ns_nens
      C(j,i)=C(i,j)
   END DO
   END DO
      
   !!-------------------------------------------------------                         
   DO irepeat=1,nrepeat      
      !!-------------------------------------------
      !! hidden part
      DO i=na1,na2
         M(i)=sum(S(:,1:ns,i))/ns_nens
         dS(:,1:ns,i) = S(:,1:ns,i)-M(i)
      END DO
         
      DO i=1,na
      DO j=na1,na2
         C(i,j)=sum(dS(:,1:ns,i)*dS(:,1:ns,j))/ns_nens
         C(j,i)=C(i,j)
      END DO
      END DO

      DO i=na1,na2
      DO j=i,na2
         C(i,j)=sum(dS(:,1:ns,i)*dS(:,1:ns,j))/ns_nens
         C(j,i)=C(i,j)
      END DO
      END DO

      CALL pseudo_inverse(C,na2,na2,C_inv)
      CALL coupling_prediction()

      cost_repeat1(irepeat)=sum(cost(1:na))/na
      cost_repeat(irepeat)=sum(cost(1:na2))/na
            
      MSE(irepeat)=sum((W_all(1:na,1:na)-W0_all(1:na,1:na))**2)/na/na
      slope((irepeat))=sum(W0_all(1:na,1:na)*W_all(1:na,1:na))/sum(W0_all(1:na,1:na)**2)
            
      W_repeat(:,:,irepeat)=W_all(:,:)
      H0_repeat(:,irepeat)=H0_all(:)
      
      WRITE(21,'(I6,5F15.10)')irepeat,cost_repeat(irepeat),cost_repeat1(irepeat),&
      cost_repeat(irepeat)*na/na2,MSE(irepeat),slope(irepeat)

      IF (mod(irepeat,10)==0) THEN 
      WRITE(*,'(I6,5F15.10)')irepeat,cost_repeat(irepeat),cost_repeat1(irepeat),&
      cost_repeat(irepeat)*na/na2,MSE(irepeat),slope(irepeat)
      END IF  
     
      !! update hidden spin:
      IF (nb>0) THEN        
         CALL update_hidden()
      END IF
      
   END DO
               
   irepeat=minloc(cost_repeat(1:nrepeat),dim=1)             
   !irepeat=nrepeat !! take the last one
   Wpred(:,:)=W_repeat(:,:,irepeat)

   WRITE(22,'(F8.4,3I6,2X,5F13.9)')ns/real(na0**2),nh,nb,irepeat,cost_repeat(irepeat),&
   cost_repeat1(irepeat),cost_repeat(irepeat)*na/na2,MSE(irepeat),slope(irepeat)
                  
   DO i=1,na2
      DO j=1,na2
      IF (NEW_CONFIG =='YES') THEN
      WRITE(32,'(2I6,2F25.12)')j,i,W0_all(j,i),Wpred(j,i)
      ELSE
      WRITE(32,'(2I6,F25.12)')j,i,Wpred(j,i)
      END IF
      END DO
   END DO  

   DO i=1,na2
      WRITE(33,'(I6,F25.12)')i,H0_repeat(i,irepeat)
   END DO 
   
   END SUBROUTINE
!!========================================================================================
!! coupling prediction
!!========================================================================================
   SUBROUTINE coupling_prediction()
   IMPLICIT NONE
      
   DO i0=1,na2
      !WRITE(*,*)'i0:',i0
      !WRITE(26,*)'i0:',i0

      DO k=1,ns
         St1(:,k)=S(:,k+1,i0)
      END DO

      H(:,1:ns)=St1(:,1:ns)  !! initial value
      DO iloop2=1,nloop2
         iloop=1 ; cost1(:)=0. ; stop_iteration=0
         !!--------------------------------------------------------------   
         DO WHILE ((iloop<=nloop).and.(stop_iteration==0))
            H_av=sum(H(:,1:ns))/ns_nens
            DO i=1,na2
            HS(i)=sum((H(:,1:ns)-H_av)*dS(:,1:ns,i))/ns_nens
            END DO        
            W(:)=matmul(HS(:),C_inv(:,:))

            !IF (iloop2==1) THEN
               H0(iloop)=H_av - sum(W(:)*M(:)) 
            !ELSE
            !   IF (mod(iloop2,2)==1) THEN      
            !      H0(iloop)=0.
            !  ELSE
            !     H0(iloop)=H_av - sum(W(:)*M(:))            
            !  END IF
            !END IF

            DO iens=1,nens
               H(iens,1:ns)=matmul(S(iens,1:ns,:),W(:))+H0(iloop)                
            END DO

            IF (i0<=na) THEN
               cost1(iloop)=sum((St1(:,1:ns)-tanh(H(:,1:ns)))**2)/ns_nens
            ELSE
               cost1(iloop)=sum((St1(:,2:ns-1)-tanh(H(:,2:ns-1)))**2)/(ns-2)/nens
            END IF      
            !WRITE(*,*)i0,iloop,cost1(iloop)
                     
            IF ((cost1(iloop)>cost1(iloop-1)).and.(iloop>1)) THEN
               stop_iteration=1
            END IF

            H=St1*H/tanh(H)
            W1(iloop,:)=W(:)            

            iloop=iloop+1                 
         END DO
      END DO
      !!--------------------------------------------------------------    
      n_iter=iloop-2                          
      
      WRITE(23,'(3I8)')irepeat,i0,n_iter   
      
      W_all(:,i0)=W1(n_iter,:)
      cost(i0)=cost1(n_iter)
      H0_all(i0)=H0(n_iter)
   END DO
      
   END SUBROUTINE
!!!======================================================================================= 
!!!======================================================================================= 
   SUBROUTINE update_hidden()
   IMPLICIT NONE   

   DO iens=1,nens
      DO k=2,ns 
         H12(na1:na2)=H0_all(na1:na2)+matmul(S(iens,k-1,1:na2),W_all(1:na2,na1:na2))
         P12(na1:na2)= 1/(1+exp(-2*S(iens,k,na1:na2)*H12(na1:na2)))  
                  
         DO i=na1,na2
            H1(:)=H0_all(:)+matmul(S(iens,k,1:na2),W_all(1:na2,:))         
     
            P11=P12(i)/product(1+exp(-2*S(iens,k+1,:)*H1(:)))          
            P21=(1.-P12(i))/product(1+exp(-2*S(iens,k+1,:)*(H1(:)-2.*W_all(i,:)*S(iens,k,i))))

            CALL random_number(rdn_S)
            IF (P21/(P11+P21) > rdn_S) THEN
               S(iens,k,i)=-S(iens,k,i)
            END IF   
         END DO

      END DO
   END DO   
   
   END SUBROUTINE    

!!========================================================================================
!!========================================================================================
   SUBROUTINE read_Gasper_data()
   IMPLICIT NONE

   !OPEN(12,file='Gasper_sort_Wmax.txt')
   OPEN(12,file='Gasper_sort-after-clean.txt') 
             
   DO iens=1,nens                 
   DO k=1,ns0
      READ(12,*)(S0(iens,k,i),i=1,na0)
   END DO
   END DO
               
   END SUBROUTINE   

!!!======================================================================================= 
!!!======================================================================================= 
   SUBROUTINE hidden_cordinate_main()
   IMPLICIT NONE   

   OPEN(51,file='Wfinal/Wall'//trim(file_name)//'.dat')
   OPEN(52,file='Wfinal/Woo'//trim(file_name)//'.dat')
   OPEN(53,file='Wfinal/Who'//trim(file_name)//'.dat')
   OPEN(54,file='Wfinal/Woh'//trim(file_name)//'.dat')
   OPEN(55,file='Wfinal/Whh'//trim(file_name)//'.dat')  
   OPEN(56,file='Wfinal/MSE'//trim(file_name)//'.dat')            
      
   CALL hidden_cordinate(na0,na,W0_all,Wpred,Wfinal,fs_correct,MSE_final,slope_final)

   !!---------------------------------------------------------      
   !! all:
   DO i=1,na0
   DO j=1,na0
      WRITE(51,*)j,i,W0_all(j,i),Wfinal(j,i)
   END DO
   END DO
   
   !! obs --> obs:
   DO i=1,na
   DO j=1,na
      WRITE(52,*)j,i,W0_all(j,i),Wfinal(j,i)
   END DO
   END DO
   
   !! hidden --> obs:
   DO i=1,na
   DO j=na1,na0
      WRITE(53,*)j,i,W0_all(j,i),Wfinal(j,i)
   END DO
   END DO
   
   !! obs --> hidden:
   DO i=na1,na0
   DO j=1,na
      WRITE(54,*)j,i,W0_all(j,i),Wfinal(j,i)
   END DO
   END DO
   
   !! hidden --> hidden:
   DO i=na1,na0
   DO j=na1,na0
      WRITE(55,*)j,i,W0_all(j,i),Wfinal(j,i)
   END DO
   END DO   

   WRITE(56,'(F8.2,2I4,2X,F8.5,10F14.10)')ns/real(na0**2),nh,nb,fs_correct,(MSE_final(i),slope_final(i),i=1,5)

   CLOSE(51) ; CLOSE(52) ; CLOSE(53) ; CLOSE(54) ; CLOSE(55) ; CLOSE(56)

   END SUBROUTINE  
!!!=======================================================================================
!!!=======================================================================================
!!! 2017.07.05: Correct the order of hidden spin
!!!=======================================================================================
   SUBROUTINE hidden_cordinate(na0,na,W0,W,W1,fs_correct,MSE_final,slope_final)
   IMPLICIT NONE
   
   INTEGER (KIND=8) :: i,j,j1,j2,na0,na,na1
   REAL (KIND=8) :: fs_correct

   INTEGER (KIND=8),DIMENSION(na0) :: i_tab,i_sign
   REAL (KIND=8),DIMENSION(na0,na0) :: W0,W,W1
   REAL (KIND=8),DIMENSION(na+1:na0,na+1:na0,1:2) :: cost
   REAL (KIND=8),DIMENSION(na+1:na0) :: costS

   REAL (KIND=8),DIMENSION(5) :: MSE_final,slope_final
   
   !!-----------------------------------            
   !!! 1) ini MSE, slope 
   na1=na+1
         
   DO i=1,na
      i_tab(i)=i
      i_sign(i)=1   
   END DO
   !!-------------------------------------------------------------------------------------
   !! 1) find hidden cordinate 
   DO i=na1,na0
   DO j=na1,na0
      cost(i,j,1)=sum((W0(1:na,i)-W(1:na,j))**2)+sum((W0(i,1:na)-W(j,1:na))**2)
      cost(i,j,2)=sum((W0(1:na,i)+W(1:na,j))**2)+sum((W0(i,1:na)+W(j,1:na))**2)            
   END DO
   END DO
   
   !!---------------------
   DO i=na1,na0   
      j1=minloc(cost(i,na1:na0,1),dim=1)+na
      j2=minloc(cost(i,na1:na0,2),dim=1)+na
      
      IF (cost(i,j1,1)<cost(i,j2,2)) THEN
         i_tab(i)=j1
         i_sign(i)=1
         !WRITE(*,*)i,i_tab(i),i_sign(i),cost(i,j1,1)         
      ELSE   
         i_tab(i)=j2
         i_sign(i)=-1
         !WRITE(*,*)i,i_tab(i),i_sign(i),cost(i,j2,2)
      END IF
      WRITE(*,*)i,i_tab(i),i_sign(i)      
   END DO
            
   !!---------------------------------------------------------
   !! 2) WRITE coupling predicted:
   DO i=1,na0
   DO j=1,na0
      W1(j,i)=W(i_tab(j),i_tab(i))*i_sign(i)*i_sign(j)
   END DO
   END DO   
   !!---------------------------------------------------
   !! 3) percent correct hidden spin:
   
   DO i=na1,na0
      costS(i)=sum(abs(S0(:,1:ns,i)-S(:,1:ns,i_tab(i))*i_sign(i)))/ns/2/nens    
   END DO
   fs_correct=1.-sum(costS(na1:na0))/(na0-na)
      
   !!-----------------------------------------------------
   !! all:
   MSE_final(1)=sum((W0(1:na0,1:na0)-W1(1:na0,1:na0))**2)/(na0**2)
   slope_final(1)=sum(W0(1:na0,1:na0)*W1(1:na0,1:na0))/sum(W0(1:na0,1:na0)**2)     
   !! obs --> obs:                                      
   MSE_final(2)=sum((W0(1:na,1:na)-W1(1:na,1:na))**2)/(na**2)
   slope_final(2)=sum(W0(1:na,1:na)*W1(1:na,1:na))/sum(W0(1:na,1:na)**2)                                            
   !! hidden --> obs:                    
   MSE_final(3)=sum((W0(na1:na0,1:na)-W1(na1:na0,1:na))**2)/((na0-na)*na)  
   slope_final(3)=sum(W0(na1:na0,1:na)*W1(na1:na0,1:na))/sum(W0(na1:na0,1:na)**2) 
   !! obs --> hidden:                                                                  
   MSE_final(4)=sum((W0(1:na,na1:na0)-W1(1:na,na1:na0))**2)/((na0-na)*na)
   slope_final(4)=sum(W0(1:na,na1:na0)*W1(1:na,na1:na0))/sum(W0(1:na,na1:na0)**2) 
   !! hidden --> hidden:                                           
   MSE_final(5)=sum((W0(na1:na0,na1:na0)-W1(na1:na0,na1:na0))**2)/((na0-na)**2) 
   slope_final(5)=sum(W0(na1:na0,na1:na0)*W1(na1:na0,na1:na0))/sum(W0(na1:na0,na1:na0)**2) 
   !!-----------------------------------------------------

   END SUBROUTINE         
 
!!!=======================================================================================
!!!=======================================================================================
   SUBROUTINE generate_coupling(na,g,W_all) 
   IMPLICIT NONE

   INTEGER (KIND=8) :: na,i,j
   REAL (KIND=8) :: g,W_mean,sigma     
   REAL (KIND=8) :: beta(na**2),W_all(na,na)

   !!!====================================================================================
   W_mean=0.
   sigma=g/sqrt(real(na))
      
   CALL generate_gaussian_random(na**2,sigma,W_mean,beta)
   
   !OPEN(14,file='W0_all.dat')   
   
   DO i=1,na
      DO j=1,na
         W_all(j,i)=beta((j-1)*na+i)
         !WRITE(14,*)j,i,W_all(j,i)
      END DO
   END DO
   !CLOSE(14)
    
   END SUBROUTINE

!!!=======================================================================================
!!! 2017.03.14: Generate coupling strengths as checkerboard
!!!=======================================================================================
   SUBROUTINE generate_coupling_checkerboard(na,Wmax,W) 
   IMPLICIT NONE

   INTEGER (KIND=8),PARAMETER :: nb=2
   INTEGER (KIND=8) :: i,j,k,kj,na
   REAL (KIND=8) :: Wmax,db,rdn,Wnoise=1.0
   REAL (KIND=8) :: W(na,na),ib(0:na)

   db=real(na)/nb

   DO i=0,nb  
      ib(i)=db*i
      !WRITE(*,*)i,ib(i)
   END DO
   
   W=Wmax   
   !DO k=0,nb
   DO i=1,na
      DO j=1,na
         CALL random_number(rdn)
         W(i,j)=Wmax*(1.+(rdn-0.5)*Wnoise)      
      END DO
   END DO
   !END DO
   
   !!------------------------------------------------
   DO k=0,nb,2 
      DO i=1,na
         IF ((ib(k)<i).and.(i<=ib(k+1))) THEN
         DO kj=0,nb,2
            DO j=1,na
            IF ((ib(kj)<j).and.(j<=ib(kj+1))) THEN     

            CALL random_number(rdn)
            W(i,j)=-Wmax*(1.+(rdn-0.5)*Wnoise)
 
            ENDIF
            END DO
         END DO
         END IF
      END DO
   END DO

   DO k=1,nb,2 
      DO i=1,na
         IF ((ib(k)<i).and.(i<=ib(k+1))) THEN
         DO kj=1,nb,2
            DO j=1,na
            IF ((ib(kj)<j).and.(j<=ib(kj+1))) THEN     
               CALL random_number(rdn)
               W(i,j)=-Wmax*(1.+(rdn-0.5)*Wnoise)

            ENDIF
            END DO
         END DO
         END IF
      END DO
   END DO

   END SUBROUTINE
!!!=======================================================================================
!!! 2017.03.1: Generate coupling strengths as circles
!!!=======================================================================================
   SUBROUTINE generate_coupling_circle(na,Wmax,W) 
   IMPLICIT NONE

   INTEGER (KIND=8),PARAMETER :: nr=6 !! for N=40
   !INTEGER (KIND=8),PARAMETER :: nr=6   !! for N=100
   INTEGER (KIND=8) :: na,i,j,k,nr2
   REAL (KIND=8) :: Wmax,na1
   REAL (KIND=8) :: W(na,na),R0(na,na)
   REAL (KIND=8) :: R(2*nr+1)
   !!-----------------------------------------
   na1=(na+1)/2. ; nr2=2*nr   
   DO i=1,nr2+1
      R(i)=i*na/2./nr
      !WRITE(*,*)i,R(i)
   END DO
   
   W=Wmax
   DO i=1,na
      DO j=1,na
      R0(i,j)=sqrt((i-na1)**2.+(j-na1)**2.)
      
      DO k=1,nr2,2         
         IF ((R(k)<R0(i,j)).and.(R0(i,j)<=R(k+1))) THEN 
            W(i,j)=-Wmax/log(R(k))
         ENDIF
      END DO

      DO k=2,nr2,2         
         IF ((R(k)<R0(i,j)).and.(R0(i,j)<=R(k+1))) THEN 
            W(i,j)=Wmax/log(R(k))
         ENDIF
      END DO
      
      END DO
   END DO
   
   END SUBROUTINE  

!!========================================================================================
!!========================================================================================   
   SUBROUTINE generate_config(na,ns,W,S)
   IMPLICIT NONE   
  
   INTEGER (KIND=8) :: na,ns,i,k
   REAL (KIND=8) :: H0,P,rdn_S
   REAL (KIND=8) :: W(na,na),S(nens,ns+1,na)
   !CHARACTER (LEN=50) :: sequence
!!========================================================================================   
   !!! Set initial S
   DO iens=1,nens

      DO i=1,na 
         CALL random_number(rdn_S)
         IF (rdn_S>0.5) THEN
            S(iens,1,i)=1.
         ELSE
            S(iens,1,i)=-1.
         END IF
      END DO

   DO k=1,ns
      DO i=1,na
         
      H0=sum(W(1:na,i)*S(iens,k,1:na))                      
      P=exp(H0)/(exp(H0)+exp(-H0))
            
      CALL random_number(rdn_S)
      IF (P>rdn_S) THEN
         S(iens,k+1,i)=1.
      ELSE
         S(iens,k+1,i)=-1.
      END IF
     
      END DO

   END DO  

   WRITE(*,*)S(1,ns+1,na)
 
   END DO
   
   !OPEN(unit=16,file='train.dat')    
   !WRITE(sequence,'("("I4,"F5.1,2X,F24.20)")')na
   !DO k=1,ns   
   !   WRITE(16,sequence)S(k,1:na),H(k)
   !END DO
   
   END SUBROUTINE   

!!========================================================================================
!!========================================================================================   
   SUBROUTINE generate_gaussian_random(n,sigma,average,y1_tab)
   IMPLICIT NONE
   
   INTEGER(KIND=8) :: i,n
   REAL(KIND=8) :: rdn1,rdn2,x1,x2,y1,y2,w
   REAL(KIND=8) :: y1_av,y2_av,y12_av,y22_av,y1_dev,y2_dev,sigma,average
      
   REAL(KIND=8),DIMENSION(n) :: y1_tab,y2_tab
   
   !!!-----------------------------------
   !OPEN(unit=21,file='y.dat')
     
   DO i=1,n
      w=10.
      DO WHILE ((w.ge.1.).or.(w.eq.0.))
         CALL random_number(rdn1)
         CALL random_number(rdn2)
         x1=2.*rdn1-1.
         x2=2.*rdn2-1.
         w=x1**2.+x2**2.
      END DO
      
      w=sigma*sqrt(-2.*log(w)/w)
      y1=x1*w
      y2=x2*w

      y1_tab(i)=y1 + average
      y2_tab(i)=y2 + average
      
   END DO   

   !!!----------------
   !DO i=1,n
      !WRITE(21,*)i,y1_tab(i),y2_tab(i)
   !END DO

   !!!-------------------------------
   !!! Checking:

   y1_av=sum(y1_tab(1:n))/n ; y2_av=sum(y2_tab(1:n))/n
   y12_av=sum(y1_tab(1:n)**2)/n ; y22_av=sum(y2_tab(1:n)**2)/n
   y1_dev=sqrt(y12_av-y1_av**2.) ; y2_dev=sqrt(y22_av-y2_av**2.)
   
   !WRITE(*,*)y1_av,y2_av,y1_dev,y2_dev 
   
   !CLOSE(21)
   
   END SUBROUTINE      
   
!!========================================================================================
!!========================================================================================
!! Pseudo Inverse Matrix.
!! ALL the below program related to pseudo_inverse are downloaded at 
!!! https://people.sc.fsu.edu/~jburkardt/f_src/svd_demo/svd_demo.html
!  Licensing: This code is distributed under the GNU LGPL license.
!  Modififed: 19 June 2012
!  Author: John Burkardt
!    input: A(M,N) : the matrix whose singular value decomposition we are investigating.
!    outpout:
!    S(M,N): the diagonal factor in the singular value decomposition of A.
!    U(M,M), the first orthogonal factor in the singular value decomposition of A.
!    V(N,N), the second orthogonal factor in the singular value decomposition of A.
!
!! PSEUDO_INVERSE computes the pseudoinverse.
!  Discussion:
!    Given the singular value decomposition of a real MxN matrix A:
!      A = U * S * V'
!    where
!      U is MxM orthogonal,
!      S is MxN, and entirely zero except for the diagonal;
!      V is NxN orthogonal.
!
!    the pseudo inverse is the NxM matrix A+ with the form
!      A+ = V * S+ * U'
!    where
!      S+ is the NxM matrix whose nonzero diagonal elements are
!      the inverses of the corresponding diagonal elements of S.
!!------------------------------------------------------------------
   SUBROUTINE pseudo_inverse(a,m,n,a_pseudo)
   implicit none

   INTEGER(kind = 8) :: m,n,i,ncut
   REAL(kind = 8):: s(m,n),sp(n,m),u(m,m),v(n,n),a(m,n),a_pseudo(n,m)!,a_pseudobig(n,m)
   
   !!-----------------------
   ncut=n
   !WRITE(*,*)'ncut:',ncut
   
   CALL SVD_lapack (a,m,n,u,s,v)

   sp(1:n,1:m) = 0.
   DO i=1,min(m,n)
      IF (s(i,i) /=0.) THEN
         sp(i,i) = 1./s(i,i)
      END IF
   END DO

   !a_pseudo(1:n,1:m) = matmul(v(1:n,1:n),matmul(sp(1:n,1:m),transpose(u(1:m,1:m))))
   
   !! 2016.05.16: Tai: use only large value of s (or small of sp): -----
   u=transpose(u)
   !ncut=20
   a_pseudo(1:n,1:m) = matmul(v(1:n,1:ncut),matmul(sp(1:ncut,1:ncut),u(1:ncut,1:m)))
   
  ! a_pseudobig(1:n,1:m) = matmul(v(1:n,1:ncut1),matmul(sp(1:ncut1,1:ncut1),u(1:ncut1,1:m)))   
   !!---------------------------
   
   END SUBROUTINE 
!!-------------------------------------------------------------
   SUBROUTINE SVD_lapack (a,m,n,u,s,v)
   implicit none

   integer (kind = 8):: m,n,i,lda,ldu,ldv,lwork,info
   character jobu
   character jobv
   real(kind = 8) :: a(m,n),a_copy(m,n),sdiag(min(m,n)),s(m,n),u(m,m),v(n,n)
   real(kind = 8), allocatable, DIMENSION (:) :: work

   lwork = max(3 *min(m,n) + max(m,n),5*min(m,n))

   allocate (work(1:lwork))
   !!Compute the eigenvalues and eigenvectors:
   jobu = 'A' ; jobv = 'A' ; lda = m ; ldu = m ; ldv = n

   !! The input matrix is destroyed by the routine. Since we need to keep
   !! it around, we only pass a copy to the routine:
   a_copy(1:m,1:n) = a(1:m,1:n)

   !!! Call function dgesvd() from LAPACK library
   CALL dgesvd(jobu,jobv,m,n,a_copy,lda,sdiag,u,ldu,v,ldv,work,lwork,info)

   !IF (info /= 0) then
   ! write ( *, '(a)' ) ' '
   ! write ( *, '(a)' ) 'R8MAT_SVD_LAPACK - Failure!'
   ! write ( *, '(a)' ) '  The SVD could not be calculated.'
   ! write ( *, '(a)' ) '  LAPACK routine DGESVD returned a nonzero'
   ! write ( *, '(a,i8)' ) '  value of the error flag, INFO = ', info
   ! return
   !END IF
!
!  Make the MxN matrix S from the diagonal values in SDIAG.

   !!16.05.13: Tai : write eigenvalue
   !WRITE(31,*)'eigenvalue:'
   !WRITE(31,*)sdiag
         
   s(1:m,1:n) = 0.0D+00
   DO i = 1, min(m,n)
      s(i,i) = sdiag(i)
   END DO

   !Transpose V:
   v=transpose(v)

   !deallocate(work)

   END SUBROUTINE

!!=======================================================================================
!! Computing singular values of a matrix using LAPACK
!!!=======================================================================================
   SUBROUTINE rank_matrix(m,n,a,rank)
   IMPLICIT NONE
   INTEGER (KIND=8) :: i, m, n, lda, lds, lwork, lwmax, info, rank
   REAL (KIND=8),PARAMETER :: eps= 5e-10
   REAL (KIND=8) :: a(m,n),S(max(m,n)),U(m,m),VT(n)
   REAL (KIND=8),DIMENSION (:),ALLOCATABLE :: work
   !!--------------------------------------------------------
   
   lwmax = max(3*min(m,n)+max(m,n), 5*min(m,n))
   !    = max(3*5 + 8, 25) = 25
   lda = max(1,m)
   lds = max(1,m,n)  
 
   ALLOCATE(work(lwmax))   
   OPEN(21,file='singular_values.dat') 
 
  ! workspace query: calculates the optimal size of the WORK array
  lwork = -1
  CALL DGESVD('N', 'N', m, n, A, lda, S, U, m, VT, n, WORK, lwork, info)
  
  !write(*,*) 'lwmax   = ', lwmax
  !write(*,*) 'work(1) = ', work(1)
  ! compute workspace size
  lwork = min(lwmax, int(work(1)))
  !write(*,*) 'lwork = ', lwork

  ! solve
  CALL DGESVD('N', 'N', m, n, A, lda, S, U, m, VT, n, WORK, lwork, info)

  write(21,*) 'info = ', info
  if (info .eq. 0) then
    write(21,*) 'Singular values were succesfully computed:'
    ! print the solution x
    do i=1, lds
      write(21,*) i, S(i)
    end do  
  else
    write(21,*) '* Error computing singula values!'
  end if
  !write(*,*)

  ! The rank of the matrix is the number of singular values that are not zero

  write(21,*) 'epsilon = ', eps
  rank = 0
  do i=1, lds
    if (abs(S(i)).gt.eps) then
      rank = rank + 1
    end if
  end do  
  !write(*,*) 'Rank of the Matrix =', rank

   !9 format('s[', i1, ']= ', f27.20)  

   !CLOSE(21)
   DEALLOCATE(work)
   
   END SUBROUTINE
!====================================================================
!! http://ww2.odu.edu/~agodunov/computing/programs/book2/Ch06/Inverse.f90
!! Computing Inverse matrix C = A^{-1} 
!! Method: Based on the Doolittle LU method
!! Input: a(n,n) - array of coefficients for matrix A 
!! Output: c(n,n) - inverse matrix of A
!! Comments: the original matrix a(n,n) will be destroyed during the calculation
   !====================================================================
   SUBROUTINE Inverse_Matrix(a,c,n)
   IMPLICIT NONE

   !INTEGER (KIND=8),PARAMETER:: n=3
   INTEGER (KIND=8) :: n
   INTEGER (KIND=8) :: i,j,k
   REAL    (KIND=8) :: coeff
   REAL    (KIND=8),DIMENSION(n,n) :: a,c,L,U
   REAL    (KIND=8),DIMENSION(n) :: b,d,x

   ! Input: matrix A
   !  a(1,1)=3. ; a(1,2)=2. ;  a(1,3)=4. 
   !  a(2,1)=2. ; a(2,2)=-3. ; a(2,3)=1.
   !  a(3,1)=1. ; a(3,2)=1. ;  a(3,3)=2.
   !====================================================================
   ! step 0: initialization for matrices L and U and b
   ! Fortran 90/95 aloows such operations on matrices
   L=0.0
   U=0.0
   b=0.0

   !! step 1: forward elimination -------------------------
   DO k=1, n-1
      DO i=k+1,n
         coeff=a(i,k)/a(k,k)
         L(i,k) = coeff
         DO j=k+1,n
            a(i,j) = a(i,j)-coeff*a(k,j)
         END DO
      END DO
   END DO

   !! Step 2: prepare L and U matrices ---------------------
   !! L matrix is a matrix of the elimination coefficient; the diagonal elements are 1.0
   DO i=1,n
     L(i,i) = 1.0
   END DO
   !! U matrix is the upper triangular part of A
   DO j=1,n
     DO i=1,j
       U(i,j) = a(i,j)
     END DO
   END DO

   !! Step 3: compute columns of the inverse matrix C--------
   DO k=1,n
      b(k)=1.0
      d(1) = b(1)
   !! Step 3a: Solve Ld=b using the forward substitution
      DO i=2,n
      d(i)=b(i)
         DO j=1,i-1
         d(i) = d(i) - L(i,j)*d(j)
      END DO
      END DO
   !! Step 3b: Solve Ux=d using the back substitution
      x(n)=d(n)/U(n,n)
      DO i = n-1,1,-1
      x(i) = d(i)
         DO j=n,i+1,-1
            x(i)=x(i)-U(i,j)*x(j)
         END DO
      x(i) = x(i)/u(i,i)
      END DO
   ! Step 3c: fill the solutions x(n) into column k of C
      DO i=1,n
         c(i,k) = x(i)
      END DO
     b(k)=0.0
   END DO

   !! Out put: ---------
   !WRITE(*,*)c(1,1),c(1,2),c(1,3)
   !WRITE(*,*)c(2,1),c(2,2),c(2,3)
   !WRITE(*,*)c(3,1),c(3,2),c(3,3) 

   END SUBROUTINE
!!========================================================================================
!!========================================================================================
   SUBROUTINE initial_CPU_time()
   IMPLICIT NONE

   INTEGER (KIND=8),DIMENSION(8) :: time

   CALL DATE_AND_TIME(values=time)     ! Get the current time

   time_2=time(2) ; time_3=time(3) ; time_5=time(5) ; time_6=time(6) ; time_7=time(7)

   END SUBROUTINE
   
!!!=======================================================================================  
   SUBROUTINE computation_time()
   IMPLICIT NONE

   INTEGER (KIND=8),DIMENSION(8) :: time1
   INTEGER (KIND=8) :: run_date,run_hour,run_minute,run_second
   REAL    (KIND=8) :: run_time

   !OPEN (24,file='time_run.dat')

   CALL DATE_AND_TIME(values=time1)     ! Get the current time
      
   run_time = (time1(2)-time_2)*1296000.+(time1(3)-time_3)*86400.&
             +(time1(5)-time_5)*3600.+(time1(6)-time_6)*60.+(time1(7)-time_7) !! second (s)
   run_date=int(run_time/86400.)
   run_hour=int(run_time/3600.-run_date*24.)
   run_minute=int(run_time/60.-run_date*1440.-run_hour*60.)
   run_second=int(run_time-run_date*86400.-run_hour*3600.-run_minute*60.)

   WRITE(24,*)run_time,run_date,run_hour,run_minute,run_second

   !WRITE(24,*)'run_date  :',run_date
   !WRITE(24,*)'run_hour  :',run_hour
   !WRITE(24,*)'run_minute:',run_minute
   !WRITE(24,*)'run_second:',run_second

   CLOSE(24)

   END SUBROUTINE computation_time      
          
   END PROGRAM
