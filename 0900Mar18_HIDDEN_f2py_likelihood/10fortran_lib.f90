!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!! 2018.03.07: fortran lib for python
!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!==============================================================================
!!==============================================================================
    SUBROUTINE coupling_prediction(N,L,S,St1,dS,M,C_inv,W2,H02,cost2)
    IMPLICIT NONE
    !! W here = W in python code (the same order: from j--> i: W(i,j))
    INTEGER,PARAMETER ::  nloop=1000
    INTEGER :: N,L,i,iloop,stop_iteration,n_iter,i0

    REAL (KIND=8) :: H_av
    REAL (KIND=8) :: S(L,N),St1(L,N),dS(L,N),W(N),H0(nloop),cost(0:nloop),H(L)
    REAL (KIND=8) :: HS_av(N),C_inv(N,N),M(N),W1(nloop,N),W2(N,N),H02(N),cost2(N)
          
   !f2py intent(in) :: S,St1,dS,M,C_inv
   !f2py intent(hide), depend(S) :: L=shape(S,0),N=shape(S,1)
   !f2py intent(out) :: W2,H02,cost2
    
    DO i0=1,N
    H(1:L)=St1(1:L,i0) ; iloop=1 ; cost(0)=100. ; stop_iteration=0
    !!--------------------------------------------------------------   
    DO WHILE ((iloop<=nloop).and.(stop_iteration==0))
        H_av=sum(H(1:L))/L
        
        DO i=1,N
        HS_av(i)=sum((H(1:L)-H_av)*dS(1:L,i))/L       
        END DO

        W(:)=matmul(HS_av(:),C_inv(:,:))
        H0(iloop)=H_av - sum(W(:)*M(:)) 

        H(1:L)=matmul(S(1:L,:),W(:))+H0(iloop)                

        cost(iloop)=sum((St1(1:L,i0)-tanh(H(1:L)))**2)/L     
        !WRITE(*,*)iloop,cost(iloop)
                  
        IF (cost(iloop)>=cost(iloop-1)) THEN
         stop_iteration=1
        END IF

        H(:)=St1(:,i0)*H(:)/tanh(H(:))
        W1(iloop,:)=W(:)            

        iloop=iloop+1                 
    END DO
    !!--------------------------------------------------------------    
    n_iter=iloop-2

    W2(i0,:)=W1(n_iter,:)
    cost2(i0)=cost(n_iter)
    H02(i0)=H0(n_iter)

    !WRITE(*,*)i0,cost2(i0)
    END DO    

    END SUBROUTINE
!!==============================================================================
!!==============================================================================
    SUBROUTINE update_hidden_ens(N1,N2,L0,W,H0,S,S_out,nens) 
    IMPLICIT NONE
    !! W here = W.T in python code
    INTEGER :: N1,N2,t,i,L0,L01,nens
    REAL (KIND=8) :: W(N2,N2),H0(N2),S(L0,N2),S_out(L0,N1:N2)
    REAL (KIND=8) :: H11(N1:N2),P11(N1:N2),H12(N2),P1,P2,rdn_S(L0,N1:N2)

    !f2py intent(in) :: N1,W,H0,S,nens
    !f2py intent(hide), depend(S) :: L0=shape(S,0),N2=shape(W,0)
    !f2py intent(out) :: S_out

    L01=L0/nens
    CALL random_number(rdn_S)

    DO t=1,L0
        IF (mod(t,L01) > 1) THEN  !! t = 2, .., L01-1      
            H11(N1:N2)=H0(N1:N2)+matmul(S(t-1,:),W(:,N1:N2))
            P11(N1:N2)=1/(1+exp(-2*H11(N1:N2)))  
                   
            DO i=N1,N2
                S(t,i)=1.0
                H12(:)=H0(:)+matmul(S(t,:),W(:,:)) 
                P1=P11(i)/product(1+exp(-2*S(t+1,:)*H12(:)))          
                P2=(1-P11(i))/product(1+exp(-2*S(t+1,:)*(H12(:)-2.*W(i,:))))
                S(t,i) = sign(1.d+00,(P1/(P1+P2)-rdn_S(t,i)))
            END DO
        END IF
    END DO

    !DO t=1,L0
    !    IF (mod(t,L01) > 1) THEN  !! t = 2, .., L01-1      
    !        H11(N1:N2)=H0(N1:N2)+matmul(S(t-1,:),W(:,N1:N2))
    !        P11(N1:N2)=1/(1+exp(-2*H11(N1:N2)))  
                   
     !       DO i=N1,N2
     !           S(t,i)=1.0
     !           H12(:)=H0(:)+matmul(S(t,:),W(:,:)) 
     !           P1=P11(i)/product(1+exp(-2*S(t+1,:)*H12(:)))          
     !           P2=(1-P11(i))/product(1+exp(-2*S(t+1,:)*(H12(:)-2.*W(i,:))))
     !           S(t,i) = sign(1.d+00,(P1/(P1+P2)-rdn_S(t,i)))
      !      END DO
      !  ELSE

        !IF (mod(t,L01) == 0) THEN  !! t = L01 
        !    H11(N1:N2)=H0(N1:N2)+matmul(S(t-1,:),W(:,N1:N2))
        !    P11(N1:N2)=1/(1+exp(-2*H11(N1:N2)))                     
        !    S(t,N1:N2) = sign(1.d+00,(P11(N1:N2)-rdn_S(t,N1:N2)))
        !ELSE  !! t=1,L01+1,... 
        !    DO i=N1,N2
        !        S(t,i)=1.0
        !        H12(:)=H0(:)+matmul(S(t,:),W(:,:)) 
        !        P1=1/product(1+exp(-2*S(t+1,:)*H12(:)))          
        !        P2=1/product(1+exp(-2*S(t+1,:)*(H12(:)-2.*W(i,:))))
        !        S(t,i) = sign(1.d+00,(P1/(P1+P2)-rdn_S(t,i)))
        !    END DO
        !END IF
        !END IF

    !END DO 

    S_out(:,N1:N2)=S(:,N1:N2)

    END SUBROUTINE 
!!==============================================================================
   SUBROUTINE update_hidden(N1,N2,L,W,H0,S,S_out) 
   IMPLICIT NONE
   !! W here = W.T in python code
   
   INTEGER :: N1,N2,t,L,i
   REAL (KIND=8) :: W(N2,N2),H0(N2),S(L,N2),S_out(L,N1:N2)
   REAL (KIND=8) :: H11(N1:N2),P11(N1:N2),H12(N2),P1,P2,rdn_S(L,N1:N2)

   !f2py intent(in) :: N1,W,H0,S
   !f2py intent(hide), depend(S) :: L=shape(S,0),N2=shape(W,0)
   !f2py intent(out) :: S_out

   CALL random_number(rdn_S)

   DO t=2,L-1

      H11(N1:N2)=H0(N1:N2)+matmul(S(t-1,:),W(:,N1:N2))
      P11(N1:N2)=1/(1+exp(-2*H11(N1:N2)))  
               
      DO i=N1,N2
         S(t,i)=1.0
         H12(:)=H0(:)+matmul(S(t,:),W(:,:)) 
         P1=P11(i)/product(1+exp(-2*S(t+1,:)*H12(:)))          
         P2=(1-P11(i))/product(1+exp(-2*S(t+1,:)*(H12(:)-2.*W(i,:))))

         S(t,i) = sign(1.d+00,(P1/(P1+P2)-rdn_S(t,i)))

      END DO

   END DO 

   S_out(:,N1:N2)=S(:,N1:N2)
   
   END SUBROUTINE    
