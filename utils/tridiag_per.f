      SUBROUTINE tridiag_per(A, B, C, F, M)
C%%%%% A is lower diagonal %%%%%
C%%%%% B is middle diagonal %%%%
C%%%%% C is upper diagonal %%%%%
C%%%%% F is right-hand side %%%%
C%%%%% M is dimension of b %%%%%


      INTEGER M,j

      REAL A(M), B(M), C(M), F(M)
      REAL Q(M), S(M)
      REAL p,fm

      fm = F(M)

C   forward elimination sweep
      Q(1) = -C(1)/B(1)
      F(1) = F(1)/B(1)
      S(1) = -A(1)/B(1)

      do j=2,M
         p=1./(B(j)+A(j)*Q(j-1))
         Q(j)=-C(j)*p
         F(j)=(F(j)-A(j)*F(j-1))*p
         S(j)=-A(j)*S(j-1)*p
      enddo

C   backward pass
      Q(M)=0.
      S(M)=1.
 
      do j=M-1,1,-1
         S(j)=S(j)+Q(j)*S(j+1)
         Q(j)=F(j)+Q(j)*Q(j+1)
      enddo

C   final pass

      F(M) = (fm-C(M)*Q(1)-A(M)*Q(M-1))/
     $	     (C(M)*S(1)+A(M)*S(M-1)+B(M))
      do j=1,M-1
         F(j)=F(M)*S(j)+Q(j)
      enddo

      RETURN
      END

