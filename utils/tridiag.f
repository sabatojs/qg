      SUBROUTINE tridiag(A, B, C, F, M)
C%%%%% a is lower diagonal %%%%%
C%%%%% b is middle diagonal %%%%
C%%%%% c is upper diagonal %%%%%
C%%%%% f is right-hand side %%%%
C%%%%% M is dimension of b %%%%%


      INTEGER M,j

      COMPLEX A(M-1), B(M), C(M-1), F(M), Q(M)
      COMPLEX p

      Q(1)=-C(1)/B(1);
      F(1)=F(1)/B(1);

      do j=2,M-1
	p=1/(B(j)+A(j-1)*Q(j-1))
	Q(j)=-C(j)*p
	F(j)=(F(j)-A(j-1)*F(j-1))*p;
      end do

	p=1/(B(M)+A(M-1)*Q(M-1));
	F(M)=(F(M)-A(M-1)*F(M-1))*p;

      do j=M-1,1,-1
	F(j)=F(j)+Q(j)*F(j+1);
      end do

      RETURN
      END

