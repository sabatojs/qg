      subroutine crank_nicol_mean(q,dy,dt,nu,M)

      integer M,j
      real q(M),dy,dt,nu

      real a
      complex sub(M-1),dia(M),sup(M-1)

      complex rhs(M)

      a=nu*dt/(dy*dy)

         rhs(1)=q(1)
         dia(1)=1.
         sup(1)=0.
      do j=2,M-1
         rhs(j) = q(j)
         sub(j-1)=-a
         dia(j)=1.+2.*a
         sup(j)=-a
      enddo
         rhs(M)=q(M)
         sub(M-1)=0.
         dia(M)=1.

         call tridiag(sub,dia,sup,rhs,M)

      do j=1,M
         q(j)=real(rhs(j))
      enddo


      return
      end
