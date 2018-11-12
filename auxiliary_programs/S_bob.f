       program SRJ_func
c   This program uses the so called S-function in order to calculate a natural
c  cubic spline interpolation through 
c   a given set of n points y(i),U(i).
c
c   The S-functions are piesewise polynomials defined through linear
c  combination of the functions A,B,C and D
c   (see the corresponding defintions). The coefficients rKL are calculated in
c  the routine Lkoef_ne(n,y,rKL). 
c   
c   The natural cubic spline is then defined as U(y)=sum_{i=1}^{n} U(i)S_i. So
c  a total of n S-functions are defined.
c
c   The S-function is defined in the double precision function
c  Scalc(x,m,n,y,rKL). Here x is the argument of the
c   S-function and m is the index of the S-function.
c
       implicit double precision (a-h,o-z)
       implicit integer (i,j,k,l,m,N)

       PARAMETER (LMAX=102,NMAX=10000)
       dimension  y(1:LMAX),rKL(1:LMAX,1:LMAX),U(1:LMAX)
       write(*,*)
       write(*,*) '--------------   S_func   ---------------'
       write(*,*) '      Last modified on 7.07.2008     '
       write(*,*)
       do i=1,LMAX
        y(i)=0.0
        U(i)=0.0
       end do

      do i=1,LMAX
       do m=1,LMAX
        rKL(m,i)=0.0
        end do
      end do

c      --------------reading general information------------
       write(*,*) 'reading s_func.inf'
       
       open(2,File='s_func.inf')
c
c      Here the input information from s_func.inf is red.
c     
c      This is the total number of points n

       read(2,*) n
c
c     And this is the values of y and U
c
       do i=1,n
        read(2,*) y(i),U(i)
       end do

       close(2)
c    Here the rKL coefficients are calculated given the grid y.
        call Lkoef(n,y,rKL)
      do i=1,n
          write(6,600) (rkl(i,m),m=1,n)
  600 format(8(1Pd10.2)/(10x,7d10.2))
          end do
        write(*,*) 'L-coef - OK!'
c    The spline function is calculated in 1000 points between y(1) and y(n)
        open(1,file='s_func.out')
        step=(y(n)-y(1))/(1000-1)
        do i=-550,1450
        x=y(1)+(i-1)*step
        Uspl=0.0
c    Here for each of the 1000 points x, the value of the spline is calculated
c  as a sum
c    oven the n values of U(m) multipled by the corresponding S_m function
         do m=1,n
          Uspl=Uspl+U(m)*Scalc(x,m,n,y,rKL)
         end do
         write(1,*) x,Uspl
        end do
        end

c***********************************************************************
      double precision function Scalc(x,m,n,y,rKL)
      INTEGER  LMAX,I,K,KK,M,N
      PARAMETER (LMAX=102)
      REAL*8  x,y1,y2,y(1:LMAX),rKL(1:LMAX,1:LMAX)
      k= 0
      kk= 0
      do i=2,n
c... select interval
          if ((x.gt.y(i-1)).and.(x.le.y(i)))  k=i
          end do
      if (x.lt.y(1)) then
          k=2
          kk=1
          end if
      if (x.gt.y(n)) then
          k=n
          kk=1
          end if
      if(x.eq.y(1)) k=2
      y1=y(k-1)
      y2=y(k)
      Scalc= 0.d0
      IF(kk.eq.0) 
     1    Scalc= rKL(m,k)*((y1-x)*(((y1-x)/(y1-y2))**2-1)/6)*(y1-y2)
     2         + rKL(m,k-1)*((x-y2)*(((x-y2)/(y1-y2))**2-1)/6)*(y1-y2)
      IF(k.EQ.m) Scalc= Scalc + (y1-x)/(y1-y2)
      IF(k-1.EQ.m) Scalc= Scalc + (x-y2)/(y1-y2)
      end
c23456789 123456789 123456789 123456789 123456789 123456789 123456789 12

c***********************************************************************
      subroutine Lkoef(n,x,A)   
c*** Based on nespl subroutine          
      INTEGER LMAX
      PARAMETER (LMAX=102)
      INTEGER I,J,N,INDX(1:LMAX)
      REAL*8 X(1:LMAX),A(1:LMAX,1:LMAX),B(1:LMAX,1:LMAX), d
c
      A(1,1)=(x(3)-x(1))/3
      A(1,2)=(x(3)-x(2))/6
      do i=2,n-3
          A(i,i-1)=(x(i+1)-x(i))/6
          A(i,i)=(x(i+2)-x(i))/3
          A(i,i+1)=(x(i+2)-x(i+1))/6
          end do
      A(n-2,n-3)=(x(n-1)-x(n-2))/6
      A(n-2,n-2)=(x(n)-x(n-2))/3  
      do i=1,n-2
          B(i,i)=1/(x(i+1)-x(i))
          B(i,i+1)=-1/(x(i+2)-x(i+1))-1/(x(i+1)-x(i))
          B(i,i+2)=1/(x(i+2)-x(i+1))
          end do  
      call ludcmp(A,n-2,LMAX,indx,d)
      do i=1,n 
          call lubksb(A,n-2,LMAX,indx,B(1,i))
          end do 
      do i=1,n-2
          do j=1,n
              A(j,i+1)=B(i,j)
              end do
          end do 
      do i=1,n
          A(i,1)=0.0
          A(i,n)=0.0
          end do
      end
c23456789 123456789 123456789 123456789 123456789 123456789 123456789 12

c***********************************************************************
      SUBROUTINE ludcmp(a,n,np,indx,d)
      INTEGER n,np,indx(n),NMAX
      double precision d,a(np,np),TINY
      PARAMETER (NMAX=500,TINY=1.0e-20)
      INTEGER i,imax,j,k
      double precision aamax,dum,sum,vv(NMAX)
      d=1.d0
      do  i=1,n
          aamax=0.
          do  j=1,n
              if (abs(a(i,j)).gt.aamax) aamax=abs(a(i,j))
              enddo
          if (aamax.eq.0.) pause 'singular matrix in ludcmp'
          vv(i)=1.d0/aamax
          enddo
      do  j=1,n
          do  i=1,j-1
              sum=a(i,j)
              do  k=1,i-1
                  sum=sum-a(i,k)*a(k,j)
                  enddo
              a(i,j)=sum
              enddo
          aamax=0.
          do  i=j,n
              sum=a(i,j)
              do  k=1,j-1
                  sum=sum-a(i,k)*a(k,j)
                  enddo
              a(i,j)=sum
              dum=vv(i)*abs(sum)
              if (dum.ge.aamax) then
                  imax=i
                  aamax=dum
                  endif
              enddo
          if(j.ne.imax)then
              do  k=1,n
                  dum=a(imax,k)
                  a(imax,k)=a(j,k)
                  a(j,k)=dum
                  enddo
              d=-d
              vv(imax)=vv(j)
              endif
          indx(j)=imax
          if(a(j,j).eq.0.)a(j,j)=TINY
              if(j.ne.n)then
                  dum=1./a(j,j)
                  do  i=j+1,n
                      a(i,j)=a(i,j)*dum
                      enddo
                  endif
          enddo
      return
      END
c23456789 123456789 123456789 123456789 123456789 123456789 123456789 12

c***********************************************************************
      SUBROUTINE lubksb(a,n,np,indx,b)
      INTEGER i,ii,j,ll, n,np,indx(n)
      double precision a(np,np),b(n), sum
      ii=0
      do  i=1,n
          ll=indx(i)
          sum=b(ll)
          b(ll)=b(i)
          if (ii.ne.0)then
              do  j=ii,i-1
                  sum=sum-a(i,j)*b(j)
                  enddo
            else if (sum.ne.0.) then
              ii=i
            endif
          b(i)=sum
          enddo
      do  i=n,1,-1
          sum=b(i)
          do  j=i+1,n
              sum=sum-a(i,j)*b(j)
              enddo
          b(i)=sum/a(i,i)
          enddo
      return
      END
c23456789 123456789 123456789 123456789 123456789 123456789 123456789 12
