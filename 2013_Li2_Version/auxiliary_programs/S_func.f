       program S_func
c
c
c**This program uses the so called S-function to calculate a natural 
c  cubic spline interpolation through  a given set of n points y(i),U(i).  
c  The S-functions are piesewise polynomials defined through linear 
c  combination of the functions A,B,C and D (see corresponding defintions).
c  The coefficients rKL are calculated in the routine Lkoef_ne(n,y,rKL).
c  The natural cubic spline is then defined as U(y)=sum_{i=1}^{n} U(i)S_i,
c  so a total of n S-functions are defined.
c
c** The S-function is the double precision function Scalc(x,m,n,y,rKL),
c  where x is the argument and m the index of the S-function.
c
       implicit double precision (a-h,o-z)
       implicit integer (i,j,k,l,m,N)
       PARAMETER (LMAX=102,NMAX=10000)
       dimension  y(1:LMAX),rKL(1:LMAX,1:LMAX),U(1:LMAX)

       write(*,*)
       write(*,*) '--------------   S_func   ---------------'
       write(*,*) '      Last modified on 7.07.2008     '
       
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
c  Here the input information from s_func.inf is red. 
c  This is the total number of points n
       read(2,*) n
c   
c     And this is the values of y and U
c       
      do i=1,n
          read(2,*) y(i),U(i)
          end do  
      close(2)
c
c    Here the rKL coefficients are calculated given the grid y.
c
      call Lkoef_ne(n,y,rKL)        
      write(*,*) 'L-coef - OK!'
c
c** The spline function is calculated at 1000 points between y(1) and y(n)
c       
      open(1,file='s_func.out')
       step=(y(n)-y(1))/(1000-1)
      do i=1,1000
          x=y(1)+(i-1)*step
          Uspl=0.0 
c** For each of 1000 points x, the value of the spline is calculated as a sum
c   oven the n values of U(m) multipled by the corresponding S_m function
          do m=1,n
              Uspl=Uspl+U(m)*Scalc(x,m,n,y,rKL)
              end do
          write(1,*) x,Uspl
          end do
      end
c23456789 123456789 123456789 123456789 123456789 123456789 123456789 12

c***********************************************************************
      double precision function Scalc(x,m,n,y,rKL)
      implicit double precision (a-h,o-z)
      implicit integer (i,j,k,l,N,m)
      PARAMETER (LMAX=102)
      dimension  y(1:LMAX),rKL(1:LMAX,1:LMAX)
      k=0
      do i=2,n
c... select interval
          if ((x.gt.y(i-1)).and.(x.le.y(i)))  k=i
          end do
      if (x.eq.y(1)) k=2
      if (k.eq.0) goto 110
       y1=y(k-1)
      y2=y(k)
      Scalc=ndirac(k,m)*A(x,y1,y2)+ndirac(k-1,m)*B(x,y1,y2)+
     + C(x,y1,y2)*rKL(m,k)+D(x,y1,y2)*rKL(m,k-1)        
  110       continue  	
      end
c23456789 123456789 123456789 123456789 123456789 123456789 123456789 12

c***********************************************************************
      double precision function A(z,x1,x2)
      double precision z,x1,x2
      A=(x1-z)/(x1-x2)
      end
c23456789 123456789 123456789 123456789 123456789 123456789 123456789 12

c***********************************************************************
      double precision function B(z,x1,x2)
      double precision z,x1,x2
      B=(z-x2)/(x1-x2)
      end
c23456789 123456789 123456789 123456789 123456789 123456789 123456789 12

c***********************************************************************
      double precision function C(z,x1,x2)
      double precision z,x1,x2
      C=((x1-z)*(((x1-z)/(x1-x2))**2-1)/6)*(x1-x2)
      end
c23456789 123456789 123456789 123456789 123456789 123456789 123456789 12

c***********************************************************************
      double precision function D(z,x1,x2)
      double precision z,x1,x2
      D=((z-x2)*(((z-x2)/(x1-x2))**2-1)/6)*(x1-x2)
      end
c23456789 123456789 123456789 123456789 123456789 123456789 123456789 12

c***********************************************************************
      subroutine Lkoef_ne(n,x,A)   
c*** Based on nespl subroutine          
      implicit double precision (a-h,m,p-z)
      implicit integer (i,j,n,l,k) 
      PARAMETER (LMAX=102)
      dimension X(1:LMAX),A(1:LMAX,1:LMAX),
     +          B(1:LMAX,1:LMAX),indx(1:LMAX)
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
      d=1.
      do  i=1,n
          aamax=0.
          do  j=1,n
              if (abs(a(i,j)).gt.aamax) aamax=abs(a(i,j))
              enddo
          if (aamax.eq.0.) pause 'singular matrix in ludcmp'
          vv(i)=1./aamax
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
      INTEGER n,np,indx(n)
      double precision a(np,np),b(n)
      INTEGER i,ii,j,ll
      double precision sum
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

c***********************************************************************
      integer function ndirac(n1,n2)
      integer n1,n2
      ndirac=0
      if (n1.eq.n2) then
          ndirac=1
          end if
      end
c23456789 123456789 123456789 123456789 123456789 123456789 123456789 12
