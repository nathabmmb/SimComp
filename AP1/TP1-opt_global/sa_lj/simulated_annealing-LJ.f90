!  LJ.f90 
!
!  FUNCTIONS:
!	LJ         - Entry point of console application.
!

!****************************************************************************
!
!  PROGRAM: LJ
!
!  PURPOSE:  Entry point for the console application.
!
!****************************************************************************

	program LJ

	implicit none

	! Variables
	integer i,j,k, idum, imax, nat !este �ltimo es el n�mero de �tomos en el cluster
	real*8 rx(150),ry(150),rz(150),xn(150),yn(150),zn(150),dv,f,t,t0
	real*8 xlow(150),ylow(150),zlow(150),flow !mejor configuraci�n a cada temperatura
	real*8 drmax !desplazamiento maximo
	real*8 ran1,fold,fnew
	real*8 alpha !factor de descenso de T
	integer itrial, iaccept
	integer i_T, T_steps
	real*8 ratio
	
	real*8 tmin,drmax_min
		
	open(file='salida.txt',unit=16,status='unknown')
	open(file='salida.xyz',unit=15,status='unknown')


!	parametros corrida
	t0=10.d0   !define temperatura inicial
	t=t0
	imax=10000  !cantidad de pasos de descenso de temperatura
	drmax=0.1d0 !desplazamiento maximo permitido a cada atomo
	
	tmin=1.d-4  !temperatura a la cual se detiene (criterio de convergencia)
    drmax_min=1.d-5
	
	
	alpha=0.995d0  ! T[n+1] = alpha * T[n]

	write(*,*) 'ingrese random seed y numero de atomos'
	read(*,*) idum, nat
	
	T_steps=1000*nat !pasos de MC a cada temperatura


!   generacion de configuracion inicial...

	do i=1,nat
	 rx(i)=(ran1(idum)-0.5d0)*2*dble(nat)**0.3333333333333d0
	 ry(i)=(ran1(idum)-0.5d0)*2*dble(nat)**0.3333333333333d0
	 rz(i)=(ran1(idum)-0.5d0)*2*dble(nat)**0.3333333333333d0
	enddo

!   inicializa la mejor configuracion
	xlow=rx
	ylow=ry
	zlow=rz
	flow=f(xlow,ylow,zlow,nat)


	write(16,*) 'coordenadas atomicas iniciales'
	do j=1,nat
	 write(16,'(3(f10.5,2x))') rx(j),ry(j),rz(j)
	enddo
	write(16,*)
	write(16,'(a13,f40.25)') 'energia     =',f(rx,ry,rz,nat)



	write(*,*)

!   inicializa algunos parametros	
	ratio = 1.d0
	itrial = 0
	iaccept = 0


!	comienzo de la corrida

	do i=1,imax       !pasos de descenso de temperatura
	
	    	
	do i_T=1,T_steps  !pasos de MC a cada temperatura
	
	do k=1,nat    !movida de todos los �tomos en el cluster al mismo tiempo
	
	call center(rx,ry,rz,nat) !centra el el centro de masa en el origen
	
	xn(k)=rx(k)+drmax*(2.d0*ran1(idum)-1.d0)
	yn(k)=ry(k)+drmax*(2.d0*ran1(idum)-1.d0)
	zn(k)=rz(k)+drmax*(2.d0*ran1(idum)-1.d0)

	if(abs(xn(k)).gt.dble(nat)**0.3333333333333d0) rx(k)=dble(nat)**0.3333333333333d0
	if(abs(yn(k)).gt.dble(nat)**0.3333333333333d0) ry(k)=dble(nat)**0.3333333333333d0
	if(abs(zn(k)).gt.dble(nat)**0.3333333333333d0) rz(k)=dble(nat)**0.3333333333333d0

	! si la siguiente linea esta comentada, mueve un atomo a la vez
	enddo !fin movida secuencial
	
	fnew=f(xn,yn,zn,nat)
	fold=f(rx,ry,rz,nat)
	dv=fnew-fold

	if(fnew.lt.flow) then
	 xlow=xn
	 ylow=yn
	 zlow=zn
	 flow=fnew
	endif
	
		

	if(dv.le.0.d0) then
	rx=xn
	ry=yn
	rz=zn

	iaccept=iaccept+1


	else if(exp(-dv/t).gt.ran1(idum))then
	rx=xn
	ry=yn
	rz=zn

	iaccept=iaccept+1

	
	endif
	

100	itrial = itrial +1

! si la siguiente linea esta comentada, mueve todos los atomos a la vez
 !enddo !fin movida secuencial



!   regulacion de drmax
	if(itrial.eq.500) then
	 ratio=dble(iaccept)/dble(itrial)
	 if(ratio.lt.0.5d0) then
      drmax = drmax / 1.05d0
	 else
	  drmax = drmax * 1.05d0
	 endif
	 if(drmax.gt.1.d0) drmax = 1.d0

	itrial = 0
	iaccept = 0
	endif

	enddo !a cada temperatura
	
	
!************************************************
!           para que arranque a cada temperatura a la mejor configuraci�n
!	rx=xlow
!	ry=ylow
!	rz=zlow
!************************************************

    !calcula siguiente temperatura

	t=t*alpha
	if(t.lt.tmin) then
     write(*,*)
     write(*,*) 'Temp < Tmin : convergencia alcanzada'
     write(*,*)
     goto 200
    endif

	if(drmax.lt.drmax_min) then
     write(*,*)
     write(*,*) 'drmax < drmax_min : convergencia alcanzada'
     write(*,*)
     goto 200
    endif
!***************************
!   escribir salida cada cierto numero de pasos

    if(mod(i,10).eq.0)then
    
	write(*,'(a13,i5)') 'paso        =',i
	write(*,'(a13,f15.10)') 'temperatura =',t
	write(*,'(a13,f30.25)') 'energia     =',f(rx,ry,rz,nat)
	write(*,*) 'mejor energia=',flow
	write(*,*) 'drmax =',drmax
	write(*,'(a7,f30.25)') 'ratio =',ratio
	write(*,*)



	write(16,'(a13,i5)') 'paso        =',i
	write(16,'(a13,f15.10)') 'temperatura =',t
	write(16,*) 'drmax       =',drmax  !'(a13,f18.7)'
	write(16,'(a13,f30.25)') 'energia     =',f(rx,ry,rz,nat)
	write(16,*) 'mejor energia=',flow
	write(16,*)
	write(16,*) 'coordenadas atomicas'

	do j=1,nat
	 write(16,'(3(f10.5,2x))') rx(j),ry(j),rz(j)
	enddo

	write(16,*)
	write(16,*)
	


	write(15,*) nat
	write(15,*) 'coordenadas'
	do j=1,nat
	 write(15,'(a1,3(2x,f10.5))') 'H',rx(j),ry(j),rz(j)
	enddo


	endif
!***************************

    enddo !nuevo paso de descenso de T

200	close(16)
    
    call center(rx,ry,rz,nat) !centra el el centro de masa en el origen
    
	open(file='final.xyz',unit=16,status='unknown')

	write(16,*) nat
	write(16,*) 'coordenadas'
	do j=1,nat
	 write(16,'(a1,3(2x,f10.5))') 'H',rx(j),ry(j),rz(j)
	enddo


	end program LJ

!***********************************************************

	subroutine center(rx,ry,rz,nat) !para volver a desplazar el centro de masa al origen

	real*8  xmean,ymean,zmean

	real*8  rx(150),ry(150),rz(150)
	integer i

	xmean=0.d0
	ymean=0.d0
	zmean=0.d0

	do i=1,nat
	 xmean=xmean+rx(i)
	 ymean=ymean+ry(i)
	 zmean=zmean+rz(i)
	enddo

	xmean=xmean/dble(nat)
	ymean=ymean/dble(nat)
	zmean=zmean/dble(nat)

	do i=1,nat
	 rx(i)=rx(i)-xmean
	 ry(i)=ry(i)-ymean
	 rz(i)=rz(i)-zmean
	enddo	

	end subroutine
	
!***********************************************************

	function f(x,y,z,n)

	integer i,j
	real*8 f,x(150),y(150),z(150)
	real*8 dsq
		
	f=0.d0

	do i=1,n-1
	 do j=i+1,n
	
	 dsq=((x(i)-x(j))*(x(i)-x(j))+(y(i)-y(j))*(y(i)-y(j))+(z(i)-z(j))*(z(i)-z(j)))
	 
	 f=f+(dsq**(-6.d0)-dsq**(-3.d0))

	 enddo
	enddo
	
	f=4.d0*f


	end function

!***********************************************************

!     generacion de un numero aleatorio ran1 

      real*8 function ran1(idum)
      implicit real*8(a-h,o-z)
      dimension r(97)
	integer idum
      parameter (m1=259200,ia1=7141,ic1=54773,rm1=3.8580247e-6)
      parameter (m2=134456,ia2=8121,ic2=28411,rm2=7.4373773e-6)
      parameter (m3=243000,ia3=4561,ic3=51349)
      data iff /0/
      save
      if (idum.lt.0.or.iff.eq.0) then
	iff=1
	ix1=mod(ic1-idum,m1)
	ix1=mod(ia1*ix1+ic1,m1)
	ix2=mod(ix1,m2)
	ix1=mod(ia1*ix1+ic1,m1)
	ix3=mod(ix1,m3)
	do 11 j=1,97
	  ix1=mod(ia1*ix1+ic1,m1)
	  ix2=mod(ia2*ix2+ic2,m2)
	  r(j)=(dfloat(ix1)+dfloat(ix2)*rm2)*rm1
11      continue
	idum=1
      endif
      ix1=mod(ia1*ix1+ic1,m1)
      ix2=mod(ia2*ix2+ic2,m2)
      ix3=mod(ia3*ix3+ic3,m3)
      j=1+(97*ix3)/m3
      if(j.gt.97.or.j.lt.1)pause
      ran1=r(j)
      r(j)=(dfloat(ix1)+dfloat(ix2)*rm2)*rm1
      return
      end
