	program SA

	implicit none

	! Variables
	integer i,j,k, idum, imax, ndim !este ultimo es la dimension del problema
	real*8 p(150),pn(150),dv,f,t,t0
	real*8 plow(150),flow !mejor configuraci�n a cada temperatura
	real*8 drmax !desplazamiento maximo
	real*8 ran1,fold,fnew
	real*8 alpha !factor de descenso de T
	integer itrial, iaccept
	integer i_T, T_steps
	real*8 ratio
	
	real*8 tmin,drmax_min
		
	open(file='salida.txt',unit=16,status='unknown')


!	parametros corrida
	t0=10.d0   !define temperatura inicial
	t=t0
	imax=10000  !cantidad de pasos de descenso de temperatura
	drmax=0.1d0 !desplazamiento maximo permitido a cada atomo
	
	tmin=1.d-5  !temperatura a la cual se detiene (criterio de convergencia)
        drmax_min=1.d-3
	
	
	alpha=0.995d0  ! T[n+1] = alpha * T[n]

	write(*,*) 'ingrese random seed y dimension'
	!call flush(*)
	read(*,*) idum, ndim
	
	T_steps=1000*ndim !pasos de MC a cada temperatura


!   generaci�n de configuraci�n inicial...
!	do i=1,ndim
!	 p(i)=(2.d0*ran1(idum)-1.d0)*5.d0
!	enddo

	write(*,*) 'configuracion inicial: '
	read(*,*) p(1:ndim)


!   inicializa la mejor configuracion
	plow=p
	flow=f(p,ndim)


	write(16,*) 'punto de partida'

    do i=1,ndim
	 write(16,*) 'x('  , i,')=',p(i)
	enddo
	write(16,*)
	write(16,'(a13,f40.25)') 'valor     =',f(p,ndim)


	write(16,*)

!   inicializa algunos parametros	
	ratio = 1.d0
	itrial = 0
	iaccept = 0


!	comienzo de la corrida

	do i=1,imax       !pasos de descenso de temperatura
	
	    	
	do i_T=1,T_steps  !pasos de MC a cada temperatura
	
	
!	pn=rx
!	yn=ry
!	zn=rz

	do k=1,ndim    !movida de todos los �tomos en el cluster al mismo tiempo
	

	pn(k)=p(k)+drmax*(2.d0*ran1(idum)-1.d0)

	! si la siguiente linea esta comentada, mueve una coordenada a la vez
	enddo !fin movida secuencial
	
	fnew=f(pn,ndim)
	fold=f(p,ndim)
	dv=fnew-fold

	if(fnew.lt.flow) then
	 plow=pn
	 flow=fnew
	endif
	
		

	if(dv.le.0.d0) then
	p=pn

	iaccept=iaccept+1


	else if(exp(-dv/t).gt.ran1(idum))then
	p=pn

	iaccept=iaccept+1
	
	endif
	

100	itrial = itrial +1

! si la siguiente linea esta comentada, mueve todas las coordenadas a la vez
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
!	p=plow
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
	write(*,'(a13,f30.25)') 'valor f     =',f(p,ndim)
	write(*,*) 'mejor valor=',flow
	write(*,*) 'drmax =',drmax
	write(*,'(a7,f30.25)') 'ratio =',ratio
	write(*,*)



	write(16,'(a13,i5)') 'paso        =',i
	write(16,'(a13,f15.10)') 'temperatura =',t
	write(16,*) 'drmax       =',drmax  !'(a13,f18.7)'
	write(16,'(a13,f30.25)') 'valor f     =',f(p,ndim)
	write(16,*) 'mejor valor=',flow
	write(16,*)
	write(16,*) 'valor punto'

    do j=1,ndim
	 write(16,*) 'x('  , j,')=',p(j)
	enddo
	write(16,*)

	write(16,*)
	write(16,*)
	
	endif
!***************************

    enddo !nuevo paso de descenso de T

200	close(6)
    
	open(file='final.txt',unit=16,status='unknown')

	write(*,*)
	write(16,*) ndim
	write(16,*) 'mejor punto'
    do i=1,ndim
	 write(16,*) 'x('  , i,')=',p(i)
	enddo
	write(16,*)

	write(*,*) ndim
	write(*,*) 'mejor punto'
    do i=1,ndim
	 write(*,*) 'x('  , i,')=',p(i)
	enddo


	end program SA

!***********************************************************


!***********************************************************

	function f(p,n)

	integer i
	real*8 f,p(150)
		
	f=0.d0
	do i=1,n
	f = f+(p(i)**4.d0-16.d0*p(i)**2+5.d0*p(i))/2.d0
	enddo

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
      if(j.gt.97.or.j.lt.1) read(*,*)
      ran1=r(j)
      r(j)=(dfloat(ix1)+dfloat(ix2)*rm2)*rm1
      return
      end
