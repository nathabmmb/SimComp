      program gradiente_conjugado
	
	integer, parameter :: n = 8  !dimensionalidad del problema
	real*8, dimension(n) :: p,df
	integer iter
	integer i
	real*8 ftol,fret

	ftol=1.d-10  !tolerancia como criterio de convergencia
	
	!define punto de partida para la busqueda
	!p=(/-1.d0,-1.d0, 0.d0, 0.d0/)
	
	write(*,*) 'posicion inicial: '
	read(*,*) p
      
      write(*,*) 'punto de partida:'
      do i=1,n
	write(*,*) 'x('  , i,')=',p(i)
	enddo
	

!     llama al procedimiento de gradiente conjugado
	call frprmn(p,n,ftol,iter,fret)
!     **********************************************	


      !calcula la funcion y el gradiente en el minimo encontrado
      call dfunc(n,p,df)
      
!     escribe el resultado 
	write(*,*)     
	write(*,*)'ubicacion del minimo encontrado' 
      do i=1,n
	write(*,*) 'x('  , i,')=',p(i)
	enddo
	write(*,*)
	write(*,*)'valor =' ,fret
	
	write(*,*)
	write(*,*)'grad(f) en el minimo' 
      do i=1,n
	write(*,*) 'grad(f)('  , i,')=',df(i)
	enddo
	write(*,*)
	
	write(*,*)'iteraciones =' ,iter
	write(*,*)

	end program gradiente_conjugado

!******************************************	
!***** funcion a minimizar
!******************************************	
	subroutine func(n,p,ff)
	real*8 ff
	real*8, dimension(n) :: p
	integer i

	ff=0
	do i=1,n
	ff = ff+(p(i)**4.d0-16.d0*p(i)**2+5.d0*p(i))/2.d0
	enddo
	end subroutine func
!******************************************	
!******************************************	
!******************************************	



!******************************************	
!***** gradiente de la funcion a minimizar
!******************************************	
	subroutine dfunc(n,p,df)
	real*8, dimension(n) :: p, df
	
	integer i
	do i=1,n
	df(i)=4.d0*p(i)**3-32.d0*p(i)+5.d0
	enddo

	end subroutine dfunc
!******************************************	
!******************************************	
!******************************************	






