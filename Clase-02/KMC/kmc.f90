 Program GNi100
 Implicit real(a-h,o-z)
 common/nn/nx,ny,nc
 common/Eact/P, Temp, Ea0, Eint
 integer step, hpic, P(4,2)
 integer, allocatable, dimension(:,:) :: C, GRID, N1
 real, allocatable, dimension(:,:) :: VELO


 open(1,file="Tar")
 open(100,file="surf.xyz")
 open(101,file="atoms.xyz")
 
 read(1,*)nx, ny, Conc
 read(1,*)Ea0, Eint
 read(1,*)nrun, hpic
 read(1,*)Temp, idum
 
 nc=nint(Conc*float(nx*ny))
 print*
 print*,"  KMC Simulation with",nc,"particles"
 
 allocate(C(nc,2))
 allocate(GRID(nx,ny))
 allocate(N1(nx,ny))
 
 allocate(VELO(nc,4))
 
 C=0; GRID=0
 
  print*,"  Setting up the surface..."
  do ic=1,nc
  100 i=int(ran2(idum)*float(nx))+1
      j=int(ran2(idum)*float(ny))+1
      if(GRID(i,j).gt.0)goto 100

	GRID(i,j)=1
	C(ic,1)=i
	C(ic,2)=j
  enddo
 
 call picSurf(nx,ny)
 call picAtom(nc,C)
 call ProcessList

 
  print*,"  Getting neighbors..."
  N1=0
  do i=1,nx
  do j=1,ny
  if(GRID(i,j).gt.0)then  
            do m=1,4          
            ix=i+P(m,1)
            iy=j+P(m,2)
            call correctIJ(ix,iy)
            N1(ix,iy)=N1(ix,iy)+1
            enddo
  endif
  enddo
  enddo


!_____________KMC loop
 print*,"  Running KMC..."
 time=0.
 DO step=1,nrun
 
  call getVel(step,i,j,GRID,C,N1,  VELO)

  vt=sum(VELO)
  vr=ran2(idum)*vt
  vsum=0.
  do ic=1,nc
    do ip=1,4
    vsum=vsum+VELO(ic,ip)
    if(vsum.gt.vr)goto 110
    enddo
  enddo
  
  if(ic.gt.nc)ic=nc
  
!_________________________ Get initial-final states
! 110 dt=-1./vt*log(ran2(idum))
110 dt=-1./VELO(ic,ip)*log(ran2(idum))
    time=time+dt
    
    i=C(ic,1)
    j=C(ic,2)
    
    in=i+P(ip,1)
    jn=j+P(ip,2)
    call correctIJ(in,jn)
    
    C(ic,1)=in
    C(ic,2)=jn 
!_________________________  Update Ngbs                      
            do m=1,4          
            ix=i+P(m,1)
            iy=j+P(m,2)
            call correctIJ(ix,iy)
            N1(ix,iy)=N1(ix,iy)-1
            enddo
            
            do m=1,4          
            ix=in+P(m,1)
            iy=jn+P(m,2)
            call correctIJ(ix,iy)
            N1(ix,iy)=N1(ix,iy)+1
            enddo            
  
!_________________________  Update positions

  GRID(i,j)=0
  GRID(in,jn)=1

  ii=mod(step,hpic)   
  if(ii.eq.0)call picAtom(nc,C)
  
  ENDDO
  
  print*
  print*,"-----------------------------------------------------"
  print*, "       Simulation time =", time, "s"
  print*,"-----------------------------------------------------"
  print*
 END
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 SUBROUTINE picSurf(nx,ny)
 implicit real(a-h,o-z)
 
 Nat=nx*ny !*2
 write(100,*)Nat
 write(100,*)
 dz=0.8
 z1=10.
 z2=z1-dz
 
 do i=1,nx

	do j=1,ny
	x=float(i)-0.5
	y=float(j)-0.5
	write(100,*)"Ni", x,y,z1
! 	x=x+0.5; y=y+0.5
! 	write(100,*)"Ni", x,y,z2
	enddo
 enddo
 End
 
 
 SUBROUTINE picAtom(nc,C)
 implicit real(a-h,o-z)
 integer C(nc,2)
 
 Nat=nx*ny*2
 write(101,*)nc
 write(101,*)
 dz=0.8
 z1=10.
 
 do i=1,nc
 write(101,*)"C", C(i,:), z1+dz
 enddo
 End
 

 ! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX 
 SUBROUTINE ProcessList
 Implicit real(a-h,o-z)
 integer P(4,2)
 real Ea(2,2,0:3,0:3), Erep(2,2), dEa(5,5,3)
 common/Eact/P, Temp, Ea0, Eint

!______Surface Diffusion
 P(1,1)= 1 ; P(1,2)= 0  !right
 P(2,1)= 0 ; P(2,2)=-1  !down
 P(3,1)=-1 ; P(3,2)= 0  !left
 P(4,1)= 0 ; P(4,2)= 1  !up
 End
 
 
 
 SUBROUTINE correctIJ(i,j)
 implicit real(a-h,o-z)
 common/nn/nx,ny,nc

 if(i.gt.nx)i=i-nx
 if(i.lt.1)i=nx+i
 if(j.gt.ny)j=j-ny
 if(j.lt.1)j=ny+j
 End 
 

 SUBROUTINE getVel(step,i,j,GRID,C,N1,  VELO)
 implicit real(a-h,o-z)
 integer P(4,2)
 integer GRID(nx,ny), C(nc,2), N1(nx,ny), step
 real VELO(nc,4)
 common/nn/nx,ny,nc
 common/Eact/P, Temp, Ea0, Eint
 
 A0=1.E13
 BkT=8.617332478E-5*Temp

 
 VELO=0.
 do ic=1,nc
    i=C(ic,1) ; j=C(ic,2)
    Nini=N1(i,j)

            ! get repulsion from initial site
            Rep0=0.
            do n=1,4
            in=i+P(n,1)
            jn=j+P(n,2)
            call correctIJ(in,jn)
               
                if(GRID(in,jn).eq.0)then
                Nfin=N1(in,jn)-1
                Ea=Ea0+Eint*float(Nfin-Nini)
                VELO(ic,n)=A0*exp(-Ea/BkT)
                endif
                
            enddo

 enddo
 END 
 
 
 
 
  FUNCTION ran2(idum)
  INTEGER idum,IM1,IM2,IMM1,IA1,IA2,IQ1,IQ2,IR1,IR2,NTAB,NDIV
  REAL ran2,AM,EPS,RNMX
!   DOUBLE PRECISION ran2,AM,EPS,RNMX
  PARAMETER (IM1=2147483563,IM2=2147483399,AM=1.d0/IM1,IMM1=IM1-1)
  PARAMETER (IA1=40014,IA2=40692,IQ1=53668,IQ2=52774,IR1=12211,IR2=3791)
  PARAMETER (NTAB=32,NDIV=1+IMM1/NTAB,EPS=1.2d-16,RNMX=1.d0-EPS)
  INTEGER idum2,j,k,iv(NTAB),iy
  SAVE iv,iy,idum2
  DATA idum2/123456789/, iv/NTAB*0/, iy/0/
  if (idum.le.0) then
    idum=max(-idum,1)
    idum2=idum
    do 11 j=NTAB+8,1,-1
      k=idum/IQ1
      idum=IA1*(idum-k*IQ1)-k*IR1
      if (idum.lt.0) idum=idum+IM1
      if (j.le.NTAB) iv(j)=idum
11      continue
    iy=iv(1)
  endif
  k=idum/IQ1
  idum=IA1*(idum-k*IQ1)-k*IR1
  if (idum.lt.0) idum=idum+IM1
  k=idum2/IQ2
  idum2=IA2*(idum2-k*IQ2)-k*IR2
  if (idum2.lt.0) idum2=idum2+IM2
  j=1+iy/NDIV
  iy=iv(j)-idum2
  iv(j)=idum
  if(iy.lt.1)iy=iy+IMM1
  ran2=min(AM*iy,RNMX)
  return
  END
  