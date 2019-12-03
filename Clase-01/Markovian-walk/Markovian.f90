  Program QuantumSurfing
  implicit real*8 (a-h,o-z)
  real*8 kT
!   real, dimension(:,:), allocatable :: At
  
  open(1,file="Tar")
  open(10,file="PES.xyz")
  open(11,file="random-walk.xyz")
  open(12,file="markov-walk.xyz")
  
! ----------------------------Tarjeta-------------------------------
  read(1,*)idum0
  read(1,*)x1,y1
  read(1,*)x2,y2
  read(1,*)Temp
  kT=8.6173324E-5*Temp
  Ac=1. !E14
  pi=dacos(-1.d0)

!                             print*,Ac*exp(-0.1/kT)
!                             print*,Ac*exp(-0.2/kT)
!                             print*,Ac*exp(-0.5/kT)
  idum=idum0
  x=x1; y=y1
  dxy_i=1.          ! distance threshold for pts_i--pnt_B (initial)
  dxy_f=0.01
  it=0
  nnn=1
  

! ---------------------------- RANDOM WALK  
100 it=it+1 
                                    if(it.gt.8000)then
                                    print*,"No se encontró el punto B por Random-walk"
                                    goto 101
                                    endif
    
    ang=ran2(idum)*2.d0*pi
    dx=dxy_i*dcos(ang)
    dy=dxy_i*dsin(ang)
    
    x=x+dx
    y=y+dy
    
    call picMe(x1,y1,x2,y2,x,y)
  
    
    r=sqrt((x2-x)**2+(y2-y)**2)
    if(r.gt.1.)goto 100
    
    
    print*,"Point B reached in step", it

! ---------------------------- MC WALK    
! ---------- Make PES
101 continue
    Npts=160
    xi=-1. ; xf=2.*x2+1. ; dx=(xf-xi)/float(Npts)
    yi=-1. ; yf=2.*y2+1. ; dy=(yf-yi)/float(Npts)
    
    write(10,*)(Npts+1)**2
    write(10,*)
    
    do i=0,Npts
    x=xi+float(i)*dx
        do j=0,Npts
        y=yi+float(j)*dy
        
        f=fun(x2,y2,x,y)
        write(10,*)"H", x, y, f 
        
        enddo
    enddo
    
! ---------- Make Markovian Walk
  idum=idum0
  x=x1; y=y1
  it=0
  E0=fun(x2,y2,x1,y1)
  
200 it=it+1  
                                    if(it.gt.80000)then
                                    print*,"No se encontró el punto B por Metropoli-walk"                                    
                                    stop
                                    endif
      
      ang=ran2(idum)*2.d0*pi
      dx=dxy_i*dcos(ang)
      dy=dxy_i*dsin(ang)
!         write(8,*)it,dx,dy
    
    x=x+dx
    y=y+dy
    
    E=fun(x2,y2,x,y)
    dE=(E-E0)*0.1

    if(dE.lt.0.)then
        call picMe2(x1,y1,x2,y2,x,y,E)
        E0=E
    else
        
!         write(8,*)dE
        P=Ac*exp(-dE/kT)
        r=ran2(idum)
        
        if(P.gt.r)then
            call picMe2(x1,y1,x2,y2,x,y,E)
            E0=E
        else
            x=x-dx
            y=y-dy
            goto 200
        endif
        
    endif

    r=sqrt((x2-x)**2+(y2-y)**2)
   
    if(r.gt.dxy_i)then
        goto 200
    else
        if(nnn.eq.1)print*,"Point B reached in step", it
        nnn=2
        dxy_i=0.5*dxy_i
        if(dxy_i.le.dxy_f)stop
        goto 200
    endif
!     goto 200
    
    
!     print*,"Point B reached in step", it

  END


     
! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  Function fun(x2,y2,x,y)
  implicit real*8 (a-h,o-z)
  
  f1=-10.*exp(-0.02*(x-x2)**2-0.02*(y-y2)**2)
  
!   hole=-12.*exp(-0.2*(x-x2)**2-0.2*(y-y2)**2)
  hole=-12.*exp(-0.2*(x-x2)**2-0.2*(y-y2)**2)
  bord=0.04*(x-x2)**2+0.04*(y-y2)**2
  
!   fun=f1+bord
  fun=-f1+bord+hole
  
  END

  
  
  SUBROUTINE picMe(x1,y1,x2,y2,x,y)
  implicit real*8 (a-h,o-z)
  
  write(11,*)1; write(11,*)
!   write(11,*)"N", x1, y1, "0."
!   write(11,*)"O", x2, y2, "0."
  write(11,*)"C", x , y , "20."
  END
  
  
  
  SUBROUTINE picMe2(x1,y1,x2,y2,x,y,E)
  implicit real*8 (a-h,o-z)
  
  write(12,*)1; write(12,*)
!   write(12,*)"N", x1, y1, "0."
!   write(12,*)"O", x2, y2, "0."
  write(12,*)"C", x , y , E+1.6
  END 
  
  
  
  FUNCTION ran2(idum)
  INTEGER idum,IM1,IM2,IMM1,IA1,IA2,IQ1,IQ2,IR1,IR2,NTAB,NDIV
  DOUBLE PRECISION ran2,AM,EPS,RNMX
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
   
