!******************************************	
!******************************************	
!******************************************	
!     aca comienza el calculo
!******************************************	
!******************************************	
!******************************************	
      SUBROUTINE frprmn(p,n,ftol,iter,fret)
      INTEGER iter,n,NMAX,ITMAX
      DOUBLE PRECISION fret,ftol,p(n),EPS,ff !func
      EXTERNAL func
      PARAMETER (NMAX=50,ITMAX=200,EPS=1.d-10)
CU    USES dfunc,func,linmin
      INTEGER its,j
      DOUBLE PRECISION dgg,fp,gam,gg,g(NMAX),h(NMAX),xi(NMAX)
      call func(n,p,ff)
      fp=ff
      call dfunc(n,p,xi)
      do 11 j=1,n
        g(j)=-xi(j)
        h(j)=g(j)
        xi(j)=h(j)
11    continue
      do 14 its=1,ITMAX
        iter=its
        call linmin(p,xi,n,fret)
        if(2.d0*abs(fret-fp).le.ftol*(abs(fret)+abs(fp)+EPS))return
      call func(n,p,ff)
        fp=ff
        call dfunc(n,p,xi)
        gg=0.d0
        dgg=0.d0
        do 12 j=1,n
          gg=gg+g(j)**2
C         dgg=dgg+xi(j)**2
          dgg=dgg+(xi(j)+g(j))*xi(j)
12      continue
        if(gg.eq.0.d0)return
        gam=dgg/gg
        do 13 j=1,n
          g(j)=-xi(j)
          h(j)=g(j)+gam*h(j)
          xi(j)=h(j)
13      continue
14    continue
      pause 'frprmn maximum iterations exceeded'
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software 63$&.





      SUBROUTINE linmin(p,xi,n,fret)
      INTEGER n,NMAX
      DOUBLE PRECISION fret,p(n),xi(n),TOL
      real*8 dbrent
      PARAMETER (NMAX=50,TOL=1.d-4)
CU    USES brent,f1dim,mnbrak
      INTEGER j,ncom
      DOUBLE PRECISION ax,bx,fa,fb,fx,xmin,xx,pcom(NMAX),xicom(NMAX)
     *,brent
      real*8 f1dim,df1dim
      COMMON /f1com/ pcom,xicom,ncom
      EXTERNAL f1dim,df1dim
      ncom=n
      do 11 j=1,n
        pcom(j)=p(j)
        xicom(j)=xi(j)
11    continue
      ax=0.d0
      xx=1.d0
      call mnbrak(ax,xx,bx,fa,fx,fb,f1dim,n)
      fret=dbrent(ax,xx,bx,f1dim,df1dim,TOL,xmin,n)
      do 12 j=1,n
        xi(j)=xmin*xi(j)
        p(j)=p(j)+xi(j)
12    continue
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software 63$&.




      FUNCTION f1dim(x,n)
      INTEGER n,NMAX
      DOUBLE PRECISION f1dim,ff,x !func
      PARAMETER (NMAX=50)
CU    USES func
      INTEGER j,ncom
      DOUBLE PRECISION pcom(NMAX),xicom(NMAX),xt(NMAX)
      COMMON /f1com/ pcom,xicom,ncom
      do 11 j=1,ncom
        xt(j)=pcom(j)+x*xicom(j)
11    continue
      call func(n,xt,ff)
      f1dim=ff
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software 63$&.






      FUNCTION df1dim(x,n)
      INTEGER n,NMAX
      DOUBLE PRECISION df1dim,x
      PARAMETER (NMAX=50)
CU    USES dfunc
      INTEGER j,ncom
      DOUBLE PRECISION df(NMAX),pcom(NMAX),xicom(NMAX),xt(NMAX)
      COMMON /f1com/ pcom,xicom,ncom

      do 11 j=1,ncom
        xt(j)=pcom(j)+x*xicom(j)
11    continue
      call dfunc(n,xt,df)
      df1dim=0.d0
      do 12 j=1,ncom
        df1dim=df1dim+df(j)*xicom(j)
12    continue
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software 63$&.




      SUBROUTINE mnbrak(ax,bx,cx,fa,fb,fc,func,n)
      DOUBLE PRECISION ax,bx,cx,fa,fb,fc,func,GOLD,GLIMIT,TINY
      EXTERNAL func
      PARAMETER (GOLD=1.618034d0, GLIMIT=100.d0, TINY=1.d-20)
      DOUBLE PRECISION dum,fu,q,r,u,ulim
      integer n
  
        fa=func(ax,n)
      fb=func(bx,n)
        if(fb.gt.fa)then
          dum=ax
          ax=bx
          bx=dum
          dum=fb
          fb=fa
          fa=dum
        endif
        cx=bx+GOLD*(bx-ax)
      fc=func(cx,n)
1       if(fb.ge.fc)then
        r=(bx-ax)*(fb-fc)
        q=(bx-cx)*(fb-fa)
        u=bx-((bx-cx)*q-(bx-ax)*r)/(2.d0*sign(max(abs(q-r),TINY),q-r))
        ulim=bx+GLIMIT*(cx-bx)
        if((bx-u)*(u-cx).gt.0.d0)then
        fu=func(u,n)
          if(fu.lt.fc)then
            ax=bx
            fa=fb
            bx=u
            fb=fu
            return
          else if(fu.gt.fb)then
            cx=u
            fc=fu
            return
          endif
          u=cx+GOLD*(cx-bx)
        fu=func(u,n)
        else if((cx-u)*(u-ulim).gt.0.d0)then
        fu=func(u,n)
          if(fu.lt.fc)then
            bx=cx
            cx=u
            u=cx+GOLD*(cx-bx)
            fb=fc
            fc=fu
          fu=func(u,n)
          endif
        else if((u-ulim)*(ulim-cx).ge.0.d0)then
          u=ulim
        fu=func(u,n)
        else
          u=cx+GOLD*(cx-bx)
        fu=func(u,n)
        endif
        ax=bx
        bx=cx
        cx=u
        fa=fb
        fb=fc
        fc=fu
        goto 1
      endif
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software 63$&.



      FUNCTION dbrent(ax,bx,cx,f,df,tol,xmin,n)
      INTEGER ITMAX
      DOUBLE PRECISION dbrent,ax,bx,cx,tol,xmin,df,f,ZEPS
      EXTERNAL df,f
      PARAMETER (ITMAX=100,ZEPS=1.0d-10)
      INTEGER iter
      DOUBLE PRECISION a,b,d,d1,d2,du,dv,dw,dx,e,fu,fv,fw,fx,olde,tol1
     *,tol2,u,u1,u2,
     *v,w,x,xm
      LOGICAL ok1,ok2
      integer n

      a=min(ax,cx)
      b=max(ax,cx)
      v=bx
      w=v
      x=v
      e=0.d0
      fx=f(x,n)
      fv=fx
      fw=fx
      dx=df(x,n)
      dv=dx
      dw=dx
      do 11 iter=1,ITMAX
        xm=0.5d0*(a+b)
        tol1=tol*abs(x)+ZEPS
        tol2=2.d0*tol1
        if(abs(x-xm).le.(tol2-.5d0*(b-a))) goto 3
        if(abs(e).gt.tol1) then
          d1=2.d0*(b-a)
          d2=d1
          if(dw.ne.dx) d1=(w-x)*dx/(dx-dw)
          if(dv.ne.dx) d2=(v-x)*dx/(dx-dv)
          u1=x+d1
          u2=x+d2
          ok1=((a-u1)*(u1-b).gt.0.d0).and.(dx*d1.le.0.d0)
          ok2=((a-u2)*(u2-b).gt.0.d0).and.(dx*d2.le.0.d0)
          olde=e
          e=d
          if(.not.(ok1.or.ok2))then
            goto 1
          else if (ok1.and.ok2)then
            if(abs(d1).lt.abs(d2))then
              d=d1
            else
              d=d2
            endif
          else if (ok1)then
            d=d1
          else
            d=d2
          endif
          if(abs(d).gt.abs(0.5d0*olde))goto 1
          u=x+d
          if(u-a.lt.tol2 .or. b-u.lt.tol2) d=sign(tol1,xm-x)
          goto 2
        endif
1       if(dx.ge.0.d0) then
          e=a-x
        else
          e=b-x
        endif
        d=0.5d0*e
2       if(abs(d).ge.tol1) then
          u=x+d
          fu=f(u,n)
        else
          u=x+sign(tol1,d)
          fu=f(u,n)
          if(fu.gt.fx)goto 3
        endif
        du=df(u,n)
        if(fu.le.fx) then
          if(u.ge.x) then
            a=x
          else
            b=x
          endif
          v=w
          fv=fw
          dv=dw
          w=x
          fw=fx
          dw=dx
          x=u
          fx=fu
          dx=du
        else
          if(u.lt.x) then
            a=u
          else
            b=u
          endif
          if(fu.le.fw .or. w.eq.x) then
            v=w
            fv=fw
            dv=dw
            w=u
            fw=fu
            dw=du
          else if(fu.le.fv .or. v.eq.x .or. v.eq.w) then
            v=u
            fv=fu
            dv=du
          endif
        endif
11    continue
      pause 'dbrent exceeded maximum iterations'
3     xmin=x
      dbrent=fx
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software 63$&.
  