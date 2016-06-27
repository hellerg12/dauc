c    TEST FOR DIFFERENCE IN NESTED EMPIRICAL AUCs

      program dauc
      
      INTEGER  px,pz,lda,np,n0,n1,p,id(20)
      real norrv,runif,f(20),fsm(20),frsm(20),frmrcsm(20),
     1 xx(1000,1000),step(20),stepsm(20),v(20,20),steprsm(20),
     2 phi(1000,1000),norcdf,info(20,20),xpi(1000),d1(20),dd1(20,20),
     3 dd(20,20),wr(20),wi(20),tempa(20,20),fr(20),lccsq,
     4 u(1000,1000),Iinfo(9,9),sig(20,20),zi(1000),frmrc(20),
     5 d(20),s(1000,1000,9),stepr(20),ddbg(20,20),var(20,20),
     6 Idd(20,20),Iddt(20,20),sig2(20,20),xbr(1000),xxbr(1000,1000),
     7 fmrc(20),fmrcsm(20),dv(20,20),ddbb(20,20),Iddbb(20,20),
     8 wt(1000),wwt(1000,1000),ww(1000,1000,20),sig1(20,20),
     9 ddgb(20,20),ddgg(20,20),tempb(20,20),tempc(20,20),dinv(20,20),
     1 vv(20,20),wmat(1000,20),eu(1000,1000),xlcf(1000),xlcr(1000)
      character*30 infile
      external FUNCTN,FUNCTNR,FUNCTNSM,FUNCTNRSM
      COMMON  jn,ip,ipx,y(1000),w(1000,20),hndf

      write(*,*) 'sample size'
      read(*,*) jn
      jn2=int(jn/2) 
      xn=jn
      xn2=xn**2
      rxn=sqrt(xn)
      write(*,*) 'number of existing covariates'
      read(*,*) px
      write(*,*) 'number of new covariates'
      read(*,*) pz
      p=px+pz
      ip=p
      ipx=px
      write(*,*) ' Columns of covariate matrix for existing factors '
      write(*,*) ' (enter anchor variable first)  '
      do 201 j=1,px
      write(*,*) 'Factor column ', j
      read(*,*) id(j)
 201  continue
      write(*,*) ' Columns of covariate matrix for new factors '
      do 202 j=1,pz
      write(*,*) 'Factor column ', j
      read(*,*) id(px+j)
 202  continue
      
      np=20
      jb=50000  
      write(*,*) 'seed '
      read(*,*) iseed
      write(*,*) 'File where data is kept'
      read(*,*) infile
      write(*,*) 
      write(*,*)
      const=1.0/2.5066
      conv=0.00001
      iter=500     
      lda=20
      pval=0.0
      pvalsm=0.0
      
      open (UNIT=28, FILE=INFILE, STATUS='UNKNOWN')
      do 900 i=1,jn
      read(28,*) y(i),(wmat(i,j),j=1,p)
 900  continue
    
 
      do 11 i=1,jn
      do 12 j=1,p
      w(i,j) = wmat(i,id(j))
 12   continue
 11   continue     
 
      xn0=0.0
      xn1=0.0
      do 15 i=1,jn
      if(y(i).eq.1.0)  xn1=xn1+1.0
      if(y(i).eq.0.0)  xn0=xn0+1.0
 15   continue     

c#######################################################################    
c  NELDER-MEAD ALGORITHM USED TO FIND PARAMETERS THAT MAXIMIZE THE NP AUC_f
 
        do 21 k=1,p-1
        f(k)=0.0
        step(k)=0.1
 21     continue       
        nop=p-1 
        nloop=1 
        iprint=-1 
        stopcr=conv 
        iquad=0 
        simp=0 
        max=iter 
        call minim(f,step,nop,func,ymin,max,iprint,stopcr,nloop, 
     &   iquad,simp,var,functn,ifault)  
     
        fmrc(1)=1.0
        do 22 k=1,p-1
        fmrc(k+1)=f(k)  
 22     continue       


c  NELDER-MEAD ALGORITHM USED TO FIND PARAMETERS THAT MAXIMIZE THE NP AUC_r
        if(px.gt.1)  then
        do 23 k=1,px-1
        fr(k)=0.0
        stepr(k)=0.1
 23     continue       
        nop=px-1 
        nloop=1 
        iprint=-1 
        stopcr=conv 
        iquad=0 
        simp=0 
        max=iter 
        call minim(fr,step,nop,funcr,ymin,max,iprint,stopcr,nloop, 
     &   iquad,simp,var,functnr,ifault)  
     
        frmrc(1)=1.0
        do 24 k=1,px-1
        frmrc(k+1)=fr(k)  
 24     continue   
        
        else
        frmrc(1)=1.0
        endif    

c#######################################################################
c   BANDWIDTH CALCULATION
      do 31 i=1,jn
      wt(i)=0.0
      do 32 k=1,p
      wt(i)=wt(i)+(w(i,k)*fmrc(k))
 32   continue
 31   continue
 
 
      xu_1=0.0
      xu_2=0.0
      do 33 i=1,jn
      do 34 j=1,jn
      u(i,j)=wt(i)-wt(j)
      xu_1=xu_1+u(i,j)
      xu_2=xu_2+(u(i,j)**2)
 34   continue     
 33   continue
      seu=sqrt((xu_2/xn2)-((xu_1/xn2)**2))
      hnpdf=seu/(float(jn)**0.20)
      hndf=seu/(float(jn)**0.333)
c#######################################################################  
c  NELDER-MEAD ALGORITHM USED TO FIND PARAMETERS THAT MAXIMIZE THE SM AUCS_F
 
        do 41 k=1,p-1
        fsm(k)=0.0
        stepsm(k)=0.1
 41     continue       
        nop=p-1 
        nloop=1 
        iprint=-1 
        stopcr=conv 
        iquad=0 
        simp=0 
        max=iter 
        call minim(fsm,stepsm,nop,funcsm,ymin,max,iprint,stopcr,nloop, 
     &   iquad,simp,var,functnsm,ifault)    
     
        fmrcsm(1)=1.0
        do 42 k=1,p-1
        fmrcsm(k+1)=fsm(k)  
 42     continue    

 
c  NELDER-MEAD ALGORITHM USED TO FIND PARAMETERS THAT MAXIMIZE THE SM AUCS_R
        if(px.gt.1)  then
        do 43 k=1,px-1
        frsm(k)=0.0
        steprsm(k)=0.1
 43     continue       
        nop=px-1 
        nloop=1 
        iprint=-1 
        stopcr=conv 
        iquad=0 
        simp=0 
        max=iter 
        call minim(frsm,stepsm,nop,funcrsm,ymin,max,iprint,stopcr,nloop, 
     &   iquad,simp,var,functnrsm,ifault)    
     
        frmrcsm(1)=1.0
        do 44 k=1,px-1
        frmrcsm(k+1)=frsm(k)  
 44     continue     

        else
        frmrcsm(1)=1.0
        endif        
        
        write(*,*)
        write(*,*) 'Minimum AUC regression coefficients - full model'
        write(*,*)
        write(*,*) 'Coefficients for existing (X) covariates'
        do 45 k=1,px
        write(*,801)  k,fmrc(k)
 45     continue
        write(*,*) 'Coefficients for new (Z) covariates'
        do 46 k=px+1,p
        write(*,801)  k,fmrc(k)
 46     continue
        write(*,*)
        write(*,*)
        write(*,*) 'Minimum AUC regression coefficients - reduced model'
        write(*,*)
        write(*,*) 'Coefficients for existing (X) covariates'
        do 47 k=1,px
        write(*,801)  k,frmrc(k)
 47     continue
        write(*,*)
        write(*,*)
 801    format(i2,5x,f6.3)
c#######################################################################  
c  FIRST AND SECOND DERIVATIVES OF AUC     
      aucf=0.0
      aucr=0.0 
      smaucf=0.0
      smaucr=0.0
      xind=0.0
      do 51 k=1,p 
      d(k)=0.0
      do 52 l=1,p
      dd(k,l)=0.0
 52   continue
 51   continue    
  
c   PHI is NORMAL DENSITY (FOR SMOOTHING) 
      do 61 i=1,jn
      xbr(i)=0.0
      wt(i)=0.0
      do 62 k=1,px
      xbr(i)=xbr(i)+(w(i,k)*frmrcsm(k))
      wt(i)=wt(i)+(w(i,k)*fmrcsm(k))
 62   continue
      do 63 k=px+1,p
      wt(i)=wt(i)+(w(i,k)*fmrcsm(k))
 63   continue
      xlcf(i) = exp(wt(i))/(1+exp(wt(i)))
      xlcr(i) = exp(xbr(i))/(1+exp(xbr(i)))
 61   continue

      do 64 i=1,jn
      do 65 j=1,jn 
      xxbr(i,j)=xbr(i)-xbr(j)
      wwt(i,j)=wt(i)-wt(j)
      do 66 k=1,p
      ww(i,j,k)=w(i,k)-w(j,k)
 66   continue
 65   continue
 64   continue
c#######################################################################  
c  	AUC DIFFERENCE aucf-aucr
c  	SMOOTHED AUC DIFFERENCE smaucf-smaucr
c   XIND IS SUM OF INDICATORS Y_i > Y_j  [n_0*n_1]

      do 71 i=1,jn
      do 72 j=1,jn 
      phi(i,j)=const*exp(-((((wwt(i,j)/hnpdf)**2)/2)))  
      do 73 k=1,p    
      s(i,j,k) = (phi(i,j)*(ww(i,j,k)/hnpdf)) 
 73   continue
      s(i,j,1) =0.0

      if(y(i).gt.y(j))  then
      do 74 k=1,p
      d(k)=d(k)+s(i,j,k)
      do 75 l=1,p
      dd(k,l)=dd(k,l)-
     &  (phi(i,j)*(wwt(i,j)/hnpdf)*(ww(i,j,k)/hnpdf)*(ww(i,j,l)/hnpdf))
 75   continue
 74   continue     

      smaucf=smaucf+norcdf(wwt(i,j)/hndf)
      smaucr=smaucr+norcdf(xxbr(i,j)/hndf)
      if(wwt(i,j).gt.0.0)  aucf=aucf+1
      if(xxbr(i,j).gt.0.0)  aucr=aucr+1
      xind=xind+1.0    
      endif     
 72   continue
 71   continue
      
      do 70 k=1,p
      dd(1,k)=0.0
      dd(k,1)=0.0
 70   continue     

      aucf=aucf/xind
      aucr=aucr/xind
      smaucf=smaucf/xind
      smaucr=smaucr/xind   

      tstat=2.0*float(jn)*(aucf-aucr)
      tstatsm=2.0*float(jn)*(smaucf-smaucr)
      
      write(*,*) 'AUCs for the full and reduced models '
      write(*,802) aucf,aucr
 802  format(5x,f5.3,3x,f5.3)      
c#######################################################################     
c  ACCOUNT FOR ANCHORING BY FIRST VARIABLE
      do 76 k=1,p-1
      d1(k)=d(k+1)/xind
      do 77 l=1,p-1    
      dd1(k,l)=dd(k+1,l+1)/xind
      Idd(k,l)=dd1(k,l)
      if((k.le.px-1).and.(l.le.px-1))  ddbb(k,l)=dd1(k,l)
      if((k.le.px-1).and.(l.le.px-1))  Iddbb(k,l)=dd1(k,l)
      if((k.le.px-1).and.(l.gt.px-1))  ddbg(k,l-(px-1))=dd1(k,l)
      if((k.gt.px-1).and.(l.le.px-1))  ddgb(k-(px-1),l)=dd1(k,l)
      if((k.gt.px-1).and.(l.gt.px-1))  ddgg(k-(px-1),l-(px-1))=dd1(k,l)      
 77   continue
 76   continue    
      call MATINV(Idd,lda,p-1,iflag)  
      call MATINV(Iddbb,lda,px-1,iflag)
      
      do 131 k=1,pz
      do 132 l=1,px-1
      tempb(k,l)=0.0
      do 133 j=1,px-1
      tempb(k,l)=tempb(k,l)+(ddgb(k,j)*Iddbb(j,l))
 133  continue
 132  continue
 131  continue
 
      do 134 k=1,pz
      do 135 l=1,pz
      tempc(k,l)=0.0
      do 136 j=1,px-1
      tempc(k,l)=tempc(k,l)+(tempb(k,j)*ddbg(j,l))
 136  continue
 135  continue
 134  continue
 
 
      do 137 k=1,pz
      do 138 l=1,pz
      dinv(k,l)=ddgg(k,l)-tempc(k,l)
 138  continue
 137  continue 
      
c#######################################################################  
      do 81 k=1,p-1
      do 82 l=1,p-1
      sig1(k,l)=0.0
      sig2(k,l)=0.0
      do 83 i=1,jn
      do 84 j=1,jn
      do 85 ij=1,jn
      if((y(i).eq.1).and.(y(j).eq.0).and.(y(ij).eq.0).and.(ij.ne.j))then
      sig1(k,l)=sig1(k,l)+(s(i,j,k+1)*s(i,ij,l+1))
      endif
      if((y(i).eq.1).and.(y(j).eq.0).and.(y(ij).eq.1).and.(ij.ne.i))then
      sig2(k,l)=sig2(k,l)+(s(i,j,k+1)*s(ij,j,l+1))
      endif      
 85   continue
 84   continue
 83   continue 
      sig1(k,l)=sig1(k,l)/(xn1*xn0*(xn0-1))
      sig2(k,l)=sig2(k,l)/(xn1*xn0*(xn1-1))
      sig(k,l)=((xn/xn1)*sig1(k,l))+((xn/xn0)*sig2(k,l)) 
 82   continue
 81   continue  
c#######################################################################   
      do 91 k=1,p-1
      do 92 l=1,p-1
      tempa(k,l)=0.0
      Iddt(k,l)=Idd(l,k)
      do 93 j=1,p-1
      tempa(k,l)=tempa(k,l)+(Idd(k,j)*sig(j,l))
 93   continue
 92   continue
 91   continue
 
      do 94 k=1,p-1
      do 95 l=1,p-1
      v(k,l)=0.0
      do 96 j=1,p-1
      v(k,l)=v(k,l)+(tempa(k,j)*Iddt(j,l))   
 96   continue
      if((k.gt.px-1).and.(l.gt.px-1))  vv(k-(px-1),l-(px-1))=v(k,l)      
 95   continue
 94   continue
 
c######################################################################       
      do 101 k=1,p-px
      do 102 l=1,p-px
      dv(k,l)=0.0
      do 103 kl=1,p-px
      dv(k,l)=dv(k,l)+(dinv(k,kl)*vv(kl,l))
 103  continue
 102  continue
 101  continue               

c  COMPUTE EIGENVALUES (WR) OF DV MATRIX
      call elmhes(dv,p-px,np)
      call hqr(dv,p-px,np,wr,wi)
      
c   EIGENVALUES (WR) ARE NEGATIVE BECAUSE THE EXPANSION IS AROUND
c   THETA_EST RATHER THAN THETA_TRUE
      do 111 ib=1,jb
      lccsq=0.0
      do 112 k=1,p-px
      lccsq=lccsq-(wr(k)*((norrv(iseed))**2)) 
 112  continue
      if(tstat.lt.lccsq)    pval=pval+1
      if(tstatsm.lt.lccsq)  pvalsm=pvalsm+1
 111  continue     
      write(*,*)
      write(*,*)
      write(*,*)'p-value for diff in AUCs '
      if((pvalsm/float(jb)).lt.0.0001)  then
      write(*,*) '   < 0.0001 '
      else
      write(*,803) pvalsm/float(jb)
      endif
 803  format(5x,f6.4)      
c######################################################################       
c  COMPUTE A CONFIDENCE INTERVAL WHEN DELTA > 0

      if((pvalsm/float(jb)).le.0.05)  then
      eubar=0.0
      do 121 i=1,jn
      do 122 j=1,jn
      eu(i,j) = norcdf(wwt(i,j)/hndf)-norcdf(xxbr(i,j)/hndf)
      if(y(i).gt.y(j))  eubar=eubar+eu(i,j)
 122  continue
 121  continue
      eubar=eubar/xind  
      
      sigci1=0.0
      sigci2=0.0
      do 141 i=1,jn
      do 142 j=1,jn
      do 143 k=1,jn
      if((y(i).eq.1).and.(y(j).eq.0).and.(y(k).eq.0).and.(k.ne.j)) then
      sigci1=sigci1+((eu(i,j)-eubar)*(eu(i,k)-eubar))
      endif
      if((y(i).eq.1).and.(y(j).eq.0).and.(y(k).eq.1).and.(k.ne.i)) then
      sigci2=sigci2+((eu(i,j)-eubar)*(eu(k,j)-eubar))
      endif      
 143   continue  
 142   continue
 141   continue 
      sigci1=sigci1/(xn1*xn0*(xn0-1))
      sigci2=sigci2/(xn1*xn0*(xn1-1))
      vstat=(((xn/xn1)*sigci1)+((xn/xn0)*sigci2))/xn  
      sestat = sqrt(vstat)
         
      if((smaucf-smaucr).lt.1.0E-10)  then
      srdiff=1.0E-10
      else
      srdiff=smaucf-smaucr
      endif
      tsr=sqrt(srdiff)
      tsrse=sestat/(sqrt(4*(srdiff)))  
      xlcb = (tsr - (1.96*tsrse))**2
      if((tsr - (1.96*tsrse)).lt.0)  xlcb = 0         
      xucb = (tsr + (1.96*tsrse))**2
      if(xucb.gt.1.0)  xucb=1.0   
      write(*,*)
      write(*,*)   
      write(*,*) '95% back transformed confidence interval '
      write(*,*) 'for the difference in the AUC parameters'
      write(*,804) xlcb, xucb  
      endif 
 804  format(5x,f6.4,3x,f6.4)      
c######################################################################       
c 800   format(3x,3(f6.4,2x))          
c       open(21,file=outfile,status='unknown',access='append')
c#######################################################################     
      END 
c#######################################################################     
      SUBROUTINE hpsort(n,ra,ind) 
      INTEGER n,ind(n) 
      real ra(n) 
      INTEGER i,ir,j,l,ira 
      real rra 
      if (n.lt.2) return 
      l=n/2+1 
      ir=n 
 10   continue 
      if(l.gt.1)then 
       l=l-1 
       rra=ra(l) 
       ira=ind(l) 
      else 
       rra=ra(ir) 
       ira=ind(ir) 
       ra(ir)=ra(1) 
       ind(ir)=ind(1) 
       ir=ir-1 
       if(ir.eq.1)then 
         ra(1)=rra 
         ind(1)=ira 
         return 
       endif 
      endif 
      i=l 
      j=l+l 
 20    if(j.le.ir)then 
       if(j.lt.ir)then 
         if(ra(j).lt.ra(j+1))j=j+1 
       endif 
       if(rra.lt.ra(j))then 
         ra(i)=ra(j) 
         ind(i)=ind(j) 
         i=j 
         j=j+j 
       else 
         j=ir+1 
       endif 
      goto 20 
      endif 
      ra(i)=rra 
      ind(i)=ira 
      goto 10 
      END 



      real FUNCTION norcdf(x) 
      real x 
      real dx,t,P,B1,B2,B3,B4,B5,PI 
      PARAMETER (P=0.2316419,B1=0.31938153,B2=-0.356563782, 
     1 B3=1.781477937,B4=-1.821255978,B5=1.330274429,PI=3.14159265359) 
      t = 1/(1+(P*abs(x))) 
      dx = exp(-x**2/2)/sqrt(2*PI) 
      if (x.ge.0) then 
         norcdf = 1.0-dx*(B1*t+B2*t**2+B3*t**3+B4*t**4+B5*t**5) 
      else 
         norcdf = dx*(B1*t+B2*t**2+B3*t**3+B4*t**4+B5*t**5) 
      endif 
      RETURN 
      END 
 
 
      FUNCTION runif(idum)
      INTEGER idum,IM1,IM2,IMM1,IA1,IA2,IQ1,IQ2,IR1,IR2,NTAB,NDIV
      real runif,AM,EPS,RNMX
      PARAMETER (IM1=2147483563,IM2=2147483399,AM=1./IM1,IMM1=IM1-1,
     *IA1=40014,IA2=40692,IQ1=53668,IQ2=52774,IR1=12211,IR2=3791,
     *NTAB=32,NDIV=1+IMM1/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
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
 11            continue
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
      runif=min(AM*iy,RNMX)
      return
      END

      real FUNCTION norrv(idum)
      INTEGER idum
c      real norrv
CU    USES runif
      INTEGER iset
      real gset,v1,v2,runif
      SAVE iset,gset
      DATA iset/0/
      if (iset.eq.0) then
 1    v1=6.28318530718*runif(idum)
      v2=sqrt(-2.0*alog(runif(idum)))
      gset=v2*sin(v1)
      norrv=v2*cos(v1)
      iset=1
      else
      norrv=gset
      iset=0
      endif
      return
      END


       SUBROUTINE FUNCTN(f,func) 
 
       real f(20),func 
       COMMON  jn,ip,ipx,y(1000),w(1000,20),hn
       
       real  wt(1000),wwt(1000,1000),norcdf
 
      func=0.0
      xind=0.0
      
      do 101 i=1,jn
      wt(i)=w(i,1)
      do 102 k=2,ip
      wt(i)=wt(i)+(w(i,k)*f(k-1))
 102  continue
 101  continue
      
      do 111 i=1,jn
      do 112 j=1,jn  
      wwt(i,j)=wt(i)-wt(j)   
      if(y(i).gt.y(j))  then
      xind=xind+1.0
      if(wwt(i,j).gt.0)  func=func-1.0
      endif
 112  continue     
 111  continue
      func=func/xind
      
      
       RETURN
       END





       SUBROUTINE FUNCTNSM(fsm,funcsm) 
 
       real fsm(20),funcsm 
       COMMON  jn,ip,ipx,y(1000),w(1000,20),hndf
       
       real  wt(1000),wwt(1000,1000),norcdf
 
      funcsm=0.0
      xind=0.0
      
      do 101 i=1,jn
      wt(i)=w(i,1)
      do 102 k=2,ip
      wt(i)=wt(i)+(w(i,k)*fsm(k-1))
 102  continue
 101  continue
      
      do 111 i=1,jn
      do 112 j=1,jn  
      wwt(i,j)=wt(i)-wt(j)   
      if(y(i).gt.y(j))  then
      xind=xind+1.0
      funcsm=funcsm-norcdf(wwt(i,j)/hndf)
      endif
 112  continue     
 111  continue
      funcsm=funcsm/xind
      
  
       RETURN
       END


       SUBROUTINE FUNCTNR(fr,funcr) 
 
       real fr(20),funcr 
       COMMON  jn,ip,ipx,y(1000),w(1000,20),hndf
       
       real  xbr(1000),xxbr(1000,1000),norcdf
 
      funcr=0.0
      xind=0.0
      
      do 101 i=1,jn
      xbr(i)=w(i,1)
      do 102 k=2,ipx
      xbr(i)=xbr(i)+(w(i,k)*fr(k-1))
 102  continue
 101  continue
 
      do 111 i=1,jn
      do 112 j=1,jn  
      xxbr(i,j)=xbr(i)-xbr(j)    
      if(y(i).gt.y(j))  then
      xind=xind+1.0
      if(xxbr(i,j).gt.0)  funcr=funcr-1.0
      endif
 112  continue     
 111  continue
      funcr=funcr/xind
     
  
       RETURN
       END
 


       SUBROUTINE FUNCTNRSM(frsm,funcrsm) 
 
       real frsm(20),funcrsm 
       COMMON  jn,ip,ipx,y(1000),w(1000,20),hndf
       
       real  xbr(1000),xxbr(1000,1000),norcdf
 
      funcrsm=0.0
      xind=0.0
      
      do 101 i=1,jn
      xbr(i)=w(i,1)
      do 102 k=2,ipx
      xbr(i)=xbr(i)+(w(i,k)*frsm(k-1))
 102  continue
 101  continue
 
      do 111 i=1,jn
      do 112 j=1,jn  
      xxbr(i,j)=xbr(i)-xbr(j)   
      if(y(i).gt.y(j))  then
      xind=xind+1.0
      funcrsm=funcrsm-norcdf(xxbr(i,j)/hndf)
      endif
 112  continue     
 111  continue
      funcrsm=funcrsm/xind
      
  
       RETURN
       END
 
 
      subroutine minim(f,step,nop,func,ymin,max,iprint,stopcr,nloop, 
     *iquad,simp,var,functn,ifault) 
c 
c 
c   the parameter list is as follows - 
c   f         on entry, the starting values of the parameters 
c             on exit, the parameter values specifying the minimum point 
c   step      the initial step sizes 
c   nop       the number of parameters(including those to be held fixed) 
c   func      on exit, the function value at the minimum found by iteration 
c   ymin      on exit, minimum from fitted quadratic surface (jbd modification) 
c   max       the maximum number of function evaluations to be allowed 
c   iprint    parameter to control output from minim 
c             less than 0    no output 
c             equal to 0     reports of initial evidence of convergence, 
c                            final convergence ,with parameter and 
c                            function values there), overrunning of max, 
c                            and the fitting of the quadratic surface if 
c                            called for. 
c             greater than 0 as for iprint = 0, plus a progress report 
c                            on the minimisation every iprint function 
c                            evaluations. 
c   stopcr    stopping criterion 
c   nloop     convergence is tested for every nloop times the process 
c             changes the simplex. after initial convergence, nloop 
c             further changes are allowed before testing for final 
c             convergence 
c   iquad     =1  if fitting of quadratic surface required 
c             =0  if not 
c   simp      criterion for expanding simplex to overcome rounding 
c             errors before fitting quadratic surface 
c   var       on exit, contains all (jbd modification) elements of the inverse  
c             of the information matrix - dimension is nop*(nop+1)/2  
c   functn    a subroutine functn(f,func) which given an array f of 
c             parameter values returns the function value func 
c   ifault    on exit - 
c                  ifault=1 if no. of function evaluations exceeds max, 
c                  ifault=2 if information matrix is not positive 
c                           semi-definite 
c                  ifault=3 if nop<1 
c                  ifault=4 if nloop<1 
c                  ifault=0 otherwise 
c 
c 
c     f,step,nop,iprint,stopcr,nloop,iquad and simp (see ** below) 
c     must be set by the calling program 
c 
c 
c **  if iquad = 0 , simp is not used 
c 
c     f,step ( and var if iquad = 1 ) must have dimension at least nop 
c     in the calling program 
c 
c   as the program is currently set up, it will deal with up to 20 
c   parameters 
c 
      implicit real (a-h,o-z) 
      dimension f(20),step(20) 
      dimension g(21,20),h(21),pbar(20),pstar(20),pstst(20) 
      dimension aval(20),bmat(210),pmin(20),vc(210),var(210),temp(20) 
      equivalence(pmin,pstar),(aval,pbar) 
      external functn
c 
c     a is reflection coefficient, b is contraction coefficient, and 
c     c is expansion coefficient 
      data a/1.00/ 
      data b/0.50/ 
      data c/2.00/ 
      data zero/0.00/ 
      data one/1.00/ 
      data two/2.00/ 
c 
c     if progress reports desired, print heading for them 
      iw=7 
      if(iprint)5,5,55 
 55       write(iw,100)iprint 
c 100 
 100  format(22h progress report every,i4,21h function evaluations// 
     *24h eval. no.  func. value ,10x,10hparameters/)
                                                                      
c     *24h eval. no.  func. value ,10x,10hparameters/) 
c 
 5                   ifault = 0 
      if(nop.le.0) ifault=3 
      if(nloop.le.0) ifault = 4 
      if (ifault.ne.0) return 
c 
c     nap is the number of parameters to be varied,i.e. with step not 0. 
      nap=0 
      loop=0 
      iflag=0 
      do 1 i=1,nop 
      if(step(i).ne.zero) nap=nap+1 
 1         continue 
c 
c     if nap=0 evaluate function at starting point and return 
      if(nap)3,4,3 
 4         call functn(f,func) 
      return 
c 
c     set up initial simplex 
 3         do 6 i=1,nop 
 6                     g(1,i)=f(i) 
      irow=2 
      do 7 i=1,nop 
      if(step(i).eq.zero) go to 7 
      do 9 j = 1,nop 
 9               g(irow,j)=f(j) 
      g(irow,i)=g(irow,i)+step(i) 
      irow=irow+1 
 7         continue 
      np1=nap+1 
      neval=0 
      do 10 i=1,np1 
      do 11 j=1,nop 
 11             f(j)=g(i,j) 
      call functn(f,h(i)) 
      neval=neval+1 
c 
c     all points of the initial simplex are output if iprint> 0 
      if(iprint)10,10,12 
 12       write(iw,101)neval,h(i),(f(j),j=1,nop) 
 101          format(/3x,i4,3x,7g14.6/(24x,6g14.6)/) 
 10                 continue 
c 
c     now follows the basic loop, i.e. given a simplex, to determine 
c     new simplex and test for convergence as required (following the 
c     flow chart given in nelder and mead) 
c 
c     to statement 13 , determine maximum and minimum points of 
c     current simplex (function values are hmax and hmin) 
 45                        loop=loop+1 
      imax=1 
      imin=1 
      do 13 i=2,np1 
      if(h(i)-h(imax))15,15,14 
 14       imax=i 
 15            if(h(i)-h(imin))16,13,13 
 16                  imin = i 
 13                         continue 
      hmax = h(imax) 
      hmin = h(imin) 
c 
c     to statement 18 , find centroid of all vertices, excluding the 
c     maximum 
      do 17 i=1,nop 
 17             pbar(i)=zero 
      do 18 i=1,np1 
      if(i-imax)19,18,19 
 19       do 20 j=1,nop 
 20                  pbar(j) = pbar(j)+g(i,j) 
 18                           continue 
      do 602 j=1,nop 
 602           pbar(j) = pbar(j)/nap 
c 
c     reflect maximum through pbar to pstar, and evaluate function at 
c     pstar (giving hstar) 
      do 21 i=1,nop 
 21             pstar(i)=a*(pbar(i)-g(imax,i))+pbar(i) 
      call functn (pstar,hstar) 
c 
c     next 3 statements test if progress report is required at present, 
c     and if so, provide one 
c     this procedure occurs frequently in the program 
      neval=neval+1 
      if(iprint)57,57,56 
 56       if(mod(neval,iprint).eq.0) write(iw,101)neval,hstar,(pstar(j), 
     1    j=1,nop) 
 57            if(hstar-hmin)22,23,23 
c 
c     if hstar less than hmin reflect pbar through pstar to give pstst, 
c     and evaluate function there (giving hstst) 
 22                  do 24 i=1,nop 
 24                               pstst(i)=c*(pstar(i)-pbar(i))+pbar(i) 
      call functn (pstst,hstst) 
      neval=neval+1 
      if(iprint)60,60,59 
 59       if(mod(neval,iprint).eq.0) write(iw,101)neval,hstst,(pstst(j), 
     1     j=1,nop) 
 60            if(hstst-hmin)25,26,26 
c 
c     if hstst less than hmin replace maximum point of current simplex 
c     by pstst and hmax by hstar, then test (statement 26 onward) 
 25                  do 27 i=1,nop 
      if(step(i).ne.0.00) g(imax,i)=pstst(i) 
 27       continue 
      h(imax)=hstst 
      go to 41 
c 
c     if hstar not less than hmin, test if hstar greater than function 
c     value at all vertices other than maximum one 
 23       do 28 i=1,np1 
      if(i-imax)29,28,29 
 29       if(hstar-h(i))26,28,28 
 28            continue 
c 
c     if it is less than (at least) one of the vertices, replace maximum 
c     point of current simplex by pstar and hmax by hstar, then test 
c     (statement 26 onward) 
c 
c     if hstar greater than all function values excluding the maximum, 
c     test if hstar greater than hmax 
c     if not, replace maximum point by pstar and hmax by hstar for 
c     whichever simplex now in store (i.e. depending on whether hstar 
c     was greater or less than hmax), calculate the contracted point 
c     pstst and the function value there, hstst 
      if(hstar-hmax)30,30,31 
 30       do 32 i=1,nop 
      if(step(i).ne.zero) g(imax,i)=pstar(i) 
 32       continue 
      hmax=hstar 
      h(imax)=hstar 
 31       do 33 i=1,nop 
 33                  pstst(i)=b*g(imax,i)+(one-b)*pbar(i) 
      call functn(pstst,hstst) 
      neval=neval+1 
      if(iprint) 63,63,62 
 62       if(mod(neval,iprint).eq.0) write(iw,101)neval,hstst,(pstst(j), 
     1    j=1,nop) 
 63            if(hstst-hmax)35,35,34 
c 
c     if hstst less than hmax, replace maximum point by pstst and hmax 
c     by hstst, then test (statement 41 onward) 
 35                  do 36 i=1,nop 
      if(step(i).ne.zero) g(imax,i)=pstst(i) 
 36       continue 
      h(imax)=hstst 
      go to 41 
c 
c     if hstst not less than hmax, replace each point of the current 
c     simplex by a point midway between its current position and the 
c     position of the minimum point of the current simplex. evaluate 
c     function at each new vertex then test (statement 41 onward) 
 34       do 38 i=1,np1 
      if(i.eq.imin) go to 38 
      do 39 j=1,nop 
      if(step(j).ne.zero) g(i,j)=(g(i,j)+g(imin,j))/2.0 
 39       f(j) = g(i,j) 
      call functn (f,h(i)) 
      neval=neval+1 
      if(iprint)38,38,65 
c 65       if(mod(neval,iprint).eq.0)write(iw,101)neval,h(i),(f(j),j=1,nop) 
 65   if(mod(neval,iprint).eq.0)write(iw,101)neval,h(i),(f(j),j=1,nop)
 38            continue 
      go to 41 
 26       do 40 i=1,nop 
      if(step(i).ne.zero) g(imax,i)=pstar(i) 
 40       continue 
      h(imax)=hstar 
c 
c     if loop = nloop, begin tests for convergence 
c     otherwise, go back to beginning of basic loop 
 41       if(loop-nloop)45,46,45 
c 
c     test for convergence 
c     calculate mean and standard deviation of function values of 
c     current simplex 
 46            hstd=zero 
      hmean=zero 
      do 42 i=1,np1 
 42             hmean=hmean+h(i) 
      hmean=hmean/np1 
      do 601 i=1,np1 
 601           hstd = hstd+(h(i)-hmean)**2 
      hstd = sqrt(hstd/np1) 
c 
c     calculate centroid of current simplex, f, and function value 
c     there, func 
      do 53 i=1,nop 
      if(step(i).eq.zero) go to 53 
      f(i)=zero 
      do 54 j=1,np1 
 54             f(i)=f(i)+g(j,i) 
      f(i)=f(i)/np1 
 53       continue 
      call functn (f,func) 
      neval=neval+1 
      if(iprint.le.0) go to 700 
      if(mod(neval,iprint).eq.0)write(iw,101)neval,func,(f(j),j=1,nop) 
c 
c     if number of function evaluations to date has overrun the limit 
c     set (max), set ifault = 1 and return 
 700     if(neval-max) 44,44,43 
 43           if(iprint)68,67,67 
 67                 write(iw,102)max 
c 102                      
 102     format(40h number of function evaluations exceeds ,i4/)            
      write(iw,103)hstd 
 103  format(51h standard error of function values of last simplex ,g13.   
     1     6/) 
c      write(iw,104)(f(i),i=1,nop) 
c 104  format(28h centroid of last simplex  ,6g14.6/(28x,6g14.6)/)
                                                               
      write(iw,105)func 
 105     format(31h  function value at centroid   ,g13.6/) 
 68           ifault = 1 
      return 
 44       if(hstd-stopcr)72,48,48 
c 
c     if the standard deviation calculated above is not less than the 
c     criterion (stopcr), set iflag and loop to zero and begin basic 
c     loop again 
 48            iflag=0 
      loop=0 
      go to 45 
 72       if(iprint)47,70,70 
 70            write(iw,106) 
 106                format(2h *,/33h  initial evidence of convergence/) 
c      write(iw,104)(f(i),i=1,nop) 
      write(iw,105)func 
c 
c     if the standard deviation is less than the criterion test iflag 
c     iflag=0 there was no evidence of convergence on last test 
c     iflag=1 there was evidence of convergence on last test 
 47       if(iflag)49,50,49 
c 
c     if iflag=0, set iflag=1 and save mean of function values of 
c     current simplex. go to beginning of basic loop. 
 50            iflag=1 
 51                  savemn = hmean 
      loop=0 
      go to 45 
c 
c     if iflag=1, test if change in mean is less than the criterion 
c     (stopcr). if it is, process has converged 
 49       if(abs(savemn-hmean).ge.stopcr)go to 51 
      if(iprint)74,73,73 
 73       write(iw,107)neval 
 107          format(4(/),36h process converges on minimum after ,i4, 
     1         21h function evaluations//) 
      write(iw,108)(f(i),i=1,nop) 
 108     format(14h minimum at   ,6g14.6/,(14x,6g14.6)/) 
      write(iw,109)func 
 109     format(/26h minimum function value   ,g13.6/) 
      write(iw,110) 
 110     format(//16h end  of  search/1x,15(1h*)/) 
 74           continue 
 
       RETURN 
       END 



 
       SUBROUTINE MATINV (A,LDA,N,IFLAG) 
C
C-----------------------------------------------------------------------
C   MATINV   WRITTEN BY CHARLES P. REEVE, STATISTICAL ENGINEERING
C            DIVISION, NATIONAL BUREAU OF STANDARDS, GAITHERSBURG,
C            MARYLAND  20899
C
C   FOR: COMPUTING THE INVERSE OF A GENERAL N BY N MATRIX IN PLACE,
C        I.E., THE INVERSE OVERWRITES THE ORIGINAL MATRIX.  THE STEPS 
C        OF THE ALGORITHM ARE DESCRIBED BELOW AS THEY OCCUR.  ROW
C        INTERCHANGES ARE DONE AS NEEDED IN ORDER TO INCREASE THE
C        ACCURACY OF THE INVERSE MATRIX.  WITHOUT INTERCHANGES THIS
C        ALGORITHM WILL FAIL WHEN ANY OF THE LEADING PRINCIPAL
C        SUBMATRICES ARE SINGULAR OR WHEN THE MATRIX ITSELF IS
C        SINGULAR.  WITH INTERCHANGES THIS ALGORITHM WILL FAIL ONLY
C        WHEN THE MATRIX ITSELF IS SINGULAR.  THE LEADING PRINCIPAL
C
C                                   [A B C]
C        SUBMATRICES OF THE MATRIX  [D E F]  ARE  [A]  AND  [A B] .
C                                   [G H I]                 [D E]
C
C   SUBPROGRAMS CALLED: -NONE-
C
C   CURRENT VERSION COMPLETED JANUARY 15, 1987
C
C   REFERENCE: STEWART, G.W., 'INTRODUCTION TO MATRIX COMPUTATIONS',
C              ACADEMIC PRESS, INC., 1973
C-----------------------------------------------------------------------
C   DEFINITION OF PASSED PARAMETERS
C
C     * A = MATRIX (SIZE NXN) TO BE INVERTED (REAL)
C
C   * LDA = LEADING DIMENSION OF MATRIX A [LDA>=N] (INTEGER)
C
C     * N = NUMBER OF ROWS AND COLUMNS OF MATRIX A (INTEGER)
C
C   IFLAG = ERROR INDICATOR ON OUTPUT (INTEGER)   INTERPRETATION: 
C           -2 -> TOO MANY ROW INTERCHANGES NEEDED - INCREASE MX
C           -1 -> N>LDA
C            0 -> NO ERRORS DETECTED
C            K -> MATRIX A FOUND TO BE SINGULAR AT THE KTH STEP OF
C                 THE CROUT REDUCTION (1<=K<=N)
C
C   * INDICATES PARAMETERS REQUIRING INPUT VALUES 
C-----------------------------------------------------------------------
C
      implicit real (a-h,o-z)    
      PARAMETER (MX=100)
      DIMENSION A(LDA,*),IEX(MX,2)
      IFLAG = 0
C
C--- CHECK CONSISTENCY OF PASSED PARAMETERS
C
      IF (N.GT.LDA) THEN
         IFLAG = -1 
         RETURN
      ENDIF
C
C--- COMPUTE A = LU BY THE CROUT REDUCTION WHERE L IS LOWER TRIANGULAR
C--- AND U IS UNIT UPPER TRIANGULAR (ALGORITHM 3.4, P. 138 OF THE
C--- REFERENCE)
C
      NEX = 0
      DO 70 K = 1, N
         DO 20 I = K, N
            S = A(I,K)
            DO 10 L = 1, K-1
               S = S-A(I,L)*A(L,K)
 10                               CONTINUE
            A(I,K) = S
 20                      CONTINUE
C
C--- INTERCHANGE ROWS IF NECESSARY
C
         Q = 0.0
         L = 0
         DO 30 I = K, N
            R = ABS(A(I,K))
            IF (R.GT.Q) THEN
               Q = R
               L = I
            ENDIF
 30                      CONTINUE
         IF (L.EQ.0) THEN
            IFLAG = K
            RETURN
         ENDIF
         IF (L.NE.K) THEN
            NEX = NEX+1
            IF (NEX.GT.MX) THEN
               IFLAG = -2
               RETURN
            ENDIF
            IEX(NEX,1) = K
            IEX(NEX,2) = L
            DO 40 J = 1, N
               Q = A(K,J)
               A(K,J) = A(L,J)
               A(L,J) = Q
 40                               CONTINUE
         ENDIF
C
C--- END ROW INTERCHANGE SECTION
C
         DO 60 J = K+1, N
            S = A(K,J)
            DO 50 L = 1, K-1
               S = S-A(K,L)*A(L,J)
 50                               CONTINUE
            A(K,J) = S/A(K,K) 
 60                      CONTINUE
 70                                    CONTINUE
C
C--- INVERT THE LOWER TRIANGLE L IN PLACE (SIMILAR TO ALGORITHM 1.5,
C--- P. 110 OF THE REFERENCE) 
C
      DO 100 K = N, 1, -1
         A(K,K) = 1.0/A(K,K)
         DO 90 I = K-1, 1, -1 
            S = 0.0 
            DO 80 J = I+1, K
               S = S+A(J,I)*A(K,J)
 80                               CONTINUE
            A(K,I) = -S/A(I,I)
 90                      CONTINUE
 100                                  CONTINUE
C
C--- INVERT THE UPPER TRIANGLE U IN PLACE (ALGORITHM 1.5, P. 110 OF
C--- THE REFERENCE) 
C
      DO 130 K = N, 1, -1
         DO 120 I = K-1, 1, -1
            S = A(I,K)
            DO 110 J = I+1, K-1
               S = S+A(I,J)*A(J,K)
 110                             CONTINUE
            A(I,K) = -S
 120                    CONTINUE
 130                                 CONTINUE
C
C--- COMPUTE INV(A) = INV(U)*INV(L)
C
      DO 160 I = 1, N
         DO 150 J = 1, N
            IF (J.GT.I) THEN
               S = 0.0
               L = J
            ELSE
               S = A(I,J)
               L = I+1
            ENDIF
            DO 140 K = L, N
               S = S+A(I,K)*A(K,J)
 140                             CONTINUE
            A(I,J) = S
 150                    CONTINUE
 160                                 CONTINUE
C
C--- INTERCHANGE COLUMNS OF INV(A) TO REVERSE EFFECT OF ROW 
C--- INTERCHANGES OF A
C
      DO 180 I = NEX, 1, -1
         K = IEX(I,1)
         L = IEX(I,2)
         DO 170 J = 1, N
            Q = A(J,K)
            A(J,K) = A(J,L)
            A(J,L) = Q
 170                    CONTINUE
 180                                 CONTINUE
      RETURN
      END 




       	SUBROUTINE elmhes(a,n,np) 
      	INTEGER n,np
      	real a(np,np)
c  Reduction to Hessenberg form  by  the elimination method.  The  real,  nonsymmetric, n by
c  n  matrix a, stored in  an  array of  physical dimensions np  by  np, is replaced by  an  upper
c  Hessenberg matrix with  identical eigenvalues. Recommended, but  not required, is that this routine be preceded  
c  balanc. On  output, the  Hessenberg matrix is in elements a(i,j) with  i < j+1.  Elements with  i > j+1 are  to   
c  be thought of as  zero,  but are  returned with random  values.
      	INTEGER i,j,m
      	real x,y
      	do 17 m=2,n-1	
      	x=0. 
      	i=m
	      do 11 j=m,n                           
	      if(abs(a(j,m-1)).gt.abs(x))then 
	      x=a(j,m-1)
	      i=j
	      endif 
 11     continue 
	      if(i.ne.m)then                       
	      do 12 j=m-1,n 
	      y=a(i,j) 
	      a(i,j)=a(m,j) 
	      a(m,j)=y
 12     continue 
	      do 13 j=1,n 
	      y=a(j,i) 
	      a(j,i)=a(j,m) 
	      a(j,m)=y
 13     continue 
	      endif
	      if(x.ne.0.)then                    
	      do 16 i=m+1,n 
	      y=a(i,m-1) 
	      if(y.ne.0.)then
	      y=y/x
	      a(i,m-1)=y 
	      do 14 j=m,n
	      a(i,j)=a(i,j)-y*a(m,j)
 14     continue 
	      do 15 j=1,n
	      a(j,m)=a(j,m)+y*a(j,i)
 15     continue 
	      endif 
 16     continue 
	      endif 
 17     continue 
	      return 
	      END
			
			


	     SUBROUTINE hqr(a,n,np,wr,wi) 
	     INTEGER n,np
	     real a(np,np),wi(np),wr(np)
c  Finds all eigenvalues  of an n by n upper  Hessenberg matrix  a that is stored  in an np by np
c  array.  On input  a can  be exactly  as output from elmhes x11.5;  on output it  is destroyed. The  real and     
c  imaginary  parts  of the  eigenvalues  are returned in wr  and  wi, respectively.
	      INTEGER i,its,j,k,l,m,nn
	      real anorm,p,q,r,s,t,u,v,w,x,y,z
	      anorm=0.                                                                             
	      do 12   i=1,n
	      do 11   j=max(i-1,1),n 
	      anorm=anorm+abs(a(i,j))
 11     continue 
 12     continue 
	      nn=n
	      t=0	
 1	    if(nn.ge.1)then 	
	      its=0
 2	    do 13 l=nn,2,-1	
	      s=abs(a(l-1,l-1))+abs(a(l,l))
	      if((a(l,l-1)+s).eq.s)goto 3  
 13     continue 
	      l=1
 3	    x=a(nn,nn)
	      if(l.eq.nn)then	
	      wr(nn)=x+t 
	      wi(nn)=0 
	      nn=nn-1
	      else
	      y=a(nn-1,nn-1)
	      w=a(nn,nn-1)*a(nn-1,nn)
	      if(l.eq.nn-1)then	
	      p=0.5*(y-x) 
	      q=p**2+w 
	      z=sqrt(abs(q)) 
        x=x+t
	      if(q.ge.0.)then	
	      z=p+sign(z,p)
	      wr(nn)=x+z
	      wr(nn-1)=wr(nn) 
	      if(z.ne.0.)wr(nn)=x-w/z 
	      wi(nn)=0.
	      wi(nn-1)=0.
	      else	
	      wr(nn)=x+p
	      wr(nn-1)=wr(nn)
	      wi(nn)=z 
	      wi(nn-1)=-z
	      endif 
	      nn=nn-2
	      else	
c	if(its.eq.30)pause  'too  many  iterations in  hqr' 
	      if(its.eq.10.or.its.eq.20)then	
	      t=t+x
	      do 14   i=1,nn 
	      a(i,i)=a(i,i)-x
 14     continue 
	      s=abs(a(nn,nn-1))+abs(a(nn-1,nn-2))
	      x=0.75*s 
        y=x
	      w=-0.4375*s**2
	      endif 
	      its=its+1
	      do 15 m=nn-2,l,-1	
	      z=a(m,m)
	      r=x-z 
	      s=y-z 
	      p=(r*s-w)/a(m+1,m)+a(m,m+1)      
	      q=a(m+1,m+1)-z-r-s 
	      r=a(m+2,m+1)
	      s=abs(p)+abs(q)+abs(r)                  
	      p=p/s 
	      q=q/s 
	      r=r/s
	      if(m.eq.l)goto 4
	      u=abs(a(m,m-1))*(abs(q)+abs(r))
	      v=abs(p)*(abs(a(m-1,m-1))+abs(z)+abs(a(m+1,m+1)))
	      if(u+v.eq.v)goto 4                           
 15     continue 
 4	    do 16 i=m+2,nn 
 	      a(i,i-2)=0.
	      if (i.ne.m+2) a(i,i-3)=0. 
 16     continue 
	      do 19 k=m,nn-1	
	      if(k.ne.m)then 
	      p=a(k,k-1)	
	      q=a(k+1,k-1)
	      r=0.
	      if(k.ne.nn-1)r=a(k+2,k-1) 
	      x=abs(p)+abs(q)+abs(r) 
	      if(x.ne.0.)then
	      p=p/x	
	      q=q/x 
	      r=r/x
	      endif
	      endif
	      s=sign(sqrt(p**2+q**2+r**2),p)
	      if(s.ne.0.)then 
	      if(k.eq.m)then
	      if(l.ne.m)a(k,k-1)=-a(k,k-1)
	      else
	      a(k,k-1)=-s*x 
	      endif
	      p=p+s                                           
	      x=p/s 
	      y=q/s 
	      z=r/s 
	      q=q/p 
	      r=r/p
	      do 17 j=k,nn                         
	      p=a(k,j)+q*a(k+1,j)
	      if(k.ne.nn-1)then 
	      p=p+r*a(k+2,j)
	      a(k+2,j)=a(k+2,j)-p*z
	      endif
	      a(k+1,j)=a(k+1,j)-p*y 
	      a(k,j)=a(k,j)-p*x
 17     continue 
	      do 18 i=l,min(nn,k+3)          
	      p=x*a(i,k)+y*a(i,k+1)
	      if(k.ne.nn-1)then 
	      p=p+z*a(i,k+2)
	      a(i,k+2)=a(i,k+2)-p*r
	      endif
	      a(i,k+1)=a(i,k+1)-p*q 
	      a(i,k)=a(i,k)-p
 18     continue 
	      endif 
 19     continue 
	      goto 2	
	      endif 
	      endif
	      goto 1 	
	      endif 
	      return 
	      END


 
 
