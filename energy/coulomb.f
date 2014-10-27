C
C Author:  Steven C. Howell
C Purpose: Calculate interaction energy between every pair of atoms
C Created: January 2014
C
C $Id: collision.f 72 2014-10-14 20:18:58Z schowell $
C
C         1         2         3         4         5         6         7
C123456789012345678901234567890123456789012345678901234567890123456789012
C
        subroutine coulomb(coor,charge,e,t,switchd,nbcutoff,n,el,el2)
      
        integer n,i,j,np
        double precision coor(n,3) 
        double precision charge(n)
        double precision e,t,switchd,nbcutoff
        double precision el,el2

        double precision pi,qconv,qconv2,eps,ang2m,jtokjpmol
        double precision kjtokcal,conv,conv2
        double precision xi,yi,zi,xj,yj,zj
        double precision dx,dy,dz,dx2,dy2,dz2
        double precision rij,switchscale,vij
        double precision arg1,arg2,arg3,qi,qj,kb
        
cf2py intent(in) :: coor,charge,e,t,n,switchd,nbcutoff
cf2py intent(out):: el,el2
cf2py intent(hide):: i,j,np,pi,qconv,qconv2,eps,ang2m,jtokjpmol
cf2py intent(hide):: kjtokcal,conv,conv2
cf2py intent(hide):: xi,yi,zi,xj,yj,zj
cf2py intent(hide):: dx,dy,dz,dx2,dy2,dz2
cf2py intent(hide):: rij,switchscale,vij
cf2py intent(hide):: arg1,arg2,arg3,qi,qj,kb

C eV --> C                        (yes)
        qconv=1.602177E-19        
C q*q --> C^2                     (yes) 
        qconv2=qconv*qconv        
C epsilon --> C^2/(J*m)           (yes) 
        eps=8.854187817E-12       
C angstrom --> m                  (yes) 
        ang2m=1E-10
C J --> KJ/mol                    (yes) 
        jtokjpmol=6.02214129E+20 
C KJ --> kcal                     (yes)
        kjtokcal=1d0/4.184

C kB --> (J/K):
        kb=1.3806488E-23
        
        pi=dacos(-1d0)
C        conv=(qconv2/(4.0*pi*eps*ang2m*e))*jtokjpmol*kjtokcal
C        conv=332.4683
         conv=332.06/e

C       conversion to unitless energy (t: temperature in kelvin)
        conv2=qconv2/(4.0*pi*eps*ang2m*e*kb*t)        

        np=0

C  bppbi*bppbj/rij

        el = 0d0
        el2 = 0d0
        switchscale = 1d0

        do 200,i=1,n-3
            qi=charge(i)
            xi=coor(i,1)
            yi=coor(i,2)
            zi=coor(i,3)
            do 100,j=i+3,n
               qj=charge(j)
               xj=coor(j,1)
               yj=coor(j,2)
               zj=coor(j,3)

               dx = xj-xi
               dy = yj-yi
               dz = zj-zi

               dx2=dx*dx
               dy2=dy*dy
               dz2=dz*dz
               rij=sqrt(dx2+dy2+dz2)

               if(rij . LE . switchd) then
                  switchscale=1d0
               else if(rij .GT. switchd .AND. rij .LE. nbcutoff) then 
                  arg1 = (nbcutoff-rij)*(nbcutoff-rij)
                  arg2 = (nbcutoff+2.0*rij-3.0*switchd)
                  arg3 = ((nbcutoff-switchd)**3.0)
                  switchscale = (arg1*arg2)/arg3
               else if (rij .GT. nbcutoff) then
                  switchscale=0d0
               end if

               vij=((qi*qj)/(rij))*switchscale

               el=el+vij
             
C               write(*,*) 'el = ',el
 
               np=np+1

  100   continue
  200   continue

        el2 = el*conv2
        el = el*conv

        end
        
C123456789012345678901234567890123456789012345678901234567890123456789012
C
        subroutine screen_coulomb(coor,charge,e,t,ld,switchd,nbcutoff,n,
     cel,el2)
      
        integer n,i,j,np
        double precision coor(n,3) 
        double precision charge(n)
        double precision ld,e,t,switchd,nbcutoff
        double precision el,el2

        double precision pi,qconv,qconv2,eps,ang2m,jtokjpmol
        double precision kjtokcal,conv,conv2
        double precision xi,yi,zi,xj,yj,zj
        double precision dx,dy,dz,dx2,dy2,dz2
        double precision rij,switchscale,vij
        double precision arg1,arg2,arg3,qi,qj,kb
        
cf2py intent(in) :: coor,charge,ld,e,t,n,switchd,nbcutoff
cf2py intent(out):: el,el2
cf2py intent(hide):: i,j,np,pi,qconv,qconv2,eps,ang2m,jtokjpmol
cf2py intent(hide):: kjtokcal,conv,conv2
cf2py intent(hide):: xi,yi,zi,xj,yj,zj
cf2py intent(hide):: dx,dy,dz,dx2,dy2,dz2
cf2py intent(hide):: rij,switchscale,vij
cf2py intent(hide):: arg1,arg2,arg3,qi,qj,kb

C eV --> C                        (yes)
        qconv=1.602177E-19        
C q*q --> C^2                     (yes) 
        qconv2=qconv*qconv        
C epsilon --> C^2/(J*m)           (yes) 
        eps=8.854187817E-12       
C angstrom --> m                  (yes) 
        ang2m=1E-10
C J --> KJ/mol                    (yes) 
        jtokjpmol=6.02214129E+20 
C KJ --> kcal                     (yes)
        kjtokcal=1d0/4.184

C kB --> (J/K):
        kb=1.3806488E-23
        
        pi=dacos(-1d0)
C        conv=(qconv2/(4.0*pi*eps*ang2m*e))*jtokjpmol*kjtokcal
C        conv=332.4683
         conv=332.06/e

C       conversion to unitless energy (t: temperature in kelvin)
        conv2=qconv2/(4.0*pi*eps*ang2m*e*kb*t)        

        np=0

C  bppbi*bppbj/rij

        el = 0d0
        el2 = 0d0
        switchscale = 1d0

        do 200,i=1,n-3
            qi=charge(i)
            xi=coor(i,1)
            yi=coor(i,2)
            zi=coor(i,3)
            do 100,j=i+3,n
               qj=charge(j)
               xj=coor(j,1)
               yj=coor(j,2)
               zj=coor(j,3)

               dx = xj-xi
               dy = yj-yi
               dz = zj-zi

               dx2=dx*dx
               dy2=dy*dy
               dz2=dz*dz
               rij=sqrt(dx2+dy2+dz2)

               if(rij . LE . switchd) then
                  switchscale=1d0
               else if(rij .GT. switchd .AND. rij .LE. nbcutoff) then 
                  arg1 = (nbcutoff-rij)*(nbcutoff-rij)
                  arg2 = (nbcutoff+2.0*rij-3.0*switchd)
                  arg3 = ((nbcutoff-switchd)**3.0)
                  switchscale = (arg1*arg2)/arg3
               else if (rij .GT. nbcutoff) then
                  switchscale=0d0
               end if

               vij=((qi*qj)/(rij))*exp(-rij/ld)*switchscale

               el=el+vij
             
C               write(*,*) 'el = ',el
 
               np=np+1

  100   continue
  200   continue

        el2 = el*conv2
        el = el*conv

        end
        
C ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
C ~~~~~~Original Code from J. E. Curtis~~~~~
C ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
C
Csh         subroutine inter(coor1,coor2,charge1,charge2,vdwp1,vdwp2,
Csh      cswitchd,nbcutoff,slength,er,boxl,sflag,vdwflag,natoms1,natoms2,
Csh      cnonbondenergy)
Csh       
Csh         integer natoms1,natoms2,i,j,np,sflag,vdwflag
Csh         double precision coor1(natoms1,3) 
Csh         double precision coor2(natoms2,3) 
Csh         double precision charge1(natoms1)
Csh         double precision charge2(natoms2)
Csh         double precision vdwp1(natoms1,2)
Csh         double precision vdwp2(natoms2,2)
Csh 
Csh         double precision rij,pi,qconv,qconv2,eps,qi,qj,er
Csh         double precision epsi,epsj,rmi,rmj,epsv,rminij,slength
Csh         double precision ang2m,jtokjpmol,kjtokcal,conv
Csh         double precision x2,y2,z2,dx2,dy2,dz2,switchscale,vij
Csh         double precision elenergy,vdwenergy,switchd,nbcutoff
Csh         double precision arg1,arg2,arg3,rij6,rij12,boxl,invboxl
Csh         double precision nonbondenergy
Csh 
Csh C        write(*,*) 'fortran: switchd = ',switchd
Csh C        write(*,*) 'fortran: nbcutoff = ',nbcutoff
Csh 
Csh cf2py intent(in) :: coor1,charge1,natoms1,switchd,nbcutoff,er,slength
Csh cf2py intent(in) :: coor2,charge2,natoms2,vdwp1,vdwp2,boxl,sflag
Csh cf2py intent(in) :: vdwflag
Csh cf2py intent(out):: nonbondenergy
Csh cf2py intent(hide):: i,j,pi,qconv,qconv2,eps,ang2m,jtokjpmol,kjtokcal
Csh cf2py intent(hide):: conv,np,x2,y2,z2,dx2,dy2,dz2,rij,switchscale,vij
Csh cf2py intent(hide):: arg1,arg2,arg3,rij6,rij12,q1,q2,invboxl
Csh cf2py intent(hide):: epsi,epsj,rmi,rmj,epsv,rminij
Csh 
Csh C eV --> C
Csh         qconv=1.620177E-19
Csh C q*q --> C^2
Csh         qconv2=qconv*qconv
Csh C epsilon --> C^2/(J*m)
Csh         eps=8.854187816E-12
Csh C angstrom --> m
Csh         ang2m=1E-10
Csh C J --> KJ/mol
Csh         jtokjpmol=6.022137E+20
Csh C KJ --> kcal
Csh         kjtokcal=1d0/4.184
Csh 
Csh C qconv2=1.0
Csh 
Csh         pi=dacos(-1d0)
Csh         conv=(qconv2/(4.0*pi*eps*ang2m))*jtokjpmol*kjtokcal
Csh C
Csh         conv=332.4683
Csh 
Csh         np=0
Csh 
Csh C  qi*qj/(er*rij)
Csh 
Csh         invboxl = 1d0/boxl
Csh 
Csh         elenergy = 0d0
Csh         vdwenergy = 0d0
Csh         switchscale = 1d0
Csh         nonbondingenergy = 0d0
Csh 
Csh         do 200,i=1,natoms1
Csh             qi=charge1(i)
Csh             epsi = vdwp1(i,1) 
Csh             rmi = vdwp1(i,2) 
Csh             x1=coor1(i,1)
Csh             y1=coor1(i,2)
Csh             z1=coor1(i,3)
Csh             do 100,j=1,natoms2
Csh                qj=charge2(j)
Csh                epsj = vdwp2(j,1) 
Csh                rmj = vdwp2(j,2) 
Csh                x2=coor2(j,1)
Csh                y2=coor2(j,2)
Csh                z2=coor2(j,3)
Csh 
Csh                dx = x2-x1
Csh                dy = y2-y1
Csh                dz = z2-z1
Csh 
Csh                dx = dx - boxl*anint(dx*invboxl)
Csh                dy = dy - boxl*anint(dy*invboxl)
Csh                dz = dz - boxl*anint(dz*invboxl)
Csh 
Csh                dx2=dx*dx
Csh                dy2=dy*dy
Csh                dz2=dz*dz
Csh 
Csh                rij=sqrt(dx2+dy2+dz2)
Csh 
Csh                if(rij . LT. 1E-2) then
Csh C                  write(*,*) 'rij is screwed = ',rij
Csh C                  call flush()
Csh                   exit
Csh                endif
Csh 
Csh                epsv = sqrt(epsi*epsj)
Csh                rminij = rmi+rmj
Csh 
Csh                if(rij . LE . switchd) then
Csh                   switchscale=1d0
Csh                else if(rij .GT. switchd .AND. rij .LE. nbcutoff) then 
Csh                   arg1 = (nbcutoff-rij)*(nbcutoff-rij)
Csh                   arg2 = (nbcutoff+2.0*rij-3.0*switchd)
Csh                   arg3 = ((nbcutoff-switchd)**3.0)
Csh                   switchscale = (arg1*arg2)/arg3
Csh                else if (rij .GT. nbcutoff) then
Csh                   switchscale=0d0
Csh                end if
Csh 
Csh                if(sflag .EQ. 0) then
Csh                   vij=((qi*qj)/(er*rij))*switchscale
Csh                else if(slag .EQ. 1) then
Csh                   switchscale = exp(-rij/slength)
Csh                   vij=((qi*qj)/(er*rij))*switchscale
Csh                end if
Csh 
Csh                elenergy=elenergy+vij
Csh              
Csh C               write(*,*) 'el = ',elenergy
Csh  
Csh                rij6 = (rminij/rij)**6.0
Csh                rij12 = rij6 * rij6
Csh                
Csh                if(vdwflag .EQ. 1) then
Csh                   vij=epsv*(rij12-2.0*rij6)*switchscale
Csh                else
Csh                   vij=epsv*(rij12)*switchscale
Csh                endif
Csh 
Csh                vdwenergy=vdwenergy+vij
Csh C               write(*,*) 'vdw = ',vdwenergy
Csh 
Csh                np=np+1
Csh 
Csh   100   continue
Csh   200   continue
Csh 
Csh         elenergy = elenergy*conv
Csh 
Csh         nonbondenergy = elenergy + vdwenergy
Csh 
Csh         end
