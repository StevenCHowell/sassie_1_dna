C $Id: collision.f,v 1.1 2014-03-27 16:13:25 schowell Exp $
C         1         2         3         4         5         6         7
C123456789012345678901234567890123456789012345678901234567890123456789012

        subroutine wca_nbyn(coor,r0,rN,c0,cN,w,wca,n)
        integer n,i,j,r0,rN,c0,cN
        double precision coor(n,3) 
        double precision wca(n,n)
        double precision cutoff,w,rij
        double precision x2,y2,z2,dx2,dy2,dz2
        
cf2py intent(in) :: coor,r0,rN,c0,cN,w,wca,n
cf2py intent(in,out):: wca
cf2py intent(hide):: i,j,cutoff,rij
cf2py intent(hide):: x2,y2,z2,dx2,dy2,dz2,rij

C    to call this from python:
C
C    import sys ; sys.path.append('./')
C    import electrostatics
C
C    elenergy = elecrostatics.fel(coor,charge,nbcutoff,nbmargin,n)
C
C
C    to setup and incorporate into python:
C
C    sudo python setup_electrostatics.py build
C    cp build/lib*/electrostatics.o (or .so) to directory that you are
C    running in

        cutoff = 2.**(1./6.)*w

        do 200,i=r0,rN
            x1=coor(i,1)
            y1=coor(i,2)
            z1=coor(i,3)
            do 100,j=c0,cN
               x2=coor(j,1)
               y2=coor(j,2)
               z2=coor(j,3)

               dx2=(x1-x2)*(x1-x2)
               dy2=(y1-y2)*(y1-y2)
               dz2=(z1-z2)*(z1-z2)
               rij=sqrt(dx2+dy2+dz2)

               if(rij . LT . cutoff) then
                  wca(i,j)=(w/rij)**12.-(w/rij)**6.+0.25
               end if
  100   continue
  200   continue

        end
C         1         2         3         4         5         6         7
C123456789012345678901234567890123456789012345678901234567890123456789012

        subroutine wca_nbym(coor_r,coor_c,r0,rN,c0,cN,w,wca,n,m)
        integer n,m,c0,cN,r0,rN,i,j
        double precision coor_r(n,3),coor_c(m,3)
        double precision wca(n,m)
        double precision cutoff,w,rij
        double precision x2,y2,z2,dx2,dy2,dz2
        
cf2py intent(in) :: coor_r,coor_c,r0,rN,c0,cN,w,wca,n,m
cf2py intent(in,out):: wca
cf2py intent(hide):: i,j,cutoff,rij
cf2py intent(hide):: x2,y2,z2,dx2,dy2,dz2,rij

C    to call this from python:
C
C    import sys ; sys.path.append('./')
C    import electrostatics
C
C    elenergy = elecrostatics.fel(coor,charge,nbcutoff,nbmargin,natoms)
C
C
C    to setup and incorporate into python:
C
C    sudo python setup_electrostatics.py build
C    cp build/lib*/electrostatics.o (or .so) to directory that you are
C    running in

        cutoff = 2.**(1./6.)*w

        do 200,i=r0,rN
            x1=coor_r(i,1)
            y1=coor_r(i,2)
            z1=coor_r(i,3)
            do 100,j=c0,cN
               x2=coor_c(j,1)
               y2=coor_c(j,2)
               z2=coor_c(j,3)

               dx2=(x1-x2)*(x1-x2)
               dy2=(y1-y2)*(y1-y2)
               dz2=(z1-z2)*(z1-z2)
               rij=sqrt(dx2+dy2+dz2)

               if(rij . LT . cutoff) then
                  wca(i,j)=(w/rij)**12.-(w/rij)**6.+0.25
               end if
  100   continue
  200   continue

        end
C         1         2         3         4         5         6         7
C123456789012345678901234567890123456789012345678901234567890123456789012

        subroutine wca_d(coor,tbead,w,wca,n)
        integer n,tbead,i,j
        double precision coor(n,3) 
        double precision wca(n,n)
        double precision cutoff,w,rij
        double precision x2,y2,z2,dx2,dy2,dz2
        
cf2py intent(in) :: coor,tbead,w,wca,n
cf2py intent(in,out):: wca
cf2py intent(hide):: i,j,cutoff,rij
cf2py intent(hide):: x2,y2,z2,dx2,dy2,dz2,rij

C    to call this from python:
C
C    import sys ; sys.path.append('./')
C    import electrostatics
C
C    wca = elecrostatics.calcwca(coor,n,tbeads,w,wca)
C
C
C    to setup and incorporate into python:
C
C    f2py -m electrostatics -c ./electrostatics.f


        cutoff = 2.**(1./6.)*w

        do 200,i=tbead,n
            x1=coor(i,1)
            y1=coor(i,2)
            z1=coor(i,3)
            do 100,j=1,i-1
               x2=coor(j,1)
               y2=coor(j,2)
               z2=coor(j,3)

               dx2=(x1-x2)*(x1-x2)
               dy2=(y1-y2)*(y1-y2)
               dz2=(z1-z2)*(z1-z2)
               rij=sqrt(dx2+dy2+dz2)

               if(rij . LT . cutoff) then
                  wca(i,j)=(w/rij)**12.-(w/rij)**6.+0.25
               end if
  100   continue
  200   continue

        end
C         1         2         3         4         5         6         7
C123456789012345678901234567890123456789012345678901234567890123456789012