C
C Author:  Steven C. Howell
C Purpose: Calculate interaction energy between every pair of atoms
C Created: January 2014
C
C $Id$
C
C         1         2         3         4         5         6         7
C123456789012345678901234567890123456789012345678901234567890123456789012

        subroutine wca_nbyn(coor,r0,rN,c0,cN,w,wca,n)
        integer n,i,j,r0,rN,c0,cN
        double precision coor(n,3)
        double precision wca(n,n)
        double precision cutoff,w,rij
        double precision x2,y2,z2,dx2,dy2,dz2
        double precision x1,y1,z1

cf2py intent(in) :: coor,r0,rN,c0,cN,w,wca,n
cf2py intent(in,out):: wca
cf2py intent(hide):: i,j,cutoff,rij
cf2py intent(hide):: x2,y2,z2,dx2,dy2,dz2,rij
cf2py intent(hide):: x1,y1,z1

C    to call this from python:
C
C    import sys ; sys.path.append('./')
C    import electrostatics
C
C    elenergy = collision.wca_nbyn(coor,charge,nbcutoff,nbmargin,n)
C
C
C    to setup and incorporate into python:
C
C    sudo python setup_collision.py build
C    cp build/lib*/collision.o (or .so) to directory that you are
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
               else
                  wca(i,j)=0                 
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
        double precision x1,y1,z1

cf2py intent(in) :: coor_r,coor_c,r0,rN,c0,cN,w,wca,n,m
cf2py intent(in,out):: wca
cf2py intent(hide):: i,j,cutoff,rij
cf2py intent(hide):: x2,y2,z2,dx2,dy2,dz2,rij
cf2py intent(hide):: x1,y1,z1

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
               else
                  wca(i,j)=0                 
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
        double precision x1,y1,z1

cf2py intent(in) :: coor,tbead,w,wca,n
cf2py intent(in,out):: wca
cf2py intent(hide):: i,j,cutoff,rij
cf2py intent(hide):: x2,y2,z2,dx2,dy2,dz2,rij
cf2py intent(hide):: x1,y1,z1

        cutoff = 2.**(1./6.)*w

        do 200,i=tbead,n
            x1=coor(i,1)
            y1=coor(i,2)
            z1=coor(i,3)
            do 100,j=1,i-1
               x2=coor(j,1)
               y2=coor(j,2)
               z2=coor(j,3)

               dx = x2-x1
               dy = y2-y1
               dz = z2-z1

               dx2=dx*dx
               dy2=dy*dy
               dz2=dz*dz
               rij=sqrt(dx2+dy2+dz2)

               if(rij . LT . cutoff) then
                  wca(i,j)=(w/rij)**12.-(w/rij)**6.+0.25
               else
                  wca(i,j)=0
               end if
  100   continue
  200   continue

        end
C         1         2         3         4         5         6         7
C123456789012345678901234567890123456789012345678901234567890123456789012

        subroutine overlap1(coor1,natoms1,cutoff,check)
        double precision coor1(natoms1,3)
        double precision cutoff
        integer natoms1,check,count
        double precision x1,y1,z1,x2,y2,z2,diff2,dist

cf2py intent(in) :: coor1,cutoff
cf2py intent(out):: check
cf2py intent(hide)::natoms1
cf2py intent(hide)::x1,y1,z1,x2,y2,z2,diff2,dist

        count = 1
        check = 0
        do 200,i=1,natoms1-1
            x1=coor1(i,1)
            y1=coor1(i,2)
            z1=coor1(i,3)
            do 100,j=i+1,natoms1
                x2=coor1(j,1)
                y2=coor1(j,2)
                z2=coor1(j,3)
                diff2=(x1-x2)*(x1-x2)+(y1-y2)*(y1-y2)+(z1-z2)*(z1-z2)
                dist=sqrt(diff2)
                if(dist . LT . cutoff) then
                    check=1
                    exit
                endif
            count = count + 1

  100       continue

            if(check==1) then
                exit
            endif
  200   continue

        end
C         1         2         3         4         5         6         7
C123456789012345678901234567890123456789012345678901234567890123456789012

        subroutine overlap2(coor1a,coor1b,n1a,n1b,cutoff,check)
        double precision coor1a(n1a,3),coor1b(n1b,3)
        double precision cutoff
        integer n1a,n1b,check,count
        double precision x1,y1,z1,x2,y2,z2,diff2,dist

cf2py intent(in) :: coor1a,coor1b,cutoff
cf2py intent(out):: check
cf2py intent(hide)::n1a,n1b
cf2py intent(hide)::x1,y1,z1,x2,y2,z2,diff2,dist

        count = 1
        check = 0
        do 200,i=1,n1a
            x1=coor1a(i,1)
            y1=coor1a(i,2)
            z1=coor1a(i,3)
            do 100,j=1,n1b
                x2=coor1b(j,1)
                y2=coor1b(j,2)
                z2=coor1b(j,3)
                diff2=(x1-x2)*(x1-x2)+(y1-y2)*(y1-y2)+(z1-z2)*(z1-z2)
                dist=sqrt(diff2)
                if(dist . LT . cutoff) then
C                    write (*,*) "collision between beads"
                    write (*,*) dist, i, j
C                    write (*,*) "b:", j
C                    write (*,*) "dist",dist
                    check=1
                    exit
                endif
            count = count + 1

  100       continue

            if(check==1) then
                exit
            endif
  200   continue

        end
C         1         2         3         4         5         6         7
C123456789012345678901234567890123456789012345678901234567890123456789012

        subroutine overlap_dist(coor1a,coor1b,dist,n1a,n1b)
        double precision coor1a(n1a,3),coor1b(n1b,3),dist(n1a*n1b)
        integer n1a,n1b
        double precision x1,y1,z1,x2,y2,z2,diff2

cf2py intent(in) :: coor1a,coor1b,dist
cf2py intent(in,out):: dist
cf2py intent(hide)::n1a,n1b
cf2py intent(hide)::x1,y1,z1,x2,y2,z2,diff2

        do 200,i=1,n1a
            x1=coor1a(i,1)
            y1=coor1a(i,2)
            z1=coor1a(i,3)
            do 100,j=1,n1b
                x2=coor1b(j,1)
                y2=coor1b(j,2)
                z2=coor1b(j,3)
                diff2=(x1-x2)*(x1-x2)+(y1-y2)*(y1-y2)+(z1-z2)*(z1-z2)
                dist((i-1)*n1a+j)=sqrt(diff2)
  100       continue
  200   continue

        end
C         1         2         3         4         5         6         7
C123456789012345678901234567890123456789012345678901234567890123456789012
