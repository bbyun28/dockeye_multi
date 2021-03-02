      program make_conformers
c
c read pdb file, generate bonds, find movable torions, generate conformers
c
      implicit none
      integer maxbnd,maxatm
      parameter ( maxbnd = 10000 )
      parameter ( maxatm = 10000 )
      integer bndorder(maxatm)
      integer bndlist(10,maxatm)
      integer nsubset,isubset(maxatm)
      integer nbnd,natm,i
      integer nring,iring(maxatm)
      integer ntor,torlist(4,maxatm)
      real*4 atxyz(3,maxatm)
      integer lunit
      character*80 crdfil
      character*14 atlab(maxatm)
c-------------------

      lunit = 10
      if(iargc().lt.1)then
        print *,'Usage: make_conformers pdbfile'
        stop
      end if
      call getarg(1,crdfil)
c      print *,'PDB file?'
c      read(5,'(a)')crdfil
c
c generate bond list
c
      call getbnd(lunit,crdfil,natm,nbnd,bndorder,bndlist,atlab,atxyz)
c
c find # of segments/molecules
c
      call findseg(natm,nbnd,bndorder,bndlist,nsubset,isubset)
c      nring = nbnd - (natm - nsubset)
c      print *,'expected # of rings: ',nring
c
c find rings
c
      call findring(natm,nbnd,bndorder,bndlist,nring,iring,atlab)
      open(11,file='make_bond_sub.lst')
      do i = 1,natm
        write(11,'(3i8)')i,isubset(i),iring(i)
      end do
      close(11)
c
c clip one segment based on distance to another segment
c
      call findtor(natm,nbnd,bndorder,bndlist,iring,atlab,ntor,torlist)
      stop
      end


      subroutine findtor(natm,nbnd,bndorder,bndlist,iring,atlab,ntor,torlist)
      implicit none
      integer maxbnd,maxatm
      parameter ( maxbnd = 10000 )
      parameter ( maxatm = 10000 )
      integer bndorder(maxatm)
      integer bndlist(10,maxatm)
      integer nbnd,natm
      integer iring(maxatm)
      integer ntor,torlist(4,maxatm)
      integer i,j,j1,k,i1,i2,i3,i4
      character*14 atlab(maxatm)
c---------------------
      ntor = 0
      do i2 = 1,natm
        if(bndorder(i2).ge.2)then
          do j = 1,bndorder(i2)
            i3 = bndlist(j,i2)
            if((i3.gt.i2).and.(bndorder(i3).ge.2))then   ! prevent double counting
c                print *,i2,i3,atlab(i2),atlab(i3)
              if((iring(i2).eq.0).or.(iring(i3).eq.0))then
                ntor = ntor + 1
                do j1 = 1,bndorder(i2)
                  k = bndlist(j1,i2)
                  if(k.ne.i3)i1=k
                end do
                do j1 = 1,bndorder(i3)
                  k = bndlist(j1,i3)
                  if(k.ne.i2)i4=k
                end do
                torlist(1,ntor) = i1
                torlist(2,ntor) = i2
                torlist(3,ntor) = i3
                torlist(4,ntor) = i4
c                write(6,'(2i6,2a14,4i6)')i2,i3,atlab(i2),atlab(i3),(torlist(k,ntor),k=1,4)
c                write(6,'(4(a14,a))')atlab(i1),'|',atlab(i2),'|',atlab(i3),'|',atlab(i4),'|'
              end if
            end if
          end do
        end if
      end do
      print *,'# of torsions found: ',ntor
      open(11,file='makebond.tor')
      write(11,'(I8)')ntor
      do i = 1,ntor
        i1 = torlist(1,i)
        i2 = torlist(2,i)
        i3 = torlist(3,i)
        i4 = torlist(4,i)
        write(11,'(4i8,3x,4(a14,a))')(torlist(k,i),k=1,4),atlab(i1),'|',atlab(i2),'|',atlab(i3),'|',atlab(i4),'|'
      end do
      close(11)
      return
      end
c-------------------------------------------------------------------
      subroutine findring(natm,nbnd,bndorder,bndlist,nring,iring,atlab)
      implicit none
      integer maxbnd,maxatm
      parameter ( maxbnd = 10000 )
      parameter ( maxatm = 10000 )
      integer bndorder(maxatm)
      integer bndlist(10,maxatm)
      integer nlevel(maxatm),ilevel
      integer nbnd,natm
      integer nring,iring(maxatm)
      integer i,j,k
      integer ndone
      integer i1,i2,i3,i4,i5,i6,i7
      character*14 atlab(maxatm)
c---------------------
      nring = 0
c recursively set atom levels in bonding tree
      ilevel = 1
      ndone = 0
      do i = 1,natm
        if(bndorder(i).eq.1)then
          nlevel(i) = ilevel
          ndone = ndone + 1
        else
          nlevel(i) = 0
        end if
      end do
      if(ndone.eq.0)then
        print *,'ERROR: every atom is in a ring!'
        stop
      end if
      do while(ndone.lt.natm)
        print *,ndone,ilevel
        do i = 1,natm
          if(nlevel(i).eq.ilevel)then
            do j = 1,bndorder(i)
              k = bndlist(j,i)
              if(nlevel(k).eq.0)then
                ndone = ndone + 1
                nlevel(k) = ilevel + 1
              end if
            end do
          end if
        end do
        ilevel = ilevel + 1
      end do
      do i = 1,natm
        print *,i,nlevel(i),atlab(i)
      end do
c      while(ndone.lt.natm)
c        ilevel = ilevel + 1
c        end do
c      end do
      end
c-------------------------------------------------------------------
      subroutine findseg(natm,nbnd,bndorder,bndlist,nsubset,isubset)
      implicit none
      integer maxbnd,maxatm
      parameter ( maxbnd = 10000 )
      parameter ( maxatm = 10000 )
      integer bndorder(maxatm)
      integer bndlist(10,maxatm)
      integer nbnd,natm
      integer nsubset,isubset(maxatm),ndone,icurr,inext
      integer nstack,istack(maxatm)
      integer idone(maxatm)
      integer nring
      integer i,j,k
c--------
      do i = 1,natm
        idone(i) = 0
      end do
      nsubset = 1
100   i = 1
      do while((idone(i).eq.1).and.(i.le.natm))
        i = i + 1
      end do
      nstack = 1
      istack(nstack) = i
      idone(i) = 1 ! change
c
c run thru points
c
      do while(nstack.gt.0)
c
c pop next point off stack
c
        icurr = istack(nstack)
        nstack = nstack - 1
c        if(idone(icurr).eq.0)then ! change
c
c if new point
c
c          idone(icurr) = 1  ! change
          isubset(icurr) = nsubset
          ndone = ndone + 1
c         print *,'doing point: ',icurr
c
c push its neighbors onto stack if they haven't been done
c
          do k = 1,bndorder(icurr)
            inext = bndlist(k,icurr)
            if(idone(inext).eq.0)then
              nstack = nstack + 1
              istack(nstack) = inext
              idone(inext) = 1 ! change
            end if
          end do
c        end if
      end do
c      print *,'# done, # points; ',ndone,natm
      if(ndone.lt.natm)then
        nsubset = nsubset + 1
        goto 100
      end if
      print *,'# of atom subsets: ',nsubset,' > one means more than 1 molecule/segment'
      if(nsubset > 1)then
        print *,'ERROR: can only deal with 1 molecule at a time'
        stop
      end if
      nring = nbnd - (natm - nsubset) 
      print *,'expected # of rings: ',nring

      return
      end 
c-------------------------------------------------------------------
      subroutine getbnd(lunit,crdfil,natm,nbnd,bndorder,bndlist,atlab,atxyz)
      implicit none
      integer maxbnd,maxatm
      parameter ( maxbnd = 10000 )
      parameter ( maxatm = 10000 )
      integer bndorder(maxatm)
      integer bndlist(10,maxatm)
      integer nbnd,natm
      integer nhyd
      integer lunit
      character*80 crdfil,line
      character*14 atlab(maxatm)
      real*4 atrad(maxatm),atxyz(3,maxatm)
      real*4 dist, radsum,slop
c------------------------------------
      integer nrec,i,j,k,n
c------------------------------------
      character atnam*6,elem*2
c------------------------------------
      data slop / 0.2 /
c------------------------------------
c
c open pdb file
c
      print *,lunit
      print *,crdfil
      open(lunit,file=crdfil,err=900,status='old')
      i = 0
      nhyd = 0
300   read(lunit,'(a)',end=400)line
        if((line(1:6).ne.'ATOM  ').and.(line(1:6).ne.'HETATM'))goto 300
c
c read header and atom names
c
        i = i + 1
        if(i.gt.maxatm)then
          print *,'exceeded maxatm: ',maxatm
          stop
        end if
        atlab(i) = line(13:26)
        atnam = line(12:17)
        call up(atnam,6)
        call elb(atnam,6)
        read(line(31:54),'(3f8.3)')(atxyz(k,i),k=1,3)
c
c assign covalent radius
c
c R.T. Sanderson in Chemical Periodicity, Reinhold, New York, USA, 1962.
c L.E. Sutton (ed.) in Table of interatomic distances and configuration in 
c molecules and ions, Supplement 1956-1959, Special publication No. 18, Chemical Society, London, UK, 1965.
c J.E. Huheey, E.A. Keiter, and R.L. Keiter in Inorganic Chemistry : 
c Principles of Structure and Reactivity, 4th edition, HarperCollins, New York, USA, 1993.
c W.W. Porterfield in Inorganic chemistry, a unified approach, Addison 
c Wesley Publishing Co., Reading Massachusetts, USA, 1984.
c A.M. James and M.P. Lord in Macmillan's Chemical and Physical Data, Macmillan, London, UK, 1992.
c
c to eliminate H atoms, set radius to large -ve number instead of 0.32
        if(atnam(1:1).eq.'H')then
          nhyd = nhyd + 1
          atrad(i) =  0.32
        else if(atnam(1:2).eq.'1H')then
          nhyd = nhyd + 1
          atrad(i) =  0.32
        else if(atnam(1:2).eq.'2H')then
          nhyd = nhyd + 1
          atrad(i) =  0.32
        else if(atnam(1:2).eq.'3H')then
          nhyd = nhyd + 1
          atrad(i) =  0.32
        else if(atnam(1:3).eq.'CAL')then
          atrad(i) = 1.74
        else if(atnam(1:2).eq.'Cl')then
          atrad(i) = 0.99
        else if(atnam(1:2).eq.'Na')then
          atrad(i) = 1.54
        else if(atnam(1:1).eq.'C')then
          atrad(i) = 0.77
c          print *,'carbon'
        else if(atnam(1:1).eq.'N')then
          atrad(i) = 0.75
        else if(atnam(1:1).eq.'O')then
          atrad(i) = 0.73
        else if(atnam(1:1).eq.'P')then
          atrad(i) = 1.06
        else if(atnam(1:1).eq.'S')then
          atrad(i) = 1.02
        else
          print *,'WARNING I could not find radius for ',atnam
          atrad(i) = 0.
        end if
c        print *,atnam
      goto 300
400   continue
c
c close PDB file
c
      close(lunit)
      natm = i
      print *,'PDB atom records read: ',natm
      print *,'number of hydrogens: ',nhyd
      nbnd = 0
c      print *,'Atom radii: '
      do i = 1,natm
        bndorder(i) = 0
c        write(6,*)i,atrad(i)
      end do
      do i = 1,natm
        do j = i+1,natm
          dist = 0.
          do k = 1,3
            dist = dist + (atxyz(k,i) - atxyz(k,j))**2
          end do
          dist = sqrt(dist) - (atrad(i) + atrad(j) + slop)
          if(dist.lt.0.)then
            nbnd = nbnd + 1
c            print *,'i,j: ',i,j
            bndorder(i) = bndorder(i) + 1
            bndorder(j) = bndorder(j) + 1
            bndlist(bndorder(i),i) = j
            bndlist(bndorder(j),j) = i
          end if
        end do
      end do
      print *,'# of bonds: ',nbnd
c      print *,'bond orders: '
      do i = 1,natm
c        write(6,*)i,bndorder(i),atlab(i)
        do j = 1,bndorder(i)
          k = bndlist(j,i)
c          print *,'- ',atlab(k)
        end do
      end do
      return
900      print *,'cannot find crd file ',crdfil
      return
      end
