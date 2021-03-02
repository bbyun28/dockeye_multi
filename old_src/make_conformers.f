      program make_conformers
c
c read pdb file, generate bonds, find movable tosrions, generate conformers
c write new file
c
      implicit none
      integer maxatm,maxord
      parameter ( maxatm = 10000 )
      parameter ( maxord = 10 )
      integer bndorder(maxatm)
      integer bndlist(maxord,maxatm)
      integer nsubset,isubset(maxatm)
      integer torlist(2,maxatm),ntor
      integer torset(4,maxatm)
      integer nbnd,natm
      real*4 atxyz(3,maxatm)
      integer lunit
      character*80 crdfil
      character*14 atlab(maxatm)
      integer i,j,k
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
c find torsions
c
      call findtor(natm,bndorder,bndlist,atlab,ntor,torlist)
c
c generate sets of 4 atoms that define a torsion angle
c
      do i = 1,ntor
        j = torlist(1,i)
        torset(2,i) = j
        torset(3,i) = torlist(2,i)
      end do
      stop
      end
c-------------------------------------------------------------------
      subroutine findtor(natm,bndorder,bndlist,atlab,ntor,torlist)
      implicit none
      integer maxatm,maxord
      parameter ( maxatm = 10000 )
      parameter ( maxord = 10 )
      integer bndorder(maxatm),bndorder1(maxatm)
      integer bndlist(maxord,maxatm)
      integer torlist(2,maxatm),ntor
      integer nlevel(maxatm),ilevel
      integer natm
      integer ijoin(maxatm)
      integer i,j,k
      integer ndone,inew
      integer i1,i2,i3,i4,i5,i6
      character*14 atlab(maxatm)
c---------------------
      ntor = 0
c recursively prune branches in bodning tree by removing atoms that have just 1 neighbor
c if  they are not leaves (first level of pruning) then there is a torsion bond
      ilevel = 1
      do i = 1,natm
        nlevel(i) = 0
        bndorder1(i) =  bndorder(i) ! working copy of bond orders
      end do
200   continue
        ndone = 0
        do i = 1,natm
          if(nlevel(i).eq.0)then  ! if unassigned atom
            if(bndorder1(i).eq.1)then ! if now leaf of tree
              nlevel(i) = ilevel
              ndone = ndone + 1
              do j = 1,bndorder(i) ! decrease its neigbors bond order
                k = bndlist(j,i)
                if(nlevel(k).eq.0)then
                  bndorder1(k) = bndorder1(k) - 1
                  ndone = ndone + 1
                  if(ilevel.gt.1)then
                    ntor = ntor + 1
                    torlist(1,ntor) = i
                    torlist(2,ntor) = k
                  end if
                end if
              end do
            end if
          end if
        end do
        ilevel = ilevel + 1
      if(ndone.gt.0) goto 200
c
c all branches gone everything has bond order >=2 and is part of a ring
c or between two (sets of) rings
c
c      do i = 1,natm
c        print *,i,nlevel(i),bndorder1(i),atlab(i)
c      end do
c
c for every remaining bond, if it is a possible torsion it will divide atoms into two disjoint sets
c
      do i = 1,natm
        if(bndorder1(i).ge.2)then
          j = bndorder(i)
          do k = 1,j
            i1 = bndlist(k,i)
            if((i1.gt.i).and.(bndorder1(i1).ge.2))then
c              print *,'possible inter-ring torsion bond: ',atlab(i),atlab(i1)
c label one of the atoms, and keep labeling neighbors, neighbors of neighbors
c until done, but cannot cross i-i1 bond
              do i3 = 1,natm
                ijoin(i3) = 0
              end do
              ndone = 1
              inew = 1
              ijoin(i) = 1
              do while(inew.ne.0)
                inew = 0
                do i3 = 1,natm
                  i4 = bndorder(i3)
                  do i5 = 1,i4
                    i6 = bndlist(i5,i3)
                    if((i3.eq.i).and.(i6.eq.i1))then
                    else
                      if((i3.eq.i1).and.(i6.eq.i))then
                      else
                        if((ijoin(i3).eq.1).and.(ijoin(i6).eq.0))then
                          ijoin(i6) = 1
                          inew = 1
                          ndone = ndone + 1
                        end if
                        if((ijoin(i3).eq.0).and.(ijoin(i6).eq.1))then
                          ijoin(i3) = 1
                          inew = 1
                          ndone = ndone + 1
                        end if
                      end if
                    end if
                  end do
                end do
              end do
              if(ndone.lt.natm)then
                print *,'possible inter-ring torsion bond: ',atlab(i),atlab(i1)
                print *,ndone
                ntor = ntor + 1
                torlist(1,ntor) = i
                torlist(2,ntor) = i1
              end if
            end if
          end do
        end if
      end do
      print *,'# of torsions: ',ntor
      do i = 1,ntor
        j = torlist(1,i)
        k = torlist(2,i)
        print *,atlab(j),atlab(k)
      end do
c      while(ndone.lt.natm)
c        ilevel = ilevel + 1
c        end do
c      end do
      end
c-------------------------------------------------------------------
      subroutine findseg(natm,nbnd,bndorder,bndlist,nsubset,isubset)
      implicit none
      integer maxatm,maxord
      parameter ( maxatm = 10000 )
      parameter ( maxord = 10 )
      integer bndorder(maxatm)
      integer bndlist(maxord,maxatm)
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
      integer maxatm,maxord
      parameter ( maxatm = 10000 )
      parameter ( maxord = 10 )
      integer bndorder(maxatm)
      integer bndlist(maxord,maxatm)
      integer nbnd,natm
      integer nhyd
      integer lunit
      character*80 crdfil,line
      character*14 atlab(maxatm)
      real*4 atrad(maxatm),atxyz(3,maxatm)
      real*4 dist,slop
c------------------------------------
      integer i,j,k
c------------------------------------
      character atnam*6
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
