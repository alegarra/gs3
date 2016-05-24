 
      subroutine sortqr(a,n,m,ix,iy,col,icol)
c sorts matrix of reals by icol columns given by col, in ascending
c order. Matrix a is n x m, only first ix rows and iy columns are sorted
c !! same as sortqm except rows and columns switched 
      use kinds
      integer :: i,j,s,l,r,j1,j2,k,a1,a2,iy,icol,ix,m,n
      real (rh) ::a(:,:),x(20),x1
      integer :: stack(50,2),col(:)

      s=1
      stack(1,1)=1
      stack(1,2)=ix
10    l=stack(s,1)
      r=stack(s,2)
      s=s-1
20    i=l
      j=r
      j1=(l+r)/2
      do  k=1,icol
         x(k)=a(col(k),j1)
      enddo
30    do  k=1,icol
         a1=a(col(k),i)
         a2=x(k)
         if (a1.lt.a2) then
            i=i+1
            goto 30
         else if(a1.gt.a2) then
            goto 40
         endif
      enddo
40    do k=1,icol
         a1=a(col(k),j)
         a2=x(k)
         if (a1.gt.a2) then
            j=j-1
            goto 40
         else if(a1.lt.a2) then
            goto 47
         endif
      enddo
47    if (i.le.j) then
         do j2=1,iy
            x1=a(j2,i)
            a(j2,i)=a(j2,j)
            a(j2,j)=x1
         enddo
         i=i+1
         j=j-1
      endif
      if (i.le.j) goto 30
      if (i.lt.r) then
         s=s+1
         stack(s,1)=i
         stack(s,2)=r
      endif
      r=j
      if (l.lt.r) goto 20
      if (s.ne.0) goto 10
      end
 
      integer function hashvr(dat,r,ind,m,n,mode,nr)
c returns address of  element of ind()  containing real vector dat.
c Searching goes with a HASH technique. When mode=1 and dat was not
c there before, dat is written  into ind("hash",1..r). When mod=0
c and dat is not found, this function returns 0.
c rank of ind is m x n,
c dat contains r numbers.
c !! same as hashvec except that rows and columns switched !!
c !! same as hashves except that on low memory returns -1 and reals
      use kinds
      integer r,mode,nr,m,n
      real (rh) :: ind(:,:),dat(r)
      integer iaddr,prime(10),ie,iaddress,k,izer,ieq,i1,i
!AL      data ie/17/,prime/433,1039,53,761,197,31,887,389,277,97/
      data ie/17/,prime/5813,1039,6679,16603,9803,2063,887,389,277,97/ !this one is pasted from hashv1
      
c        
      iaddr=0
      do  i=1,r
         iaddr=iaddr+prime(i)*dat(i)
      enddo
      iaddress=abs(mod(iaddr,m))+1
 
      do  k=1,5000
         izer=0
         ieq=0
         do  i=1,r
            i1=ind(i,iaddress)
            if (i1.ne.dat(i)) ieq=1
            if (ind(i,iaddress).ne.0) izer=1
          enddo
 
 
         if (izer.eq.0 .or. ieq.eq.0) then
             if (izer.eq.0 .and. mode.eq.1) then
                 do  i=1,r
                    ind(i,iaddress)=dat(i)
                 enddo
                 nr=nr+1
             endif
             if (mode.eq.0 .and.izer.eq.0) then
                 hashvr=0
               else
                 hashvr=iaddress
             endif
             return
         endif
         iaddress=mod(iaddress+ie-1,m)+1
      enddo
 
      hashvr=-1
      return
      end                              
      
      subroutine clusterr(x,m1,m2,ix)
c clusters ix rows of real x containing nonzero last element
c at the beginning of the matrix
c !! same as cluster except rows and columns switched !!
      use kinds
      integer :: m1,m2,i,ix,j1,k1,k2
      real (rh):: x(:,:)


      j1=1
      k1=1
      do  i=1,ix

130      if (j1.le.m1) then
             if (x(1,j1).ne.0 ) then
                j1=j1+1
                goto 130
              endif
         endif

140      if (k1.le.m1) then
             if (x(1,k1).eq.0 ) then
                k1=k1+1
                goto 140
             endif
         endif

         if (k1.gt.j1) then
           do  k2=1,m2
             x(k2,j1)=x(k2,k1)
             x(k2,k1)=0
           enddo
         else
           k1=k1+1
         endif
      enddo
      end

      subroutine hashia1(x,nhash,n,ia,ja,a,m)
c copies data from hash to ia-ja-a form
c Entries are sorted within rows
      use kinds
      integer nhash,n,m,ia(:),ja(:),tmp(n)
      real (rh):: x(:,:)
      real (r8):: a(:)
      integer i,j,k,cut

      cut=0
c
c count the number of entries in each column of a
      do i=1,n
         tmp(i)=0
      enddo
      do i=1,nhash
         j=x(1,i)
         if (j.ne.0) then
            if (j.gt.n) then
                 x(1,i)=0
                 cut=cut+1
            elseif (x(2,i).gt.n) then
                 x(1,i)=0
                 cut=cut+1
            else
                 tmp(j)=tmp(j)+1
            endif
         endif
      enddo

c create the row count for b
      ia(1)=1
      do i=1,n
         ia(i+1)=ia(i)+tmp(i)
      enddo

      if (ia(n+1)-1.gt.m) then
         print*,'Too small parameter m in hashia:should be > ', ia(n+1)
         stop
      endif

c load a into b
      do i=1,n
         tmp(i)=0
      enddo
      do i=1,nhash
            j=x(1,i)
            if (j.ne.0) then
            k=ia(j)+tmp(j)
            ja(k)=x(2,i)
            a(k)=x(3,i)
            tmp(j)=tmp(j)+1
         endif
      enddo
c
      if (cut.gt.0) then
          print*,cut,
     +        ' elements in HASH->IJA conversion had index(es) over ',n
      endif
      end
      
      integer function hashv1(dat,r,ind,m,n,mode,nr)
c returns address of  element of ind()  containing integer vector dat.
c Searching goes with a HASH technique. When mode=1 and dat was not
c there before, dat is written  into ind("hash",1..r). When mod=0
c and dat is not found, this function returns 0.
c rank of ind is m x n,
c dat contains r numbers.
c !! same as hashvec except that rows and columns switched !!
c !! same as hashves except that on low memory returns -1 
      integer r,mode,nr,m,n
      integer :: ind(:,:),dat(r)
      integer iaddr,prime(10),ie,iaddress,k,izer,ieq,i1,i
      data ie/17/,prime/5813,1039,6679,16603,9803,2063,887,389,277,97/
c        
      iaddr=0
      do  i=1,r
         iaddr=iaddr+prime(i)*dat(i)
      enddo
      iaddress=abs(mod(iaddr,m))+1
 
      do  k=1,10000
         izer=0
         ieq=0
         do  i=1,r
            i1=ind(i,iaddress)
            if (i1.ne.dat(i)) ieq=1
            if (ind(i,iaddress).ne.0) izer=1
          enddo
 
 
         if (izer.eq.0 .or. ieq.eq.0) then
             if (izer.eq.0 .and. mode.eq.1) then
                 do  i=1,r
                    ind(i,iaddress)=dat(i)
                 enddo
                 nr=nr+1
             endif
             if (mode.eq.0 .and.izer.eq.0) then
                 hashv1=0
               else
                 hashv1=iaddress
             endif
             return
         endif
         iaddress=mod(iaddress+ie-1,m)+1
      enddo
 
      hashv1=-1
      print*,'hash table search too slow'
      stop
      return
      end                              

