c last modified Nov 23, 1994
c
c list of subroutines needed for FSPAK
c Miguel PEREZ-ENCISO
c Ignacy MISZTAL
c Mauricio ELZO
c
c second
c count_nza
c checkb
c checkst
c fastfb
c fullfb
c fgsfct
c figureout
c fratio
c fsinverse
c fsinverse1
c ival
c loadap
c loadu
c lzero
c makepath
c mkdet
c mkldet
c permute
c rowtocol
c rowtocolu
c rzero
c segpath
c spermute
c unfold
c zero
c zerout
c
      function count_nza(n,ia,ja,iout)
c     mpe, 17 January 1995
c counts no. of nze off-diagonal and checks whether there
c is any empty row.
      implicit none
      integer count_nza,n,i,iout,ia(1),ja(1),k,n0
      n0=0
      count_nza=0
      do i=1,n
         if(ia(i).eq.ia(i+1)) then
            n0=n0+1
         else
            do k=ia(i),ia(i+1)-1
               if(ja(k).ne.i) count_nza=count_nza+1
            enddo
         endif
      enddo
      if (n0.ne.0) write(iout,*) 'WARNING! ',n0,' empty rows'
      return
      end
c     -----------------
      subroutine checkb (n,ia,ja,a,d,b,x,tol,x1,perm,iout,flag)
c     -----------------
c     mpe, 18 Dec 1992 - 17 January 1995
      implicit none
      integer n,ia(1),ja(1),iout,flag,i,j,k,perm(1)
      real*8 dev,tol,a(1),b(1),d(1),x(1),x1(1)
      dev=0.
      do i=1,n
         x1(i)=0.
      enddo
      do i=1,n
         do k=ia(i),ia(i+1)-1
            j=ja(k)
            x1(i)=x1(i)+a(k)*b(j)
            if(i.lt.j) x1(j)=x1(j)+a(k)*b(i)
         enddo
      enddo
      do i=1,n
         WRITE(IOUT,*) I,X1(I),X(I)
         if(ia(i+1).gt.ia(i)) dev=dev+abs(x1(i)-x(i))/n
      enddo
      if(dev.gt.tol) then
          flag=1
          write(iout,*) 'ERROR in CHECKB: not accurate solution'
          write(iout,*) 'observed deviation = ',dev
          write(iout,*) 'maximum allowed    = ',tol
      endif
      return
      end

c     ------------------
      subroutine checkst (option,needed,maxneed,available,iout,flag)
c     ------------------
      integer option,needed,maxneed,available,iout,flag
      maxneed=max(maxneed,needed)
      if(needed.gt.available) then
c         write(iout,*) ' insufficient storage in option ',option
c         write(iout,*) ' avaliable ',available,'; needed ',needed
         flag=option*100+1
      endif
      return
      end

c     -------------------
      subroutine compress (nc,ipos,bc,b)
c     -------------------
c compreses nc positions of b into bc
      integer ipos(1)
      real*8 bc(1),b(1)
      do i=1,nc
         bc(i)=b(ipos(i))
      enddo
      return
      end

c     -----------------
      subroutine expand (nc,ipos,bc,b)
c     -----------------
c expand a compressed vector bc into b
c MPE
      integer ipos(1)
      real*8 bc(1),b(1)
      do i=1,nc
         b(ipos(i))=bc(i)
      enddo
      return
      end

c     -----------------
      subroutine fastfb
c     -----------------
     +                  (fpath,fstart,fnseg,bpath,bstart,bnseg,
     +                   iu,ju,iju,u,d,b)
c performs fast forward/backward on a triangular system iuju,iju,u
c M_P-E Thu Aug 13, 1992 14:34:41
      implicit none
      integer fpath(1),fstart(1),fnseg,bpath(1),bstart(1),bnseg,
     +        iu(1),ju(1),iju(1),i,j,k,ik,iseg
      real*8 u(1),d(1),b(1),bi
c fast forward in order within segments with reverse order btw segments
      do iseg=fnseg,1,-1
         do ik=fstart(iseg),fstart(iseg+1)-1
            i=fpath(ik)
            bi=b(i)
            do k=iu(i),iu(i+1)-1
               j=ju(iju(i)+k-iu(i))
               b(j)=b(j)+u(k)*bi
            enddo
          enddo
c fast diagonal
          do ik=fstart(iseg),fstart(iseg+1)-1
             i=fpath(ik)
             b(i)=b(i)*d(i)
          enddo
       enddo
c fast backward in reverse order within segments with reverse order btw seg
c segm
       do iseg=1,bnseg
          do ik=bstart(iseg+1)-1,bstart(iseg),-1
             i=bpath(ik)
             bi=b(i)
             do k=iu(i),iu(i+1)-1
                j=ju(iju(i)+k-iu(i))
                bi=bi+u(k)*b(j)
             enddo
             b(i)=bi
          enddo
       enddo
       return
       end

c     -----------------
      subroutine fullfb (n,iu,ju,iju,u,d,b)
c     -----------------
c performs full forward/backward on a triangular system iuju,iju,u
c M_P-E Sat Jul 11, 1992 17:54:59
      implicit none
      integer n,iu(1),ju(1),iju(1),i,j,k
      real*8 u(1),d(1),b(1),bi
c full forward
      do i=1,n
         bi=b(i)
         do k=iu(i),iu(i+1)-1
            j=ju(iju(i)+k-iu(i))
            b(j)=b(j)+u(k)*bi
         enddo
      enddo
c diagonal
       do i=1,n
          b(i)=b(i)*d(i)
       enddo
c full backward
       do i=n,1,-1
          bi=b(i)
          do k=iu(i),iu(i+1)-1
             j=ju(iju(i)+k-iu(i))
             bi=bi+u(k)*b(j)
          enddo
          b(i)=bi
       enddo
       return
       end

c     ---------------
      function fratio()
c     ---------------
c returns the ratio of storage occupied by reals and integers;
c the only two values returned are 1 and 2. Orginally written to
c facilitate easy porting between Crays (ratio=1) to other
c computers (ratio=2, if reals are real*8)
c Ignacy Misztal
      integer fratio,ix(2),const
      parameter(const=12345)
      real*8 y
      equivalence (y,ix(1))
      ix(2)=const
      y=0
      if (ix(2).eq.const) then
         fratio=1
        else
         fratio=2
      endif
      end

csppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppp
c Modified for sp semidefinite matrices (mae,1993)
c (partial inclusion of modifications done by Kachman (UNE)
c  on subr. gsfct as described by Curt Finley (UCD), 1993)
C Slightly modified and restructured
C M_P-E Sat Jul 11, 1992 16:50:38
C From George & Liu (1981)
C----- SUBROUTINE FGSFCT
C***************************************************************
C***************************************************************
C******     GSFCT ..... GENERAL SPARSE SYMMETRIC FACT     ******
C***************************************************************
C***************************************************************
C
C     PURPOSE - THIS SUBROUTINE PERFORMS THE SYMMETRIC
C        FACTORIZATION FOR A GENERAL SPARSE SYSTEM, STORED IN
C        THE COMPRESSED SUBSCRIPT DATA FORMAT.
C                                                                         1
C     INPUT PARAMETERS -                                                  1
C        NEQNS - NUMBER OF EQUATIONS.                                     1
C        XLNZ - INDEX VECTOR FOR LNZ.  XLNZ(I) POINTS TO THE              1
C               START OF NONZEROS IN COLUMN I OF FACTOR L.                1
C        (XNZSUB, NZSUB) - THE COMPRESSED SUBSCRIPT DATA                  1
C               STRUCTURE FOR FACTOR L.                                   1
C                                                                         1
C     UPDATED PARAMETERS -                                                1
C        LNZ - ON INPUT, CONTAINS NONZEROS OF A, AND ON                   1
C               RETURN, THE NONZEROS OF L.                                2
C        DIAG - THE DIAGONAL OF L OVERWRITES THAT OF A.                   2
C        IFLAG - THE ERROR FLAG.  IT IS SET TO 1 IF A ZERO OR             2
C               NEGATIVE SQUARE ROOT OCCURS DURING THE                    2
C               FACTORIZATION.                                            2
C        OPS   - A DOUBLE PRECISION COMMON PARAMETER THAT IS              2
C                INCREMENTED BY THE NUMBER OF OPERATIONS                  2
C                PERFORMED BY THE SUBROUTINE.                             2
C                                                                         2
C     WORKING PARAMETERS -                                                2
C        LINK - AT STEP J, THE LIST IN                                    3
C                  LINK(J), LINK(LINK(J)), ...........                    3
C               CONSISTS OF THOSE COLUMNS THAT WILL MODIFY                3
C               THE COLUMN L(*,J).                                        3
C        FIRST - TEMPORARY VECTOR TO POINT TO THE FIRST                   3
C               NONZERO IN EACH COLUMN THAT WILL BE USED                  3
C               NEXT FOR MODIFICATION.                                    3
C        TEMP - A TEMPORARY VECTOR TO ACCUMULATE MODIFICATIONS.           3
C                                                                         3
C***************************************************************          3
C                                                                         4
cmae
      SUBROUTINE  FGSFCT ( NEQNS, XLNZ, LNZ, XNZSUB, NZSUB, DIAG,
     1                     LINK, FIRST, TEMP, IFLAG,                      4
     2                     tol,irank)
C                                                                         4
C***************************************************************          4
C                                                                         4
cmae
         integer*4 irank
         real*8 tol
         REAL*8 COUNT, OPS
         COMMON  /SPKOPS/ OPS                                             4
         REAL*8 DIAG(1), LNZ(1), TEMP(1), DIAGJ, LJK
         INTEGER LINK(1), NZSUB(1)                                        4
         INTEGER FIRST(1), XLNZ(1), XNZSUB(1),                            5
     1           I, IFLAG, II, ISTOP, ISTRT, ISUB, J,                     5
     1           K, KFIRST, NEQNS, NEWK                                   5
C                                                                         5
C***************************************************************          5
C                                                                         5
cmae
c        ------------------------------
c        initialize irank and dmin
c        ------------------------------
         irank=0
!         dmin=1.d-08
         dmin=1d-20
! AL	 
cmae
C        ------------------------------                                   5
C        INITIALIZE WORKING VECTORS ...                                   5
C        ------------------------------                                   5
         DO I = 1, NEQNS
            LINK(I) = 0                                                   6
            TEMP(I) = 0.0                                                 6
         ENDDO
C        --------------------------------------------                     6
C        COMPUTE COLUMN L(*,J) FOR J = 1,...., NEQNS.                     6
C        --------------------------------------------                     6
         DO J = 1, NEQNS
C           -------------------------------------------                   6
C           FOR EACH COLUMN L(*,K) THAT AFFECTS L(*,J).                   6
C           -------------------------------------------                   6
            DIAGJ = 0.0                                                   7
            NEWK = LINK(J)                                                7
            K    = NEWK
            DO WHILE (K.NE.0)
               NEWK = LINK(K)                                             7
C              ---------------------------------------                    7
C              OUTER PRODUCT MODIFICATION OF L(*,J) BY                    7
C              L(*,K) STARTING AT FIRST(K) OF L(*,K).                     7
C              ---------------------------------------                    7
               KFIRST = FIRST(K)                                          7
               LJK    = LNZ(KFIRST)                                       8
               DIAGJ = DIAGJ + LJK*LJK                                    8
               OPS  = OPS + 1.00                                          8
               ISTRT = KFIRST + 1                                         8
               ISTOP = XLNZ(K+1) - 1                                      8
               IF ( ISTOP .GE. ISTRT )  THEN
C                 ------------------------------------------              8
C                 BEFORE MODIFICATION, UPDATE VECTORS FIRST,              8
C                 AND LINK FOR FUTURE MODIFICATION STEPS.                 8
C                 ------------------------------------------              8
                  FIRST(K) = ISTRT                                        9
                  I = XNZSUB(K) + (KFIRST-XLNZ(K)) + 1                    9
                  ISUB = NZSUB(I)                                         9
                  LINK(K) = LINK(ISUB)                                    9
                  LINK(ISUB) = K                                          9
C                 ---------------------------------------                 9
C                 THE ACTUAL MOD IS SAVED IN VECTOR TEMP.                 9
C                 ---------------------------------------                 9
                  DO II = ISTRT, ISTOP                                    9
                     ISUB = NZSUB(I)                                      9
                     TEMP(ISUB) = TEMP(ISUB) + LNZ(II)*LJK               10
                     I = I + 1                                           10
                  ENDDO                                                  10
                  COUNT = ISTOP - ISTRT + 1                              10
                  OPS  = OPS + COUNT
               ENDIF
               K    = NEWK
            ENDDO
C           ----------------------------------------------               10
C           APPLY THE MODIFICATIONS ACCUMULATED IN TEMP TO               10
C           COLUMN L(*,J).                                               10
C           ----------------------------------------------               10
            DIAGJ = DIAG(J) - DIAGJ                                      11
cmae
c           IF ( DIAGJ .LE. 0.0E0 )  THEN
C              ------------------------------------------------------
C              ERROR - ZERO OR NEGATIVE SQUARE ROOT IN FACTORIZATION.
C              ------------------------------------------------------
c              IFLAG = J
c              RETURN
c           ENDIF
c           DIAGJ = SQRT(DIAGJ)
cmae        DIAG(J) = DIAGJ
cmae
c  this is the core of Kachman's modification to subr. gsfct
            if (diag(j).ge.dmin) then
               if (diagj.ge.(tol*diag(j))) then
                  irank=irank+1
                  diagj = sqrt(diagj)
                  diag(j) = diagj
                  diagj=1.d0/diagj
               else
                  diag(j)=0.0d0
                  diagj=0.0d0
               endif
            else
               diag(j)=0.0d0
               diagj=0.0d0
            endif
cmae
            ISTRT = XLNZ(J)                                              11
            ISTOP = XLNZ(J+1) - 1                                        11
            IF ( ISTOP .GE. ISTRT )  THEN
               FIRST(J) = ISTRT                                          11
               I = XNZSUB(J)                                             11
               ISUB = NZSUB(I)                                           11
               LINK(J) = LINK(ISUB)                                      12
               LINK(ISUB) = J                                            12
c dir$ ivdep
               DO II = ISTRT, ISTOP
                  ISUB = NZSUB(I)                                        12
cmae              LNZ(II) = ( LNZ(II)-TEMP(ISUB) ) / DIAGJ               12
                  LNZ(II) = ( LNZ(II)-TEMP(ISUB) ) * DIAGJ               12
                  TEMP(ISUB) = 0.0E0                                     12
                  I = I + 1                                              12
               ENDDO
               COUNT = ISTOP - ISTRT + 1                                 12
               OPS  = OPS + COUNT
         ENDIF
      ENDDO
      RETURN
      END                                                                13

c     --------------------
      subroutine figureout (n,avail,iout,st,irank)
c     --------------------
c     mpe
      integer n,avail,irank
      real*8 st(*)
      write(iout,'(a,f25.6)')'               **************'
      write(iout,'(a,f25.6)')'               **** FSPAK ***'
      write(iout,'(a,f25.6)')'               **************'
      write(iout,'(a,f25.6)')'               MPE / IM / MAE'
      write(iout,'(a,f25.6)')'                   Jun 1994'
      write(iout,'(a,f25.6)')
      write(iout,'(a,f25.6)')'              SPARSE STATISTICS'
      write(iout,'(a,I25)')  '      DIMENSION OF MATRIX     =',n
      write(iout,'(a,I25)')  '      RANK                    =',irank
      write(iout,'(a,I25)')  '      STORAGE AVAILABLE       =',avail
      write(iout,'(a,I25)')  '      MAXIMUM NEEDED          ='
     +                                                  ,IVAL(ST(3))
      write(iout,'(a,I25)')  '      NZE IN UPPER TRIANGULAR ='
     +                                                ,IVAL(ST(4))+N
      write(iout,'(a,I25)')  '      NZE IN FACTOR           ='
     +                                                  ,IVAL(ST(2))
      write(iout,'(a,i25)')  '      NO. OF CALLS NUM FACT   ='
     +                                                  ,INT(ST(11))
      write(iout,'(a,i25)')  '      NO. OF CALLS SOLVE      ='
     +                                                  ,INT(ST(12))
      write(iout,'(a,i25)')  '      NO. OF CALLS SPARS SOLV ='
     +                                                  ,INT(ST(13))
      write(iout,'(a,i25)')  '      NO. OF CALLS DET / LDET ='
     +                                                  ,INT(ST(14))
      write(iout,'(a,i25)')  '      NO. OF CALLS SPARS INV  ='
     +                                                  ,INT(ST(15))
      if(st(13).ne.0) then
      write(iout,'(a,i25)')  '      MAX NO. NODES           ='
     +                                                  ,INT(ST(30))
      write(iout,'(a,i25)')  '      MAX PATH LENGTH         ='
     +                                                  ,INT(ST(31))
      write(iout,'(a,f25.3)')'      AVG NO. NODES FPATH     ='
     +                                                ,ST(32)/ST(13)
      write(iout,'(a,f25.3)')'      AVG NO. NODES BPATH     ='
     +                                                ,ST(33)/ST(13)
      write(iout,'(a,f25.3)')'      AVG FPATH LENGTH        ='
     +                                                ,ST(34)/ST(13)
      write(iout,'(a,f25.3)')'      AVG BPATH LENGTH        ='
     +                                                ,ST(35)/ST(13)
      endif
      write(iout,'(a,f25.6)')'      TOTAL CPU TIME IN FSPAK =',st(20)
      write(iout,'(a,f25.6)')'      TIME FOR FINDING ORDER  =',st(21)
      write(iout,'(a,f25.6)')'      TIME FOR SYMBOLIC FAC   =',st(22)
      write(iout,'(a,f25.6)')'      TIME FOR NUMERICAL FAC  =',st(24)
      write(iout,'(a,f25.6)')'      TIME FOR SOLVE          =',st(25)
      write(iout,'(a,f25.6)')'      TIME FOR SPARSE SOLVE   =',st(27)
      write(iout,'(a,f25.6)')'      TIME FOR SPARSE INVERSE =',st(26)
      return
      end

c     -------------
      function ival (st)
c     -------------
      integer st
      ival=st
      return
      end

c     -----------------
      subroutine loadap
c     -----------------
     +                  (n,ia,ja,a,iap,jap,ap,d,ip,tmp)
c this subroutine loads half sparse stored elements of ia,ja,a into half sp
c sparse stored
c structure iap,jap,ap for non-diagonal elements and into d for diagonal el
C elements;
c elements are reordered according to permutation vectors p,ip
c M_P-E Tue Jun 30, 1992 10:23:01
      implicit none
      integer n,ia(1),ja(1),iap(1),jap(1),ip(1),tmp(1),i,j,k,kk
      real*8 a(1),ap(1),d(1)
      do i=1,n
         tmp(i)=0
      enddo
c compute no. of entries per upper stored row
      do k=1,n
         i=ip(k)
         do kk=ia(k),ia(k+1)-1
            j=ip(ja(kk))
            if(i.eq.j .and. i.le.n) then
               d(i)=a(kk)
            elseif(i.gt.j .and. i.le.n) then
               tmp(j)=tmp(j)+1
            elseif(i.lt.j .and. j.le.n) then
               tmp(i)=tmp(i)+1
            endif
         enddo
      enddo
c fill iap
      iap(1)=1
      do i=1,n
         iap(i+1)=iap(i)+tmp(i)
         tmp(i)=0
      enddo
c fill jap,ap
      do k=1,n
         i=ip(k)
         do kk=ia(k),ia(k+1)-1
            j=ip(ja(kk))
            if(i.gt.j .and. i.le.n) then
               jap(iap(j)+tmp(j))=i
               ap(iap(j)+tmp(j))=a(kk)
               tmp(j)=tmp(j)+1
            elseif(i.lt.j .and. i.le.n) then
               jap(iap(i)+tmp(i))=j
               ap(iap(i)+tmp(i))=a(kk)
               tmp(i)=tmp(i)+1
            endif
         enddo
      enddo
      return
      end

c     ----------------
      subroutine loadu (n,iap,jap,ap,iu,ju,iju,u,tmp)
c     ----------------
c sparse stored matrix iap,jap,ap is copied into sparse compressed storage
C GE MATRIX IU,JU,IJU,U
c where u allows storage for fill-ins
c M_P-E Tue Jun 30, 1992 11:01:23
      implicit none
      integer n,iap(1),jap(1),iu(1),ju(1),iju(1),i,j,k
      real*8 ap(1),u(1),tmp(1)
      do i=1,n
         tmp(i)=0.
      enddo
      do i=1,n
c - -    scatters
         do k=iap(i),iap(i+1)-1
            j=jap(k)
            tmp(j)=ap(k)
         enddo
c - -    loads u and zero out tmp
         do k=iu(i),iu(i+1)-1
            j=ju(iju(i)+k-iu(i))
            u(k)=tmp(j)
            tmp(j)=0.
         enddo
      enddo
      return
      end

c     ----------------
      subroutine lzero (n,l)
c     ----------------
c initialize a logical vector
c MPE
      logical l(1)
      do i=1,n
         l(i)=.false.
      enddo
      return
      end

c     -------------------
      subroutine makepath
c     -------------------
     +           (n,path,iu,ju,iju)
c obtains factorization path stored as a linked list where path(i)
c is next node of the path of element i
c M. Perez-Enciso, Madison,Thu May 21,1992
      integer n,path(n),iu(n+1),ju(1),iju(n)
      do i=1,n
         path(i)=0
         if(iu(i+1).gt.iu(i)) path(i)=ju(iju(i))
      enddo
      return
      end

c     -----------------
      subroutine merlin (n,iu,u,d,tol)
c     -----------------
c MPE
c mae
c skip zero diagonals (machine zero = tol)
      integer n,iu(1)
      real*8 u(1),d(1),tol
      do i=1,n
         if(d(i).gt.tol)then
            do j=iu(i),iu(i+1)-1
               u(j)=-u(j)/d(i)
            enddo
            d(i)=1./(d(i)*d(i))
         endif
      enddo
      return
      end

c     ----------------
      subroutine mkdet (n,d,b,irank,tol,iout,flag)
c     ----------------
c computes determinant, M_P-E
c 11/22/99 IM
c mae
c skip zero diagonals
      integer i,iout,irank,flag,n
      real*8 b,d(1),tol
      b=1.
      do i=1,n
         if(d(i).ne.0)then
            b=b/d(i)
         endif
      enddo
      end

c     -----------------
      subroutine mkldet (n,d,b,irank,tol,iout,flag)
c     -----------------
c computes log determinant
c M_P-E Mon Aug 10, 1992 15:43:45
c Mon Feb 28, 1994 10:57:00 (revisited)
c 11/22/99 IM
c mae
c skip zero diagonals
      integer i,iout,irank,flag,n
      real*8 b,d(1),tol
      b=0.
      do i=1,n
         if(d(i).ne. 0)then
            b=b-log(d(i))
         endif
      enddo
      end

c     ------------------
      subroutine permute (n,vec1,vec2,p)
c     ------------------
c permutes vector vec1 according to permutation vector p into vec2
c M_P-E Wed Jul 1, 1992 11:12:01
      integer n,i,p(1)
      real*8 vec1(1),vec2(1)
      do i=1,n
         vec2(p(i))=vec1(i)
      enddo
      return
      end


c     -------------------
      subroutine rowtocol(aia,aja,aa,bia,bja,ba,tmp,n)
c     -------------------
c Matrix a, in sparse i-j-a form and stored row-wise, is converted to b,
c stored column-wise; after the move, the rows of b are sorted in
c increasing order.
c a and b are have n rows and columns, and tmp is
c is a temporary integer vector of size n.
c I. Misztal
      integer n,aia(1),aja(1),bia(1),bja(1),tmp(1),i,j,k
      real*8 aa(1),ba(1)

c count the number of entries in each row of a
      call zero(n,tmp)
      do i=1,n
         do j=aia(i),aia(i+1)-1
            tmp(aja(j))=tmp(aja(j))+1
         enddo
      enddo

c create the row count for b
      bia(1)=1
      do i=1,n
         bia(i+1)=bia(i)+tmp(i)
      enddo

c load a into b
      call zero(n,tmp)
      do i=1,n
         do j=aia(i),aia(i+1)-1
            k=bia(aja(j))+tmp(aja(j))
            bja(k)=i
            ba(k)=aa(j)
            tmp(aja(j))=tmp(aja(j))+1
         enddo
      enddo
      end

c     --------------------
      subroutine rowtocolu (aia,aju,aiju,aa,bia,bja,ba,tmp,n)
c     --------------------
c Matrix a, in sparse i-ju-iju -a form and stored row-wise, is converted to
C TO B,
c stored column-wise; after the move, the rows of b are sorted in
c increasing order.
c a and b are have n rows and columns, and tmp is
c is a temporary integer vector of size n.
c I. Misztal
      integer n,aia(1),bia(1),bja(1),tmp(1),i,j,k,
     *        aju(1),aiju(1)
      real*8 aa(1),ba(1)

c count the number of entries in each row of a
      call zero(n,tmp)
      do i=1,n
         do j=aia(i),aia(i+1)-1
            tmp(aju(aiju(i)-aia(i)+j))=tmp(aju(aiju(i)-aia(i)+j))+1
         enddo
      enddo

c create the row count for b
      bia(1)=1
      do i=1,n
         bia(i+1)=bia(i)+tmp(i)
      enddo

c load a into b
      call zero(n,tmp)
      do i=1,n
         do j=aia(i),aia(i+1)-1
            k=bia(aju(aiju(i)-aia(i)+j))+tmp(aju(aiju(i)-aia(i)+j))
            bja(k)=i
            ba(k)=aa(j)
            tmp(aju(aiju(i)-aia(i)+j))=tmp(aju(aiju(i)-aia(i)+j))+1
         enddo
      enddo
      end

c     ----------------
      subroutine rzero (n,b)
c     ----------------
c initialize a real vector
c MPE
      real*8 b(1)
      do i=1,n
         b(i)=0.
      enddo
      return
      end

c     ------------------
      subroutine segpath
c     ------------------
     +           (n,node,nseg,sstart,spath,path,found)
c computes a segmented path for the nseg nodes stored in node
c M. Perez-Enciso, Madison, Thu May 21, 1992
      implicit none
      integer n,nseg, node(nseg),sstart(nseg+1),spath(1),path(n)
     +        ,ik,k,iseg
      logical found(n)
      ik=0
      do iseg=1,nseg
         k=node(iseg)
         sstart(iseg)=ik+1
         do while(k.ne.0.and..not.found(k))
            ik=ik+1
            spath(ik)=k
            found(k)=.true.
            k=path(k)
         enddo
      enddo
      sstart(nseg+1)=ik+1
c zero out boolean working vector
      do ik=1,sstart(nseg+1)-1
         found(spath(ik))=.false.
      enddo
      return
      end

c     -------------------
      subroutine spermute (n,vect1,vect2,perm)
c     -------------------
c permutes an integer compressed vector
c MPE
      integer vect1(1),vect2(1),perm(1)
      do i=1,n
         vect2(i)=perm(vect1(i))
      enddo
      return
      end

c     -----------------
      subroutine unfold (n,uia,uja,ia,ja,tmp)
c     -----------------
c this subroutine copies a half sparse stored matrix into a sparse
c stored matrix, excluding diagonals
c option=0: no reordering performed, only ia and ja created
c       =1: reordering and sorting performed, a,d also created
c I. Misztal, M. Perez-Enciso
      implicit none
      integer n,uia(1),uja(1),ia(1),ja(1)
     +        ,i,j,k,tmp(1)
c count the number of entries above and below diagonal
      do i=1,n
         tmp(i)=0
      enddo
      do i=1,n
         do k=uia(i),uia(i+1)-1
            j=uja(k)
            if(i.ne.j .and. j.le.n .and. j.ne.0) then
               tmp(i)=tmp(i)+1
               tmp(j)=tmp(j)+1
            endif
         enddo
      enddo
c fill ia ( <--> xadj)
      ia(1)=1
      do i=1,n
         ia(i+1)=ia(i)+tmp(i)
         tmp(i)=0
      enddo
c fill ja ( <--> adjncy )
      do i=1,n
         do k=uia(i),uia(i+1)-1
            j=uja(k)
            if(j.ne.i .and. j.le.n .and. j.ne.0) then
               ja(ia(i)+tmp(i))=j
               ja(ia(j)+tmp(j))=i
               tmp(i)=tmp(i)+1
               tmp(j)=tmp(j)+1
            endif
         enddo
      enddo
      return
      end

c     ---------------
      subroutine zero (n,x)
c     ---------------
      integer n,x(n),i
      do i=1,n
         x(i)=0
      enddo
      end

c     -----------------
      subroutine zerout
c     -----------------
     +                  (b,path,sstart,nseg)
      integer path(1),sstart(1),nseg
      real*8 b(1)
      do i=1,sstart(nseg+1)-1
         b(path(i))=0.
      enddo
      return
      end
C IM 10/24/93-11/23/94
C***********************************************************************
C  Sinverse1 -- inverse OF SPARSE SYMMETRIC POSITIVE DEFINITE SYSTEM OF
C         LINEAR EQUATIONS  Z=inv(M)  GIVEN UT-D-U FACTORIZATION OF M;
C
C         Replaces a in input matrix ia-ja-a with the corresponding
c         inverse elements
C
C         Compared to Sinverse, this version has smaller memory
C         requirements and does not expand ia-ja-a to full storage
C***********************************************************************
      SUBROUTINE  fsinverse1
     *     (N, P, IP, D, IJU,JU,IU,U, TMPu,tmpz,tmpi,tmpi1,
     *      jja,ia,ja,a)
      implicit none
      INTEGER  P(1),  IP(1), IJU(1), JU(1), IU(1), TMPi1(1),
     *           tmpi(1),jja(1),ia(1),ja(1),n,j1,k1
      real*8  D(1), U(1), a(1), TMPu(1),tmpz(1)

      integer i,j,k,iiu
C
C  ADDITIONAL PARAMETERS
C
C    TMPu,tmpz   - real ONE-DIMENSIONAL WORK ARRAYs;  DIMENSION = N
C    jja - integer on-dimensional work array of order ia(n+1)-1
C    tmpi,tmpi1 - integer one-dimensional work array; dimension=N,
c                 tmpi and tmpi1 can occupy the same storage as
c                 tmpu and tmpz
C
C-----------------------------------------------------------------------
      iiu(i,j)=ju(iju(i)-iu(i)+j)
c

c  CAUTION: U is really -U

C----invert in site, z overwrites u
C            tmpu keeps unfolded current column of u;
C            tmpz keeps unfolded current column of z;
      do i=1,n
         tmpz(i)=0
         tmpu(i)=0
      enddo

      do i=n,1,-1
c dir$ ivdep
         do j=iu(i), iu(i+1)-1
            j1=iiu(i,j)
            tmpu(j1)=u(j)
            tmpz(j1)=u(j)*d(j1)
         enddo

c-------off-diagonal elements
         do j=iu(i),iu(i+1)-1
            j1=iiu(i,j)
c           tmpz(j1)=tmpz(j1)+tmpu(j1)*d(j1)
c dir$ ivdep
            do k=iu(j1),iu(j1+1)-1
               k1=iiu(j1,k)
               tmpz(j1)=tmpz(j1)+tmpu(k1)*u(k)
               tmpz(k1)=tmpz(k1)+tmpu(j1)*u(k)
            enddo
          enddo

c------diagonal element last
c dir$ ivdep
         do j=iu(i),iu(i+1)-1
            d(i)=d(i)+u(j)*tmpz(iiu(i,j))
         enddo

c store inverse and zero all nonzeroes in tmp and tmp1
c dir$ ivdep
         do j=iu(i), iu(i+1)-1
            tmpu(iiu(i,j))=0
            u(j)=tmpz(iiu(i,j))
            tmpz(iiu(i,j))=0
         enddo
       enddo

c replace a in ia-ja-a with inverse elements

      do i=1,n
         tmpi1(i)=i
      enddo
c first partly permute ia-ja-a
      call pperm(n,ip,ia,ja,a,tmpi,jja)
      do i=1,n
         do j=iu(ip(i)),iu(ip(i)+1)-1
            tmpu(p(iiu(ip(i),j)))=u(j)
         enddo
         tmpu(i)=d(ip(i))
         do j=ia(i),ia(i+1)-1
            a(j)=tmpu(ja(j))
         enddo
       enddo

c reverse the permutation
      call pperm(n,tmpi1,ia,ja,a,tmpi,jja)
      end


C IM - 10-19-1993-11/23/94
C*********************************************************************
C Partly permute upper-triangular matrix ia-ja-a using permutation vector
C p.
C The permutation is such that it does not change the indices of elements,
C but insures, that after permutation all elements would be in the
C upper diagonal.
C*********************************************************************
        subroutine pperm(n,p,ia,ja,a,tmpia,tmprow)
        implicit none
        integer n,p(1),ia(1),ja(1),tmpia(1),tmprow(1),i,j,k
        real*8 a(1),x

c extra parameters
c tmpia   - integer vector of size n
c tmprow  - integer vector of size ia(n+1)-n

c zero temporary row vector for reorder matrix
        do i=1,n
           tmpia(i)=0
        enddo
c set ja to store column entries and tmprow of row entries
        do i=1,n
           do j=ia(i),ia(i+1)-1
              if (p(i).gt.p(ja(j))) then
                  k=ja(j)
                  ja(j)=i
                  tmprow(j)=k
                 else
                   tmprow(j)=i
                   k=i
                endif
                tmpia(k)=tmpia(k)+1
           enddo
        enddo

c set ia to permuted ia
       do i=1,n
          ia(i+1)=ia(i)+tmpia(i)
          tmpia(i)=ia(i+1)
       enddo

c assign rows to appropriate addresses
        do i=ia(n+1)-1,ia(1),-1
           j=tmprow(i)
           tmpia(j)=tmpia(j)-1
           tmprow(i)=tmpia(j)
        enddo

c final permutation
        do i=1,ia(n+1)-1
10         j=tmprow(i)
           if (i.ne.j) then
c                swap entries i and j
              k=ja(i)
              ja(i)=ja(j)
              ja(j)=k
              k=tmprow(i)
              tmprow(i)=tmprow(j)
              tmprow(j)=k
              x=a(i)
              a(i)=a(j)
              a(j)=x
              goto 10
           endif
        enddo
        end
C  IM - 1/10/92 - 10/24/93
C***********************************************************************
C  Sinverse -- inverse OF SPARSE SYMMETRIC POSITIVE DEFINITE SYSTEM OF
C         LINEAR EQUATIONS  Z=inv(M)  GIVEN UT-D-U FACTORIZATION OF M;
C
C         Expands input matrix ia-ja-a to full storage and replaces
c         values of a with those of the inverse
C***********************************************************************
      SUBROUTINE  fsinverse
     *     (N, P, IP, D, IJU,JU,IU,U, TMP,za,zd,zja,
     *      utia,utja,tmpi,tmpi1,ia,ja,a)
      INTEGER  P(1),  IP(1), IJU(1), JU(1), IU(1), TMPi1(1),
     *          zja(1), tmpi(1),
     *         utia(1),utja(1),ia(1),ja(1)
      real*8  D(1), U(1),  za(1), a(1), TMP(1),zd(1)

      integer i,j,k,l,freeuz,nnzero
C
C  ADDITIONAL PARAMETERS
C
C    TMP   - real ONE-DIMENSIONAL WORK ARRAY;  DIMENSION = N
C    iu,zja,za - inverse, the lower triangular form, in the i-j-a
C                 form; diagonals in ZD
C    utia,utja - U in the i-j-a format, indices only; sorted within rows
C    tmpi,tmpi1 - integer one-dimensional work array; dimension=N
C
C-----------------------------------------------------------------------
C
c  CAUTION: U is really -U
c create indices for U transposed, sorted, za used as temporary
      call rowtocolu(iu,ju,iju,u,utia,utja,za,tmpi,n)

c create I index for Z the same as in U; J index is
c added during computations


C----invert, tmpi contains number of filled elements in each row of z,
C            tmp keeps unfolded current column of z;
c             remember: U is -U
      do i=1,n
         tmpi(i)=0
         tmp(i)=0
      enddo

      do i=n,1,-1
c------scatter existing elements of column i of z into tmp
         do j=iu(i), iu(i)+tmpi(i)-1
            tmp(zja(j))=za(j)
         enddo

c------diagonal element first
         freeuz=iu(i)+tmpi(i)
         tmp(i)=d(i)
         do j=iu(i),iu(i+1)-1
            tmp(i)=tmp(i)+u(j)*tmp(ju(iju(i)-iu(i)+j))
         enddo

         zd(i)=tmp(i)
         nnzero=1
         tmpi1(nnzero)=i

c-------off-diagonal elements
         do j=utia(i+1)-1,utia(i),-1
            freeuz=iu(utja(j))+tmpi(utja(j))
            l=utja(j)
c dir$ ivdep
            do k=iu(l),iu(l+1)-1
               tmp(l)=tmp(l)+u(k)*tmp(ju(iju(l)-iu(l)+k))
            enddo
            za(freeuz)=tmp(l)
            zja(freeuz)=i
            tmpi(utja(j))=tmpi(utja(j))+1
            nnzero=nnzero+1
            tmpi1(nnzero)=l
          enddo

c zero all nonzeroes in tmp
         do j=1,nnzero
            tmp(tmpi1(j))=0
         enddo
       enddo


c in Z, zero all entries not in input matrix A
c ----- since A is upper diagonal, create index structure for A transpose
      call rowtocol(ia,ja,a,utia,utja,a,tmpi,n)

       call zero(n,tmpi)
       call zero(n,tmpi1)
       do i=1,n
c  --------scatter row i of A (upper + lower triangle
          do j=ia(i),ia(i+1)-1
             tmpi(ja(j))=1
          enddo
          do j=utia(i),utia(i+1)-1
             tmpi(utja(j))=1
          enddo
          do j=iu(ip(i)),iu(ip(i)+1)-1
             if (tmpi(p(zja(j))).eq.0) then
                zja(j)=0
             endif
          enddo
c--------zero tmpi and count number of entries in upper+lower diagonal  a
          do j=ia(i),ia(i+1)-1
             tmpi(ja(j))=0
             k=ja(j)
             if (i.ne.k) tmpi1(i)=tmpi1(i)+1
             tmpi1(k)=tmpi1(k)+1
          enddo
          do j=utia(i),utia(i+1)-1
             tmpi(utja(j))=0
          enddo
       enddo


c create row structure for upper+lower diagonal a
      ia(1)=1
      do i=1,n
         ia(i+1)=ia(i)+tmpi1(i)
      enddo

c copy z to a, from half to full stored
      call zero(n,tmpi)
      do i=1,n
         k=p(i)
         a(ia(k)+tmpi(k))=zd(i)
         ja(ia(k)+tmpi(k))=k
         tmpi(k)=tmpi(k)+1
         do j=iu(i),iu(i+1)-1                   
           l=p(i)
           if (zja(j).ne.0) then
                k=p(zja(j)) 
                a(ia(l)+tmpi(l))=za(j)
                a(ia(k)+tmpi(k))=za(j)
                ja(ia(l)+tmpi(l))=k
                ja(ia(k)+tmpi(k))=l
                tmpi(l)=tmpi(l)+1
                tmpi(k)=tmpi(k)+1
            endif
         enddo
      enddo
      return
      end
