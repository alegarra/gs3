c      ----------------
       subroutine fspak
c      ----------------
     +                  (option,n,ia,ja,a,b,flag,iout,ioor,
     +                   available,needed,is,fnode,bnode,
     +                   fnseg,bnseg,fx,feqb,irank)
c -----------------------------------------------------------
c Public Sparse Matrix Package for Symmetric Matrices
c Authors: Miguel Perez-Enciso, Ignacy Misztal, Mauricio Elzo
c Mon Aug 24, 1992 12:39:46 - Jan 17, 1995
c -----------------------------------------------------------
      implicit none
      integer zero,order,symfac,numfac,solve,sinit,ssolve
     +       ,det,ldet,check,spin,spin1,restart,stat,equal,mx_st
      parameter (order=10,  symfac=20,  numfac=40, solve=50,
     +           sinit=51,  ssolve=52,  det=54,    ldet =55,
     +           check=56,  spin=60,  spin1=61, restart=70,
     +           stat=80,  equal=1,   mx_st=50 )
      integer available,r_avail
     +       ,delta,flag,ioor,iout,n,ia(1),ja(1),is(available)
     +       ,fnode(1),bnode(1),fnseg,bnseg,feqb,irank,ratio,fratio
     +       ,path,found,fpath,fstart,fpnode,bpath,bstart,bpnode
     +       ,maxju,maxu,needed,maxneed,i,adjncy,d,invp,iu,iju,ju,u
     +       ,maxint,nza,nzad,option,perm,iap,jap,ap,xadj,first,ios
     +       ,zja,za,zd,utia,utja,work1,work2,work3,work4,work5,work6
     +       ,even,count_nza
      real*8 a(1),b(1),fx(1),st(mx_st),t_t,t0,tol
      real*4 second
      equivalence (st(1),maxju),(st(2),maxu),(st(3),maxneed),(st(4),nza)
c     this is to avoid that the value of ratio is lost
c     with dynamic allocation memory
      save ratio
      data st/mx_st*0./,first/0/
c..   tolerance for accuracy in solution (in options numfac and check)
!     +    ,tol/1.d-8/
     +    ,tol/1.d-20/
!AL     

c function to ensure that real*8 variables are assigned an even address.
      even(i)=2*(i/2)+1

      flag=0
      t_t=second()
c     -----------------------------------------------
c     first time it checks whether ratio equal 1 or 2
c     -----------------------------------------------
      if (first.eq.0) then
         ratio=fratio()
         first=1
      endif
c     -------------------------------------------
c     no. of off-diagonal nze in upper triangular
c     -------------------------------------------
c ----comment either line; first one is faster but takes more memory
      nza=(ia(n+1)-1)
c      nza=count_nza(n,ia,ja,iout) 
c     --------------------
c     no. of nze in adjncy
c     --------------------
      nzad=2*nza
c     ----------------
c     allocate storage
c     ----------------
      perm   = 1
      invp   = perm + n
c     -------------------------
      if (option.eq.order) then
c     -------------------------
         xadj   = invp   + n
         adjncy = xadj   + n+1
         work1  = adjncy + nzad
         work2  = work1 + n
         work3  = work2 + n
         work4  = work3 + n
         needed = work4 + n
         call checkst (option,needed,maxneed,available,iout,flag)
         if (flag.ne.0) return
c        --------------------------------------------------
c        max integer (used as flag in the ordering routine)
c        --------------------------------------------------
         maxint=9999999
c        ---------------------------------------------------------
c        parameter for type of ordering -1=MD, 0=MMD, >0=slack MMD
c        ---------------------------------------------------------
         delta=0
         call unfold
     +               (n,ia,ja,is(xadj),is(adjncy),is(work1))
         t0=second()
c        -------------------------------
c        if mmd ordering routines by Liu
c        -------------------------------
         call genmmd
     +             (n,is(xadj),is(adjncy),is(invp),is(perm),delta,
     +              is(work1),is(work2),is(work3),is(work4),maxint,
     +              maxju)
c        ---------------------------
c        elseif md sparspak routines
c        ---------------------------
c        work5 = work4 + n
c        work6 = work5 + n
c        needed= work6 + n
c        call checkst (option,needed,maxneed,available,iout,flag)
c        if(flag.ne.0) return
c        call genqmd
c    +               (n,is(xadj),is(adjncy),is(perm),is(invp),
c    +                is(work1),is(work2),is(work3),is(work4),
c    +                is(work5),is(work6),maxju)
         st(21)=st(21)+second()-t0
c        ------------------------------
c        write down permutation vectors
c        ------------------------------
         rewind ioor
         write(ioor,*) n,maxju
         write(ioor,'(10i7)') (is(i),i=perm,perm+n-1),
     +                        (is(i),i=invp,invp+n-1)
c     -------------------------------
      elseif (option.eq.restart) then
c     -------------------------------
         rewind ioor
         i=0
         read(ioor,*,end=1) i,maxju
1        if(i.eq.n) then
            read(ioor,'(10i7)',end=2) (is(i),i=perm,perm+n-1),
     +                                (is(i),i=invp,invp+n-1)
2           inquire(ioor,iostat=ios)
            if (ios.ne.0) then
               write(iout,*) 'error in unit IOOR '
               flag=-2
            endif
         else
            if (i.ne.0) then
            write(iout,*) 'incorrect size of permutation vector'
            write(iout,*) 'size ',i,' found; size ',n, ' required'
            endif
            flag=-1
         endif
c     ------------------------------
      elseif (option.eq.symfac) then
c     ------------------------------
         iu     = invp   + n
         iju    = iu     + n+1
         ju     = iju    + n+1
         xadj   = ju     + maxju
         adjncy = xadj   + n+1
         work1  = adjncy + nzad
         work2  = work1  + n
         work3  = work2  + n
         needed = work3  + n
         call checkst (option,needed,maxneed,available,iout,flag)
         if (flag.ne.0) return
         t0=second()
         call unfold
     +               (n,ia,ja,is(xadj),is(adjncy),is(work1))
         call smbfct
     +               (n,is(xadj),is(adjncy),is(perm),is(invp),is(iu),
     +                maxu,is(iju),is(ju),maxju,is(work1),is(work2),
     +                is(work3),flag)
         st(22)=st(22)+second()-t0
      else
         iu        = invp  + n
         iju       = iu    + n+1
         ju        = iju   + n+1
         d         = even(ju    + maxju)
         u         = even(d     + n*ratio)
c        --------------------------
         if (option.eq.numfac) then
c        --------------------------
            iap    = u     + maxu*ratio
            jap    = iap   + n+1
            ap     = even(jap   + nza)
            work1  = even(ap    + nza*ratio)
            needed = work1 + n*ratio
            call checkst (option,needed,maxneed,available,iout,flag)
            if (flag.ne.0) return
            t0=second()
c           ---------------------------------------
c           first load numerical values of a into u
c           ---------------------------------------
            call loadap
     +                (n,ia,ja,a,is(iap),is(jap),is(ap),is(d),is(invp),
     +                 is(work1))
            call loadu
     +                (n,is(iap),is(jap),is(ap),is(iu),is(ju),is(iju),
     +                 is(u),is(work1))
            work1  = u     + maxu*ratio
            work2  = work1 + n
            work3  = even(work2 + n)
            needed = work3 + n*ratio
            call checkst (option,needed,maxneed,available,iout,flag)
            if (flag.ne.0) return
c           -----------------------
c           numerical factorization (slightly modified sparspak routine)
c           ----------------------------------
            call fgsfct
     +                 (n,is(iu),is(u),is(iju),is(ju),is(d),is(work1),
     +                  is(work2),is(work3),flag,tol,irank)
            call merlin (n,is(iu),is(u),is(d),tol)
            st(24)=st(24)+second()-t0
            st(11)=st(11)+1
c           ----------
c           zero pivot
c           ----------
            if(flag.ne.0) flag=numfac*1000000 + flag
c        -----------------------------
         elseif (option.eq.solve) then
c        -----------------------------
            work1  = even(u     + maxu*ratio)
            needed = work1 + n*ratio
            t0=second()
            call permute (n,b,is(work1),is(invp))
            call fullfb
     +                   (n,is(iu),is(ju),is(iju),is(u),is(d),is(work1))
            call permute (n,is(work1),b,is(perm))
            st(25)=st(25)+second()-t0
            st(12)=st(12)+1
c        --------------------------------------------------
         elseif (option.ge.sinit.and.option.le.ssolve) then
c        --------------------------------------------------
            path  = u     + maxu*ratio
            found = path  + n
            work1 = even(found + n)
            fstart= work1 + n*ratio
            fpnode= fstart+ fnseg+1
            needed= fpnode+ fnseg
            call checkst
c                        ------------------------------------------
c                        allows a safety margin for fpath and bpath
c                        ------------------------------------------
     +                   (option,needed+200,maxneed,available,iout,
     +                    flag)
            if(flag.ne.0) return
c           -------------------------
            if (option.eq.sinit) then
c           -------------------------
c              initialize and computes factorization table
c              -------------------------------------------
               call rzero (n,is(work1))
               call lzero (n,is(found))
               call makepath (n,is(path),is(iu),is(ju),is(iju))
c           ------------------------------
            elseif (option.eq.ssolve) then
c           ------------------------------
               t0=second()
               fpath = fpnode + fnseg
c              ------------------------------------------------------
c              obtains forward path (in the permuted vector of nodes)
c              ------------------------------------------------------
               call spermute (fnseg,fnode,is(fpnode),is(invp))
               call segpath
     +                      (n,is(fpnode),fnseg,is(fstart),is(fpath),
     +                       is(path),is(found))
               st(34)=st(34)+is(fstart+fnseg)-1
               if(feqb.eq.equal) then
                  needed= fpath + is(fstart+fnseg)-1
                  call checkst
     +                        (option,needed,maxneed,available,iout,
     +                         flag)
                  if (flag.ne.0) return
c                 -----------------------------
c                 initialize working rhs vector
c                 -----------------------------
                  call zerout (is(work1),is(fpath),is(fstart),fnseg)
c                 ----------------
c                 input rhs values
c                 ----------------
                  call expand (fnseg,is(fpnode),fx,is(work1))
c                 -----
c                 solve
c                 -----
                  call fastfb
     +                        (is(fpath),is(fstart),fnseg,is(fpath),
     +                         is(fstart),fnseg,is(iu),is(ju),is(iju),
     +                         is(u),is(d),is(work1))
c                 ---------------------
c                 stores solutions in b
c                 ---------------------
                  call compress (fnseg,is(fpnode),b,is(work1))
                  st(31)= max0(int(st(31)),is(fstart+fnseg)-1)
                  st(35)=st(35)+is(fstart+fnseg)-1
               else
                  bstart = fpath + is(fstart+fnseg)-1
                  bpnode = bstart+ bnseg+1
                  bpath  = bpnode+ bnseg
c                 ----------------------------------
c                 obtains backward path if necessary
c                 ----------------------------------
                  call spermute (bnseg,bnode,is(bpnode),is(invp))
                  call segpath
     +                         (n,is(bpnode),bnseg,is(bstart),
     +                          is(bpath),is(path),is(found))
                  needed= bpath + is(bstart+bnseg)-1
                  call checkst
     +                        (option,needed,maxneed,available,iout,
     +                         flag)
                  if (flag.ne.0) return
                  call zerout (is(work1),is(bpath),is(bstart),bnseg)
                  call zerout (is(work1),is(fpath),is(fstart),fnseg)
                  call expand (fnseg,is(fpnode),fx,is(work1))
                  call fastfb
     +                        (is(fpath),is(fstart),fnseg,is(bpath),
     +                         is(bstart),bnseg,is(iu),is(ju),is(iju),
     +                         is(u),is(d),is(work1))
                  call compress (bnseg,is(bpnode),b,is(work1))
                  st(31)=max(int(st(31)),is(bstart+bnseg)-1)
                  st(35)=st(35)+is(bstart+bnseg)-1
               endif
               st(13)=st(13)+1
               st(27)=st(27)+second()-t0
               st(30)=max(int(st(30)),fnseg)
               st(31)=max(int(st(31)),bnseg)
               st(32)=st(32)+fnseg
               st(33)=st(33)+bnseg
            endif
c        ---------------------------
         elseif (option.eq.det) then
c        ---------------------------
            call mkdet (n,is(d),b,irank,tol,iout,flag)
            st(14)=st(14)+1
c        ----------------------------
         elseif (option.eq.ldet) then
c        ----------------------------
            call mkldet (n,is(d),b,irank,tol,iout,flag)
            st(14)=st(14)+1
c        ----------------------------
         elseif (option.eq.check) then
c        ----------------------------
            work1  = u     + maxu*ratio
            call checkb 
     +                  (n,ia,ja,a,is(d),b,fx,tol,is(work1),is(invp)
     +                  ,iout,flag)
c        ----------------------------
         elseif (option.eq.spin) then
c        ----------------------------
            zja    = u     + maxu*ratio
            za     = even(zja   + maxu+n)
            zd     = even(za    + (maxu+n)*ratio)
            utia   = zd    + n*ratio
            utja   = utia  + n+1
            work1  = utja  + maxu+n
            work2  = work1 + n
            work3  = even(work2 + n)
            needed = work3 + n*ratio
            call checkst (option,needed,maxneed,available,iout,flag)
            if (flag.ne.0) return
            t0=second()
            call fsinverse
     +                    (n,is(perm),is(invp),is(d),is(iju),is(ju),
     +                     is(iu),is(u),is(work3),is(za),is(zd),
     +                     is(zja),is(utia),is(utja),
     +                     is(work1),is(work2),ia,ja,a)
            st(26)=st(26)+second()-t0
            st(15)=st(15)+1
c        ----------------------------
         elseif (option.eq.spin1) then
c        ----------------------------
            work1  = even(u     + maxu*ratio)
            work2  = even(work1 + n*ratio)
            work5  = even(work2 + n*ratio)
            needed = work5 + nza+n
            call checkst (option,needed,maxneed,available,iout,flag)
            if (flag.ne.0) return
            t0=second()
            call fsinverse1
     +                    (n,is(perm),is(invp),is(d),is(iju),is(ju),
     +                     is(iu),is(u),is(work1),is(work2),is(work1),
     +                     is(work2),is(work5),ia,ja,a)
            st(26)=st(26)+second()-t0
            st(15)=st(15)+1
c        ----------------------------
         elseif (option.eq.stat) then
c        ----------------------------
            call figureout (n,available,iout,st,irank)
c        ----
         else
c        ----
            write(iout,*) 'option no. ',option,' not implemented'
            flag=-99
         endif
      endif
      st(20)=st(20)+second()-t_t
      return
      end
