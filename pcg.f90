module pcg
! Solves the system of equations A*sol=rhs by preconditioned gradient algorithm. 
!  The default preconditioner is diagonal but can be blockdiagonal for
!  blksize>1. The block implementation uses the sparse_hashm structure and
!  fspak90, and is designed primarly for testing and not efficiency. More
!  efficient implementation would use an array of dense symmetric matrices. 

   use denseop; use sparseop
   implicit none

    interface solve_pcg
       module procedure solve_densem_pcg,&
                        solve_sparse_hashm_pcg!,&
!                        solve_sparse_ija_pcg
    end interface
    
   contains
   
   
   
   subroutine solve_densem_pcg(A,rhs,sol,blksize)
   ! solve A sol = rhs by preconditioned conjugate gradient for densem A
   !   Preconditioner can be block-diagonal if blksize>1
   type (densem)::A
   real (rh)::rhs(:),sol(:)
   integer,optional::blksize
   real (r8),allocatable::m(:),r(:),p(:),z(:),w(:)
   real (r8)::alpha,beta,tau,oldtau,conv,val
   integer::i,j,k,block
   
   !variables for blockdiagonal if any
   type (sparse_hashm) :: m_hash
   type (sparse_ija) :: m_ija
   
   if (present(blksize)) then
         block=blksize
      else
         block=1
    endif
    
    if (block == 1) then
        allocate(m(a%n))
      else
        call zerom(m_hash,a%n) !,int(a%n*block*(block+1)/2.0*1.3))
   endif       
   
   ! find the preconditioner 
   if (block == 1 ) then 
        !diagonal preconditioner
        m=0
        do i=1,a%n
!           if (a%x(i,i)/=0) m(i)=1/a%x(i,i)
            val=getm(i,i,a)
            if (val /= 0) m(i)=1/val
        enddo
     else   
        ! block preconditioner; not efficiently implemented but simple
        do i=1,a%n,block
           do j=0,block-1
              do k=0,block-1
                 call addm(getm(i+j,i+k,a),i+j,i+k,m_hash)
              enddo
           enddo
        enddo
        m_ija=m_hash
        call reset(m_hash)
        call fspak90('invert',m_ija)
        call fspak90('reset',m_ija)
     endif
   
   allocate (r(a%n),p(a%n),z(a%n),w(a%n))
   sol=0
   r=rhs
  
   do k=1,default_rounds
      if (block==1) then
            z=m*r
         else
            z=mult(m_ija,r)
      endif
      tau=dot_product(z,r)
      if (k == 1) then
         beta=0; p=z
       else
         beta=tau/oldtau
         p=z+beta*p
      end if
      w=mult(a,p)
      alpha=tau/dot_product(p,w)
      sol=sol+alpha*p
      if (mod(k,100) /= 0) then
           r=r-alpha*w
         else
           r=rhs-mult(a,real(sol,r8))
      endif
      conv=dot_product(r,r)/dot_product(rhs,rhs)
      print*,'round ',k,'   convergence=',conv
      if (conv < default_conv) exit
      oldtau=tau
   enddo
   
   if (block == 1) then
         deallocate(m)
      else
         call reset(m_ija)
   endif
   deallocate (r,p,z,w)
   print*,k,' iterations,   convergence criterion=',conv

  end subroutine
  
     
   
   subroutine solve_sparse_hashm_pcg(A,rhs,sol,blksize)
   ! solve A sol = rhs by preconditioned conjugate gradient for sparse_hashm A
   !   Preconditioner can be block-diagonal if blksize>1
   type (sparse_hashm)::A
   real (rh)::rhs(:),sol(:)
   integer,optional::blksize
   real (rh),allocatable::m(:),r(:),p(:),z(:),w(:)
   real (r8)::alpha,beta,tau,oldtau,conv,val
   integer::i,j,k,block
   
   !variables for blockdiagonal if any
   type (sparse_hashm) :: m_hash
   type (sparse_ija) :: m_ija
   
   if (present(blksize)) then
         block=blksize
      else
         block=1
    endif
    
    if (block == 1) then
        allocate(m(a%n))
      else
        call zerom(m_hash,a%n) !,int(a%n*block*(block+1)/2.0*1.3))
   endif       
   
   ! find the preconditioner 
   if (block == 1 ) then 
        !diagonal preconditioner
        m=0
        do i=1,a%n
!           if (a%x(i,i)/=0) m(i)=1/a%x(i,i)
            val=getm(i,i,a)
            if (val /= 0) m(i)=1/val
        enddo
     else   
        ! block preconditioner; not efficiently implemented but simple
        do i=1,a%n,block
           do j=0,block-1
              do k=0,block-1
                 call addm(getm(i+j,i+k,a),i+j,i+k,m_hash)
              enddo
           enddo
        enddo
        m_ija=m_hash
        call reset(m_hash)
        call fspak90('invert',m_ija)
        call fspak90('reset',m_ija)
     endif
   
   allocate (r(a%n),p(a%n),z(a%n),w(a%n))
   sol=0
   r=rhs
  
   do k=1,default_rounds
      if (block==1) then
            z=m*r
         else
            z=mult(m_ija,real(r,r8))
      endif
      tau=dot_product(z,r)
      if (k == 1) then
         beta=0; p=z
       else
         beta=tau/oldtau
         p=z+beta*p
      end if
      w=mult(a,p)
      alpha=tau/dot_product(p,w)
      sol=sol+alpha*p
      if (mod(k,100) /= 0) then
           r=r-alpha*w
         else
           r=rhs-mult(a,sol)
      endif
      conv=dot_product(r,r)/dot_product(rhs,rhs)
      print*,'round ',k,'   convergence=',conv
      if (conv < default_conv) exit
      oldtau=tau
   enddo
   
   if (block == 1) then
         deallocate(m)
      else
         call reset(m_ija)
   endif
   deallocate (r,p,z,w)
   print*,k,' iterations,   convergence criterion=',conv

  end subroutine
  

  end module
 
