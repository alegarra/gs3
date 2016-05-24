! Collection of predominantly sparse matrix structures and operations

module sparsem
use kinds
implicit none

integer::sparsem_msg=2	! message level. 0=no messages, 3=max messages
			! not yet implemented

!=====================================================================
! Type definitions
!=====================================================================

type densem     !traditional dense square matrix
   integer :: n
   real (r8) ,pointer::x(:,:)=>null()
end type densem


type dense_symm       !upper stored symmetric dense matrix
   integer ::n
   real (r8) ,pointer::x(:)=>null()
end type dense_symm


type sparse_hashm
    integer:: n,&       ! for compatibility mainly
            nel,&       ! number of elements
            filled,&    ! number of filled elements
            status      ! 1 if ready to hash, 2 if in sorted order
     real (rh) , pointer :: x(:,:)=>null()
end type sparse_hashm


type sparse_ija
     integer :: n,nel      !n is number of equations; nel is number of nonzeroes
     integer, pointer::ia(:)=>null()      !will be ia(n+1)
     integer, pointer::ja(:)=>null()      !will be ja(m)
     real (r8), pointer::a(:)=>null()     !will be a(m)
end type



type ll_element  ! basic element of linked-list
     integer::col
     real::val
     type (ll_element),pointer::next=>null()
end type    

type sparse_ll	! actual linked-list structure
    integer n
    type (ll_element),pointer::ia(:)=>null()
end type

!------------------------------------------
! Defaults and constants for subroutines
!------------------------------------------
 integer        ::     default_hash_size=5000,&
                       default_rounds = 1000  !maximum number of rounds
 real           ::     default_conv=1e-10,&   ! stopping criterion
                       default_relax=1.0,&      !relaxation factor
                       maxfill=.9,&          ! max fill of hash table
                       hash_incr=1.5          ! increase of hash table
 logical        ::     default_zerosol=.true. ! zero solutions before 
                                              ! iterative solving	       
					      
					      					      
 integer,parameter ::  conv_straight=0,&      ! straight conversion from hash
                                              !  to ija formats
                      &conv_upper_full=1      !	conversion from upper to full
		                              ! storage				      
 


!=====================================================================
!Interfaces
!=====================================================================

  interface init
      module procedure init_densem,&
                       init_dense_symm,&
                       init_sparse_hashm,&
                       init_sparse_ija
  end interface

  interface zerom
      module procedure zero_densem,&
                       zero_dense_symm,&
                       zero_sparse_hashm,&
                       zero_sparse_ija
  end interface

  interface reset
      module procedure deall_densem,&
                       deall_dense_symm,&
                       deall_sparse_hashm,&
                       deall_sparse_ija
  end interface

  interface addm
     module procedure add_densem, &
                      add_dense_symm,&
                      add_sparse_hashm
  end interface

  interface setm
     module procedure set_densem,&
                      set_dense_symm,&
                      set_sparse_hashm
  end interface

  interface getm
     module procedure get_densem,&
                      get_dense_symm,&
                      get_sparse_hashm,&
                      get_sparse_ija
  end interface


  interface solve_iterm
     module procedure solve_dense_iter,&
                      solve_dense_symm_iter,&
                      solve_sparse_hashm,&
                      solve_sparse_ija
  end interface

  interface printm
     module procedure print_densem,&
                      print_dense_symm,&
                      print_sparse_hashm,&
                      print_x,&
                      print_sparse_ija
  end interface

  interface   block
     module PROCEDURE block_densem,&
                      block_hashm,&
                      block_ija
  end interface

  interface trace
     module procedure trace_densem,&
  		      trace_hashm,&
  		      trace_ija
  end interface

  interface traceblock
     module procedure traceblock_densem,&
  		      traceblock_hashm,&
		      traceblock_ija
  end interface

  interface quadrf
     module procedure quadr_densem,&
                      quadr_hashm,&
                      quadr_ija
  end interface

  interface assignment (=)
      module procedure convert_half_dense,&
                       convert_hash_densem,&
                       convert_hash_ija,&
                       convert_ija_densem,&
                       convert_ija_hash,&
                       convert_densem_hashm,&
                       convert_densem_ija,&
		       copy_dense,&
		       copy_symm,&
		       copy_hashm,&
		       copy_ija
		       
  end interface

  interface	!external function
      function hashvr(dat,r,ind,m,n,mode,nr)
      use kinds
      integer ::hashvr,r,mode,nr,m,n
      real (rh)::ind(:,:),dat(r)
      end function
  end interface


  contains

!=====================================================================
!Module subroutines
!=====================================================================

!----------------------------------------------------------------
!  initialize; 
!  necessary on systems where pointers not initialized at startup
!----------------------------------------------------------------
subroutine init_densem(x)
 ! Initializes a dense matrix
 type (densem) :: x
!
 x%n=0
 nullify(x%x)
end subroutine


subroutine init_dense_symm(x)
 ! Initializes dense symmetric matrix
 type (dense_symm) :: x
!
 x%n=0
 nullify(x%x)
end subroutine

subroutine init_sparse_hashm(x)
 ! Initializes sparse matrix in hash form
 type (sparse_hashm) :: x
!
 x%n=0; x%filled=0; x%status=0
 nullify(x%x)
end subroutine

subroutine init_sparse_ija(x)
 ! Initializes sparse matrix in ija form
 type (sparse_ija) :: x
!
 x%n=0; x%nel=0
 nullify(x%ia,x%ja,x%a)
end subroutine

subroutine init_ll(x)
 ! Initializes sparse matrix in linked-list form
 type (sparse_ll) :: x
!
 x%n=0; nullify(x%ia)
end subroutine

!--------------
!  zero
!--------------
subroutine zero_densem(x,n)
 ! Allocates and zero a dense matrix of n x n
 integer :: n
 type (densem) :: x
!
 if (associated(x%x)) then
       if (n /= x%n) then
          deallocate(x%x)   ! different dimension, deallocate and allocate
          allocate(x%x(n,n))
       endif
    else
       allocate(x%x(n,n))   
 endif
 x%n=n
 x%x=0
end subroutine


subroutine zero_dense_symm(x,n)
 ! Allocates and zero symmetric upper-stored dense matrix of n x n
 integer :: n
 type (dense_symm) :: x
!
 if (associated(x%x)) then
       if (n /= x%n) then
          deallocate(x%x)   ! different dimension, deallocate and allocate
          allocate(x%x((n*(n+1))/2))
       endif
    else
       allocate(x%x((n*(n+1))/2))
 endif
 x%n=n
 x%x=0
end subroutine

subroutine zero_sparse_hashm(x,n,max_elem)
 ! Allocates and zero sparse matrix in hash form
 integer :: n
 integer, optional :: max_elem
 type (sparse_hashm) :: x
!
!
 if (present(max_elem)) then
        x%nel=max_elem
    else
           ! for matrix used previously apply old dimensions
        if (.not. associated(x%x)) x%nel=default_hash_size
 endif
!
 if (associated(x%x)) then
       if (n /= x%n) then
          deallocate(x%x)   ! different dimension, deallocate and allocate
          allocate(x%x(3,x%nel))
       endif
    else
       allocate(x%x(3,x%nel))
 endif
 x%n=n; x%filled=0; x%status=1
 x%x=0
end subroutine

subroutine zero_sparse_ija(x,n)
 ! Allocates and zero sparse matrix in ija form
 integer :: n
 type (sparse_ija) :: x
!
! Actual allocation for the IJA format can be done through conversion
! from other programs, so here only deallocate if allocated

 if (associated(x%ia)) then
      call reset(x)
 endif
end subroutine


subroutine zero_ll(x,n)
 ! Allocates and zero a linked-list matrix of n x n
 integer :: n,i
 type (sparse_ll) :: x
 type(ll_element),pointer::current
!
 if (associated(x%ia)) then
       
       if (n /= x%n) then
            ! different dimension, deallocate and allocate
            ! call reset(x)
             print*,'reset in sparse_ll not yet implemented'
             stop
          else
           ! zero all elements but do not deallocate
             do i=1,n
                current=>x%ia(i)
                do while (associated(current))
                   current%val=0
                   current=>current%next
                enddo  
             enddo
       endif
    else
       ! new matrix
       allocate(x%ia(n))
       nullify(x%ia)
 endif
 x%n=n
end subroutine


!----------------------------
!  reset (deallocate) matrix
!----------------------------

subroutine deall_densem(x)
 ! Deallocates  dense matrix
 type (densem) :: x
!
 x%n=0
 if (associated(x%x)) deallocate(x%x)
end subroutine

subroutine deall_dense_symm(x)
 ! deallocates  symmetric upper-stored dense matrix 
 type (dense_symm) :: x
!
 x%n=0
 if (associated(x%x)) deallocate(x%x)
end subroutine

subroutine deall_sparse_hashm(x)
 ! deallocates sparse hash matrix
 type (sparse_hashm) :: x
!
 x%n=0; x%filled=0; x%status=0
 x%nel=0
 if (associated(x%x)) deallocate(x%x)
end subroutine

subroutine deall_sparse_ija(x)
 ! deallocates sparse ija matrix
 type (sparse_ija) :: x
!
 x%n=0; x%nel=0
 if (associated(x%ia))deallocate(x%ia,x%ja,x%a)
end subroutine

!--------------
!  add to matrix
!--------------

subroutine add_densem(a,i,j,x)
! adds a(i,j) to dense matrix x
 real (rh) :: a
 type (densem):: x
 integer i,j

 x%x(i,j)=x%x(i,j)+a
end subroutine


subroutine add_dense_symm(a,i,j,x)
! adds a(i,j) to dense symmetric upper-stored matrix x
 real (rh) :: a
 type (dense_symm):: x
 integer i,j,k
 if (j >=i) then
    k=pos_symm(i,j,x%n)
    x%x(k)=x%x(k)+a
 endif
end subroutine


 function pos_symm(i,j,n) result (address)
 !finds position of (i,j) element of a n.n matrix upper-stored
 integer :: i,j,n,address
 
 if (j >= i) then
      address=(i-1)*n-(i*(i-3))/2+j-i
   else
      address=(j-1)*n-(j*(j-3))/2+i-j
 endif
 end function


recursive subroutine add_sparse_hashm(a,i,j,x,storage_type)
! adds a(i,j) to sparse hash upper-stored matrix x
! optional storage_type is 'u' for upper-trangular (deafult) and 'f' for full
 
 real (rh):: a
 type (sparse_hashm):: x,y
 integer i,j,k,full
 character,optional::storage_type
 logical::upper_storage=.true.
 
 if (present(storage_type)) then
     if (storage_type == 'f') upper_storage=.false.
 endif
     
 if (x%status /=1) then
    print*,'Structure sparse_hashm not ready; probably destroyed by solving'
    stop
 endif
 if (.not. upper_storage .or. j >=i) then
    full=hashvr((/real(i,rh),real(j,rh)/),2,x%x,x%nel,3,1,x%filled)
    if (full == -1 .or. float(x%filled)/x%nel > maxfill ) then
!     Matrix full, copy to a matrix hash_incr times larger
       call init(y)
       call zerom(y,x%n,int(hash_incr*x%nel))
       do k=1,x%nel
          if (x%x(1,k) /= 0) then
              call addm(x%x(3,k),int(x%x(1,k)),int(x%x(2,k)),y)
          endif
       enddo
       print*,'hash matrix increased from ',x%nel,' to ',int(hash_incr*x%nel)
       ! Move y to x
       call reset(x)
       x%n=y%n; x%nel=y%nel;x%filled=y%filled;x%status=y%status;x%x=>y%x;
       nullify(y%x)
       full=hashvr((/real(i,rh),real(j,rh)/),2,x%x,x%nel,3,1,x%filled)
    endif
    if (full >0) then
       x%x(3,full)=x%x(3,full)+a
       else
         call add_sparse_hashm(a,i,j,x)
    endif	 
	 
 endif
end subroutine

!--------------
!  set matrix element
!--------------

subroutine set_densem(a,i,j,x)
! x(i,j)=a for dense square matrix
 real (rh) :: a
 type (densem):: x
 integer i,j

 x%x(i,j)=a
end subroutine


subroutine set_dense_symm(a,i,j,x)
! x(i,j)=a for dense symmetric upper-stored matrix x
 real (rh) :: a
 type (dense_symm):: x
 integer i,j,k
 if (j >=i) then
    k=pos_symm(i,j,x%n)
    x%x(k)=a
 endif
end subroutine


subroutine set_sparse_hashm(a,i,j,x)
! x(i,j)=a for sparse hash upper-stored matrix x

 real (rh) :: a
 type (sparse_hashm):: x
 integer i,j,full
!
 if (x%status /=1) then
    print*,'Structure sparse_hashm not ready; probably sorted'
    stop
 endif
 if (j >=i) then
    full=hashvr((/real(i,rh),real(j,rh)/),2,x%x,x%nel,3,1,x%filled)
    if (full == -1) then
       print*,'hash matrix too small'
       stop
    endif
    x%x(3,full)=a
 endif
end subroutine


!-------------------------
!  get scalar from matrix
!-------------------------

function get_densem(i,j,x) RESULT (a)
! a=x(i,j) for dense matrix x
 real (rh) :: a
 type (densem):: x
 integer i,j

 a=x%x(i,j)
end function


function get_dense_symm(i,j,x) RESULT(a)
! a=x(i,j) for dense symmetric upper-stored matrix x
 real (rh) :: a
 type (dense_symm):: x
 integer i,j,k
 if (j >=i) then
    k=pos_symm(i,j,x%n)
    a=x%x(k)
 endif
end function


function get_sparse_hashm(i,j,x) RESULT(a)
! a=x(i,j) for sparse hash matrix

 real (rh) :: a
 type (sparse_hashm):: x
 integer i,j,full
!
 if (x%status /=1) then
    print*,'Structure sparse_hashm not ready; probably sorted'
    stop
 endif
 if (j >=i) then
    full=hashvr((/real(i,rh),real(j,rh)/),2,x%x,x%nel,3,0,x%filled)
    if (full <1) then
       a=0
     else
       a=x%x(3,full)
    endif
  else
    a=0
 endif
end function

function get_sparse_ija(i,j,x) RESULT(a)
! a=x(i,j) for sparse IJA matrix
 real (r8) :: a
 type (sparse_ija) :: x
 INTEGER,INTENT(IN) :: i,j
 integer:: k

 do k=x%ia(i),x%ia(i+1)-1
   if (x%ja(k) ==j) then
       a=x%a(k)
       return
   end if
 enddo
 a=0    ! element not found
END function
!----------------------------------
!  solve linear system of equations
!----------------------------------

subroutine solve_dense_iter(x,rs,sol)
! solves x * sol = rs for dense matrix x
 type (densem) :: x
 real (rh) :: rs(:),sol(:)

 if (default_zerosol) sol=0

 call solve_dense_gs(x%n,x%x,rs,sol)
end subroutine

subroutine solve_dense_symm_iter(x,rs,sol)
! solves x * sol = rs for dense symmetric upper-stored matrix x
 type (dense_symm) :: x
 real (rh) :: rs(:),sol(:)

if (default_zerosol) sol=0
   
 call solve_dense_symm_gs(x%n,x%x,rs,sol)
end subroutine

subroutine solve_sparse_hashm(x,rs,sol)
! solves x * sol = rs for sparse-hash matrix x assumed upper-stored 

 interface
      subroutine sortqr(a,n,m,ix,iy,col,icol)
        use kinds
        integer :: n,m,ix,iy,icol
        real (rh)::a(:,:)
        integer :: col(:)
      end subroutine
      subroutine clusterr(x,m1,m2,ix)
        use kinds
        integer m1,m2,ix
        real (rh)::x(:,:)
      end subroutine
 end interface

 type (sparse_hashm) :: x
 real (rh) :: rs(:),sol(:)
 
 if (default_zerosol) sol=0
 
 select case (x%status)
     case (1)           !hash ready, needs sorting
        call clusterr(x%x,x%nel,3,x%filled)
        call sortqr(x%x,x%nel,3,x%filled,3,(/1,2/),2)
        x%status=2
     case (2)           ! already sorted, no action
        continue
     case default
        print*,'Unknown status of hash matrix'
 end select

 call solve_sparse_hashm_gs(x%n,x%x,rs,sol,x%filled)
end subroutine

subroutine solve_sparse_ija(x,rs,sol)
! solves x * sol = rs for sparse-ija matrix x assumed upper-stored
 type (sparse_ija) :: x
 real (rh) :: rs(:),sol(:)

if (default_zerosol) sol=0   

 call solve_sparse_ija_gs(x%n,x%ia,x%ja,x%a,rs,sol,x%nel)
end subroutine

subroutine default_iter(conv,maxround,relax,zerosol)
! changes default parameters for convergence criterion,
! maximum number of rounds, relaxation factor, and zeroing solutions.
 integer, optional:: maxround
 real,optional::conv,relax
 logical,optional::zerosol

 if (present(conv)) default_conv=conv
 if (present(maxround)) default_rounds=maxround
 if (present(relax)) default_relax=relax
 if (present(zerosol)) default_zerosol=zerosol
end subroutine

!--------------
!  print matrix
!---------------------------------------------------------------------------
subroutine print_x(x,mode)
! prints regular real matrix x by rows
 real x(:,:)
 integer i
 character (*),optional:: mode

 do i=1,size(x,dim=2)
    print '(20f6.2)',x(:,i)
 enddo
end subroutine

subroutine print_densem(x,mode)
! prints dense matrix
 type (densem) :: x
 integer :: i
 character (*),optional:: mode

 do i=1,x%n
   print '(11F7.2)',x%x(i,1:x%n)
 enddo
end subroutine

subroutine print_dense_symm(x,mode)
! prints dense symmetric upper-stored matrix
 type (dense_symm) :: x
 type (densem) :: y
 character (*),optional:: mode

 y=x
 call printm(y)
end subroutine

subroutine print_sparse_hashm(x,mode)
! prints sparse hash matrix, assumed upper-stored
! mode='internal' prints internal representation
 type (sparse_hashm) :: x
 type (densem) :: y
 character (*),optional:: mode
 integer ::i

 if (present(mode)) then
    if (mode(1:3)=='int') then
       do i=1,x%nel
          if (x%x(1,i) /= 0) then
             print '(2i8,f12.4)',int(x%x(1,i)),int(x%x(2,i)),x%x(3,i)
          endif
       enddo
       return
    endif
 endif
 y=x
 call printm(y)
end subroutine


subroutine print_sparse_ija(x,mode)
! prints sparse ija matrix, assumed upper-stored
! mode='internal' prints internal representation
 type (sparse_ija) :: x
 type (densem) :: y
 character (*),optional:: mode
 integer ::i,j

 if (present(mode)) then
    if (mode(1:3)=='int') then
       do i=1,x%n
          do j=x%ia(i),x%ia(i+1)-1
             print '(2i8,f12.4)',i,x%ja(j),x%a(j)
          enddo
       enddo
       return
    endif
 endif
 y=x
 call printm(y)
end subroutine

!--------------
! block
!--------------

function block_densem(x,i1,i2,j1,j2,step) RESULT (y)
! y=x(i1:i2:step,j1:j2:step) for dense matrices
  type (densem):: x,y
  INTEGER :: i1,i2,j1,j2
  integer,optional::step
  integer :: st

  if (present(step)) then
      st=step
     else
      st=1
  endif

 !
  if (i2-i1 /= j2-j1) then
     PRINT*,'BLOCK_DENSEM: matrix not square'
     stop
  end if

  call zerom(y,(i2-i1)/st+1)
  y%x=x%x(i1:i2:st,j1:j2:st)
END function

function block_hashm(x,i1,i2,j1,j2,step) RESULT (y)
! y=x(i1:i2:step,j1:j2:step) for hash matrices
  type (sparse_hashm):: x,y
  INTEGER,INTENT(IN) :: i1,i2,j1,j2
  INTEGER :: i,ix,jy
  integer,optional::step
  integer :: st

  if (present(step)) then
      st=step
     else
      st=1
  endif

 !
  if (i2-i1 /= j2-j1) then
     PRINT*,'BLOCK_HASHM: matrix not square'
     stop
  end if

  call zerom(y,(i2-i1)/st+1)
   do i=1,x%nel
          if (x%x(1,i) /= 0) then
             ix=x%x(1,i)
             jy=x%x(2,i)
             if (ix >=i1 .and. ix<=i2 .and. mod(ix-i1,st) == 0 .and.&
                 jy >=j1 .and. jy<=j2 .and. mod(jy-j1,st) == 0) then
		   call addm(x%x(3,i),(ix-i1)/st+1,(jy-j1)/st+1,y)
             endif
          endif
   enddo
END function

function block_ija(x,i1,i2,j1,j2,step) RESULT (y)
! y=x(i1:i2:step,j1:j2:step) for ija matrices
  type (sparse_ija) :: x,y
  INTEGER,INTENT(IN) :: i1,i2,j1,j2
  INTEGER ::i,j,k
  integer,allocatable::work(:)
  integer,optional::step
  integer :: st

  if (present(step)) then
      st=step
     else
      st=1
  endif
 !
  if (i2-i1 /= j2-j1) then
     PRINT*,'BLOCK_IJA: matrix not square'
     stop
  end if

  call zerom(y,(i2-i1)/st+1);
  allocate(work((i2-i1)/st+1)); work=0
 
  do i=i1,i2,st
      do j=x%ia(i),x%ia(i+1)-1
         if (x%ja(j) >=j1 .and. x%ja(j) <=j2 .and. mod(x%ja(j)-j1,st)==0) then
             work((i-i1)/st+1)=work((i-i1)/st+1)+1
         end if
      enddo
  enddo

  y%n=(i2-i1)/st+1; allocate(y%ia(y%n+1),y%ja(sum(work)),y%a(sum(work)))
  y%ia(1)=1; 
  do i=1,y%n
     y%ia(i+1)=y%ia(i)+work(i)
  enddo
  work=0

  do i=i1,i2,st
      do j=x%ia(i),x%ia(i+1)-1
         if (x%ja(j) >=j1 .and. x%ja(j) <=j2 .and. mod(x%ja(j)-j1,st)==0) then
             k=y%ia((i-i1)/st+1)+work((i-i1)/st+1)
             y%ja(k)=(x%ja(j)-j1)/st+1
             y%a(k)=x%a(j)
             work((i-i1)/st+1)=work((i-i1)/st+1)+1
         end if
      enddo
  enddo
  deallocate(work)
END function

!--------------
! trace
!--------------

function trace_densem(x,y) RESULT(tr)
! tr=trace(x,y) for dense matrices
  TYPE(densem):: x,y
  REAL (r8) :: tr
  !
  tr=0
  if (x%n /=y%n) then
     PRINT*,'TRACE_DENSEM: Matrices of different size'
     stop
  end if
  tr=SUM(x%x*y%x)
end function


function rel_diff(x,y)
  ! calculates relative difference of x and y
  real (r8) ::x,y,rel_diff

  if (max(abs(x),abs(y)) /=0) then
       rel_diff=abs(x-y)/max(abs(x),abs(y))
     else
       rel_diff=0
  endif
  end function


function trace_hashm(x,y) RESULT(tr)
! tr=trace(x,y) for dense matrices
  TYPE(sparse_hashm),target:: x,y
  TYPE(sparse_hashm),POINTER::mat1,mat2
  REAL (r8) :: tr,trlo,trup,trdiag,a
  INTEGER :: i,j,k
  !
  tr=0; trlo=0; trup=0; trdiag=0

  if (x%n /=y%n) then
     PRINT*,'TRACE_DENSEM: Matrices of different size'
     stop
  end if
  ! the first matrix should be smaller
  if (x%filled < y%filled) then
      mat1=>x; mat2=>y
     else
      mat1=>y; mat2=>x

  end if
  do i=1,mat1%nel
     if (mat1%x(1,i) /=0) then
        j= INT(mat1%x(1,i));  k=INT(mat1%x(2,i))
        a=mat1%x(3,i)* getm(j,k,mat2)
        if (j >  k) trlo=trlo+a
        if (j == k) trdiag=trdiag+a
        if (j < k)  trup=trup+a
     end if
  end do

! handle triangular matrices
  if (trlo == 0) then
       tr=2*trup+trdiag
     elseif (rel_diff(trlo,trup) < 1d-6 ) then
       tr=trlo+trdiag+trup
     else
       print*,'TRACE_HASHM inconsistent trace value: tr[lower,diag,upper]= ',&
              trlo,trdiag,trup
  endif
end function

function trace_ija(x,y) RESULT(tr)
! tr=trace(x,y) for ija matrices
  TYPE(sparse_ija):: x,y
  REAL (r8):: tr,trlo,trup,trdiag
  REAL(r8) :: work(x%n)
  integer :: i,j,k
  !
  if (x%n /=y%n) then
     PRINT*,'TRACE_IJA: Matrices of different size'
     stop
  end if
  !
  work=0; tr=0; trlo=0; trup=0; trdiag=0

  do i=1,x%n                    !work=x(i,:)
     do j=x%ia(i),x%ia(i+1)-1
        work(x%ja(j))=x%a(j)
     end do
     do j=y%ia(i),y%ia(i+1)-1    !tr=tr+work'*y(i,:)
        k=y%ja(j)
        if (i >  k) trlo=trlo+y%a(j)*work(k)
        if (i == k) trdiag=trdiag+y%a(j)*work(k)
        if (i < k)  trup=trup+y%a(j)*work(k)
        tr=tr+y%a(j)*work(k)
     end do
     do j=x%ia(i),x%ia(i+1)-1    ! work=0
        work(x%ja(j))=0
     end do
  end do
!
! handle triangular matrices
  if (trlo == 0) then
       tr=2*trup+trdiag
     elseif (rel_diff(trlo,trup) < 1d-6 ) then
       tr=trlo+trdiag+trup
     else
       print*,'TRACE_IJA: inconsistent trace value: tr[lower,diag,upper]= ',&
              trlo,trdiag,trup
  endif
end function


!--------------
! traceblock
!--------------

function traceblock_densem(x,y,i1,i2,j1,j2,step) RESULT (tr)
! y=trace(x,y(i1:i2:step,j1:j2:step)) for dense matrices
  type (densem):: x,y
  real (r8) :: tr
  INTEGER :: i1,i2,j1,j2
  integer,optional::step
  integer :: st

  if (present(step)) then
      st=step
     else
      st=1
  endif
 !
  if (i2-i1 /= j2-j1) then
     PRINT*,'BLOCKTRACE_DENSEM: ranges not equal'
     stop
  end if

  if (x%n /=(i2-i1)/st+1) then
     PRINT*,'TRACE_DENSEM: Matrices of different size'
     stop
  end if
 
  tr=SUM(x%x*y%x(i1:i2:st,j1:j2:st))
END function


function traceblock_hashm(x,y,i1,i2,j1,j2,step) RESULT (tr)
! y=trace(x,y(i1:i2:step,j1:j2:step)) for hash matrices
  type (sparse_hashm):: x,y
  INTEGER,INTENT(IN) :: i1,i2,j1,j2
  real (r8)::tr
  INTEGER :: i,j,k,k1,k2,l,l1,l2,st
  integer,optional::step

  if (present(step)) then
      st=step
     else
      st=1
  endif  
 !y
  if (i2-i1 /= j2-j1) then
     PRINT*,'BLOCKTRACE_HASHM: ranges not equal'
     stop
  end if
  
  tr=0

  do i=1,x%nel
     if (x%x(1,i) /=0) then
        k= INT(x%x(1,i));  l=INT(x%x(2,i))
	k1=(k-1)*st+i1;  l1=(l-1)*st+j1
	k2=(l-1)*st+i1;  l2=(k-1)*st+j1
	if (k2 > l2) then
	   j=k2; k2=l2; l2=j
	endif
	   
	if ( k == l) then
            tr=tr+getm(k1,l1,y)*x%x(3,i)
          else
            tr=tr+(getm(k1,l1,y)+getm(k2,l2,y))*x%x(3,i)
	endif       
     end if
  end do  
  end function
   
 
 function traceblock_ija(x,y,i1,i2,j1,j2,step) RESULT (tr)
! y=trace(x,y(i1:i2:step,j1:j2:step)) for ija matrices
  type (sparse_ija):: x,y
  INTEGER,INTENT(IN) :: i1,i2,j1,j2
  real (r8)::tr
  INTEGER :: i,j,k,k1,k2,l,l1,l2,m,st
  integer,optional::step

  if (present(step)) then
      st=step
     else
      st=1
  endif  
 !y
  if (i2-i1 /= j2-j1) then
     PRINT*,'BLOCKTRACE_HASHM: ranges not equal'
     stop
  end if
  
  tr=0; 
 
  do i=1,x%n
     do j=x%ia(i),x%ia(i+1)-1        
        k= i;  l=x%ja(j)
	k1=(k-1)*st+i1;  l1=(l-1)*st+j1
	k2=(l-1)*st+i1;  l2=(k-1)*st+j1
	if (k2 > l2) then
	   m=k2; k2=l2; l2=m
	endif
	if ( k == l) then
            tr=tr+getm(k1,l1,y)*x%a(j)
          else
            tr=tr+(getm(k1,l1,y)+getm(k2,l2,y))*x%a(j)
	endif       
     end do
  end do  
  end function
  
  
!----------------
! quadratic form
!----------------

function quadr_densem(x,a,y) RESULT(q)
! q=x'Ay for dense matrices
  TYPE(densem):: a
  real (rh) :: x(:),y(:)
  REAL (r8) :: q
  integer :: i,j
  !
  if (size(x) /=a%n .or. size(y) /=a%n ) then
     PRINT*,'QUADR_DENSEM: incompatible dimensions'
     stop
  end if

  q=0
  do i=1,a%n
     do j=1,a%n
        q=q+a%x(i,j)*x(i)*y(j)
     enddo
  enddo
end function

function quadr_hashm(x,a,y) RESULT(q)
! q=x'Ax for hash matrices
  TYPE(sparse_hashm):: a
  real (rh) :: x(:),y(:)
  REAL (r8) :: q
  integer :: i,j,k
  !
  if (size(x) /=a%n .or. size(y) /=a%n ) then
     PRINT*,'QUADR_HASHM: incompatible dimensions'
     stop
  end if
  !
  q=0

  do j=1,a%nel
     if (a%x(1,j) /= 0) then
        i=a%x(1,j); k=a%x(2,j)
        if (i == k) q=q+a%x(3,j)*x(i)*y(k)
        if (i < k)  q=q+a%x(3,j)*(x(i)*y(k)+x(k)*y(i))
     endif
  end do
end function

function quadr_ija(x,a,y) RESULT(q)
! q=x'Ax for ija matrices
  TYPE(sparse_ija):: a
  real (rh):: x(:),y(:)
  REAL (r8):: q
  integer :: i,j,k
  !
  if (size(x) /=a%n .or. size(y) /=a%n ) then
     PRINT*,'QUADR_IJA: incompatible dimensions'
     stop
  end if
  !
  q=0

  do i=1,a%n
     do j=a%ia(i),a%ia(i+1)-1
        k=a%ja(j)
        if (i == k) q=q+a%a(j)*x(i)*y(k)
        if (i < k)  q=q+a%a(j)*(x(i)*y(k)+x(k)*y(i))
     end do
  end do
!
end function

!----------------------------------------
! lower level solving
!----------------------------------------
subroutine solve_dense_gs(n,lhs,rhs,sol)
! finds sol in the system of linear equations: lhs*sol=rhs
!     the solution is iterative by Gauss-Seidel
 integer :: n
 real (r8) :: lhs(n,n)
 real (rh) :: rhs(n),sol(n),solnew,eps,ss
 integer :: round,i
 
 round=0
 do
   eps=0; ss=0; round=round+1
    do i=1,n
       if (lhs(i,i).eq.0) cycle
       solnew=sol(i)+default_relax*(rhs(i)-sum(lhs(i,:)*sol))/lhs(i,i)
       eps=eps+ (sol(i)-solnew)**2; ss=ss+solnew*solnew
       sol(i)=solnew
    end do
    eps=eps/ss
    if (eps < default_conv .or. round >= default_rounds) exit
 end do
 if (sparsem_msg >=2) then
   print*,'solutions computed in ',round,' rounds of iteration'
   print*,'reached convergence criterion ',eps
 endif
end subroutine 

subroutine solve_dense_symm_gs(n,lhs,rhs,sol)
! As solve_dense_gs but for summetric dense upper-stored matrix
 integer :: n
 real (r8):: lhs(:)
 real (rh):: rhs(n),sol(n),solnew,eps,ss,sum
 integer :: round,i,j
 
 round=0
 do
    eps=0; ss=0; round=round+1
    do i=1,n
       if (lhs(pos_symm(i,i,n)).eq.0) cycle
       sum=0
       do j=1,n
          sum=sum+lhs(pos_symm(i,j,n))*sol(j)
       enddo
       solnew=sol(i)+default_relax*(rhs(i)-sum)/lhs(pos_symm(i,i,n))
       eps=eps+ (sol(i)-solnew)**2; ss=ss+solnew*solnew
       sol(i)=solnew
    end do
    eps=eps/ss
    if (eps < default_conv .or. round >= default_rounds) exit
 end do
 if (sparsem_msg >=2) then
     print*,'solutions computed in ',round,' rounds of iteration'
     print*,'reached convergence criterion ',eps
 endif    
end subroutine 


subroutine solve_sparse_hashm_gs(n,lhs,rhs,sol,filled)
! As solve_dense_gs but for sparse hash assumed upper-stored matrix
 use kinds
 integer :: n,filled
 real (rh) :: lhs(:,:)
 real (rh) :: rhs(n),sol(n),solnew,eps,ss,sum
 integer :: round,i,j,l1,l2,row
 integer:: ia(n+1)   !row pointer as IA in IJA format
 real (rh):: diag(n),rhsadj(n) !diagonals and adjusted right hand sides

 ia=0;diag=0
 ! create row pointers
   row=1
   ia(row)=1
   do i=1,filled
      l1=lhs(1,i)
      if (l1 == lhs(2,i)) diag(l1)=lhs(3,i)
      if (l1 /= row) then
         row=l1
         ia(row)=i
      endif
   enddo
   ia(row+1)=filled+1
 ! Add pointers to empty rows
   do i=row,2,-1
      if (ia(i) == 0) ia(i)=ia(i+1)
   enddo

 round=0;
 do
    eps=0; ss=0; round=round+1; rhsadj=rhs
    do i=1,n
       if (diag(i) == 0) cycle
       sum=0
       do j=ia(i),ia(i+1)-1
          sum=sum+lhs(3,j)*sol(int(lhs(2,j)))
       enddo
       solnew=sol(i)+default_relax*(rhsadj(i)-sum)/diag(i)
       eps=eps+ (sol(i)-solnew)**2; ss=ss+solnew*solnew
       sol(i)=solnew
      ! adjust right hand sides because of triangular storage
       do j=ia(i),ia(i+1)-1
          l2=lhs(2,j)
          rhsadj(l2) = rhsadj(l2)-lhs(3,j)*sol(i)
       enddo
    end do
    eps=eps/ss
    if (eps < default_conv .or. round >= default_rounds) exit
    print '(" round ",i6,"  convergence=",g12.4)' ,round,eps
 end do
 if (sparsem_msg >=2) then
    print*,'solutions computed in ',round,' rounds of iteration'
    print*,'reached convergence criterion ',eps
 endif   
end subroutine

subroutine solve_sparse_ija_gs(n,ia,ja,a,rhs,sol,nel)
! As solve_dense_gs but for sparse ija assumed upper-stored matrix
 integer :: n,nel,ia(:),ja(:)
 real (r8) :: a(:)
 real (rh) ::rhs(n),sol(n),solnew,eps,ss,sum
 integer :: i,j,round
 real (rh):: diag(n),rhsadj(n) !diagonals and adjusted right hand sides

 diag=0
 do i=1,n
   do j=ia(i),ia(i+1)-1
      if (i == ja(j)) diag(i)=a(j)
   enddo
 enddo
 
 round=0;
 do
    eps=0; ss=0; round=round+1; rhsadj=rhs
    do i=1,n
       if (diag(i) == 0) cycle
       sum=0
       do j=ia(i),ia(i+1)-1
          sum=sum+a(j)*sol(ja(j))
       enddo
       solnew=sol(i)+default_relax*(rhsadj(i)-sum)/diag(i)
       eps=eps+ (sol(i)-solnew)**2; ss=ss+solnew*solnew
       sol(i)=solnew
      ! adjust right hand sides because of triangular storage
       do j=ia(i),ia(i+1)-1
          rhsadj(ja(j)) = rhsadj(ja(j))-a(j)*sol(i)
       enddo
    end do
    eps=eps/ss
    if (eps < default_conv .or. round >= default_rounds) exit
 end do
 if (sparsem_msg >=2) then 
    print*,'solutions computed in ',round,' rounds of iteration'
    print*,'reached convergence criterion ',eps
 endif   
end subroutine

!=====================================================================
!Conversion routines
!=====================================================================

subroutine convert_half_dense(y,x)
! converts symmetric half-stored x into full-stored y
type (densem), intent(inout) :: y
type (dense_symm), intent(in) :: x
real (rh) :: a
integer :: i,j

call zerom(y,x%n)
do i=1,x%n
   do j=1,x%n
      a=x%x(pos_symm(i,j,x%n))
      call addm(a,i,j,y)
      if (i /= j) call addm(a,j,i,y)
   enddo
enddo
end subroutine

subroutine convert_hash_densem(y,x)
! converts hash matrix x into full-stored y
type (densem), intent(inout) :: y
type (sparse_hashm), intent(in) :: x
integer i

call zerom(y,x%n)
do i=1,x%nel
   if (x%x(1,i).ne.0) then
      call addm(x%x(3,i),int(x%x(1,i)),int(x%x(2,i)),y)
      if (int(x%x(1,i)) /= int(x%x(2,i))) &
                     call addm(x%x(3,i),int(x%x(2,i)),int(x%x(1,i)),y)
   endif
enddo
end subroutine


subroutine convert_hash_ija(y,x)
! converts symmetric hash matrix x into ija matrix y
 type (sparse_ija), intent(inout) :: y
 type (sparse_hashm), intent(in) :: x
 
 call convert_hash_ija_general(y,x)
end subroutine

 subroutine convert_hash_ija_general(y,x,conv_type)
 ! converts symmetric hash matrix x into ija matrix y, fully stored
 ! optional conv_type specifies the type of conversion:
 !    conv_upper_full - from half- upper stored to full stored
 !    conv_straight - as the original matrix (default)
 
 type (sparse_ija), intent(inout) :: y
 type (sparse_hashm), intent(in) :: x
 integer,optional::conv_type
 integer:: i,j,k,l,row(x%n),conv_type_o=0
 real (rh) :: zero=0


   if (present(conv_type)) conv_type_o=conv_type
   
!  determine number of entries per row excluding diagonals
   row=0
   do k=1,x%nel
      i= x%x(1,k)
      if (i > 0) then 
         j=x%x(2,k)
	 select case (conv_type_o)
	   case (conv_upper_full)
               if (i <j) then
	          row(i)=row(i)+1
                  row(j)=row(j)+1
               endif   
	   case default
	       if (i /= j) row(i)=row(i)+1    
	 end select      
      endif
   enddo

!  initialize ija structure
   call zerom(y,x%n)
   y%n=x%n; y%nel=sum(row)+y%n
   allocate(y%ia(y%n+1),y%ja(y%nel),y%a(y%nel))

!  set row pointers and initialize diagonal elements
   y%ia(1)=1  
   do i=1,y%n
      j=y%ia(i)
      y%ja(j)=i
      y%a(j)=0
      y%ia(i+1)=y%ia(i)+row(i)+1
   enddo

!  set values
   row=1	!diagonals are always first in row
   do k=1,x%nel
      i= x%x(1,k)
      if (i > 0) then
         j=x%x(2,k)
	 if (i == j) then
	    y%a(y%ia(i))=x%x(3,k)	!diagonal
         elseif (i < j) then
         
	    l=y%ia(i)+row(i)	!regular or upper diagonal
            y%ja(l)=j
            y%a(l)=x%x(3,k)
            row(i)=row(i)+1
	    
            if (conv_type_o == conv_upper_full) then
               l=y%ia(j)+row(j)	!lower diagonal
               y%ja(l)=i
               y%a(l)=x%x(3,k)
               row(j)=row(j)+1
	    endif   
	    
         endif
      endif
   enddo
   !print*,'x'
   !call printm(x)
   !print*,'y'
   !call printm(y)
   !print*
   !call printm(y,'internal')
   
 end subroutine




subroutine convert_ija_densem(y,x)
! converts ija matrix x into full matrix y
 type (densem), intent(inout) :: y
 type (sparse_ija), intent(in) :: x
 integer i,j
 real (rh) ::val

 call zerom(y,x%n)
 y%n=x%n
 
 do i=1,x%n
    do j=x%ia(i),x%ia(i+1)-1
       if (i <= x%ja(j)) then	!ignore lower-stored elements
          val=x%a(j)
          call addm(val,i,x%ja(j),y)
          if (i /= x%ja(j)) call addm(val,x%ja(j),i,y)
       endif  
    enddo
 enddo
 end subroutine

subroutine convert_ija_hash(y,x)
! converts ija matrix x into hashm matrix y
 type (sparse_hashm), intent(inout) :: y
 type (sparse_ija), intent(in) :: x
 integer i,j
 real (rh) ::val

 call zerom(y,x%n)
 y%n=x%n
 
 do i=1,x%n
    do j=x%ia(i),x%ia(i+1)-1
       val=x%a(j)
       call addm(val,i,x%ja(j),y)
    enddo
 enddo
 end subroutine

subroutine convert_densem_hashm(y,x)
! converts dense matrix x into hash-stored y
type (sparse_hashm), intent(inout) :: y
type (densem), intent(in) :: x
real (rh) :: a
integer :: i,j

call zerom(y,x%n)
do i=1,x%n
   do j=1,x%n
      if (x%x(i,j).ne.0) then
         a=x%x(i,j)
         call addm(a,i,j,y)
      endif
   enddo
enddo
end subroutine


subroutine convert_densem_ija(y,x)
! converts dense matrix x into ija matrix y
 type (sparse_ija), intent(inout) :: y
 type (densem), intent(in) :: x
 integer, allocatable:: tmp(:)
 integer:: i,j,k

 allocate(tmp(x%n))
 tmp=0;
 ! create row counts of nonzeros
   do i=1,x%n
      do j=i,x%n
         if (x%x(i,j) /= 0 .or. i==j) then
             tmp(i)=tmp(i)+1
         endif
      enddo
   enddo

   call zerom(y,x%n)
 
   y%n=x%n; y%nel=sum(tmp)
   allocate(y%ia(x%n+1),y%ja(y%nel),y%a(y%nel))

   y%ia(1)=1
   do i=1,x%n
      y%ia(i+1)=y%ia(i)+tmp(i)
   enddo

! Load values
   tmp=0
   do i=1,x%n
      do j=i,x%n
         if (x%x(i,j) /= 0 .or. i==j) then
             k=y%ia(i)+tmp(i)
             y%ja(k)=j
             y%a(k)=x%x(i,j)
             tmp(i)=tmp(i)+1
         endif
      enddo
   enddo
   deallocate(tmp)
end subroutine


!=====================================================================
!Copy routines
!=====================================================================

subroutine copy_dense(y,x)
! copies full-stored x into full-stored y
type (densem), intent(inout) :: y
type (densem), intent(in) :: x

call zerom(y,x%n)
y%x=x%x
end subroutine


subroutine copy_symm(y,x)
! copies half-stored x into half-stored y
type (dense_symm), intent(inout) :: y
type (dense_symm), intent(in) :: x

call zerom(y,x%n)
y%x=x%x
end subroutine


subroutine copy_hashm(y,x)
! copies hash x into hash y
type (sparse_hashm), intent(inout) :: y
type (sparse_hashm), intent(in) :: x
type (sparse_ija)::scr

select case(x%status)
   case(2)	! matrix sorted, copy by conversion to ija and back
      scr=x
      y=scr
      call reset(scr)
   case(1) 
     call zerom(y,x%n)
     deallocate(y%x)
     y%n=x%n
     y%nel=x%nel
     y%filled=x%filled
     y%status=x%status
     allocate(y%x(size(x%x,dim=1),size(x%x,dim=2)))
     y%x=x%x
   case default
     print*,'sparse_hashm matrix empty'
     stop
end select       
end subroutine


subroutine copy_ija(y,x)
! copies IJA x into IJA y
type (sparse_ija), intent(inout) :: y
type (sparse_ija), intent(in) :: x

call zerom(y,x%n)
y%n=x%n
y%nel=x%nel
allocate(y%ia(size(x%ia)),y%ja(size(x%ja)),y%a(size(x%a)))
y%ia=x%ia
y%ja=x%ja
y%a=x%a
end subroutine


end module
