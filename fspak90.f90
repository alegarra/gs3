module sparseop
! sparse matrix operations including fspak90
  use sparsem
  implicit none

   interface fspak90
      module procedure fspak90r8,fspak90r4
   end interface
   
! matrix by vector multiplication   
   interface mult
       module procedure mult_densem_vec,mult_hashm_vec,mult_ija_vec,&
                        mult_vec_densem,mult_vec_hashm,mult_vec_ija
   end interface
   
! matrix by scalar multiplication   
   interface multmatscal		
	module procedure mult_densem_scal, mult_hashm_scal,mult_ija_scal
   end interface    
 
 
  contains



subroutine fspak90r4(operation,ija,rs4,sol4,det,msglev,maxmem,rank,rs8,sol8)
! interface to fspak; most parameters above optional and memory
! allocation dynamic
  
  interface
    subroutine fspak(mode,n,ia,ja,a,sol,flag,io1,io2,memory,&
                 mem_needed,work,i1,i2,i3,i4,i5,i6,rank)
    use kinds
    integer :: mode,n,ia(1),ja(1), flag,io1,io2,memory,&
                 mem_needed,i1,i2,i3,i4,i5,i6,rank 
    real (r8) :: a(1),sol(1),work(1)                       
    end subroutine
  end interface  
  
  character (len=*)  :: operation
  type (sparse_ija):: ija
  real(r8),optional :: rs8(:),sol8(:),det
  real(r4),optional :: rs4(:),sol4(:)
  real (r8), allocatable:: sol(:)
 !optional parameters
  integer,optional::rank,msglev,maxmem
  integer :: msglev_o,maxmem_o
  integer,save :: rank_o,status=-1,last_n=0
        ! status = -1 (initial), 1 (ordering done), 2 (factorization done)
  integer,save :: memory=100        !initial memory, modified by later programs
  real (r8),save,pointer :: work(:)=>null()
  real , external:: second


  if (last_n /= ija%n) then
        status=min(status,0)     !redo everything if dimension changes
        last_n=ija%n
  endif
!
  if (.not. present(msglev)) then
        msglev_o=0   ! default message level = minimal
      else
        msglev_o=msglev
  endif
  if (.not. present(maxmem)) then
        maxmem_o=huge(maxmem)         ! infinite memory limit if maxmem missing
     else
        maxmem_o=maxmem
  endif
  select case (operation(1:3))
       case ('fac')             !factorization
            status=min(status,1)        !nullify factorization if done
            call do_fact
       case('sol')            !solution
            allocate(sol(ija%n))
            if (present(sol8) .and. present(rs8)) then
                sol=rs8
                call do_fact
                call run_fspak('solve')
                sol8=sol
            elseif (present(sol4) .and. present(rs4)) then
                sol=rs4
                call do_fact
                call run_fspak('solve') 
                sol4=sol
              else
                print*,'FSPAK90: optional parameter SOL or RS missing'
                stop
            endif
            deallocate(sol)
       case('inv')           !sparse inversion
            call do_fact
            call run_fspak('inverse')
            status=1            !factorization lost
       case('res')            !reset storage
            if (status > -1) deallocate(work)
            status=-1
            memory=100
       case('det')              !determinant
            if (.not. present(det)) then
                print*,'FSPAK90: optional parameter DET missing'
              else
                call do_fact
                call run_fspak('det')
             endif
       case('lde')             !log determinant
            if (.not. present(det)) then
                print*,'FSPAK90: optional parameter DET misiing'
              else
                call do_fact
                call run_fspak('ldet')
            endif
       case('sta')             !print statistics
            call run_fspak('statistics')
       case default
          print*,'FSPAK90 operation ',operation,' unknown'
  end select


  contains
 
 subroutine do_fact
 ! executes all stages of fspak so that the factorization is available

 if (status == -1) then                
    allocate(work(memory/2+1))  ! real(r8) needed as integer
    status=0
 endif

 if (status == 0) then
    call run_fspak('ordering')
    call run_fspak('symbolic_fact')
    status=1
 endif
 
 if (status == 1) then
    call run_fspak('numeric_fact')
    status=2
 endif
 end subroutine do_fact


 subroutine run_fspak(op)
 ! runs operation op of fspak, adding more memory if necessary
 character (len=*) :: op
 real(r8),pointer,save::work1(:)=>null()
 real(r8) :: sc(2)   !scratch
 integer :: mem_needed,i,flag

! print*,op
 do
   select case(op)
      case ('ordering')
            open(99,file='fspak90.ord')
            call fspak(10,ija%n,ija%ia,ija%ja,ija%a,sc,flag,6,99,memory,&
                 mem_needed,work,i,i,i,i,i,i,rank_o)
            close (99)
            if (msglev_o > 2) then
               print*,'FSPAK90: size=',ija%n,'  #nonzeroes=',ija%nel
            endif
            if (msglev_o > 1) print*,'ordering time =',second()
      case('symbolic_fact')
            call fspak(20,ija%n,ija%ia,ija%ja,ija%a,sc,flag,6,99,memory,&
                 mem_needed,work,i,i,i,i,i,i,rank_o)
            if (msglev_o > 1) print*,'symbolic factorization time=',second()
      case('numeric_fact')
            call fspak(40,ija%n,ija%ia,ija%ja,ija%a,sc,flag,6,99,memory,&
                    mem_needed,work,i,i,i,i,i,i,rank_o)
            if (msglev_o > 1) print*,'numerical factorization time=',second()
      case('solve')
            call fspak(50,ija%n,ija%ia,ija%ja,ija%a,sol,flag,6,99,memory,&
                    mem_needed,work,i,i,i,i,i,i,rank_o)
      case('det')
            call fspak(54,ija%n,ija%ia,ija%ja,ija%a,sc,flag,6,99,memory,&
                    mem_needed,work,i,i,i,i,i,i,rank_o)
            det=sc(1)
      case('ldet')
            call fspak(55,ija%n,ija%ia,ija%ja,ija%a,sc,flag,6,99,memory,&
                    mem_needed,work,i,i,i,i,i,i,rank_o)
            det=sc(1)
      case('inverse')
            call fspak(61,ija%n,ija%ia,ija%ja,ija%a,sc,flag,6,99,memory,&
                 mem_needed,work,i,i,i,i,i,i,rank_o)
             if (msglev_o >1) print*,'inversion time=',second()
      case('statistics')
            call fspak(80,ija%n,ija%ia,ija%ja,ija%a,sc,flag,6,99,memory,&
                    mem_needed,work,i,i,i,i,i,i,rank_o)
   end select
   if (flag == 0) then
          exit
      elseif (flag >= 4000000) then
        print*,'FSPAK: Zero or negative pivot for row ROW',&
            'during numerical factorization in row:',flag-4000000
        stop
     ! if memory insufficient, assign more
      elseif (flag > 100) then
        if (mem_needed > maxmem_o) then
            print*,'FSPAK90: memory needed:',mem_needed
            print*,'         maxmimum memory allowed:',maxmem
          else
            if (msglev_o>0) print*,'memory increased from:',memory,' to:',&
                                                          mem_needed
            allocate(work1(mem_needed/2+1))    !create new memory
            work1(1:memory/2+1)=work           !copy old data
            deallocate(work)           !release old memory
            work=>work1                !let work point to memory of work1
            memory=mem_needed
        endif
      else
        print*,'fspak: unknown error flag: ',flag
        stop
   endif
 enddo
 if (present(rank)) rank=rank_o
 
 end subroutine

end subroutine fspak90r4


!-------------------------------------------------------------------------
! Overloading of fspak90 so that it works with real*8 solutions and RHS
!-------------------------------------------------------------------------

subroutine fspak90r8(operation1,ija1,rs1,sol1) !,msglev1,maxmem1,rank1)
! overloaded fspak90 so that solutions and right hand sides can be 
! double precision without using optional parameters 
 
  character (len=*)  :: operation1
  type (sparse_ija):: ija1
  real(r8):: rs1(:),sol1(:)
 !optional parameters
 ! integer,optional::rank1,msglev1,maxmem1
  call fspak90r4(operation1,ija1,rs8=rs1,sol8=sol1)
end subroutine fspak90r8


!--------------------------------------
! Multiplication of matrix by vector
!--------------------------------------

function mult_densem_vec(A,y) result (z)
! z = A * y for A in densem format
type (densem)::A
real (r8)::y(a%n),z(a%n)
!
z=matmul(A%x,y) 
end function


function mult_vec_densem(yy,A) result (z)
! z = y' * A  for A in densem format
type (densem)::A
real (r8)::yy(a%n),z(a%n)
!
z=matmul(yy,A%x) 
end function


function mult_hashm_vec(A,y) result (z)
! z = A * y for A in sparse_hashm form.
type (sparse_hashm)::A
real (r8)::y(a%n ),z(a%n)
integer::i,j,k
!
 z=0
 do j=1,a%nel
   if (a%x(1,j) /=0) then
      i=a%x(1,j); k=a%x(2,j)
     if (i == k) z(i)=z(i)+a%x(3,j)*y(k)
     if (i < k)  then
        z(i)=z(i)+a%x(3,j)*y(k)
        z(k)=z(k)+a%x(3,j)*y(i)
     endif    
   endif
 end do
end function
 
 
 function mult_vec_hashm(yy,A) result (z)
! z = y' * A for A in sparse_hashm form.
type (sparse_hashm)::A
real (r8)::yy(a%n ),z(a%n)
integer::i,j,k
!
 z=0
 do j=1,a%nel
   if (a%x(1,j) /=0) then
     i=a%x(1,j); k=a%x(2,j)
     if (i == k) z(i)=z(i)+a%x(3,j)*yy(k)
     if (i < k)  then
        z(i)=z(i)+a%x(3,j)*yy(k)
        z(k)=z(k)+a%x(3,j)*yy(i)
     endif    
   endif
 end do
end function
 
  

function mult_ija_vec(A,y) result (z)
! z = A * y for A in sparse_ija form.
type (sparse_ija)::A
real (r8)::y(a%n),z(a%n)
integer::i,j,k
!
   z=0
   do i=1,a%n
     do j=a%ia(i),a%ia(i+1)-1
        k=a%ja(j)
        if (i == k) z(i)=z(i)+a%a(j)*y(k)
        if (i < k)  then
            z(i)=z(i)+a%a(j)*y(k)
            z(k)=z(k)+a%a(j)*y(i)
        endif    
     end do
  end do  
end function


function mult_vec_ija(yy,A) result (z)
! z = y' * A for A in sparse_ija form.
type (sparse_ija)::A
real (r8)::yy(a%n),z(a%n)
integer::i,j,k
!
   z=0
   do i=1,a%n
     do j=a%ia(i),a%ia(i+1)-1
        k=a%ja(j)
        if (i == k) z(i)=z(i)+a%a(j)*yy(k)
        if (i < k)  then
            z(i)=z(i)+a%a(j)*yy(k)
            z(k)=z(k)+a%a(j)*yy(i)
        endif    
     end do
  end do  
end function


! -----------------------------------
!  Multiplication of matrix by scalar
!-------------------------------------

subroutine mult_densem_scal(A,y)
! A = A * y for A in densem format
type (densem)::A
real (r8)::y
!
A%x=A%x*y  
end subroutine

subroutine mult_hashm_scal(A,y)
! A = A * y for A in sparse_hashm form.
type (sparse_hashm)::A
real (r8)::y
integer::i,j,k
!
 A%x(3,:)=A%x(3,:)*y
 
end subroutine
 

subroutine mult_ija_scal(A,y)
! A = A * y for A in sparse_ija form.
type (sparse_ija)::A
real (r8)::y
integer::i,j,k
!
   A%a=A%a*y
  
end subroutine

end module
