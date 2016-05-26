!    This program is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    This program is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with this program.  If not, see <http://www.gnu.org/licenses/>.






! this is a program to estimate variances in a MM including
! mean, one fixed effect (gender) and SNP's
! THIS IS THE ONLY MODEL BY NOW!!
! All SNP have the same variance
! Estimation is by Gibbs Sampling with iteration on data with or without variance component estilation
! or BLUP by Gauss Seidel
! using a "solving by effects" technique (e.g. Notes by Misztal, Janss & DeJong 1999,
! Misztal & Gianola 1987)
! includes one animal effect which assumes ONE RECORD PER ANIMAL
! A-1 is stocked in dense
! predict is read as option and produces a file 'predictions' with predicted phenotypes (with sol taken from the
! solutions file read in the par file)

! Andres Legarra
! 17/1/07
! saving for continuation implemented
! Uses:
!     -a parameter file (example: together.animal.par) This par file
!      includes an option to solve the system by Gauss Seidel (BLUP) or MCMC  (MCMCBLUP)
!      with known variance components, or to estimate these and
!      the solutions (VCE). Also options for MCMC and convergence
!     -pedigree and data files detailed in the parameter file
!     -a file continuation.txt to continue the MCMC if needed
!     -the continuation of the var samples is in varfile_cont
! Produces:
!     -a file with samples of the variance components, name detaild in the par file (if option VCE)
!     -a file with solutions
! TODO: 
!     -X and Z are redundanT!!
!     -model flexibility
! last change 8/6/07

! gs_sparse: I try to put more model flexibility. How:
! variable number of effects ?
! sparse storage of A-1 largely overpass gs_aonly
! should I stock X in sparse as well??


! this is a modification to run a mixture model, as proposed by JME
! a_i ~ (1-p)dirac(0) + p N(0,vara) for the genomic effects only!
! by classical MCMC for mixture models
! (stochastic search variable selection) 
! There is no problem of label switching here.
! ML should be attainable without much effort on the basis of BLUP solutions
! I use variable i_inmodel to say if a given SNP is in or out. This is later used for estimation of p 
! serious bug corrected 7/5/07

! 1/8/07
! This is a copy from gs_sparse_mixt, which is able to 
!   -select given rows from the data file
!   -include one permanent effect "cage" (not difficult to be extended)
!   -disable mixture approach if needed
!   -model flexibility (
! SNP have to be kept as the last column in data file for ease of reading
! ToDo:
!   -handle missing values (skipping them or doing data augmentation?)
!   verificar los resultados
! 4/3/08
! I have added writing of EBVs
! Doc
! Corrected inconsistency in predict and avinmodel
! Corrected bug in computation of length(in_data)
! 29/5/08 Corrected bug in computation of length(in_data)
! 10/08/09 Missing covariates for snp's allowed (with effect=0)
!          Corrected too small rec options
!          Introduced heterogeneous variances
!          Corrected compute_EBVs (not sure)
! 19/11/09 Corrected: it was inconsistent to compute a_hat dividing by avinmodel because it is not the posterior mean of a
!  and complicates computing
!  Predict_y rewritten (predicts with mean included; not veru useful for cross-validation unless pre-corrected phenotypes)
! 29/3/11 strong rewriting since...
!         today, including prints of vara * 2*pbar*qbar
!         also, does not check "animals in genotypes" if there are no gentypes in the model
! 31/5/11 GS3 with Bayesian Lasso
!         recycling from BLUP_lasso, only for additive SNP effects and in BL2Var mode (not in BL1Var)
!
! 21/10/11 Reads parameter file from command line: legarra@snp$ ./GS3 myfile.par
! 04/10/2012 corrected error in coputation of eps
! 1/2/2013 added checking of at least a fixed effect
!           initialized W (thanks to Vitaly)
! 12/2/13 bug corrected: read geno file when no genotypes in the model. All kinds of weird outputs   
! 19/02/2013 corrected "name" of EBVs on output as these are not EBVs
!            changed to coherent system: aa-> 0 -> -1a  ; aA -> 1-> 0a ; AA -> 2-> +1a in set up of equations
!               THIS IS CHANGE IS INCOMPATIBLE BACKWARDS (a effects will switch sign)
!            missing dominant set to average of the column
!            freq computed for "2" 
! 13/03/2013 weights=0 converted in 1d-50
! 6/8/2013 Many small improvements: corrected bug in VCE with missing data, in computation of ngeno, check for a fixed effect...
! 7/8/13 bug in reading genotypes > 500000
! 8/5/2014 bug in eps (criteria changed), bug in rhs buildup with markers
! 6/1/2016 the program does not stop (but gives a warning) if there are more phenotypes than
! genotypes, in order to allow repeated records 
! 24/5/2016 Bayes Factors for MCMCGBLUP models

! TODO
! modify guess of format so that we can include <50 SNPs
! include prevalence
!


program gs_bayesian_lasso
use kinds
use random
use stat
use sparsem
use denseop
use boots
use util
use aux_options

! we parameterize the model in terms of an effect
!   +1/2a for SNP 1
!   -1/2a for SNP 2
!   +d for SNP combination 12
! therefore: 11 -> a 12 -> d 22 -> -a


implicit none
character(20):: version= '2.5.0'
character(20):: date='24 May 2016 '
integer,parameter::maxsnp=1000000
! parameters of the model
integer:: neff,&
          nSNP
                    
integer,allocatable:: nlev(:),offset(:),pos_eff(:),type_eff(:)
character(20):: type_this_eff,type_eff_char(7)=(/'cross        ',&
                                                 'cov          ',&
                                                 'add_animal   ',&
                                                 'perm_diagonal',&
                                                 'add_SNP      ',&
                                                 'dom_SNP      ',&
                                                 'ind_SNP      '/)
integer,parameter:: cross=1,cov=2,add_animal=3,perm_diagonal=4,add_SNP=5,dom_SNP=6,ind_SNP=7

logical,allocatable:: inmodel(:),&
                      i_inmodel(:)
integer :: neq,efaSNP,efdSNP,efanim,efp,traitcol,maxrec=10000

type design_matrix 
  ! this uses TR 15581 standard, available in most compilers in particular gfortran and g95
  real(r8),allocatable:: cov(:) !value of the covariable
  integer,allocatable:: eq(:) !number of equation
  integer,allocatable:: rec(:) !number of record
  integer,allocatable:: ia_eq(:,:) !address of the first and last elements in X for the i-th eq
  logical:: weighted !whether X has been transformed (multiplying by weights) or not
end type design_matrix
  
real(r8),allocatable:: in_data(:)                    
type(design_matrix) ::X ! for "sparse" effects
real(r8),allocatable::y(:),ywiggle(:),ywiggleold(:),&
                      phen(:),& !phenotype observed (for binary trait)
                      xpx(:),& ! diagonal of LHS
!                      Z(:,:),& !add_SNP effects 
                      W(:,:),& !dom_SNP effects 
                      sol(:),& !vector of solutions
                      freq(:),& ! allelic frequencies of "2"
                      sumsol(:),& !vector of accumumated solutions
                      ssqsol(:),& !vector of sum of squares of acumumated solutions
                      diagAinv(:),& !diag of A-inverse
                      avinmodel(:),& ! vector with p of a given effect being accounted for in the model
                                     ! in fact this is done only for SNP effects
                      sumxpx(:),& !vector with product X 1 filtered with i_inmodel
                      D(:),& !diagonal vector with inverses of variances (A-1, D-1, etc) for each effect
                      weights(:),& ! diagonal matrix with weights for heterogeneous residual variances
                      tau2(:),&      ! individual variances of SNP effects (only for additive ones)
                      sumtau2(:),& ! its sum,
                      sstau2(:) ! sum of squares
                      
!integer(2),allocatable::Z(:,:) !add_SNP effects 
real(r8),allocatable::Z(:,:) !add_SNP effects 
                      
type(sparse_hashm):: Ainv_hash                      
type(sparse_ija):: Ainv_ija
real(r8):: vara,vard,vare,varg,varp,lambda2
real(r8):: svara,svard,svare,svarg,svarp,slambda2
real(r8):: ssvara,ssvard,ssvare,ssvarg,ssvarp,sslambda2
real(r8):: vara_ap,vard_ap,vare_ap,varg_ap,varp_ap ! apriori
real(r8):: dfvara,dfvare,dfvard,dfvarg,dfvarp !degrees of freedom for the a priori distribution
real(r8):: val,lhs,rhs,eps,conv_crit=1e-10,edc
!read-in variables
character(len=maxsnp):: genome 
character(maxsnp):: line
integer,allocatable:: genotype(:,:) !nanim,nsnp
integer,allocatable:: SNPoneanim(:),& !(nSNP) !SNP combination for one animal
                      id(:),& !id for the n-th record
                      id_anim(:),& !animal id for the n-th record
                      id_geno(:) !animal id for the n-th genotype
character(200):: parfile,datafile,method,varfile,solfile,pedfile,genofile,informat='',fileweights

integer:: i,j,k,l,io,ndata,pos,iterout,iterin,pos1,pos2,iter,nanim,beginouter,record_id,hetvarcol,ngeno,pos_anim
integer:: niter,burnin,thin,correction
logical:: simulation=.false.,MCMCBLUP=.false.,GaussSeidel=.false.,VCE=.false.,&
          continuation=.false.,predict=.false.,use_mixture=.false.,BL1Var=.false.,BL2Var=.false.,&
          BinaryTrait=.false.,weights_snp=.false.
logical:: yWZ_weighted=.false.

! new variables for the mixture model
real(r8):: pa(2),aprioria(2),sumpa(2) !mixture parameters for the a dditive SNP effect
real(r8):: pd(2),apriorid(2),sumpd(2) !mixture parameters for the a dditive SNP effect
real(r8):: like1,like2,aposteriori(2),v0,v1,rj
integer:: includeda,includedd

integer:: saving=10000
integer:: ip_snp
real(r8):: sumpq
real(r8)::MisVal=-9999d0
real:: aa(20)
integer:: nmis=0,colweights=0

! new variables for BayesFactor
logical:: lBF=.false.
integer:: nBF,firstBF,lastBF
! at which points do we evaluate BF (fixed size windows)
real(r8),allocatable:: covBF(:,:,:)


call print_version()
call get_command_argument(1,value=parfile,status=io)
if (parfile=="--version") then
  stop
endif  
if(io/=0) then
        print *, ('what parameter file?')
        read(*,*)parfile
endif
open(unit=10,file=parfile,status='old')
efaSNP=0
efdSNP=0 
efanim=0
efp=0

call read_par_file()

!call get_clock_seeds()

open(unit=1,file=datafile,status='old')
!open(unit=2,file='yX.txt',status='replace')
if(VCE) then
  if (continuation) then
    open(unit=3,file=trim(varfile)//'_cont',status='replace') !create new file
  else
    open(unit=3,file=varfile,status='replace')
  endif
endif
open(unit=14,file=trim(parfile)//'_EBVs',status='replace')  


!prediction  
if (predict) then
  open(unit=4,file=solfile,status='old')
  open(unit=12,file='predictions',status='replace')  
endif  


! get size
ndata=0
do
  !read(1,trim(informat),iostat=io)in_data
  read(1,*,iostat=io)in_data
  do i=1,neff
    if( (type_eff(i)==cross).or.&
        !(type_eff(i)==cov).or.& covariates are assumed 1 level (i.e. no nesting)
        (type_eff(i)==perm_diagonal)  ) then
        if (nlev(i)<nint(in_data(pos_eff(i)))) then
          print *,'BEWARE'
          print *,'number of levels for effect ',i,'type ',type_eff_char(type_eff(i)), &
                  'changed from',nlev(i),'to',nint(in_data(pos_eff(i))),'line',ndata+1
          nlev(i)=max(nlev(i),nint(in_data(pos_eff(i))))
        endif
    else if (type_eff(i)==cov) then
      nlev(i)=1
    endif        
  enddo
  
  if(io/=0) exit
  ndata=ndata+1
enddo

rewind(1)

nanim=0
if(efanim/=0) then
  open(unit=7,file=pedfile,status='old')
  nanim=0
  do
    read(7,*,iostat=io)
    if(io/=0) exit
    nanim=nanim+1
  enddo
  rewind(7)
  nlev(efanim)=nanim
endif
ngeno=0
if( efasnp/=0 .or. efdsnp /=0) then 
  open(unit=8,file=genofile,status='old',recl=maxsnp)
        ! find out format
        ! this whole section is from I Aguilar
        ! determine Number of snp and postion of first marker in SNP_file
        !
        ! it is further assumed that the id is in the middle (AL)
        ! Note that it does not work properly if the last snp is before col 50
        read(8,'(a)') line
        ip_snp=index(line(1:50)," ",back=.true.)+1
        ! this trick is to be able to read a given number of SNPs for testing the software
        ! AL
        if(nSnp>=0)then
                nSnp=index(line(ip_snp:)," ")-1
        else
                nSnp=-nSnp
        endif
        write(informat,'(a,i0,a,i0,a)') '(i',ip_snp-2,',1x,',nsnp,'i1)'
        if (nSNP==-1) then
         print '(a,i0,a)', 'Number of SNPs greater than ', maxsnp,&
                           ', or error in reading SNP file'
         stop
        endif
        print*,'Column position in file for the first marker: ',ip_snp
        print*,'Format to read SNP file: ',informat
        print*,'Number of SNPs :', nSNP
  rewind(8)
  do
    !read(8,*,iostat=io)
    read(8,informat,iostat=io)
    if(io/=0) exit
    ngeno=ngeno+1
    if(mod(ngeno,1000)==0) write(*,'(a1)',advance='NO'),'*'
  enddo
  rewind(8)
endif    


where (type_eff==add_SNP)
  nlev=nSNP
end where
where (type_eff==dom_SNP)
  nlev=nSNP
end where
! put code to guess number of levels for fixed effects

neq=sum(nlev)

call print_info()

offset(1)=0
do i=2,neff !get the address of the last level of the i-1th effect
  offset(i)=sum(nlev(1:i-1))
enddo  
print *,offset

!set up of X
allocate(X%cov(maxrec))
allocate(X%rec(maxrec))
allocate(X%eq(maxrec))
allocate(X%ia_eq(neq,2))
X%cov=0
X%rec=0
X%eq=0
X%ia_eq=0
X%weighted=.false.


allocate( y(ndata),&
          id(ndata),&
          id_geno(ngeno),&
          id_anim(ndata),&
          ywiggle(ndata),&
          ywiggleold(ndata),&
          phen(ndata),&
          weights(ndata),&
          xpx(neq),&
          D(neq),&
          i_inmodel(neq),&
          diagAinv(nanim),&
          genotype(ngeno,nsnp),&
          SNPoneanim(nSNP),&
          freq(nSNP),&
          sumsol(neq),&
          ssqsol(neq),&
          sumxpx(neq),&
          avinmodel(neq),&
          tau2(neq),&
          sumtau2(neq),&
          sstau2(neq),&
          sol(neq) ) 
if(lBF) then
	if(firstBF==0) firstBF=1
	if(lastBF==0) lastBF=nsnp
	print *,'Bayes Factor is computed from marker:',firstBF,'to marker',lastBF
	allocate(covBF(firstBF:lastBF,0:nBF-1,0:nBF-1))
	covBF=0
endif
X%cov=0d0; X%rec=0; X%eq=0
y=0d0; ywiggle=0d0 ; phen=0d0
sol=0d0
diagAinv=0d0
xpx=0d0
id_anim=0
if(efanim/=0) then 
  call init(Ainv_hash)
  call init(Ainv_ija)
  call zerom(Ainv_hash,nanim)
  call zerom(Ainv_ija,nanim)
endif
if(efaSNP/=0) then
  allocate(Z(ndata,nSNP))
  Z=0d0
endif
if(efdSNP/=0) then
  allocate(W(ndata,nSNP))
  W=0d0
endif
sumsol=0d0
ssqsol=0d0
beginouter=1
avinmodel=0d0 
sumpa=0d0
sumpd=0d0

if(efanim/=0) then
  ! set up A inverse
  call add_g_add(3)
  print *,'a-1'
  close (7)
endif

! --------------
! read genotypes
! --------------
if( efasnp/=0 .or. efdsnp /=0) then 
    ! read file
    rewind(8)
    do i=1,ngeno
            read(8,informat)id_geno(i),genotype(i,:)
            if(mod(i,1000)==0) print *,i,'-th record'
    enddo
    print *,'read genotypes'
endif

!set up design matrix for f.e. , SNPs come later
pos=0
createX: do i=1,ndata
!  if((efaSNP/=0).or.(efdSNP/=0)) then
!    read(1,trim(informat),iostat=io)in_data,genome(1:nSNP*2)
!  else
  read(1,*,iostat=io)in_data
!  endif
  if(io/=0) exit
  if(mod(i,1000)==0) print *,i,'-th record'
  y(i)=in_data(traitcol)
  phen(i)=y(i)
  if(hetvarcol/=0) then 
    if(BinaryTrait) then
            print *,'weights are not compatible with BinaryTrait'
        print *,'stop'
        print *,'-----------'
        stop
    endif
    weights(i)=sqrt(in_data(hetvarcol)) !sqrt(edc) see blup_lasso
  else
    weights(i)=1d0
  endif
  if (phen(i)==MisVal .or. weights(i)==0) then
  	weights(i)=1d-50 ! can't put 0 because later we divide by weights
	nmis=nmis+1
	y(i)=0d0
	phen(i)=0d0
  endif
  if(efanim/=0) then
    id_anim(i)=in_data(pos_eff(efanim))
  endif
  id(i)=in_data(record_id)
  do j=1,neff
    val=0d0
    !regular effects
    if(type_eff(j)==cov) then
      pos=pos+1
      if(inmodel(j)) val=in_data(pos_eff(j))
      pos1=add(j,1)      
      call add_X(cov=val,rec=i,eq=pos1,pos=pos)
    elseif( (type_eff(j)==cross).or. & 
            (type_eff(j)==add_animal).or. &
            (type_eff(j)==perm_diagonal) ) then
      pos=pos+1
      if(inmodel(j)) val=1d0
      pos1=add(j,nint(in_data(pos_eff(j))))
      call add_X(cov=val,rec=i,eq=pos1,pos=pos)
!      print *,'read',in_data(pos_eff(j)),'effect',j
!      print *,val,i,pos1,pos
!      pause
    ! SNP effects  
    elseif( (type_eff(j)==add_SNP).or. & 
            (type_eff(j)==dom_SNP) ) then
      ! find out this animal in list of genotypes
      pos_anim=find_out(id(i),id_geno)
      if(pos_anim==0) then 
              print *,'anim',id(i),'not found in genotypes'
        stop
      endif
      
      !SNPoneanim=splittoint(genome(1:nSNP*2),nSNP)
      SNPoneanim=genotype(pos_anim,:)

      do k=1,nSNP
        if (type_eff(j)==add_SNP)then
          ! create 'a' effect
          val=0d0
          select case (SNPoneanim(k)) 
              case(0) 
                val=-1d0
              case(1)
                val=0.d0
              case(2)
                      val=1d0
              case default
                val=-200d0 ! this allows to identify missing values
                !print *,'neither 1 nor 2 --> assumed missing',id(i),'animal',k,'SNP'
          end select
          if(.not.inmodel(j)) val=0d0
          pos1=k
          Z(i,pos1)=val
          !Z(pos_anim,pos1)=val WRONG because the ggod position is i
        endif
        if (type_eff(j)==dom_SNP) then
        ! create d effect if both SNPs differ
          select case (SNPoneanim(k)) 
              case(0) 
                val=0d0
              case(1)
                val=1.d0
              case(2)
                      val=0d0
              case default
                val=-200d0 ! this allows to identify missing values
                !print *,'neither 1 nor 2 --> assumed missing',id(i),'animal',k,'SNP'
          end select
          if(.not.inmodel(j)) val=0d0
          pos1=k
          W(i,pos1)=val
          !W(pos_anim,pos1)=val
        endif
      enddo
    endif
  enddo
enddo  createX

print *,'number of missing values for phenotype: ',nmis

! missing values are fixed to the average value of its column
if(efaSNP/=0) then
  do i = 1,nsnp
    val=sum(Z(:,i), mask=(Z(:,i) > -100d0) )/count(Z(:,i) > -100d0) !mean of non-missing values
    where (Z(:,i) < -100d0)
      Z(:,i)=val
    end where
  enddo
  if(nsnp<10) call printmat(Z)
endif     
! for dominance
if(efdSNP/=0) then
  do i = 1,nsnp
    val=sum(W(:,i), mask=(W(:,i) > -100d0) )/count(W(:,i) > -100d0) !mean of non-missing values
    where (W(:,i) < -100d0)
      W(:,i)=val
    end where
  enddo
endif     


sumpq=0d0
freq=0d0
do i=1,nsnp
        freq(i)= &
        ! freq of the i-th snp
        ! count of "A" (0-> aa 1->  aA 2-> AA)
        ( 2d0*count(genotype(:,i)==2)+count(genotype(:,i)==1) ) / &
        ! count of genotypes (2*genotyped individuals at that locus)
        (2d0*count(genotype(:,i)/=5))
enddo
open(unit=1001,file='freq',status="replace")
do i=1,nsnp
	write(1001,*) freq(i)
enddo
sumpq=sum(freq*(1d0-freq))

print *,'sumpq=',sumpq,'meanp= ',mean(freq),'varp=',compute_var(freq)
print *,'Number of heterozygous SNPs: ',count( freq*(1d0-freq) /=0 )
deallocate(genotype)

call increase_X(pos) !trim X
call sort_X() !sort equation-wise
call index_X() !set up index to 1st eq term
close(1)



! model stuff
call setup_inmodel() !prepare logical vector defining if eq i is in the model or not.


! call simulation
if (simulation) then
  print *,'not fully implemented yet'
  call simulate()
  !stop
endif  


! prediction
if(predict) then
  !read solutions
  read(4,*) !skip title
  do i=1,neq
    read(4,*)pos1,pos2,sol(i),ssqsol(i),avinmodel(i)
  enddo  
  call predict_y() 
  !write it out
  write(12,*) 'id true prediction'
  do i=1,ndata
    write(12,*)id(i),y(i),ywiggle(i)
  enddo  
  print *,'predictions written'
  !compute correlation  
  if(compute_var(y)/=0) print *,'correlation ',compute_corr(y,ywiggle),'vary vary^', compute_var(y),compute_var(ywiggle)

  call compute_EBV() 
  print *,'--prediction finished, end of program!--'
  stop !!end of program here

endif  


call transform_X(weights,'multiply')
call transform_yZW(weights,'multiply')

! begin of computations

if(BinaryTrait) y=0d0
ywiggle=y
ywiggleold=y

!diag of LHS for sparse effects 
do i=1,neq
        pos1=X%ia_eq(i,1)
        pos2=X%ia_eq(i,2)
        if(pos1==0)cycle !for example for animal effects w/o record
        xpx(i)=dot_product(X%cov(X%ia_eq(i,1):X%ia_eq(i,2)),X%cov(X%ia_eq(i,1):X%ia_eq(i,2)))
        ! wrong?
enddo


! id for SNP effects
if(efaSNP/=0) then
        pos1=add(efaSNP,1)-1
        do i=1,nlev(efaSNP)
                xpx(pos1+i)=dot_product(Z(:,i),Z(:,i))
        enddo
endif
if(efdSNP/=0) then
        pos1=add(efdSNP,1)-1
        do i=1,nlev(efdSNP)
                xpx(pos1+i)=dot_product(W(:,i),W(:,i))
        enddo
endif
! write [y : X ]' for checking (disabled)

close(2)
if((.not.continuation).and.(VCE))write(3,*)'vara vard varg varp vare pa_1 pd_1 2varapqpi lambda2'

!gibbs sampler
if(continuation) then
        if(VCE.or.MCMCBLUP) call restart_iter()
        beginouter=iter/thin
        call compute_residuals()
        print *,'correction'
endif  


    ! set up inverses of the variance for random effects
D=0d0
tau2=0d0
if(efaSNP/=0) then
        D(add(efaSNP,1):add(efaSNP,nlev(efaSNP)))=1d0/vara
        tau2(add(efaSNP,1):add(efaSNP,nlev(efaSNP)))=vara
endif        
if(efdSNP/=0) D(add(efdSNP,1):add(efdSNP,nlev(efdSNP)))=1d0/vard
if(efp/=0)    D(add(efp,1):add(efp,nlev(efp)))=1d0/varp

if(weights_snp) then
        call read_weights_snp()
endif        

 
outer: do iterout=beginouter,niter/thin
  inner: do iterin=1,thin
    if (BinaryTrait) call predict_liability() ! sets y and ywiggle
    eps=0d0
    iter=(iterout-1)*thin+iterin
    
    ! set up D only if in model, otherwise 0
    ! set up inverses of the variance for random effects
    if(efaSNP/=0) then
        if(.not.BL2Var) then
                D(add(efaSNP,1):add(efaSNP,nlev(efaSNP)))=1d0/vara
        else
                D(add(efaSNP,1):add(efaSNP,nlev(efaSNP)))=&
                        1d0/tau2(add(efaSNP,1):add(efaSNP,nlev(efaSNP)))
        endif
    endif


    if(efdSNP/=0) D(add(efdSNP,1):add(efdSNP,nlev(efdSNP)))=1d0/vard
    if(efp/=0)    D(add(efp,1):add(efp,nlev(efp)))=1d0/varp
  
    if  (mod(iter,correction)==0) then
      call compute_residuals()
      print *,'correction'
    endif
    ! call predict_missing_y() ! in fact I predict missing residuals, is it useful?
                               ! not implemented yet
    
    loopeffects: do j=1,neff ! all effects 
      
      if(.not.inmodel(j)) then 
!        print *,'skipping ef',j
        cycle loopeffects !skip effects out of model
      endif
      do i=add(j,1),add(j,nlev(j)) 
        !form lhs
        lhs=xpx(i)/vare+D(i)
        if(type_eff(j)==add_SNP) then 
            rhs=dot_product(Z(:,i-offset(j)),ywiggle)/vare+xpx(i)/vare*sol(i)
        elseif(type_eff(j)==dom_SNP) then
            rhs=dot_product(W(:,i-offset(j)),ywiggle)/vare+xpx(i)/vare*sol(i)
        else 
            !sampling for sparse effects
            ! form rhs with y corrected by (iter-1) computation of b(i)
            rhs=0d0
            do k=X%ia_eq(i,1),X%ia_eq(i,2) !through all records affected by 'i'
              if(k==0) exit
              rhs=rhs+X%cov(k)*ywiggle( X%rec(k) )
            enddo
            rhs=rhs/vare+xpx(i)/vare*sol(i)
            if(type_eff(j)==add_animal) then
              pos1=add(efanim,1)
              lhs=lhs+diagAinv(i-offset(j))/varg
              pos2=add(efanim,nlev(efanim))
              rhs=update_gs_one_col(i-offset(j),Ainv_ija,varg,sol(pos1:pos2),rhs)
            endif
        endif
        ! sample solution from its conditional
        if(lhs==0) cycle
        if(.not.GaussSeidel) then
          val=normal(rhs/lhs,1d0/lhs)
          ! or do Gauss Seidel
        else
          val=rhs/lhs
        endif
        
        if(use_mixture) then
        ! note that this only works for MCMC
          if ((type_eff(j)==add_SNP)) then
           ! compute loglikelihood for state 1 (i -> in model) and 0 (not in model)
                ! Notes by RLF, p 47/67
                v1=xpx(i)*vare+(xpx(i)**2)*(1d0/D(i)) ! remember that D(i) is inverted
                v0=xpx(i)*vare
                rj=rhs*vare
                !state delta=0
                like2=log_like_normal((/rj/),v0)
                !state delta=1
                like1=log_like_normal((/rj/),v1)


             ! old, incorrect version
             ! like1=log_like_normal( ywiggle - Z(:,i-offset(j))*(val-sol(i))  ,vare )
           ! compute loglikelihood for state 2 (i -> not in model)
             ! old, incorrect version
             ! like2=log_like_normal( ywiggle + Z(:,i-offset(j))*sol(i),vare)
           ! I need to use loglikelihoods to avoid overflows
            aposteriori=(/like1,like2/)+log(pa) 
          elseif ((type_eff(j)==dom_SNP)) then
           ! compute loglikelihood for state 1 (i -> in model) and 0 (not in model)
                ! Notes by RLF, p 47/67
                v1=xpx(i)*vare+(xpx(i)**2)*D(i)
                v0=xpx(i)*vare
                rj=rhs*vare
                !state delta=0
                like2=log_like_normal((/rj/),v0)
                !state delta=1
                like1=log_like_normal((/rj/),v1)
           !  like1=log_like_normal( ywiggle - W(:,i-offset(j))*(val-sol(i))  ,vare )
           ! compute loglikelihood for state 2 (i -> not in model)
           !  like2=log_like_normal( ywiggle + W(:,i-offset(j))*sol(i),vare)
           ! I need to use loglikelihoods to avoid overflows
            aposteriori=(/like1,like2/)+log(pd) 
          endif
          if ((type_eff(j)==add_SNP).or.&
              (type_eff(j)==dom_SNP)) then
            ! a posteriori is in the log scale
            ! so I use a function to pass it to 0-1
            i_inmodel(i)=importance((/.true.,.false./),log2p(aposteriori))
            if(.not.i_inmodel(i))then 
              val=0d0
            endif
          endif
        endif

        !convergence criteria for Gauss Seidel
        eps=eps+(val-sol(i))**2
        ! correct y for that effect
        if(type_eff(j)==add_SNP) then 
          ywiggle=ywiggle - Z(:,i-offset(j))*(val-sol(i))  
        elseif(type_eff(j)==dom_SNP) then
          ywiggle=ywiggle - W(:,i-offset(j))*(val-sol(i))  
        else 
          do k=X%ia_eq(i,1),X%ia_eq(i,2) !through all records affected by 'i'
            if(k==0) exit !unknowns not in data
            ywiggle(X%rec(k))=ywiggle(X%rec(k)) - X%cov(k)*(val-sol(i))
          enddo
        endif
        !update sol
        sol(i)=val
      enddo
    enddo loopeffects
    eps=eps/sum(sol**2)
    ! new,different eps based on residuals to avoid scale sizes
    eps=sum((ywiggle-ywiggleold)**2)/sum(ywiggle**2)
    !  print *,iter,'eps',eps
    ywiggleold=ywiggle
         
    if(efaSNP/=0) includeda=count(i_inmodel(add(efaSNP,1):add(efaSNP,nlev(efaSNP))))
    if(efdSNP/=0) includedd=count(i_inmodel(add(efdSNP,1):add(efdSNP,nlev(efdSNP))))
            
    if ((eps<conv_crit).and.(GaussSeidel)) then
      avinmodel=1d0
      print *,iter,'eps',eps
      exit outer
    endif
    
    if(VCE) then
      ! sample vars
      ! form sums of squares and expectations (Rao-Blackwell estimates)
      if(efaSNP/=0) then
        if(.not.BL2Var) then
                if(dfvara<1000) then
                  pos1=add(efaSNP,1)
                  pos2=add(efaSNP,nlev(efaSNP))
                  vara=(sum(sol(pos1:pos2)**2) + vara_ap*dfvara)/ &
                  (dble(includeda) + dfvara) 
                  vara=chin(vara,includeda+dfvara)
                  tau2(pos1:pos2)=vara ! same variance for each SNP
                endif
        else
                ! sample individual variances for each SNP (additive SNPs only)
                pos1=add(efaSNP,1)
                pos2=add(efaSNP,nlev(efaSNP))
                do pos =pos1,pos2
                        !pos=add(efaSNP,i)
                        tau2(pos)=1d0/rinvGauss(mu=sqrt(lambda2/(sol(pos)**2)),lambda=lambda2)
                enddo
                ! sample lambda
                pos1=add(efaSNP,1)
                pos2=add(efaSNP,nlev(efaSNP))
                lambda2= random_gamma( &
                real(nlev(efaSNP),r8), & !shape
                2d0/sum(tau2(pos1:pos2)) ) !scale= 1/rate
                vara=2d0/lambda2
        endif
      endif

      if(efdSNP/=0) then
        if(dfvard<1000) then
          vard=(sum(sol(add(efdSNP,1):add(efdSNP,nlev(efdSNP)))**2) + vard_ap*dfvard)/ &
          (dble(includedd) + dfvard)
          vard=chin(vard,dble(includedd)+dfvard)
        endif
      endif
      if(efanim/=0) then
        if(dfvarg<1000) then
          pos1=add(efanim,1)
          pos2=add(efanim,nlev(efanim))
          varg=(quadrf(sol(pos1:pos2),Ainv_ija,sol(pos1:pos2)) + varg_ap*dfvarg)/ &
          (nlev(efanim) + dfvarg)
          varg=chin(varg,dble(nlev(efanim))+dfvarg)
        endif
      endif

      if(efp/=0) then
        if(dfvarp<1000) then
          pos1=add(efp,1)
          pos2=add(efp,nlev(efp))
          varp=(sum(sol(add(efp,1):add(efp,nlev(efp)))**2) + varp_ap*dfvarp)/ &
          (nlev(efp) + dfvarp)
          varp=chin(varp,dble(nlev(efp))+dfvarp)
        endif
      endif

      if(dfvare<1000) then
!        vare=(sum(ywiggle**2) + vare_ap*dfvare )/ &
!        (ndata-nmis + dfvare) 
        ! modified to account for missing values, those with weights<1d-50
        vare=sum(ywiggle**2,mask= weights>1d-50 )/&
        (ndata-nmis + dfvare) 
        
        vare=chin(vare,dble(ndata-nmis)+dfvare)
      endif 
      
      ! sampling proportions
      if(use_mixture) then
              if(.true.) then
                if(efaSNP/=0) then
                  !print *,' sampling p from ~Beta',includeda,'+',aprioria(1),'/',nlev(efaSNP)-includeda,'+',aprioria(2)
                  pa(1)=random_beta(includeda+aprioria(1),nlev(efaSNP)-includeda+aprioria(2))
                  pa(2)=1d0-pa(1)
                endif
                if(efdSNP/=0) then
                  !print *,' sampling p from ~Beta',includedd,'+',apriorid(1),'/',nlev(efdSNP)-includedd,'+',apriorid(2)
                  pd(1)=random_beta(includedd+apriorid(1),nlev(efdSNP)-includedd+apriorid(2))
                  pd(2)=1d0-pd(1)
                endif
              endif
      endif
    endif

    if(.not.GaussSeidel) then
      if(iter>burnin)then
       svara=svara+vara
       svard=svard+vard
       svare=svare+vare
       svarp=svarp+varp
       svarg=svarp+varp
       slambda2=slambda2+lambda2
       ssvara=svara+vara**2
       ssvard=svard+vard**2
       ssvare=svare+vare**2
       ssvarp=svarp+varp**2
       ssvarg=svarg+varg**2
       sslambda2=slambda2+lambda2**2
       sumsol=sumsol+sol
       ssqsol=ssqsol+sol**2
       where(i_inmodel) avinmodel=avinmodel+1d0    
       sumpa=sumpa+pa
       sumpd=sumpd+pd
       sumtau2=sumtau2+tau2
       sstau2=sstau2+tau2**2
       if(lBF) then
                ! covBF(i,j,k) contains, at the i marker, the covariance of markers i+j, i+k 
       		do i=firstBF,lastBF-(nBF-1)
			do j=0,nBF-1
				do k=0,nBF-1
					covBF(i,j,k)=covBF(i,j,k)+sol(add(efaSNP,i+j))*sol(add(efaSNP,i+k))
				enddo
			enddo
		enddo
       endif
      endif
    endif
    
    if( (VCE.or.MCMCBLUP) .and. (mod(iter,saving)==0) ) call save_iter()
  enddo inner
  if(GaussSeidel) print *,'eps: ',eps
   
if(VCE)  write(3,'(30g14.5)') vara,vard,varg,varp,vare,pa(1),pd(1),2d0*vara*sumpq*pa(1),lambda2
    print *,iter,'ef 1 to 3',sol(1:3),'vara,vard,varg,varp,vare,pa(1),pd(1),includeda, 2varapqpi lambda2', &
        vara,vard,varg,varp,vare,pa(1),pd(1),includeda,2d0*sumpq*vara*pa(1),lambda2
  call printtime()

enddo outer

if(.not.GaussSeidel) then 
  sol=sumsol/(niter-burnin) 
  !ssqsol now contains the posterior sd. err. of the solutions
  ssqsol=sqrt( ssqsol/(niter-burnin) - sol**2 )
  !sstau2 now contains the posterior sd. err. of the tau2
  tau2=sumtau2/(niter-burnin)
  sstau2=sqrt( sstau2/(niter-burnin) - tau2**2 )
  !posterior means of variance components
  vara=svara/(niter-burnin)
  vard=svard/(niter-burnin)
  varg=svarg/(niter-burnin)
  vare=svare/(niter-burnin)
  varp=svarp/(niter-burnin)
  lambda2=slambda2/(niter-burnin)
  pa=sumpa/(niter-burnin)
  pd=sumpd/(niter-burnin)
  ! and s.d.
  sslambda2=sslambda2/(niter-burnin)-lambda2**2
  ssvara=ssvara/(niter-burnin)-vara**2
  ssvard=ssvard/(niter-burnin)-vard**2
  ssvare=ssvare/(niter-burnin)-vare**2
  ssvarp=ssvarp/(niter-burnin)-varp**2
  ssvarg=ssvarg/(niter-burnin)-varg**2
  open(unit=11,file=trim(adjustl(parfile))//'_finalEstimates',status='replace')
  write (11,*) 'vara vard varp varg vare lambda2 pa(1) pd(1)'
  write (11,*) 'posterior means'
  write (11,*) vara, vard, varp, varg, vare, lambda2, pa(1), pd(1)
  write (11,*) 'posterior s d'
  write (11,*) sqrt(ssvara), sqrt(ssvard), sqrt(ssvarp), sqrt(ssvarg), sqrt(ssvare), sqrt(sslambda2)
  close(11)

  write (*,*) '--------------------------------'
  write (*,*) 'Final Estimates'
  write (*,*) 'vara vard varp varg vare lambda2 pa(1) pd(1)'
  write (*,*) 'posterior means'
  write (*,*) vara, vard, varp, varg, vare, lambda2, pa(1), pd(1)
  write (*,*) 'posterior s d'
  write (*,*) sqrt(ssvara), sqrt(ssvard), sqrt(ssvarp), sqrt(ssvarg), sqrt(ssvare), sqrt(sslambda2)
  write (*,*) '--------------------------------'

  
  avinmodel=avinmodel/(niter-burnin)
  if(lBF) then
  	! compute average cross-product  
  	covBF=covBF/(niter-burnin)
	! covBF(i,j,k) contains, at the i marker, the covariance of markers i+j, i+k 
  	! compute covariances across loci
	do i=firstBF,lastBF-(nBF-1)
		do j=0,nBF-1
			do k=0,nBF-1
				covBF(i,j,k)=covBF(i,j,k)  -sol(add(efaSNP,i+j))*sol(add(efaSNP,i+k))
			enddo
		enddo
	enddo
	call store_BF()
  endif


endif

call store_solutions()

call transform_X(weights,"divide")
call transform_yZW(weights,"divide")
call compute_EBV()

contains 
                      
function splittoint(genome,nSNP) result(SNP)
! converts one chain of SNP alleles 1211121112 in a matrix of integers
!                                   1 1 1 1 1
!                                   2 1 2 1 2  (transposed)
implicit none
character(len=*) genome
integer:: nSNP,SNP(nSNP,2),aux(nSNP*2)
logical,save:: first=.true.,free=.false.
character(len=20),save:: fmat
! I can put here code to guess if genome has free or condensed format
!does not worj because if this is free format then genome='1       '....
if (first) then !build format
        print *,'first 10 loci: ',genome(1:20)
    write(fmat,'(i10)')nSNP*2
    fmat=trim('('//adjustl(fmat))//'i1)'
    print *,'format= ',fmat
    first=.false.
endif
read(genome,fmt=fmat)aux
SNP=reshape(aux,shape(SNP),order=(/2,1/))

end function

  subroutine add_g_add(type)
  ! IMisztal subroutine heavily modified to construct A-1 (without genetic variance!!)
! generates contributions for additive sire or animal effect
  integer :: type,k,l,m,n,io,animal,sire,dam,par_stat,p(3)
  real(r8) ::w(3),res(4),val
!
  select case (type)
     case (3) !animal model
        w=(/1., -.5, -.5/)
        res=(/2., 4/3., 1., 0./)
     case (4) ! sire model
        w=(/1., -.5, -.25/)
        res=(/16/11., 4/3., 16/15., 1./)
  end select

  do
     read(7,*,iostat=io) animal, sire, dam    
     if (io /= 0) exit
     p(1)=animal
     p(2)=sire
     p(3)=dam
     par_stat=1
     if(sire==0) par_stat=par_stat+1
     if(dam==0) par_stat=par_stat+1
     
     do k=1,3
        do l=1,3        
           if (p(k) /= 0 .and. p(l) /=0) then
             m=p(k)
             n=p(l)
             val=w(k)*w(l)*res(par_stat)
             if(m==n) diagAinv(m)=diagAinv(m)+val
!             call addm(val,m,n,Ainv_hash,'f') !this stocks the matrix in full !bug? did not work well
             call addm(val,m,n,Ainv_hash) 
           endif
        enddo
     enddo
  enddo
  Ainv_ija=ainv_hash
call convert_hash_ija_general(ainv_ija,ainv_hash,conv_upper_full)  !convert to full-stored Ainv matrix
!  call convert_hash_ija_full(ainv_ija,ainv_hash) 
  call reset(Ainv_hash)
  
  print *,'n'
  print *,ainv_ija%n
  print *,'nel'
  print *,ainv_ija%nel
  print *,'ia(1:3)'
  print *,ainv_ija%ia(1:3)
  print *,'ja(1st col)'
  !1st col
  print *,ainv_ija%ja(ainv_ija%ia(1):ainv_ija%ia(2)-1)
  
  end subroutine

function update_gs_one_col(whichcol,Ainv_ija,varg,sol,rhs) result(rhsupdated)
!corrects rhs for the whichcol element in sol by the other correlated effects
! this is, one round in one element of GaussSeidel
implicit none
type(sparse_ija) :: Ainv_ija
integer::whichcol,i
real(r8):: rhs,sol(:),rhsupdated,varg

!if(Ainv_ija%n /= size(sol)) then
!  print *,Ainv_ija%n,size(sol),'Ainv_ija%n /= size(sol)'
!  stop
!endif  
!if(whichcol>Ainv_ija%n) then
!  print *,'whichcol>Ainv_ija%n',whichcol,Ainv_ija%n
!  stop
!endif

rhsupdated=rhs
do i=Ainv_ija%ia(whichcol),Ainv_ija%ia(whichcol+1)-1

!  print *,whichcol,i,Ainv_ija%a(i),Ainv_ija%ja(i)

  if(Ainv_ija%ja(i)/=whichcol) &
    rhsupdated=rhsupdated-Ainv_ija%a(i)/varg*sol(Ainv_ija%ja(i))
enddo  

end function

subroutine printtime
INTEGER  DATE_TIME (8)
CHARACTER(LEN = 12) REAL_CLOCK (3)

CALL DATE_AND_TIME (REAL_CLOCK (1), REAL_CLOCK (2), &
                REAL_CLOCK (3), DATE_TIME)

print *,real_clock(1)(7:8),'/',real_clock(1)(5:6),&
'/',real_clock(1)(1:4),' ', &
real_clock(2)(1:2),':',real_clock(2)(3:4),':',real_clock(2)(5:6)

end subroutine

subroutine read_par_file()
implicit none
integer :: ni
character(len=200):: mes(3)
logical:: nomean


read(10,*) 
read(10,'(a)') datafile; print*,datafile
read(10,*)
read(10,'(a)') pedfile; print*,pedfile
read(10,*)
read(10,'(a)') genofile; print*,genofile
read(10,*)
read(10,*) nSNP
read(10,*)
read(10,*) method
read(10,*)
read(10,*) simulation
read(10,*) 
read(10,*)
read(10,*) niter
read(10,*)
read(10,*) burnin
read(10,*) 
read(10,*) thin
read(10,*)
read(10,*) conv_crit
read(10,*)
read(10,*) correction
read(10,*)
read(10,'(a)') varfile
read(10,*)
read(10,'(a)') solfile
read(10,*)
read(10,*) traitcol,hetvarcol
read(10,*)
read(10,*) neff 
allocate(pos_eff(neff),type_eff(neff),nlev(neff))
read(10,*)

pos_eff=0
type_eff=0
do i=1,neff
  read(10,*)pos_eff(i),type_this_eff,nlev(i)
! guess kind of effects
  do j=1,size(type_eff_char)
    if (type_this_eff==type_eff_char(j)) type_eff(i)=j
    if(type_eff(i)==add_SNP) then
        efaSNP=i
        if(pos_eff(i)/=0) then
                print *,'position of add_SNP set to 0 because it is in the genotype file'
                pos_eff(i)=0
        endif
    endif
    if(type_eff(i)==dom_SNP) then
            efdSNP=i
        if(pos_eff(i)/=0) then
                print *,'position of dom_SNP set to 0 because it is in the genotype file'
                pos_eff(i)=0
        endif
    endif
    if(type_eff(i)==add_animal) efanim=i
    if(type_eff(i)==perm_diagonal) efp=i
  enddo
enddo  

!read(10,*)
!read(10,'(a200)') informat
read(10,*)
read(10,*)
read(10,*) vara,dfvara 
read(10,*)
read(10,*) vard,dfvard
read(10,*)
read(10,*) varg,dfvarg
read(10,*)
read(10,*) varp,dfvarp
read(10,*)
read(10,*) vare,dfvare
read(10,*)
read(10,*) record_id
read(10,*)
read(10,*) continuation
!if( (any(type_eff==add_SNP)).or.(any(type_eff==dom_SNP)) ) then
!  allocate(in_data(maxval((/pos_eff-1,record_id,traitcol,hetvarcol/))))
!else  !no SNP effects
  allocate(in_data(maxval( (/pos_eff,record_id,traitcol,hetvarcol/))))
!endif

vara_ap=vara
vard_ap=vard
varg_ap=varg
varp_ap=varp
vare_ap=vare

allocate(offset(neff),&
         inmodel(neff))


read(10,*)
read(10,*) inmodel
if(.true.)then
        read(10,*)
        read(10,*) aprioria
        read(10,*)
        read(10,*) apriorid
        use_mixture=.false.
        read(10,*)
        read(10,*) use_mixture
        if(use_mixture) then
          print *,'---------------------------------'
        !  print *,'          WARNING                '
          print *,'           using BayesCPi '
        !  print *,'mixture not usable :-( , disabled'
        !  print *,'mixture not usable :-( , trying..'
          print *,'---------------------------------'
          !use_mixture=.false.
        else
                aprioria(1)=1d0
                aprioria(2)=0d0
                apriorid(1)=1d0
                apriorid(2)=0d0
        endif
endif
  

GaussSeidel=.false.; VCE=.false.; MCMCBLUP=.false.;predict=.false.
select case (adjustl(method))
  case('MCMCBLUP')
    MCMCBLUP=.true.
    if(use_mixture) then
      print *,'mixture only compatible with VCE or PREDICT, disabled'
      use_mixture=.false.
    endif
  case('BLUP')
    GaussSeidel=.true.
    if(use_mixture) then
      print *,'mixture only compatible with VCE or PREDICT, disabled'
      use_mixture=.false.
    endif
  case('VCE')
    VCE=.true.
  case('PREDICT')
    predict=.true.
  case default
    print *,'method: ',method,'not available'
    stop
end select
if(GaussSeidel.or.predict) continuation=.false.

if(VCE.or.MCMCBLUP) then
  if (dfvara==-2) then
     vara_ap=0d0
  endif  
  if (dfvard==-2) then
     vard_ap=0d0
  endif  
  if (dfvarg==-2) then
     varg_ap=0d0
  endif  
  if (dfvarp==-2) then
     varp_ap=0d0
  endif  
  if (dfvare==-2) then
     vare_ap=0d0
  endif  
endif

if(sum(abs(aprioria-0d0))<1d-50) then
  pa=0.5d0
else
  pa=dble(aprioria)/sum(aprioria)
endif
if(sum(abs(apriorid-0d0))<1d-50) then
  pd=0.5d0
else
  pd=dble(apriorid)/sum(apriorid)
endif

call getoption_unit(10,'BayesianLasso',n=ni,xc=mes)
if(ni>=1 .and. VCE) then
        if(mes(1)=='ParkCasella') BL1Var=.true.
        if(mes(1)=='Tibshirani') BL2Var=.true.
endif
call getoption_unit(10,'BinaryTrait',n=ni)
if(ni==0) then
        BinaryTrait=.true.
	correction=10000000
	vare=1d0
	vare_ap=1d0
	dfvare=1d10
else
        BinaryTrait=.false.
endif

call getoption_unit(10,'MissingValue',n=ni,x=aa)
if(ni>=1) then
	MisVal=aa(1)
	print *,'OPTION Missing value of y changed from -9999 to',MisVal
endif

call getoption_unit(10,'SNP_weights',n=ni,xc=mes,x=aa)
if(ni>=1) then
        weights_snp=.true.
        colweights=aa(2)
        fileweights=mes(1)
        print *,'OPTION Weights for SNPs in column ',colweights,' of file ',trim(adjustl(fileweights))
endif 

call getoption_unit(10,'BayesFactor',n=ni,x=aa)
if(ni>=1) then
	lBF=.true.
	nBF=aa(1)
	print *,'BayesFactors  to windows from 1 up to of n consecutive SNPs: ',nBF
	if(ni==3) then
		! compute from markers firstBF to marker lastBF
		firstBF=aa(2)
		lastBF=aa(3)
	else
		! the 0 also works as a flag
		firstBF=0
		lastBF=0
	endif
else
	lBF=.false.
	nBF=1
endif

call is_there_a_mean()

end subroutine

subroutine is_there_a_mean()
implicit none
logical:: continu=.false.
integer:: i
! we need a cross-classified effect and also inmodel

do i=1,neff
        if (type_eff(i)==cross) then
                continu=continu.or.inmodel(i)
        endif
enddo

if(.not.continu) then
        print *,'-----------------------------------'
        print *,'                                   '
        print *,'------      WARNING      ----------'
        print *,'                                   '
        print *,' your model does not contain any fixed effect'
        print *,' like an overall mean '
        print *,' are you certain  you want to continue? '
        print *,' (T/F ?)'
        print *,'                                   '
        print *,'------      WARNING      ----------'
        print *,'                                   '
        print *,'-----------------------------------'
        read *,continu
        if(.not. continu) stop
endif
end subroutine

subroutine print_version()

print *
print *,'-----------------------'
print *,'--        GS3        --'
print *,version
print *,'-- with Bayesian Lasso --'

print *,'-----------------------'
print *,'Copyright (c) A.Legarra '
print *,' A. Ricard, O. Filangi'
print *,'    INRA, FRANCE'
print *,date
print * 
print *,'    This program is distributed in the hope that it will be useful,'
print *,'    but WITHOUT ANY WARRANTY; without even the implied warranty of'
print *,'    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the'
print *,'    GNU General Public License for more details.'

print *,'    You should have received a copy of the GNU General Public License'
print *,'    along with this program.  If not, see <http://www.gnu.org/licenses/>.'
!print *                              
!print *,' 42.924959,  -1.818967' 
!print *,'-34.693109, -58.376423'
print *,'----------------------'
call printtime()
end subroutine

subroutine print_info()

!call print_version()

print *,'parameter file: ',parfile
print *,'data file: ',datafile
print *,'with: ',ndata,'records'
print *,'reading positions',pos_eff
print *,'the record id is in column',record_id
print *,'trait read in',traitcol,'with weight in col',hetvarcol

print *,'pedigree file: ',pedfile
print *,'with: ',nanim,'records read'
print *,'genotype file: ',genofile
print *,'with: ',ngeno,'records read'
if ( (ngeno<ndata) .and. ((efaSNP/=0).or.(efdSNP/=0)) ) then
        print *,'               |-------- WARNING ----------------------------|'
        print *,'               | more records than genotypes                 |'
        print *,'               ----------------------------------------------|'
        !stop
endif


print *,'model with',neff,' effects= '
do i=1,neff
  if(inmodel(i)) then
    select case(type_eff(i))
    case (add_SNP) 
      print *,'  -> additive SNP effect in position',pos_eff(i)
    case (dom_SNP)
      print *,'  -> dominant SNP effect in position',pos_eff(i)
    case (ind_SNP)
      print *,'  -> individual SNP effect in position',pos_eff(i)
    case (add_animal)
      print *,'  -> additive infinitesimal effect in position',pos_eff(i)
    case (cross)
      print *,'  -> generic cross-classified ''fixed'' effect in position',pos_eff(i)
    case (cov)
      print *,'  -> generic covariable (not nested) in position',pos_eff(i)
    case (perm_diagonal)
      print *,'  -> generic environmental random effect in position',pos_eff(i)
    case default
      print *,'not supported or incorrectly written',type_eff(i),'for effect',i
      stop
    end select 
    print *,'     with',nlev(i),'levels'
  endif
enddo  
print *,'for a total of',neq,' equations'
print *,'length(in_data)=',size(in_data)
print *,'reading format',trim(informat)
print *,'--------------------------'
if(simulation) print *,'SIMULATED DATA'
if(continuation) then
  print *,'----------------------------------------'
  print *,'Warning: continuation of a previous MCMC'
  print *,'----------------------------------------'
endif
print *,'method: ',method
print *,'variances: vara vard varg varp vare'
print *,vara,vard,varg,varp,vare
print *,'residual is updated (corrected) every',correction,'iterations'
print *,'saving for continuation every',saving,'iterations'
print *
if(GaussSeidel) then
  print *,'--Gauss Seidel parameters--'
  print *,'convergence criterion:',conv_crit
elseif(predict) then
  print *,'--predicting--'
  if(simulation) then
    print *,'predicting simulation   from solutions in ',trim(solfile),  ' to file ''predictions'' '
  else
    print *,'predicting ', trim(datafile),' from solutions in ',trim(solfile),  ' to file ''predictions'' '
  endif
else  
  print *,'--Gibbs sampling parameters--'
  print *,'total number of iterations:',niter
  print *,'discarded for sol estimates:',burnin
  print *,'thin interval:',thin
  print *,'seeds:', s1,s2,s3 !from module ecuyer_random
  if (dfvara==-2) then
    print *,'--> flat prior for vara <--'
  endif  
  if (dfvard==-2) then
    print *,'--> flat prior for vard <--'
  endif  
  if (dfvarg==-2) then
    print *,'--> flat prior for varg <--'
  endif  
  if (dfvarp==-2) then
    print *,'--> flat prior for varp <--'
  endif  
  if (dfvare==-2) then
    print *,'--> flat prior for vare <--'
  endif  
  
  print *,'a priori variances and df'
  print *,' vara ; vard; varg; varp; vare'
  print *,vara_ap,dfvara,';',vard_ap,dfvard,';',varg_ap,dfvarg,';',varp_ap,dfvarp,';',vare_ap,dfvare
  if(use_mixture) then
    print *
    print *,'-----------------------------'
    print *,'Mixture model for SNP effects'
    print *,'-----------------------------'
    print *,'apriori parameters for p(SNP~N(0,vara))'
    print *,aprioria
    if(sum(abs(aprioria-0d0))<1d-50) then
      print *,'--> flat priors for the mixture model <--'
      print *,'initial proportion set to 0.5'
    endif
    print *,'apriori parameters for p(SNP~N(0,vard))'
    print *,apriorid
    if(sum(abs(apriorid-0d0))<1d-50) then
      print *,'--> flat priors for the mixture model <--'
      print *,'initial proportion set to 0.5'
    endif
  endif
  if(BL1Var) then
          print *
        print *,'-----------------------------'
        print *,'Bayesian Lasso a la Park & Casella'
  endif
  if(BL2Var) then
          print *
        print *,'-----------------------------'
        print *,'Bayesian Lasso a la Tibshirani'
        lambda2=2d0/vara
        print *,'initial value for lambda**2',lambda2
  endif
  if(BinaryTrait) then
          print *
        print *,'-----------------------------'
        print *,'  Binary trait using threshold (probit) model'
        print *,'-----------------------------'
  endif


endif    
print *,'--------------------'
print *

end subroutine

  subroutine store_solutions
  !from blupf90
  ! store solutions in file 'solutions'
  integer e,l
  open(11,file=solfile,status='replace',recl=1000)
  write(11,'(''effect level  solution sderror p tau2 sdtau2'')')
  !write(11,'(''effect level  solution sderror p'')')
  do e=1,neff
    do l=1,nlev(e)
       write(11,'(i3,i10,30g16.8)')e,l, sol(offset(e)+l),ssqsol(offset(e)+l),avinmodel(offset(e)+l),&
         tau2(offset(e)+l),sstau2(offset(e)+l) 
       !write(11,'(i3,i10,2g16.8)')e,l, sol(offset(e)+l),avinmodel(offset(e)+l)
    enddo
  enddo
  close(11)
  print*,'solutions stored in file: ',solfile
  end subroutine

  subroutine store_BF
  integer i,j,k,l,off
  real(r8),allocatable:: a(:),S(:,:)
  real(r8):: logL1,logL0,BF
  open(11,file=trim(adjustl(parfile))//'_BF',status='replace',recl=1000)
  write(11,*)'i pos nBF BF logp0 logp1'
  do l=1,nBF
  	allocate(a(l),S(0:l-1,0:l-1))
  	! off is the approximate offset of the marker in the center of the haplotype
	! e.g., if nBF=4 -> off=2
	!          nBF=5 -> off=3
	! (I use the properties of the integer division in fortran)
  	off=nBF/2 + 1
	!p(0|0, I*vara)
  	a=0
  	S=0
	do j=0,l-1
		S(j,j)=1d0
	enddo
	S=S*vara
	logL0=log_like_normal_matrix(a,S)
  	do i=firstBF,lastBF-(nBF-1)
  		! p(0|a_hat, cov(a_hat))
  		a=sol(add(efaSNP,i):add(efaSNP,i+l))
  		S=0
		do j=0,l-1
			do k=0,l-1
				S(j,k)=covBF(i,j,k)
			enddo
		enddo
		logL1=log_like_normal_matrix(a,S)
		BF=logL0-logL1
		write(11,'(3i10,10f12.7)')i,i+off,l,BF,logL0,logL1
	enddo
	deallocate(a,S)
  enddo
  close(11)
  end subroutine

subroutine restart_iter()
! read from continuation file
implicit none
integer:: i,j
open(13,file=trim(parfile)//'_cont',status='old')
read(13,*) iter
read(13,*) s1,s2,s3
read(13,*) vara
read(13,*) vard
read(13,*) varg
read(13,*) varp
read(13,*) vare
read(13,*) pa
read(13,*) pd
read(13,*) sumpa
read(13,*) sumpd
read(13,*)
do i=1,neq
  read(13,*) j,sol(i),sumsol(i),ssqsol(i),avinmodel(i)
  if(j/=i) then
    print *,'error reading continuation',i,j
    stop
  endif
enddo
print *,'restart from iteration',iter
close(13)
end subroutine  

subroutine save_iter()
! save iterations
implicit none
integer:: i
open(13,file=trim(parfile)//'_cont',status='replace')
write(13,*) iter, 'iteration'
write(13,*) s1,s2,s3, 'seeds'
write(13,*) vara,'vara'
write(13,*) vard,'vard'
write(13,*) varg,'varg'
write(13,*) varp,'varp'
write(13,*) vare,'vare'
write(13,*) pa,'pa'
write(13,*) pd,'pd'
write(13,*) sumpa,'sumpa'
write(13,*) sumpd,'sumpd'
write(13,*) 'sol, sumsol, zscore'
do i=1,neq
  write(13,*) i,sol(i),sumsol(i),ssqsol(i),avinmodel(i)
enddo
close(13)
end subroutine  


subroutine setup_inmodel()
! for the aux variable i_inmodel fills in whether level i 
! of the equations is in model or not
implicit none
integer:: i,pos1,pos2
do i=1,neff
  pos1=add(i,1)
  pos2=add(i,nlev(i))
  i_inmodel(pos1:pos2)=inmodel(i)
enddo
end subroutine


function add(effect,level)
! address of level in effect
implicit none
integer:: effect,level,add
  add=offset(effect)+level
end function


subroutine simulate()
  ! simulation is done so that h2=0.99
  includeda=0
!  call get_clock_seeds()
  sol=0d0
  sol(1)=0d0
  sol(2)=0d0
  sol(3)=0d0
  pos1=add(efaSNP,1); pos2=add(efaSNP,nlev(efaSNP))
  kk: do i=pos1,pos2
    ! pi = 0.99
    if( importance((/.true.,.false./),(/.01d0,.99d0/)) ) then
      sol(i)=normal01()
      includeda=includeda+1
      print *,'includeda',includeda,i
    endif
  enddo kk
  ! rescale
  y=matmul(Z,sol(pos1:pos2))
  y=y*sqrt(99*vare)/(sqrt(sum(y**2)/ndata))
  print *,'var due to SNPs',sum(y**2)/ndata

  do i=1,ndata
    y(i)=100d0+y(i)+normal(0d0,vare)
  enddo
  print *,'vary',sum((y-100d0)**2)/ndata
  sol=0d0
  includeda=0
  
end subroutine

subroutine add_X(cov,rec,eq,pos)
implicit none
real(r8)::cov
integer:: rec,eq,pos
! fills in X matrix
if(pos>maxrec) then
  print *,'X too small, increasing'
  call increase_X()
endif
X%cov(pos)=cov  
X%rec(pos)=rec  
X%eq(pos)=eq  
end subroutine

subroutine increase_X(ideal)
! increase X design matrix structure by 50%
! AL 2/8/07
implicit none
type(design_matrix):: X_temp
integer:: newmaxrec,oldmaxrec
integer,optional:: ideal

! create temp storage
oldmaxrec=maxrec
allocate(X_temp%cov(oldmaxrec),&
         X_temp%rec(oldmaxrec),&
         X_temp%eq(oldmaxrec))

X_temp%cov=X%cov
X_temp%rec=X%rec
X_temp%eq=X%eq

if(present(ideal)) then 
  newmaxrec=ideal
else  
  newmaxrec=floor(1.5*oldmaxrec) !new size
endif  
print *,'changing size of X from: ', oldmaxrec, 'to: ',newmaxrec

deallocate(X%cov,&
         X%rec,&
         X%eq)

allocate(X%cov(newmaxrec),&
         X%rec(newmaxrec),&
         X%eq(newmaxrec) )
X%cov=0d0; X%rec=0; X%eq=0         

if(newmaxrec<oldmaxrec) then
  X%cov=X_temp%cov(1:newmaxrec)
  X%rec=X_temp%rec(1:newmaxrec)
  X%eq=X_temp%eq(1:newmaxrec)
else
  print *,'changing size of X from: ', oldmaxrec, 'to: ',newmaxrec
  X%cov(1:oldmaxrec)=X_temp%cov
  X%rec(1:oldmaxrec)=X_temp%rec
  X%eq(1:oldmaxrec)=X_temp%eq
endif

deallocate(X_temp%cov,&
         X_temp%rec,&
         X_temp%eq)
maxrec=newmaxrec
print *,'lengths in increase_X: ',size(X%rec),size(X%cov),size(X%eq)


end subroutine

subroutine sort_X()
! sort X so that neq is ordered in increasing order; this *should* increase speed
! first send key to be sorted
implicit none
integer:: key(maxrec),i
type(design_matrix):: temp

allocate(temp%cov(maxrec),temp%eq(maxrec),temp%rec(maxrec))
! get index
call qsortd(X%eq,key,maxrec)

do i=1,maxrec
  temp%cov(i)=X%cov(key(i))
  temp%eq(i)=X%eq(key(i))
  temp%rec(i)=X%rec(key(i))
enddo  
X%cov=temp%cov
X%rec=temp%rec
X%eq=temp%eq
  
deallocate(temp%cov,temp%eq,temp%rec)
end subroutine  

subroutine index_X()
! now index X%ia_eq(i,:) contains the row of X having the 1st & last appearance of equation i
! Note that there might be some 'holes' (for SNP effects or animals)
! Changed to return the first and *last* element
implicit none
integer:: i,oldeq

X%ia_eq(1,1)=1
oldeq=1

do i=2,maxrec
  if(X%eq(i)/=oldeq)then
    X%ia_eq(oldeq,2)=i-1 !level of eq before
    X%ia_eq(X%eq(i),1)=i !level of the eq
    oldeq=X%eq(i)
  endif  
enddo
X%ia_eq(oldeq,2)=maxrec
end subroutine  

! I need to note how this works for a simple example...

subroutine transform_X(weights,what)
! this subroutine transform X given some weights. 
! if what="multiply" then we multiply each row in X by the appropriate weight 
! for example sqrt(edc) to get a homogeneous system of equations
! if what=="divide" we make the transformation back
implicit none
integer:: i
real(r8):: weights(:)
real(r8):: ws(size(weights))
character(len=*)::what

select case(what)
case('multiply')
  if(X%weighted) then
    print *,'should not weight X if already is'
    stop
  endif
  X%weighted=.true.
  ws=weights
case('divide')
  if(.not.X%weighted) then
    print *,'should not divide X if already is'
    stop
  endif
  X%weighted=.false.
  ws=1d0/weights
end select

do i=1,ndata
  where (X%rec==i)
    X%cov=X%cov*ws(i)
  end where
enddo
print *,'transforming X -> ',what,', weighted =',X%weighted
end subroutine

subroutine transform_yZW(weights,what)
! this subroutine transform yWZ given some weights. 
! if what="multiply" then we multiply each row in X by the appropriate weight 
! for example sqrt(edc) to get a homogeneous system of equations
! if what=="divide" we make the transformation back
implicit none
integer:: i
real(r8):: weights(:)
real(r8):: ws(size(weights))
character(len=*):: what

select case(what)
case('multiply')
  if(yWZ_weighted) then
    print *,'should not weight yWZ if already is'
    stop
  endif
  yWZ_weighted=.true.
  ws=weights
case('divide')
  if(.not.yWZ_weighted) then
    print *,'should not divide yWZ if already is'
    stop
  endif
  yWZ_weighted=.false.
  ws=1d0/weights
end select

! now reparameterize Z W and y
do i=1,ndata
  y(i)=y(i)*ws(i)
  if(allocated(Z)) Z(i,:)=Z(i,:)*ws(i)
  if(allocated(W)) W(i,:)=W(i,:)*ws(i)
  ! this leaves X... 
enddo  
print *,'transforming yZW ->',what,'weighted =',yWZ_weighted

end subroutine


subroutine compute_residuals()
!compute residuals ywiggle from y and X, W, Z
! e=ywiggle=y-[X W Z] [b u a d]
implicit none
integer::i

ywiggle=y
! traverse X
do i=1,maxrec
  ywiggle(X%rec(i))=ywiggle(X%rec(i))-sol(X%eq(i))*X%cov(i)
enddo
! addSNP
if(efaSNP/=0) then
  pos1=add(efaSNP,1)
  pos2=add(efaSNP,nlev(efaSNP))
  ywiggle=ywiggle-matmul(Z,sol(pos1:pos2))
endif
! domSNP
if(efdSNP/=0) then
  pos1=add(efdSNP,1)
  pos2=add(efdSNP,nlev(efdSNP))
  ywiggle=ywiggle-matmul(W,sol(pos1:pos2))
endif
end subroutine

subroutine predict_y()
!predict y from sol and X, W, Z
! y^=[X W Z] [b u a d]
! filtered by inmodel logical matrix (so it can be 
! modified through the predict parameter file)
! y^ is stored in ywiggle
! 
! it works as follows (a bit twisted)
! e^ = y-Xb^
! y^ = y-e^ = y-(y-Xb^)=Xb^
!

implicit none
real(r8):: sol_temp(size(sol))
integer::i

sol_temp=sol
! now nullify those sols for which inmodel=F
do i =1,neff
  if(.not.inmodel(i)) then
          sol(add(i,1):add(i,nlev(i)))=0d0
  endif
enddo
call compute_residuals() !generating ywiggle
ywiggle=y-ywiggle
sol=sol_temp
end subroutine

subroutine predict_liability()
! it predicts y,  such that phen falls in the correct category,
! given y^=y-ywiggle
! afterwards, it generates ywiggle
! it works also for missing phen (==0)
real(r8):: t(size(ywiggle))
real(r8):: lb(0:2),ub(0:2),l,u
lb=(/-9d0,-9d0,0d0/)
ub=(/9d0,0d0,9d0/)
t=0d0
do i=1,ndata
	l=lb(int(phen(i)))
	u=ub(int(phen(i)))
	l=l-(y(i)-ywiggle(i))
	u=u-(y(i)-ywiggle(i))
	t(i)=gen_trunc_normal(l,u) ! actually we sample the residual
	! lb(phen(i)) < y^+ epsilon < ub(phen(i))
enddo
y=(y-ywiggle)+t
ywiggle=t ! -> y_new-y^_old 
!print *,count(y>0 .and.  phen==1),'incoherencies'
!print *,count(y<0 .and.  phen==2),'incoherencies2'
end subroutine




subroutine compute_EBV()
! this subroutine outputs ebv's for animals in data *only*
real(r8):: ebv_asnp(size(ywiggle))
real(r8):: ebv_dsnp(size(ywiggle))
real(r8):: ebv_anim(size(ywiggle))
!real(r8):: pheno(size(y))
if(yWZ_weighted .or. X%weighted) then
 print*,'yWZ_weighted .or. X%weighted',yWZ_weighted, X%weighted
 stop
endif


   ! addSNP
  ywiggle=0d0
  if(efaSNP/=0) then
     pos1=add(efaSNP,1)
     pos2=add(efaSNP,nlev(efaSNP))
     ywiggle=matmul(Z,sol(pos1:pos2))
  endif
  ebv_asnp=ywiggle
  ! domSNP
  ywiggle=0d0
  if(efdSNP/=0) then
    pos1=add(efdSNP,1)
    pos2=add(efdSNP,nlev(efdSNP))
    ywiggle=matmul(W,sol(pos1:pos2))
  endif
  ebv_dsnp=ywiggle
  ywiggle=0d0
  if(efanim/=0) then
    do i=1,ndata
      ywiggle(i)=sol(add(efanim,id_anim(i)))
    enddo
  endif
  ebv_anim=ywiggle

  write(14,*) 'id g_aSNP g_dSNP poly_anim g_overall'
  !write(14,*) 'id EBV_aSNP EBV_dSNP EBV_anim EBV_overall' ! these are not EBVs but genotypic values!!
  do i=1,ndata
!    write(14,*)id(i),pheno(i),ebv_asnp(i),ebv_dsnp(i),ebv_anim(i)
    write(14,'(i10,4g15.6)')id(i),ebv_asnp(i),ebv_dsnp(i),ebv_anim(i),ebv_asnp(i)+ebv_dsnp(i)+ebv_anim(i)
  enddo  
  print *,'EBV''s written in ',trim(parfile),'_EBVs'
  
end subroutine



SUBROUTINE qsortd(x,ind,n)
 
! Code converted using TO_F90 by Alan Miller
! Date: 2002-12-18  Time: 11:55:47

IMPLICIT NONE
INTEGER, PARAMETER  :: dp = SELECTED_REAL_KIND(12, 60)

!REAL (dp), INTENT(IN)  :: x(:)
integer, INTENT(IN)  :: x(:)
INTEGER, INTENT(OUT)   :: ind(:)
INTEGER, INTENT(IN)    :: n

!***************************************************************************

!                                                         ROBERT RENKA
!                                                 OAK RIDGE NATL. LAB.

!   THIS SUBROUTINE USES AN ORDER N*LOG(N) QUICK SORT TO SORT A REAL (dp)
! ARRAY X INTO INCREASING ORDER.  THE ALGORITHM IS AS FOLLOWS.  IND IS
! INITIALIZED TO THE ORDERED SEQUENCE OF INDICES 1,...,N, AND ALL INTERCHANGES
! ARE APPLIED TO IND.  X IS DIVIDED INTO TWO PORTIONS BY PICKING A CENTRAL
! ELEMENT T.  THE FIRST AND LAST ELEMENTS ARE COMPARED WITH T, AND
! INTERCHANGES ARE APPLIED AS NECESSARY SO THAT THE THREE VALUES ARE IN
! ASCENDING ORDER.  INTERCHANGES ARE THEN APPLIED SO THAT ALL ELEMENTS
! GREATER THAN T ARE IN THE UPPER PORTION OF THE ARRAY AND ALL ELEMENTS
! LESS THAN T ARE IN THE LOWER PORTION.  THE UPPER AND LOWER INDICES OF ONE
! OF THE PORTIONS ARE SAVED IN LOCAL ARRAYS, AND THE PROCESS IS REPEATED
! ITERATIVELY ON THE OTHER PORTION.  WHEN A PORTION IS COMPLETELY SORTED,
! THE PROCESS BEGINS AGAIN BY RETRIEVING THE INDICES BOUNDING ANOTHER
! UNSORTED PORTION.

! INPUT PARAMETERS -   N - LENGTH OF THE ARRAY X.

!                      X - VECTOR OF LENGTH N TO BE SORTED.

!                    IND - VECTOR OF LENGTH >= N.

! N AND X ARE NOT ALTERED BY THIS ROUTINE.

! OUTPUT PARAMETER - IND - SEQUENCE OF INDICES 1,...,N PERMUTED IN THE SAME
!                          FASHION AS X WOULD BE.  THUS, THE ORDERING ON
!                          X IS DEFINED BY Y(I) = X(IND(I)).

!*********************************************************************

! NOTE -- IU AND IL MUST BE DIMENSIONED >= LOG(N) WHERE LOG HAS BASE 2.

!*********************************************************************

INTEGER   :: iu(21), il(21)
INTEGER   :: m, i, j, k, l, ij, it, itt, indx
REAL      :: r
!REAL (dp) :: t
integer :: t

! LOCAL PARAMETERS -

! IU,IL =  TEMPORARY STORAGE FOR THE UPPER AND LOWER
!            INDICES OF PORTIONS OF THE ARRAY X
! M =      INDEX FOR IU AND IL
! I,J =    LOWER AND UPPER INDICES OF A PORTION OF X
! K,L =    INDICES IN THE RANGE I,...,J
! IJ =     RANDOMLY CHOSEN INDEX BETWEEN I AND J
! IT,ITT = TEMPORARY STORAGE FOR INTERCHANGES IN IND
! INDX =   TEMPORARY INDEX FOR X
! R =      PSEUDO RANDOM NUMBER FOR GENERATING IJ
! T =      CENTRAL ELEMENT OF X

IF (n <= 0) RETURN

! INITIALIZE IND, M, I, J, AND R

DO  i = 1, n
  ind(i) = i
END DO
m = 1
i = 1
j = n
r = .375

! TOP OF LOOP

20 IF (i >= j) GO TO 70
IF (r <= .5898437) THEN
  r = r + .0390625
ELSE
  r = r - .21875
END IF

! INITIALIZE K

30 k = i

! SELECT A CENTRAL ELEMENT OF X AND SAVE IT IN T

ij = i + r*(j-i)
it = ind(ij)
t = x(it)

! IF THE FIRST ELEMENT OF THE ARRAY IS GREATER THAN T,
!   INTERCHANGE IT WITH T

indx = ind(i)
IF (x(indx) > t) THEN
  ind(ij) = indx
  ind(i) = it
  it = indx
  t = x(it)
END IF

! INITIALIZE L

l = j

! IF THE LAST ELEMENT OF THE ARRAY IS LESS THAN T,
!   INTERCHANGE IT WITH T

indx = ind(j)
IF (x(indx) >= t) GO TO 50
ind(ij) = indx
ind(j) = it
it = indx
t = x(it)

! IF THE FIRST ELEMENT OF THE ARRAY IS GREATER THAN T,
!   INTERCHANGE IT WITH T

indx = ind(i)
IF (x(indx) <= t) GO TO 50
ind(ij) = indx
ind(i) = it
it = indx
t = x(it)
GO TO 50

! INTERCHANGE ELEMENTS K AND L

40 itt = ind(l)
ind(l) = ind(k)
ind(k) = itt

! FIND AN ELEMENT IN THE UPPER PART OF THE ARRAY WHICH IS
!   NOT LARGER THAN T

50 l = l - 1
indx = ind(l)
IF (x(indx) > t) GO TO 50

! FIND AN ELEMENT IN THE LOWER PART OF THE ARRAY WHCIH IS NOT SMALLER THAN T

60 k = k + 1
indx = ind(k)
IF (x(indx) < t) GO TO 60

! IF K <= L, INTERCHANGE ELEMENTS K AND L

IF (k <= l) GO TO 40

! SAVE THE UPPER AND LOWER SUBSCRIPTS OF THE PORTION OF THE
!   ARRAY YET TO BE SORTED

IF (l-i > j-k) THEN
  il(m) = i
  iu(m) = l
  i = k
  m = m + 1
  GO TO 80
END IF

il(m) = k
iu(m) = j
j = l
m = m + 1
GO TO 80

! BEGIN AGAIN ON ANOTHER UNSORTED PORTION OF THE ARRAY

70 m = m - 1
IF (m == 0) RETURN
i = il(m)
j = iu(m)

80 IF (j-i >= 11) GO TO 30
IF (i == 1) GO TO 20
i = i - 1

! SORT ELEMENTS I+1,...,J.  NOTE THAT 1 <= I < J AND J-I < 11.

90 i = i + 1
IF (i == j) GO TO 70
indx = ind(i+1)
t = x(indx)
it = indx
indx = ind(i)
IF (x(indx) <= t) GO TO 90
k = i

100 ind(k+1) = ind(k)
k = k - 1
indx = ind(k)
IF (t < x(indx)) GO TO 100

ind(k+1) = it
GO TO 90
END SUBROUTINE qsortd

subroutine read_weights_snp()
implicit none
! actually do nothing
continue

end subroutine


end                    



