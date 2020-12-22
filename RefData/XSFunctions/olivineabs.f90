! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! OLIVINEABS
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! XSPEC local model for olivine absorption edge structure
! Contains ISMdust cross-sections for absorption with Fe-K edge for
! olivine, obtained from Rogantini et al. 2018
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine olivineabs(ear, ne, param, ifl, photar)
!
! The main routine to call all subroutines
!
implicit none
integer,parameter :: num_param = 2
integer,parameter :: nemod=13047 !Number of elements for each cross section.
integer :: ne, ifl
double precision :: moliv, rshift, emod(nemod), coemod(nemod)
double precision :: bxs(nemod), bener(nemod)
double precision :: zfac
real :: ear(0:ne), param(num_param), photar(ne)
logical :: startup=.true.
character (len=40) version
version='0.1'
if(startup)then
   call xwrite(' ',10)
   call xwrite('OlivineAbs: High resolution model for silicate with Olivine Fe K'//version, 10)
   call xwrite('Continuum xsect from ISMdust silicate from Corrales+ 2016 (MNRAS, 458, 1345)',10)
   call xwrite('Optical constants come from Draine 2003 (ApJ, 598, 1026) except',10)
   call xwrite('Fe K absorption, which comes from Rogantini+ 2018 (A&A, 609, A22)', 10)
   call xwrite('WARNING: If used in conjunction with neutral metal absorption models',10)
   call xwrite('(e.g. TBabs, TBnew), be sure to change abundances',10)
   call xwrite(' ',10)
   call read_cross_sections_olivine(nemod,bxs,bener)
   startup=.false.
endif
! Model parameters
moliv = param(1)
rshift = param(2)
zfac = 1.d0/(1.d0+dble(rshift))

call absorption_olivine(moliv, zfac, emod, nemod, coemod,bxs,bener)
!
call map_to_grid_olivine(dble(ear),ne,emod,nemod,photar,coemod)
return
end subroutine olivineabs
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
subroutine read_cross_sections_olivine(nemod,bxs,bener)
!
! This routine reads cross sections and puts them on a given grid
!
! It uses the X-Spec local variable/dictionary value ISMDUSTROOT
! to locate the data file. If this is not set then it uses
! the setting of the local_dir parameter below (which should be
! changed to match the location of the data file). The path -
! however it is given - should point to the location that contains
! the dust extinction templates, i.e. /path/to/ismdust/edge_files/
!
implicit none
integer :: nemod, status
double precision :: bener(nemod), bxs(nemod)
character (*), parameter :: fileloc = 'ismdust_edge_data.fits'
character (*), parameter :: olivinereadchat = 'olivineabs: reading from '
character (len=255 + 22) :: filename2 ! len(fileloc)
character (len=255) :: ismdust_root = ''
character (len=len(olivinereadchat)+len(filename2)) :: chatmsg = ''
character (len=512) :: contxt
integer inunit,readwrite,blocksize
integer :: hdutype,colnum
integer :: felem=1, nulld=0
logical :: anynull
character (len=255) :: fgmstr, fgmodf
external :: fgmstr, fgmodf

! Where do we look for the data?
ismdust_root = trim(fgmstr('ISMDUSTROOT'))
if (ismdust_root .EQ. '') then
   ismdust_root = trim(fgmodf())
endif
! parameters to specify the opening process
status=0
readwrite=0
blocksize=1
filename2=trim(ismdust_root) // fileloc
chatmsg=olivinereadchat // filename2
call xwrite(chatmsg, 20)
! Get an unused Logical Unit Number to use to open the FITS file.
call ftgiou(inunit,status)
! Open the FITS file
call ftopen(inunit,filename2,readwrite,blocksize,status)
if ( status .NE. 0 ) then
   contxt = 'olivineabs: Failed to open '//trim(filename2)
   call xwrite(contxt,5)
   call ftfiou(-1, status)
   return
endif
! Move to the OLIVINE extension (the binary table)
hdutype = 2
call ftmnhd(inunit,hdutype,'OLIVINE',1,status)
if ( status .NE. 0 ) then
   contxt = 'olivineabs: Failed to find OLIVINE extension in '//trim(filename2)
   call xwrite(contxt,5)
endif

!Read in one energy grid (column 1)
colnum=1
call ftgcvd(inunit,colnum,1,felem,nemod,nulld,bener,anynull,status)

!Read in the absorption cross section information (column 3)
colnum=3
call ftgcvd(inunit,colnum,1,felem,nemod,nulld,bxs,anynull,status)

! Report on errors (done before closing the file in case the error
! comes from closing the file). Unfortunately the xspec API does not
! provide a way to signal an error to the calling code, so a screen
! message is used, using the same method used to report the model
! the first time it is used.
!
! This message could be displayed only once, but it is probably worth
! repeating each time it is used.
if (status .ne. 0) then
   call xwrite('olivineabs: unable to read cross sections from '//trim(filename2),5)
endif
! Close the file and free the unit number
call ftclos(inunit, status)
call ftfiou(-1, status)
end subroutine read_cross_sections_olivine
! ======================================= !
subroutine absorption_olivine(moliv, zfac, e1, nemod, coeff, bxs, bener)
!
! This is routine that calculates the optical depth given the column densities
! Finally returns the absorption coefficient exp(-tau)
!
implicit none
integer :: nemod
integer :: i
double precision :: moliv
double precision :: bener(nemod), bxs(nemod), e1(nemod)
double precision :: tau, coeff(nemod)
double precision :: zfac
real hphoto, gphoto
external hphoto, gphoto

! Calculates the optical depth and the extinction coefficient exp(-tau)
e1(1)=(bener(1)*zfac)/1.d3
do i=2,nemod
  e1(i)=(bener(i)*zfac)/1.d3
  tau=moliv * bxs(i)
  coeff(i)=dexp(-tau)
enddo

end subroutine absorption_olivine
! ======================================= !
subroutine map_to_grid_olivine(new_en,n_ne,old_en, o_ne, new_flu, old_flu)
! This routine maps to a given grid
implicit none
integer :: i, o_ne, n_ne, bmin, bmax, btemp
double precision :: new_en(0:n_ne)
double precision :: old_en(o_ne), old_flu(o_ne)
double precision :: etemp, s
real :: new_flu(n_ne)
do i=1,n_ne
  new_flu(i)=0.
  call dbinsrch_olivine(new_en(i-1),bmin,old_en,o_ne+1)
  call dbinsrch_olivine(new_en(i),bmax,old_en,o_ne+1)
  etemp = (new_en(i)+new_en(i-1))/2
  ! Linear interpolation
  if (bmin.eq.bmax) then
    if (new_en(i).le.old_en(1)) then
      s = real(old_flu(1))
    else if (new_en(i).gt.old_en(o_ne)) then
      s = real(old_flu(o_ne))
    else
      s = old_flu(bmin) + (old_flu(bmax+1)-old_flu(bmin)) * (etemp-old_en(bmin)) / (old_en(bmax+1)-old_en(bmin))
      endif
  else
    call dbinsrch_olivine(etemp,btemp,old_en,o_ne+1)
    s = old_flu(btemp) + (old_flu(btemp+1)-old_flu(btemp)) * (etemp-old_en(btemp)) / (old_en(btemp+1)-old_en(btemp))
    endif
  new_flu(i)=real(s)
  enddo
end subroutine map_to_grid_olivine
! ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine dbinsrch_olivine(e,k,ener,n)
!
! search for energy e in array ener(1:n) and return bin k for
! which ener(k).le.e.lt.ener(k+1)
! adapted from J. Wilms's tbvabs_new.f routine
!
implicit none
double precision :: e
integer :: n,k,klo,khi
double precision :: ener(n)
klo=1
khi=n
1 if (khi-klo.gt.1) then
   k=(khi+klo)/2
   if(ener(k).gt.e)then
      khi=k
   else
      klo=k
   endif
   goto 1
endif
if (klo.eq.0) then
   call xwrite('Energy out of bounds. Should not happen',5)
   return
endif
k=klo
end subroutine dbinsrch_olivine
