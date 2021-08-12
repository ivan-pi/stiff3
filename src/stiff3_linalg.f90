! Copyright 2021 Ivan Pribec
! SPDX-License-Identifier: Apache-2.0

module stiff3_linalg

  implicit none
  private

  public :: lu, back

  !> Program for decomposing a matrix A to a lower and an upper
  !> triangular form $ A = LU$.
  interface lu
    module procedure lu_sp, lu_dp
  end interface


  !> Back substitution algorithm for solution to $LUx = b$.
  !> Used together with procedure `lu`, to calculate the solution
  !> of a system of linear equations
  interface back
    module procedure back_sp, back_dp, back1_sp, back1_dp
  end interface

  integer, parameter :: sp = kind(1.0e0)
  integer, parameter :: dp = kind(1.0d0)

  interface lapack_getrf
    pure subroutine sgetrf(m,n,a,lda,ipiv,info)
      import :: sp
      real(sp), intent(inout) :: a(lda,*)
      integer, intent(out) :: ipiv(*)
      integer, intent(out) :: info
      integer, intent(in) :: m
      integer, intent(in) :: n
      integer, intent(in) :: lda
    end subroutine sgetrf
    pure subroutine dgetrf(m,n,a,lda,ipiv,info)
      import :: dp
      real(dp), intent(inout) :: a(lda,*)
      integer, intent(out) :: ipiv(*)
      integer, intent(out) :: info
      integer, intent(in) :: m
      integer, intent(in) :: n
      integer, intent(in) :: lda
    end subroutine dgetrf
  end interface lapack_getrf

  interface lapack_getrs
    pure subroutine sgetrs(trans,n,nrhs,a,lda,ipiv,b,ldb,info)
      import :: sp
      real(sp), intent(in) :: a(lda,*)
      integer, intent(in) :: ipiv(*)
      real(sp), intent(inout) :: b(ldb,*)
      character(len=1), intent(in) :: trans
      integer, intent(out) :: info
      integer, intent(in) :: n
      integer, intent(in) :: nrhs
      integer, intent(in) :: lda
      integer, intent(in) :: ldb
    end subroutine sgetrs
    pure subroutine dgetrs(trans,n,nrhs,a,lda,ipiv,b,ldb,info)
      import :: dp
      real(dp), intent(in) :: a(lda,*)
      integer, intent(in) :: ipiv(*)
      real(dp), intent(inout) :: b(ldb,*)
      character(len=1), intent(in) :: trans
      integer, intent(out) :: info
      integer, intent(in) :: n
      integer, intent(in) :: nrhs
      integer, intent(in) :: lda
      integer, intent(in) :: ldb
    end subroutine dgetrs
  end interface lapack_getrs

contains

  subroutine lu_sp(amat,ipiv,info)
    real(sp), intent(inout) :: amat(:,:)
      !! On input a matrix, on output the LU factored matrices
    integer, intent(out) :: ipiv(:)
      !! Pivot array
    integer, intent(out), optional :: info
      !! Exit status
    integer :: m, n, lda, stat
    lda = max(1,size(amat,1))
    m = size(amat,1)
    n = size(amat,2)
    call lapack_getrf(m,n,amat,lda,ipiv,stat)
    call handle_error(info,stat,"getrf")
  end subroutine lu_sp

  subroutine lu_dp(amat,ipiv,info)
    real(dp), intent(inout) :: amat(:,:)
      !! LU-factored form of the original matrix A
    integer, intent(out) :: ipiv(:)
      !! Pivot array
    integer, intent(out), optional :: info
      !! Exit status
    integer :: m, n, lda, stat
    lda = max(1,size(amat,1))
    m = size(amat,1)
    n = size(amat,2)
    call lapack_getrf(m,n,amat,lda,ipiv,stat)
    call handle_error(info,stat,"getrf")
  end subroutine lu_dp

  subroutine back1_sp(amat,bvec,ipiv,info,trans)
    real(sp), intent(in) :: amat(:,:)
      !! LU-factored form of the original matrix A
    real(sp), intent(inout), target :: bvec(:)
      !! On input contains the right-hand side of a system of
      !! equations Ax = b, on output contains the solution vector x
    integer, intent(in) :: ipiv(:)
      !! Pivot array
    integer, intent(out), optional :: info
      !! Exit status
    character(len=1), intent(in), optional :: trans

    real(sp), pointer :: bptr(:,:)

    bptr(1:size(bvec),1:1) => bvec
    call back(amat,bptr,ipiv,info,trans)
  end subroutine back1_sp

  subroutine back_sp(amat,bmat,ipiv,info,trans)
    real(sp), intent(in) :: amat(:,:)
      !! LU-factored form of the original matrix A
    real(sp), intent(inout) :: bmat(:,:)
      !! On input contains the right-hand side of a system of
      !! equations Ax = b, on output contains the solution vector x
    integer, intent(in) :: ipiv(:)
      !! Pivot array
    integer, intent(out), optional :: info
      !! Exit status
    character(len=1), intent(in), optional :: trans
    character(len=1) :: tra
    integer :: n, nrhs, lda, ldb, stat
    if (present(trans)) then
      tra = trans
    else
      tra = 'n'
    end if
    lda = max(1,size(amat,1))
    ldb = max(1,size(bmat,1))
    n = size(amat,2)
    nrhs = size(bmat,2)
    call lapack_getrs(tra,n,nrhs,amat,lda,ipiv,bmat,ldb,stat)
    call handle_error(info,stat,"getrs")
  end subroutine back_sp

  subroutine back1_dp(amat,bvec,ipiv,info,trans)
    real(dp), intent(in) :: amat(:,:)
      !! LU-factored form of the original matrix A
    real(dp), intent(inout), target :: bvec(:)
      !! On input contains the right-hand side of a system of
      !! equations Ax = b, on output contains the solution vector x
    integer, intent(in) :: ipiv(:)
      !! Pivot array
    integer, intent(out), optional :: info
      !! Exit status
    character(len=1), intent(in), optional :: trans

    real(dp), pointer :: bptr(:,:)

    bptr(1:size(bvec),1:1) => bvec
    call back(amat,bptr,ipiv,info,trans)
  end subroutine back1_dp

  subroutine back_dp(amat,bmat,ipiv,info,trans)
    real(dp), intent(in) :: amat(:,:)
      !! LU-factored form of the original matrix A
    real(dp), intent(inout) :: bmat(:,:)
      !! On input contains the right-hand side of a system of
      !! equations Ax = b, on output contains the solution vector x
    integer, intent(in) :: ipiv(:)
      !! Pivot array
    integer, intent(out), optional :: info
    character(len=1), intent(in), optional :: trans
    character(len=1) :: tra
    integer :: n, nrhs, lda, ldb, stat
    if (present(trans)) then
      tra = trans
    else
      tra = 'n'
    end if
    lda = max(1,size(amat,1))
    ldb = max(1,size(bmat,1))
    n = size(amat,2)
    nrhs = size(bmat,2)
    call lapack_getrs(tra,n,nrhs,amat,lda,ipiv,bmat,ldb,stat)
    call handle_error(info,stat,"getrs")
  end subroutine back_dp

  subroutine handle_error(stat,istat,label)
    integer, intent(out), optional :: stat
    integer, intent(in) :: istat
    character(len=*), intent(in) :: label
    if (present(stat)) then
      stat = istat
    else if (istat /= 0) then
      print '(a, *(1x, g0))', "["//label//"]", "info =", istat
      stop 1
    end if
  end subroutine handle_error

end module
