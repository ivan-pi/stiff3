module stiff3_linalg_b

  implicit none
  private

  public :: lub, backb
  public :: lapack_gbtrf, lapack_gbtrs, lapack_gbsv

  interface lub
    module procedure lub_sp, lub_dp
  end interface

  interface backb
    module procedure back_sp, back1_sp, back_dp, back1_dp
  end interface

  integer, parameter :: sp = kind(1.0e0)
  integer, parameter :: dp = kind(1.0d0)

  interface lapack_gbtrf
    subroutine sgbtrf(m,n,kl,ku,ab,ldab,ipiv,info)
      import :: sp
      integer, intent(in) :: m, n, kl, ku, ldab
      real(sp), intent(inout) :: ab(ldab,*)
      integer, intent(out) :: ipiv(*), info
    end subroutine
    subroutine dgbtrf(m,n,kl,ku,ab,ldab,ipiv,info)
      import :: dp
      integer, intent(in) :: m, n, kl, ku, ldab
      real(dp), intent(inout) :: ab(ldab,*)
      integer, intent(out) :: ipiv(*), info
    end subroutine
  end interface lapack_gbtrf

  interface lapack_gbtrs
    subroutine sgbtrs(trans,n,kl,ku,nrhs,ab,ldab,ipiv,b,ldb,info)
      import sp
      character(len=1), intent(in) :: trans
      integer, intent(in) :: n, kl, ku, nrhs, ldab, ldb
      real(sp), intent(in) :: ab(ldab,*)
      integer, intent(in) :: ipiv(*)
      real(sp), intent(inout) :: b(ldb,*)
      integer, intent(out) :: info
    end subroutine sgbtrs
    subroutine dgbtrs(trans,n,kl,ku,nrhs,ab,ldab,ipiv,b,ldb,info)
      import dp
      character(len=1), intent(in) :: trans
      integer, intent(in) :: n, kl, ku, nrhs, ldab, ldb
      real(dp), intent(in) :: ab(ldab,*)
      integer, intent(in) :: ipiv(*)
      real(dp), intent(inout) :: b(ldb,*)
      integer, intent(out) :: info
    end subroutine dgbtrs
  end interface lapack_gbtrs

  interface lapack_gbsv
    subroutine sgbsv(n,kl,ku,nrhs,ab,ldab,ipiv,b,ldb,info)
      import sp
      integer, intent(in) :: n, kl, ku, nrhs, ldab, ldb
      real(sp), intent(inout) :: ab(ldab,*), b(ldb,*)
      integer, intent(out) :: ipiv(*), info
    end subroutine sgbsv
    subroutine dgbsv(n,kl,ku,nrhs,ab,ldab,ipiv,b,ldb,info)
      import dp
      integer, intent(in) :: n, kl, ku, nrhs, ldab, ldb
      real(dp), intent(inout) :: ab(ldab,*), b(ldb,*)
      integer, intent(out) :: ipiv(*), info
    end subroutine dgbsv
  end interface lapack_gbsv

contains

  subroutine lub_sp(kl,ku,abmat,ipiv,info)
    integer, intent(in) :: kl, ku
    real(sp), intent(in) :: abmat(:,:)
    integer, intent(out) :: ipiv(:)
    integer, intent(out) :: info

    integer :: m, n, lda, stat

    ldab = max(1,size(abmat,1))
    m = size(abmat,1)
    n = size(abmat,2)
    call lapack_gbtrf(kl,ku,m,n,abmat,ldab,ipiv,stat)

    call handle_error(info,stat,"gbtrf")

  end subroutine

  subroutine lub_dp()

  end subroutine

  subroutine back1_sp()

  end subroutine
  subroutine back_sp()

  end subroutine

  subroutine back1_dp()

  end subroutine
  subroutine back_dp()

  end subroutine
end module