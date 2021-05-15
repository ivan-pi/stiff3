module linalg

  implicit none
  private

  public :: lu, back

  !   !> LU - factorization of a matrix
  !   subroutine lu(n,ipiv,a)
  !     import dp
  !     integer, intent(in) :: n
  !       !! Size of the matrix
  !     integer, intent(out) :: ipiv(n)
  !       !! Pivot array
  !     real(dp), intent(inout) :: a(n,n)
  !       !! On input a matrix, on output the LU factored matrices
  !   end subroutine
  interface lu
    module procedure lu_sp, lu_dp
  end interface


  !   !> Performs back-substition
  !   !> Used together with procedure lu, to calculate the solution
  !   !> of a system of linear equations
  !   subroutine back(n,ipiv,a,b)
  !     import dp
  !     integer, intent(in) :: n
  !       !! Size of the linear system of equations
  !     integer, intent(in) :: ipiv(n)
  !       !! Pivot array
  !     real(dp), intent(in) :: a(n,n)
  !       !! LU-factored form of the original matrix A
  !     real(dp), intent(inout) :: b(n)
  !       !! On input contains the right-hand side of a system of
  !       !! equations Ax = b, on output contains the solution vector x
  !   end subroutine
  interface back
    module procedure back_sp, back_dp
  end interface

  integer, parameter :: sp = kind(1.0e0)
  integer, parameter :: dp = kind(1.0d0)

contains

  subroutine lu_sp(n,ipiv,a)
    integer, intent(in) :: n
    integer, intent(out) :: ipiv(n)
    real(sp), intent(inout) :: a(n,n)

    integer :: lda, info
    external :: dgetrf

    lda = n

    call sgetrf(n,n,a,lda,ipiv,info)

    if (info /= 0) then
      print *, "[dgetrf] info = ", info
    end if

  end subroutine
  subroutine lu_dp(n,ipiv,a)
    integer, intent(in) :: n
    integer, intent(out) :: ipiv(n)
    real(dp), intent(inout) :: a(n,n)

    integer :: lda, info
    external :: dgetrf

    lda = n

    call dgetrf(n,n,a,lda,ipiv,info)

    if (info /= 0) then
      print *, "[dgetrf] info = ", info
    end if

  end subroutine


  subroutine back_sp(n,ipiv,a,b)
    integer, intent(in) :: n
    integer, intent(in) :: ipiv(n)
    real(sp), intent(in) :: a(n,n)
    real(sp), intent(inout) :: b(n)

    integer :: lda, ldb, info
    external :: dgetrs

    lda = n
    ldb = n

    call sgetrs('N',n,1,a,lda,ipiv,b,ldb,info)

    if (info /= 0) then
      print *, "[dgetrs] info = ", info
    end if

  end subroutine
  subroutine back_dp(n,ipiv,a,b)
    integer, intent(in) :: n
    integer, intent(in) :: ipiv(n)
    real(dp), intent(in) :: a(n,n)
    real(dp), intent(inout) :: b(n)

    integer :: lda, ldb, info
    external :: dgetrs

    lda = n
    ldb = n

    call dgetrs('N',n,1,a,lda,ipiv,b,ldb,info)

    if (info /= 0) then
      print *, "[dgetrs] info = ", info
    end if

  end subroutine

end module