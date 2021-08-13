module stiff3_linalg_dss_helpers

  use mkl_dss, only: &
    MKL_DSS_SUCCESS, MKL_DSS_ZERO_PIVOT, MKL_DSS_OUT_OF_MEMORY, &
    MKL_DSS_FAILURE, MKL_DSS_ROW_ERR, MKL_DSS_COL_ERR, &
    MKL_DSS_TOO_FEW_VALUES, MKL_DSS_TOO_MANY_VALUES, MKL_DSS_NOT_SQUARE, &
    MKL_DSS_STATE_ERR, MKL_DSS_INVALID_OPTION, MKL_DSS_OPTION_CONFLICT, &
    MKL_DSS_MSG_LVL_ERR, MKL_DSS_TERM_LVL_ERR, MKL_DSS_STRUCTURE_ERR, &
    MKL_DSS_REORDER_ERR, MKL_DSS_VALUES_ERR, &
    MKL_DSS_STATISTICS_INVALID_MATRIX, MKL_DSS_STATISTICS_INVALID_STATE, &
    MKL_DSS_STATISTICS_INVALID_STRING,

  implicit none
  private

  public :: process_retval
  public :: retval_to_string

contains

  subroutine process_retval(func,ret,ierr,errmsg)
    use, intrinsic :: iso_fortran_env, only: error_unit
    character(len=*), intent(in) :: func
    integer, intent(in) :: ret
    integer, intent(out), optional :: ierr
    character(len=*), intent(out), optional :: errmsg

    if (present(ierr)) then
      ierr = ret
    else
      if (ret /= MKL_DSS_SUCCESS) then
        write(error_unit,'(A,*(1x,g0))') "["//func//"]", retval_to_string(ret)
        error stop
      end if
    end if

    if (present(errmsg)) then
      errmsg = retval_to_string(ret)
    end if

  end subroutine

  function retval_to_string(ret) result(str)
    integer, intent(in) :: ret
    character(len=:), allocatable :: str

    select case(ret)
    case(MKL_DSS_SUCCESS); str = "SUCCESS"
    case(MKL_DSS_ZERO_PIVOT); str = "ZERO_PIVOT"
    case(MKL_DSS_OUT_OF_MEMORY); str = "OUT_OF_MEMORY"
    case(MKL_DSS_FAILURE); str = "FAILURE"
    case(MKL_DSS_ROW_ERR); str = "ROW_ERR"
    case(MKL_DSS_COL_ERR); str = "COL_ERR"
    case(MKL_DSS_TOO_FEW_VALUES); str = "TOO_FEW_VALUES"
    case(MKL_DSS_TOO_MANY_VALUES); str = "TOO_MANY_VALUES"
    case(MKL_DSS_NOT_SQUARE); str = "NOT_SQUARE"
    case(MKL_DSS_STATE_ERR); str = "STATE_ERR"
    case(MKL_DSS_INVALID_OPTION); str = "INVALID_OPTION"
    case(MKL_DSS_OPTION_CONFLICT); str = "OPTION_CONFLICT"
    case(MKL_DSS_MSG_LVL_ERR); str = "MSG_LVL_ERR"
    case(MKL_DSS_TERM_LVL_ERR); str = "TERM_LVL_ERR"
    case(MKL_DSS_STRUCTURE_ERR); str = "STRUCTURE_ERR"
    case(MKL_DSS_REORDER_ERR); str = "REORDER_ERR"
    case(MKL_DSS_VALUES_ERR); str = "VALUES_ERR"
    case(MKL_DSS_STATISTICS_INVALID_MATRIX)
      str = "STATISTICS_INVALID_MATRIX"
    case(MKL_DSS_STATISTICS_INVALID_STATE)
      str = "STATISTICTS_INVALID_STATE"
    case(MKL_DSS_STATISTICS_INVALID_STRING)
      str = "STATISTICS_INVALID_STRING"
    end select

  end subroutine

end module

module stiff3_linalg_dss

  use mkl_dss, only: &
    mkl_dss_handle, &
    dss_create, &
    dss_define_structure, &
    dss_reorder, &
    dss_factor_real, &
    dss_solve_real, &
    dss_delete, &
    dss_statistics, &
    mkl_cvt_to_null_terminated_str

  use stiff3_linalg_dss_helpers, only: &
    process_retval, retval_to_string

  implicit none
  private

  public :: dss_solver, create, destroy

  type :: dss_solver
    private
    integer :: nrows, ncols, nnz
    type(mkl_dss_handle) :: handle
    integer, allocatable :: perm(:)

  contains
    final :: destroy
    procedure :: factorize
    procedure :: solve_v, solve_m
    generic :: solve => solve_v, solve_m
    procedure :: get_perm
  end type

  interface create
    module procedure create_dss_solver
  end interface

  interface destroy
    module procedure :: destroy_dss_solver
  end interface

contains

  subroutine create_dss_solver(self,ierr,errmsg)
    type(dss_solver), intent(out) :: self
    logical, intent(in), optional :: saveperm
    integer, intent(out), optional :: ierr
    character(len=*), intent(out), optional :: errmsg

    integer :: opt
    integer :: ret

    !
    ! create the direct sparse solver handle
    !
    !   settings available include:
    !     - MSG_LVL_WARNING and TERM_LVL_WARNING
    !     - SINGLE_PRECISION
    !     - ZERO_BASED_INDEXING
    !     - REFINEMENT_OFF
    !     - OUT OF CORE

    ret = dss_create(self%handle,opt)
    call process_retval('dss_create',ret,ierr,errmsg)

    !
    ! communicate the location of non-zero elements (sparsity pattern)
    !
    ret = dss_define_structure(self%handle,&
      opt,rowIndex,nRows,nCols,columns,nNonzeros)
    call process_retval('dss_define_structure',ret,ierr,errmsg)

    !
    ! compute permutation that minimizes fill-in
    !
    ret = dss_reorder(self%handle,opt,self%perm)
    call process_retval('dss_reorder',ret,ierr,errmsg)

  end subroutine


  subroutine get_permutation_vector(self,perm)
    type(dss_solver), intent(in) :: self
    integer, intent(out) :: perm(solver%nrows)
    integer, intent(out), optional :: ierr
    character(len=*), intent(out), optional :: errmsg
    integer, parameter :: opt = MKL_DSS_GET_ORDER

    ret = dss_reorder(self%handle,opt,perm)
    call process_retval('dss_reorder',ret,ierr,errmsg)
  end subroutine

  subroutine factorize(self,rvalues,ierr,errmsg)
    type(dss_solver), intent(inout) :: self
    real(dp), intent(in) :: rvalues(:)
      ! Size should be equal to self%nnz
    integer, intent(out), optional :: ierr
    character(len=*), intent(out), optional :: errmsg

    integer :: opt, ret

    opt = MKL_DSS_INDEFINITE

    ret = dss_factor_real(self%handle,opt,rvalues)
    call process_retval("dss_factor_real",ret,ierr,errmsg)

  end subroutine factor

  subroutine solve_v(self,rhs,sol,ierr,errmsg)
    class(dss_solver), intent(inout) :: self
    real(dp), intent(in) :: rhs(:)
    real(dp), intent(out) :: sol(:)
    integer, intent(out), optional :: ierr
    character(len=*), intent(out), optional :: errmsg

    integer :: opt, nrhs, ret

    opt = MKL_DSS_DEFAULTS ! full solution step
    nrhs = 1

    ret = dss_solve_real(self%handle,opt,&
      rrhsvalues=rhs, &
      nrhs=1, &
      rsolvalues=sol)

    call process_retval("dss_solve_real",ret,ierr,errmsg)

  end subroutine

  subroutine solve_m(self,rhs,sol,ierr,errmsg)
    class(dss_solver), intent(inout) :: self
    real(dp), intent(in), target :: rhs(:,:)
    real(dp), intent(out), target :: sol(:,:)
    integer, intent(out), optional :: ierr
    character(len=*), intent(out), optional :: errmsg

    integer :: nrhs, nptr
    real(dp) :: rhsptr(:), solptr(:)

    nrhs = size(rhs,2)
    nptr = solver%nrows * nrhs
    rhsptr(1:nptr) => rhs
    solptr(1:nptr) => sol

    opt = MKL_DSS_DEFAULTS

    ret = dss_solve_real(self%handle,opt,&
      rrhsvalues=rhsptr, &
      nrhs=1, &
      rsolvalues=solptr)

    call process_retval("dss_solve_real",ret,ierr,errmsg)

  end subroutine

  subroutine stats(solver)

    call dss_statistics(handle,opt,statArr,retValues)
  end subroutine

  subroutine destroy_dss_solver(self,ierr,errmsg)
    type(dss_solver), intent(inout) :: self
    integer, intent(out), optional :: ierr
    character(len=*), intent(out), optional :: errmsg

    integer :: opt
    opt = MKL_DSS_DEFAULTS

    ret = dss_delete(self%handle,opt)
    call process_retval("dss_delete",ret,ierr,errmsg)

  end subroutine

end module

module stiff3_solver_s

  implicit none

  abstract interface
    !> Function to evaluate the right-hand side of a system of ODEs.
    !> It is assumed the system of ODEs is autonomous, meaning that
    !> the independent variable x, does not appear explicitly.
    subroutine rhs_sub(n,y,f)
      import wp
      integer, intent(in) :: n
      real(wp), intent(in) :: y(n)
      real(wp), intent(inout) :: f(n)
    end subroutine

    !> User supplied subprogram for evaluation of the Jacobian.
    subroutine jacobian_sub(n,y,ia,ja,df)
      import wp
      integer, intent(in) :: n
      real(wp), intent(in) :: y(n)
      integer, intent(in) :: ia(n+1)
      integer, intent(in) :: ja(:)
      real(wp), intent(inout) :: df(:)
    end subroutine

  end interface

end module