module stiff3_c_mod

use, intrinsic :: iso_c_binding

implicit none
private

public :: stiff3_c
public :: stiff3_fun, stiff3_jac, stiff3_out

abstract interface
  subroutine stiff3_fun(n,y,dy,udata) bind(c)
    import c_int, c_double, c_ptr
    integer(c_int), intent(in), value :: n
    real(c_double), intent(in) :: y(n)
    real(c_double), intent(out) :: dy(n)
    type(c_ptr), intent(in), value :: udata
  end subroutine
  subroutine stiff3_jac(n,y,df,udata) bind(c)
    import c_int, c_double, c_ptr
    integer(c_int), value :: n
    real(c_double), intent(in) :: y(n)
    real(c_double), intent(out) :: dy(n,n)
    type(c_ptr), intent(in), value :: udata
  end subroutine
  subroutine stiff3_out(x,n,y,iha,qa,udata) bind(c)
    import c_int, c_double, c_ptr
    real(c_double), value :: x
    integer(c_int), value :: n
    real(c_double), intent(in) :: y(n)
    integer(c_int), value :: iha
    real(c_double), value :: qa
    type(c_ptr), value :: udata
  end subroutine
end interface

contains

  subroutine stiff3_c(n,pf,pdf,pout,nprint,x0,x1,h0,eps,w,y,udata) bind(c)

    use stiff3_solver, only: stiff3
    integer(c_int), value :: n
    type(c_funptr), value :: pf, pdf, pout
    integer(c_int), value :: nprint
    real(c_double), value :: x0, x1
    real(c_double), intent(inout) :: h0
    real(c_double), intent(in), value :: eps
    real(c_double), intent(in) :: w(n)
    real(c_double), intent(inout) :: y(n)
    type(c_ptr), value :: udata

    procedure(stiff3_rhs), pointer :: f
    procedure(stiff3_jac), pointer :: df
    procedure(stiff3_out), pointer :: out

    call c_f_funpointer (pf, f)
    call c_f_funpointer (pdf, df)
    call c_f_funpointer (pout, out)

    call stiff3 (n, ffun, fjac, fout, nprint, &
                 x0, x1, h0, eps, w, y)

  contains

    subroutine ffun(n,y,dy)
      integer, intent(in) :: n
      real(wp), intent(in) :: y(n)
      real(wp), intent(out) :: dy(n)
      call f(n,y,dy,udata)
    end subroutine

    subroutine fjac(n,y,df)
      integer, intent(in) :: n
      real(wp), intent(in) :: y(n)
      real(wp), intent(inout) :: df(n,n)
      call df(n,y,df,udata)
    end subroutine

    subroutine fout(x,y,iha,qa)
      real(wp), intent(in) :: x
      real(wp), intent(in) :: y(:)
      integer, intent(in) :: iha
      real(wp), intent(in) :: qa
      call out(x,size(y),y,iha,qa,udata)
    end subroutine

  end subroutine

end module