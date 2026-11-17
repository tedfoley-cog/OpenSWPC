contains

    !! Validation instrumentation (inserted by build_fortran_ref.sh into a scratch copy of main.f90;
    !! the upstream source tree is never modified). Controlled by SWPC_DUMP_DIR / SWPC_DUMP_EVERY.
    subroutine vdump_header(io)
        integer, intent(in) :: io
        write (io) kbeg_m, kend_m, ibeg_m, iend_m
    end subroutine vdump_header

    subroutine vdump_medium()
        character(512) :: d
        integer :: io, l
        call get_environment_variable('SWPC_DUMP_DIR', d, l)
        if (l == 0) return
        open (newunit=io, file=trim(d)//'/medium.bin', access='stream', form='unformatted', status='replace')
        call vdump_header(io)
        write (io) rho(kbeg_m:kend_m, ibeg_m:iend_m), lam(kbeg_m:kend_m, ibeg_m:iend_m), &
                   mu(kbeg_m:kend_m, ibeg_m:iend_m), taup(kbeg_m:kend_m, ibeg_m:iend_m), &
                   taus(kbeg_m:kend_m, ibeg_m:iend_m)
        close (io)
    end subroutine vdump_medium

    subroutine vdump_state(it)
        integer, intent(in) :: it
        character(512) :: d, e
        character(8) :: cit
        integer :: io, l, le, every
        call get_environment_variable('SWPC_DUMP_DIR', d, l)
        if (l == 0) return
        call get_environment_variable('SWPC_DUMP_EVERY', e, le)
        every = 1
        if (le > 0) read (e, *) every
        if (mod(it, every) /= 0) return
        write (cit, '(I8.8)') it
        open (newunit=io, file=trim(d)//'/state_'//cit//'.bin', access='stream', form='unformatted', status='replace')
        call vdump_header(io)
        write (io) Vx(kbeg_m:kend_m, ibeg_m:iend_m), Vz(kbeg_m:kend_m, ibeg_m:iend_m)
        write (io) Sxx(kbeg_m:kend_m, ibeg_m:iend_m), Szz(kbeg_m:kend_m, ibeg_m:iend_m), &
                   Sxz(kbeg_m:kend_m, ibeg_m:iend_m)
        close (io)
    end subroutine vdump_state

