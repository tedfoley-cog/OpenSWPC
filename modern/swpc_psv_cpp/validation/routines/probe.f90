!! Routine-level probe of the upstream Fortran numerical utilities. Prints every result as the raw
!! IEEE bit pattern so the C++ probe (probe.cpp) can be compared bit-for-bit.
!! build_probe.sh links it against the strict reference objects; damping_profile (private in
!! m_absorb_p) is extracted verbatim from src/swpc_psv/m_absorb_p.f90 into pml_profile.inc.
module m_probe_pml
    use m_std
    implicit none
    integer :: na
    real(SP) :: fcut, dt
contains
    include 'pml_profile.inc'
end module m_probe_pml

program probe
    use m_std
    use m_fdtool
    use m_geomap
    use m_seawater
    use m_probe_pml
    implicit none
    integer :: i, j, k
    real(SP) :: t, prm(3), mxx, myy, mzz, myz, mxz, mxy, g(4), x, y, lon, lat
    real(SP) :: ts(3), zeta
    character(8) :: stf(7) = ['kupper  ', 'texp    ', 'cosine  ', 'boxcar  ', 'triangle', 'herrmann', 'asymcos ']

    do j = 1, 7
        do i = 0, 400
            t = -0.2 + i * 0.005
            prm = [0.1, 0.7, 1.1]
            write (*, '(A,1X,A,I5,1X,Z8.8)') 'stf', trim(stf(j)), i, transfer(momentrate(t, trim(stf(j)), 3, prm), 0)
        end do
    end do
    do i = 0, 36
        do j = 0, 9
            do k = 0, 12
                call sdr2moment(i * 10.0 - 180.0, j * 10.0, k * 30.0 - 180.0, mxx, myy, mzz, myz, mxz, mxy)
                write (*, '(A,3I4,6(1X,Z8.8))') 'sdr', i, j, k, transfer(mxx, 0), transfer(myy, 0), transfer(mzz, 0), &
                    transfer(myz, 0), transfer(mxz, 0), transfer(mxy, 0)
            end do
        end do
    end do
    do i = 0, 80
        write (*, '(A,I4,2(1X,Z8.8))') 'mw', i, transfer(seismic_moment(i * 0.1), 0), &
            transfer(moment_magnitude(10.0**(10 + i * 0.137)), 0)
    end do
    do i = -50, 50
        x = i * 0.0731
        write (*, '(A,I4,2I8,2(1X,Z8.8))') 'grid', i, x2i(x, -1.3, 0.025), z2k(x, -0.7, 0.0125), &
            transfer(i2x(i, -1.3, 0.025), 0), transfer(k2z(i, -0.7, 0.0125), 0)
    end do
    do i = 1, 6
        call visco_set_relaxtime(3, ts, 0.01 * i, 2.0 * i)
        zeta = visco_constq_zeta(3, 0.01 * i, 2.0 * i, ts)
        write (*, '(A,I4,5(1X,Z8.8))') 'visco', i, transfer(ts(1), 0), transfer(ts(2), 0), transfer(ts(3), 0), &
            transfer(zeta, 0), transfer(visco_chi(3, ts, zeta, 1.0), 0)
    end do
    do i = 0, 40
        lon = 139.0 + i * 0.05; lat = 35.0 + i * 0.03
        call geomap__g2c(lon, lat, 139.7604, 35.7182, 12.5, x, y)
        call geomap__c2g(x, y, 139.7604, 35.7182, 12.5, lon, lat)
        write (*, '(A,I4,4(1X,Z8.8))') 'geo', i, transfer(x, 0), transfer(y, 0), transfer(lon, 0), transfer(lat, 0)
    end do
    call seawater__init(.true.)
    do i = 0, 60
        write (*, '(A,I4,1X,Z8.8)') 'munk', i, transfer(seawater__vel(i * 0.1), 0)
    end do
    do j = 1, 3
        na = 8 * j; fcut = 0.5 * j; dt = 0.004 * j
        do i = -5, 125
            call damping_profile(i * 0.1, na * 0.1, 0.0, 12.0, g)
            write (*, '(A,2I4,4(1X,Z8.8))') 'pml', j, i, (transfer(g(k), 0), k=1, 4)
        end do
    end do
end program probe
