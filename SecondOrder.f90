PROGRAM SecondOrderEuler
  !first order euler program of n body
  !declarations
  IMPLICIT NONE
  DOUBLE PRECISION :: r(1:3,1:10), v(1:3,1:10), a0(1:3,1:10), m(1:10), dr(1:3),rad2, dt, time, tend, G, au, Msun, yr, dist(1:10)
  DOUBLE PRECISION :: tout, toutnow, KE, POT, E, E0, dE, com(1:3), cov(1:3), mtot, velc(1:10)
  DOUBLE PRECISION :: runtime, timestep, a1(1:3,1:10), comdist(1:10), ecc(1:10), dsmax(1:10), dsmin(1:10), mstar, dstar, vstar
  INTEGER :: i, j, k, n !these can all be integers because all these are incremental variables
  !courtesy initial conditions
  au = 149597870700. !m
  Msun = 1.98847D30 !kg
  yr = 31557600. !s
  read(5,*) mstar, vstar, dstar

  !initial conditions in SI
  !1 is Sun
  !2 is Earth
  !3 is Jupiter
  !4 is Saturn
  !5 is Uranus
  !6 is Neptune
  !7 is incoming star
  m(1) = 1.98847D30
  m(2) = 5.9722D24
  m(3) = 1.8982D27
  m(4) = 5.683D26
  m(5) = 8.681D25
  m(6) = 1.024D26
  m(7) = mstar*Msun
  r(1,1) = 0.
  r(2,1) = 0.
  r(3,1) = 0.
  r(1,2) = 149.60D9
  r(2,2) = 0.
  r(3,2) = 0.
  r(1,3) = 778.57D9
  r(2,3) = 0.
  r(3,3) = 0.
  r(1,4) = 1433.53D9
  r(2,4) = 0.
  r(3,4) = 0.
  r(1,5) = 2872.46D9
  r(2,5) = 0.
  r(3,5) = 0.
  r(1,6) = 4495.06D9
  r(2,6) = 0.
  r(3,6) = 0.
  r(1,7) = dstar*au
  r(2,7) = 2000*au
  r(3,7) = 0.
  v(1,1) = 0.
  v(2,1) = 0.
  v(3,1) = 0.
  v(1,2) = 0.
  v(2,2) = 29.78D3
  v(3,2) = 0.
  v(1,3) = 0.
  v(2,3) = 13.06D3
  v(3,3) = 0.
  v(1,4) = 0.
  v(2,4) = 9.68D3
  v(3,4) = 0.
  v(1,5) = 0.
  v(2,5) = 6.80D3
  v(3,5) = 0.
  v(1,6) = 0.
  v(2,6) = 5.43D3
  v(3,6) = 0.
  v(1,7) = 0.
  v(2,7) = -vstar*1000.
  v(3,7) = 0.



  a0 = 0.
  a1 = 0.
  dr = 0.
  rad2 = 0.
  dt = 1000.
  time = 0.
  tend = 2*r(2,7)/(-v(2,7))
  i = 0
  j = 0
  k = 0
  n = 7
  G = 6.67408D-11
  tout = 604800.
  toutnow = 0.
  dsmax = 0.
  dsmin = 1D100

  !switch to COM and establish E0
  com(1:3) = 0.
  cov(1:3) = 0.
  mtot = 0.
  do i = 1, n
    mtot = mtot+m(i)
  end do
  do i = 1, n
    com(1:3) = com(1:3) + r(1:3,i)*m(i)
    cov(1:3) = cov(1:3) + v(1:3,i)*m(i)
    com = com/mtot
    cov = cov/mtot
  end do
  do i = 1, n
    r(1:3,i) = r(1:3,i) - com(1:3)
    v(1:3,i) = v(1:3,i) - cov(1:3)
    KE = KE + 0.5*m(i)*(v(1,i)*v(1,i) + v(2,i)*v(2,i) + v(3,i)*v(3,i))
  end do

  do i = 1, n-1
    do j = i+1, n
      dr(1:3) = r(1:3,i) - r(1:3,j)
      rad2 = dr(1)*dr(1)+dr(2)*dr(2)+dr(3)*dr(3)
      POT = POT - G*m(i)*m(j)/sqrt(rad2)
    end do
  end do
  E0 = KE+ POT


  do !start of double do loop
    a1 = 0.

    do i = 1, n
      do j = 1, n
        if (i==j) cycle
        dr(1:3) = r(1:3,i) - r(1:3,j)
        rad2 = dr(1)*dr(1) + dr(2)*dr(2) + dr(3)*dr(3)
        a1(1:3,i) = a1(1:3,i) - G*m(j)*dr(1:3)/(rad2*sqrt(rad2))
      end do
    end do


    !get a0, a1
    do i = 1, n
      r(1:3,i) = r(1:3,i) + v(1:3,i)*dt + 0.5*a0(1:3,i)*dt*dt
      v(1:3,i) = v(1:3,i) + 0.5*(a0(1:3,i) + a1(1:3,i))*dt
    end do
    a0 = a1

    if (toutnow > tout) then !occasional running subloop
      KE = 0.
      POT = 0.
      do i = 1, n
        KE = KE + 0.5*m(i)*(v(1,i)*v(1,i) + v(2,i)*v(2,i) + v(3,i)*v(3,i))
      end do
      do i = 1, n-1
        do j = i+1, n
          dr(1:3) = r(1:3,i)-r(1:3,j)
          rad2 = dr(1)*dr(1)+dr(2)*dr(2)+dr(3)*dr(3)
          POT = POT - G*m(i)*m(j)/sqrt(rad2)
        end do
      end do
      E = KE + POT
      dE = (E-E0)/E0

      write(40,*) dE
      do i = 1, n !writing out the orbits
        dist(i) = sqrt((r(1,i)-r(1,1))**2 + (r(2,i)-r(2,1))**2 +(r(3,i)-r(3,1))**2)
        velc(i) = sqrt((v(1,i)-v(1,1))**2 + (v(2,i)-v(2,1))**2 +(v(3,i)-v(3,1))**2)
        write(i+40,*) dist(i)/au
        !write(i+400,*) velc(i)
        comdist(i) = sqrt((r(1,1))**2 + (r(2,1))**2 +(r(3,1))**2)
        !write(4000+i,*) comdist(i)
        dsmax(i) = MAX(dist(i),dsmax(i))
        dsmin(i) = MIN(dist(i),dsmin(i))
        ecc(i) = (dsmax(i)-dsmin(i))/(dsmax(i)+dsmin(i))
        if (i == 100) then
        write(40000,*) r(1,2)
        write(40001,*) r(2,2)
        write(40002,*) r(1,3)
        write(40003,*) r(2,3)
        write(40004,*) r(1,4)
        write(40005,*) r(2,4)
        write(40006,*) r(1,5)
        write(40007,*) r(2,5)
        write(40008,*) r(1,6)
        write(40009,*) r(2,6)
        write(40010,*) r(3,2)
        write(40020,*) r(1,1)
        write(40021,*) r(2,1)
        write(41000,*) (sqrt(a1(1,2)**2+a1(2,2)**2+a1(3,2)**2))
      endif

      end do
      toutnow = 0.
    endif
    write(1,*) ecc(6)
    if (ecc(6)>0.1 .AND. time>tout) write(6,*) 'Ended Due to High Neptune Eccentricity'
    if (ecc(6)>0.1 .AND. time>tout) exit
    if (ecc(5)>0.1 .AND. time>tout) write(6,*) 'Ended Due to High Uranus Eccentricity'
    if (ecc(5)>0.1 .AND. time>tout) exit
    toutnow = toutnow + dt
    time = time + dt
    if (time > tend) exit
  end do
  !write(1,*) cov
  !do i = 1,n
    !dist(i) = sqrt((r(1,i)-r(1,1))**2 + (r(2,i)-r(2,1))**2 +(r(3,i)-r(3,1))**2)
    !velc(i) = sqrt((v(1,i)-v(1,1))**2 + (v(2,i)-v(2,1))**2 +(v(3,i)-v(3,1))**2)
    !write(i+400,*) dist(i), velc(i)
  !end do
END PROGRAM SecondOrderEuler
