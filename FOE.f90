PROGRAM FirstOrderEuler
  !first order euler program of n body
  !declarations
  IMPLICIT NONE
  DOUBLE PRECISION :: r(1:3,1:10), v(1:3,1:10), a(1:3,1:10), m(1:6), dr(1:3),rad2, dt, time, tend, G, au, Msun, yr, dist(1:10)
  DOUBLE PRECISION :: tout, toutnow, KE, POT, E, E0, dE, com(1:3), cov(1:3), mtot, velc(1:10)
  DOUBLE PRECISION :: runtime, timestep
  INTEGER :: i, j, k, n !these can all be integers because all these are incremental variables
  !courtesy initial conditions
  au = 149597870700. !m
  Msun = 1.98847D30 !kg
  yr = 31557600. !s

  !initial conditions in SI
  !1 is Sun
  !2 is Earth
  !3 is Jupiter
  !4 is Saturn
  !5 is Uranus
  !6 is Neptune
  m(1) = 1.98847D30
  m(2) = 5.9722D24
  m(3) = 1.8982D27
  m(4) = 5.683D26
  m(5) = 8.681D25
  m(6) = 1.024D26
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
  v(1,1) = 0.
  v(2,1) = 0.
  v(3,1) = 0.
  v(1,2) = 0.
  v(2,2) = 29.78D3
  v(3,2) = 0.
  v(1,3) = 0.
  v(2,3) = 13.07D3
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

  read(5,*) runtime, timestep

  a = 0.
  dr = 0.
  rad2 = 0.
  dt = timestep
  time = 0.
  tend = runtime * yr
  i = 0
  j = 0
  k = 0
  n = 6
  G = 6.67408D-11
  tout = 604800.
  toutnow = 604801.


  !switch to COM and establish E0
  com(1:3) = 0.
  cov(1:3) = 0.
  do i = 1, n
    com(1:3) = com(1:3)+r(1:3,i)*m(i)
    cov(1:3) = cov(1:3)+v(1:3,i)*m(i)
    mtot = mtot+m(i)

    com = com/mtot
    cov = cov/mtot
    r(1:3,i) = r(1:3,i)-com(1:3)
    v(1:3,i) = v(1:3,i)-cov(1:3)
    KE = KE + 0.5*m(i)*(v(1,i)*v(1,i) + v(2,i)*v(2,i) + v(3,i)*v(3,i))
  end do

  do i = 1, n-1
    do j = i+1, n
      dr(1:3) = r(1:3,i)-r(1:3,j)
      rad2 = dr(1)*dr(1)+dr(2)*dr(2)+dr(3)*dr(3)
      POT = POT - (G*m(i)*m(j))/(sqrt(rad2))
    end do
  end do
  E0 = KE+ POT


  do !start of double do loop
    a = 0.
    toutnow = toutnow + dt
    time = time + dt
    do i = 1, n
      do j = 1, n
        if (i==j) cycle
        dr(1:3) = r(1:3,i) - r(1:3,j)
        rad2 = dr(1)*dr(1) + dr(2)*dr(2) + dr(3)*dr(3)
        a(1:3,i) = a(1:3,i) - G*m(j)*dr(1:3)/(rad2*sqrt(rad2))
      end do
    end do


    !get ALL a
    do i=1, n
      r(1:3,i) = r(1:3,i) + v(1:3,i)*dt
      v(1:3,i) = v(1:3,i) + a(1:3,i)*dt
    end do


    if (toutnow > tout) then !occasional running program
      KE = 0.
      POT = 0.
      do i = 1, n
        KE = KE + 0.5*m(i)*(v(1,i)*v(1,i) + v(2,i)*v(2,i) + v(3,i)*v(3,i))
      end do
      do i = 1, n-1
        do j = i+1, n
          dr(1:3) = r(1:3,i)-r(1:3,j)
          rad2 = dr(1)*dr(1)+dr(2)*dr(2)+dr(3)*dr(3)
          POT = POT - (G*m(i)*m(j))/(sqrt(rad2))
        end do
      end do
      E = KE + POT
      dE = (E-E0)/E0
      write(30,*) dE
      do i = 1, n !writing out the orbits
        dist(i) = sqrt((r(1,i)-r(1,1))**2 + (r(2,i)-r(2,1))**2 +(r(3,i)-r(3,1))**2)
        velc(i) = sqrt((v(1,i)-v(1,1))**2 + (v(2,i)-v(2,1))**2 +(v(3,i)-v(3,1))**2)
        write(i+30,*) dist(i)
        write(i+300,*) velc(i)
      end do
    toutnow = 0.
    endif


    if (time > tend) exit
  end do
  do i = 1,n
    dist(i) = sqrt((r(1,i)-r(1,1))**2 + (r(2,i)-r(2,1))**2 +(r(3,i)-r(3,1))**2)
    velc(i) = sqrt((v(1,i)-v(1,1))**2 + (v(2,i)-v(2,1))**2 +(v(3,i)-v(3,1))**2)
    write(i+100,*) dist(i), velc(i)
  end do
END PROGRAM FirstOrderEuler
