NEURON {
  POINT_PROCESS Glu3d
  RANGE src
  GLOBAL nx, ny, nz
}

DEFINE NX 35
DEFINE NY 35
DEFINE NZ 3
DEFINE N 3675 : manually set to NX*NY*NZ

UNITS {
  (molar) = (1/liter)
  (mM) = (millimolar)
  (um) = (micron)
  : (nm) = (nanometer) already defined
}

PARAMETER {
  nx = NX
  ny = NY
  nz = NZ
  dx = 5 (nm) : Cleft width 15 nm.  Synapse area ~0.12 um2. ie 350nm x 350nm
  D = .4 (um2/ms) : Free solution value is 0.76um2/ms, reduced for tortuosity
  src = 500 (1) : Fusion pore diameter 10nm.  Glutamate in vesicle 2000.
                : located at 0,0,0 so x and y symmetry.
}

STATE {
  ves000
  glu[N] (1)
  clearance (1): glu has zero boundary condition outside of domain
}

INITIAL {
  ves000 = src
}

BREAKPOINT {
  SOLVE scheme METHOD sparse
}

KINETIC scheme {
  LOCAL c
  COMPARTMENT i, dx*dx*dx {glu}
  COMPARTMENT dx*dx*dx {clearance}
  COMPARTMENT dx*dx*dx {ves000}

  : Diffusion and clearance boundary condition
  c = D*dx*(1e6)
  FROM i = 0 TO NX-1 {
    FROM j = 0 TO NY-1 {
      FROM k = 0 TO NZ-1 {
        IF (i < NX-1) {
          ~ glu[(i)+(j)*NX+(k)*NX*NY] <-> glu[(i+1)+(j)*NX+(k)*NX*NY] (c, c)
        }else{
          ~ glu[(i)+(j)*NX+(k)*NX*NY] <-> clearance (c, 0)
        }
        IF (j < NY-1) {
          ~ glu[(i)+(j)*NX+(k)*NX*NY] <-> glu[(i)+(j+1)*NX+(k)*NX*NY] (c, c)
        }else{
          ~ glu[(i)+(j)*NX+(k)*NX*NY] <-> clearance (c, 0)
        }
        IF (k < NZ-1) {
          ~ glu[(i)+(j)*NX+(k)*NX*NY] <-> glu[(i)+(j)*NX+(k+1)*NX*NY] (c, c)
        }
      }
    }
  }

  : Point Source
  :~ glu[0] << ((1)*src)
  ~ ves000 <-> glu[0] (c, c)
}



