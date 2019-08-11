NEURON {
  POINT_PROCESS Glu3d
  RANGE src
  GLOBAL nx, ny, nz
}

DEFINE NX 25
DEFINE NY 25
DEFINE NZ 5
DEFINE N 3125 : manually set to NX*NY*NZ

UNITS {
  (molar) = (1/liter)
  (mM) = (millimolar)
  (um) = (micron)
}

PARAMETER {
  nx = NX
  ny = NY
  nz = NZ
  dx = .02 (um)
  D = .1 (um2/ms)
  src = 0 (/ms)
}

STATE {
  glu[N]
}

INITIAL {
}

BREAKPOINT {
  SOLVE scheme METHOD sparse
}

KINETIC scheme {
  LOCAL c
  COMPARTMENT i, dx*dx*dx {glu}

  : Diffusion
  c = D*dx
  FROM i = 0 TO NX-1 {
    FROM j = 0 TO NY-1 {
      FROM k = 0 TO NZ-1 {
        IF (i < NX-1) {
          ~ glu[(i)+(j)*NX+(k)*NX*NY] <-> glu[(i+1)+(j)*NX+(k)*NX*NY] (c, c)
        }
        IF (j < NY-1) {
          ~ glu[(i)+(j)*NX+(k)*NX*NY] <-> glu[(i)+(j+1)*NX+(k)*NX*NY] (c, c)
        }
        IF (k < NZ-1) {
          ~ glu[(i)+(j)*NX+(k)*NX*NY] <-> glu[(i)+(j)*NX+(k+1)*NX*NY] (c, c)
        }
      }
    }
  }

  : Point Source
  ~ glu[0] << (dx*dx*dx*src)
}



