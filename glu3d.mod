NEURON {
  POINT_PROCESS Glu3d
  RANGE vesicle_src, uncage_src
  RANGE up_rate, upglu_rate, up0, kampa, knmda, Kdampa, Kdnmda
  RANGE ijampa, ijnmda
  GLOBAL nx, ny, nz, nampa, nnmda
}

DEFINE NX 35
DEFINE NY 35
DEFINE NZ 3
DEFINE NXY 1225 : manually set to NX*NY
DEFINE NXYZ 3675 : manually set to NX*NY*NZ
DEFINE Nampa 100
DEFINE Nnmda 100

UNITS {
  (molar) = (1/liter)
  (mM) = (millimolar)
  (uM) = (micromolar)
  (nM) = (nanomolar)
  (um) = (micron)
  : (nm) = (nanometer) already defined
}

PARAMETER {
  nx = NX
  ny = NY
  nz = NZ
  dx = 5 (nm) : Cleft width 15 nm.  Synapse area ~0.12 um2. ie 350nm x 350nm
  D = .4 (um2/ms) : Free solution value is 0.76um2/ms, reduced for tortuosity
  vesicle_src = 0 (/ms) : flux from 1/4 of a vescicle
                 :Fusion pore diameter 10nm.  Glutamate in vesicle 2000.
                 : located at 0,0,0 so x and y symmetry.
  uncage_src = 0 (/ms) : flux of glu over presyn surface

  up_rate = 0 (/ms) : rate of binding of up and glu
  upglu_rate = 0 (/ms) : rate of glu unbinding into uptake
  up0 : # molecules of up per surface compartment

  nampa = Nampa
  nnmda = Nnmda
  Kdampa = 5 :(uM)
  kampa = 0 (/ms)
  Kdnmda = 25 :(nM)
  knmda = 0 (/ms)
}

ASSIGNED {
  ijampa[Nampa]
  ijnmda[Nnmda]
}

STATE {
  glu[NXYZ] (1)
  clearance (1): glu has zero boundary condition outside of domain
  upt[NXY] (1) : uptake molecule binding/unbinding generates  max uptake rate
  upglu[NXY] (1)
  uptake (1)
  ampa[Nampa]
  ampaglu[Nampa]
  nmda[Nnmda]
  nmdaglu[Nnmda]
}

INITIAL {
}

BREAKPOINT {
  SOLVE scheme METHOD sparse
}

KINETIC scheme {
  LOCAL c
  COMPARTMENT i, dx*dx*dx {glu}
  COMPARTMENT dx*dx*dx {clearance}
  COMPARTMENT dx*dx*dx {uptake}
  COMPARTMENT i, dx*dx*dx {upt upglu}
  COMPARTMENT i, dx*dx*dx {ampa ampaglu nmda nmdaglu}

  : Diffusion and clearance boundary condition
  c = D*dx*(1e6)
  FROM i = 0 TO NX-1 {
    FROM j = 0 TO NY-1 {
      FROM k = 0 TO NZ-1 {
        IF (i < NX-1) {
          ~ glu[i + j*NX + k*NXY] <-> glu[(i+1) + j*NX + k*NXY] (c, c)
        }else{
          ~ glu[i + j*NX + k*NXY] <-> clearance (c, 0)
        }
        IF (j < NY-1) {
          ~ glu[i + j*NX + k*NXY] <-> glu[i + (j+1)*NX + k*NXY] (c, c)
        }else{
          ~ glu[i + j*NX + k*NXY] <-> clearance (c, 0)
        }
        IF (k < NZ-1) {
          ~ glu[i + j*NX + k*NXY] <-> glu[i + j*NX + (k+1)*NXY] (c, c)
        }
      }
    }
  }

  : Vesicle source
  ~ glu[0] << (dx*dx*dx*vesicle_src)
  : Uncage source and uptake
  FROM i = 0 TO NX-1 {
    FROM j = 0 TO NY-1 {
      ~ glu[i + j*NX] << ((dx*dx*dx)*uncage_src)
      ~ glu[i + j*NX] + upt[i + j*NX]<-> upglu[i + j*NX] (dx*dx*dx*up_rate, 0)
      ~ upglu[i + j*NX] <-> upt[i + j*NX] + uptake (dx*dx*dx*upglu_rate, 0)
    }
  }

  : ampa
  FROM i = 0 TO Nampa-1 {
    ~ ampa[i] + glu[INT(ijampa[i])] <-> ampaglu[i] (dx*dx*dx*Kdampa*kampa, dx*dx*dx*kampa)
  }
  :nmda
  FROM i = 0 TO Nnmda-1 {
    ~ nmda[i] + glu[INT(ijnmda[i])] <-> nmdaglu[i] (dx*dx*dx*Kdnmda*knmda, dx*dx*dx*knmda)
  }
}

FUNCTION ijk(i, j, k) {
  ijk = i + j*NX + k*NXY
}

