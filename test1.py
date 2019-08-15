from neuron import h, gui
s = h.Section() # just a home for the Glu2D instance
glu = h.Glu3d(s(.5))
h.load_file("test1.ses")


nx=int(h.nx_Glu3d)
ny=int(h.ny_Glu3d)
nz=int(h.nz_Glu3d)
h('nx=%d'%nx)
h('ny=%d'%ny)
h('nz=%d'%nz)
h('''func ijk() { return $1 + nx*$2 + nx*ny*$3 }''')
dx = h.dx_Glu3d # (nanometer)

def molec2mM(x):
  vol = dx*dx*dx * (1e-24) # (liters)
  return x/vol/6.0221418e23 * 1000


sl = h.SectionList()
hin = h.PlotShape[0]
hin.size(0, nx, 0, nx)
hin.scale(0, 1)
for i in range(nx):
  for k in range(nz):
    hin.hinton(glu._ref_glu[i + k*nx*ny], float(i), float(k), .9,.9)
#h.flush_list.append(hin)
hin.exec_menu("Shape Plot")

# plot of glu[i, 0, nz-1]
gluline = h.PtrVector(nx)
gluline.label("glu[i, 0, nz-1]")
for i in range(nx):
  gluline.pset(i, glu._ref_glu[int(h.ijk(i, 0, nz-1))])
glulinegraph = h.Graph[2]
gluline.plot(glulinegraph)

gluline2 = h.PtrVector(nx)
gluline2.label("glu[i, 0, 0]")
for i in range(nx):
  gluline2.pset(i, glu._ref_glu[int(h.ijk(i, 0, 0))])
gluline2.plot(glulinegraph)
#glulinegraph2.size(0, 30, 0, 3)
h.flush_list.append(glulinegraph)

# 500 Glutamate molecules in 1/4 vesicle at 0,0,0 dump in 100us
delay = 0.01 #ms
dur = 0.01 #ms
amp = 500/dur
vsrc = h.Vector([0, 0, amp, amp, 0, 0])
tsrc = h.Vector([0, delay, delay, delay+dur, delay+dur, 1e9])
vsrc.play(glu._ref_vesicle_src, tsrc, 1)

h.stdinit()
