from neuron import h, gui
s = h.Section() # just a home for the Glu2D instance
glu = h.Glu3d(s(.5))
h.load_file("test1.ses")


nx=int(h.nx_Glu3d)
ny=int(h.ny_Glu3d)
nz=int(h.nz_Glu3d)

sl = h.SectionList()
hin = h.PlotShape[0]
hin.size(0, nx, 0, nx)
hin.scale(0, 1)
for i in range(nx):
  for k in range(nz):
    hin.hinton(glu._ref_glu[i + k*nx*ny], float(i), float(k), .9,.9)
#h.flush_list.append(hin)
hin.exec_menu("Shape Plot")

vsrc = h.Vector([0, 0, 1, 1, 0, 0])
vsrc.mul(300)
tsrc = h.Vector([0, 1, 1, 2, 2, 3])
vsrc.play(glu._ref_src, tsrc, 1)

h.stdinit()
