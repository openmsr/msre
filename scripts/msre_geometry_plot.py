import matplotlib
import openmc
from materials import *

###############################################################################
#generate geometry plot of are (all safety rods fully inserted)
###############################################################################

h5m_filepath = './h5m_files/msre.h5m'

# materials
mats = openmc.Materials([salt,BeO,inconel,insulation,coolant,helium,stainless,boron])
mats.export_to_xml()

dag_univ = openmc.DAGMCUniverse(h5m_filepath)
geom = openmc.Geometry(root=dag_univ)
geom.export_to_xml()

#plotting geometry
plots = openmc.Plots()

x_width = 300
y_width = 300

#xy plot
p1 = openmc.Plot()
p1.width = (x_width,y_width)
p1.pixels = (1000, 1000)
p1.color_by = 'material'

#xz plot (split plane)
p2 = openmc.Plot()
p2.basis = 'xz'
#p2.origin = (0,5,0)
p2.width = (x_width,y_width)
p2.pixels = (1000, 1000)
p2.color_by = 'material'

p3 = openmc.Plot()
p3.basis = 'yz'
p3.width = (x_width,y_width)
p3.pixels = (1000, 1000)
p3.color_by = 'material'

plots.append(p1)
plots.append(p2)
plots.append(p3)
plots.export_to_xml()

openmc.plot_geometry()
