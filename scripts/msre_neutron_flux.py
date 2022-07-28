import matplotlib.pyplot as plt
import openmc
from materials import *

###############################################################################
#create .png file of neutron flux (all safety rods inserted)
###############################################################################

#Geometry
h5m_filepath = 'h5m_files/msre.h5m'
graveyard = openmc.Sphere(r=10000,boundary_type='vacuum')
cad_univ = openmc.DAGMCUniverse(filename=h5m_filepath,auto_geom_ids=True)
cad_cell = openmc.Cell(region=-graveyard,fill=cad_univ)
root=openmc.Universe()
root.add_cells([cad_cell])
geometry=openmc.Geometry(root)
geometry.export_to_xml()

#materials
mats = openmc.Materials([salt,BeO,inconel,insulation,coolant,helium,stainless,boron,blanket,shield])
mats.export_to_xml()

settings = openmc.Settings()
settings.temperature = {'method':'interpolation'}
settings.batches = 100
settings.inactive = 10
settings.particles = 5000
settings.max_lost_particles = 10000
source_area = openmc.stats.Box([-200., -200., -200.],[ 200.,  200.,  200.],only_fissionable = True)
settings.source = openmc.Source(space=source_area)
settings.export_to_xml()

#tallies
tallies = openmc.Tallies()

mesh = openmc.RegularMesh()
mesh.dimension = [1000,1000]
mesh.lower_left = [-300,-300]
mesh.upper_right = [300,300]

mesh_filter = openmc.MeshFilter(mesh)

tally = openmc.Tally(name='flux')
tally.filters = [mesh_filter]
tally.scores = ['flux','fission']
tallies.append(tally)

tallies.export_to_xml()

model = openmc.model.Model(geometry, mats, settings, tallies)
sp_filename = model.run()
sp = openmc.StatePoint(sp_filename)
s_tally = sp.get_tally(scores=['flux','fission'])

flux = s_tally.get_slice(scores=['flux'])
fission = s_tally.get_slice(scores=['fission'])

flux.std_dev.shape = (1000,1000)
flux.mean.shape = (1000,1000)
fission.std_dev.shape = (1000,1000)
fission.mean.shape = (1000,1000)

fig = plt.subplot(121)
fig.axis([350,650,350,650])
fig.pixels = (2000,2000)
fig.imshow(flux.mean)
plt.savefig('neutron_flux', dpi=2000)
