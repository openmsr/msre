import os
import openmc
import openmc.deplete
import openmc.lib
import numpy as np
from math import log10, sqrt
from collections import OrderedDict
import matplotlib.pyplot as plt
import pandas as pd
import math
import numpy as np
import re
from pathlib import Path
import h5py

def get_geom_level_from_res(depletion_path, index):

    if depletion_path is None:
        depletion_path = Path(os.getcwd())

    with h5py.File(depletion_path / 'msr_results.h5','r') as f:
        items = {}
        for key in f.keys():
            if re.split(r'_', key)[0] == 'geometry':
                items[int(re.split(r'_', key)[1])] = np.array(f.get(key))
        items = OrderedDict(sorted(items.items()))

        res = list(items.items())[index][1].mean()
        return res

salt_temp = 648.9
salt = openmc.Material(name="salt", temperature = salt_temp + 273.15)
salt.add_nuclide('Li6',1.31480070E-05)
salt.add_nuclide('Li7', 0.262960140146177)
salt.add_nuclide('Be9',1.1863E-01)
salt.add_nuclide('Zr90',1.0543E-02)
salt.add_nuclide('Zr91',2.2991E-03)
salt.add_nuclide('Zr92',3.5142E-03)
salt.add_nuclide('Zr94',3.5613E-03)
salt.add_nuclide('Zr96',5.7375E-04)
salt.add_nuclide('Hf174',8.3786E-10)
salt.add_nuclide('Hf176',2.7545E-08)
salt.add_nuclide('Hf177',9.7401E-08)
salt.add_nuclide('Hf178',1.4285E-07)
salt.add_nuclide('Hf179',7.1323E-08)
salt.add_nuclide('Hf180',1.8370E-07)
salt.add_nuclide('U234',1.034276246E-05)
salt.add_nuclide('U235',1.009695816E-03)
salt.add_nuclide('U236',4.227809892E-06)
salt.add_nuclide('U238',2.168267822E-03)
salt.add_nuclide('Fe54',2.8551E-06)
salt.add_nuclide('Fe56',4.4818E-05)
salt.add_nuclide('Fe57',1.0350E-06)
salt.add_nuclide('Fe58',1.3775E-07)
salt.add_nuclide('Cr50',2.1224E-06)
salt.add_nuclide('Cr52',4.0928E-05)
salt.add_nuclide('Cr53',4.6409E-06)
salt.add_nuclide('Cr54',1.1552E-06)
salt.add_nuclide('Ni58',5.8597E-06)
salt.add_nuclide('Ni60',2.2571E-06)
salt.add_nuclide('Ni61',9.8117E-08)
salt.add_nuclide('Ni62',3.1284E-07)
salt.add_nuclide('Ni64',7.9671E-08)
salt.add_nuclide('O16',5.1437E-04)
salt.add_nuclide('O17',1.8927E-07)
salt.add_nuclide('O18',9.6440E-07)
salt.add_nuclide('F19',5.9409E-01)
salt.set_density('g/cm3', 2.32151)
#salt.volume = 4560/2.32151*1000
#salt.volume = 70 * 0.0283168 *1e6 #Circulating primary salt: 70 ft3 (ORNL-4658)
salt.volume = 1.996 * 1e6 # https://info.ornl.gov/sites/publications/Files/Pub173113.pdf
#moderator blocks170
graphite = openmc.Material(name='graphite',temperature= salt_temp + 273.15)
graphite.set_density('g/cm3',1.86)
graphite.add_nuclide('C12',1)
graphite.add_s_alpha_beta('c_Graphite')

#inor-8
inor = openmc.Material(name='inor-8',temperature= salt_temp + 273.15)
inor.set_density('g/cm3',8.7745)
inor.add_element('Ni',(66+71)/2,'wo')
inor.add_element('Mo',(15+18)/2,'wo')
inor.add_element('Cr',(6+8)/2,'wo')
inor.add_element('Fe',5,'wo')
inor.add_element('C',(0.04+0.08)/2,'wo')
inor.add_element('Al',0.25,'wo')
inor.add_element('Ti',0.25,'wo')
inor.add_element('S',0.02,'wo')
inor.add_element('Mn',1.0,'wo')
inor.add_element('Si',1.0,'wo')
inor.add_element('Cu',0.35,'wo')
inor.add_element('B',0.010,'wo')
inor.add_element('W',0.5,'wo')
inor.add_element('P',0.015,'wo')
inor.add_element('Co',0.2,'wo')

#helium
helium = openmc.Material(name='helium')
helium.add_element('He',1.0)
helium.set_density('g/cm3',1.03*(10**-4))

#Control rods inconel clad
trace = 0.01
inconel = openmc.Material(name='inconel', temperature = 65.6 + 273.15)
inconel.add_element('Ni',78.5,percent_type='wo')
inconel.add_element('Cr',14.0,percent_type='wo')
inconel.add_element('Fe',6.5,percent_type='wo')
inconel.add_element('Mn',0.25,percent_type='wo')
inconel.add_element('Si',0.25,percent_type='wo')
inconel.add_element('Cu',0.2,percent_type='wo')
inconel.add_element('Co',0.2,percent_type='wo')
inconel.add_element('Al',0.2,percent_type='wo')
inconel.add_element('Ti',0.2,percent_type='wo')
inconel.add_element('Ta',0.5,percent_type='wo')
inconel.add_element('W',0.5,percent_type='wo')
inconel.add_element('Zn',0.2,percent_type='wo')
inconel.add_element('Zr',0.1,percent_type='wo')
inconel.add_element('C',trace,percent_type='wo')
inconel.add_element('Mo',trace,percent_type='wo')
inconel.add_element('Ag',trace,percent_type='wo')
inconel.add_element('B',trace,percent_type='wo')
inconel.add_element('Ba',trace,percent_type='wo')
inconel.add_element('Be',trace,percent_type='wo')
inconel.add_element('Ca',trace,percent_type='wo')
inconel.add_element('Cd',trace,percent_type='wo')
inconel.add_element('V',trace,percent_type='wo')
inconel.add_element('Sn',trace,percent_type='wo')
inconel.add_element('Mg',trace,percent_type='wo')
inconel.set_density('g/cm3',8.5)

# SS 316 control rod flexible hose
ss316 = openmc.Material(name='ss316', temperature = 65.6 + 273.15)
ss316.add_element('C',0.026,'wo')
ss316.add_element('Si',0.37,'wo')
ss316.add_element('Mn',0.16,'wo')
ss316.add_element('Cr',16.55,'wo')
ss316.add_element('Cu',0.16,'wo')
ss316.add_element('Ni',10,'wo')
ss316.add_element('P',0.029,'wo')
ss316.add_element('S',0.027,'wo')
ss316.add_element('Mo',2.02,'wo')
ss316.add_element('N',0.036,'wo')
ss316.add_element('Fe',70.622,'wo')
ss316.set_density('g/cm3',7.99)

#Control rods bushing posion material
Gd2O3 = openmc.Material()
Gd2O3.add_element('Gd',2)
Gd2O3.add_element('O',3)
Gd2O3.set_density('g/cm3',7.41)
Al2O3 = openmc.Material()
Al2O3.add_element('Al',2)
Al2O3.add_element('O',3)
Al2O3.set_density('g/cm3',3.95)
bush = openmc.Material.mix_materials([Gd2O3,Al2O3],[0.7,0.3],'wo')
bush.name='bush'
bush.temperature = 65.6 +273.15

#Concrete block
concrete = openmc.Material(name='concrete')
concrete.add_element('H',0.005,'wo')
concrete.add_element('O',0.496,'wo')
concrete.add_element('Si',0.314,'wo')
concrete.add_element('Ca',0.083,'wo')
concrete.add_element('Na',0.017,'wo')
concrete.add_element('Mn',0.002,'wo')
concrete.add_element('Al',0.046,'wo')
concrete.add_element('S',0.001,'wo')
concrete.add_element('K',0.019,'wo')
concrete.add_element('Fe',0.012,'wo')
concrete.set_density('g/cm3',2.35)

#Thermal shielding as water and SS305 (50-50)
water = openmc.Material()
water.add_element('H',2)
water.add_element('O',1)
water.set_density('g/cm3',0.997)

#stainless steel 304
ss304 =  openmc.Material()
ss304.add_element('C',0.08,'wo')
ss304.add_element('Mn',2,'wo')
ss304.add_element('P',0.045,'wo')
ss304.add_element('S',0.03,'wo')
ss304.add_element('Si',0.75,'wo')
ss304.add_element('Cr',19,'wo')
ss304.add_element('Ni',10,'wo')
ss304.add_element('N',0.1,'wo')
ss304.add_element('Fe',67.995, 'wo')
ss304.set_density('g/cm3',7.93)
shield = openmc.Material.mix_materials([water,ss304],[0.5,0.5],'vo')
shield.temperature = 32.2 + 273.15
shield.name='steelwater'

# "Careytemp 1600" by Philip Carey Manufacturing Compamy (Cincinnati)
# http://moltensalt.org/references/static/downloads/pdf/ORNL-TM-0728.pdf
insulation=openmc.Material(name='insulation')
insulation.add_element('Si',1)
insulation.add_element('O',2)
insulation.set_density('g/cm3',0.16) #https://www.osti.gov/servlets/purl/1411211

# sand water, not sure about this material
sandwater=openmc.Material(name='sandwater')
sandwater.add_element('Fe',3)
sandwater.add_element('O',4)
sandwater.set_density('g/cm3',6)

#Vessel anular steel
steel = openmc.Material(name='steel')
steel.add_element('Fe',1)
steel.set_density('g/cm3',7.85)

mats = openmc.Materials([salt, graphite, inor, helium, inconel, shield, concrete,
                         steel, ss316, sandwater, insulation, bush])


# CAD h5m files
core_h5m = 'h5m/msre_reactor_1e-2.h5m'
control_rod1_h5m = 'h5m/msre_control_rod_1e-2.h5m'

#Geometry
core = openmc.DAGMCUniverse(filename = core_h5m, auto_geom_ids = True,
                            universe_id = 1)
control_rod1 = openmc.DAGMCUniverse(filename = control_rod1_h5m,
                            auto_geom_ids = True, universe_id=2)
core_region = core.bounding_region()
cr1_region = control_rod1.bounding_region(boundary_type = 'transmission',
                                          starting_id=20000)

# Extend control rod region1 to include upwards translations
cr1_region = cr1_region | cr1_region.translate([0, 0, 150])

# Shift control rods 1 by offset to defined control rod 2 and 3 regions
offset = 10.163255 #cm, offset between cr1, cr2 and cr3
cr2_region = cr1_region.translate([-offset, 0, 0])
cr3_region = cr1_region.translate([-offset, offset, 0])

#Define cells
core_cell = openmc.Cell(region=~(cr1_region | cr2_region | cr3_region) & core_region ,
                        fill=core)
cr1_cell = openmc.Cell(name='CR1', region=cr1_region, fill=control_rod1)
cr2_cell = openmc.Cell(name='CR2', region=cr2_region, fill=control_rod1)
cr3_cell = openmc.Cell(name='CR3', region=cr3_region, fill=control_rod1)

#Fix control rods initial positions
start_pos = 19.2 #cm, geometrical distance between lower bottom and starting point
top_pos = 51 * 2.54 # cm, initial position of control rod1 with respect to start post
init_pos = 40 * 2.54
setattr(cr1_cell, 'translation', [0, 0, start_pos + init_pos])
setattr(cr2_cell, 'translation', [-offset, 0, start_pos + top_pos])
setattr(cr3_cell, 'translation', [-offset, offset, start_pos + top_pos])
geometry = openmc.Geometry([core_cell,cr1_cell,cr2_cell,cr3_cell])

# Settings
settings = openmc.Settings()
settings.temperature = {'method':'interpolation','range':(293.15,923.15)}
settings.batches = 60
settings.inactive = 20
settings.particles = 500000
settings.photon_transport = False
source_area = openmc.stats.Box([-100., -100., 0.],[ 100.,  100.,  200.],
              only_fissionable = True)
settings.source = openmc.Source(space = source_area)

#Tallies
tally = openmc.Tally(name="heating")
if settings.photon_transport:
    heating_score = 'heating'
else:
    heating_score = 'heating-local'
tally.scores.append(heating_score)
tallies = openmc.Tallies([tally])

# Depletion settings
model = openmc.model.Model(geometry,mats,settings,tallies)
#model.export_to_xml()

#Plots
colors = {salt:'yellow', graphite:'black', inor: 'grey', helium: 'cyan', inconel: 'grey',
              bush: 'blue', ss316: 'grey', concrete: 'brown', shield: 'red', insulation: 'green',
              sandwater: 'lightgreen', steel: 'grey'}
plots = openmc.Plots()
plot = openmc.Plot()
plot.basis = 'xy'
plot.width = (150,150)
plot.pixels = (1500,1500)
plot.origin = (0,0,160)
plot.color_by = 'material'
plot.colors = colors
model.plots.append(plot)
plot = openmc.Plot()
plot.basis = 'xz'
plot.width = (150,300)
plot.pixels = (750,1500)
plot.origin = (0,-5,150)
plot.color_by = 'material'
plot.colors = colors
model.plots.append(plot)
plot = openmc.Plot()
plot.basis = 'yz'
plot.width = (150,300)
plot.pixels = (500,1000)
plot.origin = (5,0,150)
plot.color_by = 'material'
plot.colors = colors
model.plots.append(plot)

depletion_path = Path(os.getcwd())

cell = model.geometry.get_cells_by_name('CR1')[0]
last_level = get_geom_level_from_res(depletion_path, index=-1)
setattr(cell, 'translation', [0, 0, last_level])

res_path = depletion_path / 'depletion_results.h5'
res = openmc.deplete.Results(res_path)

op = openmc.deplete.CoupledOperator(model,
                                    prev_results = res,
                                    normalization_mode = "energy-deposition")

df=pd.read_excel('MSRE_235_233_Power_History_R24E.xlsx')

#U235
timesteps = df['Duration (d)'][:94].values[len(res):]
power = df['Power (MWth)'][:94].values[len(res):]*1000000

#Add one further timestep at the end of every power run
# timesteps_ext = []
# power_ext = []
# for i in range(len(timesteps)):
#     power_ext.append(power[i])
#     if power[i] != 0.0:
#         # duplicate power value
#         power_ext.append(power[i])
#         # add timestep - 1 sec
#         timesteps_ext.append(timesteps[i] - 1/3600/24)
#         # add 1 sec
#         timesteps_ext.append(1/3600/24)
#     else:
#         timesteps_ext.append(timesteps[i])



integrator = openmc.deplete.CECMIntegrator(op, timesteps=timesteps,
                                           timestep_units='d', power=power)

integrator.add_transfer_rate('salt', ['Xe','Kr'], 4.067e-5)
integrator.add_transfer_rate('salt', ['Se','Nb','Mo','Tc','Ru','Rh','Pd','Ag','Sb','Te'], 8.777e-3)

integrator.add_batchwise('trans', axis = 2, cell_id_or_name = 'CR1',
                          bracket = [-2, 5], #cm
                          bracket_limit = [0, top_pos+start_pos],
                          atom_density_limit = 1e-14, #atoms/b-cm
                          tol = 0.1)

integrator.add_batchwise('refuel', mats_id_or_name = ['salt'],
                          mat_vector = {'U235': 1},
                          bracket = [1e2,1e3], #grams
                          bracket_limit = [0,1e5], #grams
                          tol = 0.01,
                          restart_level = start_pos + init_pos)

integrator.add_batchwise_wrap('1')

integrator.integrate()
