import os
import openmc
import openmc.deplete
import openmc.lib
import numpy as np
from math import log10
from collections import OrderedDict
import matplotlib.pyplot as plt
import pandas as pd
from scipy import stats
import math

def define_materials():
    # fuel salt
    salt = openmc.Material(name="salt", temperature = 911)
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
    salt.set_density('g/cm3', 2.3275)

    #moderator blocks170
    graphite = openmc.Material(name='graphite',temperature=911)
    graphite.set_density('g/cm3',1.86)
    graphite.add_nuclide('C12',1)
    graphite.add_s_alpha_beta('c_Graphite')

    #inor-8
    inor = openmc.Material(name='inor-8',temperature=911)
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
    inconel = openmc.Material(name='inconel', temperature = 911)
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
    ss316 = openmc.Material(name='ss316')
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

    # "Careytemp 1600" by Philip Carey Manufacturing Compamy (Cincinnati) from http://moltensalt.org/references/static/downloads/pdf/ORNL-TM-0728.pdf
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

    mats = openmc.Materials([salt,graphite,inor,helium,inconel,shield,concrete,steel,ss316,sandwater,insulation,bush])
    return mats

def build(make_tally=True, plot_geom=True):
    #Clean-up
    os.system("rm *.xml *.h5 *.out")

    # CAD h5m files
    core_h5m = 'h5m/MSRE_Core_1e-2/MSRE_Core_1e-2.h5m'
    control_rod1_h5m = 'h5m/Msre_control_rod_1e-2/Msre_control_rod_1e-2.h5m'

    #Materials
    mats=define_materials()

    #Geometry
    core = openmc.DAGMCUniverse(filename=core_h5m, auto_geom_ids=True, universe_id=1)
    control_rod1 = openmc.DAGMCUniverse(filename=control_rod1_h5m, auto_geom_ids=True, universe_id=2)
    core_region = core.bounding_region()
    cr1_region = control_rod1.bounding_region(boundary_type='transmission', starting_id=20000)
    # Extend control rod region to include upwards translations
    cr1_region = cr1_region | cr1_region.translate([0,0,150])
    cr2_region = cr1_region.translate([-10.163255,0,0])
    cr3_region = cr1_region.translate([-10.163255,10.163255,0])
    core_cell = openmc.Cell(region=~(cr1_region | cr2_region | cr3_region) & core_region , fill=core)
    cr1_cell = openmc.Cell(name='CR1', region=cr1_region, fill=control_rod1)
    cr2_cell = openmc.Cell(name='CR2', region=cr2_region, fill=control_rod1)
    cr3_cell = openmc.Cell(name='CR3', region=cr3_region, fill=control_rod1)
    #translate control rod1 130 cm (51 inch) up in z direction
    setattr(cr1_cell, 'translation', [0,0,130+19.2])
    setattr(cr2_cell, 'translation', [-10.163255,0,130+19.2])
    setattr(cr3_cell, 'translation', [-10.163255,10.163255,130+19.2])
    geometry = openmc.Geometry([core_cell,cr1_cell,cr2_cell,cr3_cell])

    # Settings
    settings = openmc.Settings()
    settings.temperature = {'method':'interpolation','range':(293.15,923.15)}
    settings.batches = 50
    settings.inactive = 10
    settings.particles = 50000
    settings.photon_transport = False
    source_area = openmc.stats.Box([-100., -100., 0.],[ 100.,  100.,  200.],only_fissionable = True)
    settings.source = openmc.Source(space=source_area)

    if make_tally:
        tally = OrderedDict()
        tally['general'] = openmc.Tally(name="General")
        if settings.photon_transport:
            heating_score = 'heating'
        else:
            heating_score = 'heating-local'
        tally['general'].scores.append(heating_score)

        e_min, e_max = 1e-5, 20e6
        groups = 500
        energies = np.logspace(log10(e_min), log10(e_max), groups + 1)
        energy_filter = openmc.EnergyFilter(energies)
        particle_filter = openmc.ParticleFilter(['neutron'])
        cell_filter = openmc.MaterialFilter([mats[0]])
        # mesh = openmc.RegularMesh()
        # mesh.dimension = [1, 1, 1]
        # mesh.lower_left = [39.7, -3.1, 39.2]
        # mesh.upper_right = [41.7, -2, 195]
        # mesh_filter = openmc.MeshFilter(mesh)
        tally['flux']= openmc.Tally(name="flux")
        tally['flux'].filters = [energy_filter, particle_filter] #, cell_filter]#, mesh_filter]
        tally['flux'].scores = ['flux']

        mesh = openmc.RegularMesh()
        mesh.dimension = [500, 500, 1]
        mesh.lower_left = [-100, -100, 50]
        mesh.upper_right = [100, 100, 200]
        mesh_filter = openmc.MeshFilter(mesh)
        tally['mesh'] = openmc.Tally(name="Mesh")
        tally['mesh'].scores = ['flux','absorption','fission','scatter']
        tally['mesh'].filters = [mesh_filter]
        tally['mesh'].filters.append(particle_filter)

        tally['leak'] = openmc.Tally(name='leakage')
        mesh = openmc.RegularMesh()
        mesh.dimension = [1, 1, 1]
        mesh.lower_left = core.bounding_box[0]-10
        mesh.width = (core.bounding_box[1]+10)*2
        meshsurface_filter = openmc.MeshSurfaceFilter(mesh)
        tally['leak'].filters = [meshsurface_filter]
        tally['leak'].scores = ['current']
        tallies = openmc.Tallies(tally.values())
        model = openmc.model.Model(geometry,mats,settings,tallies)
    else:
        model = openmc.model.Model(geometry,mats,settings)

    if plot_geom:
        colors = {[i for i in mats if i.name=="salt"][0]: 'yellow',
                [i for i in mats if i.name=="graphite"][0]: 'black',
                [i for i in mats if i.name=="inor-8"][0]: 'grey',
                [i for i in mats if i.name=="helium"][0]: 'cyan',
                [i for i in mats if i.name=="inconel"][0]: 'grey',
                [i for i in mats if i.name=="bush"][0]: 'blue',
                [i for i in mats if i.name=="ss316"][0]: 'grey',
                [i for i in mats if i.name=="concrete"][0]: 'brown',
                [i for i in mats if i.name=="steelwater"][0]: 'red',
                [i for i in mats if i.name=="insulation"][0]: 'green',
                [i for i in mats if i.name=="sandwater"][0]: 'lightgreen',
                [i for i in mats if i.name=="steel"][0]: 'grey'}

        plot_file = openmc.Plots()
        plot1 = openmc.Plot()
        plot1.width = [150, 150]
        plot1.pixels = [2000, 2000]
        plot1.origin = [0,0,150]
        plot1.basis = 'xy'
        plot1.color_by = "material"
        plot1.colors = colors

        plot2 = openmc.Plot()
        plot2.width = [200, 500]
        plot2.pixels = [4000,10000]
        plot2.origin = [-5,0,150]
        plot2.basis = 'yz'
        plot2.color_by = "material"
        plot2.colors = colors

        plot3 = openmc.Plot()
        plot3.width = [200, 500]
        plot3.pixels = [4000,10000]
        plot3.origin = [0,-5,150]
        plot3.basis = 'xz'
        plot3.color_by = "material"
        plot3.colors = colors

        model.plots.append(plot1)
        model.plots.append(plot2)
        model.plots.append(plot3)
        model.plot_geometry()

    return model

def run(model, mass, power):
    results=model.run()
    salt_vol = 1.65058e6 #fuel salt core volume [cm3] from OnShape model
    #core_vol = 342.76 # core vol fuel cell [cc]
    sp = openmc.StatePoint(results)
    heating = sp.get_tally(name="General").get_pandas_dataframe()["mean"].sum()*openmc.data.JOULE_PER_EV
    fac = power/heating
    t = sp.get_tally(name="flux")
    energy_filter = t.filters[0]
    energies = energy_filter.bins[:, 0]
    mean = t.mean.ravel()
    uncertainty = t.get_values(value='std_dev').ravel()

    fig, ax = plt.subplots()
    ax.plot(energies, mean*fac/salt_vol, drawstyle='steps-post')
    ax.set_xlabel('Energy [eV]')
    ax.set_ylabel(r'Flux [neutrons/cm$^2$-s]')
    ax.set_xscale('log')
    #ax.set_yscale('log')
    ax.grid(True, which='both')
    plt.title('Neutrons spectrum in fuel salt')
    plt.savefig("spectrum",dpi=600)
    #
    values = sp.get_tally(name="Mesh").get_slice(scores=['flux']).get_pandas_dataframe()["mean"]
    values = values.values.reshape(500,500)
    fig, ax = plt.subplots()
    pos = ax.imshow(values*fac/salt_vol,
                    aspect='auto',
                    origin='lower')
    cbar = plt.colorbar(pos,ax=ax,label=r'Flux [neutrons/cm$^2$-s]')
    plt.savefig("flux",dpi=600)

    values = sp.get_tally(name="Mesh").get_slice(scores=['fission']).get_pandas_dataframe()["mean"]
    values = values.values.reshape(500,500)
    fig, ax = plt.subplots()
    pos = ax.imshow(values*fac/salt_vol,
                    aspect='auto',
                    origin='lower')
    cbar = plt.colorbar(pos,ax=ax,label=r'Fission [-/cm$^2$-s]')
    plt.savefig("fission",dpi=600)

    values = sp.get_tally(name="Mesh").get_slice(scores=['absorption']).get_pandas_dataframe()["mean"]
    values = values.values.reshape(500,500)
    fig, ax = plt.subplots()
    pos = ax.imshow(values*fac/salt_vol,
                    aspect='auto',
                    origin='lower')
    cbar = plt.colorbar(pos,ax=ax,label=r'Absorption [-/cm$^2$-s]')
    plt.savefig("abs",dpi=600)

    values = sp.get_tally(name="Mesh").get_slice(scores=['scatter']).get_pandas_dataframe()["mean"]
    values = values.values.reshape(500,500)
    fig, ax = plt.subplots()
    pos = ax.imshow(values*fac/salt_vol,
                    aspect='auto',
                    origin='lower')
    cbar = plt.colorbar(pos,ax=ax,label=r'Scatter [-/cm$^2$-s]')
    plt.savefig("scatter",dpi=600)

def depletion(model, mass, power):
    vol=mass/2.3275*1000 # total volume of fuel salt [cm3]
    model.materials[0].volume=vol
    op = openmc.deplete.CoupledOperator(model, normalization_mode = "energy-deposition", chain_file='/home/lorenzo/Documents/ca_depletion_chains/ENDF-B-VIII.0_chain_msr.xml')
    msr = openmc.deplete.msr.MsrContinuous(op,model)
    msr.set_removal_rate('salt', ['Xe','Kr'], 4.067e-5)
    msr.set_removal_rate('salt', ['Se','Nb','Mo','Tc','Ru','Rh','Pd','Ag','Sb','Te'], 8.777e-3)
    integrator = openmc.deplete.CECMIntegrator(op, [5,5,30,30,30,180,95], msr_continuous=msr,
                        timestep_units='d', power=power)
    integrator.integrate(final_step = False)

def control_rod_worth(model):
    keffs = []
    for inch in np.arange(0,51,10):
        cm = inch*2.54
        setattr(cr1_cell, 'translation', [0,0,19.2+cm])
        results=model.run()
        with openmc.StatePoint(res) as sp:
            keffs.append(sp.keff.n)
    plt.figure()
    plt.scatter(np.arange(0,51,10), keffs)
    plt.xlabel('Withdrawn of control rod n. 1 [inch]')
    plt.ylabel('keff')
    plt.savefig('rod_worth')

def triton_adder(mass):

    df = pd.read_csv('ADDER/msre_simple_fuel_mod.csv')
    res = openmc.deplete.Results(f'depletion_results.h5')
    vol=mass/2.3275*1000
    mats = dict()
    mats['fiss'] = ['U235','Pu239','U239','Pu238','H3','Cs137']
    mats['fp'] = ['Ce144','Ce141','I131','Mo99','Xe135','Kr91','Pr145','Sm149']

    fuel_id = [str(mat.id) for mat in res.export_to_materials(0) if mat.name == 'salt'][0]
    t_omc = res.get_atoms(fuel_id,'U235')[0]/3600/24
    t_add = df.loc[np.where(df[df.columns[0]] == 'times')[0][0]][1:].values.astype(float)

    for cat, nucs in mats.items():
        div = round(len(nucs)/2)
        fig,ax = plt.subplots(div,2,figsize=(12,12))
        nuc_index=0
        for i in range(div):
            for j in range(2):
                ax[i,j].plot(t_omc, res.get_atoms(fuel_id,nucs[nuc_index])[1]/vol, marker='x', markersize=9, color='r',label='openmc-msr')
                ax[i,j].plot(t_add, df.loc[np.where(df[df.columns[0]] == nucs[nuc_index])[0][0]][1:].values.astype(float), marker='o', markersize=9, color='b',label='mcnp-adder')
                ax[i,j].legend()
                nuc_index +=1

        for a,n in zip(ax.flat, nucs):
            a.set(xlabel='EFPD [d]', ylabel=f'{n} [atoms/cc]')
        plt.tight_layout()
        plt.savefig(f'msre_openmc-vs-adder_{cat}', dpi=600)

    tt_add=np.unique(t_add)
    args_add = []
    for _t in tt_add:
        args_add.append(np.argwhere(t_add == _t)[0][0])

    means = []
    lim=30
    for nuc in df[df.columns[0]][9:].values:
        try:
            if res.get_atoms('1',nuc)[1].mean()>1e15:
                omc = res.get_atoms('1',nuc)[1][1:]/vol
                add = df.loc[np.where(df[df.columns[0]] == nuc)[0][0]][1:].values.astype(float).take(args_add)[1:]
                diff = (add-omc)/add *100
                diff = diff.mean()
                if not math.isnan(diff):
                    if abs(diff) < lim:
                        means.append(abs(diff))
        except:
            continue
    pd.DataFrame(means).describe()
    params = stats.gamma.fit(means)
    x = np.linspace(0, round(max(means)), 1000)
    pdf = stats.gamma.pdf(x, *params)
    plt.figure()
    plt.plot(x, pdf, label='Gamma func. data fit')
    plt.hist(means,round(max(means))*2,density=True,label='Data')
    plt.legend()
    plt.xlabel('Relative error [%]', weight='bold')
    plt.ylabel('Probability', weight='bold')
    plt.title(f'Relative error distribution for {len(means)} most abundant nuclides in fuel salt below {lim}%',fontsize=9)
    plt.savefig(f'rel_error_below{lim}%', dpi=600)

if __name__ == '__main__':
    mass = 4590 #tot fuel salt mass [kg]
    power = 8e6 #total thermal power [W]
    #run(build(make_tally=True, plot_geom=True), mass, power)
    control_rod_worth(build(make_tally=False, plot_geom=False))
    depletion(build(make_tally=False, plot_geom=False), mass, power)
    #triton_adder(mass)
