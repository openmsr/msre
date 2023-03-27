import openmc
import os
import re
import matplotlib.pyplot as plt
import seaborn as sns
regex = re.compile(r'(\d+|\s+)')

path_to_results = 'depletion_results.h5'
save_dir = Path(os.path.realpath(path_to_results)).parent
mats = {mat.name: mat.id for mat in openmc.material.Materials.from_xml(save_dir / 'materials.xml')
                    if mat.depletable}

results = openmc.deplete.Results(path_to_results)
t, keff = results.get_keff()
n_xe = 0
n_kr = 0
for nuc,_ in openmc.data.isotopes('Xe'):
    n_xe += results.get_atoms(str(mats['salt']), nuc)[1]
for nuc,_ in openmc.data.isotopes('Kr'):
    n_kr += results.get_atoms(str(mats['salt']), nuc)[1]

# Let's convert time from sec to days
t /= (3600 * 24)

plt.figure()
ax = plt.subplot()
k1, = ax.plot(t, [k[0] for k in keff], '--', c='red', label='keff')
ax1 = ax.twinx()
n2, = ax1.plot(t, n_xe, c='green', label='Xe')
n4, = ax1.plot(t, n_kr, c='blue', label='Kr')
ax.set_xlabel('Time[d]')
ax.set_ylabel(r'$k_{eff}$', color='r')
ax.tick_params(axis='y', colors='red')
ax1.set_yscale('log')
ax1.set_ylabel('Nuclides [atoms]')
ax1.legend(handles=[k1, n2, n4])
plt.savefig(f'{save_dir}/keff', dpi=600)

# Microscopic absorption cross section at 0.0253 eV
xs_xe135 = 2664214.0
xs_u235 = 686.006994850397
_, n_xe135 = results.get_atoms(str(mats['salt']), 'Xe135')
_, n_u235 = results.get_atoms(str(mats['salt']), 'U235')
# Poison fraction
pf = (xs_xe135*n_xe135)/(xs_u235*n_u235)*100
plt.figure()
plt.plot(t, pf)
plt.xlabel('Time [d]')
plt.ylabel('Xe posion fraction [%]')
plt.savefig(f'{save_dir}/fission_products', dpi=600)

inventory = dict()
for nuc,_ in openmc.data.isotopes('U'):
    inventory[nuc] = results.get_atoms(str(mats['salt']), nuc)[1] / openmc.data.AVOGADRO * openmc.data.atomic_mass(nuc) / 1000

for nuc in ['Pu238','Pu239','Pu240','Pu241','Pu242']:
    inventory[nuc] = results.get_atoms(str(mats['salt']), nuc)[1] / openmc.data.AVOGADRO * openmc.data.atomic_mass(nuc) / 1000

plt.figure()
for nuc, mass in inventory.items():
    plt.plot(t, mass, label=nuc)
plt.xlabel('Time [y]')
plt.ylabel('Mass [g]')
plt.yscale('log')
plt.ylim(1e-5)
plt.legend()
plt.savefig(f'{save_dir}/inventory', dpi=600)

# All nuclides present in the fuel at last time-step
nucs = results.export_to_materials(-1)[0].get_nuclides()

# Let's begin by making some useful groupings
gaseos = ['H', 'He', 'Ne', 'Ar', 'Kr', 'Xe', 'Rn'] #gaseous fission products
noble_metals = ['Se','Nb','Mo','Tc','Ru','Rh','Pd','Ag','Sb','Te'] # noble metals fission products
metals = ['Cr','Mn','Fe','Co','Ni','Cu','Zn','Hf','Zr','W',]
halogens = ['F','Cl','Br','I','At']
alkali_metals = ['Li','Na','K','Rb','Cs']
alkali_earths= ['Be','Mg','Ca','Sr','Ba','Ra']
lanthanides  = ['Y','La','Ce','Pr','Nd','Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb','Lu']
m_a = ['Ac','Th','Pa','Np','Am','Cm','Bk','Cf','Es','Fm','Md','No','Lr']

# Get fissile nuclides in the fuel, based on Ronen's rule for determining fissile isotopes
fissile = []

for nuc in nucs:
    elm = regex.split(nuc)[0]
    a = round(openmc.data.atomic_mass(nuc))
    z = openmc.data.ATOMIC_NUMBER[elm]
    if 90 <= z <= 100:
        ronen = 2*z -(a-z)
        if ronen in [41,43,45]:
            fissile.append(nuc)

# Calculate totat absorption rate of fissile nuclides
tot_abs_rate = 0
for nuc in fissile:
    tot_abs_rate += results.get_reaction_rate(str(mats['salt']), nuc, 'fission')[1]
    tot_abs_rate += results.get_reaction_rate(str(mats['salt']), nuc, '(n,gamma)')[1]

nuclides_stack = dict()
groups_stack = {'Gaseos':0, 'Noble metals':0, 'Metals':0, 'Halogens':0 , 'Alkali metals':0, 'Alkali earths':0, 'Lanthanides':0, 'MA':0, 'Others':0}
for nuc in nucs:

    if regex.split(nuc)[0] in ['U','Pu']:
        nuclides_stack[nuc] = results.get_reaction_rate(str(mats['salt']), nuc, 'fission')[1]
        nuclides_stack[nuc] += results.get_reaction_rate(str(mats['salt']), nuc, '(n,gamma)')[1]
        nuclides_stack[nuc] /= tot_abs_rate


    elif regex.split(nuc)[0] in gaseos:
        groups_stack['Gaseos'] += results.get_reaction_rate(str(mats['salt']), nuc, 'fission')[1]
        groups_stack['Gaseos'] += results.get_reaction_rate(str(mats['salt']), nuc, '(n,gamma)')[1]
    elif regex.split(nuc)[0] in noble_metals:
        groups_stack['Noble metals'] += results.get_reaction_rate(str(mats['salt']), nuc, 'fission')[1]
        groups_stack['Noble metals'] += results.get_reaction_rate(str(mats['salt']), nuc, '(n,gamma)')[1]
    elif regex.split(nuc)[0] in metals:
        groups_stack['Metals'] += results.get_reaction_rate(str(mats['salt']), nuc, 'fission')[1]
        groups_stack['Metals'] += results.get_reaction_rate(str(mats['salt']), nuc, '(n,gamma)')[1]
    elif regex.split(nuc)[0] in halogens:
        groups_stack['Halogens'] += results.get_reaction_rate(str(mats['salt']), nuc, 'fission')[1]
        groups_stack['Halogens'] += results.get_reaction_rate(str(mats['salt']), nuc, '(n,gamma)')[1]
    elif regex.split(nuc)[0] in alkali_metals:
        groups_stack['Alkali metals'] += results.get_reaction_rate(str(mats['salt']), nuc, 'fission')[1]
        groups_stack['Alkali metals'] += results.get_reaction_rate(str(mats['salt']), nuc, '(n,gamma)')[1]
    elif regex.split(nuc)[0] in alkali_earths:
        groups_stack['Alkali earths'] += results.get_reaction_rate(str(mats['salt']), nuc, 'fission')[1]
        groups_stack['Alkali earths'] += results.get_reaction_rate(str(mats['salt']), nuc, '(n,gamma)')[1]
    elif regex.split(nuc)[0] in lanthanides:
        groups_stack['Lanthanides'] += results.get_reaction_rate(str(mats['salt']), nuc, 'fission')[1]
        groups_stack['Lanthanides'] += results.get_reaction_rate(str(mats['salt']), nuc, '(n,gamma)')[1]
    elif regex.split(nuc)[0] in m_a:
        groups_stack['MA'] += results.get_reaction_rate(str(mats['salt']), nuc, 'fission')[1]
        groups_stack['MA'] += results.get_reaction_rate(str(mats['salt']), nuc, '(n,gamma)')[1]
    else:
        groups_stack['Others'] += results.get_reaction_rate(str(mats['salt']), nuc, 'fission')[1]
        groups_stack['Others'] += results.get_reaction_rate(str(mats['salt']), nuc, '(n,gamma)')[1]

# Divide each array by the total absorption reaction rate of fissile nuclides
for g in groups_stack.keys():
    groups_stack[g] /= tot_abs_rate

# Sort dictionary groups
groups_stack=dict(reversed(sorted(groups_stack.items(), key=lambda item: item[1][len(item)])))

# Create red color palette for groups_stack
colors = list(reversed(sns.color_palette("Reds", len(groups_stack))))

# Order uranium series
u_series = {key:value for key,value in nuclides_stack.items() if key.startswith('U')}
u_series = dict(reversed(sorted(u_series.items(), key=lambda item: item[1][len(item)])))
# Create green color palette for Uranium isotopes
colors += list(reversed(sns.color_palette("Greens", len(u_series))))

# Order plutionium series
pu_series = {key:value for key,value in nuclides_stack.items() if key.startswith('Pu')}
pu_series = dict(reversed(sorted(pu_series.items(), key=lambda item: item[1][len(item)])))
# Create blue color palette for plutonium isotopes
colors += list(reversed(sns.color_palette("Blues", len(pu_series))))
# Add uramium and plutonium series to the stack
groups_stack.update(u_series)
groups_stack.update(pu_series)

plt.figure(figsize=(15,10))
plt.stackplot(t, groups_stack.values(), labels=groups_stack.keys(),
              edgecolor="black", linewidth=0.5,colors=colors, alpha=0.8)
handles, labels = plt.gca().get_legend_handles_labels()
legend = plt.legend([handles[idx] for idx in list(reversed(np.arange(0,len(handles),1)))],
                            [labels[idx] for idx in list(reversed(np.arange(0,len(handles),1)))],
                            bbox_to_anchor=(1.05,1), loc='upper left',
                            borderaxespad=0, ncol=2,
                            fontsize=15)
plt.xlabel('Time [d]',weight='bold',fontsize=17)
plt.title('Neutrons absorption distribution per neutron absorbed in fissile isotopes',
                   weight='bold', fontsize=17)
plt.xticks(fontsize=13)
plt.yticks(fontsize=13)
plt.tight_layout()
plt.savefig(f'{save_dir}/neutrons', dpi=600)
