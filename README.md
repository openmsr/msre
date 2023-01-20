# msre
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

detailed cad model of the [msre](https://en.wikipedia.org/wiki/Molten-Salt_Reactor_Experiment) (molten salt reactor experiment), operated by oak ridge national laboratory 1965-69.

## msre core
<img src="core/docs/msre.png" width="500" height="500"/>
[core/msrecore.pdf](core/docs/msrecore.pdf) lists reference of the msre core design, documented in the old msre reports and located [here](https://github.com/openmsr/msr-archive/blob/master/README.md).

## msre step files
.step files of entire msre assembly and control rod.

## h5m
.h5m surface mesh of the previous step files for OpenMC simulation.
**Note:** the h5m files are generated with Coreform Cubit with a surface tolerance of 1e-2 cm.

## docker msre
docker conatiner which includes not only support for OpenMC with DAGMC/MOAB, embree, and double_down libraries, but also a Jupyter notebook server.

## openmc notebooks
openmc models of the msre in form of jupyter notebooks. Examples include:
- msre cad with settable control rods
- msre isothermal temperature coefficient calculation
- msre depletion analysis with fission products removal

## prerequisites
### openmc
[openmc](https://docs.openmc.org/en/stable/) automated source installation scripts for linux can be found [here](https://github.com/openmsr/openmc_install_scripts)
### CAD_to_openMC
[CAD_to_openMC](https://github.com/openmsr/CAD_to_openMC) is an open-source package to convert CAD geometry (in the form of '.step' files) into an openmc-readable h5m file


## msre heat exchangers

open-access [master's thesis](https://ltu.diva-portal.org/smash/get/diva2:1546993/FULLTEXT01.pdf) produced by Malcolm Akner about simulations of the heat exchangers of the msre, titled:

Validating results from the Molten Salt Reactor Experiment by use of turbulent CFD simulations
- A study of a modified U-tube shell-and-tube primary heat exchanger and radiator with molten salts

### msre primary heat exchanger
![](heatexchanger/docs/phexcadmodel.png)

[onshape primary heat exchanger cad model](https://cad.onshape.com/documents/03be2f510296a2e264886390/w/8cfbca3b7b9682dd4e53a998/e/54728fd981a1b4f5594c73d6), open to copy and use freely. chapter 4.1 in the [thesis](https://ltu.diva-portal.org/smash/get/diva2:1546993/FULLTEXT01.pdf) mentioned above covers the cad construction details extensively with references to original msre reports.

![](heatexchanger/docs/phexflowpaths.png)

[simscale primary heat exchanger simulation model](https://www.simscale.com/projects/MalcolmAkner/phex_-_final_version/). simulation results for primary heat exchanger can be viewed in chapter 6.2.2.1 of the [thesis](https://ltu.diva-portal.org/smash/get/diva2:1546993/FULLTEXT01.pdf), with comparisons to msre data in chapter 7.1.1.

![](heatexchanger/docs/phexreal.png)

primary heat exchanger produced and installed in the msre.


### msre radiator
![](heatexchanger/docs/radiatorcadmodel.png)

[onshape radiator cad model](https://cad.onshape.com/documents/bf944323ed6a82e05924078c/w/2a25d73c5a3a66824d2d5fbd/e/a83d5535602a053216fedff4) open to copy and use freely. chapter 4.2 of the [thesis](https://ltu.diva-portal.org/smash/get/diva2:1546993/FULLTEXT01.pdf) covers the cad construction details with origianl references to msre reports.

![](heatexchanger/docs/radiatorreal.png)

radiator produced and installed in the msre.

---

please contact [me](https://github.com/aslakstubsgaard) if you want to contribute.
note that this work and the cad models are under the GNU General Public License v3.0

---
