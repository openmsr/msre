import numpy as np
import pandas as pd
import math
pi = math.pi

# Perturbations
# SOURCE INSERTION
# No source insertion
sourcedata = np.array([0, 0, 0])
sourcetime = np.array([0, 50, 100])
# % 1 (n/no)/s for 10 seconds
# sourcedata = np.array([0, 10, 0])
# sourcetime = np.array([0, 10, 20])
source = pd.Series(sourcedata, index=sourcetime)

# REACTIVITY INSERTION
# No reactivity insertion
simtime = 10
reactdata = np.array([0, 5E-4])
reacttime = np.array([0, 2500])
# Periodic 60 PCM for 50 seconds
# simtime = 500
# periodic = np.array([[0, 0], [50, 6e-4], [100, 0], [150, -6e-4], [200, 0], [250, 6e-4], [300, 0], [350, -6e-4], [400, 0]])
# reactdata = periodic[:, 1]
# reacttime = periodic[:, 0]
# Step up 60 pcm
# simtime = 1000
# reactdata = np.array([0, 6e-3])
# reacttime = np.array([0, 300])
# # Step down -60 pcm for 10 sec
# simtime = 100
# reactdata = np.array([0, -6e-4])
# reacttime = np.array([0, 50])
# # Pulse 600 pcm for 0.1 sec
# simtime = 30
# reactdata = np.array([0, 6e-3, 0])
# reacttime = np.array([0, 10, 10.1])

react = pd.Series(reactdata, index=reacttime)

ts_max = 1e-1  # maximum timestep (s)

# NEUTRONICS DATA
tau_l = 16.73  # ORNL-TM-0728 %16.44; % (s)
tau_c = 8.46  # ORNL-TM-0728 %8.460; % (s)
P = 8  # Thermal Power in MW ORNL-TM-1070, p.2
n_frac0 = 5  # initial fractional neutron density n/n0 (n/cm^3/s)
Lam = 2.400E-04  # mean generation time ORNL-TM-1070 p.15 U235
# Lam = 4.0E-04;  # mean generation time ORNL-TM-1070 p.15 U233
lam = np.array([1.240E-02, 3.05E-02, 1.11E-01, 3.01E-01, 1.140E+00, 3.014E+00])
beta = np.array([0.000223, 0.001457, 0.001307, 0.002628, 0.000766, 0.00023])  # U235
# beta = np.array([0.00023, 0.00079, 0.00067, 0.00073, 0.00013, 0.00009])  # U233
beta_t = np.sum(beta)  # total delayed neutron fraction MSRE
rho_0 = beta_t-sum(np.divide(beta,1+np.divide(1-np.exp(-lam*tau_l),lam*tau_c))) # reactivity change in going from stationary to circulating fuel
C0 = beta / Lam * (1.0 / (lam - (np.exp(-lam * tau_l) - 1.0) / tau_c))

# Feedback co-efficients
a_f = -8.71E-05  # U235 (drho/°C) fuel salt temperature-reactivity feedback coefficient ORNL-TM-1647 p.3 % -5.904E-05; % ORNL-TM-0728 p. 101 %
a_g = -6.66E-05  # U235 (drho/°C) graphite temperature-reactivity feedback coefficient ORNL-TM-1647 p.3 % -6.624E-05; % ORNL-TM-0728 p.101

# CORE HEAT TRANSFER PARAMETERS
# FUEL PARAMETERS - DONE
vdot_f = 7.5708E-02  # ORNL-TM-0728 % 7.571e-2; % vol. flow rate (m^3/s) ORNL-TM-1647 p.3, ORNL-TM-0728 p.12
rho_f = 2.14647E+03  # (partially enriched U-235)ORNL-TM-0728 p.8 2.243E+03; % (Th-U) density of fuel salt (kg/m^3) ORNL-TM-0728 p.8
W_f = 1.623879934566580e+02  # 1.83085e+02;%vdot_f*rho_f; % 182.78; % calcd from m_dot*cp*delT=P; vdot_f*rho_f; % fuel flow rate (kg/s)
# tau_f_c  = tau_c; % ORNL-TM-0728 % 8.45; % transit time of fuel in core (s) ORNL-TM-1070 p.15, TDAMSRE p.5
m_f = W_f * tau_c  # fuel mass in core (kg)
nn_f = 2  # number of fuel nodes in core model
mn_f = m_f / nn_f  # fuel mass per node (kg)
# cp_f     = 4.2*9/5; % (MJ/deg-C) total fuel heat capacity TDAMSRE p.5
scp_f = 1.9665E-3  # specific heat capacity of fuel salt (MJ/kg-C) ORNL-TM-0728 p.8

# Core Upflow - DONE
v_g = 1.95386  # graphite volume(m^3) ORNL-TM-0728 p. 101
rho_g = 1.860E3  # graphite density (kg/m^3) ORNL-3812 p.77, ORNL-TM-0728 p.87
m_g = v_g * rho_g  # graphite mass (kg)
cp_g = 3.6 * 9 / 5  # TDAMSRE p.5 graphite total heat capacity (MW-s/C) ORNL-TM-1647 p.3
scp_g = 1.773E-3  # cp_g/m_g; % graphite specific heat capacity (MW-s/kg-C) ORNL-TM-1647 p.3
mcp_g1 = m_g * scp_g  # (mass of material x heat capacity of material) of graphite per lump (MW-s/°C)
mcp_f1 = mn_f * scp_f  # (mass of material x heat capacity of material) of fuel salt per lump (MW-s/°C)
mcp_f2 = mn_f * scp_f  # (mass of material x heat capacity of material) of fuel salt per lump (MW-s/°C)
hA_fg = 0.02 * 9 / 5  # (fuel to graphite heat transfer coeff x heat transfer area) (MW/°C) ORNL-TM-1647 p.3, TDAMSRE p.5
k_g = 0.07  # fraction of total power generated in the graphite  ORNL-TM-0728 p.9
k_1 = 0.5  # fraction of heat transferred from graphite which goes to the first fuel lump
k_2 = 0.5  # fraction of heat transferred from graphite which goes to the second fuel lump
k_f = 0.93  # fraction of heat generated in fuel - that generated in the external loop ORNL-TM-0728 p.9
k_f1 = k_f / nn_f  # fraction of total power generated in lump f1
k_f2 = k_f / nn_f  # fraction of total power generated in lump f2

# New node for power deposited in fuel outside the core
k_out = 1 - (k_g + k_f)  # fraction of power generated in fuel in external loop ORNL-TM-0728 p.9
m_out = W_f  # (kg) Mass of node such that resident time is 1 sec (W_f needs to be defined)

# Initial conditions - DONE
Tf_in = 6.3222E+02  # in °C ORNL-TM-1647 p.2
T0_f2 = 6.5727E+02  # 6.5444E+02; % in °C 6.461904761904777e+02; ORNL-TM-1647 p.2
T0_f1 = Tf_in + (T0_f2 - Tf_in) / 2  # 6.405952380952389e+02; in °C
T0_g1 = T0_f1 + (k_g * P / hA_fg)  # 6.589285714285924e+02; in °C
# T0_out = k_out * P / m_out / scp_f + T0_f2  # in °C (scp_f needs to be defined)


# Heat Exchanger - DONE
# Geometry
d_he = 16  # (in) he diameter ORNL-TM-0728 p. 164
h_he = 72  # (in) active height % 96; %(in) he height ORNL-TM-0728 p. 164
od_tube = 0.5  # (in) coolant tube OD ORNL-TM-0728 p. 164
id_tube = od_tube - 2 * 0.042  # (in) coolant tube ID ORNL-TM-0728 p. 164
n_tube = 159  # number of coolant tubes ORNL-TM-0728 p. 164
a_tube = 254 * 144  # (in^2) total area of tubes ORNL-TM-0728 p. 164
l_tube = a_tube / n_tube / (np.pi * od_tube)  # (in) tube length
v_tube = n_tube * np.pi * (od_tube / 2) ** 2 * l_tube  # (in^3) hx shell volume occupied by tubes
v_cool = n_tube * np.pi * (id_tube / 2) ** 2 * l_tube  # (in^3) hx volume occupied by coolant
v_he = (d_he / 2) ** 2 * np.pi * h_he  # (in^3) volume of heat exchanger shell
v_he_fuel = v_he - v_tube  # (in^3) volume available to fuel in shell

# Unit conversions
in_m = 1.63871e-5  # 1 cubic inch = 1.63871e-5 cubic meters

# PRIMARY FLOW PARAMETERS - DONE
W_p = W_f  # fuel flow rate (kg/s)

m_p = v_he_fuel * in_m * rho_f  # fuel mass in PHE (kg)
nn_p = 4  # number of fuel nodes in PHE
mn_p = m_p / nn_p  # fuel mass per node (kg)
cp_p = scp_f  # fuel heat capacity (MJ/(kg-C))

# SECONDARY FLOW PARAMETERS - DONE
vdot_s = 5.36265E-02  # ORNL-TM-0728 p. 164 % 5.236E-02; % coolant volume flow rate (m^3/s) ORNL-TM-1647 p.3
rho_s = 1.922e3  # coolant salt density (kg/m^3) ORNL-TM-0728 p.8
W_s = 1.005793369810108e+02  # vdot_s*rho_s; % calcd from mdot*cp*delT; vdot_s*rho_s; % coolant flow rate (kg/s) ORNL-TM-1647 p.3

m_s = v_cool * in_m * rho_s  # coolant mass in PHE (kg)
nn_s = 4  # number of coolant nodes in PHE
mn_s = m_s / nn_s  # coolant mass per node (kg)
scp_s = 2.39E-3  # cp_s/m_s; % specific heat capacity of coolant (MJ/(kg-C) ORNL-TM-0728 p.8

A_phe = 2.359E+01  # effective area for heat transfer (primary and secondary, m^2) ORNL-TM-0728 p.164

ha_p = 6.480E-01  # heat transfer*area coefficient from primary to tubes (MW/C) ORNL-TM-1647 p.3
ha_s = 3.060E-01  # heat transfer*area coefficient from tubes to secondary (MW/C) ORNL-TM-1647 p.3

# Primary Side
mcp_pn = mn_p * cp_p  # (mass of material x heat capacity of material) of fuel salt per lump in MW-s/°C
hA_pn = ha_p / nn_s  # 3.030; % (primary to tube heat transfer coeff x heat transfer area) in MW/°C

# Tubes - DONE
nn_t = 2  # number of nodes of tubes in the model
rho_tube = 8.7745E+03  # (kg/m^3) density of INOR-8 ORNL-TM-0728 p.20
m_tn = (v_tube - v_cool) * in_m * rho_tube / nn_t  # mass of tubes (kg)
scp_t = 5.778E-04  # specific heat capacity of tubes (MJ/(kg-C)) ORNL-TM-0728 p.20
mcp_tn = m_tn * scp_t  # mass*(heat capacity) of tubes per lump in MW-s/°C

# Secondary Side - DONE
mcp_sn = mn_s * scp_s  # (mass of material x heat capacity of material) of coolant salt per lump in MW-s/°C
hA_sn = ha_s / nn_s  # (tube to secondary heat transfer coeff x heat transfer area) in MW/°C

# Initial conditions - DONE
# Primary nodes
Tp_in = T0_f2  # in °C ORNL-TM-1647 p.2
T0_p4 = Tf_in  # 6.5444E+02; % in °C 6.461904761904777e+02; ORNL-TM-1647 p.2
T0_p1 = Tp_in + (T0_p4 - Tp_in) / 4  # in °C
T0_p2 = Tp_in + 2 * (T0_p4 - Tp_in) / 4  # in °C
T0_p3 = Tp_in + 3 * (T0_p4 - Tp_in) / 4  # in °C

# Secondary nodes
Ts_in = 5.4611E+02  # in °C ORNL-TM-1647 p.2
T0_s4 = 5.7939E+02  # in °C ORNL-TM-1647 p.2
T0_s1 = Ts_in + (T0_s4 - Ts_in) / nn_s  # in °C
T0_s2 = Ts_in + 2 * (T0_s4 - Ts_in) / nn_s  # in °C
T0_s3 = Ts_in + 3 * (T0_s4 - Ts_in) / nn_s  # in °C
# Tube nodes
T0_t1 = (T0_p1 * hA_pn + T0_s3 * hA_sn) / (hA_pn + hA_sn)  # in °C
T0_t2 = (T0_p3 * hA_pn + T0_s1 * hA_sn) / (hA_pn + hA_sn)  # in °C

# Radiator Parameters - DONE

# Initial conditions - DONE
# Primary nodes
Trp_in = T0_s4  # 5.933E+02; % in °C ORNL-TM-1647 p.2
T0_rp = Ts_in  # in °C ORNL-TM-1647 p.2

# Secondary nodes - DONE
Trs_in = 37.78  # (C) air inlet temperature ORNL-TM-1647 p.2
T0_rs = 148.9  # (C) air exit temperature ORNL-TM-1647 p.2

# Radiator Geometry
od_rad = 0.01905  # (m) outer diameter of tubes in the radiator ORNL-TM-0728 p.296
tube_wall_thick = 0.0018288  # (m) thickness of tubes in the radiator ORNL-TM-0728 p.296
id_rad = od_rad - 2 * tube_wall_thick
n_rtubes = 120  # number of tubes in the radiator (rows times tubes per row) ORNL-TM-0728 p.296
l_rtube = 9.144  # (m) length of tubes in the radiator ORNL-TM-0728 p.296
v_rp = pi * (id_rad / 2) ** 2 * l_rtube * n_rtubes  # volume available to salt in the radiator
# v_rtube = pi * (od_rad / 2) ** 2 * l_rtube * n_rtubes - v_rp  # volume of metal in radiator tubes *TUBES NOT MODELED

n_tpr = 12  # number of tubes per row in the radiator matrix
n_row = 10  # number rows in the radiator matrix
tube_space = 0.0381  # (m) spacing between tubes and rows of matrix
v_rs = (n_row * od_rad + (n_row - 1) * tube_space) * (n_tpr * od_rad + (n_tpr - 1) * tube_space) * l_rtube  # volume of air inside radiator

# PRIMARY FLOW PARAMETERS - DONE
W_rp = W_s  # coolant salt flow rate (kg/s)
m_rp = v_rp * rho_s  # coolant salt mass in rad (kg)
nn_rp = 1  # number of coolant salt nodes in the radiator
mn_rp = m_rp / nn_rp  # coolant mass per node (kg)
cp_rp = scp_s  # coolant specific heat capacity (MJ/(kg-C))

# SECONDARY FLOW PARAMETERS - DONE
vdot_rs = 94.389  # ORNL-TM-0728 p. 296; 78.82; % air volume flow rate (m^3/s) ORNL-TM-1647 p.2
rho_rs = 1.1237  # air density (kg/m^3) REFPROP (310K and 0.1MPa)
W_rs = vdot_rs * rho_rs  # air flow rate (kg/s)

m_rs = v_rs * rho_rs  # coolant air mass in rad (kg)
nn_rs = 1  # number of coolant nodes in rad
mn_rs = m_rs / nn_rs  # coolant mass per node (kg)
scp_rs = 1.0085E-3  # (MJ/kg-C) specific heat capacity of air at (air_out+air_in)/2 REFPROP

A_rad = 6.503E1  # (m^2) surface area of the radiator ORNL-TM-0728 p.14
h_roverall = P / A_rad / ((T0_rp + Trp_in) / 2 - (T0_rs + Trs_in) / 2)  # cald as: P/A_rad/((T0_rp+Trp_in)/2-(T0_rs+Trs_in)/2)  3.168E-4; % (MW/m^2-C) polimi thesis

# Primary Side
mcp_rpn = mn_rp * cp_rp  # (mass of material x heat capacity of material) of fuel salt per lump in MW-s/°C
hA_rpn = h_roverall * A_rad / nn_rs  # 3.030; % (primary to secondary heat transfer coeff x heat transfer area) in MW/°C

# Secondary Side - DONE
mcp_rsn = mn_rs * scp_rs  # (mass of material x heat capacity of material) of coolant salt per lump in MW-s/°C
hA_rsn = h_roverall * A_rad / nn_rs  # (tube to secondary heat transfer coeff x heat transfer area) in MW/°C

# Pure time delays between components - DONE
tau_hx_c = 8.67  # (sec) delay from hx to core TDAMSRE p.6
tau_c_hx = 3.77  # (sec) subtracted 1 sec for external loop power generation node resident time; delay from core to fuel hx TDAMSRE p.6
tau_hx_r = 4.71  # (sec) fertile hx to core TDAMSRE p.6
tau_r_hx = 8.24  # (sec) core to fertile hx TDAMSRE p.6

first_val = (rho_0 - beta_t) * n_frac0 / Lam + lam[0] * C0[0] + lam[1] * C0[1] + lam[2] * C0[2] + lam[3] * C0[3] + lam[4] * C0[4] + lam[5] * C0[5]