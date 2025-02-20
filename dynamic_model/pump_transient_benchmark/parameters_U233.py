import numpy as np



# domain
t0 = 0.0
tf = 1000.00
T = np.arange(t0,tf,0.01)

# REACTIVITY INSERTION
inserted = 1.39e-4 # 1MW
# inserted = 1.96e-4 # 5MW
# inserted = 2.48e-4 # 8MW

# NEUTRONICS DATA
tau_l = 16.73
tau_c = 8.46
# P = 0.1
P = 1  
# P = 5
# P = 8    
n_frac0 = 1  # initial fractional neutron density n/n0
Lam = 4.0E-04
lam = np.array([1.260E-02, 3.370E-02, 1.390E-01, 3.250E-01, 1.130E+00, 2.500E+00])
beta = np.array([0.00023, 0.00079, 0.00067, 0.00073, 0.00013, 0.00009])
beta_t = np.sum(beta)  # total delayed neutron fraction MSRE
rho_0 = beta_t-sum(np.divide(beta,1+np.divide(1-np.exp(-lam*tau_l),lam*tau_c))) # reactivity change in going from stationary to circulating fuel
C0 = beta / Lam * (1.0 / (lam - (np.exp(-lam * tau_l) - 1.0) / tau_c))

# Feedback coefficients
a_f = -11.034E-5
a_g = -05.814E-5

# CORE HEAT TRANSFER PARAMETERS
vdot_f = 7.5708E-02
rho_f = 2.14647E+03
W_f = 1.623879934566580e+02
m_f = W_f * tau_c
nn_f = 2
mn_f = m_f / nn_f
scp_f = 1.9665E-3

# Core Upflow
v_g = 1.95386
rho_g = 1.860E3
m_g = v_g * rho_g
scp_g = 1.773E-3
mcp_g1 = m_g * scp_g
mcp_f1 = mn_f * scp_f
mcp_f2 = mn_f * scp_f
hA_fg = 0.02 * 9 / 5
k_g = 0.07
k_1 = 0.5
k_2 = 0.5
k_f = 0.93
k_f1 = k_f / nn_f
k_f2 = k_f / nn_f

# Heat Exchanger
d_he = 16
h_he = 72
od_tube = 0.5
id_tube = od_tube - 2 * 0.042
n_tube = 159
a_tube = 254 * 144
l_tube = a_tube / n_tube / (np.pi * od_tube)
v_tube = n_tube * np.pi * (od_tube / 2) ** 2 * l_tube
v_cool = n_tube * np.pi * (id_tube / 2) ** 2 * l_tube
v_he = (d_he / 2) ** 2 * np.pi * h_he
v_he_fuel = v_he - v_tube
in_m = 1.63871e-5
W_p = W_f
m_p = v_he_fuel * in_m * rho_f
nn_p = 4
mn_p = m_p / nn_p
cp_p = scp_f
vdot_s = 5.36265E-02
rho_s = 1.922e3
W_s = 1.005793369810108e+02
m_s = v_cool * in_m * rho_s
nn_s = 4
mn_s = m_s / nn_s
scp_s = 2.39E-3
A_phe = 2.359E+01
ha_p = 6.480E-01
ha_s = 3.060E-01
mcp_pn = mn_p * cp_p
hA_pn = ha_p / nn_s
nn_t = 2
rho_tube = 8.7745E+03
m_tn = (v_tube - v_cool) * in_m * rho_tube / nn_t
scp_t = 5.778E-04
mcp_tn = m_tn * scp_t
mcp_sn = mn_s * scp_s
hA_sn = ha_s / nn_s

# Initial conditions
Tf_in = 6.3222E+02
T0_f2 = 6.5727E+02
T0_f1 = Tf_in + (T0_f2 - Tf_in) / 2
T0_g1 = T0_f1 + (k_g * P / hA_fg)
Tp_in = T0_f2
T0_p4 = Tf_in
T0_p1 = Tp_in - (Tp_in - T0_p4) / 4
T0_p2 = Tp_in - 2 * (Tp_in - T0_p4) / 4
T0_p3 = Tp_in - 3 * (Tp_in - T0_p4) / 4
Ts_in = 5.4611E+02
T0_s4 = 5.7939E+02
T0_s1 = Ts_in + (T0_s4 - Ts_in) / nn_s
T0_s2 = Ts_in + 2 * (T0_s4 - Ts_in) / nn_s
T0_s3 = Ts_in + 3 * (T0_s4 - Ts_in) / nn_s
T0_t1 = (T0_p1 * hA_pn + T0_s3 * hA_sn) / (hA_pn + hA_sn)
T0_t2 = (T0_p3 * hA_pn + T0_s1 * hA_sn) / (hA_pn + hA_sn)

# Radiator Parameters
Trp_in = T0_s4
T0_rp = Ts_in
Trs_in = 37.78
T0_rs = 148.9
od_rad = 0.01905
tube_wall_thick = 0.0018288
id_rad = od_rad - 2 * tube_wall_thick
n_rtubes = 120
l_rtube = 9.144
v_rp = np.pi * (id_rad / 2) ** 2 * l_rtube * n_rtubes
n_tpr = 12
n_row = 10
tube_space = 0.0381
v_rs = (n_row * od_rad + (n_row - 1) * tube_space) * (n_tpr * od_rad + (n_tpr - 1) * tube_space) * l_rtube
W_rp = W_s
m_rp = v_rp * rho_s
nn_rp = 1
mn_rp = m_rp / nn_rp
cp_rp = scp_s
vdot_rs = 94.389
rho_rs = 1.1237
W_rs = vdot_rs * rho_rs
m_rs = v_rs * rho_rs
nn_rs = 1
mn_rs = m_rs / nn_rs
scp_rs = 1.0085E-3
A_rad = 6.503E1
h_roverall = P / A_rad / ((T0_rp + Trp_in) / 2 - (T0_rs + Trs_in) / 2)
mcp_rpn = mn_rp * cp_rp
hA_rpn = h_roverall * A_rad / nn_rs
mcp_rsn = mn_rs * scp_rs
hA_rsn = h_roverall * A_rad / nn_rs

# Pure time delays between components
tau_hx_c = 8.67 #+2.145 
tau_c_hx = 3.77 #+2.145
tau_hx_r = 4.71
tau_r_hx = 8.24

