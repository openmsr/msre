# Imports
from parameters_U233 import *
import numpy as np
import matplotlib.pyplot as plt
from jitcdde import t
from msrDynamics.system_objects import Node, System
import pandas as pd
import sympy as sp
from concurrent.futures import ProcessPoolExecutor
from scipy.signal import find_peaks

f_range = np.logspace(-2, 1, num=100)
# tau_l = 2 * tau_l
# tau_c = 2*tau_c
# tau_c_hx = 2 * tau_c_hx
# tau_hx_c = 2 * tau_hx_c
def process_frequency(f):

    MSRE = System()

    # radiator
    T_out_rc = Node(m = mn_rp, scp = mcp_rpn/mn_rp, W = W_rp, y0 = T0_rp)
    T_out_air = Node(m = mn_rs, scp = mcp_rsn/mn_rs, W = W_rs, y0 = T0_rs)

    # heat exchanger
    T_hf1 = Node(m = mn_p, scp = mcp_pn/mn_p, W = W_p, y0 = T0_p1)
    T_hf2 = Node(m = mn_p, scp = mcp_pn/mn_p, W = W_p, y0 = T0_p2)
    T_hf3 = Node(m = mn_p, scp = mcp_pn/mn_p, W = W_p, y0 = T0_p3)
    T_hf4 = Node(m = mn_p, scp = mcp_pn/mn_p, W = W_p, y0 = T0_p4)
    T_ht1 = Node(m = m_tn, scp = scp_t, y0 = T0_t1)
    T_ht2 = Node(m = m_tn, scp = scp_t, y0 = T0_t2)
    T_hc1 = Node(m = mn_s, scp = mcp_sn/mn_s, W = W_s, y0 = T0_s1)
    T_hc2 = Node(m = mn_s, scp = mcp_sn/mn_s, W = W_s, y0 = T0_s2)
    T_hc3 = Node(m = mn_s, scp = mcp_sn/mn_s, W = W_s, y0 = T0_s3)
    T_hc4 = Node(m = mn_s, scp = mcp_sn/mn_s, W = W_s, y0 = T0_s4)

    # core 
    n = Node(y0 = n_frac0)
    C1 = Node(y0 = C0[0])
    C2 = Node(y0 = C0[1])
    C3 = Node(y0 = C0[2])
    C4 = Node(y0 = C0[3])
    C5 = Node(y0 = C0[4])
    C6 = Node(y0 = C0[5])
    rho = Node(y0 = 0.0)

    # add reactivity input
    r = 1e-5
    def rho_insert(t):
        return r*sp.sin(f*t)

    rho_ext = MSRE.add_input(rho_insert, T)

    T_cg = Node(m = mcp_g1/scp_g, scp = scp_g, y0 = T0_g1)
    T_cf1 = Node(m = mn_f, scp = scp_f, W = W_f, y0 = T0_f1)
    T_cf2 = Node(m = mn_f, scp = scp_f, W = W_f, y0 = T0_f2)

    MSRE.add_nodes([T_out_rc,T_out_air,T_hf1,T_hf2,T_hf3,T_hf4,T_ht1,T_ht2,T_hc1,
                T_hc2,T_hc3,T_hc4,n,C1,C2,C3,C4,C5,C6,T_cg,T_cf1,T_cf2,rho])

    # dynamics 

    # radiator
    T_out_rc.set_dTdt_advective(source = T_hc4.y(t-tau_hx_r))
    T_out_rc.set_dTdt_convective(source = [T_out_air.y()], hA = [hA_rpn])

    T_out_air.set_dTdt_advective(source = Trs_in)
    T_out_air.set_dTdt_convective(source = [T_out_rc.y()], hA = [hA_rsn])

    # heat exchanger
    T_hf1.set_dTdt_advective(source = T_cf2.y(t-tau_c_hx))
    T_hf1.set_dTdt_convective(source = [T_ht1.y()], hA = [hA_pn])

    T_hf2.set_dTdt_advective(source = T_hf1.y())
    T_hf2.dTdt_convective = T_hf1.dTdt_convective

    T_hf3.set_dTdt_advective(source = T_hf2.y())
    T_hf3.set_dTdt_convective(source = [T_ht2.y()], hA = [hA_pn])

    T_hf4.set_dTdt_advective(source = T_hf3.y())
    # T_hf4.set_dTdt_convective(source = [T_ht2.y()], hA = [hA_pn])
    T_hf4.dTdt_convective = T_hf3.dTdt_convective

    # T_ht1.set_dTdt_convective(source = [T_hf1.y(),T_hf2.y(),T_hc3.y(),T_hc4.y()], hA = [hA_pn,hA_pn,hA_sn,hA_sn])
    # T_ht2.set_dTdt_convective(source = [T_hf3.y(),T_hf4.y(),T_hc1.y(),T_hc2.y()], hA = [hA_pn,hA_pn,hA_sn,hA_sn])
    T_ht1.set_dTdt_convective(source = [T_hf1.y(),T_hf1.y(),T_hc3.y(),T_hc3.y()], hA = [hA_pn,hA_pn,hA_sn,hA_sn])
    T_ht2.set_dTdt_convective(source = [T_hf3.y(),T_hf3.y(),T_hc1.y(),T_hc1.y()], hA = [hA_pn,hA_pn,hA_sn,hA_sn])

    T_hc1.set_dTdt_advective(source = T_out_rc.y(t-tau_r_hx))
    T_hc1.set_dTdt_convective(source = [T_ht2.y()], hA = [hA_sn])

    T_hc2.set_dTdt_advective(source = T_hc1.y())
    T_hc2.dTdt_convective = T_hc1.dTdt_convective

    T_hc3.set_dTdt_advective(source = T_hc2.y())
    T_hc3.set_dTdt_convective(source = [T_ht1.y()], hA = [hA_sn])

    T_hc4.set_dTdt_advective(source = T_hc3.y())
    T_hc4.dTdt_convective = T_hc3.dTdt_convective

    # core
    n.set_dndt(r = rho.y()+rho_ext, beta_eff = beta_t, Lambda = Lam, lam = lam, C = [C1.y(),C2.y(),C3.y(),C4.y(),C5.y(),C6.y()])
    C1.set_dcdt(n.y(),beta[0],Lam,lam[0],tau_c,tau_l)
    C2.set_dcdt(n.y(),beta[1],Lam,lam[1],tau_c,tau_l)
    C3.set_dcdt(n.y(),beta[2],Lam,lam[2],tau_c,tau_l)
    C4.set_dcdt(n.y(),beta[3],Lam,lam[3],tau_c,tau_l)
    C5.set_dcdt(n.y(),beta[4],Lam,lam[4],tau_c,tau_l)
    C6.set_dcdt(n.y(),beta[5],Lam,lam[5],tau_c,tau_l)

    T_cg.set_dTdt_convective(source = [T_cf1.y()], hA = [hA_fg])
    T_cg.set_dTdt_internal(source = n.y(), k = k_g*P)

    T_cf1.set_dTdt_advective(source = T_hf4.y(t-tau_hx_c))
    T_cf1.set_dTdt_convective(source = [T_cg.y()], hA = [k_1*hA_fg])
    T_cf1.set_dTdt_internal(source = n.y(), k = k_f1*P)

    T_cf2.set_dTdt_advective(source = T_cf1.y())
    T_cf2.dTdt_convective = T_cf1.dTdt_convective
    T_cf2.set_dTdt_internal(source = n.y(), k = k_f2*P)

    rho.set_drdt(sources = [T_cf1.dydt(), T_cf2.dydt(), T_cg.dydt()], coeffs = [a_f/2,a_f/2,a_g])

    MSRE.solve(T)

    i_out = [i for i in range(len(T)) if T[i] >= 500]
    n0 = n.y_out[i_out[0]-25]
    n_out = np.array(n.y_out)[i_out]

    # calculate output amplitude
    peaks, _ = find_peaks(n_out)
    troughs, _ = find_peaks(-n_out)
    amplitude = (np.mean(n_out[peaks]) - np.mean(n_out[troughs]))/2

    # calculate Gain
    input_amplitude = r  
    gain = amplitude / (input_amplitude*n0)

    # calculate Phase Shift
    peak_times = T[i_out][peaks]
    input_period = 2 * np.pi / f
    input_signal = [i[0] for i in MSRE.input.get_state(T)]
    input_peaks, _ = find_peaks(input_signal)
    time_differences = [abs(T[input_peaks[i]] - peak_times[0]) for i in range(len(input_peaks))]
    closest_peak_index = np.argmin(time_differences)
    closest_peak_time = T[input_peaks[closest_peak_index]]
    
    # calculate Phase Shift
    phase_shift = 360*(closest_peak_time-peak_times[0])/input_period

    if (phase_shift>180):
        phase_shift = -360+phase_shift
    elif (phase_shift<-180):
        phase_shift = 360-phase_shift
           
    return f, gain, phase_shift

with ProcessPoolExecutor() as executor:
    results = list(executor.map(process_frequency, f_range))

# Process the results
results_df = pd.DataFrame(results, columns=['Frequency', 'Gain', 'Phase Shift'])

# Write to CSV file
# csv_filename = f"frequency_response_results_{P}_MW_double_tau.csv"
csv_filename = f"frequency_response_results_{P}_MW.csv"
results_df.to_csv(csv_filename, index=False)

print(f"Results written to {csv_filename}")