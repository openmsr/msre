from parameters import *
import numpy as np
from scipy.integrate import ode
import matplotlib.pyplot as plt
import copy 

def dydtMSRE(t,y,delays,rho_ext):

    '''
    Returns derivative of state vector y; y' or dy/dt, of the MSRE system.

    The y vector contains the following:

        T_in_rc = Inlet temperature (°C) of radiator coolant, will be equal to 
                  the outlet temperature of the heat exchanger (T_hc4) plus 
                  relevant time delay

        T_out_rc = Outlet temperature (°C) of radiator coolant

        T_in_air = Inlet temperature (°C) of air in radiator

        T_out_air = Outlet temperature (°C) of air in radiator

        T_in_hf = Inlet temperature (°C) of heat exchanger fuel, will be equal 
                  to the outlet temperature of the core (T_cf2) plus relevant 
                  time delay

        T_hf* = Temperature (°C) of heat exchanger fuel node *

        T_ht* = Temperature (°C) of heat exchanger tube node *

        T_hc* = Temperature (°C) of heat exchanger coolant node *

        T_in_cf = Inlet temperature (°C) of core fuel, will be equal to the 
                  outlet temperature of the heat exchanger (T_hf4) plus 
                  relevant time dela

        S = neutro source perturbation term 

        rho_fb = feedback reactivity (from fuel and graphite temperatures)

        rho_ext = external reactivity (reactivity insertion)

        rho_tot = total reactivity = rho_0 + rho_fb + rho_ext (rho_0 = 
                  steady-state reactivity, constant)

        n = neutron density n(t)

        C* = precursor concentration of group *

        T_cg = Temperature (°C) of core graphite node 

        T_cf* = Temperature (°C) of core fuel node *  

        *_delay = Parameter at time t = t - delay

    Other parameters are defined in parameters.py
    '''

    # unpack state variables 
    T_out_rc, T_out_air, T_hf1, T_hf2, T_hf3, T_hf4, T_ht1, T_ht2, T_hc1,     \
    T_hc2, T_hc3, T_hc4, n, C1, C2, C3, C4, C5, C6, T_cg, T_cf1, T_cf2 = y
    
    # delay terms
    T_out_rc_delay, T_hc4_delay, T_cf2_delay, T_hf4_delay, C1_delay,          \
    C2_delay, C3_delay, C4_delay, C5_delay, C6_delay = delays

    # reactivity
    rho = (a_f/2)*((-T0_f1+T_cf1)+(-T0_f2+T_cf2)) + a_g*(-T0_g1+T_cg) + rho_ext

    # derivatives 
    dydt = [ 
    (W_rp/mn_rp)*(T_hc4_delay-T_out_rc) + (hA_rpn/mcp_rpn)*(T_out_air-T_out_rc),                  # T_out_rc                                                                                   # T_out_rc
    -((W_rs/mn_rs)+(hA_rsn/mcp_rsn))*T_out_air + (hA_rsn/mcp_rsn)*T_out_rc + (W_rs/mn_rs)*Trs_in, # T_out_air
    -((W_p/mn_p)+(hA_pn/mcp_pn))*T_hf1 + (hA_pn/mcp_pn)*T_ht1 + (W_p/mn_p)*T_cf2_delay,           # T_hf1
    (W_p/mn_p)*(T_hf1-T_hf2) + (hA_pn/mcp_pn)*(T_ht1-T_hf1),                                      # T_hf2
    -((W_p/mn_p)+(hA_pn/mcp_pn))*T_hf3 + (hA_pn/mcp_pn)*T_ht2 + (W_p/mn_p)*T_hf2,                 # T_hf3
    (W_p/mn_p)*(T_hf3-T_hf4) + (hA_pn/mcp_pn)*(T_ht2-T_hf3),                                      # T_hf4
    (2*hA_pn/mcp_tn)*(T_hf1-T_ht1) + (2*hA_sn/mcp_tn)*(T_hc3-T_ht1),                              # T_ht1
    (2*hA_pn/mcp_tn)*(T_hf3-T_ht2) + (2*hA_sn/mcp_tn)*(T_hc1-T_ht2),                              # T_ht2
    -((W_s/mn_s)+(hA_sn/mcp_sn))*T_hc1 + (hA_sn/mcp_sn)*T_ht2 + (W_s/mn_s)*T_out_rc_delay,        # T_hc1
    (W_s/mn_s)*(T_hc1-T_hc2) + (hA_sn/mcp_sn)*(T_ht2-T_hc1),                                      # T_hc2
    -((W_s/mn_s)+(hA_sn/mcp_sn))*T_hc3 + (hA_sn/mcp_sn)*T_ht1 + (W_s/mn_s)*T_hc2,                 # T_hc3
    (W_s/mn_s)*(T_hc3-T_hc4) + (hA_sn/mcp_sn)*(T_ht1-T_hc3),                                      # T_hc4
    (rho-beta_t)*n/Lam+lam[0]*C1+lam[1]*C2+lam[2]*C3+lam[3]*C4+lam[4]*C5+lam[5]*C6,               # n (no source insertion)
    n*beta[0]/Lam-lam[0]*C1-C1/tau_c+C1_delay*np.exp(-lam[0]*tau_l)/tau_c,                        # C1
    n*beta[1]/Lam-lam[1]*C2-C2/tau_c+C2_delay*np.exp(-lam[1]*tau_l)/tau_c,                        # C2
    n*beta[2]/Lam-lam[2]*C3-C3/tau_c+C3_delay*np.exp(-lam[2]*tau_l)/tau_c,                        # C3
    n*beta[3]/Lam-lam[3]*C4-C4/tau_c+C4_delay*np.exp(-lam[3]*tau_l)/tau_c,                        # C4
    n*beta[4]/Lam-lam[4]*C5-C5/tau_c+C5_delay*np.exp(-lam[4]*tau_l)/tau_c,                        # C5
    n*beta[5]/Lam-lam[5]*C6-C6/tau_c+C6_delay*np.exp(-lam[5]*tau_l)/tau_c,                        # C6
    (hA_fg/mcp_g1)*(T_cf1 - T_cg) + k_g*P*n/mcp_g1,                                               # T_cg
    W_f/mn_f*(T_hf4_delay-T_cf1) + (k_f1*P*n/mcp_f1) + (hA_fg*k_1*(T_cg - T_cf1)/mcp_f1),         # T_cf1   
    W_f/mn_f*(T_cf1 - T_cf2) + (k_f2*P*n/mcp_f2) + (hA_fg*k_2*(T_cg - T_cf1)/mcp_f2),             # T_cf2
       
    ]
    return dydt

def get_tIdx(t,tao,timeVec):
    '''
    Returns index of time t = t-tau
    '''
    td = t-tao
    diff_min = 999999.9999999
    idx = 0
    for t in enumerate(timeVec):
        diff = abs(td-t[1])
        if (diff<diff_min):
            diff_min = abs(td-t[1])
            idx = t[0]
    return idx, timeVec[idx]-td

def main():
    '''
    Sets initial conditions and calls the solver 
    '''

    # initial conditions
    y0 = [T0_rp, T0_rs, T0_p1,T0_p2, T0_p3, T0_p4, T0_t1, T0_t2, T0_s1, T0_s2, 
          T0_s3, T0_s4, n_frac0, C0[0], C0[1], C0[2], C0[3], C0[4], C0[5], 
          T0_g1, T0_f1, T0_f2]

    # initial delay terms 
    d0 = [T0_rp, T0_rs, T0_f2, T0_p4, C0[0], C0[1], C0[2], C0[3], C0[4], C0[5]] 

    # solver
    backend = 'dopri5'
    r = ode(dydtMSRE).set_integrator(backend,max_step=0.10)

    sol_interim = []
    def solout(t, y):
        sol_interim.append([t, *y])
    r.set_solout(solout)

    # timing parameters
    t0 = 0.0
    t_start = t0
    t_stop = 500.00

    # solution 
    sol = []
    
    # step-reactivity insertion
    t_insert = 2500.00
    insert = 1.0e-4
    if (t_start>=t_insert):
        rho_ext = insert
    else:
        rho_ext = 0.0

    # delay parameters
    d_terms = []
    derivs = [dydtMSRE(t0,y0,d0,rho_ext)]

    # takes one step at a time and then accounts for delay terms
    i = 0
    while (t_start < t_stop):
        # take one step
        if (i == 0):
            t_start = t0
            r.set_initial_value(y0,t0).set_f_params(d0,0.0)
            r.integrate(1.0)
            sol.append(sol_interim[0])
            sol.append(sol_interim[1])
        else:
            t_start = sol[-1][0]
            if (t_start>=t_insert):
                rho_ext = insert
            else:
                rho_ext = 0.0
            r.set_initial_value(y_next,t_start).set_f_params(d_new,rho_ext)
            r.integrate(t_start+1.0)
            sol.append(sol_interim[1])
            derivs.append(dydtMSRE(t_start,y_next,d_new,rho_ext))

        # account for delays, linear interpolation for time differences 
        # core fuel inlet
        d_new = [0]*10
        idx_cf_in = 0
        dt_cf = 0.0
        if (t_start > tau_hx_c):
            idx_cf_in, dt_cf = get_tIdx(t_start,tau_hx_c,[s[0] for s in sol])
        d_new[3] = sol[idx_cf_in][6] + dt_cf*derivs[idx_cf_in][5]

        # heat exchanger fuel inlet
        idx_hf_in = 0
        dt_hf = 0.0
        if (t_start > tau_c_hx):
            idx_hf_in, dt_hf = get_tIdx(t_start,tau_c_hx,[s[0] for s in sol])
        d_new[2] = sol[idx_hf_in][22] + dt_hf*derivs[idx_hf_in][21]

        # heat exchanger coolant inlet
        idx_hc_in = 0
        dt_hc = 0.0
        if (t_start > tau_r_hx):
            idx_hc_in, dt_hc = get_tIdx(t_start,tau_r_hx,[s[0] for s in sol])
        d_new[0] = sol[idx_hc_in][1] + dt_hc*derivs[idx_hc_in][0]

        # radiator coolant inlet
        idx_rc_in = 0
        dt_rc = 0.0
        if (t_start > tau_hx_r):
            idx_rc_in, dt_rc = get_tIdx(t_start,tau_hx_r,[s[0] for s in sol])
        d_new[1] = sol[idx_rc_in][12] + dt_rc*derivs[idx_rc_in][11]

        # precursors
        idx_c = 0
        dt_c = 0.0
        if (t_start > tau_l):
            idx_c, dt_c = get_tIdx(t_start,tau_l,[s[0] for s in sol])
        d_new[4] = sol[idx_c][14] + dt_c*derivs[idx_c][13]
        d_new[5] = sol[idx_c][15] + dt_c*derivs[idx_c][14]
        d_new[6] = sol[idx_c][16] + dt_c*derivs[idx_c][15]
        d_new[7] = sol[idx_c][17] + dt_c*derivs[idx_c][16]
        d_new[8] = sol[idx_c][18] + dt_c*derivs[idx_c][17]
        d_new[9] = sol[idx_c][19] + dt_c*derivs[idx_c][18]

        d_terms.append(d_new)

        # initial condiiton for next step
        y_next = sol[-1][1:]

        # empty interim solution
        sol_interim = []

        # display progress
        #print(f"{t_start}")

        i += 1

    # plot single parameter
    #of_interest = 13
    #ti = [s[0] for s in sol]
    #oi = [s[of_interest] for s in sol]
    #print(type(ti[0]))
    #plt.plot(ti,oi)
    #plt.show()

    # check delay behavior
    # for i in range(len(sol)-1):
    #    print(f"t: {sol[i][0]}, hf4: {sol[i][6]}, c1_delay: {d_terms[i][3]}")

    # write output data
    output_filename = f"sim_out_{t_stop}_{P}"
    results = open(output_filename,'w+')
    for k in range(len(sol)):
        for col in range(len(sol[0])):
            results.write(f"{sol[k][col]} ")
        results.write("\n")
    
    return None 

main()