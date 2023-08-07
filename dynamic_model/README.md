# Model
Non-linear dynamic model of the MSRE based on [Singh et al.](https://www.sciencedirect.com/science/article/pii/S030645491730381X). The model uses non-linear neutron kinetic equations, coupled to 1D linear heat transfer equations in order to model reactor dynamics. The neutron kinetics are expressed in terms of fractional power. As described in Singh et al., "The premise is that reactor power is proportional to neutron density, with all other parameters held fixed." The different components of the reactor, namely the core, heat exchanger and radiator are nodalized such that the dynamics of each is governed by its own system of equations, coupled to the other systems by delayed inlet/outlet terms. 

# Method
The model herein uses SciPy's ODE library. Sample implementation is shown below:

```python
backend = 'dopri5'
r = ode(dydtMSRE).set_integrator(backend,max_step=0.10)
r.integrate(100.00)
```

However, since SciPy's ODE library does not support delay differential equations, the delay terms are stored and handled manually. Since the 'dopri5' method uses adaptive time-stepping, linear interpolation is used for approximating the value of the delay terms near the closest timestep to the delay. 
