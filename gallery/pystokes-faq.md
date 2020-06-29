# PyStokes Wiki

Welcome to the PyStokes Wiki. 

1. Please read the PyStokes [Frequently Asked Questions](https://github.com/rajeshrinet/pyross/wiki/FAQ-on-PyStokes) if you want to use the library for your own research.  
2. See [PyStokes Gallery](https://github.com/rajeshrinet/pystokes/wiki/Gallery) for a selected list of research application of PyStokes.



## Phoresis and Stokesian hydrodynamics

PyStokes is a Python library for studying phoretic and hydrodynamic interactions between spherical particles when these interactions can be described by the solutions of, respectively, the Laplace and Stokes equations. The library has been specifically designed for studying these interactions in suspensions of active particles, which are distinguished by their ability to produce flow, and thus motion, in the absence of external forces or torques. Such particles are endowed with a mechanism to produce hydrodynamic flow in a thin interfacial layer, which may be due to the motion of cilia, as in microorganisms (Brennen & Winet, 1977) or osmotic flows of various kinds in response to spontaneously generated gradients of phoretic fields (Ebbens & Howse, 2010). The latter, often called autophoresis, is a generalisation of wellknown phoretic phenomena including, inter alia, electrophoresis (electric field), diffusiophoresis (chemical field) and thermophoresis (temperature field) that occur in response to externally imposed gradients of phoretic fields (Anderson, 1989).



## Phoresis



* **What is phoretic transport?** Quoting (JL Anderson, Ann Rev Fluid Mech (1989))

> Phoretic transport is defined as the movement of colloidal particles by
a field that interacts with the surface of each particle

Examples of phoresis (phoretic field) are electrophoresis (electric field), diffusiophoresis (concentration of chemical field), thermophoresis (temperature field) etc.  The phoretic motion of the colloid in this case is due to the externally imposed electric field. The phoretic transport results from an electrophoretic slip velocity on the surface of the particles as we now describe.

* **Electrophoresis: how does a neutral object (charged particle with negative counterions) propel under the effect of electric field?** Consider a particle with a charge q in a non-conducting fluid. Once an electric field E is applied, the colloid moves due to the body force qE acting on it. The motion also creates flow v, which decays as v ~1/r, where r is the distance from the colloid. On the other hand, if the fluid is electrolytic, then the flow decays as v ~ 1/r^3. Why?<br/>As described in the figure below, a charged particle in an electrolytic fluid is effectively neutral due to the cloud of counterions around it. Thus, there is no net body force and the colloid should not move! It does due to the fact that the counterions are mobile, which leads to a `slip` velocity of the surface of the particle. Solving the Stokes equation with this slip boundary condition gives v ~ 1/r^3.![Image](https://raw.githubusercontent.com/rajeshrinet/pystokes-examples/master/gallery/figs/electrophoresis.jpg)

* **What is autophoretic transport**?  The autophoretic particles self-propel due to a non-equilibrium process, such as chemical reactions on their surface. The phrase auto-phoresis is used interchangeably with self-phoresis. 

* **What is self-diffusiophoretic motion of active particles** Active particles are distinguished by the fact that the slip velocity does not result from an externally imposed field, but is spontaneously generated  in response to gradients of phoretic fields

![Image](https://raw.githubusercontent.com/rajeshrinet/pystokes-examples/master/gallery/figs/self-diffusiophoresis.jpg)
