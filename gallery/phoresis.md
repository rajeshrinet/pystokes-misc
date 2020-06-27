## Phoresis and Stokesian hydrodynamics


The PyStokes library allows to study the transport of active colloidal particles in a [Stokesian fluid](https://en.wikipedia.org/wiki/Stokes_flow). Colloidal particles are of size between few nanometers to several microns, and thus, the fluid flow is given by the Stokes equation. Active particles are distinguished by their ability to produce flow, and thus motion, in the absence of external forces or torques on them. Examples of active particles include auto-phoretic synthetic colloids [1], and microorganisms [2]. In what follows we use colloidal and particles interchangeably. 

The auto-phoretic particles self-propel due to a non-equilibrium process, such as chemical reactions on their surface. In this case, the phoretic field is the chemical concentration. In general, to quote JL Anderson [3],

> Phoretic transport is defined as the movement of colloidal particles by
a field that interacts with the surface of each particle

Examples of phoresis (phoretic field) are electrophoresis (electric field), diffusiophoresis (concentration of chemical field), thermophoresis (temperature field) etc. We now describe electrophoresis in detail. 

## Electrophoretic motion of colloidal particles

In this section, we describe the classical phoretic motion of a charged colloidal particle in a conducting fluid. The phoretic motion of the colloid in this case is due to the externally imposed electric field. The phoretic transport results from an electrophoretic slip velocity on the surface of the particles as we now describe. 

Consider a particle with a charge q in a non-conducting fluid. Once an electric field E is applied, the colloid moves due to the body force qE acting on it. The motion also creates flow v, which decays as v ~1/r, where r is the distance from the colloid. On the other hand, if the fluid is electrolytic, then the flow decays as v ~ 1/r^3. 

As described in the figure below, a charged particle in an electrolytic fluid is effectively neutral due to the cloud of counterions around it. Thus, their is no net body force and the colloid should not move! It does due to the fact that the counterions are mobile, which leads to a `slip` velocity of the surface of the particle. Solving the Stokes equation with this slip boundary condition gives v ~ 1/r^3, (see [1] for derivation)

![Image](https://raw.githubusercontent.com/rajeshrinet/pystokes-examples/master/gallery/figs/electrophoresis.jpg)

## Self-diffusiophoretic motion of active particles

Active particles are distinguished by the fact that the slip velocity does not result from an externally imposed field, but due to a non-equilibrium process on the surface of the particle. 

![Image](https://raw.githubusercontent.com/rajeshrinet/pystokes-examples/master/gallery/figs/self-diffusiophoresis.jpg)


[1] S. J. Ebbens and J. R. Howse, “In pursuit of propulsion at the nanoscale,” Soft Matter 6, 726–738 (2010).

[2] C. Brennen and H. Winet, “Fluid mechanics of propulsion by cilia and flagella,” Annu. Rev. Fluid Mech. 9, 339–398 (1977).

[3] JL Anderson, Ann Rev Fluid Mech, 21: 61-99 (1989)
