# Ristorcelli-Initial-Conditions-Isotropic-Turbulence
Initial Conditions by Ristorcelli and Blaisdell (1996) for Isotropic Turbulence.

It is important to initialize compressible Isotropic Turbulence with the appropriate initial conditions. The method developed by
Ristorcelli and Blaisdell (1996) - https://ntrs.nasa.gov/archive/nasa/casi.ntrs.nasa.gov/19960045437.pdf , initializes
the flow-field with appropriate "weakly compressible" thermodynamic and velocity fluctuations. This method requires the presence
of solenoidal (incompressible) velocity fluctuations as a starting point. There are multiple methods to develop incompressible
fluctuations, the most prominent amongst them the work by Rogallo (1981) - https://ntrs.nasa.gov/archive/nasa/casi.ntrs.nasa.gov/19810022965.pdf .

A  perturbation expansion of the thermodynamic and velocity fields are performed about their mean values, with the perturbation
parameter being the turbulent mach number (defined as the ratio root mean square velocity fluctation and mean speed of sound).

The Navier-Stokes equations are reduced to first order perturbation equations based on this expansion. A second order 
expression is used for the dilatation from which the isentropic thermodynamic fluctuations are derived.

A couple of Poisson equations (one for the pressure and the other for velocity temporal derivative) have to be solved 
to obtain the flow-field. The solver file also contains the Poisson solver. Spatial Derivatives are computed in spectral space
for high accuracy.
