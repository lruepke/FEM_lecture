Excercise: Transient case
---------------------------
The next step is to derive and implement the unsteady (time-dependent) heat diffusion equation and to solve an example problem for the cooling of the lithosphere. The unsteady heat diffusion equation looks like this:

.. math::
    :label: eq:fem_2d_transient

    \rho c_p \frac{\partial T}{\partial t} = \frac{\partial}{\partial x} k \frac{\partial T}{\partial x} + \frac{\partial}{\partial y}k\frac{\partial T}{\partial y}
    
Two differences to the previously considered steady heat diffusion are apparent. First, we now have a time derivative in addition to the spatial derivatives. Second, the material parameters (density, specific heat, thermal conductivity) are not constant (neglected) any more. Nevertheless, the diffusion part of :eq:`eq:fem_2d_transient` looks very similar to the steady-state case that we know how to solve.

We know how to handle spatial derivatives but this is the first time we encounter a time derivative in finite elements. To solve it we will use a “trick” from finite differences – we will simply write the time derivative in finite differences form.

.. math::
    :label: eq:fem_2d_transient_2

    \rho c_p \frac{T^{n+1} - T^n}{\Delta t} = \frac{\partial}{\partial x} k\frac{\partial T^{n+1}}{\partial x} + \frac{\partial}{\partial y}k\frac{\partial T^{n+1}}{\partial y}.

Re-arrange :eq:`eq:fem_2d_transient` so that all known temperatures :math:`T^n` are on the Rhs and all unknown temperatures :math:`T^{n+1}` are on the Lhs.

    #. Derive the finite element form of this equation.
    #. Implement the equation by modifying the steady-sate script
        * Note, you will get a non-zero Rhs. Consequently you will obtain a force vector for each element (containing the old temperatures) that needs to be added to the global force vector.
    #. As an example problem we will considered a simplified lithosphere cooling case. Use these model parameters:
        * Ttop=0, Tbot=1300, density=3000kg/m3, cp=1000J/Kg/K and k=3W/m/K
        * Make the box 100x100km
        * Initialize all temperatures to Ttop
        * Take a time step of 1Ma.