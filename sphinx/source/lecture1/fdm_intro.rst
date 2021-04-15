Introduction Finite Differences
===============================

Before progressing towards finite element modeling, we will learn about the Finite Difference Method (FDM, which is a somewhat easier method to solve partial differential equations. We will do so by looking at how heat conduction "works".

Background on conductive heat transport
---------------------------------------

Temperature is one of the key parameters that control and affect the geological processes that we are exploring during this class. Mantle melting, rock rheology, and metamorphism are all functions of temperature. It is therefore essential that we understand how rocks exchange heat and change their temperature.

One fundamental equation in the analysis of heat transfer is Fourier’s law for heat flow:

.. math::
    :label: eq:heat_flow
    
    \vec{q} = -k \nabla T = -k
    \begin{bmatrix}
    \frac{\partial T}{\partial x} \\
    \frac{\partial T}{\partial y} \\
    \frac{\partial T}{\partial z}
    \end{bmatrix}

It states that heat flow is directed from high to low temperatures (that’s where the minus sign comes from) and is proportional to the geothermal gradient. The proportionality constant, k, is the thermal conductivity which has units of W/m/K and describes how well a rock transfers heat. k is a typically a complex function of rock type, porosity, and temperature yet is often simplified to a constant value.

In most cases, we are interested in temperature and not heat flow so that we would like to have an equation that describes the temperature evolution in a rock. We can get this by deriving a conservation law for heat.

.. figure:: /_figures/3D_cube_heat_flow.*
   :align: center
   :name: heat_flux_box
   :figwidth: 100%

   Derivation of the energy equation. Change in internal energy is related to changes in heat fluxes into and out of the box (and a lot of other terms (e.g. advection) that are neglected here).

For simple incompressible cases, changes in internal energy can be expressed as changes in temperature timesdue density and specific heat. The change in internal energy with time can now be written as :math:`\rho c_p \frac{\partial}{\partial t} T \Delta x \Delta y \Delta z` (:math:`\rho` is density, :math:`c_p` specific heat, and :math:`T` temperature), has units of (:math:`J/s`), and must be equal to the difference between the heat flow into the box :math:`q_{in} \Delta y \Delta z \left( \frac{W}{m K}\frac{K}{m}m^2 = \frac{J}{s}\right)` and the heat flow out of the box :math:`\left( q_{in} + \frac{\partial q_{in}}{\partial x} \Delta x\right) \Delta y \Delta z` (the y and z directions are done in the same way). With these considerations, we can write a conservation equation for energy:

.. math::
    :label: eq:1D_heat_flow

    \begin{align}
    \begin{split}
    \rho c_p \frac{\partial T}{\partial t} = -\frac{\partial q_{in}}{\partial x} = \frac{\partial}{\partial x}k\frac{\partial T}{\partial x}\\
    \frac{\partial T}{\partial t} = \frac{k}{\rho c_p} \frac{\partial^2 T}{\partial x^2} = \kappa \frac{\partial^2 T}{\partial x^2}
    \end{split}
    \end{align}

:eq:`eq:1D_heat_flow` is called the heat transfer or heat diffusion equation and is one of the most fundamental equations in Earth Sciences. If the thermal conductivity is constant, we can define a thermal diffusivity, :math:`\kappa`, and write the simpler second form of the equation.

Note how the changes in heat flow can be written in terms of a divergence:

.. math::
    :label: eq:1D_heat_flow_div

    \rho c_p \frac{\partial T}{\partial t} =
    -\nabla \cdot (\vec{q}) =
    \frac{\partial}{\partial x}k\frac{\partial T}{\partial x} + \frac{\partial}{\partial y}k\frac{\partial T}{\partial y} + \frac{\partial}{\partial z}k\frac{\partial T}{\partial z}

Or in vector notation:

.. math::
    :label: eq:1D_heat_flow_vec

    \rho c_p \frac{\partial T}{\partial t} =
    -\nabla \cdot \left(\vec{q}\right) =
    \begin{bmatrix}
    \frac{\partial}{\partial x} & \frac{\partial}{\partial y} & \frac{\partial}{\partial z}
    \end{bmatrix}
    k
    \begin{bmatrix}
    \frac{\partial T}{\partial x}\\
    \frac{\partial T}{\partial y}\\
    \frac{\partial T}{\partial z}
    \end{bmatrix}

The Finite Differences Method
-----------------------------

:eq:`eq:1D_heat_flow` is a partial differential equation that describes the evolution of temperature. There are two fundamentally different ways to solve it: 1) analytically or 2) numerically. Analytical solutions have the virtue that they are exact but it is often not possible to find one for complex systems. Numerical solutions are always approximations but can be found also for very complex systems. We will first use one numerical technique called finite differences. To learn how partial differential equations are solved using finite differences, we have to go back to the definition of a derivative:

.. math::
    :label: eq:FD_heat_flow_d1

    \frac{\partial T}{\partial x} =
    \lim_{\Delta x \rightarrow 0} \frac{T(x + \Delta x) - T(x)}{\Delta x}

In our case, the above derivative describes the change in temperature with x (our space coordinate). In numerics, we always deal with discrete functions (as opposed to continuous functions), which only exist at grid points (Fig. \ref{fig:dike_setup}). We can therefore translate the above equation into computer readable form:

.. math::
    :label: eq:FD_heat_flow_d2

    \frac{\partial T}{\partial x} \approx
    \frac{T(i + 1) - T(i)}{x(i+1) - x(i)} =
    \frac{T(i + 1) - T(i)}{\Delta x}	

:math:`T(i)` is the temperature at a grid point :math:`i` and :math:`\Delta x` is the grid spacing between two grid points. Using this definition of a derivative, we can now compute the heat flow from the calculated temperature solution:

.. math::
    :label: eq:FD_heat_flow_d3

    q_x = -k \frac{\partial T}{\partial x} \approx -k \left(\frac{T(i+1)-T(i)}{x(i+1)-x(i)}\right)

This form is called the *finite differences* form and is a first step towards solving partial differential equations numerically.

Note that it actually matters in which direction you count: usually it makes life much easier if indices and physical coordinates point in the same direction, e.g. x coordinate and index increase to the right.

We have learned how we can compute derivatives numerically. The next step is to solve the heat conduction equation :eq:`eq:1D_heat_flow` completely numerically. We are interested in the temperature evolution versus time :math:`T(x, t)` which satisfies :eq:`eq:1D_heat_flow`, given an initial temperature distribution. We know already from the heat flow example how to calculate first derivatives (forward differencing):

.. math::
    :label: eq:FD_temperature

    \frac{\partial T}{\partial t} = \frac{T_i^{n+1} - T_i^n}{\Delta t}

The index $n$ corresponds to the time step and the index $i$ to the grid point (x-coordinate). Next, we need to know how to write second derivatives. A second derivative is just a derivative of a derivate. So we can write (central differencing):

.. math::
    :label: eq:FD_heat_flow_central_difference

    \kappa\frac{\partial^2 T}{\partial x^2} \approx \kappa \frac{\frac{T_{i+1}^n - T_i^n}{\Delta x}-\frac{T_{i}^n - T_{i-1}^n}{\Delta x}}{\Delta x} = \kappa \frac{T_{i+1}^n-2T_{i}^n+T_{i-1}^n}{\Delta x^2}

If we combine equation :eq:`eq:FD_temperature` and :eq:`eq:FD_heat_flow_central_difference` we get:

.. math::
    :label: eq:FD_heat_flow_explicit

    \frac{T_i^{n+1}-T_i^n}{\Delta t} = \kappa \left(\frac{T_{i+1}^n - 2T_i^n + T_{i-1}^n}{{\Delta x^2}}\right)

The last step is a rearrangement of the discretized equation, so that all known quantities (i.e. temperature at time :math:`n`) are on the right-hand side and the unknown quantities on the left-hand side (properties at :math:`n+1`). This results in:

.. math::
    :label: eq:FD_heat_flow_explicit_solution
    
    T_i^{n+1} = \frac{\kappa \Delta t}{\Delta x^2} \left( T_{i+1}^n - 2 T_i^n + T_{i-1}^n\right) + T_i^n

We have now translated the heat conduction equation :eq:`eq:1D_heat_flow` into a computer readable finite differences form.



