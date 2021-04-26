Example 1: Cooling dike
=======================

Problem description
-------------------

As a first example, we will use the finite differences form :eq:`eq:FD_heat_flow_explicit_solution` of the heat diffusion equation :eq:`eq:1D_heat_flow` to explore the cooling of a dike. We will look at a  :math:`2m` wide dike that intruded with a temperature of  :math:`1200°C` into  :math:`300°C` warm country rock. The initial conditions can then be written like this:

.. math::
    :label: eq:dike_ic

    \begin{align}
    \begin{split}
    T(x<-1 \mid x>1) = 300\\
    T(x>-1 \& x<1) = 1200
    \end{split}
    \end{align}

In addition, we assume that the temperature far away from the dike center (at :math:`\lvert L/2 \rvert`) remains at a constant temperature. The boundary conditions are thus:

.. math::
    :label: eq:dike_bc

    \begin{align}
    \begin{split}
    T(x=-\frac{L}{2}) = 300\\
    T(x=\frac{L}{2}) = 300
    \end{split}
    \end{align}

.. figure:: /_figures/1D_dike_example_stencil.*
   :align: center
   :name: dike_setup
   :figwidth: 100%
   
   Setup of the model considered here (A). A hot basaltic dike intrudes into colder country rock. Only variations in x-direction are considered; properties in the other directions are assumed to be constant. The initial temperature distribution :math:`T(x,0)` has a step-like perturbation. B) Finite difference discretization of the 1D heat equation. The finite difference method approximates the temperature at given grid points, with spacing :math:`\Delta x`. The time-evolution is also computed at given times with timestep :math:`\Delta t`.


:numref:`dike_setup` summarizes the setup of the cooling dike problem.

FDM notebook
------------

.. toctree::
    :maxdepth: 2

    jupyter/cooling_dike_1d_fdm.ipynb


Excercises
----------

The first step is to get the jupyter notebook to work. A good starting point is to make your own copy of the notebook. One way is to start-up your course python environment in a shell, start a jupyter lab or notebook, and create a new notebook. You can then copy the code blocks (cells) from the script into your notebook and complete the missing pieces. 

.. code-block:: bash
 
    conda activate py37_fem_lecture
    cd $your_working_directory$
    jupyter lab
    
Now you can explore the numerical solution.    

    * Complete the notebook and get it to work
    * Vary the parameters (e.g. use more gridpoints, a larger timestep). Notice how the numerical solution becomes unstable when the timestep is increased beyond a certain value. What does this value depend on? This is a major drawback of explicit finite difference codes such as the one presented here.
    * Record and plot the temperature evolution versus time at a distance of 5 meter from the dike/country rock contact. What is the maximum temperature the country rock experiences at this location and when is it reached?
    * Bonus question: Derive a finite-difference approximation for variable k and variable dx.


FDM notebook - T(t) at x meters from dike
------------------------------------------

.. toctree::
    :maxdepth: 2

    jupyter/cooling_dike_1d_fdm_T_at_x_m.ipynb
