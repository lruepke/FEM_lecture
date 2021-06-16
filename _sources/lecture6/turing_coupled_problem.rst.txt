Coupled diffusion problems
==========================
So far we have looked at problems that included a single unknown (e.g. temperature) at each node. Sometimes it is necessary to solve for multiple unknowns per node simultaneously. This could just be multiple displacement or velocity components (as in mechanics), or this could be necessary when equations strongly coupled.

One example, where a system of coupled equations needs to be solved are so-called Turing pattern :cite:`Turing1952`, which appear to emerge in many biological and chemical systems quasi spontaneoulsly from near homogeneous starting conditions. A modern look at it can be found in the paper by :cite:`Maini2012`. Here we follow the implementation and problem description given in the FEM book by Guy Simpson :cite:`Simpson2017`. 

Governing equations
-------------------
We will explore the evolution of two chemical compunds, :math:`A` and :math:`B`, which affect each other. One is an 'activator' and one is an 'inhibitor'. Both compounds are also produced by a background rate and the their concentrations control each other. This can be expressed by a set of coupled equations.

.. math::
    :label: eq:turing_eqs
    
    \begin{align}
    \begin{split}
    \frac{\partial A}{\partial t} &= \nabla A + \gamma (a - A +A^2 B)\\
    \frac{\partial B}{\partial t} &= d\nabla B + \gamma (b - A^2 B)
    \end{split}
    \end{align}

We can see that the equations are coupled through the source terms that contain the concentrations of both species. 

FEM discretization
------------------
Using the same FD time discretization and using an implicit scheme, we can write these equations in FEM matrix form (just like we did for temperature in the previous example):

.. math::
    :label: eq:turing_eqs_fem

    \begin{align}
    \begin{split}
    (M(1+\Delta t\gamma) + \Delta tA_A)A^{n+1} = MA^n + \Delta tF_A^{n+1}\\
    (M + \Delta tA_B)B^{n+1} = MB^n + \Delta tF_B^{n+1}
    \end{split}
    \end{align}

with the matrices defined as:

.. math::
    :label: eq:turing_eqs_fem_matrix

    \begin{align}
    \begin{split}
    M &= \int_\Omega N_i N_j  d\Omega \ \ \ \ \ \ \ i,j=1,2,...,n\\
    A_A &= \int_\Omega \left ( \frac{\partial N_i}{\partial x}\frac{\partial N_j }{\partial x} + \frac{\partial N_i}{\partial y}\frac{\partial N_j }{\partial y} d\Omega \right ) d\Omega\ \ \ \ \ \ \ i,j=1,2,...,n\\
    A_B &= \int_\Omega d \left ( \frac{\partial N_i}{\partial x}\frac{\partial N_j }{\partial x} + \frac{\partial N_i}{\partial y}\frac{\partial N_j }{\partial y} d\Omega \right ) d\Omega\ \ \ \ \ \ \ i,j=1,2,...,n\\    
    F_A &= \int_\Omega \gamma N_i \left (a + (N_jA_j)^2N_jB_j      \right )  d\Omega \ \ \ \ \ \ \ i,j=1,2,...,n\\ 
    F_B &= \int_\Omega \gamma N_i \left (b - (N_jA_j)^2N_jB_j      \right )  d\Omega \ \ \ \ \ \ \ i,j=1,2,...,n\\ 
    \end{split}
    \end{align}

We can combine the matrices on the left-hand side and get the stifness matrices :math:`K_A` and :math:`K_B`.

.. math::
    :label: eq:turing_eqs_fem_matrix_2
    
    \begin{align}
    \begin{split}
    K_A A^{n+1} = M A^n + \Delta t F_A^{n+1}\\
    K_B B^{n+1} = M B^n + \Delta t F_B^{n+1}
    \end{split}
    \end{align}

One way of solving such a coupled problem is to have two unknowns per node (two so-called degrees of freedom). In our case, the complete element stiffness matrix would look like this:

.. math::
    :label: eq:turing_eqs_fem_matrix_3
    
    \begin{bmatrix}
    {K_A}_{11} & 0 & {K_A}_{12} & 0 & {K_A}_{13} & 0 \\
    0 & {K_B}_{11} & 0 & {K_B}_{12} & 0 & {K_B}_{13} \\
    {K_A}_{21} & 0 & {K_A}_{22} & 0 & {K_A}_{23} & 0 \\
    0 & {K_B}_{21} & 0 & {K_B}_{22} & 0 & {K_B}_{23} \\
    {K_A}_{31} & 0 & {K_A}_{32} & 0 & {K_A}_{33} & 0 \\
    0 & {K_B}_{31} & 0 & {K_B}_{32} & 0 & {K_B}_{33}
    \end{bmatrix}
    \begin{bmatrix}
    A_1^{n+1}\\
    B_1^{n+1}\\
    A_2^{n+1}\\
    B_2^{n+1}\\
    A_3^{n+1}\\
    B_3^{n+1}\\
    \end{bmatrix}
    = Rhs

The unknown concentrations of :math:`A` and :math:`B` or both showing up in the solution vector. What we also note is that the coupling is not that strong, as there are no cross-terms. The equation for :math:`A_1^{n+1}`, the first row in :eq:`eq:turing_eqs_fem_matrix_3` has zeros in the columns that operate on B. The coupling comes in through the source terms that appear on the Rhs (:eq:`eq:turing_eqs`).

