MWR solutions of the steady advection-diffusion equation
========================================================

Let us explore the steady-state advection diffusion equation and find approximate solutions to it using the MWR. Find :math:`u(x)` that satisfies


.. math::
    :label: eq:ss_heat_mwr

    \begin{align}
    \begin{split}
    \frac{\partial u}{\partial t} = K \frac{\partial^2 T}{\partial x^2} - c\frac{\partial u}{\partial x}=0\\
    c\frac{\partial u}{\partial x} - K \frac{\partial^2 T}{\partial x^2} =0\\
    u(0)=0\\
    u(1)=1
    \end{split}
    \end{align}

A scalar quantity :math:`u` (e.g. temperature) is advected with velocity :math:`c` and undergoes a diffusion dependent on K. The exact solution of this problem is:

.. math::
    :label: eq:ss_heat_mwr_analytical

    u^{ex}(x)=C_1 \exp \left ( \frac{c}{K} x \right ) + C_2

By using the boundary conditions, we can find the integration constants.

.. math::
    :label: eq:ss_heat_mwr_analytical_2

    \begin{align}
    \begin{split}
    C_1 &= -C_2\\
    C_2 &= \frac{1}{1-\exp \left (\frac{c}{K} \right )}
    \end{split}
    \end{align}   

.. tip::

    Try to derive the constants yourself and check that :math:`u^{ex}` really is the exact solution!


Let’s solve by the Method of Weighted Residuals using a polynomial function as a basis. That is, let the approximating function be

.. math::
    :label: eq:mwr_example_1
    
    \tilde{u}(x)=a_0 + a_1 x + a_2 x^2

Application of the boundary conditions reveals

.. math::
    :label: eq:mwr_example_2

    \begin{align}
    \begin{split}
    \tilde{u}(0) &= 0=a_1\\
    \tilde{u}(1) &= 1=0+a_1 + a_2\\
    a_1 &= 1-a_2
    \end{split}
    \end{align}   

so that the approximating polynomial which also satisfies the boundary conditions is

.. math::
    :label: eq:mwr_example_3
    
    \tilde{u}(x)=(1-a_2)x+a_2x^2

To find the residual :math:`R(x)`, we need the first and second derivative of this function, which are

.. math::
    :label: eq:mwr_example_4

    \begin{align}
    \begin{split}
    \frac{\partial \tilde{u}}{\partial x} &= (1-a_2)+2a_2 x\\
    \frac{\partial^2 \tilde{u}}{\partial x^2} &= 2a_2
    \end{split}
    \end{align}   

So the residual is

.. math::
    :label: eq:mwr_example_5

    R(x)=c\left ( (1-a_2)+ 2a_2 x \right ) - K(2a_2)

Before continuing let us explore an example problem:
	  
.. math::
    :label: eq:mwr_example_tmp1

    \begin{align}
    \begin{split}
    c &= 1\\
    K &= 1
    \end{split}
    \end{align}   

Collocation method
------------------

For the collocation method, the residual is forced to zero at a number of discrete points. Since there is only one unknown :math:`a_2`, only one collocation point is needed. We choose (arbitrarily, but from symmetry considerations) the collocation point :math:`x = 0.5`. Thus, the equation needed to evaluate the unknown :math:`a_2` is

.. math::
    :label: eq:mwr_collo_ex_1

    R\left(\frac{1}{2}\right)=c\left ( (1-a_2)+ 2a_2 \left(\frac{1}{2}\right) \right ) - K(2a_2)=0

and the contant :math:`a_2` is

.. math::
    :label: eq:mwr_collo_ex_2

    a_2=\frac{c}{2K}

Least-squares method
--------------------

The weight function :math:`W_i`is just the derivative of :math:`R(x)` with respect to the unknown :math:`a_2`:

.. math::
    :label: eq:mwr_lsq_ex_1

    W_1(x)=\frac{dR}{da_2}=-c+2xc-2K

So the weighted residual statement becomes

.. math::
    :label: eq:mwr_lsq_ex_2

    \int_{0}^{1}W_{1}(x)\cdot R(x)dx &= 0 \\
    W_1(x)&=\frac{dR}{da_2}=c(2x-1)-2K \\
    \int_{0}^{1}\left(c(2x-1)-2K \right )\left[c \left( (1-a_2)+2a_2 x \right ) -K(2a_2)\right ] dx &= 0

the math is considerably more involved than before, but nothing more than integration of polynomial terms. Direct evaluation leads to the algebraic relation

.. math::
    :label: eq:mwr_lsq_ex_3

    a_2 = \frac{6cK}{\left (c^2 + 12K^2 \right )}

Galerkin method
---------------

In the Galerkin Method, the weight function :math:`W_1` is the derivative of the approximating function :math:`\tilde{u}` with respect to the unknown coefficient :math:`a_2`:
