2-D FEM - Heat diffusion
=========================================

Derivation of the weak form
---------------------------

Let’s move on to 2D! The strong form looks like this, we have just added the second dimension.

.. math::
    :label: eq:fem_2d_strong

    \frac{\partial}{\partial x}k\frac{\partial T_{ex}}{\partial x} + \frac{\partial}{\partial y}k\frac{\partial T_{ex}}{\partial y} =0. 


We now use 2-D shape functions in our approximate solution:

.. math::
    :label: eq:fem_aprox_funcion_2d

    T_{ex} \cong \tilde{T}= \sum_{j=1}^{n} N_j(x,y) T_j = N_j T_j, 

where the shape functions :math:`N_j(x,y)` are now function of both spatial dimensions. An example of a 2-D FEM mesh and associated shape functions are shown in :numref:`fig:shapeFunc:2D:linear`. Note that the structured quad mesh and the bi-linear shape functions are just one possibility of using the FEM in 2-D. There are many other possibility, like unstructured triangle meshes and higher order shape functions.

.. figure:: Schematic_FEM/shapeFunction_2D_Q1.svg
    :name: fig:shapeFunc:2D:linear
    :align: center

    Example of 2D linear finite element shape functions

Note that the 2-D shape functions still sum to :math:`1` everywhere:

.. math::
    :label: eq:fem_aprox_funcion_2d_2

    \sum_{j=1}^{n} N_j(x,y)  = 1. 



We proceed by substituting :eq:`eq:fem_aprox_funcion_2d` into :eq:`eq:fem_2d_strong` and by using the Galerkin method. This results in:

.. math::
    :label: eq:fem_2d_weak
    
    \int_\Omega N_i \frac{\partial}{\partial x}k\frac{\partial N_j T_j }{\partial x} + N_i \frac{\partial}{\partial y}k\frac{\partial N_j T_j }{\partial y} d\Omega=0\ \ \ \ \ \ \ i=1,2,...,n


We do the integration by parts and obtain:

.. math::
    :label: eq:fem_2d_weak_2

    \begin{align}
    \begin{split}
    -\int_\Omega \left ( \frac{\partial N_i}{\partial x}k\frac{\partial N_j T_j }{\partial x} + \frac{\partial N_i}{\partial y}k\frac{\partial N_j T_j }{\partial y} d\Omega \right ) &+ \oint_{\Gamma} \left ( N_ik\frac{\partial N_j T_j }{\partial x}  +  N_ik\frac{\partial N_j T_j }{\partial y} \right ) d\Gamma    =0\ \ \ \ \ \ \ i=1,2,...,n\\
    \Rightarrow \\
    +\int_\Omega \left ( \frac{\partial N_i}{\partial x}k\frac{\partial N_j}{\partial x} + \frac{\partial N_i}{\partial y}k\frac{\partial N_j}{\partial y} T_j d\Omega \right ) &+ \oint_{\Gamma} N_i \vec{q}\vec{n}  d\Gamma    =0\ \ \ \ \ \ \ i=1,2,...,n\\
    \end{split}
    \end{align}


Note how the signs flipped during integration by parts and by substituting :math:`\vec{q}=-k\nabla T`. The line integral has a clear physical meaning now, it's the heat flow in and out of the modeling domain.

We go on and split the modeling domain into finite elements and the first integral turns into a sum over those elements. The line integral is again obmitted for the moment.

.. math::
    :label: eq:fem_2d_weak_3

    \begin{align}
    \begin{split}
    +\int_\Omega \left ( \frac{\partial N_i}{\partial x}k\frac{\partial N_j}{\partial x} + \frac{\partial N_i}{\partial y}k\frac{\partial N_j}{\partial y} T_j \right ) d\Omega   &=0\ \ \ \ \ \ \ i=1,2,...,n\\
    \Rightarrow \\
    +\int_\Omega \left ( \frac{\partial N_i}{\partial x}k\frac{\partial N_j}{\partial x} + \frac{\partial N_i}{\partial y}k\frac{\partial N_j}{\partial y} T_j \right ) d\Omega   &= \sum_{Elements} \int_{\Omega_{e}} \left ( \frac{\partial N_i}{\partial x}k\frac{\partial N_j}{\partial x} + \frac{\partial N_i}{\partial y}k\frac{\partial N_j}{\partial y}  \right ) T_j d\Omega_{e} =0\ \ \ \ \ \ \ i=1,2,...,n\\
    \end{split}
    \end{align}

2-D FEM: connectivity
----------------------
The basic concept of the finite element method is to solve/assemble the system of equations, for e.g. heat diffusion, on each element and add all element contributions together to obtain the global matrix equation.
Every element has an element number and a certain number of nodes. We will initially use quadratic elements with four nodes. For the connectivity between elements we will need two matrices: a CGCOORD matrix that has the size [nnod,2], where 2 is the number of dimensions (x,y) and nnod is the total number of nodes in the mesh, and ELEM2NODE, which has the size [nel, nnodel] (nel is the total number of elements and nnodel is the number of nodes per element, i.e. 4), see :numref:`fig:mesh:2D:structured`. We had already used this matrices in 1-D but now their meaning becomes clear.

.. tab:: Structured mesh 1

    .. figure:: Schematic_FEM/mesh2D_structured.*
        :name: fig:mesh:2D:structured
        :align: center
    
        Structured mesh example of 2D finite element

.. tab:: Structured mesh 2

    .. figure:: Schematic_FEM/mesh2D_structured_usi.*
        :name: fig:mesh:2D:structured:usi
        :align: center
    
        :download:`Structured mesh <Schematic_FEM/mesh/structure.msh>` example of 2D finite element

.. tab:: Unstructured mesh 

    .. figure:: Schematic_FEM/mesh2D_unstructured.*
        :name: fig:mesh:2D:unstructured
        :align: center
    
        :download:`Unstructured mesh <Schematic_FEM/mesh/unstructure.msh>` example of 2D finite element

We will use a connectivity as shown in :numref:`fig:matrix:2D:structured`. Element 0 has, for example, the global nodes 0, 1, 6, 5. Note that the local element numbering is always counterclockwise (in this class). It therefore contributes to the calculations of those four temperatures. The element stiffness matrix will now be [4,4]. The contribution of all elements is added/assembled into the glob

.. tab:: Structured mesh 1

    .. figure:: Schematic_FEM/Matrix2D_structured.svg
        :name: fig:matrix:2D:structured
        :align: center

        Connectivity example of 2D finite element with structured mesh

.. tab:: Structured mesh 2

    .. figure:: Schematic_FEM/Matrix2D_structured_usi.svg
        :name: fig:matrix:2D:structured:usi
        :align: center

        Connectivity example of 2D finite element with structured mesh
 
.. tab:: Unstructured mesh

    .. figure:: Schematic_FEM/Matrix2D_unstructured.svg
        :name: fig:matrix:2D:unstructured
        :align: center

        Connectivity example of 2D finite element with unstructured mesh

Excercise
^^^^^^^^^^

#. Create the global coordinate vector GCOORD. Assume that the length of the box in x direction (lx) is 1 ([0,1]) and the length of box is y direction (ly) is also 1 ([0, 1]). The number of nodes in x is 5 (nx) and the number of nodes in y is 4 (ny).  
#. Create the connectivity matrix ELEM2NODE. The element numbering should look like this:

.. tip::

    * Hint loop over all nodes and create their respective x and y coordinates. The node numbering should be like in :numref:`fig:mesh:2D:structured` (see jupyter notebook for a discussion).
    * Functions like ceil or floor might help you with the connectivity matrix; it might also help to first compute the row we are in and then compute the index of the lower left node and move on from there.

.. tip::

    We use a function called tabulate to format the output in the notebooks. You will probably have to install it into your virtual python environment.

    .. code-block:: bash

        conda activate py37_fem_class
        conda install tabulate


.. toctree::
    :maxdepth: 2

    jupyter/2d_fem_connectivity.ipynb


Numerical integration
----------------------
In the previous session we have learned how to set up a finite element mesh and the connectivity matrix that relates elements and nodes. Next we will learn about numerical integration and how shape functions are used to interpolate on elements.

We have seen that the general finite element form, called weak form, of the heat diffusion equation can be written as (with boundary terms omitted):

.. math::
    :label: eq:fem_2d_weak_num_int

    \int_\Omega \left ( \frac{\partial N_i}{\partial x}k\frac{\partial N_j}{\partial x} + \frac{\partial N_i}{\partial y}k\frac{\partial N_j}{\partial y} T_j \right ) d\Omega   = \sum_{Elements} \int_{\Omega_{e}} \left ( \frac{\partial N_i}{\partial x}k\frac{\partial N_j}{\partial x} + \frac{\partial N_i}{\partial y}k\frac{\partial N_j}{\partial y}  \right ) T_j d\Omega_{e} =0\ \ \ \ \ \ \ i=1,2,...,n



The main problem with solving :eq:`eq:fem_2d_weak_num_int` is that the integral may be very difficult to evaluate. In finite elements we can use elements of very different geometries to mesh complex objects, which makes the 2-D integrals very hard to solve. We will therefore use numerical integration techniques like Gauss-Legendre quadrature. The problem is (or good thing  as we will see later) that basically all quadrature rules (for quads) apply to rectangular domains that are 2x2 in size and have a local coordinate system :math:`(−1 \leq (\xi,\eta) \leq 1)`. The integration of a function over a 2d domain can than be expressed as:

.. math::
    :label: eq:num_int_1

    \int_{\hat{\Omega}} F(\xi,\eta) d\xi d\eta = \int_{-1}^{1}\int_{-1}^{1}F(\xi,\eta)d\xi d\eta = \sum_{I=1}^{M} \sum_{J=1}^{N} F(\xi_I,\eta_J)W_I W_J = \sum_{ip}F(\xi_{ip},\eta_{ip})W_{ip},

where M and N denote the number of Gauss points, math:`(\xi,\eta)`denote the Gauss point coordinates, and :math:`W_I` and :math:`W_J` denote the corresponding Gauss weights. The selection of the number of Gauss points is based on the formula :math:`N=int[(p+1)/2]+1`, where p the polynomial degree to which the integrand is approximated. We will use linear shape functions so that we have 2 Gauss points (in every direction). In this case the local coordinates are :math:`\pm0.5773502692` and the weights are all :math:`1`.

So far so good, but how do we solve the integral in :eq:`eq:fem_2d_weak_num_int` using local coordinates? We will have to do a coordinate transformation for the purpose of numerical integration that relates the element :math:`\Omega_e` to the master element area :math:`\hat{\Omega}` - or, equivalently, :math:`(x,y)` to :math:`(ξ,η)`. This can be done by using interpolation functions in terms of local coordinates. We will use so-called *isoparametric* elements in which we use the same number of interpolation functions for coordinate mapping and primary variable interpolation. This implies for our four node element that the local variables vary bi-linearly over the element.
    
.. math::
    :label: eq:num_int_isoparam

    \begin{align}
    \begin{split}
    x=&\sum_{j=1}^n \Phi_j(\xi,\eta)x_j \\
    y=&\sum_{j=1}^n \Phi_j(\xi,\eta)y_j, \\
    \end{split}
    \end{align}

where :math:`x_j` are the nodal coordinates. Our primary variable (temperature) interpolation functions are still in global coordinates:

.. math::
    :label: eq:num_int_isoparam_pv

    T(x,y) = \sum_{j=1}^N N_j(x,y)T_j

We can express the shape functions N in terms of local coordinates using :eq:`eq:num_int_isoparam` and the chain rule of differentiation:

.. math::
    :label: eq:num_int_isoparam_2

    \begin{align}
    \begin{split}
    \frac{\partial N}{\partial \xi} =& \frac{\partial N}{\partial x} \frac{\partial x}{\partial \xi} + \frac{\partial N}{\partial y} \frac{\partial y}{\partial \xi} \\
    \frac{\partial N}{\partial \eta} =& \frac{\partial N}{\partial x} \frac{\partial x}{\partial \eta} + \frac{\partial N}{\partial y} \frac{\partial y}{\partial \eta} \\
    \end{split}
    \end{align}

Or in matrix form:

.. math::
    :label: eq:num_int_isoparam_3

    \begin{bmatrix}
    \frac{\partial N}{\partial \xi} \\
    \frac{\partial N}{\partial \eta}
    \end{bmatrix} = 
    \begin{bmatrix}
    \frac{\partial x}{\partial \xi} &  \frac{\partial y}{\partial \xi} \\
    \frac{\partial x}{\partial \eta} & \frac{\partial y}{\partial \eta}
    \end{bmatrix}
    \begin{bmatrix}
    \frac{\partial N}{\partial x} \\
    \frac{\partial N}{\partial y},
    \end{bmatrix} 


which gives the relation between the derivatives of :math:`N` with respect to the global and local coordinates. The matrix in :eq:`eq:num_int_isoparam_3` is called the Jacobian matrix, J, of the coordinate transformation. Unfortunately this relationship is in the wrong direction. If we have a look at our starting equation :eq:`eq:fem_2d_weak_num_int`, we realize that we need to relate the derivative with respect to global coordinates to the derivative with respect to local coordinates (which is the opposite direction). Fortunately this can be done quite easily:

.. math::
    :label: eq:num_int_isoparam_4

    \begin{bmatrix}
    \frac{\partial N}{\partial x} \\
    \frac{\partial N}{\partial y}
    \end{bmatrix} = 
    \begin{bmatrix}
    J^{-1}
    \end{bmatrix}
    \begin{bmatrix}
    \frac{\partial N}{\partial \xi} \\
    \frac{\partial N}{\partial \eta},
    \end{bmatrix} 


Using :eq:`eq:num_int_isoparam` we can spell out the different components oft the Jacobian matrix:


.. math::
    :label: eq:num_int_jacobian

    \begin{bmatrix}
    J
    \end{bmatrix} = 
    \begin{bmatrix}
    \frac{\partial x}{\partial \xi} &  \frac{\partial y}{\partial \xi} \\
    \frac{\partial x}{\partial \eta} & \frac{\partial y}{\partial \eta}
    \end{bmatrix} = 
    \begin{bmatrix}
    \sum_{j=1}^n x_j \frac{\partial \Phi_j}{\partial \xi} & \sum_{j=1}^n y_j \frac{\partial \Phi_j}{\partial \xi}\\
    \sum_{j=1}^n x_j \frac{\partial \Phi_j}{\partial \eta} & \sum_{j=1}^n y_j \frac{\partial \Phi_j}{\partial \eta},\\
    \end{bmatrix} 

Finally the element area dA=dxdy of the original element in global coordinates is transformed to:

.. math::
    :label: eq:num_int_jacobian_1

    d\Omega_e = dxdy = det(J)d\xi d\eta

We have now all the information we need to solve the integral in :eq:`eq:fem_2d_weak_num_int` numerically:

.. math::
    :label: eq:num_int_jacobian_2

    \int_{\Omega_{e}} \left ( \frac{\partial N_i}{\partial x}k\frac{\partial N_j}{\partial x} + \frac{\partial N_i}{\partial y}k\frac{\partial N_j}{\partial y}  \right ) d\Omega_{e} T_j \approx \\ \left [ \sum_{ip} \left ( \frac{N_i(\xi_{ip},\eta_{ip})}{\partial x} k \frac{N_j(\xi_{ip},\eta_{ip})}{\partial x} + \frac{N_i(\xi_{ip},\eta_{ip})}{\partial y} k \frac{N_j(\xi_{ip},\eta_{ip})}{\partial y}
    \right ) W_{ip} det(J) \right ] T_j

It should be noted that the transformation of a quadrilateral element of a mesh to a master element :math:`\hat{\Omega}` is solely for the purpose of numerically evaluating the integrals in :eq:`eq:fem_2d_weak_num_int` . No transformation of the physical domain or elements is involved in the finite element analysis.

Shape functions
----------------

By now we have learned how to solve integrals of shape functions and transform shape functions between different coordinate systems. It is really time that we learn what these shape functions are! We will use four-node rectangular element with bilinear shape functions. Remember that these shape functions look like this:

.. figure:: Schematic_FEM/shapeFunction_2D_Q1.svg
    :name: fig:shapeFunc:2D:linear_1
    :align: center

    Shape of 2-D bi-linear element shape functions

A shape function has the value 1 at its node, zero at all the others, varies linearly between them and the sum of all shape functions is always 1. A simple way to think about this is in 1D (remember the last script)! With these conventions we can now spell out a general interpolation scheme using shape functions:


.. math::
    :label: eq:shape_functions

    \begin{align}
    \begin{split}
    \tilde{T}(x,y) = \sum_{j=1}^{n} N_jT_J = N_jT_J = NT \\
    NT &= \begin{bmatrix} N_1 & N_2 & N_3 & N_4 \end{bmatrix} \begin{bmatrix} T_1 \\ T_2 \\ T_3 \\ T_4 \end{bmatrix}\\
    N_1 & = 0.25(1-\xi)(1-\eta) \\
    N_2 & = 0.25(1+\xi)(1-\eta) \\
    N_3 & = 0.25(1+\xi)(1+\eta) \\
    N_4 & = 0.25(1-\xi)(1+\eta) \\
    \end{split}
    \end{align}

where :math:`T_{j=1..4}` are the four nodal temperatures, :math:`\tilde{T}` is the temperature at an (integration) point inside the element, :math:`N` are the four shape functions, and :math:`(\xi,\eta)` are the two local coordinates between -1 and 1.