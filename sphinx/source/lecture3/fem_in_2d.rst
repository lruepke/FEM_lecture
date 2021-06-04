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

.. tab:: Four‐Node Bilinear Quadrilateral

    .. figure:: Schematic_FEM/shapeFunction_2D_Q1.*
        :name: fig:shapeFunc:2D:linear
        :align: center

        Example of 2D linear finite element shape functions
 
.. tab:: Nine‐Node Biquadratic Quadrilateral

    .. figure:: Schematic_FEM/shapeFunction_2D_Q2.*
        :name: fig:shapeFunc:2D:quad:Q2
        :align: center

        Example of 2D finite element shape functions

.. tab:: Three‐Node Triangle Element
    :new-set:

    .. figure:: Schematic_FEM/shapeFunction_2D_triangle_Q1.*
        :name: fig:shapeFunc:2D:tri:Q1
        :align: center

        Example of 2D finite element shape functions

.. tab:: Six‐Node Quadratic Triangle

    .. figure:: Schematic_FEM/shapeFunction_2D_triangle_Q2.*
        :name: fig:shapeFunc:2D:tri:Q2
        :align: center

        Example of 2D finite element shape functions

.. tab:: Ten‐Node Qubic Triangle

    .. figure:: Schematic_FEM/shapeFunction_2D_triangle_Q3.*
        :name: fig:shapeFunc:2D:tri:Q3
        :align: center

        Example of 2D finite element shape functions

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

    .. figure:: Schematic_FEM/Matrix2D_structured.*
        :name: fig:matrix:2D:structured
        :align: center

        Connectivity example of 2D finite element with structured mesh

.. tab:: Structured mesh 2

    .. figure:: Schematic_FEM/Matrix2D_structured_usi.*
        :name: fig:matrix:2D:structured:usi
        :align: center

        Connectivity example of 2D finite element with structured mesh
 
.. tab:: Unstructured mesh

    .. figure:: Schematic_FEM/Matrix2D_unstructured.*
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

2-D Shape functions
-------------------

By now we have learned how to solve integrals of shape functions and transform shape functions between different coordinate systems. It is really time that we learn what these shape functions are! We will use four-node rectangular element with bilinear shape functions. Remember that these shape functions look like this:

.. figure:: Schematic_FEM/shapeFunction_2D_Q1.*
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

Excercise
^^^^^^^^^^

Let’s try this interpolation scheme out on a single element with coordinates -1,-1 to 1,1.

.. toctree::
    :maxdepth: 2

    jupyter/FEM_2d_shapes_excercise.ipynb

2-D FEM solution of steady-state diffusion
------------------------------------------
Now we have assembled all the pieces that we need to make a finite element model and solve it. :numref:`lst:2d-fem-ss` shows the complete python code for our first solver.

Python code
^^^^^^^^^^^

.. code-block:: python 
    :linenos:
    :emphasize-lines: 60,65-85,103-110
    :name: lst:2d-fem-ss
    :caption: 2-D FEM diffusion, steady-state.

    import numpy as np
    from tabulate import tabulate
    from scipy.sparse.linalg import spsolve
    from scipy.sparse import csr_matrix
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    from shapes import shapes

    #geometry
    lx          = 1
    ly          = 1
    nx          = 51
    ny          = 51
    nnodel      = 4
    dx          = lx/(nx-1)
    dy          = ly/(ny-1)
    w_i         = 0.2 # width inclusion
    h_i         = 0.2 # heigths inclusion

    # model parameters
    k1          = 1
    k2          = 0.01
    Ttop        = 0
    Tbot        = 1
    
    nex         = nx-1
    ney         = ny-1
    nnod        = nx*ny
    nel         = nex*ney
    GCOORD      = np.zeros((nnod,2))
    T           = np.zeros(nnod) #initial T, not strictly needed

    # global coordinates
    id = 0
    for i in range(0,ny):
        for j in range(0,nx):
            GCOORD[id,0] = -lx/2 + j*dx
            GCOORD[id,1] = -ly/2 + i*dy
            id          = id + 1

    # FEM connectivity
    EL2NOD   = np.zeros((nel,nnodel), dtype=int)

    for iel in range(0,nel):
        row        = iel//nex   
        ind        = iel + row
        EL2NOD[iel,:] = [ind, ind+1, ind+nx+1, ind+nx]
        
    # Gauss integration points
    nip   = 4
    gauss = np.array([[ -np.sqrt(1/3), np.sqrt(1/3), np.sqrt(1/3), -np.sqrt(1/3)], [-np.sqrt(1/3), -np.sqrt(1/3), np.sqrt(1/3), np.sqrt(1/3)]]).T.copy()

    # Storage
    Rhs_all = np.zeros(nnod)
    I       = np.zeros((nel,nnodel*nnodel))
    J       = np.zeros((nel,nnodel*nnodel))
    K       = np.zeros((nel,nnodel*nnodel))

    for iel in range(0,nel):
        ECOORD = np.take(GCOORD, EL2NOD[iel,:], axis=0 )
        Ael    = np.zeros((nnodel,nnodel))
        Rhs_el = np.zeros(nnodel)
        
        for ip in range(0,nip):        
            # 1. update shape functions
            xi      = gauss[ip,0]
            eta     = gauss[ip,1]
            N, dNds = shapes(xi, eta)
            
            # 2. set up Jacobian, inverse of Jacobian, and determinant
            Jac     = np.matmul(dNds,ECOORD) #[2,nnodel]*[nnodel,2]
            invJ    = np.linalg.inv(Jac)     
            detJ    = np.linalg.det(Jac)
            
            # 3. get global derivatives
            dNdx    = np.matmul(invJ, dNds) # [2,2]*[2,nnodel]
            
            # 4. compute element stiffness matrix
            kel = k1
            if  abs(np.mean(ECOORD[:,0]))<w_i and abs(np.mean(ECOORD[:,1]))<h_i:
                kel=k2  
            Ael     = Ael + np.matmul(dNdx.T, dNdx)*detJ*kel # [nnodel,1]*[1,nnodel] / weights are missing, they are 1
            
            # 5. assemble right-hand side, no source terms, just here for completeness
            Rhs_el     = Rhs_el + np.zeros(nnodel)
        
        # assemble coefficients
        I[iel,:]  =  (EL2NOD[iel,:]*np.ones((nnodel,1), dtype=int)).T.reshape(nnodel*nnodel)
        J[iel,:]  =  (EL2NOD[iel,:]*np.ones((nnodel,1), dtype=int)).reshape(nnodel*nnodel)
        K[iel,:]  =  Ael.reshape(nnodel*nnodel)
        
        Rhs_all[EL2NOD[iel,:]] += Rhs_el

    A_all = csr_matrix((K.reshape(nel*nnodel*nnodel),(I.reshape(nel*nnodel*nnodel),J.reshape(nel*nnodel*nnodel))),shape=(nnod,nnod))

    # indices and values at top and bottom
    i_bot   = np.arange(0,nx, dtype=int)
    i_top   = np.arange(nx*(ny-1),nx*ny, dtype=int)
    Ind_bc  = np.concatenate((i_bot, i_top))
    Val_bc  = np.concatenate((np.ones(i_bot.shape)*Tbot, np.ones(i_top.shape)*Ttop ))

    # smart way of boundary conditions that keeps matrix symmetry
    Free    = np.arange(0,nnod)
    Free    = np.delete(Free, Ind_bc)
    TMP     = A_all[:,Ind_bc]
    Rhs_all = Rhs_all - TMP.dot(Val_bc)

    # solve reduced system
    T[Free] = spsolve(A_all[np.ix_(Free, Free)],Rhs_all[Free])
    T[Ind_bc] = Val_bc

    # postprocessing - heat flow
    Q_x  = np.zeros(nel)
    Q_y  = np.zeros(nel)
    Ec_x = np.zeros(nel)
    Ec_y = np.zeros(nel)

    for iel in range(0,nel):
        # 0. get element coordinates
        ECOORD = np.take(GCOORD, EL2NOD[iel,:], axis=0 )
    
        # 1. update shape functions
        xi      = 0
        eta     = 0
        N, dNds = shapes(xi, eta)
        
        # 2. set up Jacobian, inverse of Jacobian, and determinant
        Jac     = np.matmul(dNds,ECOORD) #[2,nnodel]*[nnodel,2]
        invJ    = np.linalg.inv(Jac)    
        detJ    = np.linalg.det(Jac)
        
        # 3. get global derivatives
        dNdx    = np.matmul(invJ, dNds) # [2,2]*[2,nnodel]
        
        # 4. heat flux per element
        kel = k1
        if  abs(np.mean(ECOORD[:,0]))<0.25 and abs(np.mean(ECOORD[:,1]))<0.25:
            kel=k2  
        Q_x[iel] = -kel*np.matmul(dNdx[0,:], np.take(T, EL2NOD[iel,:]))
        Q_y[iel] = -kel*np.matmul(dNdx[1,:], np.take(T, EL2NOD[iel,:]))
        Ec_x[iel]  = np.mean(ECOORD[:,0])
        Ec_y[iel]  = np.mean(ECOORD[:,1])

    # plotting
    fig = plt.figure()
    left, bottom, width, height = 0.1, 0.1, 0.8, 0.8
    ax = fig.add_axes([left, bottom, width, height]) 

    X = np.reshape(GCOORD[:,0], (ny,nx))
    Y = np.reshape(GCOORD[:,1], (ny,nx))
    T = T.reshape((ny,nx))

    cp = plt.contourf(X, Y, T)
    plt.colorbar(cp)
    plt.quiver(Ec_x, Ec_y, Q_x, Q_y, np.sqrt(np.square(Q_x) + np.square(Q_y)), cmap='autumn')

    ax.set_title('Temperature with heat flow vectors')
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    plt.show()import numpy as np

The key components of the code are an element loop (which corresponds to the sum over the integrals in :eq:`eq:fem_2d_weak_3` and an integration loop (which corresponds to the loop over integration points in :eq:`eq:fem_2d_weak_3`). Most of the work is done inside the integration loop:

    #. shape functions are calculated at the local coordinates of the integration point
    #. The Jacobian matrix, its determinant, and the inverse Jacobian are is calculated according to :eq:`eq:num_int_isoparam_4` 
    #. Global derivatives of shape functions at local coordinates are calculated according to :eq:`eq:num_int_jacobian`
    #. The element stiffness matrix and the element force vector are integrated.

    After the integration

    #. the local stiffness matrices and their global coefficients are stored
    #. and the element force vector is added to the global force vector


Let’s have a closer look at some of the different steps.

Jacobian:
^^^^^^^^^
The calculation of the Jacobian matrix is implemented via  :eq:`eq:num_int_isoparam_4`  as a matrix-vector multiplication:

.. math::
    :label: eq:jacobian_code

    J = 
    \begin{bmatrix}
    \frac{\partial N_1}{\partial \xi} & \frac{\partial N_2}{\partial \xi} & \frac{\partial N_3}{\partial \xi} & \frac{\partial N_4}{\partial \xi} \\
    \frac{\partial N_1}{\partial \eta} & \frac{\partial N_2}{\partial \eta} & \frac{\partial N_3}{\partial \eta} & \frac{\partial N_4}{\partial \eta}
    \end{bmatrix}
    =
    \begin{bmatrix}
    x_1 & y_1 \\
    x_2 & y_2 \\
    x_3 & y_3 \\
    x_4 & y_4 \\
    \end{bmatrix}
    =
    \begin{bmatrix}
    \frac{\partial x}{\partial \xi} & \frac{\partial y}{\partial \xi} \\
    \frac{\partial x}{\partial \eta} & \frac{\partial y}{\partial \eta}
    \end{bmatrix},

which is implemented as:

.. code-block:: python 
    
    ECOORD = np.take(GCOORD, EL2NOD[iel,:], axis=0 )
    Jac     = np.matmul(dNds,ECOORD) #[2,nnodel]*[nnodel,2]


Global derivatives of shape functions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
The global derivatives are calculated according to :eq:`eq:num_int_jacobian`, which in matrix form is:

.. math::
    :label: eq:dndx_code

    \begin{bmatrix}
    \frac{\partial N_1}{\partial x} & \frac{\partial N_2}{\partial x} & \frac{\partial N_3}{\partial x} & \frac{\partial N_4}{\partial x} \\
    \frac{\partial N_1}{\partial y} & \frac{\partial N_2}{\partial y} & \frac{\partial N_3}{\partial y} & \frac{\partial N_4}{\partial y}
    \end{bmatrix}
    =
    \begin{bmatrix}
    J^{-1}
    \end{bmatrix}
    \begin{bmatrix}
    \frac{\partial N_1}{\partial \xi} & \frac{\partial N_2}{\partial \xi} & \frac{\partial N_3}{\partial \xi} & \frac{\partial N_4}{\partial \xi} \\
    \frac{\partial N_1}{\partial \eta} & \frac{\partial N_2}{\partial \eta} & \frac{\partial N_3}{\partial \eta} & \frac{\partial N_4}{\partial \eta}
    \end{bmatrix},

which is implemented as:

.. code-block:: python 
        
    dNdx    = np.matmul(invJ, dNds) # [2,2]*[2,nnodel]

Element stiffness matrix
^^^^^^^^^^^^^^^^^^^^^^^^
The final step is to put the :eq:`eq:num_int_jacobian_2` into finite element form. Also this now quite easy:

.. math::
    :label: eq:Ael_code

    \begin{bmatrix}
    \begin{bmatrix}
    \frac{\partial N_1}{\partial x} & \frac{\partial N_1}{\partial y} \\
    \frac{\partial N_2}{\partial x} & \frac{\partial N_2}{\partial y} \\
    \frac{\partial N_3}{\partial x} & \frac{\partial N_3}{\partial y} \\
    \frac{\partial N_4}{\partial x} & \frac{\partial N_4}{\partial y}
    \end{bmatrix}
    \begin{bmatrix}
    \frac{\partial N_1}{\partial x} & \frac{\partial N_2}{\partial x} & \frac{\partial N_3}{\partial x} & \frac{\partial N_4}{\partial x} \\
    \frac{\partial N_1}{\partial y} & \frac{\partial N_2}{\partial y} & \frac{\partial N_3}{\partial y} & \frac{\partial N_4}{\partial y}
    \end{bmatrix}
    \end{bmatrix}
    \begin{bmatrix}
    T_1 \\ T_2 \\ T_3 \\ T_4
    \end{bmatrix}
    = 
    \begin{bmatrix}
    0 \\ 0 \\ 0 \\ 0
    \end{bmatrix},

where the first matrix is the element stiffness matrix, which is implemented as:

.. code-block:: python 
        
    Ael     = Ael + np.matmul(dNdx.T, dNdx)*detJ*kel # [nnodel,1]*[1,nnodel] / weights are missing, they are 1

Boundary conditions
^^^^^^^^^^^^^^^^^^^^
Before we can solve anything, we need to set boundary conditions. In the 1-D case, we had set the boundary conditions in the usual way by setting the entire row to zero, putting a 1 on the main diagonal, and setting the Rhs to the specified temperature. We had noted before that this procedure breaks matrix symmetry. We will therefore use a more elaborate method here, which was described in :cite:`Dabrowski2008` .
We add the columns that refer to nodes with boundary conditions to the righ-hand side and afterwards solve the smaller system of equations! This is implemented as

.. code-block:: python

    # indices and values at top and bottom
    i_bot   = np.arange(0,nx, dtype=int)
    i_top   = np.arange(nx*(ny-1),nx*ny, dtype=int)
    Ind_bc  = np.concatenate((i_bot, i_top))
    Val_bc  = np.concatenate((np.ones(i_bot.shape)*Tbot, np.ones(i_top.shape)*Ttop ))

    # smart way of boundary conditions that keeps matrix symmetry
    Free    = np.arange(0,nnod)
    Free    = np.delete(Free, Ind_bc)
    TMP     = A_all[:,Ind_bc]
    Rhs_all = Rhs_all - TMP.dot(Val_bc)

    # solve reduced system
    T[Free] = spsolve(A_all[np.ix_(Free, Free)],Rhs_all[Free])
    T[Ind_bc] = Val_bc

First we find all the global indices of the nodes with boundary conditions and build a vector of equal length with all the fixed values. Afterwards we move the columns to the Rhs and solve the smaller system of equations.
Note how the Rhs is a vector full of zeros at the moment. This will change when we look at the transient problem!

.. admonition:: Excercise

    You can download the python code here (:download:`2d_fem_steady_state.py <python/2d_fem_steady_state.py>` and :download:`shapes.py <python/shapes.py>`). Try to get it to work and make sure that you understand the key parts. I good way is to use the Visual Studio Code debugger and to step through the code. Make sure that you understand the various matrix shapes and the connection between the equations presented above and the matrix multiplications in the code.

Excersise - unstructured triangle meshes
----------------------------------------
In the previous section we have learned how to solve for steady-state diffusion problems in 2-D using FEM. One of strengths of the FE method is that it is quite easy to use differently shaped elements. Triangles are, for example, well suited for resolving complex geometries in 2-D.

In this excercise we will learn

    * how to use a mesh genetrator to mesh complex geometries,
    * what changes are needed to change the basic element.

The problem we will look at is basically the same as in the example given above. We will resolve steady-state heat diffusion within a rectangular region but now we will have circular inclusions of different thermal conductivity.

Theoretical background
^^^^^^^^^^^^^^^^^^^^^^^^
Let's first explore how triangle element "work". :numref:`fig:fem_triangle_coords_f`  shows a triangle in global and local coordinates. The local node numering is again counterclockwise.

.. figure:: /_figures/triangle_coords.png
    :name: fig:fem_triangle_coords_f
    :align: center
    
    Global and local coordinates of a triangle element. Black crosses mark the locations of the three integration points.

Here we look at a linear triangle with three nodes. The relation between global and local coordinates is therefore:

.. math::
    :label: eq:triangle_1

    \begin{align}
    \begin{split}
    x &= x_1 + (x_2-x_1)\xi + (x_3-x_1)\eta \\
    y &= y_1 + (y_2-y_1)\xi + (y_3-y_1)\eta \\
    \end{split}
    \end{align}


Temperature, or any other variable, will also vary linearly with the local coordinates and we can write

.. math::
    :label: eq:triangle_2

    T(\xi,\eta) = a_1 +a_2\xi + a_3\eta

We can spell this out for the three nodes:

.. math::
    :label: eq:triangle_3

    \begin{align}
    \begin{split}
    T_1 &= T(0,0) = a_1 \\
    T_2 &= T(1,0) = a_1 + a_2 \\
    T_3 &= T(0,1) = a_1 + a_3
    \end{split}
    \end{align}

And express the contants :math:`a` in terms of the nodal temperatures (or whatever our unknown is). 

.. math::
    :label: eq:triangle_4

    \begin{align}
    \begin{split}
    a_1 = T_1 \\
    a_2 = T_2 - T_1 \\
    a_3 = T_3 - T_1
    \end{split}
    \end{align}

This we put back into the definition of the approximate solution :eq:`eq:triangle_2` and get:


.. math::
    :label: eq:triangle_4

    \begin{align}
    \begin{split}
    T(\xi,\eta) &= T_1 +(T_2 - T_1)\xi + (T_3-T_1)\eta\\
    \Rightarrow \\
    T(\xi,\eta) &= (1 -\xi - \eta)T_1 + \xi T_2 + \eta T_3 \\
    \Rightarrow \\
    N_1(\xi, \eta) &=  (1 -\xi - \eta)\\
    N_2(\xi, \eta) &=  \xi \\
    N_3(\xi, \eta) &=  \eta \\
    \Rightarrow \\
    T(\xi,\eta) &= N_1 T_1 + N_2 T_2 + N_3 T_3 \\
    \end{split}
    \end{align}

Alright, now that we understand the local coordinates and the shape functions associated with linear triangle element, we are good to go!


Step 0: getting ready
^^^^^^^^^^^^^^^^^^^^^
Make a copy of the exemple solver above. Take it as a starting point and integrate the various code pieces and pieces of information below into. 

Step 1: mesh generation
^^^^^^^^^^^^^^^^^^^^^^^
We will use the mesh generator `triangle <https://www.cs.cmu.edu/~quake/triangle.html>`_ by Jonathan Shewchuk. It's one of the best 2-D mesh generators for triangle meshes. It's originally written in C but we will for convenience use a python wrapper. You can install it into your virtual environment by doing this:

.. code-block:: bash

    pip install triangle

Here is a code piece for making the the mesh:

.. code-block:: python

    import triangle as tr

    ## Create the triangle mesh
    vertices = []
    segments = []
    regions = []

    # make a box with given dims and place given attribute at its center
    def make_box(x, y, w, h, attribute):
        i = len(vertices)

        vertices.extend([[x,   y],
                        [x+w, y],
                        [x+w, y+h],
                        [x,   y+h]])

        segments.extend([(i+0, i+1),
                        (i+1, i+2),
                        (i+2, i+3),
                        (i+3, i+0)])
        
        regions.append([x+0.01*w, y+0.01*h, attribute,0.005])

    def make_inclusion(center_x, center_y, radius, points_inc, attribute):
        theta        = np.linspace(0,2*np.pi,points_inc, endpoint=False)
        xx           = np.cos(theta)
        yy           = np.sin(theta)

        i = len(vertices)   
        vertices.extend(np.array([center_x + radius*xx,center_y + radius*yy]).T)
        
        Tmp = np.array([np.arange(i, i+points_inc), np.arange(i+1, i+points_inc+1)]).T
        Tmp[-1,1] = i
        segments.extend(Tmp)
    
        regions.append([center_x, center_y, attribute,0.001])

    #geometry
    x0          = -1
    y0          = -1
    lx          = 2
    ly          = 2
    n_incl      = 5
    radius      = 0.15

    # generate input    
    make_box(x0, y0, lx, ly, 1)

    make_inclusion(-0.8, -0.3, radius, 20, 100)
    make_inclusion(-0.5, -0.75, radius, 20, 100)
    make_inclusion(-0.6, 0.5, radius, 20, 100)
    make_inclusion(-0.1, -0.3, radius, 20, 100)
    make_inclusion(0.1, 0, radius, 20, 100)
    make_inclusion(0.5, -0.2, radius, 20, 100)
    make_inclusion(0.6, .3, radius, 20, 100)
    make_inclusion(0.7, .8, radius, 20, 100)
    make_inclusion(0, .75, radius, 20, 100)
    make_inclusion(-0.5, .05, radius, 20, 100)
    make_inclusion(0.5, -.75, radius, 20, 100)

    A = dict(vertices=vertices, segments=segments, regions=regions)
    B = tr.triangulate(A, 'pq33Aa')
    #tr.compare(plt, A, B)
    #plt.show()

    # extract mesh information
    GCOORD = B.get("vertices")
    EL2NOD = B.get("triangles")
    Phases = B.get("triangle_attributes")

    nnodel = EL2NOD.shape[1]
    nel    = EL2NOD.shape[0]
    nnod   = GCOORD.shape[0]
    Phases = np.reshape(Phases,nel)


Note how the generated mesh comes back as a python dictionary from which we extract the mesh information. Note also that the array :code:`Phases` contains a marker to which region an element belongs (matrix versus inclusion).

Step 2: triangle shape functions and integration points
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
:numref:`fig:shapeFunc:2D:linear` shows shape functions for a linear triangle element. You will need to modify the function in our :code:`shapes.py` file to implement the triangular shape functions. 

.. code-block:: python

    #shape functions
    eta2 = xi
    eta3 = eta
    eta1 = 1-eta2-eta3

    N1 = eta1
    N2 = eta2
    N3 = eta3

Take this as a starting point and modify :code:`shapes.py` to return the correct three shape functions and the six derivatives!

Next we need to adapt our integration rule. Take these three integration points and weights.

.. list-table:: Triangle integration
    :header-rows: 1

    * - Integration point
      - :math:`\xi`
      - :math:`\eta`
      - weight
    * - :math:`1`
      - :math:`\frac{1}{6}`
      - :math:`\frac{1}{6}`
      - :math:`\frac{1}{6}`
    * - :math:`2`
      - :math:`\frac{2}{3}`
      - :math:`\frac{1}{6}`
      - :math:`\frac{1}{6}`
    * - :math:`3`
      - :math:`\frac{1}{6}`
      - :math:`\frac{2}{3}`
      - :math:`\frac{1}{6}`


.. tip::
    There are, of course, many different elements and associated integration rules using different numbers of integration points in FEM. Have a look at the suggested books in :ref:`Course details` .



Step 3: Boundary conditions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^
We were actually a bit "lazy" when we implemented the boundary conditions in the example above. Instead of coming up with a general solution, we identified the global node numbers on the boundaries assuming a structured quad mesh. Now we are using an unstructured triangle mesh and it is not so easy anymore to know which nodes are on the boundaries and should get boundary conditions. The "typical" way would be to use boundary markers in the mesh generation; here we use a different approach and use an ad-hoc search for the boundary nodes:

.. code-block:: python

    # indices and values at top and bottom
    tol = 1e-3
    # i_bot = np.where(abs(???) < tol)[0] #bottom nodes
    # i_top = np.where(abs(???) < tol)[0] #top nodes
    

Complete this and add the lines to the main code!

.. admonition:: Boundary conditions

    Did you notice that we are not specifying boundary conditions for lateral boundaries? Still the code seems to work. Have a look at :eq:`eq:fem_2d_weak_2` and think about which implicit assumption we are making about the line integral.


Step 4: Post-processing and plotting
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Modify the post-processing code that computes the heat fluxes, so that it works for triangles. This should just involve setting a new local coordinates for the single integration points. Use :math:`\xi=\frac{1}{3}` and  :math:`\eta=\frac{1}{3}`, which is the center of each triangle.

Now we just need to plot! Thankfully python makes it easy for us and we can use the functions :code:`triplot` and :code:`tricontourf` to plot the unstructured data.

.. code-block:: python

    # plotting
    fig = plt.figure()
    left, bottom, width, height = 0.1, 0.1, 0.8, 0.8
    ax = fig.add_axes([left, bottom, width, height]) 

    plt.triplot(GCOORD[:,0], GCOORD[:,1], EL2NOD, linewidth=0.5)
    cp = plt.tricontourf(GCOORD[:,0], GCOORD[:,1], EL2NOD, T, 10, cmap='gist_yarg')

    plt.colorbar(cp)
    plt.quiver(Ec_x, Ec_y, Q_x, Q_y, np.sqrt(np.square(Q_x) + np.square(Q_y)), cmap='hot')

    ax.set_title('Temperature with heat flow vectors')
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    plt.show()


.. figure:: /_figures/triangle_fem_diffusion.png
    :name: fig:fem_triangle_dif_f
    :align: center
    :figwidth: 80%

    2-D diffusion on unstructured mesh.
