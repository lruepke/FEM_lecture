Example 1: Cooling dike
------------------------

As a first example, we will use the finite differences form of the heat diffusion equation :eq:`eq:1D_heat_flow` to explore the cooling of a dike. We will look at a  :math:`2m` wide dike that intruded with a temperature of  :math:`1200°C` into  :math:`300°C` warm country rock. The initial conditions can then be written like this:

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

Figure \ref{fig:dike_setup} summarized the setup of the cooling dike problem:

\begin{figure}[htb]
	\center{\includegraphics[width=\textwidth]{figures/1D_dike_example_stencil.png}}
%	\vspace{3in}
	\caption{\label{fig:dike_setup} A) Setup of the model considered here. A hot basaltic dike intrudes into colder country rock. Only variations in x-direction are considered; properties in the other directions are assumed to be constant. The initial temperature distribution $T(x,0)$ has a step-like perturbation. B) Finite difference discretization of the 1D heat equation. The finite difference method approximates the temperature at given grid points, with spacing $\Delta x$. The time-evolution is also computed at given times with timestep $\Delta t$.
	}
\end{figure}

\subsection{Exercises}

\begin{enumerate}
	\item Open MATLAB and an editor and type the Matlab script in an empty file. Save the file under the name ”heat1Dexplicit.m”. Fill in the question marks and run the file by typing heat1Dexplicit in the MATLAB command window (make sure you’re in the correct directory).
	\item Vary the parameters (e.g. use more gridpoints, a larger timestep). Note that if the timestep is increased beyond a certain value (what does this value depend on?), the numerical method becomes unstable. This is a major drawback of explicit finite difference codes such as the one presented here. In the next lesson we will learn methods that do not have these limitations.
	\item Go through the rest of the handout and see how one derives finite difference approximations.
	\item Record and plot the temperature evolution versus time at a distance of $5 \unt{m}$ from the dike/country rock contact. What is the maximum temperature the country rock experiences at this location and when is it reached?
	Bonus question: Derive a finite-difference approximation for variable $k$ and variable $\Delta x$.
\end{enumerate}
\pagebreak

\lstset{language=Matlab, frame=single, basicstyle=\footnotesize} %or \small or \footnotesize etc.}
\begin{lstlisting}
%heat1Dexplicit.m
%
% Solves the 1D heat equation with an explicit finite difference scheme
clear

%Physical parameters
L       = 100; 				% Length of modeled domain [m]
Tmagma  = 1200; 			% Temperature of magma [C]
Trock   = 300; 				% Temperature of country rock [C]
kappa   = 1e-6; 			% Thermal diffusivity of rock [m2/s]
W       = 5; 				% Width of dike [m]
day     = 3600*24; 			% # seconds per day
dt      = 1*day; 			% Timestep [s]

% Numerical parameters
nx      = 201; 				% Number of gridpoints in x-direction
nt      = 500; 				% Number of timesteps to compute
dx      = L/(nx-1); 			% Spacing of grid
X       = linspace(-L/2,L/2,nx);	% Grid

% Setup initial temperature profile
T       = ones(nx)*Trock;
T(abs(X)<=W/2) = Tmagma;
time    = 0;

%Main time loop
for n=1:nt % Timestep loop
	% Compute new temperature
	Tnew = zeros(1,nx);
	for i=2:nx-1
		Tnew(i) = T(i) + ?????;
	end

	% Set boundary conditions
	Tnew(1)     = T(1);
	Tnew(nx)    = T(nx);

	% Update temperature and time
	T           = Tnew;
	time        = time+dt;

	% Plot solution
	figure(1), clf
	plot(X,Tnew);
	xlabel('x [m]')
	ylabel('Temperature [^\circ C]')
	title(['Temperature evolution after ',num2str(time/day),' days'])
	drawnow
end
\end{lstlisting}
