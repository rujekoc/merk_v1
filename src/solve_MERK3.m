function [u_n,nffast,nfslow] = solve_MERK3(A,gn,internalsolver,finalstepsolver, u0, m, tvals,h)
  % Usage: [u_n,nffast,nfslow] = solve_MERK3(A,gn,internalsolver,finalstepsolver, u0, m, tvals,h)
  %
  % Third order multirate exponential Runge Kutta solver for the  ODE problem
  %     u' = Au(t_ + g(t,u(t)), t in tvals,
  %     u(t0) = u0
  %
  % Inputs:
  %     A          =  matrix A
  %     gn         = string holding function name for g(t,u(t))
  %    internalsolver  = string holding function name for  RK method
  %                  can be one less order than MERK solver
  %                 Butcher table formatting
  %                 B = [c A;
  %                      q b]
  %                 Here, c is a vector of stage time fractions (s-by-1),
  %                    A is a matrix of Butcher coefficients (s-by-s),
  %                    q is an integer denoting the method order of accuracy,
  %                    b is a vector of solution weights (1-by-s),
  %    finalstepsolver  = string holding function name for RK method
  %                  has to be same order as MERK solver
  %     tvals      = [t0,tN] initial and final time
  %     u0         = initial value
  %     h          = slow time step
  %
  % Outputs:
  %     u_n    = solution at final time
  %     nffast = number of fast function calls
  %     nfslow = number of slow function calls
  %
  %
  % Rujeko Chinomona
  % Department of Mathematics
  % Southern Methodist University
  % April 2019
  
  % Info from butcher table for inner ODE solvers
  B = butcher(internalsolver);
  D = butcher(finalstepsolver);

  % Set problem parameters
  nffast = 0;
  nfslow = 0;

  c        = [1/2,2/3];
  c_2      = c(1);
  c_3      = c(2);
  n        = 0;
  u_n      = u0;
  t_n      = tvals(1);
  ONEMSM   = 1-sqrt(eps);

  while t_n < tvals(2)*ONEMSM

    % Set initial condition
    Y0  = u_n;

    % Determine micro time step
    h_fast = c_2*h/ceil(c_2*m);

    % Set up right hand side for modified ODE to solve for U_{n,2}
    p_n2 = gn(t_n,u_n);
    fcn  = @(t,y) A*y + p_n2;

    % Solve for U_{n,2}
    [~,Y,nflocal] = solve_ERKfast(fcn,[0,c_2*h],Y0,B,h_fast);
    U_n2 = Y(:,2);

    % Update number of fast function calls
    nffast = nffast + nflocal;

    % Define slow function contribution function
    D_ni   = @(t,c,U) gn(t+c*h,U) - p_n2;

    % Determine micro time step
    h_fast = c_3*h/ceil(c_3*m);

    % Solve for D_n2
    D_n2 = D_ni(t_n,c_2,U_n2);

    % Set up right hand side for modified ODE to solve for  U_{n,3}
    p_n3 = @(t) p_n2 + 4*t/(9*c_2*c_3^2*h)*D_n2;
    fcn = @(t,y) A*y + p_n3(t);

    % Solve for U_{n,3}
    [~,Y,nflocal] = solve_ERKfast(fcn,[0,c_3*h],Y0,B,h_fast);
    U_n3 = Y(:,2);

    % Update number of fast function calls
    nffast = nffast + nflocal;

    % Solve for D_n3
    D_n3 = D_ni(t_n,c_3,U_n3);

    % Update number of slow function calls
    % p_n2,D_n2,D_n3
    nfslow   = nfslow + 3;

    % Set up right hand side for modified ODE to solve for u_{n+1}
    q_n3 = @(t) p_n2 + (t/h)*(-2/(3*c_2*(c_2-c_3)))*D_n2 + ...
    (t/h)*(c_2/(c_3*(c_2-c_3)))*D_n3 + ...
    (t^2/(2*h^2))*(2/(c_2*(c_2-c_3)))*D_n2 - ...
    (t^2/(2*h^2))*(2/(c_3*(c_2-c_3)))*D_n3;

    fcn = @(t,y) A*y + q_n3(t);

    % Set up micro time step
    h_fast = h/m;

    % Solve for u_{n+1} on [0,h]
    [~,Y,nflocal] = solve_ERKfast(fcn,[0,h],Y0,D,h_fast);
    u_np1 = Y(:,2);

    % Update number of fast function calls
    nffast = nffast + nflocal;

    % Update time step
    t_n = t_n + h;
    n = n+1;

    % Update u value
    u_n = u_np1;

  end
end
