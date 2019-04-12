function [u_n,nffast,nfslow] = solve_MERK4(A,gn,internalsolver,finalstepsolver, u0, m, tvals,h)
  % Usage: [u_n,nffast,nfslow] = solve_MERK4(A,gn,internalsolver,finalstepsolver, u0, m, tvals,h)
  %
  % Fifth order multirate exponential Runge Kutta solver for the ODE problem
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
  nffast = 0;       % fast function calls
  nfslow = 0;       % slow function calls

  % Set problem parameters
  c       = [1/2,1/2,1/3,5/6,1/3];
  c_2      = c(1);
  c_3      = c(2);
  c_4      = c(3);
  c_5      = c(4);
  c_6      = c(5);
  n        = 0;
  u_n      = u0;
  t_n      = tvals(1);
  ONEMSM   = 1-sqrt(eps);

  while t_n < tvals(2)*ONEMSM

    % Set initial condition
    Y0 = u_n;

    % Set up right hand side for modified ODE to solve for U_{n,2}
    p_n2 = gn(t_n,u_n);
    fcn = @(t,y) A*y + p_n2;

    % Set up micro time step
    h_fast = c_2*h/ceil(c_2*m);

    % Solve for U_{n,2}
    [~,Y,nflocal] = solve_ERKfast(fcn,[0,c_2*h],Y0,B,h_fast);
    U_n2 = Y(:,2);

    % Update number of fast function calls
    nffast   = nffast + nflocal;

    % Define slow function contribution function
    D_ni   = @(t,c,U) gn(t + c*h,U) - p_n2;

    % Solve for D_n2
    D_n2 = D_ni(t_n,c_2,U_n2);

    % Set up right hand side for modified ODE to solve for U_{n,3}, U_{n,4}
    p_n3 = @(t) p_n2 + (t/(h*c_2))*D_n2;
    fcn = @(t,y) A*y + p_n3(t);

    % Determine time step, and times to evaluate at
    h_fast = c_3*h/ceil(c_3*m);
    tsteps = [0,c_4*h,c_3*h];

    %Solve for U_{n,3}, U_{n,4}
    [~,Y,nflocal] = solve_ERKfast(fcn,tsteps,Y0,B,h_fast);
    U_n4 = Y(:,2);
    U_n3 = Y(:,end);

    % Update number of fast function calls
    nffast   = nffast + nflocal;

    % Solve for D_n3 and D_n4
    D_n3 = D_ni(t_n,c_3,U_n3);
    D_n4 = D_ni(t_n,c_4,U_n4);

    % Set up right hand side for modified ODE to solve for U_{n,5}, U_{n,6}
    p_n5 = @(t) p_n2 + (t/h)*(-c_4/(c_3*(c_3-c_4)))*D_n3 +...
    (t/h)*(c_3/(c_4*(c_3-c_4)))*D_n4 +...
    (t^2/(2*h^2))*(2/(c_3*(c_3-c_4)))*D_n3 - ...
    (t^2/(2*h^2))*(2/(c_4*(c_3-c_4)))*D_n4;
    fcn = @(t,y) A*y + p_n5(t);

    % Determine time step, and times to evaluate at
    h_fast = c_5*h/ceil(c_5*m);
    tsteps = [0,c_6*h,c_5*h];

    %Solve for U_{n,5}, U_{n,6}
    [~,Y,nflocal] = solve_ERKfast(fcn,tsteps,Y0,B,h_fast);
    U_n6 = Y(:,2);
    U_n5 = Y(:,end);

    % Update number of fast function calls
    nffast   = nffast + nflocal;

    % Solve for D_n5 and D_n6
    D_n5 = D_ni(t_n,c_5,U_n5);
    D_n6 = D_ni(t_n,c_6,U_n6);

    % Update number of slow function calls
    % p_n2,D_n2,D_n3,D_n4,D_n5,D_n6
    nfslow   = nfslow + 6;

    % Set up right hand side for modified ODE to solve for  u_{n+1}
    q_n6 = @(t) p_n2 + (t/h)*(-c_6/(c_5*(c_5-c_6)))*D_n5 +...
    (t/h)*(c_5/(c_6*(c_5-c_6)))*D_n6 +...
    (t^2/(2*h^2))*(2/(c_5*(c_5-c_6)))*D_n5 - ...
    (t^2/(2*h^2))*(2/(c_6*(c_5-c_6)))*D_n6;
    fcn = @(t,y) A*y + q_n6(t);

    % Set up micro time step
    h_fast = h/m;

    % Solve for u_{n+1}
    [~,Y,nflocal] = solve_ERKfast(fcn,[0,h],Y0,D,h_fast);
    u_np1 = Y(:,2);

    % Update number of fast function calls
    nffast = nffast + nflocal;

    % Update time step
    t_n = t_n + h;
    n   = n+1;

    % Update u value
    u_n = u_np1;

  end
end
