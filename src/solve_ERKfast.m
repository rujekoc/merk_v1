function [tvals,Y,nf] = solve_ERKfast(fcn,tvals,Y0,B,hmax)
% usage: [tvals,Y,nf] = solve_ERKfast(fcn,tvals,Y0,B,hmax)
%
% THIS VERSION RETURNS SOLUTION AT MULTIPLE TIME STEPS
%
% Explicit Runge-Kutta solver for the vector-valued ODE problem
%     y' = F(t,Y), t in tvals, y in R^m,
%     Y(t0) = [y1(t0), y2(t0), ..., ym(t0)]'.
%
% Inputs:
%     fcn    = string holding function name for F(t,Y)
%     tvals  = array with times to output solution at.
%     Y0     = initial value array (column vector of length m)
%     B      = Butcher matrix for ERK coefficients, of the form
%                 B = [c A;
%                      q b]
%              Here, c is a vector of stage time fractions (s-by-1),
%                    A is a matrix of Butcher coefficients (s-by-s),
%                    q is an integer denoting the method order of accuracy,
%                    b is a vector of solution weights (1-by-s),
%     h      = time step
%
% Outputs:
%     tvals  = the same as the input array tvals
%     y      = [y(t0),y(t1),...,y(tN)], where each
%               y(t*) is a column vector of length m.
%     nf = number of function calls
%
% Rujeko Chinomona
% Department of Mathematics
% Southern Methodist University
% May 2018
% Adjusted from the code by D. Reynolds ("solve_ERK.m", Aug 2012)


% extract ERK method information from B
[~, Bcols] = size(B);
s = Bcols - 1;        % number of stages
c = B(1:s,1);         % stage time fraction array
A = B(1:s,2:s+1);     % RK coefficients
b = (B(s+1,2:s+1))';


% initialize output arrays
N = length(tvals);
m = length(Y0);
Y = zeros(m,N);
Y(:,1) = Y0;

% set the solver parameters
ONEMSM   = 1-sqrt(eps);  % coefficients to account for

% initialize temporary variables
t = tvals(1);
Ynew = Y0;

% initialize work counter
nsteps = 0;

% iterate over output time steps
for tstep = 2:length(tvals)
    h = hmax;
   % loop over internal time steps to get to desired output time
   while (t < tvals(tstep)*ONEMSM)
     h = min([h, tvals(tstep)-t]);  % stop at output time

     % initialize storage for RHS vectors
     k = zeros(length(Y0),s);

     % loop over stages
     for stage=1:s
       % construct stage solution and evaluate RHS
       %    zi = y_n + h*sum_{j=1}^{i-1} (A(i,j)*f(zj))
       z = Y0;
       for j=1:stage-1
          z = z + h*A(stage,j)*k(:,j);
       end

       % construct new stage RHS
       k(:,stage) = fcn(t+h*c(stage),z);
     end

      % compute new solution and error estimate
      %    ynew = yold + h*sum(b(j)*fj)
      Ynew  = Y0 + h*k*b;

      % increment number of internal time steps taken
      nsteps = nsteps + 1;

      % update solution and time for last successful step
      Y0 = Ynew;
      t  = t + h;

   end  % end while loop attempting to solve steps to next output time

   % store updated solution in output array
   Y(:,tstep) = Ynew;

end  % time step loop
   nf = nsteps*s;

% end solve_ERK function
end
