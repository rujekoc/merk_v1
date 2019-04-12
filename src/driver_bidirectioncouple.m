% USAGE : driver_bidirectioncouple
% driver for bi-directional coupling test problem:
%      u' = 100v + w,
%      v' =  -100u,
%      w' = -w + u,
% where u(0) = 9001/10001, v(0) = -10^5/10001,w(0) = 1000.
% We evaluate over the time interval [0,2].
% MERK methods and the MIS-KW3 method are tested and compared against a reference
% solution.
% Inner solvers can be changed to test different elements of algorithms.
% Plotting code can be commented out.
%
% Rujeko Chinomona
% Department of Mathematics
% Southern Methodist University
% April 2019

clear, close all

% Splitting for MERK methods
A   = [0,100,0;-100,0,0;1,0,0];                      % fast linear part
gn  = @(t,y) [y(3);0;-y(3)];                         % slow variable

% Splitting for MIS-KW3
ff  = @(t,y) [100*y(2);-100*y(1);y(1)];              % fast
fs  = @(t,y) [y(3); 0;-y(3)];                        % slow
fn  = @(t,y) A*y + gn(t,y);                          % full rhs


% Set time parameters for algorithm
Ti    = 0;                                         % start time
Tf    = 2.0;                                       % end time
n     = 21;
tout  = linspace(Ti,Tf,n);                         % intermediate times for solution
hmax  = (tout(2)-tout(1))/4;                       % largest macro/slow time step
h     = hmax*0.5.^(0:7);                           % macro/slow scale timesteps
Y0    = [9001/10001;-1e5/10001;1000];              % initial condition

% Read in reference solution
filename = 'bidirectionrefsolution.mat';
q        = matfile(filename);
Ytrue    = q.Yirk;

% Initialize problem variables/ allocate space
Y          = zeros(length(Y0),n);
Y(:,1)     = Y0;
err_max    = zeros(1,length(h));
nffast     = zeros(1,length(h));
nfslow     = zeros(1,length(h));
time       = zeros(1,length(h));

% Solver info
solvers      = {@solve_MERK3,@solve_MERK4,@solve_MERK5,@solve_MIS_KW3};
snames       = {'MERK3','MERK4','MERK5','MIS-KW3'};
internalsolvers  = {'ERK-3-3','ERK-4-4','Cash-Karp-ERK','Knoth-Wolke-ERK'};
finalstepsolvers  = {'ERK-3-3','ERK-4-4','Cash-Karp-ERK','Knoth-Wolke-ERK'};
orders   = [3,4,5,3];
mvalues  = [50,50,10,25];   % predetermined time scale separation factors

% Initialize parameters for plotting
colors  = {[0,0.7461,1],[1,0,0],[0.1953,0.8008,0.1953],[1,0,1]};
markers = {'^','<','>','v'};
figure(1);
hold on;
figure(2);
hold on;
figure(3);
hold on;

% Start iterating over solvers
for isol = 1:length(solvers)
  % Get information about solvers and print it out
  solver = solvers{isol};
  sname  = snames{isol};
  internalsolver  = internalsolvers{isol};
  finalstepsolver  = finalstepsolvers{isol};
  fprintf('\nRunning convergence test with %s integrator (theoretical order %i)\n',...
  snames{isol},orders(isol));
  B = butcher(internalsolver);  s = numel(B(1,:))-1;
  if strcmp(sname,'MIS-KW3')
    fprintf('\nRunning with inner ERK integrator: %s (order = %i)\n',internalsolver,B(s+1,1));
  else
    D = butcher(finalstepsolver);  s1 = numel(D(1,:))-1;
    fprintf('\nRunning with inner ERK integrators: %s (order = %i)  and %s (order = %i)\n',...
    internalsolver,B(s+1,1), finalstepsolver,D(s1+1,1));
  end

  % Set function calls
  nffast     = zeros(1,length(h));
  nfslow     = zeros(1,length(h));

  % Compute solution using multirate solver using different h values
  for j = 1:length(h)
    m = mvalues(isol);
    hfast = h(j)/m;
    if strcmp(sname,'MIS-KW3')
      % Computation for MIS
      tstart = tic;
      [t,Yout,fastcalls,slowcalls] = solver(fs, ff, tout, Y0, h(j), hfast);
      telapsed = toc(tstart);

      Y = Yout;
      nffast(j) = fastcalls;
      nfslow(j) = slowcalls;
    else
      % Computation for MERK methods
      tstart = tic;
      for i = 2:n
        Y0step     = Y(:,i-1);
        [Yout,fastcalls,slowcalls]= solver(A, gn, internalsolver,finalstepsolver,...
        Y0step,m,[tout(i-1),tout(i)],h(j));
        Y(:,i) = Yout;
        nffast(j) = nffast(j)+ fastcalls;
        nfslow(j) = nfslow(j)+ slowcalls;
      end
      telapsed = toc(tstart);
    end
    time(j) = telapsed;

    % Function calls
    fprintf('nffast = %g ,  nfslow = %g .\n',nffast(j),nfslow(j));

    % Error calculation
    err_max(j) = max(max(abs(Y - Ytrue)));
    fprintf('Accuracy/Work Results, h = %g  , m = %i:\n',h(j),m)
    fprintf('   maxerr = %.5e,    ',err_max(j));

    % Rate of convergence calculation
    if j > 1
      max_rate(j-1) = log(err_max(j-1)/err_max(j))/log(h(j-1)/h(j));
      fprintf('   maxrate = %.5e \n',max_rate(j-1));
    end
    fprintf('Time elapsed = %g.\n',telapsed);
  end

  % Best-fit convergence rate (includes stagnating pts &/floating pt round off pts)
  p = polyfit(log(h),log(err_max),1);
  fprintf('best-fit convergence rate = %g\n', p(1));


  % Plot Error vs. h
  figure(1)
  hold all
  loglog(h,err_max,[markers{isol},'-'],'color',colors{isol},...
  'LineWidth',3,'MarkerFaceColor',colors{isol})
  header1{isol} = sname;
  hold on;

  % Plot Error vs. Total function calls
  figure(2)
  hold all
  loglog(nfslow+nffast,err_max,[markers{isol},'-'],'color',colors{isol},...
  'MarkerFaceColor',colors{isol},'LineWidth',3);
  header2{isol} = sname;
  hold on ;

  % Plot Error vs. Slow function calls
  figure(3)
  hold all
  loglog(nfslow,err_max,[markers{isol},'-'],'color',colors{isol},...
  'MarkerFaceColor',colors{isol},'LineWidth',3);
  header3{isol} = sname;
  hold on

  % end for different solvers
end

% figure 1 processing
figure(1)
hold off
set(gca,'XScale','log','YScale','log','Fontsize',14);
set(gca,'xdir','reverse');
xlabel('H','Fontsize',18)
ylabel('Max Error','Fontsize',18)
legend(header1,'Location','bestoutside')

% figure 2 processing
figure(2)
hold off
set(gca,'XScale','log','YScale','log','Fontsize',14);
xlabel('Total function calls','Fontsize',18)
ylabel('Max Error','Fontsize',18)
legend(header2,'Location','BestOutside');

% figure 3 processing
figure(3)
hold off
set(gca,'XScale','log','YScale','log','Fontsize',14);
xlabel('Slow function calls','Fontsize',18)
ylabel('Max Error','Fontsize',18)
legend(header2,'Location','bestoutside');
