%% ------------------------------------------------------------------------
%   
%   Supplementary Material for "Safe-by-Design: Approximate Nonlinear Model 
%   Predictive Control with Realtime Feasibility" by 
%   Jan Olucak, Arthur Castello B. de Oliveira, and Torbjørn Cunis
%
%   Short Description: Evalaute all methods and compare in several
%   categories.
%
%
%   Needed software: - CasADi 3.6 (for CaSoS)
%                    - CaSoS
%
%
%   License: see License file in repository.   
%
% ------------------------------------------------------------------------


close all
clc
clear

% generate true to generate tikz code using matlab2tikz (matlab2tikz must
% be on MATLAB path)
genTikz = true;

% load some conversion function for attitute parameterization
addpath("helperFunc\")

% load infinitesimal MPC data
load infiniteMPC_comparison_3_runs.mat

% load full horizon MPC data (IPOPT)
load fullHor_MPC_comparison_3Runs.mat

% load full horizon MPC RTI data
load fullHor_MPC_RTI_comparison_3Runs.mat

% load data from polynomial control law
load polyControl_comparison_3_runs.mat

% load data from CBF-CLF-QP law
load CBF_CLF_QP_comparison_3_runs.mat

% load data for horizon-one
load horOneMPC_comparison_3_runs.mat

% put all methods into one cell
x_sol_all = {x_sol_vec_infMPC, x_sol_vec_fullMPC, x_sol_vec_fullMPC_RTI, x_sol_vec_poly, x_sol_vec_CBF_CLF,x_sol_vechorOneMPC};
u_sol_all = {u_sol_vec_infMPC, u_sol_vec_fullMPC, u_sol_vec_fullMPC_RTI, u_sol_vec_poly, u_sol_vec_CBF_CLF,u_sol_vechorOneMPC};
iter_all  = {iter_conv_inf, iter_all_fullMPC, iter_all_fullMPC_RTI, iter_conv_poly_inf, iter_conv_CBF_CLF,iter_convhorOne};


numRuns = length(x_sol_all);

t = linspace(0, simTime, simTime/simStepSize);

for jj = 1:3

% fprintf('------------------------------------------- \n')
% fprintf('------------------------------------------- \n')
% fprintf('Run number: %d \n',jj)
% fprintf('------------------------------------------- \n')
%% convergence time
convergenceTime_infinite(jj)    = t(iter_conv_inf{jj});
convergenceTime_full(jj)        = t(iter_all_fullMPC{jj});
convergenceTime_full_RTI(jj)    = t(iter_all_fullMPC_RTI{jj});
convergenceTime_poly(jj)        = t(iter_conv_poly_inf{jj});
convergenceTime_CBF_CLF_QP(jj)  = t(iter_conv_CBF_CLF{jj});
convergenceTime_horOne(jj)          = t(iter_convhorOne{jj});

% fprintf('------------------------------------------- \n')

% fprintf('Convergence time inf.MPC: %f s \n',convergenceTime_infinite(jj))
% fprintf('Convergence time NMPC (IPOPT): %f s \n',convergenceTime_full(jj))
% fprintf('Convergence time RTI: %f s \n',convergenceTime_full_RTI(jj))
% fprintf('Convergence time Poly.Law: %f s \n',convergenceTime_poly(jj))
% fprintf('Convergence time CBF-CLF-QP: %f s \n',convergenceTime_CBF_CLF_QP(jj))
% 
% fprintf('------------------------------------------- \n')

% %% worst-case computation time (in seconds)
% fprintf('Mean (worst-case) comp. time inf.MPC: %f (%f) ms \n',  [meanSolveTime_infMPC, maxSolveTime_infMPC])
% fprintf('Mean (worst-case) comp. time NMPC (IPOPT): %f (%f) ms \n',[meanSolveTime_full, maxSolveTime_full])
% fprintf('Mean (worst-case) comp. time RTI: %f (%f) ms \n',[meanSolveTime_full_RTI, maxSolveTime_full_RTI])
% fprintf('Mean (worst-case) comp. time Poly.Law: %f (%f) ms \n',[meanSolveTime_poly, maxSolveTime_poly])
% fprintf('Mean (worst-case) comp. time CBF-CLF-QP: %f (%f) ms \n',[meanSolveTime_CBF_CLF ,maxSolveTime_CBF_CLF])


%% evaluate performance

% -------- inf. MPC ----------

% We have to bring back MRP
D = diag([pi/180,pi/180,pi/180]);
x_sol_vec_infMPC_MRP =  x_sol_vec_infMPC{jj};

for j = 1:size(x_sol_vec_infMPC_MRP,2)
x_sol_vec_infMPC_MRP(4:6,j) = Euler1232MRP(D*x_sol_vec_infMPC_MRP(4:6,j));
end

% Define L as an anonymous function 
L = @(x, u)   u' * R_infinite * u + x'*Q_infinite*x;

% Preallocate cost array for each time step
cost_per_step = zeros(1, iter_conv_inf{jj});

% Compute L for each time step individually
for k = 1:iter_conv_inf{jj}
    cost_per_step(k) = L(x_sol_vec_infMPC_MRP(:,k), u_sol_vec_infMPC{jj}(:,k));
end

% Compute total cost using trapezoidal integration
stageCost_inf(jj) = trapz(linspace( 0,t(iter_conv_inf{jj}),length(cost_per_step)), cost_per_step);


% ---------- full-horizon (IPOPT) ----------
% to MRP
x_sol_vec_fullMPC_MRP =  x_sol_vec_fullMPC{jj};

% there was a slight mistake in the initialization of the array to store
% data not in the simulation! Array was initialized with nan values but
% only attitude was correctly overwritten!
x_sol_vec_fullMPC_MRP(1:3,1) = zeros(3,1);  

for j = 1:size(x_sol_vec_fullMPC{jj},2)
    x_sol_vec_fullMPC_MRP(4:6,j) = Euler1232MRP(D*x_sol_vec_fullMPC_MRP(4:6,j));
end


% Preallocate cost array for each time step
cost_per_step = zeros(1, iter_all_fullMPC{jj});

% Compute L for each time step individually
for k = 1:iter_all_fullMPC{jj}-1
    cost_per_step(k) = L(x_sol_vec_fullMPC_MRP(:,k), u_sol_vec_fullMPC{jj}(:,k));
end

% Compute total cost using trapezoidal integration
stageCost_full(jj) = trapz(linspace( 0,t(iter_all_fullMPC{jj}),length(cost_per_step)),cost_per_step);


% ---------- full-horizon RTI ----------

% to MRP
x_sol_vec_RTIMPC_MRP =  x_sol_vec_fullMPC_RTI{jj};


% there was a slight mistake in the initialization of the array to store
% data not in the simulation! Array was initialized with nan values but
% only attitude was correctly overwritten!
x_sol_vec_RTIMPC_MRP(1:3,1) = zeros(3,1);  

for j = 1:size(x_sol_vec_fullMPC_RTI{jj},2)
    x_sol_vec_RTIMPC_MRP(4:6,j) = Euler1232MRP(D*x_sol_vec_RTIMPC_MRP(4:6,j));
end


% Preallocate cost array for each time step
cost_per_step = zeros(1, iter_all_fullMPC_RTI{jj});

% Compute L for each time step individually
for k = 1:iter_all_fullMPC_RTI{jj}-1
    cost_per_step(k) = L(x_sol_vec_RTIMPC_MRP(:,k), u_sol_vec_fullMPC_RTI{jj}(:,k));
end

% Compute total cost using trapezoidal integration
stageCost_full_RTI(jj) = trapz(linspace( 0,t(iter_all_fullMPC_RTI{jj}),length(cost_per_step)), cost_per_step);


% ---------- Poly.Law ----------
% We have to bring back MRP
D = diag([pi/180,pi/180,pi/180]);
x_sol_vec_poly_MRP =  x_sol_vec_poly{jj};

for j = 1:size(x_sol_vec_poly_MRP,2)
    x_sol_vec_poly_MRP(4:6,j) = Euler1232MRP(D*x_sol_vec_poly_MRP(4:6,j));
end


% Preallocate cost array for each time step
cost_per_step = zeros(1, iter_conv_poly_inf{jj});

% Compute L for each time step individually
for k = 1:iter_conv_poly_inf{jj}-1
    cost_per_step(k) = L(x_sol_vec_poly_MRP(:,k), u_sol_vec_poly{jj}(:,k));
end

% Compute total cost using trapezoidal integration
stageCost_poly(jj) = trapz(linspace( 0,t(iter_conv_poly_inf{jj}),length(cost_per_step)), cost_per_step);


% ---------- CBF-CLF-QP ----------
% We have to bring back MRP
D = diag([pi/180,pi/180,pi/180]);
x_sol_vec_CBF_CLF_MRP =  x_sol_vec_CBF_CLF{jj};

for j = 1:size(x_sol_vec_CBF_CLF_MRP,2)
    x_sol_vec_CBF_CLF_MRP(4:6,j) = Euler1232MRP(D*x_sol_vec_CBF_CLF_MRP(4:6,j));
end


% Preallocate cost array for each time step
cost_per_step = zeros(1, iter_conv_CBF_CLF{jj});

% Compute L for each time step individually
for k = 1:iter_conv_CBF_CLF{jj}-1
    cost_per_step(k) = L(x_sol_vec_CBF_CLF_MRP(:,k), u_sol_vec_CBF_CLF{jj}(:,k));
end

% Compute total cost using trapezoidal integration
stageCost_CBF_CLF(jj) = trapz(linspace( 0,t(iter_conv_CBF_CLF{jj}),length(cost_per_step)), cost_per_step);





% ---------- Horizon-one ----------
% We have to bring back MRP
D = diag([pi/180,pi/180,pi/180]);
x_sol_vec_horOne_MRP =  x_sol_vechorOneMPC{jj};

for j = 1:size(x_sol_vec_horOne_MRP,2)
    x_sol_vec_horOne_MRP(4:6,j) = Euler1232MRP(D*x_sol_vec_horOne_MRP(4:6,j));
end


% Preallocate cost array for each time step
cost_per_step = zeros(1, iter_convhorOne{jj});

% Compute L for each time step individually
for k = 1:iter_convhorOne{jj}-1
    cost_per_step(k) = L(x_sol_vec_horOne_MRP(:,k), u_sol_vechorOneMPC{jj}(:,k));
end

% Compute total cost using trapezoidal integration
stageCost_horOne(jj) = trapz(linspace( 0,t(iter_convhorOne{jj}),length(cost_per_step)), cost_per_step);


% fprintf('------------------------------------------- \n')
% fprintf('Int. stage cost inf.MPC: %f  \n',stageCost_inf(jj))
% fprintf('Int. stage cost NMPC (IPOPT): %f  \n',stageCost_full(jj))
% fprintf('Int. stage cost RTI: %f  \n',stageCost_full_RTI(jj))
% fprintf('Int. stage cost Poly.Law: %f  \n',stageCost_poly(jj))
% fprintf('Int. stage cost CBF-CLF_QP: %f  \n',stageCost_CBF_CLF(jj))


%% Plotting
colors = lines(numRuns);                        % Get a colormap for different runs
Euler_names = {'\phi','\theta','\psi'};

% re-scale rate constraints (real physical constraints)
x_low =  [-Omega_bounds(1)*180/pi -Omega_bounds(2)*180/pi -Omega_bounds(3)*180/pi]';
x_up  =  -x_low;

% set up with torques in miliNetwonmeter
ulow = u_low';
uup  = u_up';

% reduce number of data points
downsample_factor = 10;

%% summarized plot for paper

iter_all_numArr = zeros(3,numRuns);
for kk = 1:numRuns 
    iter_all_numArr(:,kk) =   cell2mat(iter_all{kk});
end


figure('Name',['Motion about x-axis for Run no. ' num2str(jj)])

% ------ Roll rate ------
subplot(3,1,1)
% zero line
plot([0 t(max([iter_all_numArr(jj,:)]))], [0 0],'k-.')
hold on;
i = 1;
for j = 1:numRuns
    plot(t(1:downsample_factor:iter_all{j}{jj}), x_sol_all{j}{jj}(i,1:downsample_factor:iter_all{j}{jj}) * 180/pi, 'Color', colors(j, :));
    hold on
end
xlabel('t [s]');
ylabel(sprintf('\\omega_%c [°/s]', 'x' + (i-1)));
grid on;
      

% plot gray shadded area and dashed gray line
xLimits = [0 t(max([iter_all_numArr(jj,:)]))]; 
yDashed = x_low(i); 
plot([0 t(max([iter_all_numArr(jj,:)]))],[yDashed yDashed],'--','Color',[0.5 0.5 0.5])
miny =  x_low(i)+0.5* x_low(i);
fill([xLimits(1) xLimits(2) xLimits(2) xLimits(1)], [miny miny yDashed yDashed], [0.7 0.7 0.7],'FaceAlpha',0.4 ,'EdgeColor', 'none');

xLimits = [0 t(max([iter_all_numArr(jj,:)]))]; 
yDashed =  x_up(i);
plot([0 t(max([iter_all_numArr(jj,:)]))],[yDashed yDashed],'--','Color',[0.5 0.5 0.5])
maxy = x_up(i)+0.5*x_up(i);
fill([xLimits(1) xLimits(2) xLimits(2) xLimits(1)], [yDashed yDashed maxy maxy], [0.7 0.7 0.7],'FaceAlpha',0.4 ,'EdgeColor', 'none');

axis([0 t(max([iter_all_numArr(jj,:)])) miny maxy])

% custom legend
h = zeros(6, 1);
h(1) = plot(NaN,NaN,'-','Color',colors(1,:));
h(2) = plot(NaN,NaN,'-','Color',colors(2,:));
h(3) = plot(NaN,NaN,'-','Color',colors(3,:));
h(4) = plot(NaN,NaN,'-','Color',colors(4,:));
h(5) = plot(NaN,NaN,'-','Color',colors(5,:));
h(6) = plot(NaN,NaN,'-','Color',colors(6,:));

legend(h, 'Inf.MPC','NMPC (IPOPT)','RTI','Poly. Law','CBF-CLF','horizon-one')



% ------ Roll angle ------
subplot(312)
% zero line
plot([0 t(max([iter_all_numArr(jj,:)]))], [0 0],'k-.')
hold on;
for j = 1:numRuns
    plot(t(1:downsample_factor:iter_all{j}{jj}), x_sol_all{j}{jj}(4,1:downsample_factor:iter_all{j}{jj}), 'Color', colors(j, :));
end
% axis([0 t(max([iter_all{:}])) -180 180])
xlabel('t [s]');
ylabel([Euler_names{1} ' [deg]']);
grid on;

% custom legend
h = zeros(6, 1);
h(1) = plot(NaN,NaN,'-','Color',colors(1,:));
h(2) = plot(NaN,NaN,'-','Color',colors(2,:));
h(3) = plot(NaN,NaN,'-','Color',colors(3,:));
h(4) = plot(NaN,NaN,'-','Color',colors(4,:));
h(5) = plot(NaN,NaN,'-','Color',colors(5,:));
h(6) = plot(NaN,NaN,'-','Color',colors(6,:));

legend(h, 'Inf.MPC','NMPC (IPOPT)','RTI','Poly. Law','CBF-CLF','horizon-one')



axis([0 t(max([iter_all_numArr(jj,:)])) -40 120])

% ------ torque ------
subplot(313)
% zero line
plot([0 t(max([iter_all_numArr(jj,:)]))], [0 0],'k-.')
hold on;
for j = 1:numRuns
    plot(t(1:downsample_factor:iter_all{j}{jj}-1), u_sol_all{j}{jj}(i,(1:downsample_factor:iter_all{j}{jj}-1)), 'Color', colors(j, :));
end
xlabel('t [s]');
ylabel(sprintf('\\tau_%c [Nm]', 'x' + (i-1)));
grid on;

% plot gray shadded area and dashed gray line
xLimits = [0 t(max([iter_all_numArr(jj,:)]))]; 
yDashed = ulow(i); 
miny =  ulow(i)+0.5* ulow(i);
plot([0 t(max([iter_all_numArr(jj,:)]))],[yDashed yDashed],'--','Color',[0.5 0.5 0.5])
fill([xLimits(1) xLimits(2) xLimits(2) xLimits(1)], [miny miny yDashed yDashed], [0.7 0.7 0.7],'FaceAlpha',0.4 ,'EdgeColor', 'none');

xLimits = [0 t(max([iter_all_numArr(jj,:)]))];  
yDashed =  uup(i); 
plot([0 t(max([iter_all_numArr(jj,:)]))],[yDashed yDashed],'--','Color',[0.5 0.5 0.5]) 
maxy = uup(i)+0.5*uup(i);
fill([xLimits(1) xLimits(2) xLimits(2) xLimits(1)], [yDashed yDashed maxy maxy], [0.7 0.7 0.7],'FaceAlpha',0.4 ,'EdgeColor', 'none');
axis([0 t(max([iter_all_numArr(jj,:)])) miny maxy])

% custom legend
% custom legend
h = zeros(6, 1);
h(1) = plot(NaN,NaN,'-','Color',colors(1,:));
h(2) = plot(NaN,NaN,'-','Color',colors(2,:));
h(3) = plot(NaN,NaN,'-','Color',colors(3,:));
h(4) = plot(NaN,NaN,'-','Color',colors(4,:));
h(5) = plot(NaN,NaN,'-','Color',colors(5,:));
h(6) = plot(NaN,NaN,'-','Color',colors(6,:));

legend(h, 'Inf.MPC','NMPC (IPOPT)','RTI','Poly. Law','CBF-CLF','horizon-one')


if genTikz && jj == 1
 cleanfigure();
 matlab2tikz('comp_I.tex','width','\figW','height','\figH');
elseif genTikz && jj == 2
 cleanfigure();
 matlab2tikz('comp_II.tex','width','\figW','height','\figH'); 
elseif genTikz && jj == 3
 cleanfigure();
 matlab2tikz('comp_III.tex','width','\figW','height','\figH');
end


end


%% Evaluate aggregated table
avg_conv_Time       = [mean(convergenceTime_infinite); mean(convergenceTime_full); mean(convergenceTime_full_RTI); mean(convergenceTime_poly); mean(convergenceTime_CBF_CLF_QP);mean(convergenceTime_horOne)];          
avg_Int_stage_cost  = [mean(stageCost_inf); mean(stageCost_full); mean(stageCost_full_RTI); mean(stageCost_poly); mean(stageCost_CBF_CLF);mean(stageCost_horOne)]; 
comp_Time_mean      = [meanSolveTime_infMPC; meanSolveTime_full; meanSolveTime_full_RTI; meanSolveTime_poly; meanSolveTime_CBF_CLF;meanSolveTimehorOneMPC]; 
comp_Time_worst     = [maxSolveTime_infMPC; maxSolveTime_full; maxSolveTime_full_RTI; maxSolveTime_poly; maxSolveTime_CBF_CLF;maxSolveTimehorOneMPC]; 
% Row names
rowNames = {'Inf.MPC','NMPC (IPOPT)','RTI','Poly. Law','CBF-CLF','horizon-one'};

% Create the table
T = table(avg_conv_Time, comp_Time_mean, comp_Time_worst, avg_Int_stage_cost, 'RowNames', rowNames);

% Display the table
disp(T);





