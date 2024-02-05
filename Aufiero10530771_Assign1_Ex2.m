% Spacecraft Guidance and Navigation (2023/2024)
% Assignment #1 - Guidance - Ex 2 Impulsive guidance
% Author: Aufiero Luca
% Personal code: 10530771
% Matricola: 232972
% Total expected execution time heavily depends on initial conditions,
%   for the example solution approximately 110-130 s

%% Start

clearvars; close all; clc; format long
cspice_furnsh('ex02.tm');
ticsolve = tic;

%% Ex 2.1

timewind = {'2029 JANUARY 01 00:00:00.000 TDB'; '2029 JULY 31 23:59:59.999 TDB'};
t0 = cspice_str2et(timewind(1));
tfinal = cspice_str2et(timewind(2));
time = linspace(t0,tfinal,5000);
rSUN = cspice_spkezr('20099942', time, 'ECLIPJ2000', 'NONE', 'SUN');
rEARTH = cspice_spkezr('20099942', time, 'ECLIPJ2000', 'NONE', 'EARTH');
pMOON = cspice_spkezr('MOON', time, 'ECLIPJ2000', 'NONE', 'SUN');
pEARTH = cspice_spkezr('399', time, 'ECLIPJ2000', 'NONE', 'SUN');
rMOON = cspice_spkezr('20099942', time, 'ECLIPJ2000', 'NONE', 'MOON');
date = cspice_et2utc(time, 'C', 1);

distSUN = zeros(1,length(time));
distEARTH = zeros(1,length(time));
distMOON = zeros(1,length(time));
distES = zeros(1,length(time));
Theta = zeros(1,length(time));

for i=1:length(time)
    distSUN(i) = norm(rSUN(1:3,i));
    distEARTH(i) = norm(rEARTH(1:3,i));
    distMOON(i) = norm(rMOON(1:3,i));
    distES(i) = norm(pEARTH(1:3,i));
    u = rSUN(1:3,i);
    v = rEARTH(1:3,i);
    CosTheta = dot(u,v)/(norm(u)*norm(v));
    Theta (i) = acosd(CosTheta);
end

[minDist, timeindex_nom] = min(distEARTH);
MinDate = cspice_et2utc(time(timeindex_nom), 'C', 1);

[minDistM, timeindex_nomM] = min(distMOON);
MinDateM = cspice_et2utc(time(timeindex_nomM), 'C', 1);

figure('Name','Ex 2.1.1')
plot(time,distSUN,'color','r','LineWidth',2)
hold on
grid on
plot(time,distES,'--','color','b')
hold off
ylabel('distance [km]','FontSize',15,'FontWeight','bold')
ylim([min(distSUN)-5e6,max(distSUN)+5e6])
ax=gca; ax.FontSize = 15; ax.FontWeight = 'bold';
xticks([min(time) time(timeindex_nom) max(time)])
xticklabels(["2029-01-01","TCA: 2029-04-13","2029-07-31"]);
legend('Distance between Apophis and the Sun','Distance between the Earth and the Sun',...
    'FontSize',15,'FontWeight','bold','Location','northeast')

figure('Name','Ex 2.1.2')
plot(time,distMOON,'color','k','LineWidth',2)
grid on
ylabel('distance [km]','FontSize',15,'FontWeight','bold')
ylim([min(distMOON)-1e7,max(distMOON)+1e7])
ax=gca; ax.FontSize = 15; ax.FontWeight = 'bold';
xticks([min(time) time(timeindex_nom) max(time)])
xticklabels(["2029-01-01","TCA: 2029-04-13","2029-07-31"]);
legend('Distance between Apophis and the Moon','FontSize',15,'FontWeight','bold','Location','north')

figure('Name','Ex 2.1.3')
plot(time,distEARTH,'LineWidth',2)
grid on
ylabel('distance [km]','FontSize',15,'FontWeight','bold')
ylim([min(distEARTH)-1e7,max(distEARTH)+1e7])
ax=gca; ax.FontSize = 15; ax.FontWeight = 'bold';
xticks([min(time) time(timeindex_nom) max(time)])
xticklabels(["2029-01-01","TCA: 2029-04-13","2029-07-31"]);
legend('Distance between Apophis and the Earth','FontSize',15,'FontWeight','bold','Location','north')

figure('Name','Ex 2.1.4')
plot(time,Theta,'color','g','LineWidth',2)
grid on
ylabel('angle [deg]','FontSize',15,'FontWeight','bold')
ax=gca; ax.FontSize = 15; ax.FontWeight = 'bold';
xticks([min(time) time(timeindex_nom) max(time)])
xticklabels(["2029-01-01","TCA: 2029-04-13","2029-07-31"]);
legend('angle Earth-Apophis-Sun','FontSize',15,'FontWeight','bold','Location','northeast')

figure(); hold on; grid on; xlabel('X [km]'); ylabel('Y [km]'); zlabel('Z [km]')
plot3(rSUN(1,:),rSUN(2,:),rSUN(3,:))
plot3(pEARTH(1,:),pEARTH(2,:),pEARTH(3,:))
plot3(pMOON(1,:),pMOON(2,:),pMOON(3,:))
plot3(rSUN(1,timeindex_nom),rSUN(2,timeindex_nom),rSUN(3,timeindex_nom),'.','MarkerSize',30,'color','c')
plot3(pEARTH(1,timeindex_nom),pEARTH(2,timeindex_nom),pEARTH(3,timeindex_nom),'.','MarkerSize',60,'color','b')
plot3(rSUN(1,timeindex_nomM),rSUN(2,timeindex_nomM),rSUN(3,timeindex_nomM),'.','MarkerSize',30)
plot3(pMOON(1,timeindex_nomM),pMOON(2,timeindex_nomM),pMOON(3,timeindex_nomM),'.','MarkerSize',45,'color','r')
plot3(0,0,0,'y.','MarkerSize',100)
title('Earth, Moon and Apophis orbits','FontSize',15,'FontWeight','bold')
legend('Apophis orbit','Earth orbit','Moon orbit','Apophis at CA','Earth at CA',...
    'Apophis at CA with Moon','Moon at CA with Apophis','Sun','FontSize',15,'FontWeight','bold')

% TCA ground-track
eapo_KM = sqrt(sum(distEARTH.^2,1));
[~,i_TCA] = min(eapo_KM);
Tw_TCA = linspace(time(i_TCA)-6/24*cspice_spd, time(i_TCA)+6/24*cspice_spd,3000);
rr_TCA = cspice_spkpos('20099942',Tw_TCA,'IAU_EARTH','NONE','EARTH');
[~,lonAPO,latAPO] = cspice_reclat(rr_TCA);
lonAPO = 180/pi*lonAPO;
latAPO = 180/pi*latAPO;
figure('Name','Apophis ground-track')
title('Ground Track','FontSize',15,'FontWeight','bold')
planisphere
gt = plotGroundTrack(lonAPO,latAPO,[0.2745 0.8784 0.4745],'g',1);

%% EX 2.3

% Time strings
LWO = '2024-10-01 00:00:00 TDB';
LWC = '2025-02-01 00:00:00 TDB';

t_open_DSM = '2025-04-01 00:00:00 TDB';         % LWO + 6 Months
t_close_DSM = '2026-08-01 00:00:00 TDB';        % LWC + 18 Months

t_open_impact = '2028-08-01 00:00:00 TDB';
t_close_impact = '2029-02-28 00:00:00 TDB';

% Time conversions in et
t_LWO = (cspice_str2et(LWO));
t_LWC = (cspice_str2et(LWC));
t_O_DSM = (cspice_str2et(t_open_DSM));
t_C_DSM = (cspice_str2et(t_close_DSM));
t_O_I = (cspice_str2et(t_open_impact));
t_C_I = (cspice_str2et(t_close_impact));

% Frame and center for Spice
frame = 'ECLIPJ2000';
center = 'SSB';

% Time window for the close encounter
et0 = cspice_str2et('TDB 2029/01/01-00:00:00.00');
etf = cspice_str2et('TDB 2029/07/31-00:00:00.00');
dim = 5000;
t_span = linspace(et0,etf,dim);

% Initializations: Distances SSB-Body
rr_Apo = zeros(3,length(t_span));
rr_Sun = zeros(3,length(t_span));
rr_Earth = zeros(3,length(t_span));
rr_Moon = zeros(3,length(t_span));

% Distances Obj-Body
rr_Apo_Sun = zeros(3,dim);
rr_Apo_Earth = zeros(3,dim);

% Compute relative distances
rr_apo_moon = rr_Moon - rr_Apo;

dist_apo_sun = sqrt(sum(rr_Apo_Sun.^2,1));
dist_apo_earth = sqrt(sum(rr_Apo_Earth.^2,1));

AU = cspice_convrt(1,'AU','KM'); 
earth_radii = cspice_bodvrd('399','RADII',3);

% Lower bound for Position Vectors
lb = -40 * AU *ones(21,1);

% Lower bounds for Times
lb(1) = t_LWO;
lb(2) = t_O_DSM;
lb(3) = t_O_I;

% Lower bound for Velocities
lb(7:9) = -100;
lb(13:15) = -100;
lb(19:21) = -100;

% Upper bound for Position Vectors
ub = 40 * AU  * ones(21,1);

% Upper bound for Times
ub(1) = t_LWC;
ub(2) = t_C_DSM;
ub(3) = t_C_I;

% Upper bound for Velocities
ub(7:9) = 100;
ub(13:15) = 100;
ub(19:21) = 100;

% Initialize parameters for iterations
exflag = 0;
iter = 0;
iter_max = 50;
Re_CA = 0;

% Iterations to find a better minimum
while(Re_CA < 16)

exflag = 0;
% Iterations to find a minimum
while (exflag~=1 && iter<iter_max )

% Guess for Times
% t1 = t_LWO + (t_LWC - t_LWO)/2;
% t2 = t_O_DSM + (t_C_DSM - t_O_DSM)/2;
% t3 = t_O_I + (t_C_I - t_O_I)/2;

% Random Times
t1 = t_LWO +  rand(1) * (t_LWC - t_LWO);
t2 = t_O_DSM +  rand(1) * (t_C_DSM - t_O_DSM);
t3 = t_O_I +  rand(1) * (t_C_I - t_O_I) ;

% Random Delta v
Dv1 = rand(3,1);
Dv1 = 2 * Dv1./norm(Dv1);
Dv2 = rand(3,1);
Dv2 = 2 * Dv2./norm(Dv2);

% Guess for states 
rve = cspice_spkezr('Earth',t1,frame,'NONE',center);
x1 = [rve(1:3) ; rve(4:6) + Dv1];
[x1f, ~, ~] = keplerian_propagator(t1 ,x1 , t2, 'Sun');
x2 = [x1f(1:3) ; x1f(4:6) + Dv2];
x3 = cspice_spkezr('20099942',t3,frame,'NONE',center);

% Random initial guess for 21 variables problem
% guess = [t1; t2; t3; x1; x2; x3];

% Example solution guess
guess = [781572608.989753
827784307.513450
906878954.791016
143890498.720087
36375361.9959268
25719.7664184030
-7.14430136920452
29.6805521756719
1.61098964897752
-144237294.855872
73225319.8450688
5590789.68918057
-12.4989167052946
-22.9118929883880
-0.551506863769696
113070286.591064
30679579.4769331
1019606.72977353
-5.35923358001503
35.6865075729688
-2.03150293249896]; % 24 Re

options = optimoptions('fmincon','Display','iter-detailed','StepTolerance',1e-9,'ConstraintTolerance',1e-9);

% Find the minimum
[optimal,min_CA,exflag,output] = fmincon(@objective_fun,guess,[],[],[],[],lb,ub,@constraints, options)

% Update Iteration
iter = iter+1;

t1_opt = optimal(1); 
t2_opt = optimal(2);
t3_opt = optimal(3); 
x1_opt = optimal(4:9);
x2_opt = optimal(10:15);
x3_opt = optimal(16:21);

max_CA = -min_CA;
Re_CA = max_CA/earth_radii(1);

end
end

%% Verifications of results

% Define list of celestial bodies:
labels = {'Sun';
          'Mercury';
          'Venus';
          'Earth';
          'Moon';
          'Mars Barycenter';
          'Jupiter Barycenter';
          'Saturn Barycenter';
          'Uranus Barycenter';
          'Neptune Barycenter';
          'Pluto Barycenter'};

% Initialize propagation data (same as regular n-body)
bodies = nbody_init(labels);

rv0 = cspice_spkezr('Earth',t1_opt,frame,'NONE',center);
r0 = rv0(1:3);
v0 = rv0(4:6);

% First segment propagation
[xf1, ~, xx1] = keplerian_propagator(t1_opt, x1_opt, t2_opt, 'Sun');

% Dv1 to leave Earth Orbit
Dv1_opt = x1_opt(4:6) - v0;

% Second segment propagation
[xf2, ~, xx2] = keplerian_propagator(t2_opt, x2_opt, t3_opt, 'Sun');

% Deep Space Maneuver
Dv2_opt = x2_opt(4:6) - xf1(4:6);

% Total Delta v
Dvtot = norm(Dv1_opt) + norm(Dv2_opt);

% Integration of Apophis dynamics using a full n-body propagator
options = odeset('reltol', 1e-12, 'abstol', [1e-6*ones(1,3),1e-9*ones(1,3)],'Events',@(t,x) close_approach(t,x,frame,center,1));
[tt, xx_Apo] = ode78(@(t,x) nbody_rhs(t,x,bodies,frame), [t3_opt t3_opt*3 ], x3_opt, options);

% Time of CA
t_CA = tt(end);

% Earth position at CA
r_E_f = cspice_spkpos('Earth',t_CA,frame,'NONE',center);

% Function to minimize
d_apo_E = -norm(xx_Apo(end,1:3)' - r_E_f);

% Print the results
fprintf('\nLaunch date: %s UTC\n', cspice_et2utc(t1_opt,'C',3))
fprintf('DSM date: %s UTC\n', cspice_et2utc(t2_opt,'C',3))
fprintf('Impact date: %s UTC\n', cspice_et2utc(t3_opt,'C',3))
fprintf('Time of Closest Approach: %s UTC\n',cspice_et2utc(t_CA,'C',3))
fprintf('Delta v Launch [km/s]: %15.10f %15.10f %15.10f\n', Dv1_opt(1), Dv1_opt(2), Dv1_opt(3))
fprintf('Delta v DSM [km/s]: %15.10f %15.10f %15.10f\n', Dv2_opt(1), Dv2_opt(2), Dv2_opt(3))
fprintf('Total Delta v [km/s] (max 5 km/s): %15.10f\n', Dvtot)
fprintf('Distance of nominal Closest Approach in Earth Radii: %15.10f\n', min(distEARTH)/earth_radii(1))
fprintf('Distance of Closest Approach after impact in Earth Radii: %15.10f\n', Re_CA)

%% Retrieve Earth and Apophis states

t01 = linspace(t_LWO,t_CA,dim);
t = linspace(t_O_I,time(timeindex_nom),dim);

% Initialization
r_E = zeros(3,length(tt));
r_apo = zeros(3,dim);
r_earth = zeros(3,dim);
r_e_o = zeros(3,dim);

for i = 1 : dim

    % Earth orbit 
    r_earth(:,i) = cspice_spkpos('Earth',t01(i),frame,'NONE',center);
    
    % Earth orbit to retreive old position at CA
    r_e_o(:,i) = cspice_spkpos('Earth',t(i),frame,'NONE',center);

    % Apophis states before the impact
    r_apo(:,i) = cspice_spkpos('20099942',t(i),frame,'NONE',center);

end

%% Plots

% Multiple Shooting 
MS = figure; hold on; grid on;
plot3(0,0,0,'y.','MarkerSize',100)
plot3(r_earth(1,:),r_earth(2,:),r_earth(3,:),'--','Color','b')
plot3(x1_opt(1),x1_opt(2),x1_opt(3),'k.','MarkerSize',50)
quiver3(x1_opt(1),x1_opt(2),x1_opt(3),Dv1_opt(1)*2e7,Dv1_opt(2)*2e7,Dv1_opt(3)*2e7,'LineWidth',2)
plot3(xx1(:,1),xx1(:,2),xx1(:,3),'LineWidth',1)
plot3(xf1(1),xf1(2),xf1(3),'.','MarkerSize',50)
quiver3(xf1(1),xf1(2),xf1(3),Dv2_opt(1)*2e7,Dv2_opt(2)*2e7,Dv2_opt(3)*2e7,'LineWidth',2)
plot3(xx2(:,1),xx2(:,2),xx2(:,3),'LineWidth',1)
plot3(xf2(1),xf2(2),xf2(3),'.','MarkerSize',50)
plot3(r_apo(1,:),r_apo(2,:),r_apo(3,:),'--','Color','r')
plot3(xx_Apo(:,1),xx_Apo(:,2),xx_Apo(:,3),'m','LineWidth',3,'Color','r')
plot3(r_apo(1,end),r_apo(2,end),r_apo(3,end),'rs','MarkerSize',15)
plot3(xx_Apo(end,1),xx_Apo(end,2),xx_Apo(end,3),'r.','MarkerSize',50)
plot3(r_e_o(1,end),r_e_o(2,end),r_e_o(3,end),'bs','MarkerSize',15)
plot3(r_earth(1,end),r_earth(2,end),r_earth(3,end),'b.','MarkerSize',50)
legend('Sun','Earth orbit','Launch position','Launch Delta v','S/C journey first leg','DSM position',...
        'DSM Delta v','S/C journey second leg','Impact position','Apophis nominal orbit',...
        'Apophis orbit after the impact','Apophis at original CA','Apophis at CA after impact',...
        'Earth at original CA','Earth at CA after impact',...
        'FontSize',15,'FontWeight','bold')
xlabel('X [km]'); ylabel('Y [km]'); zlabel('Z [km]')

fprintf('\nEx 2 execution time = %.2f s\n',toc(ticsolve))

%% FUNCTIONS

function [xf, tt, xx] = keplerian_propagator(et0, x0, etf, attractor)
% KEPLERIAN_PROPAGATOR Propagate a Keplerian Orbit 

% Initialize propagation data
if isfloat(attractor)
    GM = attractor;
else
    GM = cspice_bodvrd(attractor, 'GM', 1);
end

tof = etf-et0;

options = odeset('reltol', 1e-12, 'abstol', [ones(3,1)*1e-8; ones(3,1)*1e-11]);

% Perform integration
[tt, xx] = ode78(@(t,x) keplerian_rhs(t,x,GM), [0 tof], x0, options);

% Extract state vector 
xf = xx(end,1:6)';

end

%----------------------------------------------------------------------

function [dxdt] = keplerian_rhs(~, x, GM)
% KEPLERIAN_RHS  Evaluates the right-hand-side of a 2-body (keplerian) propagator
%   Evaluates the right-hand-side of a newtonian 2-body propagator.
%
%
% Author
%   Name: ALESSANDRO 
%   Surname: MORSELLI
%   Research group: DART
%   Department: DAER
%   University: Politecnico di Milano 
%   Creation: 24/10/2021
%   Contact: alessandro.morselli@polimi.it
%   Copyright: (c) 2021 A. Morselli, Politecnico di Milano. 
%                  All rights reserved.
%
%
% Notes:
%   This material was prepared to support the course 'Satellite Guidance
%   and Navigation', AY 2021/2022.
%
%
% Inputs:
%   t   : [1, 1] epoch (unused)
%   x   : [6, 1] cartesian state vector wrt Solar-System-Barycentre and
%                 State Transition Matrix elements
%   GM  : [1, 1] gravitational constant of the body
%
% Outputs:
%   dxdt   : [6,1] RHS, newtonian gravitational acceleration only
%

% Initialize right-hand-side
dxdt = zeros(6,1);

% Extract positions
rr = x(1:3);

% Compute square distance and distance
dist2 = dot(rr, rr);
dist = sqrt(dist2);

% Position detivative is object's velocity
dxdt(1:3) = x(4:6);   
% Compute the gravitational acceleration using Newton's law
dxdt(4:6) = - GM * rr /(dist*dist2);

end

%----------------------------------------------------------------------

function d_apo_E = objective_fun(y)
% Objective function

   t3 = y(3);
   x3 = y(16:21);
   
% Define list of celestial bodies:
labels = {'Sun';
          'Mercury';
          'Venus';
          'Earth';
          'Moon';
          'Mars Barycenter';
          'Jupiter Barycenter';
          'Saturn Barycenter';
          'Uranus Barycenter';
          'Neptune Barycenter';
          'Pluto Barycenter'};

% Initialize propagation data (same as regular n-body)
bodies = nbody_init(labels);

% Select integration frame string (SPICE naming convention)
frame = 'ECLIPJ2000';
center = 'SSB';

% Integration of Apophis dynamics using a full n-body propagator
options = odeset('reltol', 1e-12, 'abstol', [1e-8*ones(1,3),1e-9*ones(1,3)],'Events',@(t,x) close_approach(t,x,frame,center,1));
[tt, xx_Apo] = ode78(@(t,x) nbody_rhs(t,x,bodies,frame), [t3 t3*3], x3, options);

% Earth position at CA
r_E_f = cspice_spkpos('Earth',tt(end),frame,'NONE',center);

% Function to minimize
d_apo_E = -norm(xx_Apo(end,1:3)'- r_E_f);

end

%----------------------------------------------------------------------

function [bodies] = nbody_init(labels)
% NBODY_INIT Initialize planetary data for n-body propagation
%   Given a set of labels of planets and/or barycentres, returns a
%   cell array populated with structures containing the body label and the
%   associated gravitational constant.
% Inputs:
%   labels : [1,n] cell-array with object labels
%
% Outputs:
%   bodies : [1,n] cell-array with struct elements containing the following
%                  fields
%                  |
%                  |--bodies{i}.name -> body label
%                  |--bodies{i}.GM   -> gravitational constant [km**3/s**2]

% Initialize cell array bodies
bodies = cell(size(labels));

% Loop over labels
for i = 1 : length(labels)

    % Store body label
    bodies{i}.name = labels{i};
   
    % Store body gravitational constant
    bodies{i}.GM = cspice_bodvrd(labels{i},'GM',1);
end

end

%----------------------------------------------------------------------

function [dxdt] = nbody_rhs(t, x, bodies, frame)
% NBODY_RHS Evaluates the right-hand-side of a N-body propagator
%   Evaluates the right-hand-side of a newtonian N-body propagator.
%   The integration centre is the Solar-System-Barycentre (SSB) and only
%   Newtonian gravitational accelerations are considered.
%
% Inputs:
%   t      : [1,1] ephemeris time (ET SPICE), seconds past J2000 (TDB)
%   x      : [6,1] cartesian state vector wrt Solar-System-Barycentre
%   bodies : [1,6] cell-array created with function nbody_init
%
% Outputs:
%   dxdt   : [6,1] RHS, newtonian gravitational acceleration only

% Initialize right-hand-side
dxdt = zeros(6,1);

% Position detivative is object's velocity
dxdt(1:3) = x(4:6);

% Extract the object position from state x
rr_ssb_obj = x(1:3);

% Loop over all bodies
for i = 1 : length(bodies)
    % Retrieve position and velocity of i-th celestial body wrt Solar
    % System Barycentre in inertial frame
    rv_ssb_body = cspice_spkezr(bodies{i}.name,t,frame,'NONE','SSB');

    % Extract object position wrt. i-th celestial body
    rr_body_obj = rr_ssb_obj - rv_ssb_body(1:3);

    % Compute square distance and distance
    dist2 = dot(rr_body_obj,rr_body_obj);
    dist = sqrt(dist2);

    % Compute the gravitational acceleration using Newton's law
    aa_grav = - bodies{i}.GM * rr_body_obj/(dist * dist2);

    % Sum up acceleration to right-hand-side
    dxdt(4:6) = dxdt(4:6) + aa_grav;

end

end

%----------------------------------------------------------------------

function [value, isterminal, direction] = close_approach(t,x,frame,center,isTerminal)
% Event Function
% This Function stops the integration at TCA

    rv_Earth = cspice_spkezr('Earth',t,frame,'NONE',center);
    dist = x(1:3) - rv_Earth(1:3);
    dot_dist = x(4:6) - rv_Earth(4:6);

    %derivative of Rho square
    value = 2 * dot(dist, dot_dist);

    isterminal = isTerminal;
    direction = 1;
end

%----------------------------------------------------------------------

function [c,ceq] = constraints(y)
% Equality and inequality constraints
% This function builds the constraints useful to fmincon function
% in order to optimize the problem

   t1 = y(1);
   t2 = y(2);
   t3 = y(3);
   x1 = y(4:9);
   x2 = y(10:15);
   x3 = y(16:21);
  
frame = 'ECLIPJ2000';
center = 'SSB';

rv0 = cspice_spkezr('Earth',t1,frame,'NONE',center);
r0 = rv0(1:3);
v0 = rv0(4:6);

% First segment propagation
[xf1, ~, ~] = keplerian_propagator(t1, x1, t2, 'Sun');

% Second segment propagation
[xf2, ~, ~] = keplerian_propagator(t2, x2, t3, 'Sun');

% Delta v
Dv1 = x1(4:6) - v0;
Dv2 = x2(4:6) - xf1(4:6);

ceq = zeros(15,1);
Dv_max = 5;

% State of Apophis before impact
rv3 = cspice_spkezr('20099942',t3,frame,'NONE',center);

% Inequality constraints
c(1) = norm(Dv1) + norm(Dv2) - Dv_max;

% Equality constraints
ceq(1:3) = x1(1:3) - r0;
ceq(4:6) = x2(1:3) - xf1(1:3);
ceq(7:9) = x3(1:3) - xf2(1:3); 
ceq(10:12) = x3(1:3) - rv3(1:3);
ceq(13:15) = x3(4:6) - xf2(4:6) * 5e-5 - rv3(4:6);
  
end

%----------------------------------------------------------------------

function p = plotGroundTrack(lon,lat,color1,color2,i)

% Plots the ground track corresponding to the given latitudes and
% longitudes. It can optionally check whether the orbit is retrograde or
% prograde to display the ground track correctly.
% 
% PROTOTYPE:
%   p = plotGroundTrack(lon,lat,color,i)
% 
% INPUT:
%   lon[n] = longitude vector [deg]
%   lat[n] = latitude vector [deg]
%   color1 = ground track color, can either be a character array or a RGB
%       decimal triplet.
%   color2 = Start/End markers color
%   i[1] = inclination of the orbit [rad]
% 
% OUTPUT:
%   p[3] = vector of initial point, final point and ground track proper
%       (line element)
% 
% -------------------------------------------------------------------------

hold on

if nargin < 4
    i = 1; % if not specified, assume prograde orbit
    if nargin < 3
        color1 = 'g'; % if not specified, color is green
        color2 = 'g'; % if not specified, color is green
    end
end

ps = plot(lon(1),lat(1),'Marker','o','Color',color2,'MarkerSize',8,'LineStyle','none','LineWidth',3,'DisplayName','Start');
pe = plot(lon(end),lat(end),'Marker','square','Color',color2,'MarkerSize',10,'LineStyle','none','LineWidth',3,'DisplayName','End');
t1 = text(lon(1),lat(1),'start','FontSize',18,'FontWeight','bold','Color',color2,'VerticalAlignment','top','HorizontalAlignment','left');
t2 = text(lon(end),lat(end)+2,'end','FontSize',18,'FontWeight','bold','Color',color2,'VerticalAlignment','bottom','HorizontalAlignment','right');

for j = 1:length(lon)-1
    if i < pi/2
        if lon(j+1)*lon(j) < 0 && lon(j+1) < 0
            lon(j) = NaN;
        end
    else
        if lon(j+1)*lon(j) < 0 && lon(j+1) > 0
            lon(j) = NaN;
        end
    end
end

pp = plot(lon,lat,'.','Color',color1,'DisplayName','Ground track','LineWidth',1.5);
legend([pp,ps,pe],'Location','northoutside','NumColumns',3,'FontSize',15)
uistack([ps,pe,t1,t2],'top')
hold off

p = [pp ps pe t1 t2];
end

%----------------------------------------------------------------------

function planisphere()

% Represents a map of the surface of Earth.

I = imread("EarthTex.jpg");
I = flip(I,1);

image([-180 180],[-90 90],I);
grid on
xticks(-180:30:180)
yticks(-90:30:90)
set(gca,'YDir','normal')
hold on
axis equal tight
xlabel('Longitude [deg]','FontSize',15)
ylabel('Latitude [deg]','FontSize',15)
end
