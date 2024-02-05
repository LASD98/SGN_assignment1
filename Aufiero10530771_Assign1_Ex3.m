% Spacecraft Guidance and Navigation (2023/2024)
% Assignment #1 - Guidance - Ex 3 Continuous guidance
% Author: Aufiero Luca
% Personal code: 10530771
% Matricola: 232972
% Total expected execution time: variable depending on tries needed to find
%   the 800 mN solution (approximately 20-30 s per 800 mN try + 60 s)

%% Load kernels

cspice_furnsh 'assignment1.tm';
clearvars; close all; clc; format long
warning('off','all')

%% Ex 3.1

ticsolve = tic;
TR = linspace(1,5/8,30); % Thrust ratio for numerical continuation
%TR = [1;15/16;14/16;13/16;12/16;11/16;10/16]; % Thrust ratio 
%TR = [1;7/8;6/8;5/8]; % Thrust ratio 
%TR = [1;5/8]; % Thrust ratio
%TR = 1; % Thrust ratio
m0 = 1000;  % [kg] s/c mass at t0

T0 = '2023-05-28-14:13:09.000 UTC';
t0 = cspice_str2et(T0); % [s]
x0Earth = cspice_spkezr('Earth',t0,'ECLIPJ2000','none','Sun'); % [km][Km/s]

km2au = cspice_convrt(1,'KM','AU');
sec2year = cspice_convrt(1,'SECONDS','YEARS');

x0Earth(1:3) = x0Earth(1:3)*km2au;
x0Earth(4:6) = x0Earth(4:6)*km2au/sec2year;
t0 = t0*sec2year;

x0 = [x0Earth;m0;t0]; % initial boundary conditions + t0

%% Ex 3.2

mu = cspice_bodvrd('Sun','GM',1); % [km^3/s^2]
Tmax = 800*1e-6; % [kN --> kg*km/s^2] Max Thrust
Isp = 3120; % [s] % Specific impulse
g0 = 9.80665*1e-3; % [km/s^2] gravitational acceleration at sea level
t0E = cspice_str2et(T0); % [s]
x0E = cspice_spkezr('Earth',t0E,'ECLIPJ2000','none','Sun'); % [km][km/s]

LU = cspice_convrt(1,'AU','KM'); % km
MU = m0; % kg
TU = sqrt(LU^3/mu); % s

mu_ad = 1;
m0_ad = m0 / MU;
g0_ad = g0 / (LU/TU^2);
Tmax_ad = Tmax / (MU*LU/TU^2);
Isp_ad = Isp / TU;
rr0_ad = x0E(1:3) / LU;
vv0_ad = x0E(4:6) / (LU/TU);

fprintf('\nADIMENSIONALIZED QUANTITIES\n')
fprintf('r0: %.10f %.10f %.10f\n', rr0_ad(1), rr0_ad(2), rr0_ad(3))
fprintf('v0: %.10f %.10f %.10f\n', vv0_ad(1), vv0_ad(2), vv0_ad(3))
fprintf('m0: %.10f\n', m0_ad)
fprintf('Isp: %.10f\n', Isp_ad)
fprintf('Tmax: %.10f\n', Tmax_ad)
fprintf('g0: %.10f\n', g0_ad)
fprintf('GM: %.10f\n', mu_ad)

%% Ex 3.3 + 3.4 FSOLVE

opt = optimoptions('fsolve', ...
                  'FunctionTolerance',1e-5, ...
                  'StepTolerance',1e-50, ...
                   'MaxFunctionEvaluations',2e3, ...
                   'Display', ...
                   'off');

e = 0;
disc = 0;

% x = a + (b-a)*rand, search x between interval (a,b)
% search Δt within interval (t1,t2)

x = -20 + (20-(-20)).*rand(7,1); % {λ0} 

t1 = 0;
t2 = 2*pi;

x(8) = t1/12 + (t2/12 - t1/12)*rand; % {Δtf}

k = 1;
n = 0;

%% 800 mN solution

ticTR = tic;
while 1
    tic
    n = n + 1;

    D = char(916);
    L = char(955);

    inguess = sprintf(['\nThrust: %d mN\n\nFIRST GUESS:\n'...  
             '%sr: %.10f %.10f %.10f; \n%sv: %.10f %.10f %.10f; \n' ...
             '%sm: %.10f; %st: %.10f[years]\n' ...
             'x = [%.32f %.32f %.32f %.32f %.32f %.32f %.32f %.32f]'';\n\n'],...
             TR(k)*800,L,x(1:3),L,x(4:6),L,x(7),D,x(8),x);
    
    fprintf('%s',inguess)

    [y,F,exitflag,output] = fsolve(@boundariesfun,x,opt,TR(k),x0);
    y(7) = abs(y(7));

    foo = output.firstorderopt;
    Tf = cspice_convrt(t0+y(8),'YEARS','SECONDS');
    Tf = cspice_timout(Tf,'YYYY-MM-DD-HR:MM:SC.### UTC');
    outguess = sprintf(['INITIAL CONDITION:\n'...  
             '%sr: %.10f %.10f %.10f; \n%sv: %.10f %.10f %.10f; \n' ...
             '%sm: %.10f; %st: %.10f [years]\nArrival Date: %s\n' ...
             'y = [%.32f %.32f %.32f %.32f %.32f %.32f %.32f %.32f]'';\n\n'],...
             L,y(1:3),L,y(4:6),L,y(7),D,y(8),Tf,y);
    
    fprintf('exit flag: %d\n\n',exitflag)
    
    fprintf('%s\n',outguess)

    % errors % [AU,YEARS]
    Er = abs(F(1:3)); % position error
    Ev = abs(F(4:6)); % velocity error
    er = norm(Er);
    ev = norm(Ev);
    r = norm(F)^2; 
    au2km = cspice_convrt(1,'AU','KM');
    year2sec = cspice_convrt(1,'YEARS','SECONDS');
    error = sprintf(['errors:\n' ...
             'er: %.10e[Km] %.10e[Km] %.10e[Km]\n', ...
             'ev: %.10e[Km/s] %.10e[Km/s] %.10e[Km/s]\n' ...
             '||er||: %.10e[Km]; ||ev||: %.10e[Km/s]\n||F||^2: %.10e\n' ...
             'First Order Optimality: %.10e\n'],Er*au2km,Ev*au2km/year2sec,er*au2km,ev*au2km/year2sec,r,foo);
    fprintf('%s',error)
    
    ANS(k).initialCondition.firstguess_x = x;
    ANS(k).initialCondition.y = y;
    ANS(k).errors.Er = Er;
    ANS(k).errors.Er_km = Er*au2km;
    ANS(k).errors.er = er;
    ANS(k).errors.Ev = Ev;
    year2sec = cspice_convrt(1,'YEARS','SECONDS');
    ANS(k).errors.Ev_kms = Ev*au2km/year2sec;
    ANS(k).errors.FirstOderOpt = foo;
    ANS(k).inguess = inguess;
    ANS(k).outguess = outguess;
    ANS(k).error = error;

    fprintf('Execution time %d mN Thrust try n. %d = %.2fs\n',TR(k)*800,n,toc)

    if exitflag == 1
        fprintf('\nExecution time %d mN Thrust = %.2fs\n\n\n',TR(k)*800,toc(ticTR))

        x = y;

        break
    elseif exitflag == 0 || exitflag == -2 || exitflag == -3
        e = e + 1; % number of mistakes

        x = -20 + (20-(-20)).*rand(7,1); % {λ0} 
        t1 = 0;
        t2 = 2*pi;
        x(8) = t1/12 + (t2/12 - t1/12)*rand; % {Δtf}

    elseif exitflag == 4 || exitflag == 3 || exitflag == 2
        disc = disc + 1; % number of discarded good guesses
    end
end

fprintf('\nExecution time fsolve (%d bad guess/es, %d good discarded) = %.2f s\n',e,disc,toc(ticsolve))


%% Numerical continuation to obtain 500 mN solution

for k = 2:length(TR)
    n = 0;

ticTR = tic;
    while 1
        tic
        n = n + 1;

        D = char(916);
        L = char(955);

        inguess = sprintf(['\nThrust: %d mN\n\nFIRST GUESS:\n'...  
                 '%sr: %.10f %.10f %.10f; \n%sv: %.10f %.10f %.10f; \n' ...
                 '%sm: %.10f; %st: %.10f[years]\n' ...
                 'x = [%.32f %.32f %.32f %.32f %.32f %.32f %.32f %.32f]'';\n\n'],...
                 TR(k)*800,L,x(1:3),L,x(4:6),L,x(7),D,x(8),x);
        
        fprintf('%s',inguess)

        [y,F,exitflag,output] = fsolve(@boundariesfun,x,opt,TR(k),x0);
        y(7) = abs(y(7));

        foo = output.firstorderopt;
        Tf = cspice_convrt(t0+y(8),'YEARS','SECONDS');
        Tf = cspice_timout(Tf,'YYYY-MM-DD-HR:MM:SC.### UTC');
        outguess = sprintf(['INITIAL CONDITION:\n'...  
                 '%sr: %.10f %.10f %.10f; \n%sv: %.10f %.10f %.10f; \n' ...
                 '%sm: %.10f; %st: %.10f [years]\nArrival Date: %s\n' ...
                 'y = [%.32f %.32f %.32f %.32f %.32f %.32f %.32f %.32f]'';\n\n'],...
                 L,y(1:3),L,y(4:6),L,y(7),D,y(8),Tf,y);
        
        fprintf('exit flag: %d\n\n',exitflag)
        
        fprintf('%s\n',outguess)

        % errors % [AU,YEARS]
        Er = abs(F(1:3)); % position error
        Ev = abs(F(4:6)); % velocity error
        er = norm(Er);
        ev = norm(Ev);
        r = norm(F)^2; 
        au2km = cspice_convrt(1,'AU','KM');
        year2sec = cspice_convrt(1,'YEARS','SECONDS');
        error = sprintf(['errors:\n' ...
                 'er: %.10e[Km] %.10e[Km] %.10e[Km]\n', ...
                 'ev: %.10e[Km/s] %.10e[Km/s] %.10e[Km/s]\n' ...
                 '||er||: %.10e[Km]; ||ev||: %.10e[Km/s]\n||F||^2: %.10e\n' ...
                 'First Order Optimality: %.10e\n'],Er*au2km,Ev*au2km/year2sec,er*au2km,ev*au2km/year2sec,r,foo);
        fprintf('%s',error)
        
        ANS(k).initialCondition.firstguess_x = x;
        ANS(k).initialCondition.y = y;
        ANS(k).errors.Er = Er;
        ANS(k).errors.Er_km = Er*au2km;
        ANS(k).errors.er = er;
        ANS(k).errors.Ev = Ev;
        year2sec = cspice_convrt(1,'YEARS','SECONDS');
        ANS(k).errors.Ev_kms = Ev*au2km/year2sec;
        ANS(k).errors.FirstOderOpt = foo;
        ANS(k).inguess = inguess;
        ANS(k).outguess = outguess;
        ANS(k).error = error;

        fprintf('Execution time %d mN Thrust try n. %d = %.2fs\n',TR(k)*800,n,toc)

        x = y;

        if exitflag == 1
            fprintf('\nExecution time %d mN Thrust = %.2fs\n\n\n',TR(k)*800,toc(ticTR))
            break
        elseif exitflag == 0 || exitflag == -2 || exitflag == -3
            e = e + 1; % number of mistakes
        elseif exitflag == 4 || exitflag == 3 || exitflag == 2
            disc = disc + 1; % number of discarded good guesses
        end
    end
end
fprintf('\nExecution time fsolve (%d bad guess/es, %d good discarded) = %.2f s\n',e,disc,toc(ticsolve))

%% ERROR TABLE

errortab = struct;
errorTab = table;
prec = '%.10e';
for k = 1:length(TR)
    errortab(k).TR = TR(k)*800;
    errortab(k).FOO = num2str(ANS(k).errors.FirstOderOpt,'%.3e');
    errortab(k).dt = num2str(ANS(k).initialCondition.y(8)*365,'%.2f');
    errortab(k).r = num2str(norm(ANS(k).errors.Er_km),prec);
    errortab(k).v = num2str(norm(ANS(k).errors.Ev_kms),prec);
    errortab(k).rx = num2str(ANS(k).errors.Er_km(1),prec);
    errortab(k).ry = num2str(ANS(k).errors.Er_km(2),prec);
    errortab(k).rz = num2str(ANS(k).errors.Er_km(3),prec);
    errortab(k).vx = num2str(ANS(k).errors.Ev_kms(1),prec);
    errortab(k).vy = num2str(ANS(k).errors.Ev_kms(2),prec);
    errortab(k).vz = num2str(ANS(k).errors.Ev_kms(3),prec);
end
errorTab = struct2table(errortab);
for i = 1 : width(errorTab)
    if height(errorTab) == 1
        break
    end
    errorTab.(i) = categorical(errorTab.(i));
end

errorTab.Properties.VariableNames = {'Thrust [mN]','First Order Optimality', ...
    [D,'t[days]'],['||',D,'r||[km]'],['||',D,'v||[km/s]'],[D,'r_x[km]'], ...
    [D,'r_y[km]'],[D,'r_z[km]'],[D,'v_x[km/s]'],[D,'v_y[km/s]'],[D,'v_z[km/s]']}

ANS(1).errortab = errorTab;
ANS = orderfields(ANS, [6,1:5]);

%% Propagation

km2au = cspice_convrt(1,'KM','AU');
sec2year = cspice_convrt(1,'SECONDS','YEARS');

for k = 1:length(TR)
    y = ANS(k).initialCondition.y;
    tf = cspice_convrt(t0+y(8),'YEARS','SECONDS');
    xM = cspice_spkezr('Venus',tf,'ECLIPJ2000','none','Sun');

    xM(1:3) = xM(1:3)*km2au;

    [~,out,t_out] = orbit_propagator([x0(1:7);y(1:7)],t0,t0+y(8),TR(k));
    
    % costates
    Lr = out(:,8:10);
    Lv = out(:,11:13);
    Lm = out(:,14);

    % norm of costates
    lr = zeros(length(t_out),1);
    lv = zeros(length(t_out),1);
    for j = 1:length(t_out)
        lr(j) = norm(Lr(j,:));
        lv(j) = norm(Lv(j,:));
    end

    % primer vector

    alpha = - Lv ./ lv;
    ANS(k).out = out;
    ANS(k).t_out = t_out;
    ANS(k).Lr = Lr;
    ANS(k).lr = lr;
    ANS(k).Lv = Lv;
    ANS(k).lv = lv;
    ANS(k).Lm = Lm;
    ANS(k).alpha = alpha;
    ANS(k).xM.xM = xM;
end
ANS = orderfields(ANS, [1:2,15,3:14]);

%% PLOT Orbit

% Environment
figure('Name','1: Ex3 Heliocentric Earth-Venus transfer','NumberTitle','off');
hold on
grid on
box on
%axis equal tight
%zlim([-1,1])
xlabel('X [AU]','FontSize',15)
ylabel('Y [AU]      ','FontSize',15,'Rotation',0,'VerticalAlignment','middle','HorizontalAlignment','center')
zlabel('Z [AU]','FontSize',15)
title('ECLIPJ2000 @Sun','FontSize',15)

plot3(0,0,0,'Marker','o','MarkerFaceColor','#EDB120','MarkerEdgeColor','k','MarkerSize',12,'LineStyle','none')

plot3(x0Earth(1),x0Earth(2),x0Earth(3),'Marker','o','MarkerFaceColor',"#0072BD",'MarkerEdgeColor','k','MarkerSize',10,'LineStyle','none')

STR_ = {'\textbf{Sun}','$\mathbf{Earth\ @t_0}$','$\mathbf{Venus\ @t_f}$'};
STR = STR_;    
str = cell(1,1);
thr = cell(1,1);

for k = [1,length(TR)]
    xM = ANS(k).xM.xM;
    out = ANS(k).out;
    thr{k} = num2str(TR(k)*800);
    
    plot3(xM(1),xM(2),xM(3),'Marker','o','MarkerFaceColor',"#D95319",'MarkerEdgeColor','k','MarkerSize',8,'LineStyle','none')
    plot3(out(:,1),out(:,2),out(:,3),'LineWidth',2)

    str(2*k-1:2*k,:) = {['\boldmath$\underline{r}\textrm{\ :\textsf{\ ',thr{k},'}\ mN thrust}$'],''};
    STR = {STR{:},str{2*k-1:2*k}};

    legend(STR,'Location','south', ...
        'Interpreter','latex',...
        'FontSize',17);
end

%% Prime Vector

STR = STR(1:end-1);
str = cell(1,1);
title('\textbf{Primer Vector / Thrust Angles} \boldmath${\underline{\hat{\alpha}}^*}$','\textbf{ECLIPJ2000 @Sun}',...
      'FontSize',17,'Interpreter','latex')
for k = [1,length(TR)]
    out = ANS(k).out;
    alpha = ANS(k).alpha;
    quiver3(out(:,1),out(:,2),out(:,3),alpha(:,1),alpha(:,2),alpha(:,3),'LineWidth',1);
    
    str(k,:) = {['\boldmath${\underline{\hat{\alpha}}^*}\textrm{:\textsf{\ ',thr{k},'}\ mN thrust}$']};
    STR = {STR{:},str{k}};

    legend(STR,'Location','south', ...
           'Interpreter','latex',...
           'FontSize',17);
end

%% Vectorial Costates

L = char(955);
figure('Name',sprintf('2: costate %sr',L),'NumberTitle','off')
hold on
grid on
box on
%axis equal tight
%zlim([-1,1])
xlabel('X [AU]','FontSize',15)
ylabel('Y [AU]      ','FontSize',15,'Rotation',0,'VerticalAlignment','middle','HorizontalAlignment','center')
zlabel('Z [AU]','FontSize',15)
title('\boldmath${\underline{\lambda}_r}$','\textbf{ECLIPJ2000 @Sun}',...
      'FontSize',17,'Interpreter','latex')

plot3(0,0,0,'Marker','o','MarkerFaceColor','#EDB120','MarkerEdgeColor','k','MarkerSize',12,'LineStyle','none')
plot3(x0Earth(1),x0Earth(2),x0Earth(3),'Marker','o','MarkerFaceColor',"#0072BD",'MarkerEdgeColor','k','MarkerSize',10,'LineStyle','none')

for k = [1,length(TR)]
    xM = ANS(k).xM.xM;
    plot3(xM(1),xM(2),xM(3),'Marker','o','MarkerFaceColor',"#D95319",'MarkerEdgeColor','k','MarkerSize',8,'LineStyle','none')
end
for k = [1,length(TR)]
    out = ANS(k).out;
    plot3(out(:,1),out(:,2),out(:,3),'LineWidth',2)
end
for k = [1,length(TR)]
    out = ANS(k).out;
    Lr = ANS(k).Lr;
    quiver3(out(:,1),out(:,2),out(:,3),Lr(:,1),Lr(:,2),Lr(:,3),'LineWidth',1);
end

STR_ = {'\textbf{Sun}','$\mathbf{Earth\ @t_0}$','$\mathbf{Venus\ @t_f}$'};
for k = [1,length(TR)]
    str0(k,:) = {''};
    str1(k,:) = {['\boldmath$\underline{r}\textrm{\ :\textsf{\ ',thr{k},'}\ mN thrust}$']};
    str2(k,:) = {['\boldmath${\underline{\lambda}_r}\textrm{:\textsf{\ ',thr{k},'}\ mN thrust}$']};
    str3(k,:) = {['\boldmath${\underline{\lambda}_v}\textrm{:\textsf{\ ',thr{k},'}\ mN thrust}$']};
end
STR = {STR_{:},'',str1{1},str1{length(TR)},str2{1},str2{length(TR)}};
legend(STR,'Location','south', ...
       'Interpreter','latex',...
       'FontSize',17);

figure('Name',sprintf('3: costate %sv',L),'NumberTitle','off')
hold on
grid on
box on
%axis equal tight
%zlim([-1,1])
xlabel('X [AU]','FontSize',15)
ylabel('Y [AU]      ','FontSize',15,'Rotation',0,'VerticalAlignment','middle','HorizontalAlignment','center')
zlabel('Z [AU]','FontSize',15)
title('\boldmath${\underline{\lambda}_v}$','\textbf{ECLIPJ2000 @Sun}', ...
      'FontSize',17,'Interpreter','latex')

plot3(0,0,0,'Marker','o','MarkerFaceColor','#EDB120','MarkerEdgeColor','k','MarkerSize',12,'LineStyle','none')
plot3(x0Earth(1),x0Earth(2),x0Earth(3),'Marker','o','MarkerFaceColor',"#0072BD",'MarkerEdgeColor','k','MarkerSize',10,'LineStyle','none')

for k = [1,length(TR)]
    xM = ANS(k).xM.xM;
    plot3(xM(1),xM(2),xM(3),'Marker','o','MarkerFaceColor',"#D95319",'MarkerEdgeColor','k','MarkerSize',8,'LineStyle','none')
end
for k = [1,length(TR)]
    out = ANS(k).out;
    plot3(out(:,1),out(:,2),out(:,3),'LineWidth',2)
end
for k = [1,length(TR)]
    out = ANS(k).out;
    Lv = ANS(k).Lv;
    quiver3(out(:,1),out(:,2),out(:,3),Lv(:,1),Lv(:,2),Lv(:,3),'LineWidth',1);
end

STR = {STR_{:},'',str1{1},str1{length(TR)},str3{1},str3{length(TR)}};

legend(STR,'Location','south', ...
        'Interpreter','latex',...
        'FontSize',17);

%% Plot Mass / Switching function & scalar costates

for k = [1,length(TR)]
    out = ANS(k).out;
    m = out(:,7)/m0; % [-] mass ratio
    Lv = ANS(k).Lv;
    Lm = ANS(k).Lm;
    St = zeros(length(m),1);
    for j = 1:length(out)
        St(j) = switchingfun(m(j),Lv(j,:),Lm(j,:));
    end
    ANS(k).m = m;
    ANS(k).St = St;
end

figure('Name','4: Thrust ratio, Mass ratio, Switching function & scalar costates','NumberTitle','off')

t = tiledlayout(6,1,'TileSpacing','tight','Padding','compact');
title(t,'Time-Optimal quantities & costates','FontSize',15,'FontWeight','bold')
xlabel(t,'time [years]','FontSize',15,'FontWeight','bold')
for i = 1:6
    nexttile(i)
    grid on
    hold on
    box on
    xlim tight
    ylim padded
end

    Tr = max(TR/4);
for k = [1,length(TR)]
    t_out = ANS(k).t_out;
    nexttile(1)
    tr = ones(length(t_out),1)/4*TR(k);
    plot(t_out-t0,tr,'LineWidth',1.5);
    title('Thrust ratio (ref. 800 mN)','FontSize',14)
    ylim([0 Tr+0.2])
    yticks([0 1])

    nexttile(2)
    m = ANS(k).m;
    plot(t_out-t0,m,'LineWidth',1.5);
    title('Mass ratio (ref. m_0)','FontSize',14)
    nexttile(3)
    St = ANS(k).St;
    plot(t_out-t0,St,'LineWidth',1.5)
    title('Switching function S_t','FontSize',14)
    nexttile(4)
    lr = ANS(k).lr;
    plot(t_out-t0,lr,'LineWidth',1.5);
    title('\boldmath$|\!|{\underline{\lambda}_r}|\!|$','Interpreter','latex','FontSize',17)
    nexttile(5)
    lv = ANS(k).lv;
    plot(t_out-t0,lv,'LineWidth',1.5);
    title('\boldmath$|\!|{\underline{\lambda}_v}|\!|$','Interpreter','latex','FontSize',17)
    nexttile(6)
    Lm = ANS(k).Lm;
    plot(t_out-t0,Lm,'LineWidth',1.5);
    title('\boldmath${\lambda_m}$','Interpreter','latex','FontSize',17)
    str(k,:) = {[thr{k},' mN thrust']};
end

nexttile(1)
legend('800 mN thrust','500 mN thrust','Location','east','FontSize',14)

%% Hamiltonian

H800 = hamiltonian(ANS(1).out,800*1e-6);
t_out800 = ANS(1).t_out;
figure()
plot(t_out800-t0,H800,'LineWidth',2)
ylabel('Hamiltonian','FontSize',15,'FontWeight','bold')
xlabel('time [years]','FontSize',15,'FontWeight','bold')
ax=gca; ax.FontSize = 15; ax.FontWeight = 'bold';

H500 = hamiltonian(ANS(end).out,500*1e-6);
t_out500 = ANS(end).t_out;
figure()
plot(t_out500-t0,H500,'LineWidth',2)
ylabel('Hamiltonian','FontSize',15,'FontWeight','bold')
xlabel('time [years]','FontSize',15,'FontWeight','bold')
ax=gca; ax.FontSize = 15; ax.FontWeight = 'bold';

%% Clear Kernel Pool

fprintf('\nEx 3 execution time = %.2f s\n',toc(ticsolve))
cspice_kclear();

%% Functions

function [mu,ustar,Tmax,Isp,g0] = constants % [AU,YEARS,KG]
mu = cspice_bodvrd('Sun','GM',1); % [km^3/s^2]
ustar = 1; % thrust throttle factor for time-optimal problem
Tmax = 800*1e-6; % [kN --> kg*Km/s^2] Max Thrust
Isp = 3120; % [s] % Specific impulse
g0 = 9.80665*1e-3; % [km/s^2] gravitational acceleration at sea level

km2au = cspice_convrt(1,'KM','AU');
sec2year = cspice_convrt(1,'SECONDS','YEARS');

mu = mu * km2au^3/sec2year^2; % [AU^3/years^2]
Tmax = Tmax * km2au/sec2year^2; % [kg*AU/years^2]
Isp = Isp * sec2year; % [years]
g0 = g0 * km2au/sec2year^2; % [AU/years^2]

end

%----------------------------------------------------------------------

function dy = TPBVP(t,x,TR) % 2 point boundary value problem

% variables
r = x(1:3);
v = x(4:6);
m = x(7);
Lr = x(8:10);
Lv = x(11:13);
Lm = x(14);

[mu,ustar,Tmax,Isp,g0] = constants;

% dynamics
rdot = v;
vdot = -mu*r/norm(r)^3 - ustar*TR*Tmax/m * Lv/norm(Lv);
mdot = -ustar*TR*Tmax/(Isp*g0);
Lrdot = -3*mu/norm(r)^5*dot(r,Lv)*r + mu/norm(r)^3*Lv;
Lvdot = -Lr;
Lmdot = -ustar*norm(Lv)*TR*Tmax/m^2;

dy = [rdot;vdot;mdot;Lrdot;Lvdot;Lmdot];

end

%----------------------------------------------------------------------

function [y,out,t_out] = orbit_propagator(x0,t0,tf,TR)
opt = odeset('AbsTol',1e-12,'RelTol',1e-12);
[t_out,out] = ode113(@TPBVP,[t0 tf],x0',opt,TR);
y = out(end,:)';
end

%----------------------------------------------------------------------

function F = boundariesfun(x,TR,x_start) 

% {λ0,tf} zero-finding variables 

Lm0 = abs(x(7)); % λm guess > 0
Dtf = abs(x(8)); % Δt guess > 0

Lx = x(1:6);

y0 = [x_start(1:7);Lx;Lm0];
t0 = x_start(8);
tf = Dtf + t0;
Y = orbit_propagator(y0,t0,tf,TR);
rf = Y(1:3);
vf = Y(4:6);
mf = Y(7);
Lrf = Y(8:10);
Lvf = Y(11:13);
Lmf = Y(14);

tf = cspice_convrt(tf,'YEARS','SECONDS');

xM = cspice_spkezr('Venus',tf,'ECLIPJ2000','none','Sun'); % [km][Km/s]

km2au = cspice_convrt(1,'KM','AU');
sec2year = cspice_convrt(1,'SECONDS','YEARS');

xM(1:3) = xM(1:3)*km2au;
xM(4:6) = xM(4:6)*km2au/sec2year;

rM = xM(1:3);
vM = xM(4:6);

F = [rf - rM; vf - vM; Lmf];

end

%----------------------------------------------------------------------

function St = switchingfun(m,Lv,Lm) % Switching function
[~,~,~,Isp,g0] = constants;
St = - norm(Lv)*Isp*g0/m - Lm;
end

%----------------------------------------------------------------------

function H = hamiltonian(out,Tmax)

r = [out(:,1) out(:,2) out(:,3)];
v = [out(:,4) out(:,5) out(:,6)];
m = out(:,7);
Lr = [out(:,8) out(:,9) out(:,10)];
Lv = [out(:,11) out(:,12) out(:,13)];
Lm = out(:,14);

[mu,ustar,~,Isp,g0] = constants; % [AU,YEARS,KG]
km2au = cspice_convrt(1,'KM','AU');
sec2year = cspice_convrt(1,'SECONDS','YEARS');
Tmax = Tmax * km2au/sec2year^2; % [kg*AU/years^2]

for i=1:length(r)
    H(i) = 1 + dot(Lr(i,:),v(i,:)) - mu/norm(r(i,:))^3*dot(r(i,:),Lv(i,:)) + Tmax/(Isp*g0)*ustar*(- norm(Lv(i,:))*Isp*g0/m(i,:) - Lm(i,:));
end

end





