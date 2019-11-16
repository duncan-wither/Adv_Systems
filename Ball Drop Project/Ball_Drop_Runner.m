%% Ball Drop Sim Runner
% this code was written to run the ball drop simulation. 
% simulation variables are as follows:
% - ball_dia
% - air_dens
% - drag_coef
% - gravity
% - coef_rest
% - mass
% - init_height
% - stop_en

% close all previous figures
close all;

%Setting Variables so they can be manipulated by function
ball_dia    = 0.067; 
air_dens    = 1.255; 
drag_coef   = 0.53; 
gravity     = -9.81;
coef_rest   = 0.75; 
mass        = 0.0577; 
init_height = 1.5; 
stop_en     = true; 

% The Settling Time Simulation:
% Setting base variables from Paper
% Running Simulation
Result = run_sim();
% Calculating Settling timeopTime', '20'
Outwith_indexes=find(Result.h>(0.02*Sim_Consts(7))); %index array of outwith bound
% index of the the step after the final outwith bound value:
Final_End_ind = Outwith_indexes(end) + 1;
%Settling Time Value
Settling_time = Result.t(Final_End_ind);
% Displaying in Console
disp("Settling Time = ")
disp(Settling_time)

%Plot for non-liniearity of time using a variable step solver
figure;plot(Result.t);xlabel('index');ylabel('Time');title('Example of Time non-linearity');

% Plot to compare Solvers
% Setting Ideal Conditions
air_dens=0;
coef_rest=1;
stop_en = false;
% Running for Variable Step
Res_var=run_sim;
% Setting Solver to Fixed Step
set_param('Ball_Drop_final', 'Solver', 'FixedStepAuto')
% Setting teh steps to match the same number as the variable step solver.
set_param('Ball_Drop_final', 'FixedStep', '20/2160')
% Running for Fixed
Res_fix=run_sim();

%Creating Height Plot
figure;plot(Res_var.t, Res_var.h,Res_fix.t,Res_fix.h);
xlabel('Time /s');ylabel('Height /m');title('Fixed vs Variable Solver');
legend('Variable Step Solver', 'Fixed Step Solver')

%Creating Total Energy Plot
figure;plot(Res_var.t, Res_var.et,Res_fix.t,Res_fix.et);
xlabel('Time /s');ylabel('Total Energy /J');title('Fixed vs Variable Solver');
legend('Variable Step Solver', 'Fixed Step Solver')

% Resetting the Solver Type
set_param('Ball_Drop_final', 'Solver', 'VariableStepAuto')

% Teminal Velovity Calculation
% Resetting  Paramiters
air_dens    = 1.255; 
coef_rest   = 0.75;
stop_en     = true;
% Setting Height
init_height = 1e3;
% setting stop time so we dont stop to early
set_param('Ball_Drop_final', 'StopTime', '40');
% running simulation
Result_term = run_sim();
% displaying result
disp("Max Speed =")
disp(max(abs(Result_term.v)))
% Reseting length of sim
set_param('Ball_Drop_final', 'StopTime', '40');

% Coef of Rest Period Estimation
%resetting height and drag
init_height = 1.5;
air_dens=0;
%set coef
coef_rest = 0.9;
% Run Sim
Res_CR_1 = run_sim();
%find times when heights are zero
times = Res_CR_1.t(find(Res_CR_1.h <= 0));
%Calculate bounce periods
bp1=[times-circshift(times,1)];
%ignore first value
bp1=bp1(2:end);

% Calculated values
t_0=bp1(1);
est_curve=t_0*(coef_rest).^((1:0.1:length(bp1))-1);

%Plotting Both
figure
plot(bp1,'*');hold on;
plot(1:0.1:length(bp1),est_curve)
title('Ball Drop Period vs. Calculation')
legend('Measured Period','Calculated Period')
xlabel('Bounce Number');ylabel('Bounce Period /s');hold off

%For Coef Rest > 1
%set coef
coef_rest = 1.1;
% Run Sim
Res_CR_2 = run_sim();
%find times when heights are zero
times = Res_CR_2.t(find(Res_CR_2.h <= 0));
%Calculate bounce periods
bp2=[times-circshift(times,1)];
%ignore first value
bp2=bp2(2:end);

% Calculated values
t_0=bp2(1);
est_curve=t_0*(coef_rest).^((1:0.1:length(bp2))-1);

%Plotting Both
figure
plot(bp2,'*');hold on;
plot(1:0.1:length(bp2),est_curve)
title('Ball Drop Period vs. Calculation')
legend('Measured Period','Calculated Period')
xlabel('Bounce Number');ylabel('Bounce Period /s');hold off

% Total Energy
figure
plot(Res_var.t,Res_var.et)
title('Total Energy from Ideal Simulation')
xlabel('Time /s'); ylabel('Energy /J');

% ITF Test
% reseting CoR
coef_rest   = 0.75; 
% done in a vaccum.
air_dens = 0;
% bounce height = 2.54m
init_height = 2.54;
Res_ITF = run_sim();
figure;plot(Res_ITF.t,Res_ITF.h);xlabel('Time .s');ylabel('Height /m');
title('ITF Bounce Test');

% Custom function was written to run the simulation and simplify the output
% structure along side doing the calculations for energy.
function output = run_sim()
mass = evalin('base', 'mass');
gravity = evalin('base', 'gravity');
sim('Ball_Drop_final')
% Pulling Results into Easy To Use Struct.
output.t    = Height.time; %Time Values
output.h    = Height.signals.values;%Height Values
output.v    = Velocity.signals.values;%Velocity Values
output.a    = Accel.signals.values;%Acceleration Values
output.da   = Derived_Accel.signals.values; %Derived Acceleration Values

% Calculating Energy for the simulation
output.ek = 0.5.*mass.*output.v.^2;     %Kinetic
output.ep = mass.*abs(gravity).*output.h;    %Potential
output.et = output.ep+output.ek;          %Total

end
