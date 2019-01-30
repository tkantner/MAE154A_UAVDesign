%MAE 154A Preliminary Aerodynamic Calculations

%This script calculates the Lift and drag coefficients at various
%velocities along with the actual forces themselves. These are then plotted
%as a function of velocity and the max L/D found at sea level and 10k feet.
%It takes the 'Master_Valid_Designs.xlsx' file as input to get the median
%values of the Weight [lbs], Wing Surface Area [ft^2], Wing span [ft],
%Aspect ratio [-], and Oswald's efficiency [-] calculated using the
%weight_calcs.m script

clear; close all; clc;

%Physical constants
rho_10k = 17.56e-4; %Atmos. Density at 10k ft [slugs/ft^3]
rho_5k = .002048;   %Atoms. Density at 5k ft [slugs/ft^3]
rho_sl = 23.77e-4; %Atmos. Density at sl [slugs/ft^3]

%Specs
endur = 2; %Endurance [hrs]
RC = 1500;  %Rate of Climb [fpm]
RC = RC/60; %Rate of Climb [fps]
ceil = 10000; %Ceiling [ft]
v_max_sl = 150;  %Max Speed @ SL [mph]
v_max_sl = v_max_sl*5280/3600; %Max speed @ SL [fps]
v_max_10k = 180;  %Max Speed @ 10k [mph]
v_max_10k = v_max_10k*5280/3600; %Max speed @ 10k [fps]
v_stall = 80; %Stall speed @ 10k [mph]
v_stall = v_stall*5280/3600; %Stall speed @ 10k [mph]

%Retrieve data from Spreadsheet
designs = xlsread('Master_Valid_Designs.xlsx');
[r, c] = size(designs); %Get size

%Design parameters
W = designs(r - 1, 2);   %Median weight from spreadsheet [lbs]
S_w = designs(r - 1, 3); %Median wing surface area from spreadsheet [ft^2]
b_w = designs(r - 1, 4); %Median wingspan from spreadsheet [ft]
A = designs(r - 1, 5);  %Median aspect ratio [-]
e = designs(r - 1, 6);  %Median Oswald's eff [-]
S_ht = designs(r - 1, 14); %Median Horizontal Tail Area [ft^2]
l_t = designs(r - 1, 15); %Median l_t [ft]
CD_0 = 0.03; %Parasitic Drag Coefficient [-]
chord = 1;  %Chord length [ft]

%Calculations at SL
v_sl = 75:1:v_max_sl; %Velocity vector [fps]
CL_sl = 2*W./(rho_sl*S_w*v_sl.^2);  %Lift coefficient [-]
CD_i_sl = CL_sl.^2/(A*e*pi);  %Induced drag coefficient [-]
CD_tot_sl = CD_i_sl + CD_0;   %Total drag coefficient [-]
D_sl = CD_tot_sl.*.5*rho_sl*S_w.*v_sl.^2;  %Total drag [lbf]

%Calculations at 10k
v_10k = v_stall:1:v_max_10k; %Velocity vector [fps]
CL_10k = 2*W./(rho_10k*S_w*v_10k.^2);  %Lift coefficient [-]
CD_i_10k = CL_10k.^2/(A*e*pi); %Induced drag coefficient [-]
CD_tot_10k = CD_i_10k + CD_0;  %Total drag coefficient [-]
D_10k = CD_tot_10k.*.5*rho_10k*S_w.*v_10k.^2;  %Total drag [lbf]

%Plotting
figure(1)
hold on;
plot(v_sl, D_sl);
plot(v_10k, D_10k);
grid on;
title({'Drag vs. Velocity';...
    [' W = ', num2str(W),' lbs, S_w = ', num2str(S_w), ' ft^2, A = ', num2str(A)]});
xlabel('Velocity [fps]'); ylabel('Drag [lbf]');
legend('Sea Level', '10,000 feet');

figure(2)
hold on;
plot(v_sl, CL_sl);
plot(v_10k, CL_10k);
grid on;
title({'CL vs. Velocity';...
    [' W = ', num2str(W),' lbs, S_w = ', num2str(S_w), ' ft^2, A = ', num2str(A)] });
xlabel('Velocity [fps]'); ylabel('CL [-]');
legend('Sea Level', '10,000 feet');


%L/D values
LD_sl_max = W/min(D_sl);
LD_10k_max = W/min(D_10k);


%Stability Derivative calculations

%Lift curve slopes are from Cl vs. Alpha graphs
a_w = 1.50/10;  %Wing lift-curve slope [deg^-1]
a_w = a_w*360/2/pi; %Wing lift-curve slope [rad^-1]
a_t = 1.50/10; %Tail lift-curve slope [deg^-1]
a_t = a_t*360/2/pi; %Tail lift-curve slope [deg^-1]

%Other Physical Parameters

h_acw = .25;  %AC of wing, wrt leading edge of wing, in proportion to chord [-]
epsilon_alpha = .15;  % Downwash efficiency loss [-]
static_margin = .15;  %Static Margin [-]
h_cg_i = 1; %COG, wrt wing leading edge, in proportion to chord inital guess[-]

iter = 50; %Number of Iterations for convergence
hcg_thresh = .01; %Threshold for convergence
i = 0;
while (i < iter)

    V_H = l_t*S_ht/(chord*S_w); %Tail volume ratio [-]
    h_act = l_t/chord + h_cg_i; %AC of tail, wrt leading edge of wing, in proportion to chord [-]

    syms h_cg_sym; %Declare symbol
    
    %Solve for center of gravity, wrt wing leading edge, prop to chord [-]
    eq = (h_acw + h_act*(V_H/(h_act - h_cg_sym))*(a_t/a_w)*(1-epsilon_alpha))...
        /(1 + (V_H/(h_act - h_cg_sym))*(a_t/a_w)*(1-epsilon_alpha))-(static_margin + h_cg_sym);
    h_cg_solution_set = vpa(solve(eq == 0, h_cg_sym));  
    h_cg = double(h_cg_solution_set(1)); %Get value
    
    h_n = static_margin + h_cg;  %Static Margin wrt wing leading edge, in proportion to chord [-]
    
    if(abs(h_cg - h_cg_i) < hcg_thresh)
        fprintf('Converged! h_cg = %d\n', h_cg);
        break;
    else
        h_cg_i = h_cg;
        i = i + 1;
    end
end
