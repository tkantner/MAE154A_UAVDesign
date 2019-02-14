%Thomas Kantner
%MAE 154A Plot Generator

%This script calculates the parasitic drag and subsequent drag polars,
%power calculations, and L/D vectors and plots them
clear; close all; clc;

%Physical constants
rho_10k = 17.56e-4; %Atmos. Density at 10k ft [slugs/ft^3]
rho_5k = .002048;   %Atoms. Density at 5k ft [slugs/ft^3]
rho_sl = 23.77e-4; %Atmos. Density at sl [slugs/ft^3]

mu_sl = 3.737e-7;  %Air dynamic viscosity at sl [lb s/ft^2]
mu_10k = 3.534e-7; %Air dynamic viscosity at 10k ft [lb s/ft^2]
nu_sl = mu_sl/rho_sl; %Air kinematic viscosity at sl [ft^2/s]
nu_10k = mu_10k/rho_10k; %Air kinematic viscosity at 10k ft [ft^2/s]

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
W_payload = 7.5 ; %Weight of the payload (Max-Viz 2300 + avionics)

params = xlsread('Valid_Designs.xlsx'); %Import Spreadsheet

%Design Parameters
W_i = params(1);  %Takeoff weight [lbs]
S_w = params(2);  %Wing Surface Area [ft^2]
b_w = params(3);  %Wingspan [ft]
e = params(5); %+ .3*rand; %Rectangular wing efficiency [-]
lam_1_4 = params(6); %Wing Quarter chord sweep [rad]
lam = params(7);   %Taper ratio [-]
N = params(9);  %Ultimate load factor [-]
L_fuse = params(10); %Length of fuselage [ft]
Wid_fuse = params(11);  %Width of fuselage [ft]
D_fuse = params(12); %Depth of fuselage [ft]
S_ht = params(13); %Horizontal tail surface area [ft^2]
l_t = params(14); %Distance from wing 1/4 MAC to tail 1/4 MAC [ft]
b_h = params(15); %Horizontal tail span [ft]
t_HR = params(16); %Horizontal tail max root thickness [in]
S_vt = params(17); %Vertical tail surface area [ft^2] 
b_v = params(18); %Vertical tail span [ft]
t_VR = params(19); %Vertical tail max root thickness [in]

%Airfoil numbers: NACA 4412 for now
x_cm = .4;  %Location of max airfoil thickness
thicc = params(8);  %Maximum thickness ratio (Last 2 digits of NACA) [-]
C_L = 1.2;  %Coeff of lift [-]

%Area parameters 
A_fuselage = 70; %Max Frontal area of fuselage [in^2]
A_wing = 72; %Max frontal area of wing [in^2]
S_wet_fuse = 1091/4; %Fuselage wetted area [in^2] -> Update with correct area values
S_wet_fuse = S_wet_fuse/144; %Fuselage wetted area [ft^2]
S_wet_wing = 491; %Wing Wetted area [in^2]
S_wet_wing = S_wet_wing/144; %Wing wetted area [ft^2]
S_wet_ht = 500/2; %Horizontal tail wetted area [in^2] -> Update with correct area values
S_wet_ht = S_wet_ht/144; %Fuselage wetted area [ft^2]
S_wet_vt = 238/2;  %Vert tail wetted area [in^2]  -> Update with correct area values
S_wet_vt = S_wet_vt/144;  %Vert tail wetted area [ft^2]

chord = S_w/b_w; %Chord Length [ft]

%Stuff for component buildup method
v_sl = linspace(50,v_max_sl); % Velocity vector at sea level [fps]
v_10k = linspace(v_stall, v_max_10k);  %Velocity vector at 10k [fps]

M_sl = v_sl/1115;  %Mach number at sl vector
M_10k = v_10k/1076; %Mach number at 10 vector
 
Re_sl = v_sl*chord/nu_sl; %Reynolds number vector at sl [-]
Re_10k = v_10k*chord/nu_10k; %Reynolds number vector at 10k [-]

C_f_sl = .455./((log10(Re_sl).^2.58).*(1 + 0.144*M_sl.^2).^.65); %Skin Friction Coeff at sl [-]
C_f_10k = .455./((log10(Re_10k).^2.58).*(1 + 0.144*M_10k.^2).^.65); %Skin Friction Coeff at 10k ft [-]

K_sl = (1 + 0.6*thicc/x_cm + 100*thicc^4)...
    *(1.34*M_sl.^0.18*cos(lam_1_4)^.28); %Form factor vector at SL [-]
K_10k = (1 + 0.6*thicc/x_cm + 100*thicc^4)...
    *(1.34*M_10k.^0.18*cos(lam_1_4)^.28); %Form factor vector at SL [-]

%Component buildup
CD0_wing_sl = K_sl*1.*C_f_sl*S_wet_wing/S_w; %Parasitic drag coeff of wing at sl [-]
CD0_wing_10k = K_10k*1.*C_f_10k*S_wet_wing/S_w; %Parasitic drag coeff of wing at 10k ft [-]

CD0_fuse_sl = K_sl*1.15.*C_f_sl*S_wet_fuse/S_w; %Parasitic drag coeff of fuselage at sl [-]
CD0_fuse_10k = K_10k*1.15.*C_f_10k*S_wet_fuse/S_w; %Parasitic drag coeff of fuselage at 10k [-]

CD0_ht_sl = K_sl*1.*C_f_sl*S_wet_ht/S_w; %Parasitic drag coeff of hor. tail at sl [-]
CD0_ht_10k = K_10k*1.*C_f_10k*S_wet_ht/S_w; %Parasitic drag coeff of hor tail at 10k ft [-]

CD0_vt_sl = K_sl*1.*C_f_sl*S_wet_vt/S_w; %Parasitic drag coeff of vert. tail at sl [-]
CD0_vt_10k = K_10k*1.*C_f_10k*S_wet_vt/S_w; %Parasitic drag coeff of vert tail at 10k ft [-]

CD0_tot_sl = CD0_wing_sl + CD0_fuse_sl + CD0_ht_sl + CD0_vt_sl;  %Total parastic Drag coeff at sl [-]
CD0_tot_10k = CD0_wing_10k + CD0_fuse_10k + CD0_ht_10k + CD0_vt_10k; %Total parastiic Drag coeff at 10k [-]

%Drag calculations
A = b_w^2/S_w; % Aspect Ratio
K = 1/(pi*A*e);

D_para_sl = .5*rho_sl*v_sl.^2*S_w.*CD0_tot_sl; %Parasitic drag sl [lbf]
D_para_10k = .5*rho_10k*v_10k.^2*S_w.*CD0_tot_10k; %Parasitic drag at 10k ft [lbf]

D_i_sl = 2*K*W_i^2./(rho_sl*v_sl.^2*S_w); %Induced Drag sl [lbf]
D_i_10k = 2*K*W_i^2./(rho_10k*v_10k.^2*S_w); %Induced drag at 10k ft[lbf]

D_tot_sl = D_para_sl + D_i_sl; %Total drag at sl [lbf]
D_tot_10k = D_para_10k + D_i_10k; %Total drag at 10k ft [lbf]

%Plots
figure(1)
plot(v_sl, D_para_sl, v_sl, D_i_sl, v_sl, D_tot_sl);
title('Drag Polar at Sea Level');
legend('Parasitic', 'Induced', 'Total')
xlabel('Velocity [fps]');
ylabel('Drag [lbs]');
grid on;

figure(2)
plot(v_10k, D_para_10k, v_10k, D_i_10k, v_10k, D_tot_10k);
title('Drag Polar at 10,000 ft');
legend('Parasitic', 'Induced', 'Total', 'Location', 'northwest')
xlabel('Velocity [fps]');
ylabel('Drag [lbs]');
grid on;

%Power calculations
P_req_sl = D_tot_sl.*v_sl; %Power required at SL [ft-lbs/s]
P_req_10k = D_tot_10k.*v_10k; %Power required at 10k ft [ft-lbs/s]
P_engine = params(21); %Power from the engine [hp]
P_avail = linspace(P_engine*550, P_engine*550)*.85; %Power available [ft-lbs/s]

%Plots
figure(3)
plot(v_sl, P_req_sl, v_sl, P_avail);
title('Power Curves at Sea Level');
legend('Power required', 'Power available')
xlabel('Velocity [fps]');
ylabel('Power [ft-lbs/s]');
grid on;

figure(4)
plot(v_10k, P_req_10k, v_10k, P_avail);
title('Power Curves at 10,000 ft');
legend('Power required', 'Power available')
xlabel('Velocity [fps]');
ylabel('Power [ft-lbs/s]');
grid on;

%L/D  Calculations
%Plots
figure(5)
plot(v_sl, W_i./D_tot_sl, v_10k, W_i./D_tot_10k);
title('L/D vs. Velocity');
legend('Sea Level', '10,000 Ft')
xlabel('Velocity [fps]');
ylabel('Lift/Drag [-]');
grid on;
