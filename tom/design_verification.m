%MAE 154A UAV Design Verification

%This master script takes the design layout input parameters (size,
%power train, etc), randomizes them, and runs calculations to determine 
%whether the design meets the mission specifications and is a stable and 
%reasonable system. Resulting data is then saved to an excel file 

clear; close all; clc;

%Import the spreadsheets for structs
engines = xlsread('Engine_Database.xlsx'); %Remember to assort by increasing
avionics = xlsread('Avionics_Weight_Budget.xlsx');
controls = xlsread('Control_Weight_Budget.xlsx');
airfoils = xlsread('Airfoil_Database.xlsx', 1);
airfoil_drag = xlsread('Airfoil_Database.xlsx', 2);

%Physical constants
rho_10k = 17.56e-4; %Atmos. Density at 10k ft [slugs/ft^3]
rho_5k = .002048;   %Atoms. Density at 5k ft [slugs/ft^3]
rho_sl = 23.77e-4; %Atmos. Density at sl [slugs/ft^3]

mu_sl = 3.737e-7;  %Air dynamic viscosity at sl [lb s/ft^2]
mu_10k = 3.534e-7; %Air dynamic viscosity at 10k ft [lb s/ft^2]
nu_sl = mu_sl/rho_sl; %Air kinematic viscosity at sl [ft^2/s]
nu_10k = mu_10k/rho_10k; %Air kinematic viscosity at 10k ft [ft^2/s]

g = 32.2; % Acceleration due to gravity [ft/s^2]

%Mission Specs -> used to verify design
endur = 2; %Endurance [hrs]
RC = 1500;  %Rate of Climb [fpm]
RC = RC/60; %Rate of Climb [fps]
R_cruise = 60; %Cruising 
h_to = 315; %Take off altitude [ft]
ceiling = 10000; %Ceiling [ft]
v_max_sl = 150;  %Max Speed @ SL [mph]
v_max_sl = v_max_sl*5280/3600; %Max speed @ SL [fps]
v_max_10k = 180;  %Max Speed @ 10k [mph]
v_max_10k = v_max_10k*5280/3600; %Max speed @ 10k [fps]
v_stall = 80; %Stall speed @ 10k [mph]
v_stall = v_stall*5280/3600; %Stall speed @ 10k [fps]
W_max = 300;   %Max weight [lbs]
W_payload = avionics(size(avionics)); %Weight of the payload [lbs]

%Propulsive Efficiencies
eta_p_climb = 0.8; %Climb prop .eff [-]
eta_p_cruise = 0.8;  %Cruise Propulsive efficiency [-]
eta_p_loit = 0.8;   %Loiter Propulsive efficiency [-]

%Number of Good and Bad Designs
n_good = 0;
n_bad = 0;

%Randomly generate designs
for n = 1:5000

%Randomly Generate a Design
af_num = ceil(rand*(size(airfoils,1) - 1)); %Aifoil Number
S_w = 3 + rand;  %Wing Surface Area [ft^2]
A = 6 + rand*2; %Aspect Ratio [-]
b_w = sqrt(A*S_w); %Wingspan [ft]
e = 0.7; %Rectangular wing efficiency [-]
lam_1_4 = rand*15*pi/180; %Wing Quarter chord sweep [rad]
lam = 1;   %Taper ratio [-]
thicc_w = airfoils(af_num,2); %Max chord thickness ratio of wing [-]
thicc_ht = airfoils(size(airfoils,1), 2); %Max chord thickness ratio of hor. tail [-]
thicc_vt = airfoils(size(airfoils,1), 2); %Max chord thickness ratio of vert. tail [-]
N = 4;  %Ultimate load factor (fixed) [-]
L_fuse = 3 + .5*rand; %Length of fuselage [ft] (70-75% of wingspan)
Wid_fuse = L_fuse*(0.1+0.1*rand);  %Width of fuselage [ft] (10-20% of fuselage length)
D_fuse = Wid_fuse; %Depth of fuselage [ft] (same as fuselage width)
l_t = 1.75 + rand*1.5; %Distance from wing 1/4 MAC to hor tail 1/4 MAC [ft]
l_v = l_t + .25*rand - .25*rand; %Distance from wing 1/4 MAC to vert tail 1/4 MAC [ft]
b_h = .5 + rand*.5; %Horizontal tail span [ft]
b_v = .5 + rand*.5; %Vertical tail span [ft]
chord_w = S_w/b_w;  %Wing Chord length [ft]
chord_f = .1 + .9*rand; %Flap chord Length [ft]
C_m = S_w/b_w; %Mean aerodynamic chord [ft]
x_cm = airfoils(af_num, 3);  %Location of max airfoil thickness
V_H = .3 + rand*.3; %Horizontal Tail Volume Ratio [-]
V_V = 0.02 + rand*.03; %Vertical Tail Volume Ratio [-]
S_ht = V_H*chord_w*S_w/l_t; %Horizontal Tail Surface Area [ft^2]
S_vt = V_V*b_w*S_w/l_v; %Vertical Tail Surface Area [ft^2]
chord_ht = S_ht/b_h; %Hor. Tail Chord [ft]
h_cg_full = .15 + rand*.1; %Center of gravity (full and empty, place fuel tank at CG) [-]
h_cg_empty = h_cg_full;
Z_t = -rand*(D_fuse/2) + rand*(D_fuse/2); %Height of the tail [ft]
W_fuel = 2 + rand*5;

v_sl = linspace(50,v_max_sl); % Velocity vector at sea level [fps]
v_10k = linspace(v_stall, v_max_10k);  %Velocity vector at 10k [fps]

%Wetted Area Estimates [ft^2]
%Wing -Raymer
if(thicc_w > 0.05)
    S_wet_w = S_w*(1.977 + 0.52*thicc_w); 
else
    S_wet_w = 2.003*S_w;
end

%Horizontal tail -Raymer
if(thicc_ht > 0.05)
    S_wet_ht = S_ht*(1.977 + 0.52*thicc_ht);
else
    S_wet_ht = 2.003*S_ht;
end

%Vertical Tail - Raymer
if(thicc_vt > 0.05)
    S_wet_vt = S_vt*(1.977 + 0.52*thicc_vt);
else
    S_wet_vt = 2.003*S_vt;
end

%Fuselage - https://onlinelibrary.wiley.com/doi/pdf/10.1002/9781118568101.app1
D_bar_fuse = .8*D_fuse; %Average fuselage width, guessing 80%
S_wet_fuse = pi*D_bar_fuse*(L_fuse - 1.3*D_bar_fuse); %Wetted area of fuselage [ft^2]
lam_fus = L_fuse/D_bar_fuse; %Fineness ratio
Vol_fuse = (pi/4)*D_bar_fuse^2*L_fuse*(1 - 2/lam_fus); %Volume of fuselage

if(chord_f > chord_w) %Flap chord can't be bigger than wing chord
    continue;
end

if(Wid_fuse < .5) %Not big enough to hold camera
    continue
end

%Initial Guesses
v_loit_i = (v_stall + v_max_10k) /2; %Loiter Velocity [fps]
v_climb_i = v_stall; %Climb Velocity [fps]
v_cruise_i = v_stall; %Cruise Velocity [fps]
W_i = 20; %Initial Weight guess [lbs]
static_margin_full_i = 0.05; %Static Margin Guess [-]
W_engine_i = 3; %Engine weight guess [lbs]
P_engine_i = 2; %Engine Power Guess [hp]
W_fuel_i = 4; %Fuel Weight [lbs]
i_t_i = -2.5*pi/180;
ind_loit_i = 1;

k = 0;
max_iter = 50;
%Thresholds for Convergence
W_thresh = 0.05; 
static_margin_full_thresh = 0.005;
i_t_thresh = 0.05;
while(k < max_iter) %Let's begin!

%-------------------------Static Weight Calculations ---------------------%

%Structure Weight Calculations (Cessna Eqs)
B = W_i*N*S_w*(1.9*A - 4)/(1 + .11*thicc_w);
W_wing = 69*(B*10^-6)^.69; %Weight of the wing [lbs]
W_fuse = .11*W_i; %Weight of the fuselage [lbs]
W_htail = 1.2*(W_i/3000)^.25*S_ht;  %Weight of Horizontal tail [lbs]
W_vtail = 1.28*S_vt; %Weight of Vertical tail [lbs]
W_booms = 2*l_t*pi*((1.375/12)^2 - (1.125/12)^2)*(2/6); %Weight of booms [lbs]
W_contsys = controls(size(controls));  %Control sys weight [lbs]

%Fuel Weights
Fuel_vol = W_fuel/6.01;  %Volume of fuel [gal]
W_eng_tot = 1.16*W_engine_i; %Total Propulsion sys weight [lbs]
W_nacelle = .175*P_engine_i;  %Nacelle Weight [lbs]
W_fuelsys = 1.25*(Fuel_vol);  %Fuel System weight (~1 lb for every gallon) [lbs]

W_tot = W_payload(1) + W_fuel + W_wing + W_fuse + W_htail + W_nacelle +...
    W_vtail + W_eng_tot + W_fuelsys + W_contsys(1) + W_booms;  %Total aircraft weight [lbs]  
    
%-----------------------------Airfoil + Lift------------------------------%

%Lift curve slopes are from Cl vs. Alpha graphs for 4412
a_w_2d = airfoils(af_num, 4);  %2-D Wing lift-curve slope [deg^-1]
a_w_2d = a_w_2d*180/pi; %2-D Wing lift-curve slope [rad^-1]
alpha_ZL = airfoils(af_num, 1); %Zero-lift AoA for NACA 4412 [deg]
alpha_ZL = alpha_ZL*pi/180; %zero lift AoA for NACA 4412 [rad]
CL_w0 = -alpha_ZL*a_w_2d;
Cl_alpha = a_w_2d; %2-D lift-curve slope [1/rad]

a_stall = 734.5; %Speed of sound @ 10k feet [mph]
a_stall = a_stall*5280/3600; %Speed of sound @ 10k feet [fps]
M_stall = v_stall/a_stall; %Mach number, vstall @10k (mph/mph)
beta_stall = (1-M_stall^2)^0.5; %Correction factor

kk_stall_w = Cl_alpha/(2*pi); %ratio between 2-d lift curve slope and elliptical lift distribution
CL_alpha_w = (2*pi*A)/(2+sqrt(((A*beta_stall)/kk_stall_w)^2+4)); %3-d lift-curve slope for wing ([-]
kk_stall_t = airfoils(size(airfoils, 1), 4)/(2*pi); %ratio between 2-d lift curve slope and elliptical lift distribution
CL_alpha_t = (2*pi*A)/(2+sqrt(((A*beta_stall)/kk_stall_t)^2+4)); %3-d lift-curve slope for wing ([-]
a_w_3d = CL_alpha_w; %3-D lift-curve slope, wing [1/rad]
a_t_3d = CL_alpha_t; %3-D lift-curve slope, tail (assume same as wing) [1/rad]
CL_0_tot = a_t_3d*(S_ht/S_w)*i_t_i;
epsilon_0 = (2*CL_0_tot)/(pi*A); 
epsilon_alpha = (2*a_w_3d)/(pi*A); % Downwash efficiency loss [-]
CL_alpha_tot = a_w_3d +...
    (S_ht/S_w)*a_t_3d*(1-epsilon_alpha); %3-D lift curve total slope for wing and tail (1/rad)

CM_acw = airfoils(af_num, 5); %Moment about the AC, (NASA Report) [-]
M_acw = CM_acw*.5*rho_10k*v_loit_i^2*S_w*chord_w; %Mom. about AC during cruise [ft-lbs]

%Find the Incidence angle for loiter
alpha_loit = (W_i/(.5*rho_10k*v_loit_i^2*S_w) + ...
    (a_t_3d*S_ht/S_w)*i_t_i)/ (a_w_3d + a_t_3d*(S_ht/S_w)*(1-epsilon_alpha)); %AOA at loiter
CL_tot_loit = CL_alpha_tot*alpha_loit + CL_0_tot; %Total Lift at loiter
CM_alpha_tot = -CL_alpha_tot*static_margin_full_i; %CM_alpha
CM_i = a_t_3d*V_H;
CL_i = -a_t_3d*(S_ht/S_w);

i_t_loit = -(CM_acw*CL_alpha_tot + CM_alpha_tot*CL_tot_loit)/...
  (CL_alpha_tot*CM_i - CM_alpha_tot*CL_i); %Incidence angle for trim at loiter [rad]

theta_f = acos(2*chord_f/chord_w - 1); %[rad]
tau = 1 - (theta_f - sin(theta_f)) / pi; % Flap effectiveness factor [-]
CM_0 = V_H*i_t_loit*a_t_3d + CM_acw; 
CL_del_e = tau*a_t_3d*(S_ht/S_w); 
CM_del_e = -tau*V_H*a_t_3d;

%Vectors for Lift/Drag Calculations
alpha_sl = (W_i./(.5*rho_sl*v_sl.^2*S_w) + (a_t_3d*S_ht/S_w)*i_t_i)./...
    (a_w_3d + a_t_3d*(S_ht/S_w)*(1-epsilon_alpha)); %AOA Vector at SL [rad]
alpha_10k = (W_i./(.5*rho_10k*v_10k.^2*S_w) + (a_t_3d*S_ht/S_w)*i_t_i)./...
    (a_w_3d + a_t_3d*(S_ht/S_w)*(1-epsilon_alpha)); %AOA

CL_tot_sl = CL_alpha_tot.*alpha_sl + CL_0_tot; %Total Lift Coeff. vector at SL [-]
CL_tot_10k = CL_alpha_tot.*alpha_10k + CL_0_tot; %Total Lift Coeff. vector at 10k [-]

del_e_sl = -(CM_0 + CM_alpha_tot.*alpha_sl)/(CM_del_e);  % Elevator deflection to trim  at SL [rad]
del_e_10k = -(CM_0 + CM_alpha_tot.*alpha_10k)/(CM_del_e); % Elevator deflection to trim  at SL [rad]

CL_w_sl = a_w_3d.*alpha_sl; %Wing Lift coeff. at SL [-]
CL_w_10k = a_w_3d.*alpha_10k; %Wing Lift coeff. at 10k [-]
CL_t_sl = a_t_3d.*alpha_sl.*(1-epsilon_alpha) - a_t_3d*i_t_loit; %Tail Lift coeff. at SL [-]
CL_t_10k = a_t_3d.*alpha_10k.*(1-epsilon_alpha) - a_t_3d*i_t_loit; %Tail Lift coeff. at 10k [-]

L_w_sl = CL_w_sl.*.5*rho_sl.*v_sl.^2*S_w; %Wing Lift Vector at SL [-]
L_w_10k = CL_w_10k.*.5*rho_10k.*v_10k.^2*S_w; %Wing Lift Vector at SL [-]
L_t_sl = CL_t_sl.*.5*rho_sl.*v_sl.^2*S_ht; %Tail Lift Vector at SL [-]
L_t_10k = CL_t_10k.*.5*rho_10k.*v_10k.^2*S_ht; %Tail Lift Vector at SL [-]

%Calculate Total Lift
L_tot_sl = .5*rho_sl*v_sl.^2.*S_w.*...
    (CL_alpha_tot.*alpha_sl + CL_i.*del_e_sl); %Total Lift at SL [lbs]
L_tot_10k = .5*rho_10k*v_10k.^2.*S_w.*...
    (CL_alpha_tot.*alpha_10k + CL_i.*del_e_10k);  %Total lift at 10k [lbs]

%-----------------------------CG/NP/SM Calculations-----------------------%

h_acw = .25;  %AC of wing, wrt leading edge of wing, in proportion to chord [-]

%From https://ocw.mit.edu/courses/aeronautics-and-astronautics/16-01
% -unified-engineering-i-ii-iii-iv-fall-2005-spring-2006/systems-labs-06/spl8.pdf
%Check tail volume ratios
if(V_H >= 0.3 && V_H <= 0.6)
    Validity.V_H = true;
else
    Validity.V_H = false;
end

if(V_V >= 0.02 && V_V <= 0.05)
    Validity.V_V = true;
else
    Validity.V_V = false;
end

h_n = h_acw + V_H*(a_t_3d/a_w_3d)*(1-epsilon_alpha); %Neutral point [-]

static_margin_full = h_n - h_cg_full;
static_margin_empty = h_n - h_cg_empty;

if(static_margin_full >= 0.05 && static_margin_full <= 0.2)
    Validity.cg_full = true;
else
    Validity.cg_full = false;
end

if(static_margin_empty >= 0.05 && static_margin_empty <= 0.2)
    Validity.cg_full = true;
else
    Validity.cg_full = false;
end

%----------------------Other Stability Calculations-----------------------%

eta = 1; %ratio of dynamic pressure at tail/dynamic pressure at wing [-]
Cm_0t = eta*V_H*a_t_3d*(epsilon_0-i_t_loit); %zero AoA moment contribution from tail
Cm_alpha_t = -eta*V_H*a_t_3d*(1-epsilon_alpha); %change in AoA moment contribution from tail [1/rad]

%Moments due to wing about CG [ft-lbs]
M_cgw_sl_full = M_acw + L_w_sl.*(h_cg_full*chord_w - h_acw*chord_w); %SL Full
M_cgw_10k_full = M_acw + L_w_10k.*(h_cg_full*chord_w - h_acw*chord_w); %10k Full
M_cgw_sl_empty = M_acw + L_w_sl.*(h_cg_empty*chord_w - h_acw*chord_w); %SL Empty
M_cgw_10k_empty= M_acw + L_w_10k.*(h_cg_empty*chord_w - h_acw*chord_w); %10k Empty

%Moment coeffs due to wing about CG [-]
CM_cgw_sl_full = M_cgw_sl_full./(.5*rho_sl*v_sl.^2*S_w*chord_w); %SL Full
CM_cgw_10k_full = M_cgw_10k_full./(.5*rho_10k*v_10k.^2*S_w*chord_w); %10k Full
CM_cgw_sl_empty = M_cgw_sl_empty./(.5*rho_sl*v_sl.^2*S_w*chord_w); %SL Empty
CM_cgw_10k_empty = M_cgw_10k_empty./(.5*rho_10k*v_10k.^2*S_w*chord_w); %10k Empty

CM_cgt_sl = -V_H*CL_t_sl; % Moment Coeff. due to tail about CG at SL[-]
CM_cgt_10k = -V_H*CL_t_10k; % Moment Coeff. due to tail about CG at 10k [-]

CM_cg_sl_full = CM_cgt_sl + CM_cgw_sl_full;
CM_cg_sl_empty = CM_cgt_sl + CM_cgw_sl_empty;
CM_cg_10k_full = CM_cgt_10k + CM_cgw_10k_full;
CM_cg_10k_empty = CM_cgt_10k + CM_cgw_10k_empty;

CL_q = 2*eta*V_H*CL_alpha_w*(1-epsilon_alpha); %lift coefficient due to pitch rate
CM_q = -(l_t/chord_w)*CL_q; %moment coefficient due to pitch rate

%Dutch Roll
B = 3; %Blaine factor, stable if > 5 but can be lower, like 3
gamma = B*b_w*CL_tot_10k(ind_loit_i); %Wing dihedral

if(gamma >= 0 && gamma <= 8)
    Validity.dihedral = true;
else
    Validity.dihedral = false;
end

%Roll Control
if(V_V*B >= .1 && V_V*B <= .2)
    Validity.roll = true;
else
    Validity.roll = false;
end

%------------------------------Drag Calculations--------------------------%

K = 1/(pi*A*e);

M_sl = v_sl/1115;  %Mach number at sl vector
M_10k = v_10k/1076; %Mach number at 10 vector
 
Re_sl = v_sl*chord_w/nu_sl; %Reynolds number vector at sl [-]
Re_10k = v_10k*chord_w/nu_10k; %Reynolds number vector at 10k [-]

C_f_sl = .455./((log10(Re_sl).^2.58).*(1 + 0.144*M_sl.^2).^.65); %Skin Friction Coeff at sl [-]
C_f_10k = .455./((log10(Re_10k).^2.58).*(1 + 0.144*M_10k.^2).^.65); %Skin Friction Coeff at 10k ft [-]

K_sl = (1 + 0.6*thicc_w/x_cm + 100*thicc_w^4)...
    *(1.34*M_sl.^0.18*cos(lam_1_4)^.28); %Form factor vector at SL [-]
K_10k = (1 + 0.6*thicc_w/x_cm + 100*thicc_w^4)...
    *(1.34*M_10k.^0.18*cos(lam_1_4)^.28); %Form factor vector at SL [-]

%Component buildup
CD0_wing_sl = K_sl*1.*C_f_sl*S_wet_w/S_w; %Parasitic drag coeff of wing at sl [-]
CD0_wing_10k = K_10k*1.*C_f_10k*S_wet_w/S_w; %Parasitic drag coeff of wing at 10k ft [-]

CD0_fuse_sl = K_sl*1.15.*C_f_sl*S_wet_fuse/S_w; %Parasitic drag coeff of fuselage at sl [-]
CD0_fuse_10k = K_10k*1.15.*C_f_10k*S_wet_fuse/S_w; %Parasitic drag coeff of fuselage at 10k [-]

CD0_ht_sl = K_sl*1.*C_f_sl*S_wet_ht/S_w; %Parasitic drag coeff of hor. tail at sl [-]
CD0_ht_10k = K_10k*1.*C_f_10k*S_wet_ht/S_w; %Parasitic drag coeff of hor tail at 10k ft [-]

CD0_vt_sl = K_sl*1.*C_f_sl*S_wet_vt/S_w; %Parasitic drag coeff of vert. tail at sl [-]
CD0_vt_10k = K_10k*1.*C_f_10k*S_wet_vt/S_w; %Parasitic drag coeff of vert tail at 10k ft [-]

CD0_tot_sl = CD0_wing_sl + CD0_fuse_sl + CD0_ht_sl + CD0_vt_sl;  %Total parastic Drag coeff at sl [-]
CD0_tot_10k = CD0_wing_10k + CD0_fuse_10k + CD0_ht_10k + CD0_vt_10k; %Total parastiic Drag coeff at 10k [-]

D_para_sl = .5*rho_sl*v_sl.^2*S_w.*CD0_tot_sl; %Parasitic drag sl [lbf]
D_para_10k = .5*rho_10k*v_10k.^2*S_w.*CD0_tot_10k; %Parasitic drag at 10k ft [lbf]

CD_iw_sl = CL_w_sl.^2*K; %Induced Drag coeff. of wing at sl [-]
CD_iw_10k = CL_w_10k.^2*K; %Induced Drag coeff of wing at 10k [-]

D_iw_sl = CD_iw_sl.*.5*rho_sl.*v_sl.^2*S_w; %Induced Wing drag at sl [lbf]
D_iw_10k = CD_iw_10k.*.5*rho_10k.*v_10k.^2*S_w; %Induced Wing drag at sl [lbf]

CD_it_sl = CL_t_sl.^2*K; %Induced Drag coeff. of tail at sl [-]
CD_it_10k = CL_t_10k.^2*K; %Induced Drag coeff of tail at 10k [-]

D_it_sl = CD_it_sl.*.5*rho_sl.*v_sl.^2*S_ht; %Induced Tail drag at sl [lbf]
D_it_10k = CD_it_10k.*.5*rho_10k.*v_10k.^2*S_ht; %Induced Tail drag at sl [lbf]

%Airfoil Drag
CD_af_sl = getAirfoilCoeffs(af_num, alpha_sl, airfoil_drag); %Get coeff from spreadsheet
CD_af_10k = getAirfoilCoeffs(af_num, alpha_10k, airfoil_drag); %Get coeff from spreadsheet
D_af_sl = CD_af_sl*.5*rho_sl.*v_sl.^2*S_w; %Total at SL [lbf]
D_af_10k = CD_af_sl*.5*rho_sl.*v_sl.^2*S_w; %Total at 10k [lbf]

%Trim Drag
A_t = b_h^2/S_ht; %Tail Aspect Ratio.
K_c = 1/(pi*A_t*.7); % [-]
CD_trim_sl = K_c*((S_ht/S_w)*CL_t_sl).^2; %Trim Drag coeff. at sl [-]
CD_trim_10k = K_c*((S_ht/S_w)*CL_t_10k).^2; %Trim Drag coeff. at 10k [-]

D_trim_sl = .5*rho_sl.*v_sl.^2.*S_ht.*CD_trim_sl; %Trim drag at sl [lbf]
D_trim_10k = .5*rho_10k.*v_10k.^2.*S_ht.*CD_trim_10k; %Trim drag at 10k [lbf]

D_tot_sl = D_para_sl + D_iw_sl + D_it_sl + D_af_sl + D_trim_sl; %Total Drag at SL [lbf]
D_tot_10k = D_para_10k + D_iw_10k + D_it_10k + D_af_10k + D_trim_10k; %Total Drag at 10k [lbf]

%---------------------------Propulsion Calcs------------------------------%

%Find Best ratios for loiter and cruise
[CL32_CD_loit, ind_loit] = max((L_tot_10k./(.5*rho_10k.*v_10k.^2*S_w)).^(3/2)./...
    (D_tot_10k./(.5*rho_10k.*v_10k.^2*S_w))); %Get max CL^(3/2)/CD and indice of loiter [-]
L_D_loit = L_tot_10k(ind_loit)/D_tot_10k(ind_loit); %Get L/D at loiter [-]
v_loit = v_10k(ind_loit); %Loiter velocity [fps]
[L_D_cr, ind_cr] = max(L_tot_10k./D_tot_10k); %Get max L/D and indice of cruise [-]
v_cruise = v_10k(ind_cr); %Cruise velocity [fps]

P_ex = RC*W_i; %Excess power [ft*lbs/s]
P_ex = P_ex/550; %Excess power [hp]

P_req_10k = D_tot_10k.*v_10k; %Power required @ 10k [ft*lbs/s]
P_req_10k = P_req_10k/550;  %Power required @ 10k[hp]

P_req_sl = D_tot_sl.*v_sl; %Power required @ SL [ft*lbs/s]
P_req_sl = P_req_sl/550;  %Power required @ SL [hp]
P_av_sl = P_ex + P_req_sl;  %Power required @ SL [hp]
[P_climb , ind_climb] = min(P_av_sl); %Get max value and indice
v_climb = v_sl(ind_climb); %Velocity of climb [fps]

%Calculate what the minimum power needed is
if(max(P_req_sl) > max(P_req_10k) && max(P_req_sl) > P_climb)
    P_needed = max(P_req_sl)/eta_p_loit;
elseif(max(P_req_10k) > max(P_req_sl) && max(P_req_10k) > P_climb)
        P_needed = max(P_req_10k)/eta_p_loit;
else
    P_needed = P_climb/eta_p_loit;
end

%Get the index of the engine that we can use
eng_index = getEngineWeight(P_needed, engines);

%Check to make sure we could get a good engine weight
if(eng_index)
    Validity.engine = true;
    W_engine = engines(eng_index,2); %Weight of the engine [lbs]
    P_engine = engines(eng_index, 1); %Engine power [hp]

else %We don't have a good engine, break out of while
    Validity.engine = false;
    break;
end %if eng_index

%------------------------------Convergence Check--------------------------%

%Weight Convergence
if(abs(W_i - W_tot) < W_thresh)
    Convergence.weight = true;
   if (W_tot <= W_max) %Check to make sure its less than max weight
       Validity.Weight = true;
   else
       Validity.Weight = false;
   end % if W_tot
else
    W_i = W_tot; %Update weight
    Convergence.weight = false;
end %if abs

%FuelWeight Convergence
if(abs(i_t_i - i_t_loit) < i_t_thresh)
    Convergence.i_t = true;
else
    i_t_i = i_t_loit; %Update weight
    Convergence.i_t = false;
end %if abs

%Loiter Velocity Convergence
if(v_loit_i == v_loit) %Discrete
    Convergence.v_loit = true;
else
    v_loit_i = v_loit;
    Convergence.v_loit = false;
end %if v_loit

%Climb Velocity Convergence
if(v_climb_i == v_climb) %Discrete
    Convergence.v_climb = true;
else
    v_climb_i = v_climb;
    Convergence.v_climb = false;
end %if v_climb

%Cruise Velocity Convergence
if(v_cruise_i == v_cruise) %Discrete
    Convergence.v_cruise = true;
else
    v_cruise_i = v_cruise;
    Convergence.v_cruise = false;
end %if v_cruise

%Engine Weight Convergence
if(W_engine_i == W_engine) %Discrete
    Convergence.W_engine = true;
else
    W_engine_i = W_engine;
    Convergence.W_engine = false;
end %if W_engine_i

%Engine Weight Convergence
if(P_engine_i == P_engine) %Discrete
    Convergence.P_engine = true;
else
    P_engine_i = P_engine;
    Convergence.P_engine = false;
end %if W_engine_i

if(ind_loit_i == ind_loit) %Discrete
    Convergence.ind_loit = true;
else
    ind_loit_i = ind_loit;
    Convergence.ind_loit = false;
end %if W_engine_i

%Check to see if all converged
converged = true;
fn_c = fieldnames(Convergence);
for j = 1:numel(fn_c)
    %If false, design is bad
    if(~Convergence.(fn_c{j}))
      converged = false;
      break; %Breakout
    end
end %for j = numel

if(converged)
    Validity.converged = true;
    break;
else
    k = k + 1;
    Validity.converged = false;
end %if converged

end %while

%-----------------------Performance Verification--------------------------%

if(Validity.converged)
%Assume Landing Fuel Consumption is Negligible, add 1.1 SF
Wf_climb = exp(-c_p_climb*(ceiling - h_to)/(eta_p_climb*...
        (1 - D_tot_sl(ind_climb)/(P_engine*550*v_climb))))*W_tot;
Wf_cr_1 = Wf_climb*(1/exp(R_cruise/2*5280*c_p_cruise...
    /(eta_p_cruise*L_D_cr))); %Fuel after cruise to fire [lbs]
Wf_loit = ((1/sqrt(Wf_cr_1)) + endur*3600*c_p_loit/...
    (eta_p_loit*CL32_CD_loit*sqrt(2*rho_10k*S_w)))^-2; %Fuel after loiter [lbs]
Wf_cr_2 = Wf_loit*(1/exp(R_cruise/2*c_p_cruise...
    /(eta_p_cruise*L_D_cr))); %Fuel after cruise from fire [lbs]

%Check to see if we have enough fuel for climb, cruise, loiter
fuel_margin = W_fuel - (W_tot - Wf_cr_2); %Fuel margin [lbf]
if(fuel_margin > 0)
    Validity.mission = true;
else
    Validity.mission = false;
end

%Check Rate of climb
RC_ac = (P_engine*550*eta_p_climb - D_tot_10k(ind_climb)*v_climb)/W_tot; %Our rate of clmb [fps]
if( RC_ac >= RC) %[ft/s]
    Validity.RC = true;
else
    Validity.RC = false;
end

%Check max speed at 10k feet
v_max_10k_ac = P_engine*550*eta_p_loit/D_tot_10k(length(D_tot_10k)); %Our max speed 10k [fps]
if( v_max_10k_ac > v_max_10k) %[ft/s]
    Validity.max_10k_speed = true;
else
    Validity.max_10k_speed = false;
end

%Check max speed at sl
v_max_sl_ac = P_engine*550*eta_p_loit/D_tot_sl(length(D_tot_sl));
if(v_max_sl_ac > v_max_sl) %[ft-lbs/s]
    Validity.max_sl_speed = true;
else
    Validity.max_sl_speed = false;
end

%Check Total Lift at Sea Level and 10k
for f =1:length(L_tot_sl)
    if(L_tot_sl(f) > W_tot)
        Validity.Lift_sl = true;
    else
        Validity.Lift_sl = false;
        break;
    end
end %for

for f =1:length(L_tot_10k)
    if(L_tot_10k(f) > W_tot)
        Validity.Lift_10k = true;
    else
        Validity.Lift_10k = false;
        break;
    end
end %for

end %if Validity

if(~converged)
    continue;
end

%-------------------------Stability Modes---------------------------------%

m = W_tot/g; %Mass of the aircraft [slugs]
X_u = -.5*rho_10k*v_loit^2*S_w/(m*v_loit)*2*CD0_tot_10k(ind_loit);
Z_u = -.5*rho_10k*v_loit^2*S_w/(m*v_loit)*2*CL_0_tot;

%Phugoid Mode
Phu = [ X_u, -g;
    -Z_u/v_loit, 0]; %Matrix
Phu_roots = eig(Phu); %Get the roots

%If negative real roots, then its stable in phugoid
if(real(Phu_roots(1)) < 0 && real(Phu_roots(2)) < 0)
    Validity.Phugoid = true;
else
    Validity.Phugoid = false;
end

%----------------------Check entire design and save-----------------------%

%Check to see if Validity Struct is good
Good_design = true;
fn = fieldnames(Validity);
for j = 1:numel(fn)
    %If false, design is bad
    if(~Validity.(fn{j}))
      Good_design = false;
      break; %Breakout
    end
end %for j = numel

if(Good_design) %If good, save the design in the struct array
    %Increase number
    n_good = n_good + 1;
    
%Plotting for sanity Checks -> Everything Appears to be normal
% figure(1)
% plot(v_sl, D_af_sl, v_sl, D_it_sl, v_sl, D_para_sl, v_sl, D_iw_sl, v_sl, D_tot_sl, v_sl, D_trim_sl);
% legend('Airfoil', 'Tail', 'Parasitic', 'Wing', 'Total', 'Trim', 'Location', 'Northwest');
% xlabel('Velocity [fps]');
% ylabel('Drag [lbf]');
% title('Drag at Sea Level');
% grid on;
% 
% figure(2)
% plot(v_10k, D_af_10k, v_10k, D_it_10k, v_10k, D_para_10k, v_10k, D_iw_10k, v_10k, D_tot_10k, v_10k, D_trim_10k);
% xlabel('Velocity [fps]');
% ylabel('Drag [lbf]');
% title('Drag at 10,000 ft');
% legend('Airfoil', 'Tail', 'Parasitic', 'Wing', 'Total', 'Trim', 'Location', 'Northwest');
% grid on;
% 
% figure(3)
% plot(v_sl, (L_tot_sl./(.5*rho_sl.*v_sl.^2*S_w))./(D_tot_sl./(.5*rho_sl.*v_sl.^2*S_w)), v_sl, ...
%     (L_tot_sl./(.5*rho_sl.*v_sl.^2*S_w)).^(3/2)./(D_tot_sl./(.5*rho_sl.*v_sl.^2*S_w)));
% xlabel('Velocity [fps]');
% ylabel('Lift/Drag Ratio [-]');
% title('Lift/Drag Ratios at Sea Level');
% legend('C_L/C_D', 'C_L^3/2/C_D');
% grid on;
% 
% figure(4)
% plot(v_10k, (L_tot_10k./(.5*rho_10k.*v_10k.^2*S_w))./(D_tot_10k./(.5*rho_10k.*v_10k.^2*S_w)),...
%     v_10k, (L_tot_10k./(.5*rho_10k.*v_10k.^2*S_w)).^(3/2)./(D_tot_10k./(.5*rho_10k.*v_10k.^2*S_w)));
% xlabel('Velocity [fps]');
% ylabel('Lift/Drag Ratio [-]');
% title('Lift/Drag Ratios at 10,000 ft');
% legend('C_L/C_D', 'C_L^3/2/C_D');
% grid on;
% 
% figure(5)
% plot(v_sl, P_engine*550.*ones(1,100)*eta_p_loit, v_sl, D_tot_sl.*v_sl);
% xlabel('Velocity [fps]');
% ylabel('Power [ft-lbs/s');
% legend('Power Available', 'Power Required', 'Location', 'Northwest');
% grid on;
% 
% figure(6)
% plot(v_10k, P_engine*550.*ones(1,100)*eta_p_loit, v_10k, D_tot_10k.*v_10k);
% xlabel('Velocity [fps]');
% ylabel('Power [ft-lbs/s');
% legend('Power Available', 'Power Required', 'Location', 'Northwest');
% grid on;

    %Design Parameters
    Good_designs(n_good).weight = W_tot;  %Total weight [lbs]
    Good_designs(n_good).S_w = S_w; %Wing Surface area [ft^2]
    Good_designs(n_good).b_w = b_w;  %Wingspan [ft]
    Good_designs(n_good).A = A;    %Aspect Ratio [-]
    Good_designs(n_good).e = e; %Rectangular wing efficiency [-]
    Good_designs(n_good).lam_1_4 = lam_1_4; %Wing Quarter chord sweep [rad]
    Good_designs(n_good).lam = lam;   %Taper ratio [-]
    Good_designs(n_good).thicc_w = thicc_w;  %Maximum thickness ratio (Last 2 digits of NACA) [-]
    Good_designs(n_good).N = N;  %Ultimate load factor [-]
    Good_designs(n_good).L_fuse = L_fuse; %Length of fuselage [ft]
    Good_designs(n_good).Wid_fuse = Wid_fuse;  %Width of fuselage [ft]
    Good_designs(n_good).D_fuse = D_fuse; %Depth of fuselage [ft]
    Good_designs(n_good).S_ht = S_ht; %Horizontal tail surface area [ft^2]  (Might have to calculate this in while loop)
    Good_designs(n_good).l_t = l_t; %Distance from wing 1/4 MAC to hor tail 1/4 MAC [ft]
    Good_designs(n_good).l_v = l_v; %Distance from wing 1/4 MAC to vert tail 1/4 MAC [ft]
    Good_designs(n_good).b_h = b_h; %Horizontal tail span [ft]
    Good_designs(n_good).S_vt = S_vt; %Vertical tail surface area [ft^2]  (Might have to calculate this in while loop)
    Good_designs(n_good).b_v = b_v; %Vertical tail span [ft]
    Good_designs(n_good).eng_ind = eng_index;  %Engine index
    Good_designs(n_good).eng_hp = engines(eng_index,1);  %Engine power [hp]
    Good_designs(n_good).W_S = W_tot/S_w;    %Wing Loading [lbs/ft^2]
    Good_designs(n_good).Preq_W = P_needed/W_tot;  %Power Loading [hp/lb]
    Good_designs(n_good).P_needed = P_needed;  %Power actually require [hp]
    Good_designs(n_good).chord_w = chord_w;  %Chord length [ft]
    Good_designs(n_good).chord_f = chord_f; %Flap chord Length [ft]
    Good_designs(n_good).C_m = C_m; %Mean aerodynamic chord [ft]
    Good_designs(n_good).chord_ht = chord_ht; %Hor. Tail chord [ft]
    Good_designs(n_good).thicc_ht = thicc_ht; %Max chord thickness ratio of hor. tail [-]
    Good_designs(n_good).thicc_vt = thicc_vt; %Max chord thickness ratio of vert. tail [-]
    Good_designs(n_good).V_H = V_H; %Horizontal Tail Volume Ratio [-]
    Good_designs(n_good).V_H = V_H; %Vertical Tail Volume Ratio [-]
    Good_designs(n_good).h_cg_full = h_cg_full; %Center of gravity full
    Good_designs(n_good).h_cg_empty = h_cg_empty; %Center of gravity empty
    Good_designs(n_good).static_margin_full = static_margin_full; % [-]
    Good_designs(n_good).static_margin_empty = static_margin_empty; % [-]
    Good_designs(n_good).del_e_sl_max = max(abs(del_e_sl)); % [rad]
    Good_designs(n_good).del_e_10k_max = max(abs(del_e_10k)); %[rad]
    Good_designs(n_good).i_t = i_t_loit; % [rad]
    Good_designs(n_good).h_n = h_n; %NP wrt LE of wing in prop to chord [-]
    Good_designs(n_good).gamma = gamma; %Dihedral angle [deg]
    Good_designs(n_good).Z_t = Z_t; %Height of the Tail from CG [ft]
    
    q = 1;
    while(alpha_sl(q)*180/pi > 8)
        q = q + 1;
        if(q == length(alpha_sl))
            break;
        end
    end
    Good_designs(n_good).v_takeoff_needed = v_sl(q);
    
    %Save the weight breakdown as well
    Good_designs(n_good).w_payload = W_payload(1);  %Weight of the payload [lbs]
    Good_designs(n_good).w_fuel = W_fuel; %Weight of the fuel [lbs]
    Good_designs(n_good).w_wing = W_wing; %Weight of the wing [lbs]
    Good_designs(n_good).w_fuselage = W_fuse;  %Weight of the fuselage [lbs]
    Good_designs(n_good).w_hor_tail = W_htail; %Weight of the horizontal tail [lbs]
    Good_designs(n_good).w_nacelle = W_nacelle; %Weight of the nacelle [lbs]
    Good_designs(n_good).w_ver_tail = W_vtail; %Weight of the vertical tail [lbs]
    Good_designs(n_good).w_engine = W_eng_tot; %Total engine weight [lbs]
    Good_designs(n_good).w_fuel_system = W_fuelsys; %Weight of the fuel system [lbs]
    Good_designs(n_good).w_control_system = W_contsys(1); %Weight of the control system [lbs]
    
    %Mission stuff
    Good_designs(n_good).v_cruise = v_cruise;
    Good_designs(n_good).v_loit = v_loit;
    Good_designs(n_good).v_climb = v_climb;
    Good_designs(n_good).L_D_loit = L_D_loit;
    Good_designs(n_good).L_D_cr= L_D_cr;
    Good_designs(n_good).fuel_margin = fuel_margin;
    Good_designs(n_good).v_max_sl = v_max_sl_ac;
    Good_designs(n_good).v_max_10k = v_max_10k_ac;
    Good_designs(n_good).RC = RC_ac;

    
    %Airfoil Stuff
    Good_designs(n_good).airfoil = af_num;
    
    %Stability Derivatives -> All at loiter
    Good_designs(n_good).CL_0 = CL_0_tot;
    Good_designs(n_good).CL_alpha = CL_alpha_tot;
    Good_designs(n_good).CD_alpha = 2*CL_tot_10k(ind_loit)*K*CL_alpha_w;
    Good_designs(n_good).CD_de = -2*CL_tot_10k(ind_loit)*K*a_t_3d*(S_ht/S_w);
    Good_designs(n_good).CL_q = CL_q;
    Good_designs(n_good).CL_del_e = CL_del_e;
    Good_designs(n_good).CD_0 = CD0_tot_10k(ind_loit);
    Good_designs(n_good).CM_0 = CM_0;
    Good_designs(n_good).CM_alpha = CM_alpha_tot;
    Good_designs(n_good).CM_alpha_dot = -2*eta*(l_t/chord_w)*V_H*a_t_3d*epsilon_alpha;
    Good_designs(n_good).CM_q = CM_q;
    Good_designs(n_good).CM_dele = CM_del_e;
    Good_designs(n_good).CY_beta = -(S_vt/S_w)*CL_alpha_t;
    Good_designs(n_good).CY_delr = 2*eta*V_V*CL_alpha_t;
    Good_designs(n_good).Cl_beta = -eta*((Z_t*S_vt)/(b_w*S_w))*CL_alpha_t;
    Good_designs(n_good).Cl_P = -a_w_3d/8;
    Good_designs(n_good).Cl_r = (CL_0_tot/4) + 2*eta*(Z_t/b_w)*V_V*CL_alpha_t;
    Good_designs(n_good).Cn_beta = eta*V_V*CL_alpha_t - (2*Vol_fuse/(S_w*b_w));
    Good_designs(n_good).Cn_P = - (CL_tot_10k(ind_loit)/8) +...
        (8/(3*pi))*(b_v/b_w)*V_V*CL_alpha_t;
    Good_designs(n_good).Cn_r = (-CD0_tot_10k(ind_loit)/4) - (2*eta*(l_v/b_w)*V_V*CL_alpha_t);

    
else
    %Increase number
    n_bad = n_bad + 1;
    
    %Design Parameters
    Bad_designs(n_bad).weight = W_tot;  %Total weight [lbs]
    Bad_designs(n_bad).S_w = S_w; %Wing Surface area [ft^2]
    Bad_designs(n_bad).b_w = b_w;  %Wingspan [ft]
    Bad_designs(n_bad).A = A;    %Aspect Ratio [-]
    Bad_designs(n_bad).e = e; %Rectangular wing efficiency [-]
    Bad_designs(n_bad).lam_1_4 = lam_1_4; %Wing Quarter chord sweep [rad]
    Bad_designs(n_bad).lam = lam;   %Taper ratio [-]
    Bad_designs(n_bad).thicc = thicc_w;  %Maximum thickness ratio (Last 2 digits of NACA) [-]
    Bad_designs(n_bad).N = N;  %Ultimate load factor [-]
    Bad_designs(n_bad).L_fuse = L_fuse; %Length of fuselage [ft]
    Bad_designs(n_bad).Wid_fuse = Wid_fuse;  %Width of fuselage [ft]
    Bad_designs(n_bad).D_fuse = D_fuse; %Depth of fuselage [ft]
    Bad_designs(n_bad).S_ht = S_ht; %Horizontal tail surface area [ft^2]  (Might have to calculate this in while loop)
    Bad_designs(n_bad).l_t = l_t; %Distance from wing 1/4 MAC to tail 1/4 MAC [ft]
    Bad_designs(n_bad).b_h = b_h; %Horizontal tail span [ft]
    Bad_designs(n_bad).S_vt = S_vt; %Vertical tail surface area [ft^2]  (Might have to calculate this in while loop)
    Bad_designs(n_bad).b_v = b_v; %Vertical tail span [ft]
    Bad_designs(n_bad).eng_ind = eng_index;  %Engine index
    Bad_designs(n_bad).eng_hp = engines(eng_index,1);  %Engine power [hp]
    Bad_designs(n_bad).W_S = W_tot/S_w;    %Wing Loading [lbs/ft^2]
    Bad_designs(n_bad).Preq_W = P_needed/W_tot;  %Power Loading [hp/lb]
    Bad_designs(n_bad).P_needed = P_needed;  %Power actually require [hp]
    Bad_designs(n_bad).chord = chord_w;  %Chord length [ft]
    Bad_designs(n_bad).chord_f = chord_f; %Flap chord Length [ft]
    Bad_designs(n_bad).C_m = C_m; %Mean aerodynamic chord [ft]
    Bad_designs(n_bad).chord_ht = chord_ht; %Hor. Tail chord [ft]
    Bad_designs(n_bad).gamma = gamma;
    Bad_designs(n_bad).Z_t = Z_t; %Height of the Tail from CG [ft]

    %Save the weight breakdown as well
    Bad_designs(n_bad).w_payload = W_payload(1);  %Weight of the payload [lbs]
    Bad_designs(n_bad).w_fuel = W_fuel; %Weight of the fuel [lbs]
    Bad_designs(n_bad).w_wing = W_wing; %Weight of the wing [lbs]
    Bad_designs(n_bad).w_fuselage = W_fuse;  %Weight of the fuselage [lbs]
    Bad_designs(n_bad).w_hor_tail = W_htail; %Weight of the horizontal tail [lbs]
    Bad_designs(n_bad).w_nacelle = W_nacelle; %Weight of the nacelle [lbs]
    Bad_designs(n_bad).w_ver_tail = W_vtail; %Weight of the vertical tail [lbs]
    Bad_designs(n_bad).w_engine = W_eng_tot; %Total engine weight [lbs]
    Bad_designs(n_bad).w_fuel_system = W_fuelsys; %Weight of the fuel system [lbs]
    Bad_designs(n_bad).w_control_system = W_contsys(1); %Weight of the control system [lbs]
    
    %Airfoil Stuff
    Bad_Designs(n_bad).airfoil_ind = af_num;
    
    %Mission stuff
    Bad_designs(n_bad).v_cruise = v_cruise;
    Bad_designs(n_bad).v_loit = v_loit;
    Bad_designs(n_bad).v_climb = v_climb;
    Bad_designs(n_bad).L_D_loit = L_D_loit;
    Bad_designs(n_bad).L_D_cr= L_D_cr;
    Bad_designs(n_bad).fuel_margin = fuel_margin;
    Bad_designs(n_bad).v_max_sl = v_max_sl_ac;
    Bad_designs(n_bad).v_max_10k = v_max_10k_ac;
    Bad_designs(n_bad).RC = RC_ac;
       
    %Validity
    Bad_designs(n_bad).Valid_V_H = Validity.V_H;
    Bad_designs(n_bad).Valid_V_V = Validity.V_V;
    Bad_designs(n_bad).Valid_cg_full = Validity.cg_full;
    Bad_designs(n_bad).Valid_dihedral = Validity.dihedral;
    Bad_designs(n_bad).Valid_engine = Validity.engine;
    Bad_designs(n_bad).Valid_converged = Validity.converged;
    Bad_designs(n_bad).Valid_Weight = Validity.Weight;
    Bad_designs(n_bad).Valid_mission = Validity.mission;
    Bad_designs(n_bad).Valid_RC = Validity.RC;
    Bad_designs(n_bad).Valid_max_10k_speed = Validity.max_10k_speed;
    Bad_designs(n_bad).Valid_max_sl_speed = Validity.max_sl_speed;
    Bad_designs(n_bad).Valid_Lift_sl = Validity.Lift_sl;
    Bad_designs(n_bad).Valid_Lift_10k = Validity.Lift_10k;
end

end %for n = 1

%Save to spreadsheet
fprintf('Number of Good Designs Found: %d\n', n_good);
fprintf('Number of Bad Designs Found: %d\n', n_bad);
fprintf('Writing Designs to Spreadsheet... ');

if(n_good)
    xlswrite('Good_Designs.xlsx',A);
    writetable(struct2table(Good_designs),'Good_Designs.xlsx');
end
if(n_bad)
    xlswrite('Bad_Designs.xlsx',A);
    writetable(struct2table(Bad_designs),'Bad_Designs.xlsx');
end
fprintf('Done!\n');
    
%TODO
%Add more airfoils to database
%TRIPLE CHECK LITERALLY EVERYTHING
