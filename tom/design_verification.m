%MAE 154A UAV Design Verification

%This master script takes the design layout input parameters (size,
%power train, etc), randomizes them, and runs calculations to determine 
%whether the design meets the mission specifications and is a stable and 
%reasonable system. If it is, then we save to an excel file 

clear; close all; clc;

%Import the spreadsheets for structs
engines = xlsread('Engine_Database.xlsx'); %Remember to assort by increasing
avionics = xlsread('Avionics_Weight_Budget.xlsx');
controls = xlsread('Control_Weight_Budget.xlsx');
airfoils = xlsread('Airfoil_Database.xlsx');

%Physical constants
rho_10k = 17.56e-4; %Atmos. Density at 10k ft [slugs/ft^3]
rho_5k = .002048;   %Atoms. Density at 5k ft [slugs/ft^3]
rho_sl = 23.77e-4; %Atmos. Density at sl [slugs/ft^3]

mu_sl = 3.737e-7;  %Air dynamic viscosity at sl [lb s/ft^2]
mu_10k = 3.534e-7; %Air dynamic viscosity at 10k ft [lb s/ft^2]
nu_sl = mu_sl/rho_sl; %Air kinematic viscosity at sl [ft^2/s]
nu_10k = mu_10k/rho_10k; %Air kinematic viscosity at 10k ft [ft^2/s]

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

%Climb conditions
eta_p_climb = 0.85; %Climb prop .eff [-]
%Cruise Conditions
eta_p_cruise = 0.85;  %Cruise Propulsive efficiency [-]
%Loiter conditions
eta_p_loit = 0.85;   %Loiter Propulsive efficiency [-]

%Number of Good and Bad Designs
n_good = 0;
n_bad = 0;

%Randomly generate designs
for n = 1:200

af_num = ceil(rand*size(airfoils,1));   
%af_num = 2;
%Randomly Generate a Design
S_w = 1 + rand*3;  %Wing Surface Area [ft^2]
A = 6 + rand*4; %Aspect Ratio [-]
b_w = sqrt(A*S_w); %Wingspan [ft]
e = 0.7; %Rectangular wing efficiency [-]
lam_1_4 = 0; %Wing Quarter chord sweep [rad]
lam = 1;   %Taper ratio [-]
thicc_w = airfoils(af_num,2); %Max chord thickness ratio of wing [-]
thicc_ht = .12; %Max chord thickness ratio of hor. tail [-]
thicc_vt = .08; %Max chord thickness ratio of vert. tail [-]
N = 4;  %Ultimate load factor (fixed) [-]
L_fuse = b_w*(0.7+.05*rand); %Length of fuselage [ft] (70-75% of wingspan)
Wid_fuse = L_fuse*(0.1+0.1*rand);  %Width of fuselage [ft] (10-20% of fuselage length)
D_fuse = Wid_fuse; %Depth of fuselage [ft] (same as fuselage width)
l_t = 1.5 + rand*2; %Distance from wing 1/4 MAC to hor tail 1/4 MAC [ft]
l_v = l_t + rand - rand; %Distance from wing 1/4 MAC to vert tail 1/4 MAC [ft]
b_h = .5 + rand*.5; %Horizontal tail span [ft]
b_v = .5 + rand*.5; %Vertical tail span [ft]
chord_w = S_w/b_w;  %Wing Chord length [ft]
l_wing = 1 + 2*rand; %Location of wing center from nose [ft]
chord_f = .1 + .9*rand; %Flap chord Length [ft]
C_m = S_w/b_w; %Mean aerodynamic chord [ft]
x_cm = airfoils(af_num, 3);  %Location of max airfoil thickness
V_H = .45; %Horizontal Tail Volume Ratio [-]
V_V = 0.03; %Vertical Tail Volume Ratio [-]
S_ht = V_H*chord_w*S_w/l_t; %Horizontal Tail Surface Area [ft^2]
S_vt = V_V*b_w*S_w/l_v; %Vertical Tail Surface Area [ft^2]
chord_ht = S_ht/b_h; %Hor. Tail Chord [ft]
l_avionics = .25 + rand*.5; %Location of avionics center from nose [ft]

%Velocity Vectors
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
S_wet_fuse = pi*D_fuse*(L_fuse - 1.3*D_fuse);

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
Fuel_vol = W_fuel_i/6.01;  %Volume of fuel [gal]
W_eng_tot = 1.16*W_engine_i; %Total Propulsion sys weight [lbs]
W_nacelle = .175*P_engine_i;  %Nacelle Weight [lbs]
W_fuelsys = 1.25*(Fuel_vol);  %Fuel System weight (~1 lb for every gallon) [lbs]

W_tot = W_payload(1) + W_fuel_i + W_wing + W_fuse + W_htail + W_nacelle +...
    W_vtail + W_eng_tot + W_fuelsys + W_contsys(1) + W_booms;  %Total aircraft weight [lbs]  
    
%-----------------------------Airfoil + Lift------------------------------%

%Lift curve slopes are from Cl vs. Alpha graphs for 4412
a_w_2d = airfoils(af_num, 4);  %2-D Wing lift-curve slope [deg^-1]
a_w_2d = a_w_2d*180/pi; %2-D Wing lift-curve slope [rad^-1]
alpha_ZL = airfoils(af_num, 1); %Zero-lift AoA for NACA 4412 [deg]
alpha_ZL = alpha_ZL*pi/180; %zero lift AoA for NACA 4412 [rad]
CL_w0 = -alpha_ZL*a_w_2d;

Cl_alpha = a_w_2d; %2-D lift-curve slope [1/rad]
Cl_0 = airfoils(af_num, 6); %2-d lift-curve slope 

a_stall = 734.5; %Speed of sound @ 10k feet [mph]
a_stall = a_stall*5280/3600; %Speed of sound @ 10k feet [fps]
M_stall = v_stall/a_stall; %Mach number, vstall @10k (mph/mph)
beta_stall = (1-M_stall^2)^0.5; %Correction factor

kk_stall = Cl_alpha/(2*pi); %ratio between 2-d lift curve slope and elliptical lift distribution
CL_alpha = (2*pi*A)/(2+sqrt(((A*beta_stall)/kk_stall)^2+4)); %3-d lift-curve slope for wing ([-]
a_w_3d = CL_alpha; %3-D lift-curve slope, wing [1/rad]
CL_0_HT = 0; %3-d CL0 for tail (Symmetric) [-]
CL_0_tot = Cl_0+(S_ht/S_w)*CL_0_HT; %3-d CL total for wing + tail [-]
epsilon_0 = (2*CL_0_tot)/(pi*A); 
epsilon_alpha = (2*a_w_3d)/(pi*A); % Downwash efficiency loss [-]
a_t_3d = CL_alpha; %3-D lift-curve slope, tail [1/rad]
CL_alpha_tot = CL_alpha +...
    (S_ht/S_w)*CL_alpha*(1-epsilon_alpha); %3-D lift curve total slope for wing and tail (1/rad)

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

%Find the Incidence angle for cruise
alpha_cruise = (W_i/(.5*rho_10k*v_cruise_i^2*S_w) + ...
    (a_t_3d*S_ht/S_w)*i_t_i)/ (a_w_3d + a_t_3d*(S_ht/S_w)*(1-epsilon_alpha)); %AOA at cruise
CL_tot_cruise = CL_alpha_tot*alpha_cruise + CL_0_tot; %Total Lift at cruise

i_t_cruise = -(CM_acw*CL_alpha_tot + CM_alpha_tot*CL_tot_cruise)/...
  (CL_alpha_tot*CM_i - CM_alpha_tot*CL_i); %Incidence angle for trim at cruise [rad]

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
CL_t_sl = a_t_3d.*alpha_sl; %Tail Lift coeff. at SL [-]
CL_t_10k = a_t_3d.*alpha_10k; %Tail Lift coeff. at 10k [-]

L_w_sl = CL_w_sl.*.5*rho_sl.*v_sl.^2*S_w; %Wing Lift Vector at SL [-]
L_w_10k = CL_w_10k.*.5*rho_10k.*v_10k.^2*S_w; %Wing Lift Vector at SL [-]
L_t_sl = CL_t_sl.*.5*rho_sl.*v_sl.^2*S_ht; %Tail Lift Vector at SL [-]
L_t_10k = CL_t_10k.*.5*rho_10k.*v_10k.^2*S_ht; %Tail Lift Vector at SL [-]


%-----------------------------CG/NP/SM Calculations-----------------------%
h_acw = .25;  %AC of wing, wrt leading edge of wing, in proportion to chord [-]

M_acw_stall = CM_acw*.5*rho_10k*v_stall^2*S_w*chord_w; %Mom. about AC during cruise [-]
M_acw_climb = CM_acw*.5*rho_sl*v_climb_i^2*S_w*chord_w; %Mom. about AC during cruise [-]

%From https://ocw.mit.edu/courses/aeronautics-and-astronautics/16-01
%-unified-engineering-i-ii-iii-iv-fall-2005-spring-2006/systems-labs-06/spl8.pdf
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

static_margin_full_vec = linspace(0.05, 0.2, 16);
h_cg_full = h_n - static_margin_full_vec;
cg_full = h_cg_full.*chord_w + l_wing;

for y = 1:length(static_margin_full_vec)
syms l_fueltank_sym;
eq1 = ((L_fuse/2)*W_fuse + (l_t + l_wing - (chord_w*.25) + ... %Take off CG wrt nose [ft]
    (chord_ht*.25))*(W_htail + W_vtail) + .9*L_fuse*(W_engine_i + W_nacelle) + ...
    + l_wing*(W_wing + W_contsys(1)) + l_fueltank_sym*(W_fuel_i + W_fuelsys) +...
    l_avionics*(W_payload(1)) + (l_t/2 + l_wing - (chord_w*.25))*W_booms)/(W_tot);

h_cg_solution_set = vpa(solve(eq1 == cg_full(y), l_fueltank_sym));
l_fueltank = double(h_cg_solution_set(1));

cg_empty = ((L_fuse/2)*W_fuse + (l_t + l_wing - (chord_w*.25) + ... %Empty CG wrt nose [ft]
     (chord_ht*.25))*(W_htail + W_vtail) + .8*L_fuse*(W_engine_i + W_nacelle) + ...
     + l_wing*(W_wing + W_contsys(1)) + l_fueltank*(W_fuelsys) +...
     l_avionics*(W_payload(1)) + (l_t/2 + l_wing - (chord_w*.25))*W_booms)/(W_tot);
 
h_cg_empty = (cg_empty - l_wing)/chord_w;
static_margin_empty = h_n - h_cg_empty;

%Check for stability
%SM must be positive, don't want it too low or high
% if(static_margin_full >= 0.05 && static_margin_full <= 0.2) 
%     Validity.CG_full = true;
% else
%     Validity.CG_full = false;
%     continue;
% end

if(static_margin_empty >= 0.05 && static_margin_empty <= 0.2)
    Validity.CG_empty = true;
else
    Validity.CG_empty = false;
    continue;
end

if(l_fueltank > .25*L_fuse && l_fueltank < .75*L_fuse)
    Validity.l_fueltank = true;
else
    Validity.l_fueltank = false;
    continue;
end

end %for

if(Validity.CG_empty)
    if(Validity.l_fueltank)
        static_margin_full = static_margin_full_vec(y);
    else
        k = k + 1;
        continue;
    end
else
    k = k + 1;
    continue;
end
    
%----------------MOMENT COEFFICIENT INDIVIDUAL COMPONENTS-----------------%

eta = 1; %ratio of dynamic pressure at tail/dynamic pressure at wing [-]
Cm_0t = eta*V_H*a_t_3d*(epsilon_0-i_t_loit); %zero AoA moment contribution from tail
Cm_alpha_t = -eta*V_H*a_t_3d*(1-epsilon_alpha); %change in AoA moment contribution from tail [1/rad]

%----------------------Other Stability Calculations-----------------------%

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

CL_q = 2*eta*V_H*CL_alpha*(1-epsilon_alpha); %lift coefficient due to pitch rate
CM_q = -(l_t/chord_w)*CL_q; %moment coefficient due to pitch rate

%------------------------------Drag Calculations--------------------------%
%This can probably be changed to a function to make it cleaner -> @TODO Later
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

%Airofoil Drag
CD_af_sl = airfoils(af_num, 7); %Airfoil Drag coeff.
CD_af_10k = CD_af_sl;
D_af_sl = CD_af_sl*.5*rho_sl.*v_sl.^2*S_w;
D_af_10k = CD_af_sl*.5*rho_sl.*v_sl.^2*S_w;

CD_tot_sl = CD0_tot_sl + CD_iw_sl + CD_it_sl + CD_af_sl; %Total Drag coeff at SL [-]
CD_tot_10k = CD0_tot_10k + CD_iw_10k + CD_it_10k + CD_af_10k; %Total Drag coeff at 10k [-]

D_tot_sl = D_para_sl + D_iw_sl + D_it_sl + D_af_sl; %Total Drag at SL [lbf]
D_tot_10k = D_para_10k + D_iw_10k + D_it_10k + D_af_10k; %Total Drag at 10k [lbf]

%---------------------------Propulsion Calcs------------------------------%

%Find Best ratios for loiter and cruise
[CL32_CD_loit, ind_loit] = max((W_i./(.5*rho_10k.*v_10k.^2*S_w))...
    .^(3/2)./(D_tot_10k./(.5*rho_10k.*v_10k.^2*S_w))); %Get max CL^(3/2)/CD and indice of loiter [-]
L_D_loit = W_i/D_tot_10k(ind_loit); %Get L/D at loiter [-]
v_loit = v_10k(ind_loit); %Loiter velocity [fps]
[L_D_cr, ind_cr] = max(W_i./D_tot_10k); %Get max L/D and indice of cruise [-]
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

    %Fuel Calculations
    c_p_climb = fuelConsumptionRate(P_climb)/(P_climb*3600*550); %Get SFC [ft^-1]
    c_p_climb = c_p_climb*(rho_5k/rho_sl -...
        (1 - (rho_5k/rho_sl))/7.55); % Get SFC, Guess at 5k [ft^-1]
    W_2 = exp(-c_p_climb*(ceiling - h_to)/(eta_p_climb*...
        (1 - D_tot_sl(ind_climb)/(P_engine*550*v_climb))))*W_i; %Fuel Weight after climb [lbs]

    c_p_cruise = fuelConsumptionRate(P_req_10k(ind_cr))...
        /(P_req_10k(ind_cr)*3600*550); %Get SFC [ft^-1]
    c_p_cruise = c_p_cruise*(rho_10k/rho_sl - (1 - (rho_10k/rho_sl))/7.55); %Change for alt [ft^-1]
    W_3 = exp(-(5280*R_cruise/2)*...
        c_p_cruise/(L_D_cr*eta_p_cruise))*W_2; %Fuel weight after cruise 1 [lbs]

    c_p_loit = fuelConsumptionRate(P_req_10k(ind_loit))...
        /(P_req_10k(ind_loit)*3600*550); %Convert units [ft^-1]
    c_p_loit = c_p_loit*(rho_10k/rho_sl - (1 - (rho_10k/rho_sl))/7.55); %Change for alt [ft^-1]
    W_4 = ((endur*3600*c_p_loit/(eta_p_loit*CL32_CD_loit*sqrt(2*rho_10k*S_w)))...
        + 1/sqrt(W_3))^-2; %Fuel after loiter [lbs]

    W_5 = exp(-(5280*R_cruise/2)*...
        c_p_cruise/(L_D_cr*eta_p_cruise))*W_4; %Fuel weight after cruise 2 [lbs]

    W_fuel = 1.1*(W_i - W_5); %Total Fuel Used (+10%) [lbs]

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
if(abs(W_fuel_i - W_fuel) < W_thresh)
    Convergence.W_fuel = true;
else
    W_fuel_i = W_fuel; %Update weight
    Convergence.W_fuel = false;
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

%Static Margin (Full) Convergence
if(abs(static_margin_full_i - static_margin_full_vec) < static_margin_full_thresh)
    Convergence.SM_full = true;
else
    static_margin_full_i = static_margin_full_vec;
    Convergence.SM_full = false;
end %if abs(static...

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
if((W_tot - Wf_cr_2) < W_fuel)
    Validity.mission = true;
else
    Validity.mission = false;
end

%Check Rate of climb
if((P_engine*550 - D_tot_10k(ind_climb)*v_climb) >= RC) %[ft-lbs/s]
    Validity.RC = true;
else
    Validity.RC = false;
end

%Check max speed at 10k feet
if(P_engine*550 > D_tot_10k(100)*v_max_10k) %[ft-lbs/s]
    Validity.max_10k_speed = true;
else
    Validity.max_10k_speed = false;
end

%Check max speed at sl
if(P_engine*550 > D_tot_sl(100)*v_max_sl) %[ft-lbs/s]
    Validity.max_sl_speed = true;
else
    Validity.max_sl_speed = false;
end

%Calculate Total Lift
L_tot_sl = .5*rho_sl*v_sl.^2.*S_w.*(CL_alpha_tot.*alpha_sl + CL_i.*del_e_sl); %Total Lift at SL [lbs]
L_tot_10k = .5*rho_10k*v_10k.^2.*S_w.*(CL_alpha_tot.*alpha_10k + CL_i.*del_e_10k);  %Total lift at 10k [lbs]

for f =1:length(L_tot_sl)
    %Check Total Lift
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

%----------------------Check entire design and save-----------------------%

if(~converged)
    continue;
end

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
    
    %Design Parameters
    Good_designs(n_good).weight = W_tot;  %Total weight [lbs]
    Good_designs(n_good).S_w = S_w; %Wing Surface area [ft^2]
    Good_designs(n_good).b_w = b_w;  %Wingspan [ft]
    Good_designs(n_good).A = A;    %Aspect Ratio [-]
    Good_designs(n_good).e = e; %Rectangular wing efficiency [-]
    Good_designs(n_good).lam_1_4 = lam_1_4; %Wing Quarter chord sweep [rad]
    Good_designs(n_good).lam = lam;   %Taper ratio [-]
    Good_designs(n_good).thicc = thicc_w;  %Maximum thickness ratio (Last 2 digits of NACA) [-]
    Good_designs(n_good).N = N;  %Ultimate load factor [-]
    Good_designs(n_good).L_fuse = L_fuse; %Length of fuselage [ft]
    Good_designs(n_good).Wid_fuse = Wid_fuse;  %Width of fuselage [ft]
    Good_designs(n_good).D_fuse = D_fuse; %Depth of fuselage [ft]
    Good_designs(n_good).S_ht = S_ht; %Horizontal tail surface area [ft^2]  (Might have to calculate this in while loop)
    Good_designs(n_good).l_t = l_t; %Distance from wing 1/4 MAC to tail 1/4 MAC [ft]
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
    Good_designs(n_good).l_wing = l_wing; %Location of wing center from nose [ft]
    Good_designs(n_good).l_fueltank = l_fueltank; %location fuel tank from nose [ft]
    Good_designs(n_good).l_avionics = l_avionics; %Location of avionics from nose [ft]
 

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
    Good_designs.(n_good).v_cruise = v_cruise;
    Good_designs.(n_good).v_loit = v_loit;
    Good_designs.(n_good).v_climb = v_climb;
    Good_designs.(n_good).L_D_loit = L_D_loit;
    Good_designs.(n_good).L_D_cr= L_D_cr;
    
    %Airfoil Stuff
    Good_designs(n_good).CL_tot = CL_tot; %3-D lift coefficient for wing and tail [-]
    Good_designs(n_good).CL_stall = CL_stall; %CL at Stall condition [-]
    Good_designs(n_good).CL_loit = CL_loit; %CL @ Vloit, 10k ft
    Good_designs(n_good).CL_cruise = CL_cruise;
    Good_designs(n_good).alpha_stall = alpha_stall; %AoA @ Vstall, 10k ft [rad]
    Good_designs(n_good).alpha_loit = alpha_loit; %AoA @ Vloit, 10k ft [rad]
    Good_designs(n_good).alpha_cr = alpha_cr; %AoA @ Vcruise, 10k ft [rad]
    
    %CG Stuff
    Good_designs(n_good).tau = tau; % Flap effectiveness factor [-]
    Good_designs(n_good).M_acw = M_acw; %Moment about the AC, [ft-lbs] -> HOW TO CALCULATE THIS
    Good_designs(n_good).CM_acw_cr = CM_acw_cr; %Mom. Coeff about AC during cruise [-]
    Good_designs(n_good).CM_acw_loit = CM_acw_loit; %Mom. Coeff about AC during cruise [-]
    Good_designs(n_good).V_H = V_H; %Tail volume ratio [-]
    Good_designs(n_good).h_n = h_n; %NP wrt LE of wing in prop to chord [-]
    Good_designs(n_good).cg_full = cg_full; %Location of CG from nose [ft]
    Good_designs(n_good).cg_empty = cg_empty; %Location of CG from nose [ft]
    Good_designs(n_good).static_margin_full = static_margin_full_vec; %Static Margin during cruise [-]
    Good_designs(n_good).static_margin_empty = static_margin_empty; %Static Margin during loiter [-]
    Good_designs(n_good).h_act_cr = h_act_cr; %AC of tail, wrt leading edge of wing, in proportion to chord [-]
    Good_designs(n_good).h_act_loit = h_act_loit; %AC of tail, wrt leading edge of wing, in proportion to chord [-]
    
    %Stability
    Good_designs(n_good).Cm_0t = Cm_0t; %zero AoA moment contribution from tail
    Good_designs(n_good).Cm_alphat = Cm_alphat; %change in AoA moment contribution from tail [1/rad]
    Good_designs(n_good).L_w_10k_loit = L_w_10k_loit; %Lift from wing during loiter [lbs]
    Good_designs(n_good).L_w_10k_cr = L_w_10k_cr; %Lift from wing during loiter [lbs]
    Good_designs(n_good).M_cgw_loit = M_cgw_loit_full; %Loiter
    Good_designs(n_good).M_cgw_cr = M_cgw_cr_full; %Cruise
    Good_designs(n_good).CM_cgw_loit = CM_cgw_loit; %Loiter
    Good_designs(n_good).CM_cgw_cr = CM_cgw_cr; %Cruise
    Good_designs(n_good).alpha_t_loit = alpha_t_loit; %Tail Eff. Angle of Attack @ loit [rad]
    Good_designs(n_good).alpha_t_cr = alpha_t_cr; %Tail Eff. Angle of Attack @ cruise [rad]
    Good_designs(n_good).CL_t_loit = CL_t_loit; %Tail coeff. of lift at loiter [-]
    Good_designs(n_good).CL_t_cr = CL_t_cr; %Tail coeff. of lift at cruise [-]
    Good_designs(n_good).CM_cgt_loit = CM_cgt_loit; % Moment Coeff. due to tail about CG at loiter [-]
    Good_designs(n_good).CM_cgt_cr = CM_cgt_cr; % Moment Coeff. due to tail about CG at cruise [-]
    Good_designs(n_good).CM_0_loit = CM_0_loit;
    Good_designs(n_good).CM_0_cr = CM_0_cr;
    Good_designs(n_good).CL_alpha = CL_alpha;
    Good_designs(n_good).CM_alpha_cr = CM_alpha_cr;
    Good_designs(n_good).CM_alpha_loit = CM_alpha_loit;
    Good_designs(n_good).CL_0 = CL_0;
    Good_designs(n_good).CL_del_e = CL_del_e; 
    Good_designs(n_good).CM_del_e = CM_del_e;
    Good_designs(n_good).delta_e_loit = delta_e_loit; %Elevator to trim [rad?]
    Good_designs(n_good).delta_e_cr = delta_e_cr; %Elevator to trim [rad?]
    Good_designs(n_good).CL_q = CL_q; %lift coefficient due to pitch rate
    Good_designs(n_good).CM_q = CM_q; %moment coefficient due to pitch rate

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
    Bad_designs(n_bad).l_wing = l_wing; %Location of wing center from nose [ft]
    Bad_designs(n_bad).l_fueltank = l_fueltank; %location fuel tank from nose [ft]
    Bad_designs(n_bad).l_avionics = l_avionics; %Location of avionics from nose [ft]
 

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
    Bad_designs(n_bad).CL_tot = CL_tot; %3-D lift coefficient for wing and tail [-]
    Bad_designs(n_bad).CL_stall = CL_stall; %CL at Stall condition [-]
    Bad_designs(n_bad).CL_loit = CL_loit; %CL @ Vloit, 10k ft
    Bad_designs(n_bad).CL_cruise = CL_cruise;
    Bad_designs(n_bad).alpha_stall = alpha_stall; %AoA @ Vstall, 10k ft [rad]
    Bad_designs(n_bad).alpha_loit = alpha_loit; %AoA @ Vloit, 10k ft [rad]
    Bad_designs(n_bad).alpha_cr = alpha_cr; %AoA @ Vcruise, 10k ft [rad]
    
    %Mission stuff
    Bad_designs.(n_bad).v_cruise = v_cruise;
    Bad_designs.(n_bad).v_loit = v_loit;
    Bad_designs.(n_bad).v_climb = v_climb;
    Bad_designs.(n_bad).L_D_loit = L_D_loit;
    Bad_designs.(n_bad).L_D_cr= L_D_cr;
    
    %CG Stuff
    Bad_designs(n_bad).tau = tau; % Flap effectiveness factor [-]
    Bad_designs(n_bad).M_acw = M_acw; %Moment about the AC, [ft-lbs] -> HOW TO CALCULATE THIS
    Bad_designs(n_bad).CM_acw_cr = CM_acw_cr; %Mom. Coeff about AC during cruise [-]
    Bad_designs(n_bad).CM_acw_loit = CM_acw_loit; %Mom. Coeff about AC during cruise [-]
    Bad_designs(n_bad).V_H = V_H; %Tail volume ratio [-]
    Bad_designs(n_bad).h_n = h_n; %NP wrt LE of wing in prop to chord [-]
    Bad_designs(n_bad).cg_full = cg_full; %Location of CG from nose [ft]
    Bad_designs(n_bad).cg_empty = cg_empty; %Location of CG from nose [ft]
    Bad_designs(n_bad).static_margin_full = static_margin_full_vec; %Static Margin during cruise [-]
    Bad_designs(n_bad).static_margin_empty = static_margin_empty; %Static Margin during loiter [-]
    Bad_designs(n_bad).h_act_cr = h_act_cr; %AC of tail, wrt leading edge of wing, in proportion to chord [-]
    Bad_designs(n_bad).h_act_loit = h_act_loit; %AC of tail, wrt leading edge of wing, in proportion to chord [-]
    
    %Stability
    Bad_designs(n_bad).Cm_0t = Cm_0t; %zero AoA moment contribution from tail
    Bad_designs(n_bad).Cm_alphat = Cm_alphat; %change in AoA moment contribution from tail [1/rad]
    Bad_designs(n_bad).L_w_10k_loit = L_w_10k_loit; %Lift from wing during loiter [lbs]
    Bad_designs(n_bad).L_w_10k_cr = L_w_10k_cr; %Lift from wing during loiter [lbs]
    Bad_designs(n_bad).M_cgw_loit = M_cgw_loit_full; %Loiter
    Bad_designs(n_bad).M_cgw_cr = M_cgw_cr_full; %Cruise
    Bad_designs(n_bad).CM_cgw_loit = CM_cgw_loit; %Loiter
    Bad_designs(n_bad).CM_cgw_cr = CM_cgw_cr; %Cruise
    Bad_designs(n_bad).alpha_t_loit = alpha_t_loit; %Tail Eff. Angle of Attack @ loit [rad]
    Bad_designs(n_bad).alpha_t_cr = alpha_t_cr; %Tail Eff. Angle of Attack @ cruise [rad]
    Bad_designs(n_bad).CL_t_loit = CL_t_loit; %Tail coeff. of lift at loiter [-]
    Bad_designs(n_bad).CL_t_cr = CL_t_cr; %Tail coeff. of lift at cruise [-]
    Bad_designs(n_bad).CM_cgt_loit = CM_cgt_loit; % Moment Coeff. due to tail about CG at loiter [-]
    Bad_designs(n_bad).CM_cgt_cr = CM_cgt_cr; % Moment Coeff. due to tail about CG at cruise [-]
    Bad_designs(n_bad).CM_0_loit = CM_0_loit;
    Bad_designs(n_bad).CM_0_cr = CM_0_cr;
    Bad_designs(n_bad).CL_alpha = CL_alpha;
    Bad_designs(n_bad).CM_alpha_cr = CM_alpha_cr;
    Bad_designs(n_bad).CM_alpha_loit = CM_alpha_loit;
    Bad_designs(n_bad).CL_0 = CL_0;
    Bad_designs(n_bad).CL_del_e = CL_del_e; 
    Bad_designs(n_bad).CM_del_e = CM_del_e;
    Bad_designs(n_bad).delta_e_loit = delta_e_loit; %Elevator to trim [rad?]
    Bad_designs(n_bad).delta_e_cr = delta_e_cr; %Elevator to trim [rad?]
    Bad_designs(n_bad).CL_q = CL_q; %lift coefficient due to pitch rate
    Bad_designs(n_bad).CM_q = CM_q; %moment coefficient due to pitch rate
    
    %Validity
    Bad_designs(n_bad).Valid_Weight = Validity.Weight;
    Bad_designs(n_bad).Valid_CG_cr = Validity.CG_cr;
    Bad_designs(n_bad).Valid_CG_loit = Validity.CG_loit;
    Bad_designs(n_bad).Valid_lift = Validity.Lift;
    Bad_designs(n_bad).Valid_mission = Validity.mission;
    Bad_designs(n_bad).Valid_RC = Validity.RC;
    Bad_designs(n_bad).Valid_max_10k_speed = Validity.max_10k_speed;
    Bad_designs(n_bad).Valid_max_sl_speed = Validity.max_sl_speed;
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

%-------------------------FUSELAGE SECTIONS-------------------------------%

% %perimeters 1-6 represent mid-point perimeter of sectioned fuselage starting from nose
% perimeter_1=18.89; % mid-section perimeters [in]
% perimeter_2=27.87;
% perimeter_5=25.15;
% perimeter_6=38.72-10.31*1.08;
% wf_1=perimeter_1/(pi); %estimation of average width of each section [in]
% wf_2=perimeter_2/pi;
% wf_5=25.15/pi;
% wf_6=perimeter_6/pi;
% wf_6p=wf_6+0.5*10.31*1.08;
% dx1=4.16; %length of each section, 1-6 (NEEDS TO BE UPDATED) [in]
% dx2=4;
% dx3=4;
% dx4=4;
% dx5=4;
% dx6=4;
% 
% xi5=2; %distance from fuselage mid-section 5 to wing T.E. (NEEDS TO BE UPDATED) [in]
% xi6=6; %distance from fuselage mid-section 6 to wing T.E. (NEEDS TO BE UPDATED) [in]
% eps_u1=1.3; %graphically determined upwash for sections 1&2 [1/rad]
% eps_u2=1.5; %graphically determined upwash for sections 1&2 [1/rad]
% eps_u3=0; %sections 3&4 are approximated to have no upwash or downwash [1/rad]
% eps_u4=0;
% eps_u5=(xi5/11.21)*(1-epsilon_alpha); %assumed linear progression in downwash from T.E to tail; estimation of downwash based 
%                                       %on linear model
% eps_u6=(xi6/11.21)*(1-epsilon_alpha);
% 
% Cm_alphaf=1/(36.5*S_w*144*11.21)*((wf_1^2)*eps_u1*dx1+(wf_2^2)...
%     *eps_u2*dx2+(wf_5^2)*eps_u5*dx5+(wf_6^6)*eps_u6*dx6); %estimation of contribution of fuselage to moment coeff from change in AoA [1/rad]
% k2_k1=0.6; %graphically determined coefficient to determine zero AoA fuselage contribution to moment
% Cm_0f=((k2_k1)/(36.5*S_w*144*11.21))*((wf_1^2)*dx1+(wf_2^2)...
%     *dx2+(wf_3^2)*dx3+(wf_4^2)*dx4+(wf_5^2)*dx5+(wf_6^2)*wf_6)*alpha_ZL; %estimation of contribution of fuselage to moment coeff at zero AoA [1/rad]
    
%TODO
%Why are we having issues with total lift?
%Add stability calculations
%    - Longitudinal Control
%    - Directional Control
%More checks -> How do we show its stable?
%Finish labeling everything and cleaning up code
%Maybe an exact parasitic drag?
%TRIPLE CHECK LITERALLY EVERYTHING
