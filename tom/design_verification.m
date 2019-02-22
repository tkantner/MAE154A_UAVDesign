%Thomas Kantner
%MAE 154A UAV Design Verification

%This master script takes the design layout input parameters (size,
%power train, etc) and runs calculations to determine whether the design
%meets the mission specifications and is a stable and reasonable system

clear; close all; clc;

%Import the spreadsheets for structs
engines = xlsread('Engine_Database.xlsx'); %Remember to assort by increasing
avionics = xlsread('Avionics_Weight_Budget.xlsx');
controls = xlsread('Control_Weight_Budget.xlsx');

%Physical constants
rho_10k = 17.56e-4; %Atmos. Density at 10k ft [slugs/ft^3]
rho_sl = 23.77e-4; %Atmos. Density at sl [slugs/ft^3]

%Mission Specs -> used to verify design
endur = 2; %Endurance [hrs]
RC = 1500;  %Rate of Climb [fpm]
RC = RC/60; %Rate of Climb [fps]
R_cruise = 60; %Cruising 
ceil = 10000; %Ceiling [ft]
v_max_sl = 150;  %Max Speed @ SL [mph]
v_max_sl = v_max_sl*5280/3600; %Max speed @ SL [fps]
v_max_10k = 180;  %Max Speed @ 10k [mph]
v_max_10k = v_max_10k*5280/3600; %Max speed @ 10k [fps]
v_stall = 80; %Stall speed @ 10k [mph]
v_stall = v_stall*5280/3600; %Stall speed @ 10k [mph]
W_max = 300;   %Max weight [lbs]
W_payload = avionics(size(avionics)); %Weight of the payload [lbs]
CD_0 = 0.04; %Estimate for now, refine later 

%Assign parameters in array to values
W_i = 20; %Initial Weight guess [lbs]
S_w = 1 + rand*3;  %Wing Surface Area [ft^2]
b_w = 2 + rand*4;  %Wingspan [ft]
e = 0.7; %Rectangular wing efficiency [-]
lam_1_4 = 0; %Wing Quarter chord sweep [rad]
lam = 1;   %Taper ratio [-]
thicc = .12; %Max chord thickness ratio [-]
N = 4;  %Ultimate load factor (fixed) [-]
L_fuse = 2.5 + rand*2; %Length of fuselage [ft]
Wid_fuse = .5 + rand*.5;  %Width of fuselage [ft]
D_fuse = .5 + rand*.5; %Depth of fuselage [ft]
S_ht = .5 + rand*.5; %Horizontal tail surface area [ft^2]
l_t = 1.5 + rand*1.5; %Distance from wing 1/4 MAC to tail 1/4 MAC [ft]
b_h = .5 + rand*.5; %Horizontal tail span [ft]
S_vt = .5 + rand*.5; %Vertical tail surface area [ft^2] 
b_v = .5 + rand*.5; %Vertical tail span [ft]
chord = S_w/b_w;  %Chord length [ft]
i_t = rand*4;  %Tail Incidence Angle [deg]
i_t = i_t*pi/180; %Tail Incidence angle [rad]
A = b_w^2/S_w; % Aspect Ratio
C_m = S_w/b_w; %Mean aerodynamic chord

%Cruise Conditions - On our way to the fire
eta_pr_cruise = 0.85;  %Cruise Propulsive efficiency [-]
c_p_cruise = 0.7;   %Cruise Propeller efficiency [lbs/hp/hr]

%Loiter conditions
eta_pr_loit = 0.7;   %Loiter Propulsive efficiency [-]
c_p_loit = 0.5;   %Loiter propeller efficiency [lbs/hp/hr]

%Takeoff and Landing Fractions
W_fuel_to_frac = 0.998;  %Takeoff fuel weight fraction [-]
W_fuel_climb_frac = 0.99; %Climb fuel weight fraction [-]
W_fuel_desc_frac = 0.995;  %Descent fuel weight fraction [-]
W_fuel_land_frac = 0.995;  %Landing fuel weight fraction [-]

%---------------------------Weight Calculations --------------------------%

k = 0;
max_iter = 50;
while(k < max_iter)
    %Drag calculations
    K = 1/(pi*A*e);
    v_sl = linspace(50,v_max_sl); % Velocity vector at sea level [fps]
    v_10k = linspace(v_stall, v_max_10k);  %Velocity vector at 10k [fps]
    
    D_para_sl = .5*rho_sl*v_sl.^2*S_w.*CD0; %Parasitic drag sl [lbf]
    D_para_10k = .5*rho_10k*v_10k.^2*S_w.*CD0; %Parasitic drag at 10k ft [lbf]

    D_i_sl = 2*K*W_i^2./(rho_sl*v_sl.^2*S_w); %Induced Drag sl [lbf]
    D_i_10k = 2*K*W_i^2./(rho_10k*v_10k.^2*S_w); %Induced drag at 10k ft[lbf]

    D_tot_sl = D_para_sl + D_i_sl; %Total drag at sl [lbf]
    D_tot_10k = D_para_10k + D_i_10k; %Total drag at 10k ft [lbf]
    
    [L_D_10k, i_10k] = max(W_i./D_tot_10k); %Get max L/D and indice of cruise [-]
    v_loit = v_10k(i_10k); %Loiter velocity [fps]
    v_cruise = v_loit; %Cruise velocity [fps]
    
    %Fuel Calculations
    W_fuel_cruise_frac = 1/(exp(R_cruise*c_p_cruise*...
        (375*eta_pr_cruise*L_D_10k)^-1));  %Fuel-Weight fraction used in cruise [-]
    W_fuel_loit_frac = 1/(exp((endur*v_loit)/...
        (375*eta_pr_loit/c_p_loit*L_D_10k)));   %Fuel-Weight fraction used in loiter [-]
    Misn_fuel_frac = W_fuel_to_frac * W_fuel_climb_frac * W_fuel_cruise_frac *...
        W_fuel_loit_frac * W_fuel_desc_frac * W_fuel_land_frac; %Total fuel-weight fraction [-]
    W_fuel = (1-Misn_fuel_frac)*W_i; %Total weight of fuel used [lbs]
    
    %Structure Weight Calculations (Cessna Eqs)
    B = W_i*N*S_w*(1.9*A - 4)/(1 + .11*thicc);
    W_wing = 69*(B*10^-6)^.69; %Weight of the wing [lbs]
    W_fuse = .11*W_i; %Weight of the fuselage [lbs]
    W_htail = 1.2*(W_i/3000)^.25*S_ht;  %Weight of Horizontal tail [lbs]
    W_vtail = 1.28*S_vt; %Weight of Vertical tail [lbs]
    
    %Power calculations
    P_ex = RC*W_i; %Excess power [ft*lbs/s]
    P_ex = P_ex/550; %Excess power [hp]
    
    P_req_10k = D_10k.*v_10k; %Power required @ 10k [ft*lbs/s]
    P_req_10k = P_req_10k/550;  %Power required @ 10k[hp]
    
    P_req_sl = D_sl.*v_sl; %Power required @ SL [ft*lbs/s]
    P_req_sl = P_req_sl/550;  %Power required @ SL [hp]
    P_av_sl = P_ex + P_req_sl;  %Power required @ SL [hp]
    [P_engine_sl , I_sl] = min(P_av_sl); %Get max value and indice
    P_climb = P_engine_sl; %Power the engine needs to produce @ SL [hp]
    v_climb = v_sl(I_sl); %Velocity of climb [fps]
    
    %Calculate what the minimum power needed is
    if(max(P_req_sl) > max(P_req_10k) && max(P_req_sl) > P_climb)
        P_needed = max(P_req_sl)/eta_pr_loit;
    elseif(max(P_req_10k) > max(P_req_sl) && max(P_req_10k) > P_climb)
            P_needed = max(P_req_10k)/eta_pr_loit;
    else
        P_needed = P_climb/eta_pr_loit;
    end
    
    %Get the index of the engine that we can use
    index = getEngineWeight(P_needed, engines);
    
    %Check to make sure we could get a good engine weight
    if(index ~= 0)
        W_engine = engines(index,2);
    
    %Additional weight values
    Fuel_vol = W_fuel/6.01;  %Volume of fuel [gal]
    W_eng_tot = 1.16*W_engine; %Total Propulsion sys weight [lbs]
    W_nacelle = .175*engines(index,1);  %Nacelle Weight [lbs]
    W_contsys = controls(size(controls));  %Control sys weight [lbs]
    W_fuelsys = 1.25*(114/454);  %Fuel System weight (1000 ml tank) [lbs] 
    
    W_tot = W_payload(1) + W_fuel + W_wing + W_fuse + W_htail + W_nacelle +...
        W_vtail + W_eng_tot + W_fuelsys + W_contsys(1);  %Total aircraft weight [lbs]
    
    %Check for convergence
    if(abs(W_i - W_tot) < W_thresh)
       if (W_tot <= W_max) %Check to make sure its less than max weight
           Validity.Weight = true;
       else
           Validity.Weight = false;
       end % if W_tot
       break; %Break out of loop
    else
        W_i = W_tot; %Update weight
        k = k + 1;
    end
        k = k + 1; 
    end %if index

end %while

%-----------------------------CG/NP/SM Calculations-----------------------%

%Lift curve slopes are from Cl vs. Alpha graphs for 4412
a_w = 1.50/10;  %Wing lift-curve slope [deg^-1]
a_w = a_w*360/2/pi; %Wing lift-curve slope [rad^-1]
a_t = 1.50/10; %Tail lift-curve slope [deg^-1]
a_t = a_t*360/2/pi; %Tail lift-curve slope [deg^-1]

%alpha = linspace(0, 10, 10); %Angle of Attack Vector [deg]
%alpha = alpha*pi/180; %Angle of Attack Vector [rad]
alpha = 4*pi/180; %Angle of attack [rad] -> Fix

h_acw = .25;  %AC of wing, wrt leading edge of wing, in proportion to chord [-]
epsilon_alpha = .15;  % Downwash efficiency loss [-] -> HOW TO CALCULATE THIS
M_acw = 0; %Moment about the AC, [ft-lbs] -> HOW TO CALCULATE THIS
CM_acw_cruise = M_acw/(.5*rho_10k*v_cruise^2*S_w*chord); %Mom. Coeff about AC during cruise [-]
CM_acw_loiter = M_acw/(.5*rho_10k*v_loit^2*S_w*chord); %Mom. Coeff about AC during cruise [-]

V_H = l_t*S_ht/(chord*S_w); %Tail volume ratio [-]
h_act = l_t/chord + h_cg_i; %AC of tail, wrt leading edge of wing, in proportion to chord [-]

%Declare symbols
syms h_cg_sym_cr;
syms h_cg_sym_loit; 

h_n = h_acw + V_H*(a_t/a_w)*(1-epsilon_alpha); %Neutral point [-]

%Solve for centers of gravity, wrt wing leading edge, prop to chord [-]
eq_cr = a_w*((h_cg_sym - h_acw) - V_H*(a_t/a_w)*(1-epsilon_alpha))*alpha...
    + CM_acw_cruise + V_H*a_t*i_t; %At cruise
eq_loit = a_w*((h_cg_sym - h_acw) - V_H*(a_t/a_w)*(1-epsilon_alpha))*alpha...
    + CM_acw_cruise + V_H*a_t*i_t; %At loiter

h_cg_solution_set_cr = vpa(solve(eq_cr == 0, h_cg_sym_cr)); %Define solution
h_cg_cr = double(h_cg_solution_set_cr(1)); %Solve Eq
h_cg_solution_set_loit = vpa(solve(eq_loit == 0, h_cg_sym_loit)); %Define solution
h_cg_loit = double(h_cg_solution_set_loit(1)); %Solve Eq

static_margin_cr = h_n - h_cg_cr; %Static Margin during cruise [-]
static_margin_loit = h_n - h_cg_loit; %Static Margin during loiter [-]

%Check for stability
%CG must be greater than neutral point
if(h_cg_cr < h_n) 
    Validity.CG_cr = true;
else
    Validity.CG_cr = false;
end

if(h_cg_loit < h_n)
    Validity.CG_loit = true;
else
    Validity.CG_loit = false;
end

%----------------------Other Stability Calculations-----------------------%

%Moments + Coefficients due to wing about CG [-]
M_cgw_loit = M_acw + L_w_10k*(h_cg_loit*chord - h_acw*chord); %Loiter
M_cgw_cr = M_acw + L_w_10k*(h_cg_cr*chord - h_acw*chord); %Cruise

CM_cgw_loit = M_cgw_loit/(.5*rho_10k*v_loit^2*S_w*chord); %Loiter
CM_cgw_cr = M_cgw_cr/(.5*rho_10k*v_cruise^2*S_w*chord); %Cruise

alpha_t = (1-epsilon_alpha)*alpha - i_t; %Tail Eff. Angle of Attack [rad]
CL_t = a_t*alpha_t; %Tail coeff. of lift [-]
CM_cgt = V_H*CL_t; % Moment Coeff. due to tail about CG

CM_0_loit = CM_acw_loit + V_H*a_t*i_t;
CM_0_cr = CM_acw_cr + V_H*a_t*i_t;

CL_alpha = a_w + a_t(S_t/S_w)*(1-epsilon_alpha);
CL_0 = -a_t*(S_t/S_w)*i_t;

%------------------------------Lift Calculations--------------------------%

%Check Lift at stall
L_tot_stall = .5*rho_10k*v_stall^2*((a_w + a_t*(S_t/S_w)*...
    (1 - epsilon_alpha))*alpha - a_t*(S_t/S_w)*i_t); %Total Lift [lbs]
if(L_tot_stall > W_i) %If more lift than weight
    Validity.Lift = true; %Mark as valid
else
    Validity.Lift = false; %Mark as invalid
end

%TODO
%Update power calculations with correct curves and efficiencies
%Add stability calculations
%    - Longitudinal Control
%    - Directional Control
%Elevator to Trim
%Stick Free Neutral Point
%More checks
%Add outer loop
%Good design validation + Saving to Excel