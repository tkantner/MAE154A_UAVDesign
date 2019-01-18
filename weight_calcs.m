%MAE 154A Preliminary Weight Calculations

%This script uses randomly generated specs to determine designs that meet
%our specs. It uses an xls (Excel) file known as 'Engine_Database.xlsx'
%to pull data about engines, and uses a function called getEngineWeight to
%determine the weight from this database. Results are saved to a xls file
%called 'Valid_Designs.xls'

clear all; close all; clc;

design_num = 0;

%Physical constants
rho_10k = 17.56e-4; %Atmos. Density at 10k ft [slugs/ft^3]
rho_sl = 23.77e-4; %Atmos. Density at sl [slugs/ft^3]

%Import the engines spreadsheet, remember to sort by ascending power
engines = xlsread('Engine_Database.xlsx');

%Specs
endur = 2; %Endurance [hrs]
RC = 3000;  %Rate of Climb [fpm]
RC = RC/60; %Rate of Climb [fps]
ceil = 10000; %Ceiling [ft]
v_max = 180;  %Max Speed [mph]
v_max = v_max*5280/3600; %Max speed [fps]
W_payload = 20; %Weight of the payload

%Cruise Conditions - On our way to the fire
L_D_cruise = 12;  %L/D Ratio for cruising [-]
R_cruise = 60;    %Cruise range [mi]
eta_pr_cruise = 0.85;  %Cruise Propulsive efficiency [-]
c_p_cruise = 0.7;   %Cruise Propeller efficiency [lbs/hp/hr]

%Loiter conditions
L_D_loit = 14;   %Loiter L/D ratio [-]
v_loit = 100; %Loiter speed [mph]
eta_pr_loit = 0.7;   %Loiter Propulsive efficiency [-]
c_p_loit = 0.5;   %Loiter propeller efficiency [lbs/hp/hr]

%Takeoff and Landing 
W_fuel_to_frac = 0.998;  %Takeoff fuel weight fraction [-]
W_fuel_climb_frac = 0.99; %Climb fuel weight fraction [-]
W_fuel_desc_frac = 0.995;  %Descent fuel weight fraction [-]
W_fuel_land_frac = 0.995;  %Landing fuel weight fraction [-]

%Design Parameters (These should be changed as necessary)
W_i = 100; %Initial Weight estimate [lbs]
C_L = 1.2; %Coefficient of Lift guess [-]
C_D0 = 0.030; %Parasitic drag coeff, typical value [-]

outer_iter = 2000; 
for k = 0:outer_iter

%Random design parameter guesses
b_w = 7*rand;  %Wingspan [ft]
e = .7 + .3*rand; %Rectangular wing efficiency [-]
lam_1_4 = 0; %Wing Quarter chord sweep [rad]
lam = 1;   %Taper ratio [-]
thicc = .15*rand;  %Maximum thickness ratio (Last 2 digits of NACA) [-]
N = 1 + rand;  %Ultimate load factor [-]
L_fuse = rand*5; %Length of fuselage [ft]
Wid_fuse = rand*2;  %Width of fuselage [ft]
D_fuse = rand*2; %Depth of fuselage [ft]
S_ht = 5*rand; %Horizontal tail surface area [ft^2]
l_t = 5*rand; %Distance from wing 1/4 MAC to tail 1/4 MAC [ft]
b_h = 2*rand; %Horizontal tail span [ft]
t_HR = 2*rand; %Horizontal tail max root thickness [in]
S_vt = 5*rand; %Vertical tail surface area [ft^2] 
b_v = 2*rand; %Vertical tail span [ft]
t_VR = 2*rand; %Vertical tail max root thickness [in]

%Iteration stuff
i = 0;
max_iter = 50;
W_thresh = 10;  %Threshold of weight difference for convergence [lbs]


while(i < max_iter)
    %Wing Area Calculations
    S_w = 2*W_i/(rho_10k*C_L*v_max^2);  %Wing Surface Area [ft^2]
   
    %Drag calculations
    A = b_w^2/S_w; 
    K = 1/(pi*A*e);
    v = linspace(0,264);
    D_10k = .5*rho_10k*v.^2*S_w*C_D0 + 2*K*W_i./(rho_10k*v.^2*S_w); %Drag Polar [lbf]
    D_sl = .5*rho_sl*v.^2*S_w*C_D0 + 2*K*W_i./(rho_sl*v.^2*S_w); %Drag Polar [lbf]
    
    %Fuel Calculations
    W_fuel_cruise_frac = 1/(exp(R_cruise*c_p_cruise*...
        (375*eta_pr_cruise*L_D_cruise)^-1));  %Fuel-Weight fraction used in cruise [-]
    W_fuel_loit_frac = 1/(exp((endur*v_loit)/...
        (375*eta_pr_loit/c_p_loit*L_D_loit)));   %Fuel-Weight fraction used in loiter [-]
    Misn_fuel_frac = W_fuel_to_frac * W_fuel_climb_frac * W_fuel_cruise_frac *...
        W_fuel_loit_frac * W_fuel_desc_frac * W_fuel_land_frac; %Total fuel-weight fraction [-]
    W_fuel = (1-Misn_fuel_frac)*W_i; %Total weight of fuel used [lbs]
    
    %Structure Weight Calculations (Nicolai Eqs)
    W_wing = 96.948*((W_i^N/10^5)^0.65*(A/cos(lam_1_4))^.57*(S_w/100)^.61*...
        ((1+lam)/(2*thicc))^.36*(1+(v_max*.592))^.5)^.993; %Weight of the wing [lbs]
    W_fuse = 200*((W_i^N/10^5)^.286*(L_fuse/10)^.857*...
        ((Wid_fuse+D_fuse)/10)*(v_max*.592/100)^.338)^1.1;  %Weight of the fuselage [lbs]
    W_htail = 127*((W_i^N/10^5)^.87*(S_ht/100)^1.2*...
        (l_t/10)^.483*(b_h/t_HR)^.5)^.458;  %Weight of Horizontal tail [lbs]
    W_vtail = 98.5*((W_i^N/10^5)^.87*(S_vt/100)^1.2*...
        (b_v/t_VR)^.5)^.458; %Weight of Vertical tail [lbs]
    
    %Power calculations
    P_ex = RC*W_i; %Excess power [ft*lbs/s]
    P_ex = P_ex/550; %Excess power [hp]
    
    P_req_10k = D_10k.*v; %Power required @ 10k [ft*lbs/s]
    P_req_10k = P_req_10k/550;  %Power required @ 10k[hp]
    P_av_10k = P_ex + P_req_10k;  %Power required @ 10k [hp]
    P_engine_10k = max(P_av_10k); %Power the engine needs to produce @ 10k [hp]
    
    P_req_sl = D_sl.*v; %Power required @ SL [ft*lbs/s]
    P_req_sl = P_req_sl/550;  %Power required @ SL [hp]
    P_av_sl = P_ex + P_req_sl;  %Power required @ SL [hp]
    P_engine_sl = max(P_av_sl); %Power the engine needs to produce @ SL [hp]
    
    % Find worst case
    if(P_engine_sl > P_engine_10k)
        P_needed = P_engine_sl;
    else
        P_needed = P_engine_10k;
    end
    
    %Get the index of the engine that we can use
    index = getEngineWeight(P_needed, engines);
    
    %Check to make sure we could get a good engine weight
    if(index ~= 0)
        W_engine = engines(index,2);
    
    
    %Additional weight values
    Fuel_vol = W_fuel/6.01;  %Volume of fuel [gal]
    W_eng_tot = 2.575*(W_engine)^.922; %Total Propulsion sys weight [lbs]
    W_fuelsys = 2.49*((Fuel_vol)^.6*(1/2)^.3*(1^.2)*(1^.13))^1.21;  %Fuel system weight [lbs]
    W_contsys= 1.08*W_i^.7;  %Control sys weight [lbs]
    
    W_tot = W_payload + W_fuel + W_wing + W_fuse + W_htail +...
        W_vtail + W_eng_tot + W_fuelsys + W_contsys;  %Total aircraft weight [lbs]
    
    %Check the weight to the compared weight
    if(abs(W_i - W_tot) < W_thresh)
        
        design_num = design_num + 1;
        %If good, save the design in the struct array
        Good_designs(design_num).weight = W_tot;
        Good_designs(design_num).S_w = S_w;
        Good_designs(design_num).b_w = b_w;  %Wingspan [ft]
        Good_designs(design_num).e = e; %Rectangular wing efficiency [-]
        Good_designs(design_num).lam_1_4 = lam_1_4; %Wing Quarter chord sweep [rad]
        Good_designs(design_num).lam = lam;   %Taper ratio [-]
        Good_designs(design_num).thicc = thicc;  %Maximum thickness ratio (Last 2 digits of NACA) [-]
        Good_designs(design_num).N = N;  %Ultimate load factor [-]
        Good_designs(design_num).L_fuse = L_fuse; %Length of fuselage [ft]
        Good_designs(design_num).Wid_fuse = Wid_fuse;  %Width of fuselage [ft]
        Good_designs(design_num).D_fuse = D_fuse; %Depth of fuselage [ft]
        Good_designs(design_num).S_ht = S_ht; %Horizontal tail surface area [ft^2]  (Might have to calculate this in while loop)
        Good_designs(design_num).l_t = l_t; %Distance from wing 1/4 MAC to tail 1/4 MAC [ft]
        Good_designs(design_num).b_h = b_h; %Horizontal tail span [ft]
        Good_designs(design_num).t_HR = t_HR; %Horizontal tail max root thickness [in]
        Good_designs(design_num).S_vt = S_vt; %Vertical tail surface area [ft^2]  (Might have to calculate this in while loop)
        Good_designs(design_num).b_v = b_v; %Vertical tail span [ft]
        Good_designs(design_num).t_VR = t_VR; %Vertical tail max root thickness [in]
        Good_designs(design_num).eng_ind = index;  %Engine index
        
        %Increase iteration index to break out
        i = max_iter;
    else
        W_i = W_tot;
        i = i + 1;
    end %if abs
    
    else
        i = i + 1;
    end %if index
    
end %while

end %for

fprintf('%d Good Designs found\n', design_num);

if(design_num > 0)
    xlswrite('Valid_Designs.xlsx',A);
    writetable(struct2table(Good_designs),'Valid_Designs.xlsx');
end
