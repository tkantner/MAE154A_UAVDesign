MAE 154A Preliminary Weight Calculations

%This script uses randomly generated parameters to determine designs that meet
%our specs. It uses an xls (Excel) file known as 'Engine_Database.xlsx'
%to pull data about engines, and uses a function called getEngineWeight to
%determine the weight from this database. Results are saved to a xls file
%called 'Valid_Designs.xls'

clear; close all; clc;

design_num = 0;

%Physical constants
rho_10k = 17.56e-4; %Atmos. Density at 10k ft [slugs/ft^3]
rho_5k = .002048;   %Atoms. Density at 5k ft [slugs/ft^3]
rho_sl = 23.77e-4; %Atmos. Density at sl [slugs/ft^3]

%Import the engines spreadsheet, remember to sort by ascending power
engines = xlsread('Engine_Database.xlsx');
avionics = xlsread('Avionics_Weight_Budget.xlsx');
controls = xlsread('Control_Weight_Budget.xlsx');

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
W_max = 300;  %Max TO weight [lbs] 
W_payload = avionics(size(avionics));

%Cruise Conditions - On our way to the fire
L_D_cruise = 11;  %L/D Ratio for cruising [-]
R_cruise = 60;    %Cruise range [mi]
eta_pr_cruise = 0.85;  %Cruise Propulsive efficiency [-]
c_p_cruise = 0.7;   %Cruise Propeller efficiency [lbs/hp/hr]

%Loiter conditions
L_D_loit = 11;   %Loiter L/D ratio [-]
v_loit = 110; %Loiter speed [mph]
eta_pr_loit = 0.7;   %Loiter Propulsive efficiency [-]
c_p_loit = 0.5;   %Loiter propeller efficiency [lbs/hp/hr]

%Takeoff and Landing 
W_fuel_to_frac = 0.998;  %Takeoff fuel weight fraction [-]
W_fuel_climb_frac = 0.99; %Climb fuel weight fraction [-]
W_fuel_desc_frac = 0.995;  %Descent fuel weight fraction [-]
W_fuel_land_frac = 0.995;  %Landing fuel weight fraction [-]

%Design Parameters (These should be changed as necessary)
C_L = 1.2; %Coefficient of Lift guess [-]
C_D0 = 0.040; %Parasitic drag coeff, typical value [-]

outer_iter = 1; 
for k = 1:outer_iter
W_i = 50; %Initial Weight estimate [lbs]

%Random design parameter guesses
b_w = 3;  %Wingspan [ft]
e = .7; %+ .3*rand; %Rectangular wing efficiency [-]
lam_1_4 = 0; %Wing Quarter chord sweep [rad]
lam = 1;   %Taper ratio [-]
thicc = .12;  %Maximum thickness ratio (Last 2 digits of NACA) [-]
LF = 1.2;  %Load Factor [-]
N = 4;  %Ultimate load factor [-]
L_fuse = 3.5; %Length of fuselage [ft]
Wid_fuse = .5;  %Width of fuselage [ft]
D_fuse = .5; %Depth of fuselage [ft]
S_ht = .5; %Horizontal tail surface area [ft^2]
l_t = 1.5; %Distance from wing 1/4 MAC to tail 1/4 MAC [ft]
b_h = 1; %Horizontal tail span [ft]
t_HR = 1.5; %Horizontal tail max root thickness [in]
S_vt = .5; %Vertical tail surface area [ft^2] 
b_v = .5; %Vertical tail span [ft]
t_VR = 1.5; %Vertical tail max root thickness [in]

if(l_t > L_fuse) %Check to see if l_t is greater - not possible
    continue; %Then continue to next iteration
end

%Iteration stuff
i = 0;
max_iter = 50; %Convergence Iterations
W_thresh = .1;  %Threshold of weight difference for convergence [lbs]


while(i < max_iter)
    %Wing Area Calculations
    S_w = 2*W_i/(rho_10k*C_L*v_stall^2);  %Wing Surface Area [ft^2]
    chord = S_w/b_w; %Chord Length
   
    %Drag calculations
    A = b_w^2/S_w; % Aspect Ratio
    K = 1/(pi*A*e);
    C_m = S_w/b_w; %Mean aerodynamic chord
    v_sl = linspace(50,v_max_sl); % Velocity vector at sea level [fps]
    v_10k = linspace(v_stall, v_max_10k);  %Velocity vector at 10k [fps]
    D_10k = .5*rho_10k*v_10k.^2*S_w*C_D0 + 2*K*W_i^2./(rho_10k*v_10k.^2*S_w); %Drag Polar [lbf]
    D_sl = .5*rho_sl*v_sl.^2*S_w*C_D0 + 2*K*W_i^2./(rho_sl*v_sl.^2*S_w); %Drag Polar [lbf]
    
    %Fuel Calculations
    W_fuel_cruise_frac = 1/(exp(R_cruise*c_p_cruise*...
        (375*eta_pr_cruise*L_D_cruise)^-1));  %Fuel-Weight fraction used in cruise [-]
    W_fuel_loit_frac = 1/(exp((endur*v_loit)/...
        (375*eta_pr_loit/c_p_loit*L_D_loit)));   %Fuel-Weight fraction used in loiter [-]
    Misn_fuel_frac = W_fuel_to_frac * W_fuel_climb_frac * W_fuel_cruise_frac *...
        W_fuel_loit_frac * W_fuel_desc_frac * W_fuel_land_frac; %Total fuel-weight fraction [-]
    W_fuel = (1-Misn_fuel_frac)*W_i; %Total weight of fuel used [lbs]
    
    %Structure Weight Calculations (Nicolai Eqs)
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
    v_climb = v_sl(I_sl);
    
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
    
    %Check the weight to the compared weight
    if(abs(W_i - W_tot) < W_thresh)
        if (W_tot <= W_max) %Checl to make sure its less than max weight
        design_num = design_num + 1;
        %If good, save the design in the struct array
        Good_designs(design_num).weight = W_tot;
        Good_designs(design_num).S_w = S_w;
        Good_designs(design_num).b_w = b_w;  %Wingspan [ft]
        Good_designs(design_num).A = A;    %Aspect Ratio [-]
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
        Good_designs(design_num).eng_hp = engines(index,1);  %Engine power [hp]
        Good_designs(design_num).W_S = W_tot/S_w;    %Wing Loading [lbs/ft^2]
        Good_designs(design_num).Preq_W = P_needed/W_tot;  %Power Loading [hp/lb]
        Good_designs(design_num).P_needed = P_needed;  %Power actually require [hp]
        
        %Break out
        break;
        end % if (W_tot <= W_max)
        
    else
        W_i = W_tot;
        i = i + 1;
    end %if abs
    
    else
        i = i + 1;
    end %if index
    
end %while

end %for

if(design_num > 0)
    fprintf('%d Good Designs found\n', design_num);
    xlswrite('Valid_Designs.xlsx',A);
    writetable(struct2table(Good_designs),'Valid_Designs.xlsx');
else
    fprintf('No Designs Found!');
end
