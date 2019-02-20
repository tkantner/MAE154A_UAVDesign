%Thomas Kantner
%MAE 154A UAV Design Verification

%This master script takes the design layout input parameters (size,
%power train, etc) and runs calculations to determine whether the design
%meets the mission specifications and is a stable and reasonable system

clear; close all; clc;

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

%Import the design parameters from the Excel file
params = xlsread('Valid_Designs.xlsx');

%Assign parameters in array to values
W_i = params(1);  %Takeoff weight [lbs]
S_w = params(2);  %Wing Surface Area [ft^2]
b_w = params(3);  %Wingspan [ft]
A = params(4);  %Wing aspect ratio
e = params(5); %Rectangular wing efficiency [-]
lam_1_4 = params(6); %Wing Quarter chord sweep [rad]
lam = params(7);   %Taper ratio [-]
thicc = params(8);  %Max chord thickness ratio [-]
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
P_engine = params(21); %Engine power [hp]
P_engine = P_engine*550; %Engine power [ft-lbs/s]
chord = S_w/b_2;  %Chord length [ft]
CL = 1.2;  %Coefficient of Lift [-]
h_cg = 0;   %TO DO: Add correct value from other spreadsheet

%Mission planning - ADD CORRECT VALUES
v_climb = 0; %Climb velocity [fps]
v_cruise = 0; %Cruising velocity [fps]
LD_cruise = 0; %Crusing L/D [-]
v_sl = linspace(50,v_max_sl); % Velocity vector at sea level [fps]
v_10k = linspace(v_stall, v_max_10k);  %Velocity vector at 10k [fps]
%Add more as needed...

%Check Lift at stall
L_10k = .5*rho_10k*v_stall^2*CL*S_w; %Lift force at stall [lbs]
if(L_10k < W_i)
    Validity.Lift = true; %Mark as valid
else
    Validity.Lift = false; %Mark as invalid
end


%TODO
%Check Rate of Climb- use drag calculations and power avail curve
%Check Endurance - figure out how much fuel is used for loiter
%Check Cruise range - figure out how much fuel for cruise
%Check stability derivatives - wait for those calculations

%Look at flow chart when I get my computer back for the rest
   


