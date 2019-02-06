close all; clear; clc;
design_num=0;
rho_10k=0.001756; %density, 10k (slug/ft^3)
rho_SL=0.00238; %density, SL (slug/ft^3)
v_max_10k=180*5280/3600; %max velocity at 10k (fps)
v_min_10k=80*5280/3600; %min velocity at 10k (fps)
v_max_SL=150*5280/3600; %max velocity at SL (fps)
W_eng=3.5+0.15+0.29; %engine weight (lbf)
W_PL=7; %payload weight (lbf)
for ii=1:50000
    W_S=13+4.*rand; %wing loading from 13-17 (lbf/ft^2) 
    AR=3+6.*rand; %aspect ratio from 3-9, (unitless)
    S_w=7+2.*rand; %wing size from 7-9 (ft^2)
    S_ht=1+3.*rand; %horizontal tail size, 1-4 (ft^2)
    S_vt=1+3*rand; %vertical tail size (ft^2)
    b_v=1+3.*rand; %vertical tail size, 1-4 (ft^2)
    b_h=1.3+0.7*rand; %horizontal tail span (ft)
    W_fuel=0+20*rand; %fuel weight
    l_t=1+3.*rand; %distance from 1/4W to 1/4T (ft)
    L_fuse=4+2.*rand; %length of fuselage, 4-6 (ft)
    Wid_fuse=1+rand; %fuselage width, 1-2 (ft)
    D_fuse=1+rand; %fuselage depth, 1-2 (ft)
    W_tot=W_S*(S_ht+S_w); %total weight based on wing loading
    thicc=0.11+0.08*rand; %max thickness ratio
    t_HR=0.02+1.98*rand; %max tail root thickness (in)
    t_VR=3*rand;
    LF=1.2+0.4*rand; %load factor
    N=1.5*LF;%ultimate load factor
    lam_1_4=0; %wing quarter chord sweep (rad)
    c_r=1+4*rand; %chord root (ft)
    b_w=sqrt(AR*S_w);
    c_t=(2*S_w)/b_w-c_r;
    lam=c_t/c_r;
    W_wing = 96.948*((W_tot^N/10^5)^0.65*(AR/cos(lam_1_4))^.57*(S_w/100)^.61*...
        ((1+lam)/(2*thicc))^.36*(1+(v_max_SL*.592)/500)^.5)^.993; %Weight of the wing [lbs]
    W_fuse = 200*((W_tot^N/10^5)^.286*(L_fuse/10)^.857*...
        ((Wid_fuse+D_fuse)/10)*(v_max_SL*.592/100)^.338)^1.1;  %Weight of the fuselage [lbs]
    W_htail = 127*((W_tot^N/10^5)^.87*(S_ht/100)^1.2*...
        (l_t/10)^.483*(b_h/t_HR)^.5)^.458;  %Weight of Horizontal tail [lbs]
    W_vtail = 98.5*((W_tot^N/10^5)^.87*(S_vt/100)^1.2*...
        (b_v/t_VR)^.5)^.458; %Weight of Vertical tail [lbs]
    W_contsys= 1.08*W_tot^.7;  %Control sys weight [lbs]
    W_landgear = 5;
    W_buildup=W_wing+W_fuse+W_htail+W_vtail+W_fuel+W_eng+W_contsys+W_landgear+W_PL; 
    if abs(W_buildup-W_tot)<0.1 && l_t<L_fuse && lam>0
        design_num = design_num + 1;
        %If good, save the design in the struct array
        Good_designs(design_num).weight = W_buildup;
        Good_designs(design_num).S_w = S_w;
        Good_designs(design_num).b_w = b_w;  %Wingspan [ft]
        Good_designs(design_num).AR = AR;    %Aspect Ratio [-]
        Good_designs(design_num).lam = lam;   %Taper ratio [-]
        Good_designs(design_num).thicc = thicc;  %Maximum thickness ratio (Last 2 digits of NACA) [-]
        Good_designs(design_num).N = N;  %Ultimate load factor [-]
        Good_designs(design_num).L_fuse = L_fuse; %Length of fuselage [ft]
        Good_designs(design_num).Wid_fuse = Wid_fuse;  %Width of fuselage [ft]
        Good_designs(design_num).D_fuse= D_fuse; %Depth of fuselage [ft]
        Good_designs(design_num).S_ht = S_ht; %Horizontal tail surface area [ft^2]  (Might have to calculate this in while loop)
        Good_designs(design_num).l_t = l_t; %Distance from wing 1/4 MAC to tail 1/4 MAC [ft]
        Good_designs(design_num).b_h = b_h; %Horizontal tail span [ft]
        Good_designs(design_num).t_HR = t_HR; %Horizontal tail max root thickness [in]
        Good_designs(design_num).S_vt = S_vt; %Vertical tail surface area [ft^2]  (Might have to calculate this in while loop)
        Good_designs(design_num).b_v = b_v; %Vertical tail span [ft]
        Good_designs(design_num).t_VR = t_VR; %Vertical tail max root thickness [in]
        Good_designs(design_num).W_S = W_buildup/S_w;    %Wing Loading [lbs/ft^2]
        
     end % if (W_tot <= W_max)
    
end
%% Calculating 3-D lift coefficient with NACA 4412 airfoil
c_avg=0.5*(c_r+c_t); %average chord length (ft)
c_mgc=(2/3)*c_r*(1+lam+lam^2)/(1+lam); %mean geometric chord (ft)
c_mac=c_mgc;
y_mgc=(b_w/6)*(1+2*lam)/(1+lam); %y location of mgc (ft)
tan_LE=tan(lam_1_4)-0.5*c_r*(lam-1)/b_w;
lam_LE=atan(tan_LE)*180/pi; %leading edge angle (deg)
x_mgc=y_mgc*tan_LE; %x location of mgc (ft)
AR_vt=(b_v^2)/S_vt; %vertical tail AR (unitless)
C_lalpha=0.5/(5*pi/180); %estimating lift-curve slope for wing (1/rad)
kk=C_lalpha/(2*pi); %ratio between lift curve slope of wing and 2pi (ideal)
alpha_ZL=-5*pi/180; %Cl of wing for 2-d airfoil at alpha=0 (rad)
C_l0=-alpha_ZL*C_lalpha;% Cl0 of wing for 2-d airfoil (unitless)
M=180/734.5; %mach number, vmax @10k (mph/mph)
beta=(1-M^2)^0.5; %correction factor
C_Lalpha=(2*pi*AR)/(2+sqrt(((AR*beta)/kk)^2+4)); %3-d lift-curve slope for wing (unitless)
C_L0=-alpha_ZL*C_Lalpha; %3-d CL0 for wing (unitless)
C_L0_HT=C_L0; %3-d CL0 for tail (unitless)
C_L0_tot=C_L0+(S_ht/S_w)*C_L0_HT; %3-d CL total for wing + tail (unitless)
d_epsilon=(2*C_Lalpha)/(pi*AR); % downwash estimate (unitless)
C_Lalpha_tot=C_Lalpha+(S_ht/S_w)*C_Lalpha*(1-d_epsilon); %3-D lift curve total slope for wing and tail (1/rad)
alpha=(-5:10).*pi./180; %AoA (radians)
C_L_tot=C_L0_tot+C_Lalpha_tot.*alpha; %3-D lift coefficient for wing and tail (unitless)
figure
plot(alpha,C_L_tot)
