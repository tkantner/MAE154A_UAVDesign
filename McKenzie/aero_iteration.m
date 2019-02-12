lam=0; %no taper (for now)
lam_1_4=0; %no wing sweep
c_r=0.78; %chord length at root (ft)
c_t=0.78; %chord length at tip (ft)
c_avg=0.5*(c_r+c_t); %average chord length (ft)
b_w=5; %wing span (ft)
b_v=1; %vertical tail span (ft)
S_vt=1; %vertical tail area (ft^2)
S_w=3.9; %wing area
AR=(b_w^2)/S_w; %aspect ratio of wing
S_ht=1; %horizontal tail area (ft^2)
c_mgc=(2/3)*c_r*(1+lam+lam^2)/(1+lam); %mean geometric chord (ft)
c_mac=c_mgc; 
y_mgc=(b_w/6)*(1+2*lam)/(1+lam); %y location of mgc (ft)
tan_LE=tan(lam_1_4)-0.5*c_r*(lam-1)/b_w;
lam_LE=atan(tan_LE)*180/pi; %leading edge angle (deg)
x_mgc=y_mgc*tan_LE; %x location of mgc (ft)
AR_vt=(b_v^2)/S_vt; %vertical tail AR (unitless)
T=readtable('/Volumes/KINGSTON/MAE 154A/airfoil_data.xlsx'); %read in airfoil data table
C_lalpha=(T{2,4}-T{1,4}*(180/pi))/((T{2,3}-T{1,3})); %estimating lift-curve slope for wing (1/rad)
kk=C_lalpha/(2*pi); %ratio between lift curve slope of wing and 2pi (ideal)
alpha_ZL=-5*pi/180; %alpha at zero lift (rad)
C_l0=-T{1,4}*C_lalpha;% Cl0 of wing for 2-d airfoil (unitless)
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
plot(alpha.*(180/pi),C_L_tot)
xlabel('Alpha (deg)')
ylabel('C_L')
title('NACA 4412 3-D lift curve slope @ Vmax, 10k')
