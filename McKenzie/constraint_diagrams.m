close all; clear; clc;
W_S=10:40;%wing loading
C_D0=0.04; %parasite drag coefficient
V_max=180*5280/(60*60); % ft/s
V_cruise=140*5280/(60*60); %ft/s
rho_SL=0.00238; %sea level density (slug/ft^3)
rho_service=0.001756; (slug/ft^3)
rho_5k=0.002048;
AR=7; %aspect ratio
e=0.7; %oswald efficiency
C_L=(2.*W_S)./(rho_SL.*(V_max.^2)); %
C_D=C_D0+(C_L.^2)./(pi.*AR.*e); %unitless drag
RC_max=50; %rate of climb (fps)
RC_service=1000/60; % (fps)
T_W=0.5.*rho_5k.*(V_max.^2).*C_D./(W_S); %thrust/weight required at 5k ft
T_W_service=RC_service./((2./0.001756).*(W_S).*sqrt(1./(pi.*AR.*e.*3.*C_D0)))+4*sqrt(C_D0/(pi*AR*e*3)); %thrust required at service ceiling, RC=1000fpm
T_W_cruise=0.5.*rho_service.*(V_cruise^2).*C_D0.*(1./W_S)+(1/(0.5*rho_service*(V_cruise^2)*pi*AR*e))*W_S; %thrust required at cruise @10k
T_W_SL=RC_max./((2./rho_SL).*(W_S).*sqrt(1./(pi.*AR.*e.*3.*C_D0)))+4*sqrt(C_D0/(pi*AR*e*3)); %thrust required at sea level at max RC
figure
plot(W_S,T_W)
hold on;
plot(W_S,T_W_service)
plot(W_S,T_W_cruise)
plot(W_S,T_W_SL)
xlabel('W/S (lbf/ft^2)'); ylabel('T/W (lbf/lbf)');
legend('Vmax, 5k ft','10k ft, RC=1k fpm','Vcruise, 10k ft','SL, RC=3k fpm')
title('Required Thrust/weight against wing loading')
grid on
%% Power curve, un-normalized
eta_p=0.8;
P_W=(T_W.*V_max)./(eta_p.*550); %power/weight required at Vmax, 5k feet (BHP)
P_W_service=(T_W_service.*V_max)./(eta_p.*550); %power/weight required at service ceiling & max velocity, with RC=1000fpm (BHP)
P_W_cruise=(T_W_cruise.*V_cruise)./(eta_p.*550); %power/weight required at cruise @10k ft (BHP)
P_W_SL=(T_W_SL.*V_cruise)./(eta_p.*550); %power/weight required at max RC, sea level, Vcruise (BHP)
figure
plot(W_S,P_W)
hold on;
plot(W_S,P_W_service)
plot(W_S,P_W_cruise)
plot(W_S,P_W_SL)
xlabel('W/S (lbf/ft^2)'); ylabel('P/W (BHP/lbf)');
legend('Vmax, 5k ft','Vmax, 10k ft, RC=1000 fpm','Vcruise, 10k ft', 'Vcruise, SL, RC=3000 fpm')
title('Required power/weight against wing loading')
grid on
%% Normalized power
sigma_5k=rho_5k/rho_SL;
sigma_service=rho_service/rho_SL;
P_W_norm=P_W/(1.132*sigma_5k-0.132);
P_W_service_norm=P_W_service/(1.132*sigma_service-0.132);
P_W_cruise_norm=P_W_cruise/(1.132*sigma_service-0.132);
figure
plot(W_S,P_W_norm)
hold on
plot(W_S,P_W_service_norm)
plot(W_S,P_W_cruise_norm)
plot(W_S,P_W_SL)
xlabel('W/S (lbf/ft^2)'); ylabel('P_required/W (BHP/lbf)');
legend('Vmax, 5k ft','Vmax, 10k ft, RC=1000 fpm','Vcruise, 10k ft', 'Vcruise, SL, RC=3000 fpm')
title('Normalized power @ sea level')
grid on




        


            

      

