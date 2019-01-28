close all; clear; clc;

W_S=10:40;%wing loading
C_D0=0.03; %parasite drag coefficient
V_max_10k=180*5280/(60*60); % ft/s
V_max_SL=130*5280/(60*60); %ft/s
V_min_10k=100*5280/(60*60);
V_min_SL=60*5280/(60*60);
rho_SL=0.00238; %sea level density (slug/ft^3)
rho_service=0.001756; %(slug/ft^3)
AR=7; %aspect ratio
e=0.7; %oswald efficiency

%Retrieve data from Spreadsheet
designs = xlsread('Master_Valid_Designs.xlsx');
[r, c] = size(designs); %Get size

C_L_SL=(2.*W_S)./(rho_service.*(V_max_SL.^2)); %
C_L_10k=(2.*W_S)./(rho_service.*(V_max_10k.^2));
C_D_SL=C_D0+(C_L_SL.^2)./(pi.*AR.*e); %unitless drag
C_D_10k=C_D0+(C_L_10k.^2)./(pi.*AR.*e);
RC_max=1500/60; %rate of climb @ SL (fps)
T_W_10k=0.5.*rho_service.*(V_max_10k.^2).*C_D_10k./(W_S); %thrust/weight required at 10k ft @V=180 mph
T_W_stall=0.5.*rho_service.*(V_min_10k.^2).*C_D_10k./(W_S);
T_W_SL=RC_max./((2./rho_SL).*(W_S).*sqrt(1./(pi.*AR.*e.*3.*C_D0)))+4*sqrt(C_D0/(pi*AR*e*3)); %thrust required at sea level at max RC
T_W_stall_SL=0.5.*rho_SL.*(V_min_SL.^2).*C_D_SL./(W_S);

figure(1)
plot(W_S,T_W_10k);
hold on;
plot(W_S,T_W_stall);
plot(W_S,T_W_SL);
plot(W_S,T_W_stall_SL);
xlabel('W/S (lbf/ft^2)'); ylabel('T/W (lbf/lbf)');
legend('Vmax, 10k ft','Stall speed @ 10k ft','Max RC @ Sea level', 'stall speed at seal level');
title('Required Thrust/weight against wing loading');
grid on;

% Power curve, un-normalized
eta_p=0.8;
P_W_10k=(T_W_10k.*V_max_10k)./(eta_p.*550); %power/weight required at Vmax, 5k feet (BHP)
P_W_stall=(T_W_stall.*V_min_10k)./(eta_p.*550); %power/weight required at cruise @10k ft (BHP)
P_W_SL=(T_W_SL.*V_max_SL)./(eta_p.*550); %power/weight required at max RC, sea level, Vcruise (BHP)
P_W_stall_SL=(T_W_stall_SL.*V_min_SL)./(eta_p.*550);

figure(2)
plot(W_S,P_W_10k);
hold on;
plot(W_S,P_W_stall);
plot(W_S,P_W_SL);
plot(W_S,P_W_stall_SL);

%Plot Individual data points
for i = 1:r
    if(i == r - 1)
        plot(designs(i, 23), designs(i, 24), 's', 'MarkerEdgeColor','red',...
            'MarkerFaceColor', 'red', 'MarkerSize', 10);
    elseif(i == r - 2)
        plot(designs(i, 23), designs(i, 24), 'x', 'MarkerEdgeColor','blue',...
            'MarkerFaceColor', 'blue', 'MarkerSize', 10);
    else
        plot(designs(i, 23), designs(i, 24), '.', 'MarkerEdgeColor','green',...
            'MarkerFaceColor', 'green', 'MarkerSize', 10);
    end
end

xlabel('W/S (lbf/ft^2)'); ylabel('P/W (BHP/lbf)');
legend('Vmax, 10k ft,','Stall speed @ 10k', 'Max RC @ sea level', 'stall speed sea level');
title('Required power/weight against wing loading');
grid on;

% Normalized power
sigma_service=rho_service/rho_SL;
P_W_10k_norm=P_W_10k/(1.132*sigma_service-0.132);
P_W_stall_norm=P_W_stall/(1.132*sigma_service-0.132);

figure(3)
plot(W_S,P_W_10k_norm);
hold on;
plot(W_S,P_W_stall_norm);
plot(W_S,P_W_SL);
plot(W_S,P_W_stall_SL);

%Plot Individual data points
for i = 1:r
    if(i == r - 1)
        plot(designs(i, 23), designs(i, 24), 's', 'MarkerEdgeColor','red',...
            'MarkerFaceColor', 'red', 'MarkerSize', 10);
    elseif(i == r - 2)
        plot(designs(i, 23), designs(i, 24), 'x', 'MarkerEdgeColor','blue',...
            'MarkerFaceColor', 'blue', 'MarkerSize', 10);
    else
        plot(designs(i, 23), designs(i, 24), '.', 'MarkerEdgeColor','green',...
            'MarkerFaceColor', 'green', 'MarkerSize', 10);
    end
end

xlabel('W/S (lbf/ft^2)'); ylabel('P_required/W (BHP/lbf)');
legend('Vmax, 10k ft','Stall speed @ 10k', 'Max RC @ sea level','stall speed sea level');
title('Normalized power @ sea level');

grid on;
