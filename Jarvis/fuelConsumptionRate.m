% Outputs the fuel consumption rate [lb/hr] given a power requirement input [hp]
function FC_rate = fuelConsumptionRate(power_required)
min_power = 0.3; %[hp]
max_power = 15; %[hp]
power = min_power:0.01:max_power; % Power Range [hp]

% Returns error message if input not within min/max range
if power_required < min_power || power_required > max_power
    error('Power Required must be between %d and %d hp',min_power,max_power)
end

for i=1:length(power)
    %% 2 Stroke Gas Engine Fuel Consumption [lb/hr] %%
    if power(i) >= 2
        FC_2gas(i) = 1.3887*power(i) - 1.3954;
    else
        FC_2gas(i) = 1E5; % Acts as +inf if its out of range
    end
    %% 4 Stroke Gas Engine Fuel Consumption [lb/hr] %%
    if power(i) >= 2.5
        FC_4gas(i) = 0.6768*power(i) - 0.1258;
    else
        FC_4gas(i) = 1E5;
    end
    %% 2 Stroke Glow Engine Fuel Consumption [lb/hr] %%
    if power(i) >= 0.3 && power(i) <= 4
        FC_2glow(i) = 2.3665*power(i) + 0.9131;
    else
        FC_2glow(i) = 1E5;
    end
    %% 4 Stroke Glow Engine Fuel Consumption [lb/hr] %%
    if power(i) >= 0.8 && power(i) <= 5
        FC_4glow(i) = 4.9123*power(i) - 1.8983;
    else
        FC_4glow(i) = 1E5;
    end
    %% Finds the lowest FC at a specific power value for the function %%
    FC_min(i) = FC_2gas(i);
    if FC_min(i) > FC_4gas(i)
        FC_min(i) = FC_4gas(i);
    end
    if FC_min(i) > FC_2glow(i)
        FC_min(i) = FC_2glow(i);
    end
    if FC_min(i) > FC_4glow(i)
        FC_min(i) = FC_4glow(i);
    end
    if round(power(i),2) == round(power_required,2)
        FC_rate = FC_min(i);
    end
end
%{
%% Plot for illustration %%
figure(1);
plot(power,FC_2gas,'o',power,FC_4gas,'o',power,FC_2glow,'o',power,FC_4glow,'o',power,FC_min,'+');
legend('2 Stroke Gas','4 Stroke Gas','2 Stroke Glow','4 Stroke Glow','Min');
xlabel('Power Required [hp]');
ylabel('Fuel Consumption [lb/hr]');
ylim([0 25]);
%}
end
