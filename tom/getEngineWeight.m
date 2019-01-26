%This Function takes a power requirement (in hp) and engine table as input and returns the
%struct index of the engine that meets the power requirement
function index = getEngineWeight(Power_hp, engines)
   
    %If less than the smallest engine, just return that one
    if(Power_hp < engines(1,1))
        index = 1;
        return;
    end

    %Loop through and find the minimum engine and return the index
    for i = 1:(size(engines,1) - 1)
        if(Power_hp > engines(i,1) && Power_hp <= engines(i + 1,1))
            index = i + 1;
            return;
        end
        
    end
    
    % None meet this return 0 and check in main script
    index = 0;
    return;
end
