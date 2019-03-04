function CD_airfoil = getAirfoilCoeffs(af_num, alpha, airfoil_drag)
    CD_airfoil = zeros(1,length(alpha)); %Allocate Space
    
    for i = 1:length(CD_airfoil) %Loop through coeff. vector
        for j = 1:size(airfoil_drag,1) %Loop through sheet
            if (alpha < airfoil_drag(1, j)*pi/180) %Find the right alpha
                break;
            end
        end
        CD_airfoil(i) = airfoil_drag(af_num, j); %Add it to the vector
    end
end