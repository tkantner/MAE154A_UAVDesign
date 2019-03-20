function CD_airfoil = getAirfoilCoeffs(af_num, alpha, airfoil_drag)
    CD_airfoil = zeros(1,length(alpha)); %Allocate Space
    
    for i = 1:length(CD_airfoil) %Loop through coeff. vector
        for j = 1:size(airfoil_drag,1) %Loop through sheet
            if (alpha(i) < airfoil_drag(j, 1)*pi/180) %Find the right alpha
                break;
            end
        end

        CD_airfoil(i) = airfoil_drag(j, af_num + 1);
    end
    
    %Fix if the AOA is too high for lower speeds at SL
    k = 1;
    while(isnan(CD_airfoil(k)))
        k = k + 1;
    end
    
    for l = 1:(k - 1)
        CD_airfoil(l) = CD_airfoil(k);
    end
end
