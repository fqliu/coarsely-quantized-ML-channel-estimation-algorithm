function gamma = func_newton(T_rr, T_lr, T_ri, T_li, gamma_pre, Temp_R, Temp_I)
%% Newton's method to solve the problem in Eq. (8)
dif = 1;
iter = 0;
while(dif > 1e-3)
       
    CC1 = T_rr * gamma_pre - Temp_R;
    CC2 = T_lr * gamma_pre - Temp_R;
    CC3 = T_ri * gamma_pre - Temp_I;
    CC4 = T_li * gamma_pre - Temp_I;
    
    
    temp1 = normpdf(CC1);
    temp2 = normpdf(CC2);
    
    temp3 = normpdf(CC3);
    temp4 = normpdf(CC4);
    
    temp5 = 0.5 * erfc(-1/sqrt(2) * CC1) - 0.5 * erfc(-1/sqrt(2) * CC2);
    temp6 = 0.5 * erfc(-1/sqrt(2) * CC3) - 0.5 * erfc(-1/sqrt(2) * CC4);
    temp7 = T_rr .* temp1 - T_lr .* temp2;
    temp8 = T_ri .* temp3 - T_li .* temp4;
    
    
    TT = temp7 ./ temp5 + temp8 ./ temp6;
    
    grad_1 = -sum(sum(TT));
    
    grad_2 = - sum(sum((-T_rr.^2 .* CC1 .* temp1 + T_lr.^2 .* CC2.*temp2) ./ temp5));
    
    
    grad_2 = grad_2 - sum(sum((-T_ri.^2 .* CC3 .* temp3 + T_li.^2 .* CC4 .* temp4)./temp6));
    grad_2 = grad_2 + sum(sum(TT.^2));
    
    gamma = gamma_pre - grad_1 / grad_2;
    
    dif = abs(gamma - gamma_pre) / abs(gamma_pre);
    
    gamma_pre = gamma;
    iter = iter +1;
end

    
    
    