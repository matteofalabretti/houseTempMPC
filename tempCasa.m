function x_dot = tempCasa(t, x, k, C, tau, T_ext, k_ext, u)

k_esima = zeros(3,3);
x_dot = zeros(6,1);

for i = 1:3
    for j = 1:3
        if i ~= j  
            k_esima(i,j) = k(i) + 4/(1+ exp(-0.5 * sqrt((x(i) - x(j))^2)));
        end
    end
    x_dot(i+3) = (u(i+3) - x(i+3))/tau(i);
end

x_dot(1) = (x(4) - k_esima(1,2) * (x(1) - x(2)) - k_esima(1,3) * (x(1) - x(3)) - k_ext *(x(1) - T_ext))/C(1) ;
x_dot(2) = (x(5) - k_esima(1,2) * (x(1) - x(2)) - k_esima(2,3) * (x(2) - x(3)) - k_ext *(x(2) - T_ext))/C(2);
x_dot(3) = (x(6) - k_esima(1,3) * (x(1) - x(3)) + k_esima(2,3) * (x(2) - x(3)) - k_ext *(x(3) - T_ext))/C(3);

end