function T = matriz_T(rho, h)
    T = rho * diag([h, h^3/12, h^3/12]); % Para W, ?x, ?y
end