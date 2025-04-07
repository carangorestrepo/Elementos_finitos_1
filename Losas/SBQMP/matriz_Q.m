function Q = matriz_Q(xi, eta)
    % xi, eta: Coordenadas naturales del punto de integración
    Q = zeros(5, 12); % 5 deformaciones × 12 coeficientes ?
    
    Q(1, 4) = 1;       % ?_x = ?4 + ?5*y + ?7*(x/2)
    Q(1, 5) = eta;      % y ? ? en coordenadas naturales
    Q(1, 7) = xi/2;     % x ? ?
    
    Q(2, 6) = 1;        % ?_y = ?6 + ?7*x + ?5*(y/2)
    Q(2, 7) = xi;
    Q(2, 5) = eta/2;
    
    Q(3, 8) = 1;        % ?_xy = ?8 + 2?5*x + 2?7*y
    Q(3, 5) = 2*xi;
    Q(3, 7) = 2*eta;
    
    Q(4, 9) = 1;        % ?_xz = ?9 + ?10*y
    Q(4, 10) = eta;
    
    Q(5, 11) = 1;       % ?_yz = ?11 + ?12*x
    Q(5, 12) = xi;
end