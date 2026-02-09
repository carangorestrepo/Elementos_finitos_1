
      dydx(v_)   = y(t_)-y(V_)/Ac;      % Dy = t-V/Ac
      dydx(t_)   = y(M_)/(EI);          % Dt = M/EI
      dydx(M_)   = y(V_) - (P)*dydx(v_);% DM = V-P*Dy
      dydx(V_)   = 0-k*y(v_);           % Dv = qyloc
      dydx(u_)   = y(fax_)/(AE);        % Du = P/AE
      dydx(fax_) = 0;                   % DP = 0
      
      

      dydx(v_)   = dydx(M_)/Ac-y(t_);   % Dy = DM/Ac-t (5A)
      dydx(t_)   = y(M_)/(EI);          % Dt = M/EI (4A)
      dydx(M_)   = y(V_) + (P)*dydx(v_);% DM = V + P*Dy (A1)
      dydx(V_)   = 0-k*y(v_);           % DV = qyloc  (A3)
      dydx(u_)   = y(fax_)/(AE);        % Du = P/AE
      dydx(fax_) = 0;                   % DP = 0
      
      
      
      
      
      
      

      
      
      
      
      