


function L_k=L_ka(k,xi,xj,xn,yn,xm,ym)    
    %          Figures 3 and 5
    if k == 5
        L_k=((xj - xi)^2 + (yj - yi)^2)^0.5;
    elseif k == 6
         L_k=((xm- xj)^2 + (ym - yj)^2)^0.5;
    elseif k == 7
         L_k=((xn - xm)^2 + (yn - ym)^2)^0.5;
    elseif k == 8
        L_k=((xi - xn)^2 + (yi - yn)^2)^0.5;
    end
end

function [C,S]=dir_cos( k,xi,xj,xn,yn,xm,ym)
    %% Figures 3 and 5
    if k == 5
        C = (xj - xi)/L_k;
        S = (yj - yi)/L_k;
    elseif k == 6
        C = (xm- xj)/L_k;
        S = (ym - yj)/L_k;
    elseif k == 7
        C = (xn - xm)/L_k;
        S = (yn - ym)/L_k;
    elseif k == 8
        C = (xi - xn)/L_k;
        S = (yi - yn)/L_k;
    end
end
    
function phi_k=phi_ka(k,nu,t,L_k)
        kappa = 5/6;
        % Equation 74
        phi_k= 2/(kappa*(1-nu))*(t/L_k(k))^2;
end

function N_i=N_ia(i, xi, eta)
        %%'''
        %%Returns the interpolation function for any given coordinate in the natural coordinate system
        %%'''
        if i == 1
            N_i=1/4*(1-eta)*(1-xi);
        elseif i == 2
            N_i=1/4*(1+eta)*(1-xi);
        elseif i == 3
            N_i= 1/4*(1+eta)*(1+xi);
        elseif i == 4
            N_i=1/4*(1-eta)*(1+xi);
        end
end
function P_k=P_ka(k, xi, eta)
        if k == 5
            P_k=1/2*(1-eta^2)*(1-xi);
        elseif k == 6
            P_k=1/2*(1+eta)*(1-xi^2);
        elseif k == 7
            P_k=1/2*(1-eta^2)*(1+xi);
        elseif k == 8
            P_k=1/2*(1-eta)*(1-xi^2);
        end
end

function J=Ja(xi, eta,x1, y1, x2, y2, x3, y3, x4, y4)
        %%'''
        %%Returns the Jacobian matrix for the element
        %%'''
        
        %% Get the local coordinates for the element
        %x1, y1, x2, y2, x3, y3, x4, y4 = x1, y1, x2, y2, x3, y3, x4, y4

        %% Return the Jacobian matrix
        J=1/4*[[x1*(eta - 1) - x2*(eta - 1) + x3*(eta + 1) - x4*(eta + 1), y1*(eta - 1) - y2*(eta - 1) + y3*(eta + 1) - y4*(eta + 1)];
               [x1*(xi - 1)  - x2*(xi + 1)  + x3*(xi + 1)  - x4*(xi - 1),  y1*(xi - 1)  - y2*(xi + 1)  + y3*(xi + 1)  - y4*(xi - 1)]];
end

function N_gamma=N_gammaa(xi, eta)

        %% Equation 44
        N_gamma=[[1/2*(1-eta),        0,1/2*(1+eta),      0];
                 [          0,1/2*(1+xi),         0,      1/2*(1-xi)]];
end
function A_gamma=A_gammaa(L5,L6,L7,L8)

        %L5 = L_k(5);
        %L6 = L_k(6);
        %L7 = L_k(7);
        %L8 = L_k(8);

        %% Equation 46
        A_gamma=[[L5/2,   0,     0,     0];
                 [  0,  L6/2,    0,     0];
                 [  0,    0,  -L7/2,    0];
                 [  0,    0,     0,  -L8/2]];
end


function Aun=A_u(L5,L6,L7,L8,C5,C6,C7,C8)
        %L5 = L_k(5);
        %L6 = L_k(6);
        %L7 = L_k(7);
        %L8 = L_k(8);

        %[C5, S5] = dir_cos(5);
        %[C6, S6] = dir_cos(6);
        %[C7, S7] = dir_cos(7);
        [%C8, S8] = dir_cos(8);

        A_un=[[-2/L5, C5, S5,  2/L5, C5, S5,   0,    0,  0,   0,    0,  0];
             [  0,    0,  0, -2/L6, C6, S6, 2/L6,  C6, S6,   0,    0,  0];
             [  0,    0,  0,   0,    0,  0, -2/L7, C7, S7,  2/L7, C7, S7];
             [ 2/L8, C8, S8,   0,    0,  0,   0,    0,  0, -2/L8, C8, S8]];
end

function A_Delta_inv_DKMQ(self)

        phi5 = phi_k(5);
        phi6 = phi_k(6);
        phi7 = phi_k(7);
        phi8 = phi_k(8);

        A_Delta_inv_DKMQ=-3/2 * [[1/(1+phi5),     0,          0,          0     ];
                                 [   0,       1/(1+phi6),     0,          0     ];
                                 [   0,           0,      1/(1+phi7),     0     ];
                                 [   0,           0,          0,      1/(1+phi8)]];
end
 

function  A_phi_Delta(self)

        phi5 = phi_k(5)
        phi6 = phi_k(6)
        phi7 = phi_k(7)
        phi8 = phi_k(8)

        A_phi_Delta=[[phi5/(1+phi5),       0,             0,             0      ];
                    [      0,       phi6/(1+phi6),       0,             0      ];
                    [      0,             0,       phi7/(1+phi7),       0      ];
                    [      0,             0,             0,       phi8/(1+phi8)]];
end
function  B_b_beta( xi, eta)

        %% Get the inverse of the Jacobian matrix
        J_inv = inv(J(xi, eta))

        %% Get the individual terms for the Jacobian inverse
        j11 = J_inv(0, 0)
        j12 = J_inv(0, 1)
        j21 = J_inv(1, 0)
        j22 = J_inv(1, 1)

        %% Derivatives of the bilinear interpolation functions
        N1_xi = 0.25*eta - 0.25
        N2_xi = 0.25 - 0.25*eta
        N3_xi = 0.25*eta + 0.25
        N4_xi = - 0.25*eta - 0.25
        N1_eta = 0.25*xi - 0.25
        N2_eta = - 0.25*xi - 0.25
        N3_eta = 0.25*xi + 0.25
        N4_eta = 0.25 - 0.25*xi

        N1x = j11*N1_xi + j12*N1_eta
        N1y = j21*N1_xi + j22*N1_eta
        N2x = j11*N2_xi + j12*N2_eta
        N2y = j21*N2_xi + j22*N2_eta
        N3x = j11*N3_xi + j12*N3_eta
        N3y = j21*N3_xi + j22*N3_eta
        N4x = j11*N4_xi + j12*N4_eta
        N4y = j21*N4_xi + j22*N4_eta

        B_b_beta=[[0, N1x,  0,  0, N2x,  0,  0, N3x,  0,  0, N4x,  0 ];
                  [0,  0,  N1y, 0,  0,  N2y, 0,  0,  N3y, 0,  0,  N4y];
                  [0, N1y, N1x, 0, N2y, N2x, 0, N3y, N3x, 0, N4y, N4x]];
end
    
function  B_b_Delta_beta=B_b_Delta_betaa( xi, eta)

        %% Get the inverse of the Jacobian matrix
        J_inv = inv(J(xi, eta))

        %% Get the individual terms for the Jacobian inverse
        j11 = J_inv(0, 0)
        j12 = J_inv(0, 1)
        j21 = J_inv(1, 0)
        j22 = J_inv(1, 1)

        %% Derivatives of the quadratic interpolation functions
        P5_xi = xi*(eta - 1)
        P6_xi = -0.5*(eta - 1)*(eta + 1)
        P7_xi = -xi*(eta + 1)
        P8_xi = 0.5*(eta - 1)*(eta + 1)
        P5_eta = 0.5*(xi - 1)*(xi + 1)
        P6_eta = -eta*(xi + 1)
        P7_eta = -0.5*(xi - 1)*(xi + 1)
        P8_eta = eta*(xi - 1)

        P5x = j11*P5_xi + j12*P5_eta
        P5y = j21*P5_xi + j22*P5_eta
        P6x = j11*P6_xi + j12*P6_eta
        P6y = j21*P6_xi + j22*P6_eta
        P7x = j11*P7_xi + j12*P7_eta
        P7y = j21*P7_xi + j22*P7_eta
        P8x = j11*P8_xi + j12*P8_eta
        P8y = j21*P8_xi + j22*P8_eta

        C5, S5 = dir_cos(5)
        C6, S6 = dir_cos(6)
        C7, S7 = dir_cos(7)
        C8, S8 = dir_cos(8)

        B_b_Delta_beta=[[    P5x*C5,          P6x*C6,          P7x*C7,          P8x*C8     ];
                        [    P5y*S5,          P6y*S6,          P7y*S7,          P8y*S8,    ];
                        [P5y*C5 + P5x*S5, P6y*C6 + P6x*S6, P7y*C7 + P7x*S7, P8y*C8 + P8x*S8]];
    
function B_b= B_ba( xi, eta)
        %%'''
        %%Returns the [B_b] matrix for bending
        %%'''

        %% Return the [B] matrix for bending
        B_b= add(B_b_beta(), B_b_Delta_beta()*A_Delta_inv_DKMQ()*Au())
end
function  B_m( r, s)

        %% Differentiate the interpolation functions
        %% Row 1 = interpolation functions differentiated with respect to x
        %% Row 2 = interpolation functions differentiated with respect to y
        %% Note that the inverse of the Jacobian converts from derivatives with
        %% respect to r and s to derivatives with respect to x and y
        dH = np.matmul(inv(J(r, s)), 1/4*[[s + 1, -s - 1, s - 1, -s + 1];                
                                          [r + 1, -r + 1, r - 1, -r - 1]]);

        % Reference 2, Example 5.5 (page 353)
        B_m =[[dH(0, 0),    0,     dH(0, 1),    0,     dH(0, 2),    0,     dH(0, 3),    0    ];
              [   0,     dH(1, 0),    0,     dH(1, 1),    0,     dH(1, 2),    0,     dH(1, 3)];
              [dH(1, 0), dH(0, 0), dH(1, 1), dH(0, 1), dH(1, 2), dH(0, 2), dH(1, 3), dH(0, 3)]];


end
function  Hb(self)
        %%'''
        %%Returns the stress-strain matrix for plate bending.
        %%'''

        %% Referemce 1, Table 4.3, page 194
        nu = nu
        E = E
        h = t

        Hb = E*h^3/(12*(1 - nu^2))*[[1,  nu,      0    ];
                                    [nu, 1,       0    ];
                                    [0,  0,  (1 - nu)/2]];
        
end

function Hs(self)
        %%'''
        %%Returns the stress-strain matrix for shear.
        %%'''
        %% Reference 2, Equations (5.97), page 422
        k = 5/6
        E = E
        h = t
        nu = nu

        Hs = E*h*k/(2*(1 + nu))*[[1, 0];
                                 [0, 1]];
end

function Cm(self)
        %%'''
        %%Returns the stress-strain matrix for an isotropic or orthotropic plane stress element
        %%'''
        
        % Apply the stiffness modification factors for each direction to obtain orthotropic
        % behavior. Stiffness modification factors of 1.0 in each direction (the default) will
        % model isotropic behavior. Orthotropic behavior is limited to the element's local
        % coordinate system.
        Ex = E*kx_mod
        Ey = E*ky_mod
        nu_xy = nu
        nu_yx = nu

        % The shear modulus will be unafected by orthotropic behavior
        % Logan, Appendix C.3, page 750
        G = E/(2*(1 + nu))

        % Gallagher, Equation 9.3, page 251
        Cm = 1/(1 - nu_xy*nu_yx)*np.array([[   Ex,    nu_yx*Ex,           0         ],
                                        [nu_xy*Ey,    Ey,              0         ],
                                        [    0,        0,     (1 - nu_xy*nu_yx)*G]])

end

function  k_b(self)
        %%'''
        %%Returns the local stiffness matrix for bending stresses
        %%'''

        Hb = Hb()
        Hs = Hs()

        % Define the gauss point for numerical integration
        gp = 1/3^0.5

        % Get the determinant of the Jacobian matrix for each gauss pointing. Doing this now will save us from doing it twice below.
        J1 = det(J(-gp, -gp))
        J2 = det(J(gp, -gp))
        J3 = det(J(gp, gp))
        J4 = det(J(-gp, gp))

        % Get the bending B matrices for each gauss point
        B1 = B_b(-gp, -gp)
        B2 = B_b(gp, -gp)
        B3 = B_b(gp, gp)
        B4 = B_b(-gp, gp)

        % Create the stiffness matrix with bending stiffness terms
        % See 2, Equation 5.94
        k = (np.matmul(B1.T, np.matmul(Hb, B1))*J1 +...
             np.matmul(B2.T, np.matmul(Hb, B2))*J2 +...
             np.matmul(B3.T, np.matmul(Hb, B3))*J3 +...
             np.matmul(B4.T, np.matmul(Hb, B4))*J4)

        % Get the shear [B_s] matrices for each gauss point
        B1 = B_s(-gp, -gp)
        B2 = B_s(gp, -gp)
        B3 = B_s(gp, gp)
        B4 = B_s(-gp, gp)

        % Add shear stiffness terms to the stiffness matrix
        %k += (np.matmul(B1.T, np.matmul(Hs, B1))*J1 +...
        %      np.matmul(B2.T, np.matmul(Hs, B2))*J2 +...
        %      np.matmul(B3.T, np.matmul(Hs, B3))*J3 +...
        %      np.matmul(B4.T, np.matmul(Hs, B4))*J4)...
        
        % Following Bathe's recommendation for the drilling degree of freedom
        % from Example 4.19 in "Finite Element Procedures, 2nd Ed.", calculate
        % the drilling stiffness as 1/1000 of the smallest diagonal term in
        % the element's stiffness matrix. This is not theoretically correct,
        % but it allows the model to solve without singularities, and should
        % have a minimal effect on the final solution. Bathe recommends 1/1000
        % as a value that is weak enough but not so small that it affect the
        % results. Bathe recommends looking at all the diagonals in the
        % combined bending plus membrane stiffness matrix. Some of those terms
        % relate to translational stiffness. It seems more rational to only
        % look at the terms relating to rotational stiffness. That will be
        % Pynite's approach.
        k_rz = min(abs(k(1, 1)), abs(k(2, 2)), abs(k(4, 4)), abs(k(5, 5)),...
                   abs(k(7, 7)), abs(k(8, 8)), abs(k(10, 10)), abs(k(11, 11))...
                   )/1000
        
        % Initialize the expanded stiffness matrix to all zeros
        k_exp = zeros(24, 24)

        % Step through each term in the unexpanded stiffness matrix
        % i = Unexpanded matrix row
        for i= 1: 12

            % j = Unexpanded matrix column
            for j =1: 12
                
                % Find the corresponding term in the expanded stiffness
                % matrix

                % m = Expanded matrix row
                reca=[0, 3, 6, 9];
                if i==reca   % indices associated with deflection in z
                    m = 2*i + 2
                end
                
                if i in [1, 4, 7, 10]  % indices associated with rotation about x
                    m = 2*i + 1
                end
                
                if i in [2, 5, 8, 11]  % indices associated with rotation about y
                    m = 2*i
                end

                % n = Expanded matrix column
                
                if j in [0, 3, 6, 9]  % indices associated with deflection in z
                    n = 2*j + 2
                end
                if j in [1, 4, 7, 10]  % indices associated with rotation about x
                    n = 2*j + 1
                end
                if j in [2, 5, 8, 11]  % indices associated with rotation about y
                    n = 2*j
                end
                
                % Ensure the indices are integers rather than floats
                m, n = round(m), round(n)

                % Add the term from the unexpanded matrix into the expanded
                % matrix
                k_exp[m, n] = k[i, j]
            end
        end

        % Add the drilling degree of freedom's weak spring
        k_exp[5, 5] = k_rz
        k_exp[11, 11] = k_rz
        k_exp[17, 17] = k_rz
        k_exp[23, 23] = k_rz


                
 end
%%%
function k_m(self)
        %%'''
        %%Returns the local stiffness matrix for membrane (in-plane) stresses.
        %Plane stress is assumed
        %%'''
        t = t
        Cm = Cm()

        % Define the gauss point for numerical integration
        gp = 1/3^0.5

        % Get the membrane B matrices for each gauss point
        % Doing this now will save us from doing it twice below
        B1 = B_m(gp, gp)
        B2 = B_m(-gp, gp)
        B3 = B_m(-gp, -gp)
        B4 = B_m(gp, -gp)

        % See reference 2 at the bottom of page 353, and reference 2 page 466
        k = t*(np.matmul(B1.T, np.matmul(Cm, B1))*det(J(gp, gp)) +
               np.matmul(B2.T, np.matmul(Cm, B2))*det(J(-gp, gp)) +
               np.matmul(B3.T, np.matmul(Cm, B3))*det(J(-gp, -gp)) +
               np.matmul(B4.T, np.matmul(Cm, B4))*det(J(gp, -gp)))
        k_exp = np.zeros((24, 24))

        % Step through each term in the unexpanded stiffness matrix
        % i = Unexpanded matrix row
        for i in range(8)

            % j = Unexpanded matrix column
            for j in range(8)
                
                % Find the corresponding term in the expanded stiffness
                % matrix

                % m = Expanded matrix row
                if i in [0, 2, 4, 6]  % indices associated with displacement in x
                    m = i*3
                end
                if i in [1, 3, 5, 7]  % indices associated with displacement in y
                    m = i*3 - 2
                end

                % n = Expanded matrix column
                if j in [0, 2, 4, 6]  % indices associated with displacement in x
                    n = j*3
                end
                if j in [1, 3, 5, 7]  % indices associated with displacement in y
                    n = j*3 - 2
                end
                
                % Ensure the indices are integers rather than floats
                m, n = round(m), round(n)

                % Add the term from the unexpanded matrix into the expanded matrix
                k_exp[m, n] = k[i, j]
            end
        end
end

%%%
    function k(self)
        %'''
        %Returns the quad element's local stiffness matrix.
        %'''

        % Recalculate the local coordinate system
        _local_coords()

        % Sum the bending and membrane stiffness matrices
        return np.add(k_b(), k_m())

%%%   
    function f( combo_name='Combo 1')
        %'''
        %Returns the quad element's local end force vector
        %'''
        
        % Calculate and return the plate's local end force vector
        return np.add(np.matmul(k(), d(combo_name)), fer(combo_name))

%%%
    function fer( combo_name='Combo 1')
        %'''
        %Returns the quadrilateral's local fixed end reaction vector.

        Parameters
        ----------
        combo_name  string
            The name of the load combination to get the consistent load vector for.
        %'''
        
        Hw = lambda r, s  1/4*np.array([[(1 - r)*(1 - s), 0, 0, (1 + r)*(1 - s), 0, 0, (1 + r)*(1 + s), 0, 0, (1 - r)*(1 + s), 0, 0]])

        % Initialize the fixed end reaction vector
        fer = np.zeros((12,1))

        % Get the requested load combination
        combo = model.LoadCombos[combo_name]

        % Define the gauss point used for numerical integration
        gp = 1/3^0.5

        % Initialize the element's surface pressure to zero
        p = 0

        % Loop through each load case and factor in the load combination
        for case, factor in combo.factors.items()

            % Sum the pressures
            for pressure in pressures

                % Check if the current pressure corresponds to the current load case
                if pressure[1] == case

                    % Sum the pressures
                    p -= factor*pressure[0]
        
        fer = (Hw(-gp, -gp).T*p*det(J(-gp, -gp))
             + Hw(gp, -gp).T*p*det(J(gp, -gp))
             + Hw(gp, gp).T*p*det(J(gp, gp))
             + Hw(-gp, gp).T*p*det(J(-gp, gp)))

        % Initialize the expanded vector to all zeros
        fer_exp = np.zeros((24, 1))

        % Step through each term in the unexpanded vector
        % i = Unexpanded vector row
        for i in range(12)
                
            % Find the corresponding term in the expanded vector

            % m = Expanded vector row
            if i in [0, 3, 6, 9]   % indices associated with deflection in z
                m = 2*i + 2
            if i in [1, 4, 7, 10]  % indices associated with rotation about x
                m = 2*i + 1
            if i in [2, 5, 8, 11]  % indices associated with rotation about y
                m = 2*i
                
            % Ensure the index is an integer rather than a float
            m = round(m)

            % Add the term from the unexpanded vector into the expanded vector
            fer_exp[m, 0] = fer[i, 0]

        return fer_exp

%%%
    def d( combo_name='Combo 1')
       %'''
       %Returns the quad element's local displacement vector
       %'''

       % Calculate and return the local displacement vector
       return np.matmul(T(), D(combo_name))

%%%
    def F( combo_name='Combo 1')
        %'''
        %Returns the quad element's global force vector

        Parameters
        ----------
        combo_name  string
            The load combination to get results for.
        %'''
        
        % Calculate and return the global force vector
        return np.matmul(inv(T()), f(combo_name))

%%%
    def D( combo_name='Combo 1')
        %'''
        %Returns the quad element's global displacement vector for the given
        load combination.

        Parameters
        ----------
        combo_name  string
            The name of the load combination to get the displacement vector
            for (not the load combination itself).
        %'''
        
        % Initialize the displacement vector
        D = np.zeros((24, 1))
        
        % Read in the global displacements from the nodes
        D[0, 0] = i_node.DX[combo_name]
        D[1, 0] = i_node.DY[combo_name]
        D[2, 0] = i_node.DZ[combo_name]
        D[3, 0] = i_node.RX[combo_name]
        D[4, 0] = i_node.RY[combo_name]
        D[5, 0] = i_node.RZ[combo_name]

        D[6, 0] = j_node.DX[combo_name]
        D[7, 0] = j_node.DY[combo_name]
        D[8, 0] = j_node.DZ[combo_name]
        D[9, 0] = j_node.RX[combo_name]
        D[10, 0] = j_node.RY[combo_name]
        D[11, 0] = j_node.RZ[combo_name]

        D[12, 0] = m_node.DX[combo_name]
        D[13, 0] = m_node.DY[combo_name]
        D[14, 0] = m_node.DZ[combo_name]
        D[15, 0] = m_node.RX[combo_name]
        D[16, 0] = m_node.RY[combo_name]
        D[17, 0] = m_node.RZ[combo_name]

        D[18, 0] = n_node.DX[combo_name]
        D[19, 0] = n_node.DY[combo_name]
        D[20, 0] = n_node.DZ[combo_name]
        D[21, 0] = n_node.RX[combo_name]
        D[22, 0] = n_node.RY[combo_name]
        D[23, 0] = n_node.RZ[combo_name]
        
        % Return the global displacement vector
        return D

%%%
    function K(self)
        %'''
        %Returns the quad element's global stiffness matrix
        %'''

        % Get the transformation matrix
        T = T()

        % Calculate and return the stiffness matrix in global coordinates
        return np.matmul(np.matmul(inv(T), k()), T)

%%% 
    % Global fixed end reaction vector
function FER( combo_name='Combo 1')
        %'''
        %Returns the global fixed end reaction vector.

        Parameters
        ----------
        combo_name  string
            The name of the load combination to calculate the fixed end
            reaction vector for (not the load combination itself).
        %'''
        
        % Calculate and return the fixed end reaction vector
        return np.matmul(inv(T()), fer(combo_name))
end

%%%  
function T(self)
        %'''
        %Returns the coordinate transformation matrix for the quad element.
        %'''

        xi = xi
        xj = xj
        yi = yi
        yj = yj
        zi = i_node.Z
        zj = j_node.Z

        % Calculate the direction cosines for the local x-axis.The local x-axis will run from
        % the i-node to the j-node
        x = [xj - xi, yj - yi, zj - zi]

        % Divide the vector by its magnitude to produce a unit x-vector of
        % direction cosines
        mag = (x[0]^2 + x[1]^2 + x[2]^2)^0.5
        x = [x[0]/mag, x[1]/mag, x[2]/mag]
        
        % The local y-axis will be in the plane of the plate. Find a vector in
        % the plate's local xy plane.
        xm = xm
        ym = ym
        zm = m_node.Z
        xy = [xm - xi, ym - yi, zm - zi]

        % Find a vector perpendicular to the plate surface to get the
        % orientation of the local z-axis.
        z = np.cross(x, xy)
        
        % Divide the z-vector by its magnitude to produce a unit z-vector of
        % direction cosines.
        mag = (z[0]^2 + z[1]^2 + z[2]^2)^0.5
        z = [z[0]/mag, z[1]/mag, z[2]/mag]

        % Calculate the local y-axis as a vector perpendicular to the local z
        % and x-axes.
        y = np.cross(z, x)
        
        % Divide the y-vector by its magnitude to produce a unit vector of
        % direction cosines.
        mag = (y[0]^2 + y[1]^2 + y[2]^2)^0.5
        y = [y[0]/mag, y[1]/mag, y[2]/mag]

        % Create the direction cosines matrix.
        dirCos = np.array([x,
                           y,
                           z])
        
        % Build the transformation matrix.
        T = np.zeros((24, 24))
        T[03, 03] = dirCos
        T[36, 36] = dirCos
        T[69, 69] = dirCos
        T[912, 912] = dirCos
        T[1215, 1215] = dirCos
        T[1518, 1518] = dirCos
        T[1821, 1821] = dirCos
        T[2124, 2124] = dirCos
        
        % Return the transformation matrix.
end

%%%
function shear( xi=0, eta=0, local=True, combo_name='Combo 1')
        %'''
        %Returns the interal shears at any point in the quad element.

        Internal shears are reported as a 2D array [[Qx], [Qy]] at the
        specified location in the (xi, eta) natural coordinate system.

        Parameters
        ----------
        xi  number
            The xi-coordinate. Default is 0.
        eta  number
            The eta-coordinate. Default is 0.
        
        %Returns
        -------
        Internal shear force per unit length of the quad element [[Qx], [Qy]]
        %'''

        % Get the plate's local displacement vector
        % Slice out terms not related to plate bending
        d = d(combo_name)[[2, 3, 4, 8, 9, 10, 14, 15, 16, 20, 21, 22], ]

        % Define the gauss point used for numerical integration
        gp = 1/3^0.5

        % Define extrapolated r and s points
        xi_ex = xi/gp
        eta_ex = eta/gp

        % Define the interpolation functions
        H = 1/4*np.array([(1 - xi_ex)*(1 - eta_ex), (1 + xi_ex)*(1 - eta_ex), (1 + xi_ex)*(1 + eta_ex), (1 - xi_ex)*(1 + eta_ex)])

        % Get the stress-strain matrix
        Hs = Hs()

        % Calculate the internal shears [Qx, Qy] at each gauss point
        q1 = np.matmul(Hs, np.matmul(B_s(-gp, -gp), d))
        q2 = np.matmul(Hs, np.matmul(B_s(gp, -gp), d))
        q3 = np.matmul(Hs, np.matmul(B_s(gp, gp), d))
        q4 = np.matmul(Hs, np.matmul(B_s(-gp, gp), d))

        % Extrapolate to get the value at the requested location
        Qx = H[0]*q1[0] + H[1]*q2[0] + H[2]*q3[0] + H[3]*q4[0]
        Qy = H[0]*q1[1] + H[1]*q2[1] + H[2]*q3[1] + H[3]*q4[1]

        if local
            
            return np.array([Qx,
                             Qy])
        
        else
            
            % Get the direction cosines for the plate's local coordinate system
            dir_cos = T()[3, 3]

            return np.matmul(dir_cos.T, np.array([Qx,
                                                  Qy,
                                                  [0]]))
        end
end

%%%   
function moment( xi=0, eta=0, local=True, combo_name='Combo 1')
        %'''
        %Returns the interal moments at any point in the quad element.

        %Internal moments are reported as a 2D array [[Mx], [My], [Mxy]] at the
        %specified location in the (xi, eta) natural coordinate system.

        %Parameters
        ----------
        xi  number
            The xi-coordinate. Default is 0.
        eta  number
            The eta-coordinate. Default is 0.
        
        %Returns
        -------
        Internal moment per unit length of the quad element [[Mx], [My], [Mxy]]
        %'''

        % Get the plate's local displacement vector
        % Slice out terms not related to plate bending
        d = d(combo_name)[[2, 3, 4, 8, 9, 10, 14, 15, 16, 20, 21, 22], ]

        % Define the gauss point used for numerical integration
        gp = 1/3^0.5

        % % Define extrapolated r and s points
        xi_ex = xi/gp
        eta_ex = eta/gp

        % Define the interpolation functions
        H = 1/4*np.array([(1 - xi_ex)*(1 - eta_ex), (1 + xi_ex)*(1 - eta_ex), (1 + xi_ex)*(1 + eta_ex), (1 - xi_ex)*(1 + eta_ex)])

        % Get the stress-strain matrix
        Hb = Hb()

        % Calculate the internal moments [Mx, My, Mxy] at each gauss point
        m1 = np.matmul(Hb, np.matmul(B_kappa(-gp, -gp), d))
        m2 = np.matmul(Hb, np.matmul(B_kappa(gp, -gp), d))
        m3 = np.matmul(Hb, np.matmul(B_kappa(gp, gp), d))
        m4 = np.matmul(Hb, np.matmul(B_kappa(-gp, gp), d))

        % Extrapolate to get the value at the requested location
        Mx = H[0]*m1[0] + H[1]*m2[0] + H[2]*m3[0] + H[3]*m4[0]
        My = H[0]*m1[1] + H[1]*m2[1] + H[2]*m3[1] + H[3]*m4[1]
        Mxy = H[0]*m1[2] + H[1]*m2[2] + H[2]*m3[2] + H[3]*m4[2]
        
        if local

            return np.array([Mx,
                             My,
                             Mxy])
        
        else
            
            % Get the direction cosines for the plate's local coordinate system
            dir_cos = T()[3, 3]

            % Convert the plate flexural stresses to global coordinates
            return np.matmul(dir_cos.T, np.array([Mx,
                                                  My,
                                                  [0]]))

%%%
    function membrane( xi=0, eta=0, local=True, combo_name='Combo 1')
        
        % Get the plate's local displacement vector
        % Slice out terms not related to membrane forces
        d = d(combo_name)[[0, 1, 6, 7, 12, 13, 18, 19], ]

        % Define the gauss point used for numerical integration
        gp = 1/3^0.5

        % Define extrapolated r and s points
        xi_ex = xi/gp
        eta_ex = eta/gp

        % Define the interpolation functions
        H = 1/4*np.array([(1 - xi_ex)*(1 - eta_ex), (1 + xi_ex)*(1 - eta_ex), (1 + xi_ex)*(1 + eta_ex), (1 - xi_ex)*(1 + eta_ex)])

        % Get the stress-strain matrix
        Cm = Cm()
        
        % Calculate the internal stresses [Sx, Sy, Txy] at each gauss point
        s1 = np.matmul(Cm, np.matmul(B_m(-gp, -gp), d))
        s2 = np.matmul(Cm, np.matmul(B_m(gp, -gp), d))
        s3 = np.matmul(Cm, np.matmul(B_m(gp, gp), d))
        s4 = np.matmul(Cm, np.matmul(B_m(-gp, gp), d))

        % Extrapolate to get the value at the requested location
        Sx = H[0]*s1[0] + H[1]*s2[0] + H[2]*s3[0] + H[3]*s4[0]
        Sy = H[0]*s1[1] + H[1]*s2[1] + H[2]*s3[1] + H[3]*s4[1]
        Txy = H[0]*s1[2] + H[1]*s2[2] + H[2]*s3[2] + H[3]*s4[2]

        if local

            return np.array([Sx,
                             Sy,
                             Txy])
        
        else

            % Get the direction cosines for the plate's local coordinate system
            dir_cos = T()[3, 3]

            % Convert the plate membrane stresses to global coordinates
            return np.matmul(dir_cos.T, np.array([Sx,
                                                  Sy,
                                                  [0]]))

