function [Ke, D_Ke_D_E] = get_Ke(E, delta_E, x, L, theta)
    %% This function calculates the bar's element stiffness matrix
    %
    %% Input:
    %       E: The elasticity modulus
    %       delta_E: The random perturbation on E
    %       x: The sectional area
    %       L: The bar length
    %       theta: The rotation angles of the bars
    %
    %% Output:
    %       Ke: The element stiffness matrix
    %       D_Ke_D_E: The 1st order derivative of Ke with respect to the random elasticity modulus E

    %% Local stiffness matrix
    Ke_local = (E + delta_E) * x / L * [...
                            1, 0, -1, 0; ...
                            0, 0, 0, 0; ...
                            - 1, 0, 1, 0; ...
                            0, 0, 0, 0];

    D_Ke_D_E_local = x / L * [...
            1, 0, -1, 0; ...
            0, 0, 0, 0; ...
            - 1, 0, 1, 0; ...
            0, 0, 0, 0];

    %% The transformation matrix
    T_1 = [...
            cos(theta), sin(theta); ...
            -sin(theta), cos(theta)];

    T_4 = T_1;

    T = blkdiag(T_1, T_4);

    %% Calculate the Ke and its derivative
    Ke = T.' * Ke_local * T;
    D_Ke_D_E = T.' * D_Ke_D_E_local * T;

end
