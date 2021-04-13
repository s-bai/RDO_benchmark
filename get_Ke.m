function [Ke, D_Ke_D_E] = get_Ke(E, delta_E, x, L, theta, kk)
    %% This function calculates the bar's element stiffness matrix
    %
    %% Input:
    %       E: The elasticity modulus
    %       delta_E: The random perturbation on E
    %       x: The sectional area
    %       L: The bar length
    %       theta: The rotation angles of the bars
    %       kk: The bar sequence number
    %
    %% Output:
    %       Ke: The element stiffness matrix
    %       D_Ke_D_E: The 1st order derivative of Ke with respect to the random elasticity modulus E

    %% Calculate the elasticity modulus
    if kk <= 6
        Ee = E(1) + delta_E(1);
    else
        Ee = E(2) + delta_E(2);
    end

    %% Local stiffness matrix (in the local coordinate system)
    Ke_local = Ee * x(kk) / L(kk) * [...
                                    1, 0, -1, 0; ...
                                    0, 0, 0, 0; ...
                                    - 1, 0, 1, 0; ...
                                    0, 0, 0, 0];

    %% The derivative of Ke_local with respect to the perturbation of the elasticity moduli
    D_Ke_D_E_local = zeros(4, 4, 2);

    if kk <= 6
        D_Ke_D_E_local(:, :, 1) = x(kk) / L(kk) * [...
                                            1, 0, -1, 0; ...
                                            0, 0, 0, 0; ...
                                            - 1, 0, 1, 0; ...
                                            0, 0, 0, 0];
    else
        D_Ke_D_E_local(:, :, 2) = x(kk) / L(kk) * [...
                                            1, 0, -1, 0; ...
                                            0, 0, 0, 0; ...
                                            - 1, 0, 1, 0; ...
                                            0, 0, 0, 0];
    end

    %% The transformation matrix
    T_1 = [...
            cos(theta(kk)), sin(theta(kk)); ...
            -sin(theta(kk)), cos(theta(kk))];

    T_4 = T_1;

    T = blkdiag(T_1, T_4);

    %% Calculate the Ke and its derivative (in the global coordinate system)
    Ke = T.' * Ke_local * T;

    D_Ke_D_E = zeros(4, 4, 2);

    for ii = 1:2
        D_Ke_D_E(:, :, ii) = T.' * D_Ke_D_E_local(:, :, ii) * T;
    end

end
