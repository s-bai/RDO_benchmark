function [U, D_U_D_E] = FEA(x, E, E_deviation, delta_E, L, theta, node_info)
    %% This function implements the FEA of the truss
    %
    %% Input:
    %       x: The sectional areas
    %       E: The elasticity moduli
    %       E_deviation: The deviations of E
    %       delta_E: The perturbation on E
    %       L: The bar lengths
    %       theta: The rotation angles of the bars
    %       node_info: The local-to-global nodal information
    %
    %% Output:
    %       U: The displacements
    %       D_U_D_E: The 1st derivative of U with respect to the random perturbation of E

    %%
    K = zeros(12,12);
    D_K_D_E = zeros(12,12);

    %%
    node_DOFs_info = zeros(10, 4);

    for ii = 1:2
        node_DOFs_info(:, 2 * ii - 1) = 2 * node_info(:, ii) - 1;
    end

    for ii = 1:2
        node_DOFs_info(:, 2 * ii) = 2 * node_info(:, ii);
    end

    %%
    for ii = 1:10

        if ii >= 1 && ii <= 6
            delta_E_temp = delta_E(1);
        else
            delta_E_temp = delta_E(2);
        end

        [Ke, D_Ke_D_E] = get_Ke(E(ii), delta_E_temp, x(ii), L(ii), theta(ii));

        row = repmat(node_DOFs_info(ii, :), 1, 4);
        column = kron(node_DOFs_info(ii, :), [1, 1, 1, 1]);

        idx = sub2ind(size(Ke), row, column);

        K(idx) = Ke(:);
        D_K_D_E(idx) = D_Ke_D_E(:);

    end

end
