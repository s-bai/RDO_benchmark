function [U, D_U_D_E] = FEA(BCs, x, E, delta_E, L, theta, node_info, F)
    %% This function implements the FEA of the truss
    %
    %% Input:
    %       BCs: The boundary conditions
    %       x: The sectional areas
    %       E: The elasticity moduli
    %       delta_E: The perturbation on E
    %       L: The bar lengths
    %       theta: The rotation angles of the bars
    %       node_info: The local-to-global nodal information
    %       F: The loads vector
    %
    %% Output:
    %       U: The displacements
    %       D_U_D_E: The 1st derivative of U with respect to the random perturbation of E

    %%
    free_DOFs = BCs.free_DOFs;

    %%
    K = zeros(12, 12);
    K_temp = zeros(12, 12);

    D_K_D_E = zeros(12, 12, 2);
    D_K_D_E_temp = zeros(12, 12);

    U = zeros(12, 1);

    %%
    node_DOFs_info = zeros(10, 4);

    for ii = 1:2
        node_DOFs_info(:, 2 * ii - 1) = 2 * node_info(:, ii) - 1;
        node_DOFs_info(:, 2 * ii) = 2 * node_info(:, ii);
    end

    %%
    for ii = 1:10

        [Ke, D_Ke_D_E] = get_Ke(E, delta_E, x, L, theta, ii);

        row = repmat(node_DOFs_info(ii, :), 1, 4);
        column = kron(node_DOFs_info(ii, :), [1, 1, 1, 1]);

        idx = sub2ind(size(K), row, column);

        K_temp(idx) = Ke(:);

        K = K + K_temp;

        K_temp = 0 * K_temp;

        for jj = 1:2
            D_Ke_D_E_temp = D_Ke_D_E(:, :, jj);

            D_K_D_E_temp(idx) = D_Ke_D_E_temp(:);
            D_K_D_E(:, :, jj) = D_K_D_E(:, :, jj) + D_K_D_E_temp;

            D_K_D_E_temp = 0 * D_K_D_E_temp;
        end

    end

    %%
    U(free_DOFs) = K(free_DOFs, free_DOFs) \ F(free_DOFs);

    %%
    D_U_D_E = zeros(12, 2);

    for ii = 1:2
        D_U_D_E(free_DOFs, ii) = -inv(K(free_DOFs, free_DOFs)) * D_K_D_E(free_DOFs, free_DOFs, ii) * U(free_DOFs);
    end

end
