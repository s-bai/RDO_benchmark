function rdo_benchmark()
    %% This function implements the RDO of a 2D three-bar truss
    % Note:
    %       Remember to set all one dimensional data arrays as column vector.

    %% Load input data from the JSON file
    JSON_data = jsondecode(fileread('input_data.json'));

    % The initial design variables (transformed to column vector)
    % x = JSON_data.x0;

    % The lower and upper limits of the design variables
    x_min = JSON_data.x_min;
    x_max = JSON_data.x_max;

    % The bar lengths
    L = JSON_data.L;

    % The nominal E
    E = JSON_data.E;

    % The DOFs of the loads
    F_DOFs = JSON_data.F_DOFs;

    % The values of the loads
    F_values = JSON_data.F_values;

    % Calculate the loads vector
    F = zeros(12, 1);
    F(F_DOFs) = F_values;

    % The fixed nodes
    fixed_nodes = JSON_data.fixed_nodes;

    %% The rotation angles of the bars
    theta = [0, 0, 0, 0, pi / 2, pi / 2, pi / 4, pi * 3/4, pi / 4, pi * 3/4];

    %% The nodes information
    node_info = JSON_data.node_info;

    % The fixed DOFs
    fixed_DOFs = sort(cat(1, 2 * fixed_nodes - 1, 2 * fixed_nodes));
    BCs.free_DOFs = setdiff(1:1:12, fixed_DOFs);

    % The DOF of the concerned displacement
    i_U = JSON_data.i_U;

    %% The standard deviations of the elasticity moduli
    E_sigma = JSON_data.E_sigma;

    % The nominal perturbation on the elasticity moduli (zero mean)
    delta_E = 0 * E_sigma;

    %% The upper limit on the volume
    V_max = JSON_data.V_max;

    %% The weighting factor of the RDO objective
    obj_beta = 1;

    %% The optimization parameters for the GA
    A = L.';
    b = V_max;
    Aeq = [];
    beq = [];
    nonlcon = [];

    lb = x_min;
    ub = x_max;

    objective_fun = @(x) get_RDO_obj(BCs, x, E, delta_E, L, theta, node_info, F, E_sigma, i_U, obj_beta);

    options = optimoptions('ga', 'PlotFcn', @gaplotbestf);

    %% Implement the RDO using the GA
    x_GA = ga(objective_fun, 10, A, b, Aeq, beq, lb, ub, nonlcon, options);

    %%
    disp('   1          2          3          4          5          6          7          8          9          10');
    disp(x_GA);

    %%
    [~, Ui_mean, Ui_sigma] = get_RDO_obj(BCs, x_GA, E, delta_E, L, theta, node_info, F, E_sigma, i_U, obj_beta);

    fprintf('Ui_mean = %7.2f; Ui_sigma = %7.2f\n', Ui_mean, Ui_sigma);

end

%% ------------------------------------------------- Local functions -------------------------------------------------
function [RDO_obj, Ui_mean, Ui_sigma] = get_RDO_obj(BCs, x, E, delta_E, L, theta, node_info, F, E_sigma, i_U, obj_beta)
    %% This function calculates the objective of the RDO
    % Input:
    %       See the corresponding comments in the main function and the subroutines
    % Output:
    %       RDO_obj: The objective value of the RDO

    %%
    [U, D_U_D_E] = FEA(BCs, x, E, delta_E, L, theta, node_info, F);

    %%
    Ui_mean = U(i_U);

    %%
    D_Ui_D_E = D_U_D_E(i_U, :);

    Ui_square_sigma = D_Ui_D_E .* D_Ui_D_E * E_sigma.^2;

    %%
    RDO_obj = Ui_mean + obj_beta * Ui_square_sigma;

    %%
    Ui_sigma = sqrt(Ui_square_sigma);
end
