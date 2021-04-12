function rdo_benchmark()
    %% This function implements the RDO of a 2D three-bar truss

    %% Load input data from the JSON file
    JSON_data = jsondecode(fileread('input_data.json'));

    % The initial design variables (transformed to column vector)
    x = (JSON_data.x0).';

    % The bar lengths
    L = (JSON_data.L).';

    % The nominal E
    E = (JSON_data.E).';

    % The DOFs of the loads
    F_DOFs = (JSON_data.F_DOFs).';

    % The values of the loads
    F_value = (JSON_data.F_value).';

    %% The rotation angles of the bars
    theta = [0, 0, 0, 0, pi / 2, pi / 2, pi / 4, pi * 3/4, pi / 4, pi * 3/4];

    %% The local-to-global nodal information
    node_info = JSON_data.node_info;

    % x = [0.003, 0.003, 0.003];
    E = [2e6, 2e6, 2e6];
    F = [100, 100];

    %% The deviations of the elasticity moduli
    E_deviation = [2e5; 2e5; 2e5];

    %% The weighting factor of the RDO objective
    obj_beta = 1;

    %% The optimization parameters for the GA
    A = [5, 5, 7.2];
    b = 0.06;
    Aeq = [];
    beq = [];
    nonlcon = [];

    lb = [0.001 0.001 0.001];
    ub = [0.01 0.01 0.01];

    objective_fun = @(x) get_RDO_obj(x, E, E_deviation, F, obj_beta);

    options = optimoptions('ga', 'PlotFcn', @gaplotbestf);

    %% Implement the RDO using the GA
    x_GA = ga(objective_fun, 3, A, b, Aeq, beq, lb, ub, nonlcon, options);

    %%
    disp(x_GA);

end

%% ------------------------------------------------- Local functions -------------------------------------------------
function [u, D1_u] = get_U(x, E, F)
    %% This local function calculates the concerned displacement
    %
    % Input:
    %       x: The design variables vector (the sectional areas)
    %       E: The elasticity moduli vector
    %       F: The loads vector
    %
    % Output:
    %       u: The value of the concerned displacement
    %       D1_u: The first derivative of the concerned displacement with respect to the random parameters

    %%
    u = (F(1) * (-2.604166667 * x(1) * E(1) + 2.604166667 * x(2) * E(2) + 1.708984375 * x(3) * E(3)) + ...
        F(2) * (1.953125 * x(1) * E(1) + 1.953125 * x(2) * E(2) + 2.604166667 * x(3) * E(3))) / ...
        (x(1) * x(2) * E(1) * E(2) + 0.7434895833 * x(1) * x(3) * E(1) * E(3) + 0.08723958333 * x(2) * x(3) * E(2) * E(3) + 0.001708984375 * x(3)^2 * E(3)^2);

    if nargout > 1
        D1_u = zeros(3, 1);

        D1_u(1) = (-2.604166667 * F(1) * x(1) + 1.953125 * F(2) * x(1)) / (x(1) * x(2) * E(1) * E(2) + 0.7434895833 * x(1) * x(3) * E(1) * E(3) + 0.08723958333 * x(2) * x(3) * E(2) * E(3) + 0.001708984375 * x(3)^2 * E(3)^2) - ...
            ((x(1) * x(2) * E(2) + 0.7434895833 * x(1) * x(3) * E(3)) * (F(1) * (-2.604166667 * x(1) * E(1) + 2.604166667 * x(2) * E(2) + 1.708984375 * x(3) * E(3)) + F(2) * (1.953125 * x(1) * E(1) + 1.953125 * x(2) * E(2) + 2.604166667 * x(3) * E(3)))) / ...
            (x(1) * x(2) * E(1) * E(2) + 0.7434895833 * x(1) * x(3) * E(1) * E(3) + 0.08723958333 * x(2) * x(3) * E(2) * E(3) + 0.001708984375 * x(3)^2 * E(3)^2)^2;

        D1_u(2) = (2.604166667 * F(1) * x(2) + 1.953125 * F(2) * x(2)) / (x(1) * x(2) * E(1) * E(2) + 0.7434895833 * x(1) * x(3) * E(1) * E(3) + 0.08723958333 * x(2) * x(3) * E(2) * E(3) + 0.001708984375 * x(3)^2 * E(3)^2) - ...
            ((x(1) * x(2) * E(1) + 0.08723958333 * x(2) * x(3) * E(3)) * (F(1) * (-2.604166667 * x(1) * E(1) + 2.604166667 * x(2) * E(2) + 1.708984375 * x(3) * E(3)) + F(2) * (1.953125 * x(1) * E(1) + 1.953125 * x(2) * E(2) + 2.604166667 * x(3) * E(3)))) / ...
            (x(1) * x(2) * E(1) * E(2) + 0.7434895833 * x(1) * x(3) * E(1) * E(3) + 0.08723958333 * x(2) * x(3) * E(2) * E(3) + 0.001708984375 * x(3)^2 * E(3)^2)^2;

        D1_u(3) = (1.708984375 * F(1) * x(3) + 2.604166667 * F(2) * x(3)) / (x(1) * x(2) * E(1) * E(2) + 0.7434895833 * x(1) * x(3) * E(1) * E(3) + 0.08723958333 * x(2) * x(3) * E(2) * E(3) + 0.001708984375 * x(3)^2 * E(3)^2) - ...
            ((0.7434895833 * x(1) * x(3) * E(1) + 0.08723958333 * x(2) * x(3) * E(2) + 0.00341796875 * x(3)^2 * E(3)) * (F(1) * (-2.604166667 * x(1) * E(1) + 2.604166667 * x(2) * E(2) + 1.708984375 * x(3) * E(3)) + F(2) * (1.953125 * x(1) * E(1) + 1.953125 * x(2) * E(2) + 2.604166667 * x(3) * E(3)))) / ...
            (x(1) * x(2) * E(1) * E(2) + 0.7434895833 * x(1) * x(3) * E(1) * E(3) + 0.08723958333 * x(2) * x(3) * E(2) * E(3) + 0.001708984375 * x(3)^2 * E(3)^2)^2;
    end

end

function RDO_obj = get_RDO_obj(x, E, E_deviation, F, obj_beta)
    %% This local function calculates the objective of the RDO

    %% Calculate the objective function
    [u, D1_u] = get_U(x, E, F);

    u_mean = u;

    u_deviation = sum(D1_u .* D1_u .* E_deviation);

    %%
    RDO_obj = u_mean + obj_beta * u_deviation;

end
