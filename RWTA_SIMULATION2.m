% Define Node coordinates
nodes = [100, 500;   % Node 0
         300, 1000;  % Node 1
         250, 300;   % Node 2
         400, 550;   % Node 3
         500, 500;   % Node 4
         550, 200;   % Node 5
         600, 550;   % Node 6
         700, 700;   % Node 7
         750, 1100;  % Node 8
         800, 800;   % Node 9
         700, 100;   % Node 10
         1200, 1250;% Node 11
         1000, 150;  % Node 12
         1300, 700]; % Node 13

% Define Edges (connections between nodes) with distances
edges = [1, 2, 1100;    % Edge between Node 0 and Node 1
         2, 3, 1600;    % Intermediate Edge between Node 1 and Node 2
         1, 3, 600;     % Edge between Node 0 and Node 2
         2, 8, 1500;    % Edge between Node 1 and Node 7
         1, 4, 1000;    % Edge between Node 0 and Node 3
         4, 5, 600;     % Edge between Node 3 and Node 4
         5, 7, 800;     % Intermediate edge between Node 4 and Node 6
         7, 8, 700;     % Intermediate edge between Node 6 and Node 7
         4, 9, 1500;    % Edge Between Node 3 and Node 8
         9, 12, 800;    % Edge Between Node 8 and Node 11
         8, 10, 700;    % Edge Between Node 7 and Node 9
         10, 12, 500;   % Edge Between Node 9 and Node 11
         9, 14, 800;    % Edge Between Node 8 and Node 13
         10, 14, 500;   % Edge Between Node 9 and Node 13
         12, 13, 300;   % Edge Between Node 11 and Node 12
         13, 14, 300;   % Edge Between Node 12 and Node 13
         3, 6, 1000;    % Edge Between Node 2 and Node 5
         6, 5, 1100;    % Edge Between Node 5 and Node 4
         6, 13, 2000;   % Edge Between Node 5 and Node 12
         6, 11, 1200;   % Edge Between Node 5 and Node 10
         10, 11, 900];  % Edge Between Node 9 and Node 10


% Generate node identifiers
num_nodes = size(nodes, 1);
node_names = cell(num_nodes, 1);
for i = 1:num_nodes
    node_names{i} = ['Node ', num2str(i-1)]; % Adjust index to start from 0
end

% Create graph using node identifiers and edge information
G = graph(edges(:, 1), edges(:, 2), edges(:, 3), node_names);

% Optionally, you can set the node coordinates as node attributes
G.Nodes.Coordinates = nodes;

% Initialize Wavelengths attribute for edges
G.Edges.Wavelengths = zeros(size(edges, 1), 1); % Assuming all wavelengths are initially set to 0

% Adjust source and destination nodes
source_node = 'Node 0';   % Adjusted source node
destination_node = 'Node 13';  % Adjusted destination node

% Find the shortest path using Dijkstra's algorithm
shortest_path = shortestpath(G, source_node, destination_node);

% Display the shortest path
disp('Shortest path:');
disp(shortest_path);

% Define parameters
t = 1; % Time for each request to complete key transmission and synchronization

% Case 1
n_values_case1 = 4;
Tr_values_case1 = [6*t, 8*t, 10*t, 12*t];

% Case 2
n_values_case2 = [2, 4, 6];
Tr_values_case2 = 9*t;

% Tr values for Solution 2
Tr_values_solution2 = {[8*t, 8*t, 10*t, 10*t], [6*t, 8*t, 10*t, 12*t]}; 

num_requests = [100, 200, 300]; % Number of connection requests
m = 40; % Number of wavelengths for TDCh

% Preallocate variables for case 1 (Solution 1 - Case 1)
wavelength_utilization_TDCh_case1 = zeros(1, length(Tr_values_case1), length(num_requests));
time_slot_utilization_case1 = zeros(1, length(Tr_values_case1), length(num_requests));
CSR_case1 = zeros(1, length(Tr_values_case1), length(num_requests));

% Preallocate variables for case 2 (Solution 1 - Case 2)
wavelength_utilization_TDCh_case2 = zeros(length(n_values_case2), length(Tr_values_case2), length(num_requests));
time_slot_utilization_case2 = zeros(length(n_values_case2), length(Tr_values_case2), length(num_requests));
CSR_case2 = zeros(length(n_values_case2), length(Tr_values_case2), length(num_requests));

% Preallocate variables for solution 2
num_methods_solution2 = length(Tr_values_solution2);
wavelength_utilization_TDCh_solution2 = zeros(num_methods_solution2, length(Tr_values_solution2{1}), length(num_requests));
time_slot_utilization_solution2 = zeros(num_methods_solution2, length(Tr_values_solution2{1}), length(num_requests));
CSR_solution2 = zeros(num_methods_solution2, length(Tr_values_solution2{1}), length(num_requests));

% Simulation for all cases and solutions
for num_idx = 1:length(num_requests)
    num = num_requests(num_idx);
    
    % Simulation for Solution 1 - Case 1
    for Tr_index = 1:length(Tr_values_case1)
        Tr = Tr_values_case1(Tr_index);
        disp(['Solution 1 - Case 1 - Tr:', num2str(Tr)]);
        [wavelength_utilization_TDCh_case1(1, Tr_index, num_idx),...
         time_slot_utilization_case1(1, Tr_index, num_idx),...
         CSR_case1(1, Tr_index, num_idx)] = Simulation(n_values_case1, Tr, num, m);
    end

    % Simulation for Solution 1 - Case 2
    for n_index = 1:length(n_values_case2)
        n = n_values_case2(n_index);
        disp(['Solution 1 - Case 2 - n:', num2str(n)]);
        [wavelength_utilization_TDCh_case2(n_index, 1, num_idx),...
         time_slot_utilization_case2(n_index, 1, num_idx),...
         CSR_case2(n_index, 1, num_idx)] = Simulation(n, Tr_values_case2, num, m);
    end

    % Simulation for Solution 2 - Method 1
    Tr_method1 = [8*t, 8*t, 10*t, 10*t];
    for Tr_index = 1:length(Tr_method1)
        Tr = Tr_method1(Tr_index);
        disp(['Solution 2 - Method 1 - Tr:', num2str(Tr)]);
        [wavelength_utilization_TDCh_solution2(1, Tr_index, num_idx), ...
         time_slot_utilization_solution2(1, Tr_index, num_idx), ...
         CSR_solution2(1, Tr_index, num_idx)] = Simulation(4, Tr, num, m);
    end

    % Simulation for Solution 2 - Method 2
    Tr_method2 = [6*t, 8*t, 10*t, 12*t];
    for Tr_index = 1:length(Tr_method2)
        Tr = Tr_method2(Tr_index);
        disp(['Solution 2 - Method 2 - Tr:', num2str(Tr)]);
        [wavelength_utilization_TDCh_solution2(2, Tr_index, num_idx), ...
         time_slot_utilization_solution2(2, Tr_index, num_idx), ...
         CSR_solution2(2, Tr_index, num_idx)] = Simulation(4, Tr, num, m);
    end
end 

% Plotting simulation results
figure; % Create a new figure

% Wavelength Utilization
subplot(1, 3, 1);
plot_utilization_histograms(Tr_values_case1, wavelength_utilization_TDCh_case1, 'Wavelength Utilization (Case 1)');
hold on;
plot_utilization_histograms(Tr_values_case2, wavelength_utilization_TDCh_case2, 'Wavelength Utilization (Case 2)');
plot_solution_histograms(Tr_values_solution2, wavelength_utilization_TDCh_solution2, 'Wavelength Utilization (Solution 2)');
xlabel('Number of Connection Requests (Tr)');
ylabel('Wavelength Utilization');
legend('Case 1', 'Case 2', 'Solution 2 - Method 1', 'Solution 2 - Method 2', 'Location', 'best');
grid on;

% Time-Slot Utilization
subplot(1, 3, 2);
plot_utilization_histograms(Tr_values_case1, time_slot_utilization_case1, 'Time-Slot Utilization (Case 1)');
hold on;
plot_utilization_histograms(Tr_values_case2, time_slot_utilization_case2, 'Time-Slot Utilization (Case 2)');
plot_solution_histograms(Tr_values_solution2, time_slot_utilization_solution2, 'Time-Slot Utilization (Solution 2)');
xlabel('Number of Connection Requests (Tr)');
ylabel('Time-Slot Utilization');
legend('Case 1', 'Case 2', 'Solution 2 - Method 1', 'Solution 2 - Method 2', 'Location', 'best');
grid on;

% CSR (Connection Security Ratio)
subplot(1, 3, 3);
plot_csr_histograms(Tr_values_case1, CSR_case1, 'CSR (Case 1)');
hold on;
plot_csr_histograms(Tr_values_case2, CSR_case2, 'CSR (Case 2)');
plot_solution_histograms(Tr_values_solution2, CSR_solution2, 'CSR (Solution 2)');
xlabel('Number of Connection Requests (Tr)');
ylabel('CSR');
legend('Case 1', 'Case 2', 'Solution 2 - Method 1', 'Solution 2 - Method 2', 'Location', 'best');
grid on;

function plot_utilization_histograms(Tr_values, utilization_data, title_str)
    % Plot histograms for wavelength or time-slot utilization
    for Tr_index = 1:numel(Tr_values)
        if size(utilization_data, 2) >= Tr_index
            histogram(utilization_data(:, Tr_index), 'DisplayName', ['Tr = ', num2str(Tr_values(Tr_index))]);
            hold on; % Add hold on to overlay histograms
        else
            disp(['Warning: Not enough data for Tr_index = ', num2str(Tr_index)]);
        end
    end
    title(title_str);
end

function plot_csr_histograms(Tr_values, csr_data, title_str)
    % Plot CSR histograms
    for Tr_index = 1:numel(Tr_values)
        plot(Tr_values, csr_data(:, Tr_index), '-o', 'DisplayName', ['Tr = ', num2str(Tr_values(Tr_index))]);
        hold on; % Add hold on to overlay plots
    end
    title(title_str);
end

function plot_solution_histograms(Tr_values, utilization_data, title_str)
    % Plot histograms for wavelength or time-slot utilization for Solution 2
    for Tr_index = 1:numel(Tr_values)
        if size(utilization_data, 2) >= Tr_index
            histogram(utilization_data(:, Tr_index), 'DisplayName', ['Tr = ', num2str(Tr_values(Tr_index))]);
            hold on; % Add hold on to overlay histograms
        else
            disp(['Warning: Not enough data for Tr_index = ', num2str(Tr_index)]);
        end
    end
    title(title_str);
end

% Simulation function definition
function [wavelength_utilization, time_slot_utilization, CSR] = Simulation(n, Tr, num_requests, m)
    disp(['Simulation - Tr:', num2str(Tr), ', Num:', num2str(num_requests)]);
    
    % Initialize matrices for results
    wavelength_utilization = zeros(1, 1, num_requests);
    time_slot_utilization = zeros(1, 1, num_requests);
    CSR = zeros(1, 1, num_requests);

    disp(['Size of wavelength_utilization before assignment: ', num2str(size(wavelength_utilization))]);
    disp(['Size of time_slot_utilization before assignment: ', num2str(size(time_slot_utilization))]);
    disp(['Size of CSR before assignment: ', num2str(size(CSR))]);

    % Loop over number of requests
    for req = 1:num_requests
        % Define the maximum expected traffic load per request
        max_traffic_load = 300; % Adjust this value based on your network's characteristics

        % Generate a random traffic volume for the request
        D = randi([1, max_traffic_load]);

        % Perform wavelength allocation using FirstFit algorithm
        wavelength_assigned = FirstFit(m, n);

        % Calculate time slot utilization
        time_slot_utilization(1, 1, req) = D / Tr;

        % Calculate Connection Security Ratio (CSR)
        if wavelength_assigned <= n
            CSR(1, 1, req) = (m - wavelength_assigned) / m;
        else
            CSR(1, 1, req) = 0; % Request not allocated a secure wavelength
        end

        % Update wavelength utilization
        wavelength_utilization(1, 1, req) = wavelength_assigned / m;
    end
    
    % Ensure that the matrices have the correct size
    if ~isequal(size(wavelength_utilization), [1, 1, num_requests]) || ...
       ~isequal(size(time_slot_utilization), [1, 1, num_requests]) || ...
       ~isequal(size(CSR), [1, 1, num_requests])
        error('Invalid size of output matrices');
    end 
end

% FirstFit algorithm for wavelength allocation
function wavelength_assigned = FirstFit(m, n)
    % Assume FirstFit algorithm without considering network topology
    % Assign the first available wavelength if it satisfies the security requirements
    for w = 1:m
        if w <= n
            wavelength_assigned = w;
            return;
        end
    end
    % If no wavelength satisfies the security requirements, assign 0
    wavelength_assigned = 0;
end