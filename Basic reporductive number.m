days = 1:10; % Time period (1 to 10 days)
beta_time_dependent = [0.2 0.25 0.3 0.35 0.4 0.45 0.5 0.55 0.6 0.65]; % Time-dependent transmission rate beta
gamma_time_dependent = [0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.66 0.1 0.1]; % Time-dependent recovery rate gamma

% Initialize parameters
N = 1000; % Total population (hypothetical)
S0 = N - 10; % Initial susceptible population
I0 = 10; % Initial infected population
R0 = 0; % Initial recovered population

% Time span
tspan = days;

% Function for SIR model with time-dependent beta and gamma
function dydt = sir_model(t, y, beta, gamma)
    S = y(1);
    I = y(2);
    R = y(3);
    N = 1000; % Total population
    dSdt = -beta(t) * S * I / N;
    dIdt = beta(t) * S * I / N - gamma(t) * I;
    dRdt = gamma(t) * I;
    dydt = [dSdt; dIdt; dRdt];
end

% Calculate R0 for each day
R_values = zeros(1, 10); % Array to store R0 values for each day

for day = 1:10
    % Define beta and gamma as functions of time
    beta_function = @(t) beta_time_dependent(day);
    gamma_function = @(t) gamma_time_dependent(day);
    
    % Solve ODE for current day
    [t, y] = ode45(@(t, y) sir_model(t, y, beta_function, gamma_function), tspan, [S0; I0; R0]);
    
    % Calculate R0 using final values
    S_final = y(end, 1);
    I_final = y(end, 2);
    R_final = y(end, 3);
    R_values(day) = beta_time_dependent(day) / gamma_time_dependent(day);
end

% Display R0 values
disp('Basic Reproductive Ratio (R0) for each day:');
disp(R_values);
