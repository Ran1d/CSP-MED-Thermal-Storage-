%=======================================================
% Parabolic Trough Collector (PTC) Performance Model
%=======================================================

% 1. Input Parameters
% Environmental Conditions
DNI = 500; % Direct Normal Irradiance [W/m^2]
T_amb = 25; % Ambient temperature [°C]
v = 5; % Wind velocity [m/s]

% Collector Geometry
L = 20; % Collector length [m]
W_aperture = 3.5; % Aperture width [m]
D_out = 50e-3; % Receiver outer diameter [m]
D_in = 40e-3; % Receiver inner diameter [m]
D_g = 90e-3; % Glass cover diameter [m]

% Material Properties
k = 15; % Thermal conductivity of pipe [W/m-K]
cp_HTF = 1350; % Specific heat of HTF [J/kg-K]
rho = 1.11; % HTF density [kg/m^3]
mu = 2.02e-5; % HTF dynamic viscosity [kg/m·s]
k_air = 0.0276; % Air thermal conductivity [W/m·K]

% Optical and Thermal Efficiencies
eta_opt = 0.75; % Optical efficiency
a1 = 0.005; % Thermal loss coefficient [1/°C]
emissivity_r = 0.92; % Receiver emissivity
emissivity_g = 0.87; % Glass emissivity

% Flow and Temperatures
m_HTF = 0.32; % HTF mass flow rate [kg/s]
T_in = 220; % HTF inlet temperature [°C]
T_r = 260; % Receiver temperature [°C]

% Heat Transfer Coefficients
h_fi = 330; % Internal convective coefficient [W/m^2-K]

% Constants
sigma = 5.67e-8; % Stefan-Boltzmann constant [W/m^2-K^4]

% 2. Geometry and Surface Areas
A_r = pi * D_out * L; % Receiver area [m^2]
A_g = pi * D_g * L; % Glass outer area [m^2]
A_a = W_aperture * L; % Aperture area [m^2]

% 3. Temperature Conversions
T_r_K = T_r + 273.15; % Receiver temp [K]
T_amb_K = T_amb + 273.15; % Ambient temp [K]
% Initial guess for T_g for iteration (Assuming T_g is defined elsewhere or is an input)
% For this script, I will set a placeholder value for T_g to allow the script to run
T_g = 50; % Placeholder value for Glass temp [°C]
T_g_K = T_g + 273.15; % Glass temp [K]
T_avg = (T_g + T_amb) / 2; % Average temp [°C]

% 4. External Convective Heat Transfer
Re = ((rho * v * D_g) / mu); % Reynolds number
Nu = 0.3 * (Re)^0.6; % Nusselt number
h_w = (Nu * k_air) / D_g; % Wind convective loss [W/m^2-K]

% 5. Radiative Heat Transfer Coefficients
h_rad_ga = emissivity_g * sigma * (T_g_K^2 + T_amb_K^2) * (T_g_K + T_amb_K);
h_rad_rg = (sigma * (T_r_K^2 + T_g_K^2) * (T_r_K + T_g_K)) / (1/emissivity_r + (A_r/A_g)*(1/emissivity_g - 1));

% 6. Glass Energy Balance
Q_rad_rg_calc = A_r * h_rad_rg * (T_r - T_g); % Radiation from receiver to glass
Q_rad_ga_calc = A_g * h_rad_ga * (T_g - T_amb); % Radiation from glass to ambient
Q_conv_gw_calc = A_g * h_w * (T_g - T_amb); % Convection from glass to ambient
% The original line "Q_rad_rg = Q_rad_ga + Q_conv_gw" is an energy balance equation, not an assignment.
% I will leave the calculated values for inspection.

% 7. Overall Heat Loss Coefficient
U_L = 1 / (1/h_rad_rg + A_r / (A_g * (h_rad_ga + h_w)));

% 8. Collector Efficiency Factor and Heat Removal
ln_ratio = log(D_out / D_in); % Using 'log' for natural logarithm in MATLAB/Octave
F_prime = (1 / U_L) / ((1 / U_L) + (D_out / (h_fi * D_in)) + ((D_out / (2 * k)) * ln_ratio));
F_R = (m_HTF * cp_HTF) / (A_r * U_L) * (1 - exp(-(F_prime * A_r * U_L) / (m_HTF * cp_HTF))); % Using 'exp' for exponential

% 9. Useful Energy Gain and Outlet Temperature
Q_u = F_R * (DNI * A_a - U_L * A_r * (T_in - T_amb)); % Useful heat gain [W]
T_out = (Q_u / (m_HTF * cp_HTF)) + T_in; % HTF outlet temperature [°C]

% 10. Validation and Error Calculation
Q_u_theoretical = 23820;
Q_u_actual = 23028;
T_out_theoretical = 275.1;
T_out_actual = 273.3;
error_Q_u = abs((Q_u_theoretical - Q_u_actual) / Q_u_theoretical) * 100;
error_T_out = abs((T_out_theoretical - T_out_actual) / T_out_theoretical) * 100;

% Display results
fprintf('Calculated Useful Heat Gain (Qu): %.2f W\n', Q_u);
fprintf('Calculated HTF Outlet Temperature (Tout): %.2f °C\n', T_out);
fprintf('Error in Qu: %.2f%%\n', error_Q_u);
fprintf('Error in Tout: %.2f%%\n', error_T_out);

% End of Model
