% Simulation of cylindrical bore drainage with continuous lateral and bottom leaks
% using a simplified advection-dispersion approach and animation

clear 
clc
% Parameters
r = 0.35;         % Bore r (m)
h0 = 0.9;         % Bore height (m)
A_bottom = pi * r^2; % Bottom area (m^2)

h202 = h0;     % Initial water h0 (m)

dt = 1;              % Time step (s)
t_final = 3600;       % Final simulation time (s)

L_lateral = h0;      % Length of the lateral leak area (m)

k_lateral = 7e-4;    % Lateral leak rate constant (m^2.5/s) - Reduced
k_bottom = k_lateral;     % Bottom leak rate constant (m^2.5/s) - Reduced

g = 9.81;            % Acceleration due to gravity (m/s^2)

% Dispersion coefficient (simplified, constant)
D = 0.0001; % m^2/s

% Number of radial segments for simplified advection-dispersion
nr = 10;
dr = r / nr;
rad = dr/2:dr:r-dr/2; %radial positions

% Initialization
t = 0:dt:t_final;
h = zeros(size(t));
h(1) = h0;
conc = zeros(nr, length(t)); %concentration over the r
conc(:,1)=1; %initial concentration is 1

% Set up the figure and animation
figure;
axis equal;
xlim([-r r]);
ylim([0 (h0)]);
hold on;


% Initialize the water level patch (rectangle approximation for circle)
water_patch = patch([-r r r -r],[h(1) h(1) 0 0],'b');

% Add a text object for displaying the current time
time_text = text(-0.3*r, h0 - 0.3, 'Time: 0 s', 'FontSize', 12, 'Color', 'k');

hold off;
xlabel('X (m)');
ylabel('Y (m)');
title('Cross section of Cylindrical Bore Drainage Simulation (Animated)');
grid on;

% Simulation loop
for i = 1:length(t)-1
    % Lateral Leak Calculation (same principle, adjusted constants)
    Q_lateral_total = 0;
    h_int_limit = min(h(i), L_lateral);
    if h_int_limit > 0
        Q_lateral_total = k_lateral *2*pi*r*h(i);  %Q = -KA \frac{\Delta h}{\Delta l} 
    end

    % Bottom Leak Calculation (same principle, adjusted constants) 
    Q_bottom = k_bottom * pi*r*r*h(i); %Q = -A * K \frac{\Delta h}{\Delta l} 

    % Total outflow
    Q_total = Q_lateral_total + Q_bottom;

    % Calculate change in water height
    dh = - (Q_total / A_bottom) * dt;

    % Update water height
    h(i+1) = h(i) + dh;

    % Prevent negative water height
    if h(i+1) < 0
        h(i+1) = 0;
        break;
    end
    
    %Simplified Advection-Dispersion
    for j=1:nr
        %Advection (simplified as a fraction of the total outflow)
        adv = -(Q_total/(A_bottom*h0))*conc(j,i)*dt;
        %Diffusion (simplified radial diffusion)
        if j>1 && j<nr
            diff = D*(conc(j+1,i)-2*conc(j,i)+conc(j-1,i))/(dr^2)*dt;
        elseif j==1
            diff = D*(conc(j+1,i)-conc(j,i))/(dr^2)*dt; %no flux boundary at r=0
        elseif j==nr
            diff = D*(-conc(j,i)+conc(j-1,i))/(dr^2)*dt; %outflow boundary
        end
        conc(j,i+1) = conc(j,i) + adv + diff;
        if conc(j,i+1)<0
            conc(j,i+1)=0;
        end
    end

    % Update the plot for animation (simplified)
    set(water_patch, 'YData', [h(i+1) h(i+1) 0 0]);
    
    % Update the time text
    set(time_text, 'String', sprintf('Time: %d s', t(i+1)));
    
    drawnow;
    pause(0.005);
end
