% N-body gravity simulation with adaptive multilevel time stepping
% ------------------------------------------------------------------------
%
% Particle initial conditions file format:
% id mass x y z xdot ydot zdot
%
% Input:    filename =  Initial condition file: ex. '2body'
%           tmax     =  Runtime in seconds: ex. 3
%           n        =  Number of levels/bins: ex. 2
%           eps      =  Softening: ex. 0.0001
%           method   =  Integrator: ex. 'Euler', 'DKD', 'KDK', 'RK4'
%
% Output:      id    =  Particle ids
%              mass  =  Particle masses
%              pos   =  Position 3-vector matrix
%              vel   =  Velocity 3-vector matrix
%              dE_E  =  Energy difference
%              E     =  Energy components
%              t     =  Time in seconds
%         timesteps  =  Timestep sizes per evaluation
%
%
% - Performance tweaks (memory pre-allocation etc.) TBD
%
% Computational Astrophysics Course Project, ETH Zurich, 2012
% mikael.mieskolainen@cern.ch

function [id,mass,pos,vel, dE_E, E, t, timesteps] = nbody(filename, tmax, n, eps, method)

% Read the input file
Input = dlmread(filename);

% Number of particles
[N,~] = size(Input);

% Time vectors
t = zeros(1,1);
t(1) = 0;
count = 2;

id = Input(:,1);                    % Particle id vector
mass = Input(:,2);                  % Particle mass vector
dts = ones(size(mass));             % Individual time steps

% Position and velocity matrices
pos = zeros(N(1), 3, 1);
vel = zeros(N(1), 3, 1);            % 3rd dimension is iteration

pos(:,:,1) = Input(:,3:5,1);        % Initialize position matrix
vel(:,:,1) = Input(:,6:8,1);        % Initialize velocity matrix

% Initial energies
Ep = potential(pos(:,:,1),mass);    % Initial potential energy
Ek = kinetic(vel(:,:,1),mass);      % Initial kinetic energy
E_old = sum(Ep) + sum(Ek);

E(count-1,1) = E_old;
E(count-1,2) = sum(Ek);
E(count-1,3) = sum(Ep);

dt_n = zeros(n,1);
level_id = zeros(size(dts));    

% Timestep for plotting
timesteps = zeros(N,1);


% Integration loop
while (t(count-1) < tmax)
    
    
    % --------------------------------------------------------------------
    % Adaptive timesteps for each particle based on current positions
    %
    a = acc(pos(:,:,count-1), mass, eps);
    
    eta = 0.1; % Tuning factor
    for i = 1:length(dts)
        dts(i) = eta*(eps/norm(a(i,:)))^0.5;
    end
    
    % New timestep t_1
    dt1 = max(dts);
    
    % --------------------------------------------------------------------

    % Calculate dt for each bin (level)
    for i = 1:n
        dt_n(i) = dt1 / (2^(i-1));
    end
    
    % Find bin assignments
    for i = 1:length(dts)
        [~,level_id(i)] = min( abs(dt_n - repmat(dts(i), n, 1)));
    end
    
    % Save timesteps for debugging
    for i = 1:N
       timesteps(i,count-1) = dt_n(level_id(i)); 
    end
    
    fprintf('\n\n');
    fprintf('Particles in bins: ');
    for i = 1:n
       fprintf('[%d|%0.8f|%d] ', i, dt_n(i), sum(level_id == i));
    end
    fprintf('\n');
    
    % Now go through every level i = 1...end
    for i = 1:n
                    
        % Create time vector with zeros for non-updating particles
        dt = ones(size(dts))*dt_n(i);
        dt(level_id ~= i) = 0;
        dt = repmat(dt, 1, 3); % replicate to x,y,z components
        
        % Now go through every step of this level
        for j = 1:2^(i-1)
            
            % Do we have the final step
            if (j == 2^(i-1) && i == n), f = 0; else, f = 1; end
            
            % ------------------------------------------------------------
            % Forward Euler
            %
            if (strcmp(method, 'Euler'))
                
                % x_n+1
                pos(:,:,count-f) = pos(:,:,count-1) + vel(:,:,count-1).*(dt/2);
                
                % v_n+1
                vel(:,:,count-f) = vel(:,:,count-1) + acc(pos(:,:,count-1), mass, eps).*dt;
            end
            %}
            
            % ------------------------------------------------------------
            % DKD (drift-kick-drift) Leap frog
            %
            if (strcmp(method, 'DKD'))
                
                % x_n+1/2
                pos(:,:,count-f) = pos(:,:,count-1) + vel(:,:,count-1).*(dt/2);

                % v_n+1
                vel(:,:,count-f) = vel(:,:,count-1) + acc(pos(:,:,count-f), mass, eps).*dt;

                % x_n+1
                pos(:,:,count-f) = pos(:,:,count-f) + vel(:,:,count-f).*(dt/2);
            end
            %}

            % ------------------------------------------------------------
            % KDK (kick-drift-kick) Leap frog
            %
            if (strcmp(method, 'KDK'))
                
                % v_n+1/2
                vel(:,:,count-f) = vel(:,:,count-1) + acc(pos(:,:,count-1),mass,eps).*(dt/2);

                % x_n+1
                pos(:,:,count-f) = pos(:,:,count-1) + vel(:,:,count-f).*dt;

                % v_n+1
                vel(:,:,count-f) = vel(:,:,count-f) + acc(pos(:,:,count-f),mass,eps).*(dt/2);
            end
            %}

            % ------------------------------------------------------------
            % Runge-Kutta 4th order (RK4)
            %
            if (strcmp(method, 'RK4'))

                kr1 = vel(:,:,count-1);

                kv1 = acc(pos(:,:,count-1), mass, eps);
                kr2 = vel(:,:,count-1) + kv1 .* (dt/2);

                kv2 = acc(pos(:,:,count-1) + kr2 .* (dt/2), mass, eps);
                kr3 = vel(:,:,count-1) + kv2 .* (dt/2);

                kv3 = acc(pos(:,:,count-1) + kr3 .* (dt/2), mass, eps);
                kr4 = vel(:,:,count-1) + kv3 .* dt;

                kv4 = acc(pos(:,:,count-1) + kr4 .* dt, mass, eps);

                % v_n+1
                vel(:,:,count-f) = vel(:,:,count-1) + (dt/6) .* (kv1 + 2*(kv2 + kv3) + kv4);
                
                % x_n+1
                pos(:,:,count-f) = pos(:,:,count-1) + (dt/6) .* (kr1 + 2*(kr2 + kr3) + kr4);

            end
            % ------------------------------------------------------------
            %}
            
        end
    end
    
    % Calculate energies
    Ep = potential(pos(:,:,count),mass);
    Ek = kinetic(vel(:,:,count),mass);
    E_tot = sum(Ep) + sum(Ek);
    
    % Energy difference
    dE = (E_tot - E_old) / E_old;
    E_old = E_tot;
    
    % Energy vectors for debugging/analysis etc.
    dE_E(count) = dE;
    E(count,1) = E_tot;
    E(count,2) = sum(Ek);
    E(count,3) = sum(Ep);
    %}
    
    % Update counter, and global time with the time step
    t(count) = t(count-1) + dt_n(1);
    count = count + 1;

    fprintf('t: %0.3f / %0.3f (sec) \n', t(count-1), tmax);
end
end


% Newtonian acceleration update function
function a = acc(pos,mass,eps)
    
    a = zeros(size(pos));  % Set acceleration to zero
    N = size(pos,1);       % Number of particles
    
    for i = 1:N
        for j = 1:N
            if (i ~= j)    % No self interaction
                
                % Vector r_ij
                vect = pos(i, 1:3) - pos(j, 1:3);
                
                % Denominator with softening
                rd3 = (eps^2 + sum(vect.^2) )^(3/2);
                
                % Update acceleration
                a(i,1:3) = a(i,1:3) - mass(j) / rd3 * vect;
            end
        end
    end
    %}

end


% Potential energy function for each particle
function Ep = potential(pos, m)

N = size(pos,1);           % Number of particles
Ep = zeros(N,1);           % Potential energy vector

for i = 1:N
    for j = 1:N
        if (i ~= j)        % No self interaction
            Ep(i) = Ep(i) - m(i)*m(j) / norm(pos(i, 1:3) - pos(j, 1:3));
        end
    end
end

end

% Kinetic energy function for each particle
function Ek = kinetic(v, m)

N = size(v,1);
Ek = zeros(N,1);

for i = 1:N
    Ek(i) = 0.5*m(i).*norm(v(i,:)).^2;
end

end
