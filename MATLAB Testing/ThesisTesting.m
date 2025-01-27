
% Thomas Key
% AA700 - Overall Simulation Concepts
%
%
% Notes:
% 
%
% This simulation uses QUATERNIONS. Matrices are used for visualization.
% Quaternions have the VECTOR as the FIRST THREE entries, scalar after.
% ALL UNITS are SI. Conversions to imperial are SEPARATE VARIABLES.
% ALL ANGLES are going to be in RADIANS. Always specify if not radians.
% Orbital distances are in KM unless otherwise specified. Always specify.
% Control distances are in M unless otherwise specified. Always specify.
%
%
% Frames of Reference
%
% B - Body Frame
% M - Center of Mass
% I - Inertial Frame

close all
clear all
clc



% % Global Orbital Parameters
% 
% gG = 9.81; % Gravity for Earth
% rE = 6371; % Earth radius in KM
% 
% 
% 
% % Situtional Orbital Parameters
% 
% initAlt = 100; 
% initPos = [rE + initAlt,0,0];
% initVel = [0,5,0]; % A particularly fast circular orbit, for now. KM/s



%% Initial Conditions & Physical Parameters

% Initial angles
initXang = 0;
initYang = 0;
initZang = 0;
initAngVec = [initXang;initYang;initZang];

% Initial angular velocity
initXrot = 0;
initYrot = 0;
initZrot = 0.1; % rad/s
initRotVec = [initXrot;initYrot;initZrot];

% Moments of inertia
mIx = 1;
mIy = 1;
mIz = 1; % kg * m^2
mImat = diag([mIx mIy mIz]);
inv_mImat = inv(mImat);

% Attitude
initQuat = angle2quat(initZang,initYang,initXang)';
% Adjust from MATLAB's [SVVV] to [VVVS] format
initQuat = [initQuat(2), initQuat(3), initQuat(4), initQuat(1)];
initQuat = initQuat/norm(initQuat);


%% Function Definitions

function out_qOmega = qOmega(w)

    out_qOmega =  [    0,  w(3), -w(2), w(1); ...
                   -w(3),     0,  w(1), w(2); ...
                    w(2), -w(1),     0, w(3); ...
                   -w(1), -w(2), -w(3),   0];

end

function outQ_der = qDeriv45(~,y)

    outQ_der = [[0;0;0];(0.5)*(qOmega(y(1:3)))*(y(4:end))];

end

qDeriv = @(omeg, qin) (0.5)*(omeg)*(qin);
wDeriv = @(w_in) [0;0;0];


%% Simulation parameters

% Step considerations
tStep     = 0.04;                          % Time-step (0.04~24fps)
timeTotal = 60;                            % Overall sim length
stepTot   = ceil(timeTotal/tStep);         % Calculate amount of steps
tVals     = linspace(0,timeTotal,stepTot); % Time values

% Containers
er_RKM_q  = zeros(4, stepTot); % Quaternion error from RKM 
er_RKM_w  = zeros(4, stepTot); % Angular vel. error from RKM
q_Val_RK  = zeros(4, stepTot); % Quaternion for RK
q_Val_RKM = zeros(4, stepTot); % Quaternion for RKM
q_Val_o45 = zeros(4, stepTot); % Quaternion for ode45
q_Dot     = zeros(4, stepTot); % Quaternion deriv.
yaw_Val   = zeros(1, stepTot); % Yaw 
pitch_Val = zeros(1, stepTot); % Pitch 
roll_Val  = zeros(1, stepTot); % Roll
w_Val_o45 = zeros(3, stepTot); % Omega
w_Val_RKM = zeros(3, stepTot); % Omega
w_Val_RK  = zeros(3, stepTot); % Omega
w_Dot     = zeros(3, stepTot); % Angular vel. vec.
tq_Vec    = zeros(3, stepTot); % Input


% Initializations
q_Val_RK(:,1)  = initQuat;   % Quaternion init.
q_Val_RKM(:,1) = initQuat;   % Quaternion init.
q_Val_o45(:,1) = initQuat;   % Quaternion init.
w_Val_o45(:,1) = initRotVec; % Omega init.
w_Val_RKM(:,1) = initRotVec; % Omega init.
w_Val_RK(:,1)  = initRotVec; % Omega init.
tq_Vec(:,1)    = [0,10,1];   % Input init.





%% Stepping Methods

% These need to be updated to follow RDKP, or the Runga Kutta Dormand
% Prince method. This is faster, and is what is technically running behind
% the scenes with ode45. While not necessary in MATLAB, this will need to
% be written explicity in Python/C++/C later on. This is for my own
% learning and understanding of a more practical implementation. 



% RKM Stepping Constants
alphas_RKM = [(1/3), (1/3), (1/2),     1];
betas_RKM  = [(1/3), (1/6), (1/6), (1/8),     0, (3/8), (1/2), 0, (-3/2), 2];
gammas_RKM = [(1/6),     0,     0, (2/3), (1/6)];

% RK4 Stepping Constants
alphas_RK = [(1/2), (1/2),     1];
betas_RK  = [(1/2),     0, (1/2),     0, 0, 1];
gammas_RK = [(1/6), (1/3), (1/3), (1/6)];

% Runga-Kutta-Merson (RKM) Stepping
for ii = 1:(stepTot)

    % Quaternion Derivatives from M. Sidi, 1997, page 104

    k1   = tStep*( qDeriv( qOmega( w_Val_RKM(:, ii) ), q_Val_RKM(:, ii) ) );
    k2   = tStep*( qDeriv( qOmega( w_Val_RKM(:, ii) ), q_Val_RKM(:, ii) ) + k1*betas_RKM(1) );
    k3   = tStep*( qDeriv( qOmega( w_Val_RKM(:, ii) ), q_Val_RKM(:, ii) ) + k1*betas_RKM(2) + k2*betas_RKM(3) );
    k4   = tStep*( qDeriv( qOmega( w_Val_RKM(:, ii) ), q_Val_RKM(:, ii) ) + k1*betas_RKM(4) + k2*betas_RKM(5) + k3*betas_RKM(6) );
    k5   = tStep*( qDeriv( qOmega( w_Val_RKM(:, ii) ), q_Val_RKM(:, ii) ) + k1*betas_RKM(7) + k2*betas_RKM(8) + k3*betas_RKM(9) + k4*betas_RKM(10) );

    k1_w = tStep*( wDeriv( w_Val_RKM(:, ii) ) );
    k2_w = tStep*( wDeriv( w_Val_RKM(:, ii) ) + k1_w*betas_RKM(1) );
    k3_w = tStep*( wDeriv( w_Val_RKM(:, ii) ) + k1_w*betas_RKM(2) + k2_w*betas_RKM(3) );
    k4_w = tStep*( wDeriv( w_Val_RKM(:, ii) ) + k1_w*betas_RKM(4) + k2_w*betas_RKM(5) + k3_w*betas_RKM(6) );
    k5_w = tStep*( wDeriv( w_Val_RKM(:, ii) ) + k1_w*betas_RKM(7) + k2_w*betas_RKM(8) + k3_w*betas_RKM(9) + k4_w*betas_RKM(10) );
    

    q_Val_RKM(:, ii + 1) = q_Val_RKM(:, ii) + k1*gammas_RKM(1)   + k2*gammas_RKM(2)   + k3*gammas_RKM(3)   + k4*gammas_RKM(4)   + k5*gammas_RKM(5);
    w_Val_RKM(:, ii + 1) = w_Val_RKM(:, ii) + k1_w*gammas_RKM(1) + k2_w*gammas_RKM(2) + k3_w*gammas_RKM(3) + k4_w*gammas_RKM(4) + k5_w*gammas_RKM(5);

    er_RKM_q(:, ii) = (1/30)*( 2*k1 - 9*k3 + 8*k4 - k5 );
 
end

% Regular RK Stepping
for ii = 1:(stepTot)

    % Quaternion Derivatives from M. Sidi, 1997, page 104

    k1   = tStep*( qDeriv( qOmega( w_Val_RK(:, ii) ), q_Val_RK(:, ii) ) );
    k2   = tStep*( qDeriv( qOmega( w_Val_RK(:, ii) ), q_Val_RK(:, ii) ) + k1*betas_RK(1) );
    k3   = tStep*( qDeriv( qOmega( w_Val_RK(:, ii) ), q_Val_RK(:, ii) ) + k1*betas_RK(2) + k2*betas_RK(3) );
    k4   = tStep*( qDeriv( qOmega( w_Val_RK(:, ii) ), q_Val_RK(:, ii) ) + k1*betas_RK(4) + k2*betas_RK(5) + k3*betas_RK(6));
    
    k1_w = tStep*( wDeriv( w_Val_RK(:, ii) ) );
    k2_w = tStep*( wDeriv( w_Val_RK(:, ii) ) + k1_w*betas_RK(1) );
    k3_w = tStep*( wDeriv( w_Val_RK(:, ii) ) + k1_w*betas_RK(2) + k2_w*betas_RK(3));
    k4_w = tStep*( wDeriv( w_Val_RK(:, ii) ) + k1_w*betas_RK(4) + k2_w*betas_RK(5) + k3_w*betas_RK(6));

    q_Val_RK(:, ii + 1) = q_Val_RK(:, ii) + k1*gammas_RK(1)   + k2*gammas_RK(2)   + k3*gammas_RK(3)   + k4*gammas_RK(4);
    w_Val_RK(:, ii + 1) = w_Val_RK(:, ii) + k1_w*gammas_RK(1) + k2_w*gammas_RK(2) + k3_w*gammas_RK(3) + k4_w*gammas_RK(4);
 
end

% RKDP or ode45
[t_o45, y_o45] = ode45(@qDeriv45, tVals, [w_Val_o45(:, 1);q_Val_o45(:, 1)]);
y_o45 = y_o45';
w_Val_o45 = y_o45(1:3,:);
q_Val_o45 = y_o45(4:end,:);


%% Visualization

% Euler angle containers
z_q_RKM = zeros(1, stepTot);
y_q_RKM = zeros(1, stepTot);
x_q_RKM = zeros(1, stepTot);

z_q_RK = zeros(1, stepTot);
y_q_RK = zeros(1, stepTot);
x_q_RK = zeros(1, stepTot);

z_q_o45 = zeros(1, stepTot);
y_q_o45 = zeros(1, stepTot);
x_q_o45 = zeros(1, stepTot);

% Convert from quaternions
for ii = 1:(stepTot)

    [x_q_RKM(ii), y_q_RKM(ii), z_q_RKM(ii)] = quat2angle(q_Val_RKM(:, ii)');
    [ x_q_RK(ii),  y_q_RK(ii),  z_q_RK(ii)] = quat2angle( q_Val_RK(:, ii)');
    [x_q_o45(ii), y_q_o45(ii), z_q_o45(ii)] = quat2angle(q_Val_o45(:, ii)');

end

% Show results

showPlots = false;

if (showPlots == true)

    figure(8)
    hold on
    grid on
    plot(tVals, rad2deg(x_q_RKM))
    plot(tVals, rad2deg(y_q_RKM))
    plot(tVals, rad2deg(z_q_RKM))
    title('Rotation angles RKM')
    legend('x rot.','y rot.','z rot.')
    xlabel('time')
    ylabel('Angles in rad')
    
    figure(9)
    hold on
    grid on
    plot(tVals, rad2deg(x_q_RK))
    plot(tVals, rad2deg(y_q_RK))
    plot(tVals, rad2deg(z_q_RK))
    title('Rotation angles RK')
    legend('x rot.','y rot.','z rot.')
    xlabel('time')
    ylabel('Angles in rad')
    
    figure(10)
    hold on
    grid on
    plot(tVals, rad2deg(x_q_o45))
    plot(tVals, rad2deg(y_q_o45))
    plot(tVals, rad2deg(z_q_o45))
    title('Rotation angles ode45')
    legend('x rot.','y rot.','z rot.')
    xlabel('time')
    ylabel('Angles in rad')
    
    figure(11)
    hold on
    grid on
    plot(tVals, er_RKM_q)
    legend('q1','q2','q3','qs')
    title('RKM Error')
    xlabel('Time')
    ylabel('Error')

end

%% 3D Representation

% Generate world with update function and output function
world = sim3d.World('Output',@outputImpl,'Update',@updateImpl);

% Choose view of the world
viewport = createViewport(world, Translation=[-4 2 1],Rotation=[0 0 -pi/8]);

% Define a coordinate syste,
coordSys = 'MATLAB';

% User-generated function for making axes.
createAxes(world,coordSys);

% Include actors
box = sim3d.Actor(ActorName='Box');               % Generate actor
createShape(box,'box',[0.25 0.25 0.25]);          % Add box mesh to actor
box.Color = [0.1 0.1 0.1];                        % Set color of actor
box.Mobility = sim3d.utils.MobilityTypes.Movable; % Make actor mobile
box.CoordinateSystem = coordSys;                  % Set actor coord. system
add(world,box);                                   % Include actor in world

% User defined step data
world.UserData.Step = 1;

% Set simulation parameters
sampletime = 0.01;
stoptime = 30;

% Run world
run(world,sampletime,stoptime)

% Delete world upon end of sim
delete(world);

% This is a wild way of running this code, but it works.
function outputImpl(world)
% Sets the actor outputs (e.g. an actor position to follow a path)
% Create actors in scene
    timePeriod = uint32(floor(world.UserData.Step/500))  % From 0->6000
    movementTime = world.UserData.Step - (timePeriod*500); % From 0-500 (12 times)
    deltaRotation = 60*(pi/180)/250;
    box = world.Actors.Box;
    switch timePeriod
        case 0
            box.Color = [1 0 0];
            if movementTime < 250
                box.Translation(1) = box.Translation(1) + 0.01;
            else
                box.Translation(1) = box.Translation(1) - 0.01;
            end
        case 1
            box.Color = [0 1 0];
            if movementTime < 250
                box.Translation(2) = box.Translation(2) + 0.01;
            else
                box.Translation(2) = box.Translation(2) - 0.01;
            end
        case 2
            box.Color = [0 0 1];
            if movementTime < 250
                box.Translation(3) = box.Translation(3) + 0.01;
            else
                box.Translation(3) = box.Translation(3) - 0.01;
            end
        case 3
            box.Color = [1 0 0];
            if movementTime < 250
                box.Rotation(1) = box.Rotation(1) + deltaRotation;
            else
                box.Rotation(1) = box.Rotation(1) - deltaRotation;
            end
        case 4 
            box.Color = [0 1 0];
            if movementTime < 250
                box.Rotation(2) = box.Rotation(2) + deltaRotation;
            else
                box.Rotation(2) = box.Rotation(2) - deltaRotation;
            end
        case 5
            box.Color = [0 0 1];
            if movementTime < 250
                box.Rotation(3) = box.Rotation(3) + deltaRotation;
            else
                box.Rotation(3) = box.Rotation(3) - deltaRotation;
            end
    end
end

function updateImpl(world)

    % Updates the step in the function
    world.UserData.Step = world.UserData.Step + 1;

end

function createAxes(world, csys)

    % Generates visible axes
    % Axes X, Y, Z represented by colors R, G, B respectively

    % Define actors
    xA = sim3d.Actor(ActorName='XAxis');
    yA = sim3d.Actor(ActorName='YAxis');
    zA = sim3d.Actor(ActorName='ZAxis');

    % Define coordinate systems
    xA.CoordinateSystem = csys;
    yA.CoordinateSystem = csys;
    zA.CoordinateSystem = csys;

    % Make thin boxes for the arrows
    createShape(xA,'box',[3 0.02 0.02]);
    createShape(yA,'box',[0.02 3 0.02]);
    createShape(zA,'box',[0.02 0.02 3]);

    % Shift centerpoints to have tips at origin
    xA.Translation = [1.5 0.01 0.01];
    yA.Translation = [0.01 1.5 0.01];
    zA.Translation = [0.01 0.01 1.5];

    xA.Color = [1 0 0]; % Red
    yA.Color = [0 1 0]; % Green
    zA.Color = [0 0 1]; % Blue

    % Include actors
    add(world,xA);
    add(world,yA);
    add(world,zA);
end