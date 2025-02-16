
% Thomas Key
% AA700 - Overall Simulation Concepts
%
%
% Notes:
% 
%
% This simulation uses QUATERNIONS. Matrices are used for visualization.
% Quaternions have the VECTOR as the LAST THREE entries, scalar FIRST.
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


%% Satellite Initial Conditions & Physical Parameters

% Define initial attitude and angular velocity
initXang = 0*(pi/4);
initYang = 0*(pi/4);
initZang = 0*(pi/4);

initXrot = 0*(pi/64);
initYrot = 10*(pi/64);
initZrot = 0*(pi/64); % rad/s

initRotVec = [initXrot;initYrot;initZrot];


% Define mass distribution
mIx = 1;
mIy = 1;
mIz = 10; % kg * m^2
mJmat = diag([mIx mIy mIz]);
inv_mJmat = inv(mJmat);


% Generate initial quaternion
initQuat = angle2quat(initZang,initYang,initXang)';
initQuat = initQuat/norm(initQuat); % Normalize


%% Orbital Initial Conditions & Physical Parameters

rEarth = 6378.137; % (meters) WGS84 Ellipsoid approximation
mEarth = 5.9722e24; % (kilogram) IAU Selected Astronomical Constants
gravConst = 6.67430e-11; % (m^2 kg^-1 s^-2) NIST Reference 2022 Codata
muEarth = (398600.4418); % (km^3 s^-2) Ries, J. C. 1992 "...Determination of the Gravitational Constant of the Earth" pp. 529-531

% ISS Orbital Parameters 2025, Feb. 14
orbitIncl = deg2rad(51.6371);   % rad
orbitEcct = 0.0003865;          % unitless
orbitRaan = deg2rad(193.6981);  % degrees
orbitPaps = 414 + rEarth;       % kilometers
orbitApog = 419 + rEarth;       % kilometers
orbitArgp = deg2rad(315.5499);  % rad

orbitSemi = (0.5)*(orbitPaps + orbitApog);              % kilometers
orbitAngm = sqrt(muEarth*orbitSemi*(1 - orbitEcct.^2)); % kg m^2 s^-1
orbitPeri = 2*pi*sqrt( (orbitSemi.^3)/(muEarth) );      % s


% Aerodynamic Drag to be added later

% Oblateness of the Earth

% Zonal Harmonics from Curtis, Orbital Mechanics for Engineering Students
% Better approximations Shaub and Junkins (2009) p. 553 (for future)
earthJ2 = 0.00108263;            % First Zonal Harmonic



%% Orbital Function Definitions

% Convert orbital parameters to RV (See Curtis p 191)
function out_RV = params2rv(muEarth, orbitParams)
    % orbitParams:
    % 
    % angm - Angular momentum
    % incl - Inclination
    % raan - Right Ascencion of the Ascending Node
    % ecct - Eccentricity
    % argp - Argument of Perigee
    % tano - True Anomaly

    % This should be converted to use quaternions in the future.

    % Position Perifocal
    rVal_p = ((orbitParams(1).^2)./muEarth).*(1./(1 + orbitParams(4).*cos(orbitParams(6))));
    rVec_p = rVal_p*[cos(orbitParams(6));sin(orbitParams(6));0];

    % Velocity Perifocal
    vVal_p = muEarth./orbitParams(1);
    vVec_p = vVal_p.*[-sin(orbitParams(6));orbitParams(4) + cos(orbitParams(6));0];

    % Annoying Rotation Matrices
    argp_c = cos(orbitParams(5));
    argp_s = sin(orbitParams(5));
    raan_c = cos(orbitParams(3));
    raan_s = sin(orbitParams(3));
    incl_c = cos(orbitParams(2));
    incl_s = sin(orbitParams(2));

    Q_z_argp = [ argp_c, argp_s, 0;
                -argp_s, argp_c, 0;
                 0,      0,      1];
    
    Q_z_raan = [ raan_c, raan_s, 0;
                -raan_s, raan_c, 0;
                 0,      0,      1];

    Q_x_incl = [ 1,      0,      0;
                 0,  incl_c, incl_s;
                 0, -incl_s, incl_c];

    Q_tot = (Q_z_argp*Q_x_incl*Q_z_raan)';

    
    rVec = Q_tot*rVec_p;
    vVec = Q_tot*vVec_p;


    out_RV = [rVec;vVec];

end

% Convert RV to orbital parameters
function out_orbitParams = rv2params(muEarth, rVec,vVec)
    % rVec & vVec are X, Y, Z in inertial frame

    % rVec - kilometers
    % vVec - kilometers/second

    rVal = sqrt(dot(rVec,rVec));

    vRad = dot(rVec,vVec)/rVal;

    angmVec = cross(rVec,vVec);
    angmVal = sqrt(dot(angmVec,angmVec));

    incl = acos(angmVec(3)/angmVal);

    nodeVec = cross([0;0;1],angmVec);
    nodeVal = sqrt(dot(nodeVec,nodeVec));

    if nodeVec(2) >= 0 

        raan = acos(nodeVec(1)/nodeVal);

    elseif nodeVec(2) < 0

        raan = 2*pi - acos(nodeVec(1)/nodeVal);

    end

    ecctVec = (1./muEarth).*(cross(vVec,angmVec) - rVec.*(muEarth./rVal));
    ecctVal = sqrt(dot(ecctVec,ecctVec));

    if ecctVec(3) >= 0

        argp = acos((dot(nodeVec,ecctVec))./(nodeVal.*ecctVal));

    elseif ecctVec(3) < 0

        argp = 2*pi - acos((dot(nodeVec,ecctVec))./(nodeVal.*ecctVal));

    end

    if (vRad >= 0)

        tano = acos(dot(ecctVec,rVec)./(ecctVal.*rVal));

    elseif (vRad < 0)

        tano = 2.*pi - acos(dot(ecctVec,rVec)./(ecctVal.*rVal));

    end

    out_orbitParams = [angmVal;incl;raan;ecctVal;argp;tano];

end


function out_orbitDeriv = orbitDeriv(~, rv_Input, muEarth)

    rVec = rv_Input(1:3);
    vVec = rv_Input(4:6);
    rVal = norm(rVec);
    oParam = (muEarth./(rVal.^3));
    
    % rhoDensity = 10e-13; % kg m^-3
    % C_d = 2.2;
    % ArVal = pi.*(1); % m^2
    % massVal = 100; % kg
    % v_Rel_vec = rv_Vecs(4:end) - cross([0;0;(72.9211e-6)],rv_Vecs(1:3)); % m/s
    % p_Aero = -(0.5).*(rhoDensity).*(norm(v_Rel_vec).^2).*(C_d).*(ArVal).*(1./massVal).*v_Rel_vec; % N

    distForce = [0;0;0];

    rdot = vVec;
    vdot = -oParam.*rVec + distForce;

    out_orbitDeriv = [rdot;vdot];


end
% Tests to confirm accuracy, based on examples in Curtis Chapter 4
%orbitParams_results = rv2params(muEarth,[-6045000;-3490000;2500000],[-3457;6618;2533]); % m^3 s^-2, meters, m/s
%rv_Results  = params2rv(muEarth,[80000*(1000^2);deg2rad(30);deg2rad(40);1.4;deg2rad(60);deg2rad(30)]);




%% Attitude Function Definitions

% For ode45, RK, and RKM.
function out_qOmega = qOmega(w)

    out_qOmega =  [    0,  w(3), -w(2), w(1); ...
                   -w(3),     0,  w(1), w(2); ...
                    w(2), -w(1),     0, w(3); ...
                   -w(1), -w(2), -w(3),   0];

end

% Quaternion multiplication
function qM = quatMult(q1,q2)

    % Quaternion derivative
    % See: "Quaternions and Dynamics" by Basile Graf
    % February 2007, found on Arxiv.org
    % Page 7, eq. 12 and Page 4, eq. 3

    if (norm(q1) == 0)
        
        % Prevent division by zero by forcing unity quaternion in omega
        q1(1) = 1;

    end

    if (norm(q2) == 0)
        
        % Prevent division by zero by forcing unity quaternion.
        q2(1) = 1;

    end
    
    % Don't normalize omega
    q2 = q2./norm(q2);

    pScal = q2(1);
    pVec = q2(2:4);
    qScal = q1(1);
    qVec = q1(2:4);

    outQ1 = pScal*qScal - dot(pVec,qVec);
    outQ2 = qScal*pVec + pScal*qVec + cross(qVec,pVec);
    qM = [outQ1;outQ2];

end

% For ode45
function outQ_der = qDeriv45(t,y,J,Ji)

    % Ensure normalized quaternions
    normQin = y(4:end)./norm(y(4:end));
    
    % Angular velocity vector -> Angular velocity quaternion
    wVec = y(1:3);
    wScal = 0;

    % Quaternion kinematics
    qM = 0.5*quatMult([wScal;wVec],normQin); % Derivative

    % Torque testing
    t_Test = -5/(t+2); % Torque slowly decreases
    t_Test2 = sin(t);

    % Body kinetics
    wTerm1 = Ji*([t_Test2;0;t_Test]);
    wTerm2 = -Ji*(cross((y(1:3)),J*y(1:3)));
    outW_der = wTerm1 + wTerm2;

    % Output
    outQ_der = [outW_der; qM];

end

% For RK and RKM
qDeriv = @(omeg, qin) (0.5)*(omeg)*(qin);
wDeriv = @(w_in) [0;0;0];

%% Simulation Pacing/Timing/Discretization

hours = 3600;

tStep     = 0.04;                          % Time-step (0.04~24fps)
timeTotal = 20;                            % Overall sim length
stepTot   = ceil(timeTotal/tStep);         % Calculate amount of steps
tVals     = linspace(0,timeTotal,stepTot); % Time values

% tVals_orbit = linspace(0, 4*hours, 1000);

%% Simulation parameters



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


%% ode45 (RKDP) for Attitude and Orbital Propagation 


[t_o45, y_o45] = ode45(@(t,y) qDeriv45(t, y, mJmat, inv_mJmat), tVals, [w_Val_o45(:, 1); q_Val_o45(:, 1)]);
y_o45 = y_o45';
w_Val_o45 = y_o45(1:3,:); % Outputs are reasonable
q_Val_o45 = y_o45(4:end,:);

initR = [8000; 0; 6000];
initV = [0; 7; 0];
[t_orbit,y_orbit] = ode45(@(t,y) orbitDeriv(t, y, muEarth), [0 6*hours], [initR;initV]);
y_orbit = y_orbit';

figure(90)
hold on
grid on
plot3(y_orbit(1,1),y_orbit(2,1),y_orbit(3,1),'ro')
plot3(y_orbit(1,:),y_orbit(2,:),y_orbit(3,:),'k')
axis equal
view([1,1,.4])
title('Testing Orbital Path')
xlabel('x')
ylabel('y')
zlabel('z')



%% Visualization

% Euler angle containers
z_q_RKM = zeros(1, stepTot);
y_q_RKM = zeros(1, stepTot);
x_q_RKM = zeros(1, stepTot);

z_q_RK  = zeros(1, stepTot);
y_q_RK  = zeros(1, stepTot);
x_q_RK  = zeros(1, stepTot);

z_q_o45 = zeros(1, stepTot);
y_q_o45 = zeros(1, stepTot);
x_q_o45 = zeros(1, stepTot);

% Convert from quaternions
for ii = 1:(stepTot)

    [z_q_RKM(ii), y_q_RKM(ii), x_q_RKM(ii)] = quat2angle(q_Val_RKM(:, ii)');
    [ z_q_RK(ii),  y_q_RK(ii),  x_q_RK(ii)] = quat2angle( q_Val_RK(:, ii)');
    [z_q_o45(ii), y_q_o45(ii), x_q_o45(ii)] = quat2angle(q_Val_o45(:, ii)');

end

% Fix discontinuities at 2pi increments
z_q_RKM = unwrap(z_q_RKM);
z_q_RK = unwrap(z_q_RK);
z_q_o45 = unwrap(z_q_o45);

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

run3D = true;

if (run3D == true)

    % NOTES:
    % 
    % The sim3d software uses its own coordinate system. 
    % I am using ISO8855, which is a right handed system.
    % ISO8855 and sim3D have y-axes that are inverted from one another.
    % The sim3d axis has left handed rotation about z.
    % The viewport will always use the sim3D origin and axes. 
    %
    % sim3d +x = ISO8855
    % sim3d -x = ISO8855
    % sim3d +y -> -y ISO8855
    % sim3d -y -> +y ISO8855
    % sim3d +z = ISO8855
    % sim3d -z = ISO8855 (positive rotation)

    % Generate world with update function and output function
    world = sim3d.World('Output',@outputImpl,'Update',@updateImpl);
    
    % Choose view of the world.
    % Viewport starts automatically at origin, facing ISO5588 and Default +x.
    viewport = createViewport(world, Translation=[2 -2 0.5],Rotation=[0 0 -(5/4)*pi]);

    pointingAngle = 3*(pi/4);

    % Axis label texts
    texts = sim3d.graphics.Text(ActorName='Texts', ...
        Translation=[1.25 0 0; 0 -1.25 0; 0 0 1.25], ...
        Rotation = [0 0 pointingAngle; 0 0 pointingAngle; 0 0 pointingAngle]);
    texts.Color = [1 0 0; 0 1 0; 0 0 1];
    texts.String = ["X";"Y";"Z"];
    add(world,texts);

    % Timestamp text
    textTime = sim3d.graphics.Text(ActorName='TimeText', Translation = [0.75 -0.75 -0.25], Rotation = [0 0 pointingAngle]);
    textTime.Color = [0 0 0];
    formatSpec = '%.3f';
    textTime.String = strcat("Time: ", num2str(world.SimulationTime,formatSpec), " s");
    add(world,textTime);
    
    % Define a coordinate syste,
    coordSys = 'ISO8855';
    
    % User-generated function for making axes.
    createAxes(world,coordSys);
    
    % Include actors
    box = sim3d.Actor(ActorName='Box');               % Generate actor
    createShape(box,'box',[0.25 0.5 0.75]);           % Add box mesh to actor
    box.Color = [0.1 0.1 0.1];                        % Set color of actor
    box.Mobility = sim3d.utils.MobilityTypes.Movable; % Make actor mobile
    box.CoordinateSystem = coordSys;                  % Set actor coord. system
    add(world,box);                                   % Include actor in world
    
    % User defined step data
    world.UserData.Step = 1;
    rot_Vals_o45 = zeros(size(q_Val_o45,2),3);
    for ii = 1:size(q_Val_o45,2)
        
        % Consider vertcat
        % sim3d takes values in degrees
        rot_Vals_o45(ii,:) =  [rad2deg(x_q_o45(ii)), rad2deg(y_q_o45(ii)), rad2deg(z_q_o45(ii))];
    
    end
    
    world.UserData.RotationValues = rot_Vals_o45;
    
    % Set simulation parameters
    sampletime = tStep;
    no_of_steps = stepTot;
    stoptime = no_of_steps*sampletime;
    
    % % Not sure how to get pacing to actually work. Bug?
    world.EnablePacing = false;
    world.PacingRate = 0.5;
    
    % Run world
    run(world,sampletime,stoptime)
    
    % Delete world upon end of sim
    % delete(world);

end

% This is a wild way of running this code, but it works.
    function outputImpl(world)
    
        box = world.Actors.Box;
        currentStep = world.UserData.Step;
        display(strcat('Current step:', num2str(currentStep)));
        zRot = world.UserData.RotationValues(currentStep, 3);
        yRot = world.UserData.RotationValues(currentStep, 2);
        xRot = world.UserData.RotationValues(currentStep, 1);
        box.Rotation = [xRot, yRot, zRot];
        cTime = num2str(world.SimulationTime,'%.2f');
        world.Actors.TimeText.String = strcat("Time: ", cTime , " s");
        
    end
    
    function updateImpl(world)
    
        % Updates the step in the function
        if (world.UserData.Step == size(world.UserData.RotationValues,1))
            
            % Prevent indexing issues with final step
            world.UserData.Step = world.UserData.Step; 
    
        else
            
            % Step forward normally
            world.UserData.Step = world.UserData.Step + 1; 
    
        end
        
        % Cheap way to pace simulation
        % pause(world.SampleTime)
    
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
        createShape(xA,'box',[1 0.02 0.02]);
        createShape(yA,'box',[0.02 1 0.02]);
        createShape(zA,'box',[0.02 0.02 1]);
    
        % Shift centerpoints to have tips at origin
        xA.Translation = [0.5 0.01 0.01];
        yA.Translation = [0.01 0.5 0.01];
        zA.Translation = [0.01 0.01 0.5];
    
        xA.Color = [1 0 0]; % Red
        yA.Color = [0 1 0]; % Green
        zA.Color = [0 0 1]; % Blue
    
        % Include actors
        add(world,xA);
        add(world,yA);
        add(world,zA);
    
    end