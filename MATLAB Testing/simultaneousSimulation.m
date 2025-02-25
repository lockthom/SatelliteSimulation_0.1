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
%
% Frames of Reference
%
% B - Body Frame
% M - Center of Mass
% I - Inertial Frame

%% Satellite Initial Conditions & Physical Parameters

% Define initial attitude and angular velocity
initXang = 0*(pi/4); % radians
initYang = 0*(pi/4);
initZang = 0*(pi/4);

initXrot = 0*(pi/64); % radians/s
initYrot = 0*(pi/64);
initZrot = 0*(pi/64); 

initRotVec = [initXrot;
              initYrot;
              initZrot];

% Define mass distribution
mIx = 1;
mIy = 1;
mIz = 10; % kg * m^2
mJmat = diag([mIx mIy mIz]);
inv_mJmat = inv(mJmat);
satCD = 2.2;
satA  = 1; % Should adjust for quaternion, but isn't super important now.

satParams = [satCD; satA];


% Generate initial quaternion
initQuat = angle2quat(initZang,initYang,initXang)';
initQuat = initQuat/norm(initQuat); % Normalize

%% Orbital Initial Conditions & Physical Parameters

rBody     = 6371.230;      % Frédéric Chambat; Bernard Valette (2001) "Physics of the Earth and Planetary Interiors"
rBody_eq  = 6378.1370;     % (kilometers) WGS84 Ellipsoid approximation
mBody     = 5.9722e24;     % (kilogram) IAU Selected Astronomical Constants
gravConst = (6.67430e-20); % (km^3 kg^-1 s^-2) NIST Reference 2022 Codata
muBody    = (398600.4418); % (km^3 s^-2) Ries, J. C. 1992 "...Determination of the Gravitational Constant of the Earth" pp. 529-531
rotBody   = 7.2921150e-5;  % (rad/s) IERS Numerical Standards 1999, mean ang. vel.

% Initial orbital parameters: (Good comparison example with Curtis book)
alt_p0    = 300;    % Perigee altitude (km)
alp_a0    = 3062;   % Apogee altitude (km)
orbitRaan = 45*deg; % Right ascension of the node (radians)
orbitIncl = 28*deg; % Inclination (radians)
orbitArgp = 30*deg; % Argument of perigee (radians)
orbitTano = 40*deg; % True anomaly (radians)
 
% Calculated values
rp0        = rBody + alt_p0;           % Perigee radius (km)
ra0        = rBody + alp_a0;           % Apogee radius (km)
orbitEcct  = (ra0 - rp0)/(ra0 + rp0);  % Eccentricity
a0         = (ra0 + rp0)/2;            % Semimajor axis (km)
orbitAngm  = sqrt(rp0*muBody*(1 + orbitEcct));  % Angular momentum (km^2/s)
orbitPeri  = 2*pi/sqrt(muBody)*a0^1.5; % Period (s)


% Zonal Harmonics from Curtis, Orbital Mechanics for Engineering Students
% Better approximations Shaub and Junkins (2009) p. 553 (for future)
bodyJ2 = 0.00108263;            % First Zonal Harmonic

bodyParams = [muBody;
               bodyJ2;
               rBody;
               rBody_eq;
               rotBody];

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

%% Orbital Function Definitions

% Convert orbital parameters to RV (See Curtis p 191)
function out_RV = params2rv(bodyParams, orbitParams)
    % orbitParams:
    % 
    % angm - Angular momentum
    % incl - Inclination
    % raan - Right Ascencion of the Ascending Node
    % ecct - Eccentricity
    % argp - Argument of Perigee
    % tano - True Anomaly

    % This should be converted to use quaternions in the future.

    muBody = bodyParams(1);

    % Position Perifocal
    rVal_p = ((orbitParams(1).^2)./muBody).*(1./(1 + orbitParams(4).*cos(orbitParams(6))));
    rVec_p = rVal_p*[cos(orbitParams(6));sin(orbitParams(6));0];

    % Velocity Perifocal
    vVal_p = muBody./orbitParams(1);
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
function out_orbitParams = rv2params(bodyParams, rVec,vVec)
    % rVec & vVec are X, Y, Z in inertial frame

    % rVec - kilometers
    % vVec - kilometers/second

    muBody = bodyParams(1);

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

    ecctVec = (1./muBody).*(cross(vVec,angmVec) - rVec.*(muBody./rVal));
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

% Get the Keplerian orbit position and velocity after time del_t
function out_RVkeplerian = orbitUniVar(bodyParams, rvVec, del_t, configVals)
    % Inputs: 
    % 1: bodyParams -> Planetary Body Information
    % 2: rvVec      -> Starting position and velocity
    % 3: del_t      -> Length of timestep
    % 4: configVals -> Tolerance for universal anomaly root-finding
    %
    % Outputs:
    % 1: out_RVkeplerian -> Position and velocity after del_t

    % Based off of Curtis. Neat concept, and useful for simulation.
    % uano -> Universal Anomaly
    % cStumpff -> "C" Stumpff Function (Karl Stumpff, 1895-1970)
    % sStumpff -> "S" Stumpff Function
    %
    % Needs functionality for failed condition with tolerances and step
    % sizes when calculating the universal anomaly. 

    % Pull out relevant info for stepping
    stepMax  = configVals(1);
    tolVal   = configVals(2);
    muBody  = bodyParams(1);

    % Format inputs for use
    rVal = norm(rvVec(1:3));
    vVal = norm(rvVec(4:6));
    vRad = dot(rvVec(1:3),rvVec(4:6))./rVal;

    % alphaVal < 0 --> Hyperbola
    % alphaVal = 0 --> Parabola
    % alphaVal > 0 --> Ellipse
    alphaVal = 2./rVal - ((vVal.^2)./muBody);

    % Good estimate for initial uano (Chobotov, 2002) citation from Curtis
    uano = sqrt(muBody).*abs(alphaVal).*del_t;
    uano_ratio = 1;
    stepVal = 0;
    while abs(uano_ratio) > tolVal && stepVal <= stepMax

        uniParam = (uano^2)*alphaVal;

        stepVal = stepVal + 1;

        if (uniParam > 0)
    
            cStumpff = ( 1./uniParam ).*( 1 - cos(sqrt(uniParam)) );
            sStumpff = ( uniParam.^(-3./2) ).*( sqrt(uniParam) - sin(sqrt(uniParam)) );
    
        elseif (uniParam < 0)
    
            cStumpff = ( -1./uniParam ).*( cosh(sqrt(-uniParam)) - 1 );
            sStumpff = ( (-uniParam).^(-3./2) ).*( sinh(sqrt(-uniParam)) - sqrt(-uniParam) );
    
        else
    
            cStumpff = 1./2;
            sStumpff = 1./6;
    
        end

        % See Curtis P. 171

        uano_c0 = (rVal).*(vRad).*(1./sqrt(muBody));
        uano_c1 = (1 - alphaVal.*rVal);

        f_uano   = uano_c0.*(uano.^2).*cStumpff + uano_c1.*(uano.^3).*sStumpff + rVal.*uano - sqrt(muBody).*del_t;
        df_duano = uano_c0.*uano.*(1 - alphaVal.*(uano.^2).*sStumpff) + uano_c1.*(uano.^2).*cStumpff + rVal;

        uano_ratio = f_uano./df_duano;

        uano = uano - uano_ratio;

    end

    % At the end of this while loop the Universal Anomaly is available in
    % the uano variable. The stumpff functions don't need to have any
    % memory so the most recent value is already stored. The universal
    % anomaly is placed into Lagrange coefficients for the calculation of
    % the osculating orbit.

    f_L  = 1 - ((uano.^2)./rVal).*cStumpff;
    g_L  = del_t - (1./sqrt(muBody)).*(uano.^3).*sStumpff;

    rNew = f_L.*rvVec(1:3) + g_L.*rvVec(4:6);

    df_L = (sqrt(muBody)./(rVal.*norm(rNew))).*(alphaVal.*(uano.^3).*sStumpff - uano);
    dg_L = 1 - ((uano.^2)./norm(rNew)).*cStumpff;

    vNew = df_L.*rvVec(1:3) + dg_L.*rvVec(4:6);

    % This output is the osculating orbit after some time del_t. It does
    % not account for any perturbations.
    out_RVkeplerian = [rNew; vNew];


end

% Determine the force vector on the orbiter due to oblateness of the body
function out_oblateForce = oblateForce(bodyParams, rvVec)
    % Inputs:
    % 1: bodyParams -> Planetary body information
    % 2: rVec_osc   -> Position after time delta_t in keplerian orbit
    %
    % Outputs:
    % 1: Force due to planet oblateness in 3x1 inertial vector.

    % Later on, it may be interesting to look at Shaub and Junkins (2009)
    % on page 553 for accelerations due to J3 -> J6. 

    muBody = bodyParams(1);
    bodyJ2 = bodyParams(2);
    rBody_equatorial = bodyParams(3);
    
    rVal = norm(rvVec);

    force_constant1 = (3./2).*(bodyJ2).*(muBody).*(rBody_equatorial^2).*(rVal^(-4));
    force_constant2 = 5.*((rvVec(3)./rVal).^2);

    force_x = force_constant1.*(rvVec(1)./rVal).*(force_constant2 - 1);
    force_y = force_constant1.*(rvVec(2)./rVal).*(force_constant2 - 1);
    force_z = force_constant1.*(rvVec(3)./rVal).*(force_constant2 - 3);

    out_oblateForce = [force_x;
                       force_y;
                       force_z];

end

function out_aeroDrag = aeroDrag(bodyParams, satParams, rvVec)

    altVal = norm(rvVec(1:3));
    C_D = satParams(1); 
    A = satParams(2);

    h = [0, 25, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 180, 200, 250, 300, 350, 400, 450, 500, 600, 700, 800, 900, 1000];
     
    %...Corresponding densities (kg/m^3) from USSA76:  
    r = [1.225, 4.008e-2, 1.841e-2, 3.996e-3, 1.027e-3, 3.097e-4, 8.283e-5, 1.846e-5, 3.416e-6, 5.606e-7, 9.708e-8, 2.222e-8, 8.152e-9, 3.831e-9, 2.076e-9, 5.194e-10, 2.541e-10, 6.073e-11, 1.916e-11, 7.014e-12, 2.803e-12, 1.184e-12, 5.215e-13, 1.137e-13, 3.070e-14, 1.136e-14, 5.759e-15, 3.561e-15];     
      
    %...Scale heights (km):
    H = [7.310, 6.427, 6.546, 7.360, 8.342, 7.583, 6.661, 5.927, 5.533, 5.703, 6.782, 9.973, 13.243, 16.322, 21.652, 27.974, 34.934, 43.342, 49.755, 54.513, 58.019, 60.980, 65.654, 76.377, 100.587, 147.203, 208.020]; 
     
    %...Handle altitudes outside of the range:
    if altVal > 1000
        altVal = 1000;
    elseif altVal < 0
        altVal = 0;
    end
     
    %...Determine the interpolation interval:
    for j = 1:27
        if altVal >= h(j) && altVal < h(j+1)
            i = j;
        end
    end
    if altVal == 1000
        i = 27;
    end
     
    %...Exponential interpolation:
    density = r(i)*exp(-(altVal - h(i))/H(i));

    vRel = rvVec(4:6) - cross([0;0;bodyParams(5)],rvVec(1:3));

    out_aeroDrag = (0.5)*(density)*(norm(vRel)^2)*C_D*A*vRel;

end

% Sum forces and generate deviations
function out_orbitDisp = orbitDisp(t, pVec, rv_Current, tp, bodyParams, satParams, configVals)
    % Inputs:
    % 1: t           -> Time in seconds
    % 2: rv_Input    -> Position and Velocity in 6x1 vector (km)
    % 3: bodyParams  -> Planetary body information
    %
    % Outputs:
    % 1: out_orbitDeriv -> del_a and del_v in 6x1 vector (km)

    oscVals = orbitUniVar(bodyParams, rv_Current, t - tp, configVals);

    rv_Pert = oscVals + pVec;

    % Aerodynamic Drag
    force_aeroDrag   = aeroDrag(bodyParams, satParams, rv_Pert);

    % Oblateness Effects
    force_oblateness = oblateForce(bodyParams, rv_Pert);

    % Solar Pressure
    % Not yet included

    % Combine all active forces
    distForce = force_aeroDrag; % + others;

    % Determine del_a and del_v
    r_K = norm(rv_Current(1:3));
    r_P = r_K + norm(pVec(1:3));

    F     = 1 - (r_K/r_P)^3; % Adjust from the book, when you have time
    del_a = -(bodyParams(1)/(r_K^3))*(pVec(1:3) - F*rv_Pert(1:3)) + distForce;
     
    out_orbitDisp  = [pVec(4:6); del_a]; 

end


%% Simulation Pacing/Timing/Discretization for Attitude

hours = 3600;
days  = 24*hours;

tStep     = 0.04;                          % Time-step (0.04~24fps)
timeTotal = 20;                            % Overall sim length
stepTot   = ceil(timeTotal/tStep);         % Calculate amount of steps
tVals     = linspace(0,timeTotal,stepTot); % Time values


%% Orbital Simulation: Encke Integration + ode45

% Test for orbitUnivar from Curtis Example 3.7
% orbitUniVar(bodyParams, rvVec, del_t, configVals)
%[testOrbitUniVar] = orbitUniVar(bodyParams, [7000;-12124;0;2.6679;4.6210;0],60*60,[1000;1.e-8]);

rv_Inits = params2rv(bodyParams, orbitParams);
 
t0      = 0; 
tF      = 100*days;                 % Initial and final time (s)
tS      = t0;
del_t   = orbitPeri/10;            % Time step for Encke procedure
options = odeset('maxstep', del_t); % Tolerance


% Determine the osculating orbit at the end of the timestep
maxSteps    = 1000;
tolVal      = 1.e-8;
uano_Params = [maxSteps; tolVal];

% Initial Perturbation
pVec = [0;0;0;0;0;0];

% Container for rvOutputs and times
actualOrbit = rv_Inits;

% First timestep
tp = t0;
tC = t0 + del_t;

% Consider Sperling Burdet method
while tC <= tF + del_t/2
     
    % orbitDisp = rates
    [~,disp_derivs] = ode45(@(t,y) orbitDisp(t, y, rv_Inits, tp, bodyParams, satParams, uano_Params), [tp tC], pVec);

    % At the starting point, propagate the orbit as if there were no
    % perturbations at all, to determine where the orbiter would have been.
    [oscVals] = orbitUniVar(bodyParams, rv_Inits, del_t, uano_Params);

    % Add the changes in position and velocity to the ideal location, using
    % data from the previous step (or initial conditions)
    rvP = oscVals + disp_derivs(end,:)';

    % Save time
    tS = [tS,tC];

    % Start from new location
    rv_Inits = rvP;

    % Store data
    actualOrbit = [actualOrbit,rvP];

    % New timestep
    tp = tC;
    tC = tC + del_t;

    pVec = zeros(6,1); % Needed?


end

% [X_S, Y_S, Z_S] = sphere;
% 
% figure(5)
% hold on
% grid on
% plot3(actualOrbit(1,:), actualOrbit(2,:), actualOrbit(3,:), 'k')
% surf(X_S * rBody, Y_S * rBody, Z_S * rBody)
% axis equal
% title("Orbit Propagation")
% xlabel("X")
% ylabel("Y")
% zlabel("Z")
% view(3)

altVals = zeros(1,length(actualOrbit));
for ii = 1:length(actualOrbit)
    
    altVals(ii) = norm(actualOrbit(1:3)) - bodyParams(3);

end

figure(5)
hold on
grid on
plot(tS/days, altVals, 'k')
title('Altitude over 100 Days')
xlabel('Time (days)')
ylabel('Altitude (km)')
 

%% Attitude Simulation parameters



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

%% ode45 (RKDP) for Attitude Propagation 


[t_o45, y_o45] = ode45(@(t,y) qDeriv45(t, y, mJmat, inv_mJmat), tVals, [w_Val_o45(:, 1); q_Val_o45(:, 1)]);
y_o45 = y_o45';
w_Val_o45 = y_o45(1:3,:); % Outputs are reasonable
q_Val_o45 = y_o45(4:end,:);


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

    %% 3D Visualization Functions

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