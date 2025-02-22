
close all
clear all
clc

% random comment for testing

%% Overall Parameters

hours = 3600;     % Hours to seconds
days  = 24*hours; % Days to seconds
deg   = pi/180;   % Degrees to radians


%% Orbital Initial Conditions & Physical Parameters

rEarth = 6371.230; % Frédéric Chambat; Bernard Valette (2001) "Physics of the Earth and Planetary Interiors"
rEarth_equatorial = 6378.1370; % (kilometers) WGS84 Ellipsoid approximation
mEarth = 5.9722e24; % (kilogram) IAU Selected Astronomical Constants
gravConst = (6.67430e-20); % (km^3 kg^-1 s^-2) NIST Reference 2022 Codata
muEarth = (398600.4418); % (km^3 s^-2) Ries, J. C. 1992 "...Determination of the Gravitational Constant of the Earth" pp. 529-531

% % ISS Orbital Parameters 2025, Feb. 14
% orbitIncl = deg2rad(51.6371);   % radians
% orbitEcct = 0.0003865;          % unitless
% orbitRaan = deg2rad(193.6981);  % degrees
% orbitPaps = 414 + rEarth;       % kilometers
% orbitApog = 419 + rEarth;       % kilometers
% orbitArgp = deg2rad(315.5499);  % radians
% 
% orbitSemi = (0.5)*(orbitPaps + orbitApog);              % kilometers
% orbitAngm = sqrt(muEarth*orbitSemi*(1 - orbitEcct.^2)); % kg m^2 s^-1
% orbitPeri = 2*pi*sqr t( (orbitSemi.^3)/(muEarth) );      % seconds


% Zonal Harmonics from Curtis, Orbital Mechanics for Engineering Students
% Better approximations Shaub and Junkins (2009) p. 553 (for future)
earthJ2 = 0.00108263;            % First Zonal Harmonic


earthParams = [muEarth;
               earthJ2;
               rEarth;
               rEarth_equatorial];

%% Orbital Function Definitions

function out_orbitUniVar = orbitUniVar(earthParams, rVec, vVec, del_t, configVals)

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
    muEarth  = earthParams(1);

    % Format inputs for use
    rVal = norm(rVec);
    vVal = norm(vVec);
    vRad = dot(rVec,vVec)./rVal;

    % alphaVal < 0 --> Hyperbola
    % alphaVal = 0 --> Parabola
    % alphaVal > 0 --> Ellipse
    alphaVal = 2./rVal - ((vVal.^2)./muEarth);

    % Good estimate for initial uano (Chobotov, 2002) citation from Curtis
    uano = sqrt(muEarth).*abs(alphaVal).*del_t;
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
            sStumpff = ( uniParam.^(-3./2) ).*( sinh(sqrt(-uniParam)) - sqrt(-uniParam) );
    
        else
    
            cStumpff = 1./2;
            sStumpff = 1./6;
    
        end

        % See Curtis P. 171

        uano_c0 = (rVal).*(vRad).*(1./sqrt(muEarth));
        uano_c1 = (1 - alphaVal.*rVal);

        f_uano   = uano_c0.*(uniParam.^2).*cStumpff + uano_c1.*(uniParam.^3).*sStumpff + rVal.*uniParam - sqrt(muEarth).*del_t;
        df_duano = uano_c0.*uniParam.*(1 - alphaVal.*(uniParam.^2).*sStumpff) + uano_c1.*(uniParam.^2).*cStumpff + rVal;

        uano_ratio = f_uano./df_duano;

        uano = uano - uano_ratio;

    end

    % At the end of this while loop the Universal Anomaly is available in
    % the uano variable. The stumpff functions don't need to have any
    % memory so the most recent value is already stored. The universal
    % anomaly is placed into Lagrange coefficients for the calculation of
    % the osculating orbit.

    f_L  = 1 - ((uano.^2)./rVal).*cStumpff;
    g_L  = del_t - (1./sqrt(muEarth)).*(uano.^3).*sStumpff;

    rNew = f_L.*rVal + g_L.*vVal;

    df_L = (sqrt(muEarth)./(rVal.*rNew)).*(alphaVal.*(uano.^3).*sStumpff - uano);
    dg_L = 1 - ((uano.^2)./rNew).*cStumpff;

    vNew = df_L.*rVal + dg_L.*vVal;

    % This output is the osculating orbit after some time del_t. It does
    % not account for any perturbations.
    out_orbitUniVar = [rNew; vNew];


end

% Convert orbital parameters to RV (See Curtis p 191)
function out_RV = params2rv(earthParams, orbitParams)
    % orbitParams:
    % 
    % angm - Angular momentum
    % incl - Inclination
    % raan - Right Ascencion of the Ascending Node
    % ecct - Eccentricity
    % argp - Argument of Perigee
    % tano - True Anomaly

    % This should be converted to use quaternions in the future.

    muEarth = earthParams(1);

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
function out_orbitParams = rv2params(earthParams, rVec,vVec)
    % rVec & vVec are X, Y, Z in inertial frame

    % rVec - kilometers
    % vVec - kilometers/second

    muEarth = earthParams(1);

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

% Generate forces on orbiter due to earths oblateness
function out_oblateForce = oblateForce(earthParams, rVec)

    % Later on, it may be interesting to look at Shaub and Junkins (2009)
    % on page 553 for accelerations due to J3 -> J6. 

    muEarth = earthParams(1);
    earthJ2 = earthParams(2);
    rEarth_equatorial = earthParams(3);
    
    rVal = norm(rVec);

    force_constant1 = (3./2).*(earthJ2).*(muEarth).*(rEarth_equatorial^2).*(rVal^(-4));
    force_constant2 = 5.*((rVec(3)./rVal).^2);

    force_x = force_constant1.*(rVec(1)./rVal).*(force_constant2 - 1);
    force_y = force_constant1.*(rVec(2)./rVal).*(force_constant2 - 1);
    force_z = force_constant1.*(rVec(3)./rVal).*(force_constant2 - 3);

    F     = 1 - (rosc/rpp)^3;
    del_a = -mu/rosc^3*(del_r - F*Rpp) + ap;

    out_oblateForce = [force_x;
                       force_y;
                       force_z];

end

function out_orbitDeriv = orbitDeriv(~, rv_Input, earthParams)

    muEarth   = earthParams(1);

    rVec = rv_Input(1:3);
    vVec = rv_Input(4:6);
    rVal = norm(rVec);
    oParam = (muEarth./(rVal.^3));

    % Aerodynamic Drag
    % Not yet included

    % Oblateness Effects
    force_oblateness = oblateForce(earthParams, rVec);

    % Solar Pressure
    % Not yet included

    distForce = force_oblateness;

    rdot = vVec;
    vdot = -oParam.*rVec + distForce;

    out_orbitDeriv = [rdot;vdot];


end

% Tests to confirm accuracy, based on examples in Curtis Chapter 4
%orbitParams_results = rv2params(muEarth,[-6045000;-3490000;2500000],[-3457;6618;2533]); % m^3 s^-2, meters, m/s
%rv_Results  = params2rv(muEarth,[80000*(1000^2);deg2rad(30);deg2rad(40);1.4;deg2rad(60);deg2rad(30)]);


%% Initial Conditions and Parameters

% Initial orbital parameters:
zp0 = 300;                       % Perigee altitude (km)
za0 = 3062;                      % Apogee altitude (km)
RA0 = 45*deg;                    % Right ascension of the node (radians)
i0  = 28*deg;                    % Inclination (radians)
w0  = 30*deg;                    % Argument of perigee (radians)
TA0 = 40*deg;                    % True anomaly (radians)
 
% Calculated values
rp0 = rEarth + zp0;              % Perigee radius (km)
ra0 = rEarth + za0;              % Apogee radius (km)
e0  = (ra0 - rp0)/(ra0 + rp0);   % Eccentricity
a0  = (ra0 + rp0)/2;             % Semimajor axis (km)
h0  = sqrt(rp0*muEarth*(1+e0));  % Angular momentum (km^2/s)
T0  = 2*pi/sqrt(muEarth)*a0^1.5; % Period (s)


%% Method: ode45 only

% initR = [-2384.46; 5729.01; 3050.46];
% initV = [-7.36138; -2.98997; 1.64354];
% 
% oParam_init = rv2params(earthParams, initR, initV);
% [t_orbit,y_orbit] = ode45(@(t,y) orbitDeriv(t, y, earthParams), linspace(0, 48*hours, 1000), [initR;initV]);
% y_orbit = y_orbit';
% 
% figure(90)
% hold on
% grid on
% plot3(y_orbit(1,1),y_orbit(2,1),y_orbit(3,1),'ro')
% plot3(y_orbit(1,:),y_orbit(2,:),y_orbit(3,:),'k')
% axis equal
% view([1,1,.4])
% title('Testing Orbital Path')
% xlabel('x')
% ylabel('y')
% zlabel('z')
% 
% orbitParams_mat = zeros(size(y_orbit));
% for ind1 = 1:length(t_orbit)
% 
%     orbitParams_mat(:,ind1) = rv2params(earthParams, y_orbit(1:3,ind1), y_orbit(4:6,ind1)); 
% 
% end
% 
% figure(91)
% hold on
% grid on
% plot(t_orbit./hours,rad2deg(orbitParams_mat(3,:))-45,'k')
% title('Change in RAAN over time')
% xlabel('Hours')
% ylabel('RAAN in Deg')
% 
% figure(92)
% hold on
% grid on
% plot(t_orbit./hours,rad2deg(orbitParams_mat(5,:))-30,'k')
% title('Change in RAAN over time')
% xlabel('Hours')
% ylabel('RAAN in Deg')

 
%% Method: Encke Integration + ode45

% THIS MATCHES EXPECTATIONS IN THE OTHER PARTS, SO UP UNTIL HERE IT IS ACCURATE. 
orbitParams = [h0, i0, RA0, e0, w0, TA0];

rv_Inits = params2rv(earthParams, orbitParams);
 
rVec_i = rv_Inits(1:3);
vVec_i = rv_Inits(4:6);
                                
rVal_i = norm(rVec_i);
vVal_i = norm(vVec_i);    
 
t0      = 0; 
tf      = 2*days;                   % Initial and final time (s)
del_t   = T0/100;                   % Time step for Encke procedure
options = odeset('maxstep', del_t); % Tolerance
t       = t0;                       % Initialize the time scalar
tsave   = t0;                       % Initialize the vector of solution times
y       = [rVec_i; vVec_i];         % Initialize the state vector
del_r   = [0;0;0];                        % Initialize position deviation
del_v   = [0;0;0];                        % Initialize velocity deviation
del_y0  = [del_r, del_v];               % Initialize the state vector perturbation

 
t = t + del_t;                      %First time step
 
%   Loop over the time interval [t0, tf] with equal increments del_t:
while t <= tf + del_t/2
    
    % Simulate the total perturbation by the end of the timestep
     [~,z] = ode45(@(t, y) rates(t, y, earthParams), [t0 t], del_y0, options);
    
    % Determine the osculating orbit at the end of the timestep
    maxSteps = 1000;
    tolVal = 1.e-8;
    uano_Params = [maxSteps;
                   tolVal];
    [oscVals] = orbit_UniVar(earthParams, rVec_i, vVec_i, t-t0, uano_Params);
 
    % At the new point, add effects of orbital perturbations (rectify)

    rOsc = oscVals(1:3);
    vOsc = oscVals(4:6);

    rVec_i     = rOsc + z(end,1:3);
    vVec_i     = vOsc + z(end,4:6);

    % Prepare for next timestep
    t0         = t;
    tsave  = [tsave;t];
    y      = [y, [rVec_i; vVec_i]];
    t      = t + del_t;
    del_y0 = zeros(6,1); 
end

t = tsave;   % Set aside all the timesteps
 
%% Plotting

%...At each solution time extract the orbital elements from the state
%   vector using Algorithm 4.2:
n_times = length(t);   %n_times is the number of solution times
for j = 1:n_times    
    r_t     = [y(j,1:3)];
    v_t     = [y(j,4:6)];
    r(j)    = norm(r_t);
    v(j)    = norm(v_t);
    coe     = rv2params(earthParams, r_t, v_t);
    angm(j) = coe(1);
    incl(j) = coe(2);
    raan(j) = coe(3);
    ecct(j) = coe(4);
    argp(j) = coe(5);
    tano(j) = coe(6);
end
 
%...Plot selected osculating elements:
 
figure(1)
subplot(2,1,1)
plot(t/3600,(raan - RA0)/deg)
title('Variation of Right Ascension')
xlabel('hours')
ylabel('{\it\Delta\Omega} (deg)')
grid on
grid minor
axis tight
 
subplot(2,1,2)
plot(t/3600,(argp - w0)/deg)
title('Variation of Argument of Perigee')
xlabel('hours')
ylabel('{\it\Delta\omega} (deg)')
grid on
grid minor
axis tight
 
figure(2)
subplot(3,1,1)
plot(t/3600,angm - h0)
title('Variation of Angular Momentum')
xlabel('hours')
ylabel('{\it\Deltah} (km^2/s)')
grid on
grid minor
axis tight
 
subplot(3,1,2)
plot(t/3600,ecct - e0)
title('Variation of Eccentricity')
xlabel('hours')
ylabel('\it\Deltae')
grid on
grid minor
axis tight
 
subplot(3,1,3)
plot(t/3600,(incl - i0)/deg)
title('Variation of Inclination')
xlabel('hours')
ylabel('{\it\Deltai} (deg)')
grid on
grid minor
axis tight
 
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function dfdt = rates(t,f,earthParams)
%~~~~~~~~~~~~~~~~~~~~~~~~~
%
% This function calculates the time rates of Encke's deviation in position
% del_r and velocity del_v.
% -----------------------------------------------------------------------
 
del_r = f(1:3)';     %Position deviation
del_v = f(4:6)';     %Velocity deviation
 
%...Compute the state vector on the osculating orbit at time t
%   (Equation 12.5) using Algorithm 3.4:
[oscVals] = orbitUniVar(earthParams, del_r, del_v, t, [1000;1.e-8]);

rOsc = oscVals(1:3);
vOsc = oscVals(4:6);
 
%...Calculate the components of the state vector on the perturbed orbit
%   and their magnitudes:
Rpp   = rOsc + del_r;
Vpp   = vOsc + del_v;
rosc  = norm(rOsc);
rpp   = norm(Rpp);
 
%...Compute the J2 perturbing acceleration from Equation 12.30:
xx    = Rpp(1); 
yy    = Rpp(2); 
zz    = Rpp(3);
 
fac   =   3/2*J2*(mu/rpp^2)*(RE/rpp)^2;
ap    =   -fac*[(1 - 5*(zz/rpp)^2)*(xx/rpp) ...
                (1 - 5*(zz/rpp)^2)*(yy/rpp) ...
                (3 - 5*(zz/rpp)^2)*(zz/rpp)];
 
%...Compute the total perturbing acceleration from Equation 12.7:            
F     = 1 - (rosc/rpp)^3;
del_a = -mu/rosc^3*(del_r - F*Rpp) + ap;
 
dfdt  = [del_v(1) del_v(2) del_v(3) del_a(1) del_a(2) del_a(3)]';
dfdt  = [del_v del_a]'; %Return the deviative velocity and acceleration
                        %to ode45.
 
end 

 

