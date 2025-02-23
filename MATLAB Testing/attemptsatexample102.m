
close all
clear all
clc

% random comment for testing

%% Overall Parameters

hours = 3600;     % Hours to seconds
days  = 24*hours; % Days to seconds
deg   = pi/180;   % Degrees to radians


%% Orbital Initial Conditions & Physical Parameters

muBody = 398600;
rBody = 6378;

bodyJ2 = 1082.63e-6;   

bodyParams = [muBody;
               bodyJ2;
               rBody];%               rBody_equatorial];

%% Initial Conditions and Parameters

% Initial orbital parameters:
zp0 = 300;                       % Perigee altitude (km)
za0 = 3062;                      % Apogee altitude (km)
RA0 = 45*deg;                    % Right ascension of the node (radians)
i0  = 28*deg;                    % Inclination (radians)
w0  = 30*deg;                    % Argument of perigee (radians)
TA0 = 40*deg;                    % True anomaly (radians)
 
% Calculated values
rp0 = rBody + zp0;              % Perigee radius (km)
ra0 = rBody + za0;              % Apogee radius (km)
e0  = (ra0 - rp0)/(ra0 + rp0);  % Eccentricity
a0  = (ra0 + rp0)/2;            % Semimajor axis (km)
h0  = sqrt(rp0*muBody*(1+e0));  % Angular momentum (km^2/s)
T0  = 2*pi/sqrt(muBody)*a0^1.5; % Period (s)
 

%% Orbital Function Definitions

% Get the orbital parameters from the position and velocity
function out_RV = params2rv(bodyParams, orbitParams)

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

% Convert position and velocity to orbital parameters
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

    else

        A = 1;

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
    
    rVal = norm(rvVec(1:3));

    force_constant1 = (3./2).*(bodyJ2).*(muBody).*(rBody_equatorial^2).*(rVal^(-4));
    force_constant2 = 5.*((rvVec(3)./rVal).^2);

    force_x = force_constant1.*(rvVec(1)./rVal).*(force_constant2 - 1);
    force_y = force_constant1.*(rvVec(2)./rVal).*(force_constant2 - 1);
    force_z = force_constant1.*(rvVec(3)./rVal).*(force_constant2 - 3);

    out_oblateForce = [force_x;
                       force_y;
                       force_z];

end

% Sum forces and generate deviations
function out_orbitDisp = orbitDisp(t, pVec, rv_Current, tp, bodyParams, configVals)
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
    % Not yet included

    % Oblateness Effects
    force_oblateness = oblateForce(bodyParams, rv_Pert);

    % Solar Pressure
    % Not yet included

    % Combine all active forces
    distForce = force_oblateness; % + others;

    % Determine del_a and del_v
    r_K = norm(rv_Current(1:3));
    r_P = r_K + norm(pVec(1:3));

    F     = 1 - (r_K/r_P)^3;
    del_a = -(bodyParams(1)/(r_K^3))*(pVec(1:3) - F*rv_Pert(1:3)) + distForce;
     
    out_orbitDisp  = [pVec(4:6); del_a]; 

end



%% Method: Encke Integration + ode45

orbitParams = [h0, i0, RA0, e0, w0, TA0];

rv_Inits = params2rv(bodyParams, orbitParams);
 
t0      = 0; 
tF      = 2*days;                   % Initial and final time (s)
tS      = t0;                       % Output vector of times
del_t   = T0/100;                   % Time step for Encke procedure
options = odeset('maxstep', del_t); % Tolerance


% Determine the osculating orbit at the end of the timestep
maxSteps = 100;
tolVal = 1.e-5;
uano_Params = [maxSteps; tolVal];

% Initial Perturbation
pVec = [0;0;0;0;0;0];

% Container for rvOutputs and times
actualOrbit = rv_Inits;

% First timestep
tp = t0;
tC = t0 + del_t;
 
% Consider Sperling Burdet
while tC <= tF + del_t/2
     
    % orbitDisp = rates
    [~,disp_derivs] = ode45(@(t,y) orbitDisp(t, y, rv_Inits, tp, bodyParams, uano_Params), [tp tC], pVec);

    % At the starting point, propagate the orbit as if there were no
    % perturbations at all, to determine where the orbiter would have been.
    [oscVals] = orbitUniVar(bodyParams, rv_Inits, del_t, uano_Params);

    % Add the changes in position and velocity to the ideal location, using
    % data from the previous step (or initial conditions)
    rvP = oscVals + disp_derivs(end,:)';

    % Save time
    tS = [tS,tC];

    if length(tS) >= 967

        A = 1;

    end

    % Start from new location
    rv_Inits = rvP;

    % Store data
    actualOrbit = [actualOrbit,rvP];

    % New timestep
    tp = tC;
    tC = tC + del_t;

    pVec = zeros(6,1);

end

 
%% Plotting

%...At each solution time extract the orbital elements from the state
%   vector using Algorithm 4.2:
n_times = length(tS);   %n_times is the number of solution times
for j = 1:n_times    
    r_t     = [actualOrbit(1:3,j)];
    v_t     = [actualOrbit(4:6,j)];
    r(j)    = norm(r_t);
    v(j)    = norm(v_t);
    coe     = rv2params(bodyParams, r_t, v_t);
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
plot(tS/3600,(raan - RA0)/deg)
title('Variation of Right Ascension')
xlabel('hours')
ylabel('{\it\Delta\Omega} (deg)')
grid on
grid minor
axis tight
 
subplot(2,1,2)
plot(tS/3600,(argp - w0)/deg)
title('Variation of Argument of Perigee')
xlabel('hours')
ylabel('{\it\Delta\omega} (deg)')
grid on
grid minor
axis tight
 