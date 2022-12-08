%% ENAE441 Project - Group 15

%% Setup

clear all; clc
load('opt2satDset3');
load('opt3satDset3');


%% Initial Orbit Determination

% Assuming observation site stays constant
lat = opt2satDset3.site_latitude_deg(1);
lon = opt2satDset3.site_longitude_deg(1);
alt = opt2satDset3.site_altitude_m(1);
lla_site = latlonalt_deg(lat, lon, alt);

idxs = [2   3   4
        6   8   10
        200 350 400
        500 600 700
        800 802 804
        840 841 842];

rv = []; % to keep track of calculated orbits
for setnum = 1:height(idxs)
    L = [];
    R = [];
    for obsnum = 1:3
        i = idxs(setnum, obsnum);
    
        % Get ECI vector to observation site
        lla_site.epoch = opt2satDset3.datetime(i);
        Rvec = eci(lla_site).position_m';
        
        % Get AER to satellite from observation site
        az = opt2satDset3.azimuth_deg(i);
        el = opt2satDset3.elevation_deg(i);
        aerhat = azelrn_deg(az, el, 1);
        
        % Get rhohat
        aerhat.epoch = opt2satDset3.datetime(i);
        eci_dir = eci(aerhat, lla_site).position_m';
        rhohat = eci_dir - Rvec;
        
        L(:, obsnum) = rhohat;
        R(:, obsnum) = Rvec;
        dt(obsnum) = opt2satDset3.datetime(i);
    end
    
    M = L\R;
    
    % Calculate a's
    tau1 = seconds(dt(1) - dt(2));
    tau3 = seconds(dt(3) - dt(2));
    
    % eq 7.83 in text
    a1 = tau3/(tau3-tau1);
    a3 = -tau1/(tau3-tau1);
    a1u = tau3*((tau3-tau1)^2 - tau3^2) / (6*(tau3-tau1));
    a3u = -tau1*((tau3-tau1)^2 - tau1^2) / (6*(tau3-tau1));
    
    % eq. 7.88 in text
    A = M(2, 1)*a1 - M(2,2) + M(2, 3)*a3;
    B = M(2, 1)*a1u + M(2, 3)*a3u;
    
    E = dot(L(:, 2), R(:, 2));
    
    mu = 3.986e14; % m^3/s^2
    rE = 6378e3; % radius of the Earth, m
    
    % eq. 7.90
    C8 = 1;
    C6 = -(A^2 + 2*A*E + norm(R(:, 2))^2);
    C3 = -2*mu*B*(A + E);
    C0 = -mu^2*B^2;;
    p = [C8 0 C6 0 0 C3 0 0 C0];
    rts = roots(p);
    
    r2norm = rts(imag(rts) == 0 & rts > rE);
    
    u = mu / r2norm^3;
    
    % eq 7.101
    c1 = a1 + a1u*u;
    c2 = -1;
    c3 = a3 + a3u*u;
    
    % eq. 7.102
    csandrhos = M*[-c1; -c2; -c3];
    
    rho1 = csandrhos(1)/c1;
    rho2 = csandrhos(2)/c2;
    rho3 = csandrhos(3)/c3;
    
    r1 = R(:, 1) + rho1*L(:, 1);
    r2 = R(:, 2) + rho2*L(:, 2);
    r3 = R(:, 3) + rho3*L(:, 3);
    
    [v1, v2, v3] = gibbs(r1, r2, r3, mu);

    rv(:, end+1:end+3) = [r1 r2 r3 v1 v2 v3]';
end

disp(rv)

%% Functions

function [v1s, v2s, v3s] = gibbs(r1s, r2s, r3s, mu)
    r1 = norm(r1s);
    r2 = norm(r2s);
    r3 = norm(r3s);
    
    D = cross(r2s, r3s) + cross(r3s, r1s) + cross(r1s, r2s);
    N = r1*cross(r2s, r3s) + r2*cross(r3s, r1s) + r3*cross(r1s, r2s);
    
    p = norm(N)/norm(D);
    
    S = (r1-r2)*r3s + (r3-r1)*r2s + (r2-r3)*r1s;
    
    e = norm(S)/norm(D);
    
    w = N/norm(N);
    q = S/norm(S);
    p = cross(q, w);
    
    % finding v1
    B = cross(D, r1s);
    L = sqrt(mu/(norm(D)*norm(N)));
    v1s = (L/r1)*B + L*S;

    % finding v2
    B = cross(D, r2s);
    v2s = (L/r2)*B + L*S;

    % finding v3
    B = cross(D, r3s);
    v3s = (L/r3)*B + L*S;
end