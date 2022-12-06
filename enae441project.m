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

setnum = 1;

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
end

M = L\R;

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