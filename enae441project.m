%% ENAE441 Project - Group 15

%% Setup

clear; clc
load('opt2satDset3.mat')
load('opt3satDset3.mat')

% Observation data set (row indices)
idxs = [
        2   3   4
        6   8   10
        200 350 400
        500 600 700
        800 802 804
        840 841 842
];

%% Initial Orbit Determination

% Assuming observation site stays constant
lat = opt2satDset3.site_latitude_deg(1);
lon = opt2satDset3.site_longitude_deg(1);
alt = opt2satDset3.site_altitude_m(1);
lla_site = latlonalt_deg(lat, lon, alt);

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
    C0 = -mu^2*B^2;
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
    
    rv1 = [r1;v1];
    rv2 = [r2;v2];
    rv3 = [r3;v3];
    
    rv(:, end+1:end+3) = [r1 r2 r3; v1 v2 v3];
end
%% Deciding The Best Initial Estimate
% Select rows from night 2
night2obs_idxs = 1:4:100;
% Force model for this section
fm = force_model(4, 4, 0, 0, 1, 1, 1000);
% Create 1d vector from idxs
idxs_1d = idxs';
idxs_1d = idxs_1d(:);
RMS_list = [];
for i = 1:length(idxs_1d)
    rownum = idxs_1d(i);
    
    % Initial estimate
    initial_est_epoch = opt2satDset3.datetime(rownum);
    initial_est = pvt(initial_est_epoch, rv(1:3,i), rv(4:6,i));
    % Propagate the orbit, skip if orbit radius too small (likely because points for IOD too close together)
    try
        eph = propagate(initial_est, initial_est_epoch, initial_est_epoch+hours(36), minutes(30), fm);
    catch exception
        RMS_list(end+1) = NaN;
        continue
    end
    % Select rows from night 2
    cols = ["right_ascension_deg", "declination_deg",...
    "site_latitude_deg", "site_longitude_deg", "site_altitude_m"];
    night2obs = select_columns(opt3satDset3, cols, true);
    night2obs = night2obs(night2obs_idxs, :);
    % Interpolate the propagated orbits to selected observation times
    eph_night2est = ephemeris_interp(eph, datetime_iso8601(night2obs.datetime));
    
    % convert to aer
    aer_day2est = aer(eph_night2est, lla_site);
    % Calculate RMS
    RMS_temp = 0;
    N = height(night2obs);
    for j = 1:N
        az_pred = aer_day2est.azimuth_deg(j);
        az = night2obs.azimuth_deg(j);
        az_diff = min(abs(az_pred-az), 360-abs(az_pred-az));
        el_pred = aer_day2est.elevation_deg(j);
        el = night2obs.elevation_deg(j);
        el_diff = min(abs(el_pred-el), 360-abs(el_pred-el));
        RMS_temp = RMS_temp + az_diff^2 + el_diff^2;
    end
    RMS_list(end+1) = sqrt(RMS_temp / N);
end
% Get index of min RMS
[~, i] = min(RMS_list);
rv_best = rv(:, i)
rownum_best = idxs_1d(i);
% Time to use for initial estimate
dt_best = opt2satDset3.datetime(rownum_best);
%% Estimation Setup

cols = ["right_ascension_deg", "declination_deg",...
"site_latitude_deg", "site_longitude_deg", "site_altitude_m"];
night1obs = select_columns(opt2satDset3, cols, true);
night2obs = select_columns(opt3satDset3, cols, true);

% site location
oapchile = make_station("OAP-Chile", lat, lon, alt);

% Initial estimate
initial_est = pvt(dt_best, rv_best(1:3), rv_best(4:6));

% Define a force model to use initially
fm = force_model(4, 4, 0, 0, 1, 1, 1000);

%% Varying data points
% Make subsets of the original data to test over
night1samples = {};
night1samples{1} = night1obs(1:4:100,:);
night1samples{2} = night1obs(371:4:470,:);
night1samples{3} = night1obs(747:4:847,:);
night1samples{4} = night1obs(1:16:800,:);
% Night 2 sample to use for RMS calculations
night2sample = night2obs(1:4:100,:);

RMS_list = [];
sampling_est = initial_est;
sampling_epoch = dt_best;
for i = 1:length(night1samples)
    % Find the orbit
    out = determine_orbit(sampling_est, oapchile, night1samples{i}, fm);
    %disp(out.details)
    % Initial estimate
    sampling_epoch = out.estimated.epoch;
    sampling_est = pvt(sampling_epoch, out.estimated.position_m, out.estimated.velocity_ms);
    % Propagate the orbit, skip if orbit radius too small (likely because points for IOD too close together)
    try
        eph = propagate(sampling_est, sampling_epoch, sampling_epoch+hours(36), minutes(30), fm);
    catch exception
        RMS_list(end+1) = NaN;
        continue
    end
    % Interpolate the propagated orbits to selected observation times
    eph_night2est = ephemeris_interp(eph, datetime_iso8601(night2sample.datetime));
    
    % convert to aer
    aer_day2est = aer(eph_night2est, lla_site);
    % Calculate RMS
    RMS_temp = 0;
    N = height(night2sample);
    for j = 1:N
        az_pred = aer_day2est.azimuth_deg(j);
        az = night2obs.azimuth_deg(j);
        az_diff = min(abs(az_pred-az), 360-abs(az_pred-az));
        el_pred = aer_day2est.elevation_deg(j);
        el = night2obs.elevation_deg(j);
        el_diff = min(abs(el_pred-el), 360-abs(el_pred-el));
        RMS_temp = RMS_temp + az_diff^2 + el_diff^2;
    end
    RMS_list(end+1) = sqrt(RMS_temp / N);
end
[~, i] = min(RMS_list);
% sample to use in other tests
best_sample = night1samples{i};
%% Varying force model degree/order

% initial force models to determine which order/degree works best
OD_force_models = [constants.force_twobody, force_model(2, 0, 0, 0, 1, 1, 1000), force_model(2, 2, 0, 0, 1, 1, 1000), force_model(20, 20, 0, 0, 1, 1, 1000)];

RMS_list = [];
OD_fm_est = initial_est;
OD_fm_epoch = dt_best;
for k = 1:length(OD_force_models)
    % find the orbit
    orb = determine_orbit(OD_fm_est, oapchile, best_sample, OD_force_models(k));
    % Initial estimate
    OD_fm_epoch = out.estimated.epoch;
    OD_fm_est = pvt(OD_fm_epoch, out.estimated.position_m, out.estimated.velocity_ms);
    % Propagate the orbit, skip if orbit radius too small (likely because points for IOD too close together)
    try
        eph = propagate(OD_fm_est, OD_fm_epoch, OD_fm_epoch+hours(36), minutes(30), fm);
    catch exception
        RMS_list(end+1) = NaN;
        continue
    end
    % Interpolate the propagated orbits to selected observation times
    eph_night2est = ephemeris_interp(eph, datetime_iso8601(night2sample.datetime));
    
    % convert to aer
    aer_day2est = aer(eph_night2est, lla_site);
    % Calculate RMS
    RMS_temp = 0;
    N = height(night2sample);
    for j = 1:N
        az_pred = aer_day2est.azimuth_deg(j);
        az = night2obs.azimuth_deg(j);
        az_diff = min(abs(az_pred-az), 360-abs(az_pred-az));
        el_pred = aer_day2est.elevation_deg(j);
        el = night2obs.elevation_deg(j);
        el_diff = min(abs(el_pred-el), 360-abs(el_pred-el));
        RMS_temp = RMS_temp + az_diff^2 + el_diff^2;
    end
    RMS_list(end+1) = sqrt(RMS_temp / N);
end
[~, i] = min(RMS_list);
% best force model
best_fm = OD_force_models(i);
%% Try force models based on other params

% TODO - run for different luni-solar, then try 3 sets of solar radiatoin
% pressure parameters

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
