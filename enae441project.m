%% ENAE441 Project - Group 15

%% Setup

clear; clc
load('opt2satDset3.mat')
load('opt3satDset3.mat')

%% Initial Orbit Determination

% Assuming observation site stays constant
lat = opt2satDset3.site_latitude_deg(1);
lon = opt2satDset3.site_longitude_deg(1);
alt = opt2satDset3.site_altitude_m(1);
lla_site = latlonalt_deg(lat, lon, alt);

idxs = [2   3   4
        200 500 750
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
    
    rv(:, end+1:end+3) = [rv1 rv2 rv3];
     
    % rv(:, end+1:end+3) = [r1 r2 r3 v1 v2 v3]';
end

disp(rv)

%%Deciding The Best Orbit

N = [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32].*31; %This is array of indicies that are being picked out from the second dataset


%The following for loop is picking out the datetimes from the secoond data set and placing them in their own array

for idx=1:width(N)

    datetimes_of_obs(idx,1)=datetime_iso8601(opt3satDset3.datetime(N(idx))); %used the datetime_iso8601() function because of a time zone error


end


%The following for loops are to propogated through each estimated state vector to find out the RMS for each state vector


for k=1:9

    r=[rv(1,k);rv(2,k);rv(3,k)]; %position vector for a specific index
    v=[rv(4,k);rv(5,k);rv(6,k)]; %velocity vector for a specific index

    pvt_of_init_orbit=pvt(opt2satDset3.datetime(idxs(k)),r,v); %pvt of the initial orbit at estimation time


    for w=1:32


        if opt3satDset3.azimuth_deg(N(w))>180

        obs_azimuth(w,1)= opt3satDset3.azimuth_deg(N(w))-180; %observation azimuth array of each N index

        else

        obs_azimuth(w,1)= opt3satDset3.azimuth_deg(N(w));

        end


        if opt3satDset3.elevation_deg(N(w))>180

        obs_elevation(w,1)=opt3satDset3.elevation_deg(N(w))-180; %observation elevation array of each N index

        else

        obs_elevation(w,1)=opt3satDset3.elevation_deg(N(w));

        end

        %insert force model in propogation when we Healy tells us what it is

        propogated_orbit=propagate(pvt_of_init_orbit,opt2satDset3.datetime(idxs(k)),hours(1),hours(1)); %Propogating initial orbit of the specific N index 

       
        interp_state_eci=ephemeris_interp(propogated_orbit,datetimes_of_obs(w)); %interpolating initial orbit to the specified datetimes

        interp_state_aer=aer(interp_state_eci,lla_site); %turning the interpolated state to aer format


  if interp_state_aer.azimuth_deg>180

        pred_azimuth(w,1)=interp_state_aer.azimuth_deg-180; %placing interpolated azimuth in the predicted azimuth array

  else

        pred_azimuth(w,1)=interp_state_aer.azimuth_deg; 


  end



  if interp_state_aer.elevation_deg>180

        pred_elevation(w,1)=interp_state_aer.elevation_deg-180; %placing interpolated elevation in the predicted elevation aray
  else

      pred_elevation(w,1)=interp_state_aer.elevation_deg;

  end

    end

RMS(k,1)=sqrt((1/width(N))*sum((obs_azimuth-pred_azimuth).^2+(obs_elevation-pred_elevation).^2)); %calculating the RMS


end


[value,index]=min(RMS); %finding the lowest RMS

best_rv=rv(:,index); %setting the rv variable equal to the orbit with the lowest RMS



%% Estimation
cols = ["right_ascension_deg", "declination_deg",...
"site_latitude_deg", "site_longitude_deg", "site_altitude_m"];
myobs = select_columns(opt2satDset3, cols, true);
% Define a force model
force_model = force_model(4, 4, 0, 0, 1, 1, 1000);
% site location
oapchile = make_station("OAP-Chile", lat, lon, alt);
% Make a subset of the original data
night1_early_25pts = myobs(1:4:100,:);

% ***** TEMP VALUES *****
% Time to use for initial estimate
initial_est_epoch = opt2satDset3.datetime(4);
% Initial estimate, currently based on the 4th entry of the estimates
initial_est = pvt(initial_est_epoch, rv(1:3,4), rv(4:6,4));

% find the orbit
runnum = determine_orbit(initial_est, oapchile, night1_early_25pts, force_model);
% propagate the orbit
orbProp = propagate(initial_est,initial_est_epoch,initial_est_epoch+hours(1),10,force_model);
% interpolate the propagated orbits to selected observation times
ei = ephemeris_interp(orbProp,orbProp.epoch); % NOTE: used the same times as before, this is not what we want to do

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
