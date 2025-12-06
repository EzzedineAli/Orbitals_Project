clc;
clear;
close all;

mu = 398600;
Re = 6378 ;
tle = "tle.txt";
nilesat_new = "nilesat_new.txt";
nilesat = "NILESAT_201.txt";
galileo = "GSAT0101 (GALILEO-PFM).txt";
GPS = "GPS.txt";

filename = galileo;

[lon_cont, lat_cont, time_vec,r_vec,perigee_idx, perigee_time] = tle_to_latlong(filename);

% lon_cont =0; lat_cont=0;r_vec = [0;0;0];perigee_idx= 1; perigee_time = 0;

[sat_name, epoch_time,i, Omega_epoch, ecc, omega_epoch, M_epoch, n,a,T] = readsat(filename);
load("earth.mat")


now_time = juliandate(datetime('now', 'TimeZone', 'UTC'));
t_final = now_time*24*60*60 + T;
dt = 100;

[r_i,lambda_deg, phi_deg,t,Omega,omega ] = get_orbit(epoch_time,i, Omega_epoch, ecc, omega_epoch, M_epoch,n,a, t_final, dt);
Omega_deg = wrapTo360( rad2deg(Omega));
omega_deg = wrapTo360( rad2deg(omega));

figure

hold on
plot3(r_i(1,:),r_i(2,:),r_i(3,:))
plot3(r_vec(1,:),r_vec(2,:),r_vec(3,:),'-.')
plot3(r_i(1,1),r_i(2,1),r_i(3,1),'o')
plot3(r_vec(1,1),r_vec(2,1),r_vec(3,1),'o')
plot3(r_vec(1,perigee_idx),r_vec(2,perigee_idx),r_vec(3,perigee_idx),'s')
[xs, ys, zs ] = sphere;
xsr = xs*Re;
ysr = ys*Re;
zsr = zs*Re;
surf(xsr,ysr,zsr)
grid on
xlabel('X')
ylabel('Y')
zlabel('Z')
title("ECI Frame")
legend("Orbit of sat "+sat_name,"Matlab toolbox","Epoch Orbit of sat "+sat_name,"Epoch Matlab toolbox","Perigee MATLAB toolbox")
hold off


figure
hold on
plot(long,lat,'.')
plot(lon_cont,lat_cont,'.')
plot(lambda_deg,phi_deg,'.')
xlabel('long (deg)')
ylabel('lat (deg)')
title("ECI Frame")
legend("Earth","Matlab toolbox","Orbit of sat "+sat_name)
hold off



function [r_i,lambda_deg, phi_deg,t,Omega,omega] = get_orbit(epoch_time,i, Omega_epoch, ecc, omega_epoch, M_epoch,n,a, t_final, dt)
    mu = 398600;
    Re = 6378 ;
    j2 = 1.08263e-3;
    omega_earth = 2*pi/24/60/60;

    t_epoch_tp = M_epoch/n;
    tp = epoch_time - t_epoch_tp;
    n_t = (t_final-epoch_time)/dt;
    t = linspace(epoch_time,t_final,n_t);

    

    E = zeros(size(t));

    for j = 1:length(t)
        M = n*(t(j)-tp);
        kepler = @(E_val) E_val - ecc*sin(E_val) - M;
        E(j) = fzero(kepler,M);
    end
    
    nu = 2.*atan2(sqrt(1+ecc).*sin(E/2),sqrt(1-ecc)*cos(E/2));
    r = a.*(1-ecc^2)./(1+ecc.*cos(nu));
    r_p = [r.*cos(nu); r.*sin(nu); zeros(size(t))];
    
    Omegadot = -(3/2*sqrt(mu)*j2*Re^2/(1-ecc^2)^2/a^(7/2))*cos(i);
    omegadot = -(3/2*sqrt(mu)*j2*Re^2/(1-ecc^2)^2/a^(7/2))*(5/2*sin(i)^2-2);
    Omega = Omegadot*(t-epoch_time) +Omega_epoch;
    omega = omegadot*(t-epoch_time) +omega_epoch;
    r_i = zeros(3,length(t));
    for j = 1:length(t)
    
        T3O = [cos(Omega(j)), sin(Omega(j)), 0; ...
            -sin(Omega(j)), cos(Omega(j)), 0; ...
            0               0           1;];
        T1i = [ 1,   0,      0; ...
            0,  cos(i), sin(i);....
            0,  -sin(i), cos(i);];
    
        T3o = [cos(omega(j)), sin(omega(j)), 0; ...
            -sin(omega(j)), cos(omega(j)), 0; ...
            0               0           1;];
    
        Tip = (T3o * T1i * T3O).';
        r_i(:,j) = Tip * r_p(:,j);
    end
    theta = gmst_from_jd(t/24/60/60);
    x_e =  r_i(1,:).*cos(theta) + r_i(2,:).*sin(theta);
    y_e = -r_i(1,:).*sin(theta) + r_i(2,:).*cos(theta);
    z_e =  r_i(3,:);
    lambda = atan2(y_e, x_e);
    phi    = atan2(z_e, sqrt(x_e.^2 + y_e.^2));

    lambda_deg = rad2deg(lambda);
    phi_deg    = rad2deg(phi);
    % lambda = atan2(r_i(2,:),r_i(1,:));
    % phi = atan2(r_i(3,:),sqrt(r_i(1,:).^2+r_i(2,:).^2));
    % lambda_ecef= wrapToPi(lambda - omega_earth *(t-epoch_time));
    % lambda_deg = rad2deg(lambda_ecef);
    % phi_deg = rad2deg(phi);


end


function [sat_name, epoch_time,i, Omega_epoch, ecc, omega_epoch, M_epoch, n_sec,a,T] = readsat(filename)
    
    mu = 398600;
  
    sat_name = filename;
    orbit_vars = readmatrix(filename);
    epoch = orbit_vars(1,4);
    temp = num2cell(orbit_vars(2, 3:end-1));
    [i, Omega_epoch, ecc, omega_epoch, M_epoch, n] = temp{:};
    i = deg2rad(i);
    omega_epoch = deg2rad(omega_epoch);
    Omega_epoch = deg2rad(Omega_epoch);
    M_epoch = deg2rad(M_epoch);
    n_rad = 2*pi*n;
    n_sec = n_rad/24/60/60;
    ecc = ecc*1e-7;
    a = (mu/n_sec^2)^(1/3);
    T = 2*pi/n_sec;
        
    epoch_year = num2str(epoch, '%.0f');
    epoch_year = str2double(epoch_year(1:2));
    
    if epoch_year >= 57
        epoch_year = 1900 + epoch_year;
    else
        epoch_year = 2000 + epoch_year;
    end
    
    epoch_day = num2str(epoch, '%.10f');
    epoch_day = str2double(epoch_day(3:end));
    %% New part
    epoch_day_int = floor(epoch_day);
    epoch_day_frac = epoch_day - epoch_day_int;
    epoch_date = datetime(epoch_year,1,1,0,0,0,'TimeZone','UTC') + days(epoch_day_int-1) + seconds(epoch_day_frac*86400);
    epoch_jd = juliandate(epoch_date);
    epoch_time = epoch_jd*24*60*60;
    
    
    %% end of new part
    % epoch_time = (epoch_year*365+epoch_day-1)*24*60*60;

end

function theta = gmst_from_jd(JD)
    T = (JD - 2451545.0)/36525;

    theta = 280.46061837 ...
          + 360.98564736629 * (JD - 2451545.0) ...
          + 0.000387933 * T.^2 ...
          - (T.^3)/38710000;

    theta = deg2rad(mod(theta, 360));   % wrap to [0, 2pi)
end