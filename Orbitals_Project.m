clc; clear vars; close all;
function Create_Orbit(Orbit_in)
    RE = 6738;
    mu = 398600;
    load("earth.mat", 'long', 'lat');
    a = Orbit_in{1};
    e = Orbit_in{2};
    TrueAnomaly_0 = Orbit_in{3};
    Inclination_0 = Orbit_in{4};
    RightAscension_0 = Orbit_in{5};
    ArgumentOfPerigee_0 = Orbit_in{6};
    Orbit_Name = Orbit_in{7};
    j2 = 1.08263e-3;
    p = a * (1 - e^2);
    T = 2 * pi * sqrt(a^3 / mu);
    n = 2 * pi / T;
    E_0 = 2 * atan( sqrt((1 - e) / (1 + e)) * tan(TrueAnomaly_0 / 2) );
    t_0 = (1 / n) * (E_0 - e * sin(E_0));
    t_vec = linspace(0, 2 * T, 1e3).';
    N = length(t_vec);
    E_vec = zeros(N,1);
    r_vec = zeros(N,1);
    x_vec = zeros(N,1);
    y_vec = zeros(N,1);
    z_vec = zeros(N,1);
    DCMbi_vec = cell(N,1);
    r_inertial_vec = zeros(N,3);
    RightAscension_vec = zeros(N,1);
    ArgumentOfPerigee_vec = zeros(N,1);
    for i = 1:N
        if i == 1
            E_vec(i) = fzero(@(E) getEQN(E, t_vec(i), t_0, e, n), 0);
        else
            E_vec(i) = fzero(@(E) getEQN(E, t_vec(i), t_0, e, n), E_vec(i - 1));
        end
    end
    function EQN_out = getEQN(E, t_1, t_0, e, n)
        EQN_out = E - e * sin(E) - n * (t_1 - t_0);
    end
    TrueAnomaly_vec = 2 * atan( sqrt((1+e)/(1-e)) * tan(E_vec./2) );
    RightAscension_dot = -(3/2) * n * j2 * (RE/p)^2 * cos(Inclination_0);
    ArgumentOfPerigee_dot = (3/4) * n * j2 * (RE/p)^2 * (4 - 5*sin(Inclination_0)^2);
    DaysPerYear = 365.25;
    Omega_Earth = 2*pi/(24*60*60*DaysPerYear/366.25);
    for i = 1:N
        RightAscension_vec(i) = RightAscension_0 + RightAscension_dot * t_vec(i) - Omega_Earth * t_vec(i);
        ArgumentOfPerigee_vec(i) = ArgumentOfPerigee_0 + ArgumentOfPerigee_dot * t_vec(i);
        RA = RightAscension_vec(i);
        AP = ArgumentOfPerigee_vec(i);
        IN = Inclination_0;
        DCMbi_vec{i} = [ cos(RA)*cos(AP) - sin(RA)*sin(AP)*cos(IN), -cos(RA)*sin(AP) - sin(RA)*cos(AP)*cos(IN), sin(RA)*sin(IN);
            sin(RA)*cos(AP) + cos(RA)*sin(AP)*cos(IN), -sin(RA)*sin(AP) + cos(RA)*cos(AP)*cos(IN), -cos(RA)*sin(IN);
            sin(AP)*sin(IN), cos(AP)*sin(IN), cos(IN) ];
    end
    for i = 1:N
        r_vec(i) = a*(1-e^2)/(1 + e*cos(TrueAnomaly_vec(i)));
        x_vec(i) = r_vec(i) * cos(TrueAnomaly_vec(i));
        y_vec(i) = r_vec(i) * sin(TrueAnomaly_vec(i));
        z_vec(i) = 0;
        r_inertial_vec(i,1:3) = DCMbi_vec{i} * [x_vec(i); y_vec(i); z_vec(i)];
    end
    x_inertial_vec = r_inertial_vec(:,1);
    y_inertial_vec = r_inertial_vec(:,2);
    z_inertial_vec = r_inertial_vec(:,3);
    R_inertial = sqrt(x_inertial_vec.^2 + y_inertial_vec.^2 + z_inertial_vec.^2);
    phi = asin(z_inertial_vec ./ R_inertial);
    gamma = atan2(y_inertial_vec, x_inertial_vec);
    figure(Name=Orbit_Name);
    tiledlayout(3,2,"TileSpacing","compact","Padding","compact");
    nexttile(1);
    plot(t_vec, E_vec);
    xlabel('t (s)'); ylabel('E (rad)');
    title('Eccentric Anomaly E');
    nexttile(3);
    plot(t_vec, TrueAnomaly_vec);
    xlabel('t (s)'); ylabel('$\theta$ (rad)', 'Interpreter','latex');
    title('True Anomaly $\theta$ (rad)', 'Interpreter','latex');
    nexttile(5);
    plot(x_vec, y_vec);
    xlabel('X (km)'); ylabel('Y (km)');
    title('2D Orbit Position');
    nexttile(6);
    plot3(x_inertial_vec, y_inertial_vec, z_inertial_vec);
    xlabel('X (km)'); ylabel('Y (km)'); zlabel('Z (km)');
    title('3D Orbit (ECI)');
    grid on;
    nexttile([2 1]);
    hold on;
    plot(long, lat, '.r');
    plot(gamma*180/pi, phi*180/pi, '.b');
    xlabel('Longitude (deg)');
    ylabel('Latitude (deg)');
    title(Orbit_Name);
    axis equal;
    grid on;
end

function Orbit_Out = Parse_TLE(Orbit_File)
    Orbit_File_lines = splitlines(Orbit_File);
    Orbit_Name = Orbit_File_lines{1};
    Orbit_TLE_1 = Orbit_File_lines{2};
    Orbit_TLE_2 = Orbit_File_lines{3};
    Inclination_0 = str2double(Orbit_TLE_2(9:16));
    RightAscension_0 = str2double(Orbit_TLE_2(18:25));
    e = str2double(Orbit_TLE_2(27:33)) * 1e-7;
    ArgumentOfPerigee_0 = str2double(Orbit_TLE_2(35:42));
    MeanAnomaly_deg = str2double(Orbit_TLE_2(44:51));
    n_rev_per_day = str2double(Orbit_TLE_2(53:63));
    MeanAnomaly = deg2rad(MeanAnomaly_deg);
    n = n_rev_per_day * 2*pi / 86400;
    TrueAnomaly_0 = 2*atan(sqrt((1+e)/(1-e)) * tan(MeanAnomaly/2));
    mu = 398600;
    a = (mu / n^2)^(1/3);
    Orbit_Out{1} = a;
    Orbit_Out{2} = e;
    Orbit_Out{3} = TrueAnomaly_0;
    Orbit_Out{4} = deg2rad(Inclination_0);
    Orbit_Out{5} = deg2rad(RightAscension_0);
    Orbit_Out{6} = deg2rad(ArgumentOfPerigee_0);
    Orbit_Out{7} = Orbit_Name;
end

function Orbit_Out = RV_Orbit(Orbit_File)
    Orbit_File_lines = splitlines(Orbit_File);
    Orbit_Name = Orbit_File_lines{1};
    Orbit_TLE_1 = py.str(Orbit_File_lines{2});
    Orbit_TLE_2 = py.str(Orbit_File_lines{3});
    Satalite = py.sgp4.api.Satrec();
    Satalite = Satalite.twoline2rv(Orbit_TLE_1, Orbit_TLE_2);
    Julian_Date = Satalite.jdsatepoch;
    Julian_Date_Fraction = Satalite.jdsatepochF;
    Output_sgp4 = Satalite.sgp4(Julian_Date, Julian_Date_Fraction);
    r_sgp4 = Output_sgp4{2};
    v_sgp4 = Output_sgp4{3};
    r = double([r_sgp4{1}, r_sgp4{2}, r_sgp4{3}]);
    v = double([v_sgp4{1}, v_sgp4{2}, v_sgp4{3}]);
    mu = 398600;
    r_mag = norm(r);
    v_mag = norm(v);
    h = cross(r, v);
    h_mag = norm(h);
    Inclination_0 = acos(h(3)/h_mag);
    k = [0 0 1];
    n = cross(k, h);
    n_mag = norm(n);
    RightAscension_0 = acos(n(1)/n_mag);
    if n(2) < 0
        RightAscension_0 = 2*pi - RightAscension_0;
    end
    e = (cross(v, h)/mu) - r/r_mag;
    e_mag = norm(e);
    ArgumentOfPerigee_0 = acos(dot(n, e) / (n_mag * e_mag));
    if e(3) < 0
        ArgumentOfPerigee_0 = 2*pi - ArgumentOfPerigee_0;
    end
    TrueAnomaly_0 = acos(dot(e, r) / (e_mag * r_mag));
    if dot(r, v) < 0
        TrueAnomaly_0 = 2*pi - TrueAnomaly_0;
    end
    a = 1 / (2/r_mag - v_mag^2/mu);
    Orbit_Out{1} = a;
    Orbit_Out{2} = e_mag;
    Orbit_Out{3} = TrueAnomaly_0;
    Orbit_Out{4} = Inclination_0;
    Orbit_Out{5} = RightAscension_0;
    Orbit_Out{6} = ArgumentOfPerigee_0;
    Orbit_Out{7} = Orbit_Name;
end
Orbit_File = fileread('GSAT0101_(GALILEO-PFM).txt');
Orbit = Parse_TLE(Orbit_File);
Create_Orbit(Orbit);
Orbit = RV_Orbit(Orbit_File);
Create_Orbit(Orbit);
Orbit_File = fileread('NILESAT_201.txt');
Orbit = Parse_TLE(Orbit_File);
Create_Orbit(Orbit);
Orbit = RV_Orbit(Orbit_File);
Create_Orbit(Orbit);

