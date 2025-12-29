% ============== ecef2enu.m ==============
function [e, n, u] = ecef2enu(dx, dy, dz, lat, lon, ~)
% ECEF2ENU - 将ECEF坐标系下的向量转换为本地东北天(ENU)坐标系下的向量。
lat = deg2rad(lat);
lon = deg2rad(lon);

sin_lat = sin(lat);
cos_lat = cos(lat);
sin_lon = sin(lon);
cos_lon = cos(lon);

e = -sin_lon * dx + cos_lon * dy;
n = -sin_lat * cos_lon * dx - sin_lat * sin_lon * dy + cos_lat * dz;
u =  cos_lat * cos_lon * dx + cos_lat * sin_lon * dy + sin_lat * dz;
end