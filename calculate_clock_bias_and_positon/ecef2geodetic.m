% ============== ecef2geodetic.m ==============
function [lat, lon, h] = ecef2geodetic(x, y, z)
% ECEF2GEODETIC - 将地心地固(ECEF)坐标转换为大地坐标(纬度, 经度, 高程)。
% WGS-84椭球体参数
a = 6378137.0;        % 赤道半径
f = 1/298.257223563;  % 扁率
e_sq = f * (2 - f);   % 偏心率的平方

lon = atan2(y, x);
p = sqrt(x^2 + y^2);

% 迭代计算纬度和高程
lat = atan2(z, p * (1 - e_sq));
h = 0;
N = a;
for i = 1:5 % 迭代5次足够
    sin_lat = sin(lat);
    N = a / sqrt(1 - e_sq * sin_lat^2);
    h = p / cos(lat) - N;
    lat = atan2(z, p * (1 - e_sq * N / (N + h)));
end

lat = rad2deg(lat);
lon = rad2deg(lon);
end