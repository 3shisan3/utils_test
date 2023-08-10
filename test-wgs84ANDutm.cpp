#include <iostream>
#include <cmath>

// 定义一些常量
const double PI = 3.14159265358979;
const double a = 6378137.0; // WGS84椭球体长半轴
const double b = 6356752.3142; // WGS84椭球体短半轴
const double k0 = 0.9996; // UTM投影的比例因子
const double FE = 500000.0; // False Easting
const double FN = 0.0; // False Northing
const double e = sqrt(1 - pow(b, 2) / pow(a, 2)); // 椭球体第一偏心率
const double e2 = pow(e, 2); // 椭球体第二偏心率的平方

// 定义一些辅助函数
double degToRad(double deg) {
    return deg * PI / 180.0;
}

double radToDeg(double rad) {
    return rad * 180.0 / PI;
}

// WGS84坐标系转UTM坐标系
void wgs84ToUtm(double latitude, double longitude, int& zone, double& easting, double& northing) {
    double latRad = degToRad(latitude);
    double lonRad = degToRad(longitude);

    // 计算UTM投影中央经线的带号
    zone = static_cast<int>((longitude + 180.0) / 6.0) + 1;

    // 计算相关参数
    double lonOrigin = degToRad((zone - 1) * 6 - 180 + 3); // 中央经线的起始经度
    double N = a / sqrt(1 - e2 * pow(sin(latRad), 2)); // 卯酉圈曲率半径
    double T = pow(tan(latRad), 2);
    double C = e2 / pow(1 - e2, 2) * pow(cos(latRad), 2);
    double A = (lonRad - lonOrigin) * cos(latRad);
    double M = a * ((1 - e2 / 4 - 3 * pow(e2, 2) / 64 - 5 * pow(e2, 3) / 256) * latRad -
                    (3 * e2 / 8 + 3 * pow(e2, 2) / 32 + 45 * pow(e2, 3) / 1024) * sin(2 * latRad) +
                    (15 * pow(e2, 2) / 256 + 45 * pow(e2, 3) / 1024) * sin(4 * latRad) -
                    (35 * pow(e2, 3) / 3072) * sin(6 * latRad));

    // 计算UTM坐标
    easting = k0 * N * (A + (1 - T + C) * pow(A, 3) / 6 +
                       (5 - 18 * T + pow(T, 2) + 72 * C - 58 * e2) * pow(A, 5) / 120) + FE;
    northing = k0 * (M + N * tan(latRad) * (pow(A, 2) / 2 +
                                            (5 - T + 9 * C + 4 * pow(C, 2)) * pow(A, 4) / 24 +
                                            (61 - 58 * T + pow(T, 2) + 600 * C - 330 * e2) * pow(A, 6) / 720)) + FN;
}

// UTM坐标系转WGS84坐标系
void utmToWgs84(int zone, double easting, double northing, double& latitude, double& longitude) {
    double x = easting - FE;
    double y = northing - FN;

    // 计算相关参数
    double lonOrigin = degToRad((zone - 1) * 6 - 180 + 3); // 中央经线的起始经度
    double M = y / k0;
    double mu = M / (a * (1 - e2 / 4 - 3 * pow(e2, 2) / 64 - 5 * pow(e2, 3) / 256));
    double e1 = (1 - sqrt(1 - e2)) / (1 + sqrt(1 - e2));
    double C1 = e2 / (1 - e2) * pow(cos(mu), 2);
    double T1 = pow(tan(mu), 2);
    double N1 = a / sqrt(1 - e2 * pow(sin(mu), 2));
    double R1 = a * (1 - e2) / pow(1 - e2 * pow(sin(mu), 2), 1.5);
    double D = x / (N1 * k0);

    // 计算纬度
    latitude = mu - (N1 * tan(mu) / R1) *
               (pow(D, 2) / 2 - (5 + 3 * T1 + 10 * C1 - 4 * pow(C1, 2) - 9 * e2) * pow(D, 4) / 24 +
                (61 + 90 * T1 + 298 * C1 + 45 * pow(T1, 2) - 252 * e2 - 3 * pow(C1, 2)) * pow(D, 6) / 720);

    // 计算经度
    longitude = lonOrigin + (D - (1 + 2 * T1 + C1) * pow(D, 3) / 6 +
                             (5 - 2 * C1 + 28 * T1 - 3 * pow(C1, 2) + 8 * e2 + 24 * pow(T1, 2)) * pow(D, 5) / 120) / cos(mu);
    latitude = radToDeg(latitude);
    longitude = radToDeg(longitude);
}

int main() {
    // 示例：将WGS84坐标(39.9087, 116.3975)转换为UTM坐标
    double latitude = 39.9087;
    double longitude = 116.3975;
    int zone;
    double easting, northing;

    wgs84ToUtm(latitude, longitude, zone, easting, northing);

    // 输出结果
    std::cout << "UTM Zone: " << zone << std::endl;
    std::cout << "Easting: " << easting << std::endl;
    std::cout << "Northing: " << northing << std::endl;

    // 示例：将UTM坐标(32649, 4421411)转换为WGS84坐标
    // int utmZone = 49;
    // double utmEasting = 32649;
    // double utmNorthing;
    double wgsLatitude, wgsLongitude;

    // utmNorthing = 4421411;

    utmToWgs84(zone, easting, northing, wgsLatitude, wgsLongitude);

    // 输出结果
    std::cout << "Latitude: " << wgsLatitude << std::endl;
    std::cout << "Longitude: " << wgsLongitude << std::endl;

    return 0;
}
