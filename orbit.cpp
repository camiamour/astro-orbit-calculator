#include <iostream>
#include <cmath>
using namespace std;

int main() {
    const double G = 6.674e-11;     // gravitational constant
    double M, r;
    cout << "Enter mass of central body (kg): ";
    cin >> M;
    cout << "Enter orbital radius (m): ";
    cin >> r;

    double v = sqrt(G * M / r);     // orbital velocity
    double T = 2 * M_PI * r / v;    // orbital period

    cout << "\nOrbital velocity: " << v << " m/s";
    cout << "\nOrbital period: " << T << " seconds\n";
    return 0;
}
