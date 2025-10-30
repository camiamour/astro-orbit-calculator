/*
  Astro Orbital Calculator
  Author: Camila Torres Juarez
  Purpose: Quick orbital calculations 
*/

#include <iostream>
#include <iomanip>
#include <cmath>
#include <string>
#include <vector>
#include <sstream>
#include <stdexcept>

using std::cin;
using std::cout;
using std::string;
using std::vector;

// ---------- Constants ----------
constexpr double G = 6.67430e-11; // m^3 kg^-1 s^-2 (CODATA 2018)
inline double PI() { return std::acos(-1.0); }
inline double TWO_PI() { return 2.0 * PI(); }

// ---------- Data Types ----------
struct Planet {
    string name;
    double mass_kg;
    double radius_m;
    double rot_period_s;
};

static const vector<Planet> PRESETS = {
    {"Earth", 5.972e24, 6.371e6, 86164.0905},
    {"Mars",  6.417e23, 3.3895e6, 88642.6848},
    {"Moon",  7.342e22, 1.7374e6, 2360584.685},
    {"Sun",   1.989e30, 6.9634e8, 2160000.0}
};

// ---------- Helpers ----------
double askNumber(const string& prompt) {
    while (true) {
        cout << prompt;
        string s;
        if (!std::getline(cin, s)) throw std::runtime_error("Input closed.");
        std::stringstream ss(s);
        double x;
        if (ss >> x) return x;
        cout << "  [Warning] Please enter a numeric value.\n";
    }
}

double km_to_m(double km) { return km * 1000.0; }
double m_to_km(double m)  { return m / 1000.0; }

// ---------- Orbital Mechanics ----------
double circSpeed(double M, double r) {
    if (M <= 0 || r <= 0) throw std::invalid_argument("Mass and radius must be > 0.");
    return std::sqrt(G * M / r);
}

double orbitPeriod(double M, double r) {
    if (M <= 0 || r <= 0) throw std::invalid_argument("Mass and radius must be > 0.");
    return TWO_PI() * std::sqrt((r*r*r) / (G * M));
}

double escapeSpeed(double M, double r) {
    if (M <= 0 || r <= 0) throw std::invalid_argument("Mass and radius must be > 0.");
    return std::sqrt(2.0 * G * M / r);
}

struct Hohmann {
    double dv1;
    double dv2;
    double total;
    double t_half;
};

Hohmann hohmannTransfer(double M, double r1, double r2) {
    if (M <= 0 || r1 <= 0 || r2 <= 0) throw std::invalid_argument("M, r1, r2 must be > 0.");
    const double mu = G * M;

    const double v1 = std::sqrt(mu / r1);
    const double v2 = std::sqrt(mu / r2);

    const double a  = 0.5 * (r1 + r2);
    const double v_peri = std::sqrt(mu * (2.0/r1 - 1.0/a));
    const double v_apo  = std::sqrt(mu * (2.0/r2 - 1.0/a));

    Hohmann out{};
    out.dv1   = std::fabs(v_peri - v1);
    out.dv2   = std::fabs(v2 - v_apo);
    out.total = out.dv1 + out.dv2;
    out.t_half = PI() * std::sqrt((a*a*a) / mu);
    return out;
}

double synchronousOrbitalRadius(double M, double rot_period_s) {
    if (M <= 0 || rot_period_s <= 0) throw std::invalid_argument("M and period must be > 0.");
    const double mu = G * M;
    return std::cbrt((rot_period_s * rot_period_s * mu) / (TWO_PI() * TWO_PI()));
}

// ---------- UI ----------
void displayHeader(const Planet& p) {
    cout << "\n=============================================\n";
    cout << "  Astro Orbital Toolkit  |  Central Body: " << p.name << "\n";
    cout << "=============================================\n\n";
}

Planet pickPlanet() {
    cout << "Choose a central body:\n";
    for (size_t i = 0; i < PRESETS.size(); ++i)
        cout << "  " << (i+1) << ") " << PRESETS[i].name << "\n";
    cout << "  " << PRESETS.size()+1 << ") Custom\n";

    int choice = 0;
    while (true) {
        double c = askNumber("Selection #: ");
        choice = static_cast<int>(c);
        if (choice >= 1 && choice <= static_cast<int>(PRESETS.size())+1) break;
        cout << "  [Warning] Invalid choice.\n";
    }

    if (choice <= static_cast<int>(PRESETS.size())) return PRESETS[choice-1];

    Planet p;
    p.name = "Custom";
    p.mass_kg = askNumber("Mass of central body (kg): ");
    p.radius_m = km_to_m(askNumber("Mean radius (km, optional): "));
    p.rot_period_s = askNumber("Rotation period (s, 0 if unknown): ");
    return p;
}

void doCircular(const Planet& p) {
    cout << "\n--- Circular Orbit ---\n";
    const double r_km = askNumber("Orbital radius from center (km): ");
    const double r = km_to_m(r_km);

    const double v = circSpeed(p.mass_kg, r);
    const double T = orbitPeriod(p.mass_kg, r);
    const double ve = escapeSpeed(p.mass_kg, r);

    cout << std::fixed << std::setprecision(3);
    cout << "\n@ r = " << r_km << " km around " << p.name << ":\n";
    cout << "  Circular speed:      " << v  << " m/s\n";
    cout << "  Orbital period:      " << T  << " s  (" << T/3600.0 << " hr)\n";
    cout << "  Escape velocity:     " << ve << " m/s\n";
}

void doHohmann(const Planet& p) {
    cout << "\n--- Hohmann Transfer ---\n";
    const double r1_km = askNumber("Initial orbit radius r1 (km): ");
    const double r2_km = askNumber("Target orbit radius r2 (km): ");
    const double r1 = km_to_m(r1_km), r2 = km_to_m(r2_km);

    const Hohmann h = hohmannTransfer(p.mass_kg, r1, r2);

    cout << std::fixed << std::setprecision(3);
    cout << "\nAround " << p.name << ":\n";
    cout << "  Δv1 (to transfer): " << h.dv1   << " m/s\n";
    cout << "  Δv2 (circularize): " << h.dv2   << " m/s\n";
    cout << "  Total Δv:          " << h.total << " m/s\n";
    cout << "  Transfer time:     " << h.t_half << " s (" << h.t_half/3600.0 << " hr)\n";
}

void doSynchronous(const Planet& p) {
    cout << "\n--- Synchronous Orbit ---\n";
    double T = p.rot_period_s;
    if (T <= 0) T = askNumber("Rotation period of body (s): ");
    else cout << "Using " << p.name << " rotation period: " << std::fixed << std::setprecision(3)
               << T << " s (" << T/3600.0 << " hr)\n";

    const double a = synchronousOrbitalRadius(p.mass_kg, T);
    const double h = (p.radius_m > 0) ? (a - p.radius_m) : -1.0;

    cout << std::fixed << std::setprecision(3);
    cout << "\nSynchronous radius (center): " << m_to_km(a) << " km\n";
    if (p.radius_m > 0)
        cout << "Altitude above surface:      " << m_to_km(h) << " km\n";
    else
        cout << "[Warning] Provide body radius to compute altitude above surface.\n";
}

int main() {
    std::ios::sync_with_stdio(false);
    cin.tie(nullptr);

    cout << "Welcome! Quick orbital calculator for demos + learning.\n";
    Planet planet = pickPlanet();
    displayHeader(planet);

    cout << std::fixed << std::setprecision(3);
    cout << "Planet info (approx):\n";
    cout << "  Mass:   " << planet.mass_kg << " kg\n";
    cout << "  Radius: " << m_to_km(planet.radius_m) << " km\n";
    if (planet.rot_period_s > 0)
        cout << "  Rot P:  " << planet.rot_period_s << " s (" << planet.rot_period_s/3600.0 << " hr)\n";
    else
        cout << "  Rot P:  (unknown)\n";

    while (true) {
        cout << "\nChoose an option:\n";
        cout << "  1) Circular orbit (v, T, v_escape)\n";
        cout << "  2) Hohmann transfer Δv (r1 -> r2)\n";
        cout << "  3) Synchronous orbit altitude\n";
        cout << "  4) Change central body\n";
        cout << "  0) Exit\n";

        int op = static_cast<int>(askNumber("Selection: "));
        try {
            if      (op == 1) doCircular(planet);
            else if (op == 2) doHohmann(planet);
            else if (op == 3) doSynchronous(planet);
            else if (op == 4) { planet = pickPlanet(); displayHeader(planet); }
            else if (op == 0) { cout << "\nGoodbye! Clear skies!\n"; break; }
            else cout << "  [Warning] Invalid choice.\n";
        } catch (const std::exception& e) {
            cout << "  [Warning] " << e.what() << "\n";
        }
    }
    return 0;
}

