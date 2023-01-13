#include <iostream>
#include <cmath>
#include <fstream>

struct vec2d {double x, y;};

double dot_product (vec2d r, vec2d v) {
    return (r.x*v.x) + (r.y*v.y);
}

double magnitude (vec2d x) {
    return sqrt(pow(x.x,2)+pow(x.y,2));
}

double dvdt (double r, double u, double G, double M) {
    return (r*pow(u,2))-(G*M/pow(r,2));
}

double dudt (double u, double v, double r) {
    return -2.*u*v/r;
}

using namespace std;

int main() {

    int a, b, c, d, e;
    int n;

    cout << "Select 1 for asteroid collision, select 2 for planet collision." << endl;
    cin >> a;
    cout << endl;

    cout << "How many time steps would you like to take for this simulation (array size)?";
    cin >> n;
    cout << endl;

    float x_array[n], y_array[n];

    // Asteroid Collision
    if (a==1) {
        double omega, r, theta; // Initial conditions: angular velocity, radius, angular displacement
        double M, m; // Masses of planet and asteroid respectively
        double vx_initial, vy_initial, x_initial, y_initial; // Initial velocity and position of asteroid
        vec2d r_initial, v_initial;

        cout << "Specify the mass of the system-dominating planet in kg located at the origin." << endl;
        cin >> M;
        cout << endl;

        cout << "Please specify the direction of travel in 2D vector form (x,y), where the magnitude of the vector is the linear velocity of the asteroid in km/h." << endl;
        cin >> b >> c;
        vx_initial = 0.27*b;
        vy_initial = 0.27*c;
        cout << endl;

        cout << "Specify the location (in km) of the asteroid in 2D vector form." << endl;
        cin >> d >> e;
        x_initial = 1000*d;
        y_initial = 1000*e;
        cout << endl;

        v_initial.x = vx_initial, v_initial.y = vy_initial; // Obtain v vector
        r_initial.x = x_initial, r_initial.y = y_initial; // Obtain r vector
        r = magnitude(r_initial); // Determine initial radius of asteroid r

        cout << "Specify the mass of the asteroid in kg." << endl;
        cin >> m;
        cout << endl;

        // Find angle between these vectors
        double angle;
        angle = acos(dot_product(r_initial, v_initial) / ((magnitude(r_initial))*magnitude(v_initial)));

        // Break velocity vector into parallel and perp components with r vector, as well as angular v and displacement
        double v_para, v_perp;
        v_para = magnitude(v_initial)*cos(angle);
        v_perp = magnitude(v_initial)*sin(angle);
        omega = v_perp/r;
        theta = atan((r_initial.y/r_initial.x));

        // Get time parameters and step size
        double tyears, tseconds, dt;
        cout << "Specify the time period in years over which you would like to run the program, and the number of steps." << endl;
        cin >> tyears;
        cout << endl;
        tseconds = tyears*31536000;
        dt = tseconds/n;

         // Declare arrays and initial conditions
         double t_array[n], u_array[n], v_array[n], dudt_array[n],  dvdt_array[n], theta_array[n], r_array[n];
         theta_array[0] = theta;
         r_array[0] = r;
         v_array[0] = v_para;
         u_array[0] = omega;
         t_array[0] = 0;
         int i;
         const double G = 6.67e-11;
         double f = 0;

         // Solve equations
         for (i=0; i < (n-1); i++) {
             dvdt_array[i] = dvdt(r_array[i], u_array[i], G, M);
             dudt_array[i] = dudt(u_array[i], v_array[i], r_array[i]);
             v_array[i+1] = v_array[i] + dt*dvdt_array[i];
             u_array[i+1] = u_array[i] + dt*dudt_array[i];
             r_array[i+1] = r_array[i] + dt*v_array[i+1];
             theta_array[i+1] = theta_array[i] + dt*u_array[i+1];
             x_array[i] = r_array[i]*cos(theta_array[i]);
             y_array[i] = r_array[i]*sin(theta_array[i]);
             t_array[i+1] = t_array[i] + dt;
         }

         // Output results
         for (i=0; i < n; i++) {
             cout << t_array[i]/31536000 << " , " << theta_array[i] << " , " << r_array[i]/1000 << " , " << endl;
         }

    }

}
