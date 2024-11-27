// This is a toy code to explore how to link BLTC to FARRSIGHT
// We want to use Cmakelist for compiling both cpu(multi-threads) and gpu verisions
// This file inlcudes FieldStructure.cpp and FieldStructure.hpp
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <iomanip> 
#include <chrono>
#include <vector>
#include <iomanip>
#include <cstring>
#include <sys/times.h>
#include <cfloat>
#include <cassert>
using namespace std;
using namespace std::chrono;

#include "FieldStructure.hpp"

int main(int argc, char** argv) {
    // Parse arguments 
    double Lx = std::stod(argv[1]);
    double greens_epsilon = std::stod(argv[2]);
    double mac = std::stod(argv[3]);
    double degree = std::stod(argv[4]);
    double max_source = std::stod(argv[5]);
    double max_target = std::stod(argv[6]);
    int numpars_s = std::stoi(argv[7]);
    int verbosity = 0;
    int use_treecode = std::stoi(argv[8]);
    int Nstep = std::stoi(argv[9]);
    double dt = std::stoi(argv[10]);

    cout << "L = " << Lx << endl;
    cout << "greens_epsilon = " << greens_epsilon << endl;
    cout << "mac = " << mac << endl;
    cout << "degree = " << degree << endl;
    cout << "numpars_s = " << numpars_s << endl;

    double pi = 3.14159265358979323846;
    ElectricField* calculate_e;
    if (use_treecode > 0) {
        calculate_e = new E_MQ_Treecode(Lx, greens_epsilon, mac, degree, max_source, max_target, verbosity);
        cout << "using tree code" << endl;
    } else {
        calculate_e = new E_MQ_DirectSum(Lx, greens_epsilon);
        cout << "using direct sum" << endl;
    }

    std::vector<double> particles_x(numpars_s);
    std::vector<double> q_ws(numpars_s);
    double dx = Lx / numpars_s;
    for (int i = 0; i < numpars_s; i++) {
        particles_x[i] = (i + 1) * dx + 0.5 * cos(2 * pi * (i + 1) * dx);
        q_ws[i] = 1.0;
    }

    std::vector<double> e_field(numpars_s);
    auto start = high_resolution_clock::now();
    // Call the calculate_e operator
    for (int i = 0; i < Nstep; i++) {
        cout << "step: " << i << endl;
        (*calculate_e)(e_field.data(), particles_x.data(), e_field.size(), particles_x.data(), q_ws.data(), particles_x.size());
        for (int j = 0; j < particles_x.size(); j++) {
            particles_x[j] = particles_x[j] + dt * e_field[j];
        }
        // for (int j = 0; j < 5; j++) {
        //     cout << "j = " << j << ",particle_x = " << particles_x[j] << endl;
        //     cout << "e_field = " << e_field[j] << endl;
        // }
    }
   
    auto end = high_resolution_clock::now();
    auto duration_us = duration_cast<milliseconds>(end - start); // microseconds
    std::cout << "Time taken: " << duration_us.count() << " milliseconds" << std::endl;


    for (int i = 0; i < numpars_s; i++) {
        particles_x[i] = (i + 1) * dx + 0.5 * cos(2 * pi * (i + 1) * dx);
    }


// Get a referenced e_field using direct sum 
    cout << "For the reference direct sum" << endl;
    std::vector<double> e_field_direct_reference(numpars_s);
    ElectricField* e_diect = new E_MQ_DirectSum(Lx, greens_epsilon);
    auto start1 = high_resolution_clock::now();
    for (int i = 0; i < Nstep; i++) {
        cout << "step: " << i << endl;
        (*e_diect)(e_field_direct_reference.data(), particles_x.data(), e_field_direct_reference.size(), particles_x.data(), q_ws.data(), particles_x.size());
        for (int j = 0; j < particles_x.size(); j++) {
            particles_x[j] = particles_x[j] + dt * e_field_direct_reference[j];
        }
        // for (int j = 0; j < 5; j++) {
        //     cout << "j = " << j << ",particle_x = " << particles_x[j] << endl;
        //     cout << "e_field = " << e_field_direct_reference[j] << endl;
        // }
    }
    auto end1 = high_resolution_clock::now();
    auto duration_us1 = duration_cast<milliseconds>(end1 - start1); // microseconds
    std::cout << "Time taken for the reference direct sum: " << duration_us1.count() << " milliseconds" << std::endl;


// get error : 
    double err2_ex = 0.0;
    double sum_d_ex = 0.0;
    double sum_n_ex = 0.0;
    for (int i = 0; i < numpars_s; i++) {
        sum_n_ex += (e_field_direct_reference[i] - e_field[i]) * (e_field_direct_reference[i] - e_field[i]);
        sum_d_ex += (e_field_direct_reference[i] * e_field_direct_reference[i]);
    }
    err2_ex = sqrt(sum_n_ex / sum_d_ex);

    cout << setprecision(25) << "Error E_2_ex is " << err2_ex << endl;
    cout << "done with this run" << endl;
    cout << endl;
    return 0;
}
