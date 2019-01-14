#pragma once

#include <cmath>
#include <vector>
#include <cstdint>
#include <fstream>

class heat_equation {
public:
    /*!
     * Constructor class of heat equation. Requires exact solutions for the entire surface.
     * @param heat_sources - The heat sources function of the hear equation.
     * @param exact_solution - Exact solution of exact equation. Required to find the boundary conditions.
     * @param diffusivity_coefficient - Diffusivity coefficient of the heat equation.
     * @param h_num - The number of spatial splits.
     * @param time_layers_num - The number of time layers.
    */
    heat_equation(double (*heat_sources)(double, double),
                  double (*exact_solution)(double, double),
                  double diffusivity_coefficient,
                  uint32_t h_num,
                  uint32_t time_layers_num);
    /*!
     * The heat equation class constructor. Requires functions of boundary conditions(initial time layer, left and right
     * bounds). Unlike the previous constructor does not require an exact solution on the entire surface.
     * @param heat_sources - The heat sources function of the hear equation.
     * @param initial_time_layer - The initial time layer function. Boundary condition.
     * @param left_bound - The left boundary condition.
     * @param right_bound - The right boundary condition.
     * @param diffusivity_coefficient - Diffusivity coefficient of the heat eaquation.
     * @param h_num - The number of spatial splits.
     * @param time_layers_num - The number of time layers.
    */
    heat_equation(double (*heat_sources)(double, double),
                  double (*initial_time_layer)(double),
                  double (*left_bound)(double),
                  double (*right_bound)(double),
                  double diffusivity_coefficient,
                  uint32_t h_num,
                  uint32_t time_layers_num);

    /*!
     * The function of finding the error of calculation. The error is total.
     * @param exact_solution - The function of exact solution of PDE.
     * @return - The total error.
     */
    double get_error(double (*exact_solution)(double, double)) const;
    void write_error_plot(double (*exact_solution)(double, double), std::string const & path = ".") const;
    void write_result(std::string const & path = ".") const;

public:
    heat_equation() = delete;
    ~heat_equation() = default;

private:
    /*!
     * Modified Thomas Algorithm. The diagonal of the matrix composed of the same elements,
     * as well upper and lower diagonal coincide.
     * @param free_part - F from Ax = F.
     * @param result - The result(x) vector.
     */
    void modified_thomas_alg(std::vector<double> const & free_part, std::vector<double> & result);

private:
    //Time and space boundaries + diffusivity coefficient(a):
    double const x_left_bound_ = 0.0;
    double const x_right_bound_ = 1.0;
    double const time_left_bound_ = 0.0;
    double const time_right_bound_ = 1.0;

    //Unknown at compile time(define in constructor):
    uint32_t const h_num_;
    uint32_t const time_layers_num_;
    double const a_;

    //Get the values during the initialization of
    //h_hum and time_layers_num.
    double const h_;
    double const tau_;
    double const matrix_above_;
    double const matrix_main_;

    std::vector<std::vector<double>> results_;
};
