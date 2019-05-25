#pragma once

#include <cmath>
#include <vector>
#include <cstdint>
#include <fstream>
#include <filesystem>

#include "../json/json.hpp"
using json = nlohmann::json;

class HeatEquation {
public:
    /*!
     * Constructor class of heat equation. Requires exact solutions for the entire surface.
     * @param heat_sources - The heat sources function of the hear equation.
     * @param exact_solution - Exact solution of exact equation. Required to find the boundary conditions.
     * @param diffusivity_coefficient - Diffusivity coefficient of the heat equation.
     * @param h_num - The number of spatial splits.
     * @param time_layers_num - The number of time layers.
    */
    HeatEquation(double_t (*heat_sources)(double_t, double_t),
                  double_t (*exact_solution)(double_t, double_t),
                  double_t diffusivity_coefficient,
                  uint32_t h_num,
                  uint32_t time_layers_num);
    /*!
     * The heat equation class constructor. Requires functions of boundary conditions(initial time layer, left and right
     * bounds). Unlike the previous constructor does not require an exact solution on the entire surface.
     * @param heat_sources - The heat sources function of the hear equation.
     * @param initial_time_layer - The initial time layer function. Boundary condition.
     * @param left_bound - The left boundary condition.
     * @param right_bound - The right boundary condition.
     * @param diffusivity_coefficient - Diffusivity coefficient of the heat equation.
     * @param h_num - The number of spatial splits.
     * @param time_layers_num - The number of time layers.
    */
    HeatEquation(double_t (*heat_sources)(double_t, double_t),
                  double_t (*initial_time_layer)(double_t),
                  double_t (*left_bound)(double_t),
                  double_t (*right_bound)(double_t),
                  double_t diffusivity_coefficient,
                  uint32_t h_num,
                  uint32_t time_layers_num);

    /*!
     * The function of finding the error of calculation. The error is total.
     * @param exact_solution - The function of exact solution of PDE.
     * @return - The total error.
     */
    double_t GetError(double_t (*exact_solution)(double_t, double_t)) const noexcept;
    void WriteExactSolution(double_t (*exact_solution)(double_t, double_t), std::filesystem::path const & path = ".") const;
    void WriteExactSolutionJSON(double_t (*exact_solution)(double_t, double_t), std::filesystem::path const & path = ".") const;
    void WriteError(double_t (*exact_solution)(double_t, double_t), std::filesystem::path const & path = ".") const;
    void WriteErrorJSON(double_t (*exact_solution)(double_t, double_t), std::filesystem::path const & path = ".") const;
    void WriteResult(std::filesystem::path const & path = ".") const;
    void WriteResultJSON(std::filesystem::path const & path = ".") const;

public:
    HeatEquation() = delete;
    ~HeatEquation() = default;

private:
    /*!
     * Modified Thomas Algorithm. The diagonal of the matrix composed of the same elements,
     * as well upper and lower diagonal coincide.
     * @param free_part - F from Ax = F.
     * @param result - The result(x) vector.
     */
    void ModifiedThomasAlg(std::vector<double_t> const & free_part, std::vector<double_t> & result) const noexcept;
    std::ofstream GenerateOfstream(std::filesystem::path const & path, std::string const & type) const;
    json GenerateJSON(std::string const & type) const;

private:
    //Time and space boundaries + diffusivity coefficient(a):
    double_t const x_left_bound_      {0.0};
    double_t const x_right_bound_     {1.0};
    double_t const time_left_bound_   {0.0};
    double_t const time_right_bound_  {1.0};

    //Unknown at compile time(define in constructor):
    uint32_t const h_num_;
    uint32_t const time_layers_num_;
    double_t const a_;

    //Get the values during the initialization of
    //h_hum and time_layers_num.
    double_t const h_;
    double_t const tau_;
    double_t const matrix_above_;
    double_t const matrix_main_;

    std::vector<std::vector<double_t>> results_;
};
