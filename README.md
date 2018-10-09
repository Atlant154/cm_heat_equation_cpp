# Heat equation. C++ implementation.

![Logo](docs/logo.jpg)

## Problem

We need to solve a [partial differential equation](https://en.wikipedia.org/wiki/Partial_differential_equation) of [heat equation](https://en.wikipedia.org/wiki/Heat_equation) using [implicit Euler method](https://en.wikipedia.org/wiki/Backward_Euler_method).  

## Requirements

* CMake.
* G++ compiler.
* Python 3* with installed matplotlib.

## How to start

1. Clone repo: `git clone https://github.com/Atlant154/cm_heat_equation_cpp.git`
2. Move to directory: `cd cm_heat_equation_cpp`
3. Build the project: `cmake -DCMAKE_BUILD_TYPE=Release .. && make`
4. Run the program: `./cm_heat_equation_cpp`

## Visualization

After running the program in the `./result` directory, the `result.txt` file is generated.  
To visualize the data obtained during the execution of the program, do the following command: `cd results && python3 vizualization.py`

![Visualization](docs/vis.png)

Also in the `./result` directory there is a `check.py` script, which is easily configured to get an image of the exact solution, if there is one.