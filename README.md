# Уравнение теплопроводности. Heat equation. C++.

![Logo](docs/logo.jpg)

## Build status

[![Build Status](https://travis-ci.org/Atlant154/cm_heat_equation_cpp.svg?branch=master)](https://travis-ci.org/Atlant154/cm_heat_equation_cpp)

## Problem

Нужно было решить [PDE(уравнение в частных производных)](https://en.wikipedia.org/wiki/Partial_differential_equation) [heat equation(уравнения теплопроводности)](https://en.wikipedia.org/wiki/Heat_equation) исользуя [implicit Euler method(неявный метод Эйлера)](https://en.wikipedia.org/wiki/Backward_Euler_method).  
Для решения была написана библиотека `heateq`, а также приложение, которое позволяет её использовать. 

## Requirements

* [CMake](https://cmake.org/) v3.10 performance guaranteed.
* [G++ compiler](https://gcc.gnu.org/) v8.2 performance guaranteed.
* [Python](https://www.python.org/) 3* with installed [matplotlib](https://matplotlib.org/) for visualization.

## How to start

1. Клонировать репозиторий с сабмодулями: `git clone --recurse-submodules https://github.com/Atlant154/cm_heat_equation_cpp.git`
2. Переместится в директорию приложения: `cd cm_heat_equation_cpp`
3. Создать директорию для сборки: `mkdir build && cd build`
3. Собрать проект: `cmake -DCMAKE_BUILD_TYPE=Release .. && make`
4. Запустить программу: `./cm_heat_equation_cpp --tau 1500 --splits 1500 --write`

Программа имеет следующие флаги/опции запуска:

```
  -h,--help                   Print this help message and exit  
  -t,--tau UINT REQUIRED      The number of time layers  
  -s,--splits UINT REQUIRED   The number of of spatial splits  
  -w,--write                  Write results to file  
  -o,--output-path TEXT       The output directory path, i.e. ~/projects/visualization. Default: ./
```

## Visualization

После выполнения программы с `write`-флагом, в указанной директории появится файлы требуемые для визуализации. Построить графики приближенного решения, ошибки и точного решения можно запустив
скрипт `heateq/visualization/visualization.py` в директории с сохранёнными файлами: `python3 visualization.py` 

![Visualization](docs/vis.png)

## Solution and hacks

### Thomas Algorithm

Для решения трёхдиагональной матрицы, получаемой в ходе решения, обычно используется
[Thomas Algorithm(алгоритм прогонки)](https://en.wikipedia.org/wiki/Tridiagonal_matrix_algorithm).
При решении конкретно этой задачи мы получаем трёхдиагональную матрицу специального вида: 
элементы верхней и нижней диагонали равны одному значению, потому для наиболее эффективного использования
использовался видоизменённый алгоритм прогонки:  
```C++
void HeatEquation::ModifiedThomasAlg(std::vector<double_t> const & free_part, std::vector<double_t> & result) const noexcept {
    assert(free_part.size() == result.size());
    std::size_t const n = free_part.size();
    std::vector<double_t> alpha(n - 1), beta(n - 1);

    double_t common_factor;

    alpha[n - 2] = -m_matrixAboveDiagonalElement / m_matrixMainDiagonalElement;
    beta[n - 2] = free_part.back() / m_matrixMainDiagonalElement;

    for (auto iter {n - 2}; iter > 0; --iter) {
        common_factor   = 1. / (m_matrixMainDiagonalElement + m_matrixAboveDiagonalElement * alpha[iter]);
        alpha[iter - 1] = -m_matrixAboveDiagonalElement * common_factor;
        beta[iter - 1]  = (free_part[iter] - beta[iter] * m_matrixAboveDiagonalElement) * common_factor;
    }

    result[0] = (free_part[0] - m_matrixAboveDiagonalElement * beta[0]) / (m_matrixMainDiagonalElement + m_matrixAboveDiagonalElement * alpha[0]);

    for (std::size_t iter{1}; iter < n; ++iter)
        result[iter] = alpha[iter - 1] * result[iter - 1] + beta[iter - 1];
}
```
Такое изменение даёт существенный прирост производительности, но использовать этот алгоритм в других задачах не представляется возможным.

### Finding solution

Для хранения временных результатов и решения используется `std::vector`. Я постарался максимально оптимизировать работу с данными, потому просадка должна быть не более нескольких процентов в сравнении
с сырыми указателями.

## Benchmarks

Тестовые результаты получены при помощи google-benchmark при сотне временных слоёв.

|      Splits:      |     128    |     512    | 1'024      | 8'192      | 32'768    | 131'072   |
|:-----------------:|:----------:|:----------:|------------|------------|-----------|-----------|
| Cpp(unoptimized): | 361.917 μs | 1439.70 μs | 2897.31 μs | 23783.5 μs | 94867  μs | 386435 μs |

Процессор: I7-6700.

## Calculation error

Скорость аппроксимации: ![error](docs/error.png). Тесты скорости можно найти в `heateq/tests`, для проведения тестирования был использован google-test.
