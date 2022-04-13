# Authors
- [Adam Piszczek](https://github.com/Adam1904)
- [Wiktoria Zych](https://github.com/polarrev)

# Rosenbrock function optimization
Rosenbrock function optimization with four different methods (unconstrained optimization). What is noticeable in the figure below, it is a convex function with a hard-to-find minimum.

<p align="center">
  <img width="561" height="421" src="https://github.com/Adam1904/rosenbrock-function-optimization/blob/main/Rosenbrockfunction.png">
</p>

## How to Run

```sh
> main.m
```
## Setup

Before the first run, you need to make sure that you have installed the required dependency.

## Dependiencies 
The application requires [Optimization Toolbox](https://www.mathworks.com/products/optimization.html) to be installed; besides, it uses functions from the standard library.

## About
This script aims to find the minimum at the lowest number of function calls; the algorithm should also have the best computational cost. The main problem with reaching the minimum is overshooting the searched point. Gradient algorithms perform additional iterations around the optimal point. In this case, you should choose the optimization parameters carefully so as not to perform redundant iterations.
