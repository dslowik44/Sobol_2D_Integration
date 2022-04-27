# Sobol_2D_Integration
Estimate a 2-D integral on unit square with Sobol sequence.
Provides functions, `sobol_2D_integral` and `Sobol_2D.integrate`, in C++ and C#, to integrate the provided function on the unit square, [0,1] x [0,1].
The code shows how to do this using the function `circle4Indicator` providing estimates of Pi.

The C# function, `Sobol_2D.integrate`, has third argument bool `prll`, to toggle Parallel processing of the integrand evaluation loop. 
There are Stopwatchs in the code to decide if this is worthwhile. 
It is worthwhile when the time spent evaluating the integrand is relatively long. 
For example, removing `Thread.Sleep(1);` from `circle4Indicator` gives faster times with `prll = false`.
The `Parallel.For` has a `MaxDegreeOfParallelism` options parameter to set max number of threads that can also be adjusted.
