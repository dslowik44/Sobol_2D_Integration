using System.Diagnostics; // for Stopwatch

using Sobol2dInteg;

// There are a few Stopwatches timing blocks of code within Sobol_2D.integrate.
// When the time spent in f(x,y) is relatively large,
// total time spent by Sobol_2D.integrate is:
// Outer time = [time spent evaluating integrand, func(x,y),] times [the number of points
// used in the integration,] divided by [the number of concurrent threads.]
//
// In the regime sleep_ms >= 1 the time spent in func(x, y) dominates.
// 
class Program
{
    static double circle4Indicator(double x, double y)
    {  // Integrand is 4 * Indicator function, == 4 inside the unit circle
       // with (x, y) on [0,1]x[0,1].
        int sleep_ms = 1;  // Used to illustrate effects of a slow integrand function.
        Thread.Sleep(sleep_ms); 
        return x * x + y * y < 1.0 ? 4.0 : 0.0;
    }

    Func_delegate ff = new Func_delegate( circle4Indicator );

    static void Main(string[] args)
    {
        uint[] a = new uint[] { 100, 1000, 10000 }; // Number of Sobol points to use.

        Console.WriteLine("MC simulate fraction of (x,y) 2-D Sobol sequence points falling " +
                          "inside unit circle:" );
        
        var watch = Stopwatch.StartNew();
        double[] Pi_calc = Sobol_2D.integrate(a, circle4Indicator, true);
        watch.Stop();
	    Console.WriteLine($"Outer: {watch.ElapsedMilliseconds}ms");
        
        for (int i=0; i < Pi_calc.Length; ++i) {
            Console.WriteLine("N: " + $"{a[i]}".PadLeft(10) + $", Calculated Pi: {Pi_calc[i]:0.000000}, " +
                              $"Diff: {Pi_calc[i] - Math.PI:0.000000}");
        }
    }
}
