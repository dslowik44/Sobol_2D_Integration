using System.Diagnostics; // for Stopwatch
using System.Collections.Concurrent;

namespace Sobol2dInteg
{
    /* Code inspired by sobol.cc available at:
       https://web.maths.unsw.edu.au/~fkuo/sobol/
    */

    public delegate double Func_delegate(double x, double y);

    public class Sobol_2D
    {

        public static double[] integrate(uint[] a, Func_delegate func, bool prll)
        {
            /* Estimate the 2 dimensional integral over the unit square of function provided as second
            argument. Returns this estimate as the Expectation (sample mean) for each of the number of 
            points specified in a[]. Assumes a[] is strictly increasing.
            Last argument specifies whether to use Parallel.For over func(x,y) integration loop 
            - Set to true when time spent is dominated by time to call the integrand func(x,y).
            */

            double TWO_m32 = Math.Pow(2.0, -32);

            // est_means will hold estimated integral values corresponding to each element of a.
            double[] est_means = new double[a.Length];
            uint N = a[a.Length - 1];  // Largest number of Sobol points that will be needed.

            // L = max number of bits needed for this many points:
            uint L = (uint)(Math.Ceiling(Math.Log((double)N) / Math.Log(2.0)));

            // C[i] = index from the right of the first zero bit of i < N:
            uint[] C = new uint[N];
            C[0] = 1;
            for (uint i = 1; i < N; ++i)
            {
                C[i] = 1;
                uint value = i;
                while ((value & 1) != 0)
                {
                    value >>= 1;
                    C[i]++;
                }
            }

            // Compute direction numbers V_x[1] to V_x[L], scaled by pow(2,32)
            uint[] V_x = new uint[L + 1];
            for (int i = 1; i <= L; ++i)
                V_x[i] = (uint)(1u << (32 - i));


            // Compute direction numbers V_y[1] to V_y[L]:
            uint[] V_y = new uint[L + 1];
            V_y[1] = (1u << 31);
            for (uint i = 2; i <= L; i++)
            {
                V_y[i] = V_y[i - 1] ^ (V_y[i - 1] >> 1);
            }

            if (prll) 
            {
                // Use V_x, V_y and C to get X[] and Y[]:
                var watch = Stopwatch.StartNew();
                uint[] X = new uint[N];
                uint[] Y = new uint[N];
                double[] values = new double[N];
                values[0] = func(0, 0);
                X[0] = 0;
                Y[0] = 0;
                for (uint i = 1; i < N; i++)
                {
                    X[i] = X[i - 1] ^ V_x[C[i - 1]];  // Note*: C[i] never 0.
                    Y[i] = Y[i - 1] ^ V_y[C[i - 1]];
                }
                watch.Stop();
                Console.WriteLine($"Set X, Y: {watch.ElapsedMilliseconds}ms"); 

                var watch2 = Stopwatch.StartNew();
                var threadIDs = new ConcurrentDictionary<int, int>();
                Parallel.For(1, (long)N, new ParallelOptions { MaxDegreeOfParallelism = 20 }, i =>
                {
                    threadIDs.AddOrUpdate(Thread.CurrentThread.ManagedThreadId, 1, (key, oldValue) => oldValue + 1);
                    double x = X[i] * TWO_m32;
                    double y = Y[i] * TWO_m32;
                    values[i] = func(x, y);  
                });


                double cumsum = 0.0;
                int a_idx = 0;
                for (int i = 0; i < N; ++i)
                {
                    cumsum += values[i];
                    if (i == a[a_idx] - 1)
                    {
                        est_means[a_idx] = cumsum / a[a_idx];
                        ++a_idx;
                    }
                }
                watch2.Stop();
                Console.WriteLine($"parallel integration: {watch2.ElapsedMilliseconds}ms");

                // Print info on threads used:
                Console.WriteLine($"Used {threadIDs.Keys.Count()} threads:");
                foreach (var key in threadIDs.Keys) {
                    Console.WriteLine($"Thread ID: {key} called func {threadIDs[key]} times.");
                }

            }
            else // not parallel:
            {
                // Use V_x, V_y and C to get (x, y) Sobol points.
                // Accumulate f(x, y) into accum for mean estimate:
                var watch2 = Stopwatch.StartNew();
                uint X = 0;
                uint Y = 0;
                int a_idx = 0;
                double cumsum = func(0, 0);
                for (uint i = 1; i < N; i++) // sequential scan to largest N == a[a.Length]:
                {
                    X = X ^ V_x[C[i - 1]];  // Note*: C[i] never 0.
                    Y = Y ^ V_y[C[i - 1]];
                    double x = (double)X * TWO_m32;
                    double y = (double)Y * TWO_m32;
                    cumsum += func(x, y);
                    if (i == a[a_idx] - 1)
                    {
                        est_means[a_idx] = cumsum / a[a_idx];
                        ++a_idx;
                    }
                }
                watch2.Stop();
                Console.WriteLine($"{watch2.ElapsedMilliseconds}ms");
            }

            return est_means;
        }
    }

}


/* // Failed parallel speedup: To much contention for shared accumulator among threads.

               for (int a_idx = 0; a_idx < a.Length; ++a_idx)
                {
                    var watch2 = Stopwatch.StartNew();
                    // Calculate accum for this a[a_idx]:
                    uint startidx = a_idx == 0 ? 1 : a[a_idx - 1];
                    if (a[a_idx] > 2000) // Parallel.For for large enough set of points:
                    {
                        // Accumulate f(x, y) into accum for mean estimate:
                        Parallel.For(startidx, (long)a[a_idx], i =>
                        {
                            double x = X[i] * TWO_m32;
                            double y = Y[i] * TWO_m32; 
                            Thread.Sleep(1); 
                            double val = x * x + y * y < 1.0 ? 4.0 : 0.0; //func(x, y);  5% speed increase.
                            // Atomically increment accum:
                            double new_current_accum = accum;
                            while (true)
                            { // This loops on avg 9 times => Too much contention for accum among the threads.
                                double current_accum = new_current_accum;
                                double new_accum = current_accum + val;
                                new_current_accum = Interlocked.CompareExchange(ref accum, new_accum, current_accum);
                                if (new_current_accum == current_accum)
                                    break;
                            }
                        });
                    }
                    else // Non-parallel calulation of accum for this a[a_idx]:
                    {
                        for (uint i = startidx; i < a[a_idx]; i++)
                        {
                            double x = X[i] * TWO_m32; 
                            double y = Y[i] * TWO_m32;
                            //double val = func(x, y);
                            accum += x * x + y * y < 1.0 ? 4.0 : 0.0; //val;
                        }
                    }
                    est_means[a_idx] = accum / a[a_idx];
                    watch2.Stop();
                    Console.WriteLine($"a_idx: {a_idx}: {watch2.ElapsedMilliseconds}ms");
                }
            }
*/
