using System;
using System.Collections.Generic;
using System.Linq;
using System.Numerics;
using System.Text;
using System.Threading.Tasks;

namespace Utilities
{
    public interface GeneralSpline
    {
        Complex Evaluate(double tx, int index);
    }

    public class PeriodicCubicSpline : GeneralSpline
    {
        private List<double> x;
        private List<Complex> y;
        private Complex[] a, b, c, d;
        private int n;

        /// <summary>
        /// list of initial points is (x0,y0) (x1,y1) ... (xn, yn) - total n + 1 points
        /// http://www.math.ou.edu/~npetrov/project-5093-s11.pdf
        /// </summary>
        /// <param name="X"></param>
        /// <param name="Y"></param>

        public PeriodicCubicSpline(IEnumerable<double> X, IEnumerable<Complex> Y)
        {
            this.x = new List<double>(X);
            this.y = new List<Complex>(Y);

            this.n = this.x.Count;
            if (n <= 1)
            {
                return;
            }

            this.x.Add(x[x.Count - 1] * 2 - x[x.Count - 2]);

            double[] h = new double[n];

            a = new Complex[n + 1];
            b = new Complex[n];
            c = new Complex[n + 1];
            d = new Complex[n];

            double[] diag = new double[n];
            double[] coln = new double[n];

            //step 1. Init
            for (int i = 0; i < n; i++)
            {
                a[i] = y[i];
                h[i] = x[i + 1] - x[i];
            }
            a[n] = a[0];

            Complex[] alpha = new Complex[n];
            alpha[0] = ((a[1] - a[0]) / h[0] - (a[0] - a[n - 1]) / h[n - 1]) * 3;
            diag[0] = (h[n - 1] + h[0]) * 2;
            coln[0] = h[n - 1];
            coln[n - 2] = h[n - 2];
            for (int i = 1; i < n; i++)
            {
                alpha[i] = ((a[i + 1] - a[i]) / h[i] - (a[i] - a[i - 1]) / h[i - 1]) * 3;
                diag[i] = (x[i + 1] - x[i - 1]) * 2;    // = (h[i-1]+h[i])*2
            }

            //step 2. Gaussian elimination
            for (int i = 0; i + 2 < n; i++)
            {
                //operation: row[i + 1] -= row[i] * (h[i] / diag[i])
                double factor = h[i] / diag[i];
                //lower[i] = 0;
                diag[i + 1] -= h[i] * factor;
                coln[i + 1] -= factor * coln[i];
                alpha[i + 1] -= factor * alpha[i];
            }
            //step 3. Back substitution
            for (int i = n - 2; --i >= 0;)
            {
                //operation: row[i] -= row[i + 1] * (h[i] / diag[i + 1])
                double factor = h[i] / diag[i + 1];
                coln[i] -= factor * coln[i + 1];
                alpha[i] -= factor * alpha[i + 1];
            }

            {
                double factor = h[n - 1] / diag[0];
                //row[n - 1] -= row[0] * factor 
                diag[n - 1] -= factor * coln[0];
                alpha[n - 1] -= factor * alpha[0];
            }
            {
                double factor = h[n - 2] / diag[n - 2];
                //row[n - 1] -= row[n - 2] * factor 
                diag[n - 1] -= factor * coln[n - 2];
                alpha[n - 1] -= factor * alpha[n - 2];
            }

            {
                c[n - 1] = alpha[n - 1] / diag[n - 1];
                for (int i = n - 1; --i >= 0;)
                {
                    c[i] = (alpha[i] - c[n - 1] * coln[i]) / diag[i];
                }
            }
            c[n] = c[0];
            for (int j = n; --j >= 0;)
            {
                b[j] = (a[j + 1] - a[j]) / h[j] - h[j] * (c[j + 1] + c[j] * 2) / 3;
                d[j] = (c[j + 1] - c[j]) / (h[j] * 3);
            }

            /*
            Complex[] l = new Complex[n + 1];
            Complex[] miu = new Complex[n + 1];
            Complex[] z = new Complex[n + 1];
            l[0] = Complex.One; miu[0] = z[0] = Complex.Zero;
            for (int i = 1; i < n; i++)
            {
                l[i] = (x[i + 1] - x[i - 1]) * 2 - h[i - 1] * miu[i - 1];
                miu[i] = h[i] / l[i];
                z[i] = (alpha[i] - h[i - 1] * z[i - 1]) / l[i];
            }
            l[n] = Complex.One; c[n] = z[n] = Complex.Zero;

            for (int j = n; --j >= 0; )
            {
                c[j] = z[j] - miu[j] * c[j + 1];
                b[j] = (a[j + 1] - a[j]) / h[j] - h[j] * (c[j + 1] + c[j] * 2) / 3;
                d[j] = (c[j + 1] - c[j]) / (h[j] * 3);
            } */
        }

        public Complex Evaluate(double tx, int index)
        {
            if (n <= 0)
            {
                return Complex.Zero;
            }
            if (n == 1)
            {
                return y[0];
            }
            index = Math.Min(Math.Max(index, 0), n - 1);
            tx -= x[index];
            return a[index] + tx * (b[index] + tx * (c[index] + tx * d[index]));
            //return (b[index] + (2 * c[index] + 3 * d[index] * tx) * tx) * 0.10 + 0.5; //first derivative
            //return 2 * (c[index] + 3 * d[index] * tx) * 0.01 + 0.5;   //second derivative
        }
    }
}
