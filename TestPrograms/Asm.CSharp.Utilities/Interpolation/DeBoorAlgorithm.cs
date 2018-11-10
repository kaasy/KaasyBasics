using System;
using System.Collections.Generic;
using System.Linq;
using System.Numerics;
using System.Text;

namespace Utilities
{
    //DeBoor Evaluated Spline
    public class DeBoorAlgorithm : GeneralSpline
    {
        private double[] knots;
        private Complex[] controlPoints;
        private int degree;

        private Complex controlPoint(int index)
        {
            //index += 1;   //index addition works from 1 to this.degree
            index += this.degree;
            index = Math.Min(Math.Max(index, 0), this.controlPoints.Length - 1);
            return this.controlPoints[index];
        }

        private double knot(int index)
        {
            index = Math.Min(Math.Max(index, 0), this.knots.Length - 1);
            return this.knots[index];
        }

        //http://www2.math.ou.edu/~npetrov/project-5093-s11.pdf
        public DeBoorAlgorithm(IEnumerable<double> knots, IEnumerable<Complex> controlPoints, int degree)
        {
            this.knots = knots.ToArray();
            this.controlPoints = controlPoints.ToArray();
            this.degree = degree;
        }

        public Complex Evaluate(double tx, int index)
        {
            Complex[] dk = new Complex[degree + 1];
            int shift = index - degree;
            for (int i = degree + 1; --i >= 0;)
            {   //"minimum index: -degree"; "maximum index: index"
                dk[i] = controlPoint(i + shift);
            }
            for (int k = 1; k <= degree; k++)
            {
                for (int i = degree; i >= k; i--)
                {
                    double alpha = (tx - knot(i + shift)) / (knot(i + 1 + index - k) - knot(i + shift));
                    alpha = double.IsInfinity(alpha) || double.IsNaN(alpha) ? 1.0 : alpha;
                    dk[i] = dk[i - 1] + alpha * (dk[i] - dk[i - 1]);
                }
            }
            return dk[degree];
        }
    }

    public static class ComplexFunctionExtensions
    {
        public static double Distance(this Complex p1, Complex p2)
        {
            return (p1 - p2).Magnitude;
        }

        public static double Square(this double x)
        {
            return x * x;
        }
    }
}
