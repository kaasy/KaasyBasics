using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Utilities
{
    //http://code.activestate.com/recipes/393090/
    public static class AccurateSummation
    {
        private static double unoptimizedSubtract(double x, double y)
        {
            return x - y;
        }

        private static List<double> getPartials(this IEnumerable<double> listing)
        {
            List<double> partials = new List<double>();               // sorted, non-overlapping partial sums
            foreach (double lx in listing)
            {
                int i = 0;
                double x = lx;
                for (int k = 0; k < partials.Count; k++)
                {
                    double y = partials[k];
                    if (Math.Abs(x) < Math.Abs(y))
                    {
                        double auxiliary = x;
                        x = y;
                        y = auxiliary;
                    }
                    double hi = x + y;
                    double lo = y - unoptimizedSubtract(hi, x);
                    if (lo != 0)
                    {
                        partials[i++] = lo;
                    }
                    x = hi;
                }
                partials.RemoveRange(i, partials.Count - i);
                partials.Add(x);
            }
            double result = partials.Sum();
            return partials;
        }

        private class DoubleAbsComparerDescending : IComparer<double>
        {
            public int Compare(double x, double y)
            {
                return Math.Abs(y).CompareTo(Math.Abs(x));
            }
        }

        private static double getAccurateSumSlow(this IEnumerable<double> numbers)
        {
            List<double> list = new List<double>(numbers);
            var comparer = new DoubleAbsComparerDescending();
            bool modified = true;
            while (modified)
            {
                modified = false;
                list.Sort(comparer);
                for (int i = 0; i + 1 < list.Count; i += 2)
                {
                    double sum = list[i] + list[i + 1];
                    double rest = list[i + 1] - unoptimizedSubtract(sum, list[i]);
                    modified = modified || list[i] != sum || list[i + 1] != rest;
                    list[i] = sum;
                    list[i + 1] = rest;
                }
                list.RemoveAll(x => x == 0);
            }
            return list.Count == 0 ? 0 : list.Sum();
        }

        public static double GetAccurateSum(this IEnumerable<double> listing)
        {
            List<double> partials = listing.getPartials();
            double result = partials.Sum();
            return result;
        }

        public static bool UnitTest()
        {
            List<double> list = new List<double>();
            Random random = new Random(1001);
            int n = 1000 * 10;
            for (int i = n; --i >= 0;)
            {
                double number = (random.Next(2) == 0 ? -1 : 1) * Math.Pow(10.0, random.Next(301)) * random.NextDouble();
                list.Add(number);
            }
            var partials = list.getPartials();
            list.AddRange(partials.Select(p => -p));
            list.Add(1);
            double sum1 = list.GetAccurateSum();
            double sum2 = list.Sum();
            double sum3 = getAccurateSumSlow(list);
            return sum1 == sum3;
        }
    }
}
