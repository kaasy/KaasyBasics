using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Diagnostics;

namespace Utilities
{
    public static class FourierReal
    {
        public static void FFT(this IList<ComplexNumber> list, int precisionDigits)
        {
            int n = list.Count;
            if ((n & (n - 1)) != 0)
            {
                throw new InvalidOperationException("This FFT version is only supported on power of 2 lengths.");
            }
            ComplexNumber[] powers = new ComplexNumber[n];
            int m = n >> 1;
            List<ComplexNumber> roots = new List<ComplexNumber>();
            RealNumber zero = new RealNumber(0L, precisionDigits);
            RealNumber one = new RealNumber(1L, precisionDigits);
            roots.Add(new ComplexNumber(zero, one));
            for (int k = 4; k < n; k *= 2)
            {
                RealNumber lastCos = roots[roots.Count - 1].Real;
                RealNumber cos = ((one + lastCos) >> 1).GetSqrt();
                RealNumber sin = ((one - lastCos) >> 1).GetSqrt();
                roots.Add(new ComplexNumber(cos, sin));
            }
            fill(powers, roots, precisionDigits, m, m);
            for (int i = m; --i > 0;)
            {   //iterate from up to down to initialize correctly.
                powers[i] = powers[i * 2];
            }
            list.binaryReverseSwap();
            list.fft(powers, 0, n);
        }

        public static void IFFT(this IList<ComplexNumber> list, int precisionDigits)
        {
            list.FFT(precisionDigits);
            for (int i = 1, j = list.Count - 1; i < j; i++, j--)
            {
                ComplexNumber auxiliary = list[i];
                list[i] = list[j];
                list[j] = auxiliary;
            }
        }

        private static void fill(ComplexNumber[] roots, List<ComplexNumber> smallRoots, int precisionDigits, int start, int n)
        {
            RealNumber zero = new RealNumber(0L, precisionDigits);
            RealNumber one = new RealNumber(1L, precisionDigits);
            roots[start] = new ComplexNumber(one, zero);
            int index = smallRoots.Count;
            for (int k = 1; k < n; k *= 2)
            {
                index--;
                for (int j = k; --j >= 0;)
                {
                    roots[start + k + j] = roots[start + j] * smallRoots[index];
                }
            }
        }

        private static void binaryReverseSwap(this IList<ComplexNumber> list)
        {
            int n = list.Count;
            for (int i = 1, j = 0; i < n; i++)
            {
                int bit = n;
                do
                {
                    bit >>= 1;
                } while ((j & bit) != 0);
                j ^= n - bit;
                if (i < j)
                {
                    ComplexNumber auxiliary = list[i];
                    list[i] = list[j];
                    list[j] = auxiliary;
                }
            }
        }

        private static void fft(this IList<ComplexNumber> list, ComplexNumber[] powers, int start, int n)
        {
            n >>= 1;
            if (n <= 0)
            {
                return;
            }
            int middle = start + n;
            list.fft(powers, start, n);    //recursive on even indices.
            list.fft(powers, middle, n);   //recursive on odd indices.
            for (int i = n; --i >= 0;)
            {
                ComplexNumber number1 = list[start + i];
                ComplexNumber number2 = list[middle + i] * powers[n + i];
                list[start + i] = number1 + number2;
                list[middle + i] = number1 - number2;
            }
        }

        public static bool UnitTest(int seed)
        {
            Random random = new Random(seed);
            int n = 512;
            ComplexNumber[] a0 = new ComplexNumber[n];
            ComplexNumber[] a1 = new ComplexNumber[n];
            int precisionDigits = 2048 / 64;
            for (int i = n; --i >= 0;)
            {
                a0[i] = a1[i] = new ComplexNumber(
                    new RealNumber(random.NextDouble() * 2 - 1, precisionDigits),
                    new RealNumber(random.NextDouble() * 2 - 1, precisionDigits));
            }
            a0.FFT(precisionDigits);
            a0.IFFT(precisionDigits);
            RealNumber scale = new RealNumber(1.0 / n, 1);
            for (int i = n; --i >= 0;)
            {
                a0[i] *= scale;
            }
            RealNumber error = a0.GetMaximumError(a1);
            RealNumber tolerance = RealNumber.One >> (precisionDigits * 64 - 16);
            bool errorOK = error < tolerance;
            return errorOK;
        }

        public static RealNumber GetMaximumError(this IList<ComplexNumber> list1, IList<ComplexNumber> list2)
        {
            RealNumber maximumErrorSquared = 0;
            for (int i = Math.Max(list1.Count, list2.Count); --i >= 0;)
            {
                maximumErrorSquared = RealNumber.Max(maximumErrorSquared, (list1[i] - list2[i]).Energy);
            }
            return maximumErrorSquared.GetSqrt();
        }
    }
}
