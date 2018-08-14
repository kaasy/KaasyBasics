using System;
using System.Collections.Generic;
using System.Linq;
using System.Numerics;
using System.Runtime.InteropServices;
using System.Text;
using System.Threading.Tasks;

namespace Utilities
{
    public static class FourierTransform235
    {
        public static void FullForwardFFT(this ComplexFloat[] data, float[] temporaryArray = null)
        {
            int N = data.Length;
            Stack<int> factors = new Stack<int>(GetSlowFactorization(N).OrderBy(x => x));
            ComplexFloat[] roots = null;
            ComplexFloat[] temp = null;
            ComplexFloat[] simpleRoots = null;
            ReversedPowersFFTSwapAtEnd(data, 0, N, factors, 0, 1, ref roots, ref temp, ref simpleRoots);

            if (temporaryArray == null || temporaryArray.Length < N * 2)
            {
                Array.Resize(ref temporaryArray, N * 2);
            }
            for (int i = N; --i >= 0;)
            {
                var number = data[i];
                temporaryArray[i * 2 + 0] = number.Real;
                temporaryArray[i * 2 + 1] = number.Imaginary;
            }
            ReversedBaseIterate(factors, (i, reversed) => data[reversed] = new ComplexFloat(temporaryArray[i * 2 + 0], temporaryArray[i * 2 + 1]));
        }

        public static void FullForwardIFFT(this ComplexFloat[] data, float[] temporaryArray = null)
        {
            FullForwardFFT(data, temporaryArray);
            for (int i = 1, j = data.Length - 1; i < j; i++, j--)
            {
                var auxiliary = data[i];
                data[i] = data[j];
                data[j] = auxiliary;
            }
        }

        public static List<int> GetSlowFactorization(int N)
        {
            var result = new List<int>();
            if (N <= 3)
            {
                result.Add(N);
                return result;
            }
            while (N % 2 == 0)
            {
                result.Add(2);
                N /= 2;
            }
            while (N % 3 == 0)
            {
                result.Add(3);
                N /= 3;
            }

            int p = 5;
            int add = 2;
            while (N / p >= p)
            {
                while (N % p == 0)
                {
                    result.Add(p);
                    N /= p;
                }
                p += add;
                add ^= 6;
            }
            if (N > 1)
            {
                result.Add(N);
            }
            return result;
        }

        private static void ReversedPowersFFTSwapAtEnd(ComplexFloat[] data, int index, int N, Stack<int> factors, int wpower, int npower,
            ref ComplexFloat[] roots, ref ComplexFloat[] temp, ref ComplexFloat[] simpleRoots)
        {
            if (N <= 1 || factors.Count <= 0)
            {
                return;
            }

            int factor = factors.Pop();
            int N1 = N / factor;
            int npf = npower * factor;
            DirectFFTWithStartTwiddle(data, index, N1, N1, wpower, npf, factor, ref roots, ref temp, ref simpleRoots);
            for (int i = factor; --i >= 0;)
            {
                ReversedPowersFFTSwapAtEnd(data, index + i * N1, N1, factors, wpower + npower * i, npf, ref roots, ref temp, ref simpleRoots);
            }
            factors.Push(factor);
        }

        private const float M_SQRT_3_4 = 0.86602540378443864676372317075294f;   // == sqrt(3.0/4) = SIN(PI/3);
        private const float M_ROOT5_C1 = 0.25f;
        private const float M_ROOT5_C2 = 0.55901699437494742410229341718282f;   // M_SQRT5 / 4;
        private const float M_ROOT5_C3 = 0.58778525229247312916870595463907f;   // = sqrt((5 - M_SQRT5) / 8);   // = Im[w^2] = SIN(PI*2/10)
        private const float M_ROOT5_C4 = 0.95105651629515357211643933337938f;   //sqrt((5 + M_SQRT5) / 8);      // = Im[w] = SIN(2*PI/5)

        private static void DirectFFTWithStartTwiddle(ComplexFloat[] data, int index, int multiplier, int length, int wpower, int npower, int prime,
            ref ComplexFloat[] roots, ref ComplexFloat[] temp, ref ComplexFloat[] simpleRoots)
        {
            if (roots == null || roots.Length < prime)
            {
                Array.Resize(ref roots, prime);
            }
            Complex expectedRoot = Complex.FromPolarCoordinates(1.0, (Math.PI * 2) * wpower / npower);
            Complex power = Complex.One;
            for (int i = 1; i < prime; i++)
            {
                power *= expectedRoot;
                roots[i] = power.ToComplexFloat();
            }

            for (; --length >= 0; index++)
            {
                switch (prime)
                {
                    case 2:
                        {
                            ComplexFloat v0 = data[index];
                            ComplexFloat v1 = data[index + multiplier] * roots[1];
                            data[index] = v0 + v1;
                            data[index + multiplier] = v0 - v1;
                        }
                        break;
                    case 3:
                        //w^2 = -1-w
                        {
                            ComplexFloat
                                z0 = data[index],
                                z1 = data[index + multiplier] * roots[1],
                                z2 = data[index + 2 * multiplier] * roots[2];
                            ComplexFloat
                                t1 = z1 + z2,
                                t2 = z0 - t1 * 0.5f,
                                t3 = (z1 - z2).MulI() * M_SQRT_3_4;
                            data[index] = z0 + t1;
                            data[index + multiplier] = t2 + t3;
                            data[index + 2 * multiplier] = t2 - t3;
                        }
                        break;
                    case 5:
                        //http://www2.itap.physik.uni-stuttgart.de/lehre/vorlesungen/SS08/simmeth/fftalgorithms.pdf
                        {
                            ComplexFloat
                                z0 = data[index],
                                z1 = data[index + multiplier] * roots[1],
                                z2 = data[index + 2 * multiplier] * roots[2],
                                z3 = data[index + 3 * multiplier] * roots[3],
                                z4 = data[index + 4 * multiplier] * roots[4];
                            ComplexFloat
                                t1 = z1 + z4,
                                t2 = z2 + z3,
                                t3 = z1 - z4,
                                t4 = z2 - z3,
                                t5 = t1 + t2,
                                t6 = (t1 - t2) * M_ROOT5_C2,
                                t7 = z0 - t5 * M_ROOT5_C1,
                                t8 = t7 + t6,
                                t9 = t7 - t6,
                                t10 = (t3 * M_ROOT5_C4 + t4 * M_ROOT5_C3).MulI(),
                                t11 = (t3 * M_ROOT5_C3 - t4 * M_ROOT5_C4).MulI();
                            data[index] = z0 + t5;
                            data[index + multiplier] = t8 + t10;
                            data[index + 2 * multiplier] = t9 + t11;
                            data[index + 3 * multiplier] = t9 - t11;
                            data[index + 4 * multiplier] = t8 - t10;
                        }
                        break;
                    default:
                        if (simpleRoots == null || simpleRoots.Length != prime)
                        {
                            Array.Resize(ref simpleRoots, prime);
                            for (int i = prime; --i > 0;)
                            {
                                simpleRoots[i] = ComplexFloat.Root((Math.PI * 2) * i / prime);
                            }
                        }
                        if (temp == null || temp.Length < prime)
                        {
                            Array.Resize(ref temp, prime);
                        }
                        ComplexFloat sum = temp[0] = data[index];
                        for (int i = prime; --i > 0;)
                        {
                            temp[i] = data[index + i * multiplier] * roots[i];
                            sum += temp[i];
                        }
                        data[index] = sum;
                        for (int i = prime; --i > 0;)
                        {
                            ComplexFloat result = temp[0];
                            for (int j = prime, k = 0; --j > 0;)
                            {
                                k -= i;
                                k = k < 0 ? k + prime : k;  //k = i * j mod prime.
                                result += temp[j] * simpleRoots[k];
                            }
                            data[index + i * multiplier] = result;
                        }
                        break;
                }
            }
        }

        public static void ReversedBaseIterate(IEnumerable<int> inputFactors, Action<int, int> action)
        {
            List<int> factors = new List<int>();
            List<int> product = new List<int>();
            int n = 1;
            product.Add(n);
            foreach (int factor in inputFactors)
            {
                factors.Add(factor);
                n *= factor;
                product.Add(n);
            }

            int[] values = new int[factors.Count];
            int reversed = 0;
            for (int i = 0; i < n; i++)
            {
                action(i, reversed);

                for (int k = factors.Count; --k >= 0;)
                {
                    values[k]++;
                    reversed += product[k];
                    if (values[k] < factors[k])
                    {
                        break;
                    }
                    values[k] = 0;
                    reversed -= product[k + 1];
                }
            }
        }

        private static int map1(int i, int n)
        {
            i <<= 1;
            return i < n ? i : n * 2 - i - 1;
        }

        private static readonly double INV_SQRT_2 = 0.70710678118654752440084436210485;    // Math.Sqrt(0.5);
        private static readonly float INV_SQRT_2f = 0.70710678118654752440084436210485f;    // Math.Sqrt(0.5);

        public static void EvenDCTType2(this float[] data)
        {
            int N = data.Length;
            if ((N & 1) != 0)
            {
                throw new Exception("DCT operates only on even data sizes.");
            }
            N >>= 1;    //N = number of complex numbers.
            var v = new ComplexFloat[N];
            for (int i = N; --i >= 0;)
            {
                v[i] = new ComplexFloat(data[map1(i * 2, N * 2)], data[map1(i * 2 + 1, N * 2)]);
            }

            FullForwardFFT(v, data);

            var scalej = new ComplexFloat(INV_SQRT_2f, -INV_SQRT_2f);
            Complex firstRoot = Complex.FromPolarCoordinates(1.0, (Math.PI / 4) / N);
            Complex secondRoot = Complex.FromPolarCoordinates(1.0, (Math.PI * 5 / 4) / N);
            Complex firstIterator = INV_SQRT_2 * firstRoot;
            Complex secondIterator = INV_SQRT_2 * secondRoot.MulI();
            for (int i = 1, j = N - 1; i <= j; i++, j--)
            {
                var vi = v[i];
                var vj = v[j].Conjugate;
                var t1 = (vi + vj) * firstIterator.ToComplexFloat();
                var t2 = (vi - vj) * secondIterator.ToComplexFloat();
                var a = (t1 - t2);
                var b = (t1 + t2) * scalej;
                data[i] = a.Real;
                data[N * 2 - i] = a.Imaginary;
                data[j] = b.Real;
                data[N * 2 - j] = -b.Imaginary;
                firstIterator *= firstRoot;
                secondIterator *= secondRoot;
            }

            var number = v[0];
            data[N] = number.Real - number.Imaginary;
            data[0] = number.Real + number.Imaginary;
        }

        public static void EvenDCTType3(this float[] data)
        {
            int N = data.Length;
            if ((N & 1) != 0)
            {
                throw new Exception("Inverse DCT operates only on even data sizes.");
            }
            N >>= 1;
            var v = new ComplexFloat[N];
            v[0] = new ComplexFloat(data[0] + data[N], data[0] - data[N]);
            //A.Conjug * B = (A * B.Conjug).Conjug

            var scalej = new ComplexFloat(INV_SQRT_2f, INV_SQRT_2f);
            Complex firstRoot = Complex.FromPolarCoordinates(1.0, (-Math.PI / 4) / N);
            Complex secondRoot = Complex.FromPolarCoordinates(1.0, (-Math.PI * 5 / 4) / N);
            Complex firstIterator = INV_SQRT_2 * firstRoot;
            Complex secondIterator = INV_SQRT_2 * secondRoot.DivI();
            for (int i = 1, j = N - 1; i <= j; i++, j--)
            {
                var a = new ComplexFloat(data[i], data[N * 2 - i]);
                var b = new ComplexFloat(data[j], -data[N * 2 - j]) * scalej;
                var t1 = (a + b) * firstIterator.ToComplexFloat();
                var t2 = (a - b) * secondIterator.ToComplexFloat();
                v[i] = (t1 - t2);
                v[j] = (t1 + t2).Conjugate;
                firstIterator *= firstRoot;
                secondIterator *= secondRoot;
            }

            FullForwardIFFT(v, data);

            for (int i = N; --i >= 0;)
            {
                var number = v[i];
                data[map1(i * 2, N * 2)] = number.Real;
                data[map1(i * 2 + 1, N * 2)] = number.Imaginary;
            }
        }
    }

    public static class Fourier235UnitTest
    {
        private static void SlowDCT2(float[] data)
        {
            int n = data.Length;
            double[] cos = new double[n * 4];
            for (int i = n * 4; --i >= 0;)
            {
                cos[i] = Math.Cos(Math.PI * 2 / (n * 4) * i);
            }
            double[] result = new double[n];
            double[] array = new double[n];
            for (int i = n; --i >= 0;)
            {
                for (int j = n; --j >= 0;)
                {
                    array[j] = data[j] * cos[(2 * j + 1) * i % (n * 4)];
                }
                result[i] = AccurateSummation.GetAccurateSum(array);
            }
            double scale = Math.Sqrt(2.0);
            for (int i = n; --i > 0;)
            {
                data[i] = (float)(result[i] * scale);
            }
            data[0] = (float)result[0];
        }

        private static void SlowFFT(ComplexFloat[] data)
        {
            int n = data.Length;
            var result = new Complex[n];
            var powers = new Complex[n];
            var input = new Complex[n];
            var array = new Complex[n];
            for (int i = n; --i >= 0;)
            {
                input[i] = new Complex(data[i].Real, data[i].Imaginary);
                powers[i] = Complex.FromPolarCoordinates(1.0, Math.PI * 2 * i / n);
            }
            for (int i = n; --i >= 0;)
            {
                for (int j = n; --j >= 0;)
                {
                    array[j] = input[j] * powers[i * j % n];
                }
                result[i] = new Complex(
                    AccurateSummation.GetAccurateSum(array.Select(x => x.Real)),
                    AccurateSummation.GetAccurateSum(array.Select(x => x.Imaginary)));
            }
            for (int i = n; --i >= 0;)
            {
                data[i] = new ComplexFloat((float)result[i].Real, (float)result[i].Imaginary);
            }
        }

        public static bool DCTUnitTest()
        {
            Random random = new Random(1001);
            int n = 15 * 21 * 4;
            var dcti1 = new float[n];
            var dcti2 = new float[n];
            var dcti3 = new float[n];
            for (int i = n; --i >= 0;)
            {
                dcti1[i] = dcti2[i] = dcti3[i] = (float)(random.NextDouble() * 256 - 128);
            }
            dcti2.EvenDCTType2();
            SlowDCT2(dcti3);
            double maximumError1 = dcti2.GetMaximumError(dcti3) * Math.Sqrt(1.0 / n);
            dcti2.EvenDCTType3();
            float scale = 1.0f / n;
            for (int i = n; --i >= 0;)
            {
                dcti2[i] *= scale;
            }
            double maximumError2 = dcti1.GetMaximumError(dcti2);
            return maximumError1 < 1E-4 && maximumError2 < 1E-4;
        }

        public static bool FFTUnitTest()
        {
            Random random = new Random(1001);
            int n = 15 * 21 * 2;
            var ffti1 = new ComplexFloat[n];
            var ffti2 = new ComplexFloat[n];
            var ffti3 = new ComplexFloat[n];
            for (int i = n; --i >= 0;)
            {
                ffti1[i] = ffti2[i] = ffti3[i] = new ComplexFloat(
                    (float)(random.NextDouble() * 256 - 128),
                    (float)(random.NextDouble() * 256 - 128));
            }
            ffti2.FullForwardFFT();
            SlowFFT(ffti3);
            double maximumError1 = ffti2.GetMaximumError(ffti3) * Math.Sqrt(1.0 / n);
            ffti2.FullForwardIFFT();
            float scale = 1.0f / n;
            for (int i = n; --i >= 0;)
            {
                ffti2[i] *= scale;
            }
            double maximumError2 = ffti1.GetMaximumError(ffti2);
            return maximumError1 < 1E-4 && maximumError2 < 1E-4;
        }
    }
}
