using System;
using System.Collections.Generic;
using System.Linq;
using System.Numerics;
using System.Runtime.InteropServices;
using System.Text;
using System.Threading.Tasks;

namespace Utilities
{
    public static class FourierTransform235Double
    {
        public static void FullForwardFFT(this Complex[] data, double[] temporaryArray = null)
        {
            int N = data.Length;
            Stack<int> factors = new Stack<int>(FourierTransform235.GetSlowFactorization(N).OrderBy(x => x));
            Complex[] roots = null;
            Complex[] temp = null;
            Complex[] simpleRoots = null;
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
            FourierTransform235.ReversedBaseIterate(factors, (i, reversed) => data[reversed] = new Complex(temporaryArray[i * 2 + 0], temporaryArray[i * 2 + 1]));
        }

        public static void FullForwardIFFT(this Complex[] data, double[] temporaryArray = null)
        {
            FullForwardFFT(data, temporaryArray);
            for (int i = 1, j = data.Length - 1; i < j; i++, j--)
            {
                var auxiliary = data[i];
                data[i] = data[j];
                data[j] = auxiliary;
            }
        }

        private static void ReversedPowersFFTSwapAtEnd(Complex[] data, int index, int N, Stack<int> factors, int wpower, int npower,
            ref Complex[] roots, ref Complex[] temp, ref Complex[] simpleRoots)
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

        private const double M_SQRT_3_4 = 0.86602540378443864676372317075294;   // == sqrt(3.0/4) = SIN(PI/3);
        private const double M_ROOT5_C1 = 0.25;
        private const double M_ROOT5_C2 = 0.55901699437494742410229341718282;   // M_SQRT5 / 4;
        private const double M_ROOT5_C3 = 0.58778525229247312916870595463907;   // = sqrt((5 - M_SQRT5) / 8);   // = Im[w^2] = SIN(PI*2/10)
        private const double M_ROOT5_C4 = 0.95105651629515357211643933337938;   //sqrt((5 + M_SQRT5) / 8);      // = Im[w] = SIN(2*PI/5)

        public static ComplexFloat ToComplexFloat(this Complex x)
        {
            return new ComplexFloat((float)x.Real, (float)x.Imaginary);
        }

        public static Complex MulI(this Complex x)
        {
            return new Complex(-x.Imaginary, x.Real);
        }

        public static Complex DivI(this Complex x)
        {
            return new Complex(x.Imaginary, -x.Real);
        }

        private static void DirectFFTWithStartTwiddle(Complex[] data, int index, int multiplier, int length, int wpower, int npower, int prime,
            ref Complex[] roots, ref Complex[] temp, ref Complex[] simpleRoots)
        {
            if (roots == null || roots.Length < prime)
            {
                Array.Resize(ref roots, prime);
            }
            for (int i = prime; --i > 0;)
            {
                roots[i] = Complex.FromPolarCoordinates(1.0, (Math.PI * 2) * (i * wpower) / npower);
            }

            for (; --length >= 0; index++)
            {
                switch (prime)
                {
                    case 2:
                        {
                            var v0 = data[index];
                            var v1 = data[index + multiplier] * roots[1];
                            data[index] = v0 + v1;
                            data[index + multiplier] = v0 - v1;
                        }
                        break;
                    case 3:
                        //w^2 = -1-w
                        {
                            Complex
                                z0 = data[index],
                                z1 = data[index + multiplier] * roots[1],
                                z2 = data[index + 2 * multiplier] * roots[2];
                            Complex
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
                            Complex
                                z0 = data[index],
                                z1 = data[index + multiplier] * roots[1],
                                z2 = data[index + 2 * multiplier] * roots[2],
                                z3 = data[index + 3 * multiplier] * roots[3],
                                z4 = data[index + 4 * multiplier] * roots[4];
                            Complex
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
                                simpleRoots[i] = Complex.FromPolarCoordinates(1.0, (Math.PI * 2) * i / prime);
                            }
                        }
                        if (temp == null || temp.Length < prime)
                        {
                            Array.Resize(ref temp, prime);
                        }
                        var sum = temp[0] = data[index];
                        for (int i = prime; --i > 0;)
                        {
                            temp[i] = data[index + i * multiplier] * roots[i];
                            sum += temp[i];
                        }
                        data[index] = sum;
                        for (int i = prime; --i > 0;)
                        {
                            var result = temp[0];
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

        private static int map1(int i, int n)
        {
            i <<= 1;
            return i < n ? i : n * 2 - i - 1;
        }

        private static readonly double INV_SQRT_2 = 0.70710678118654752440084436210485;    // Math.Sqrt(0.5);

        public static void EvenDCTType2(this double[] data)
        {
            int N = data.Length;
            if ((N & 1) != 0)
            {
                throw new Exception("DCT operates only on even data sizes.");
            }
            N >>= 1;    //N = number of complex numbers.
            var v = new Complex[N];
            for (int i = N; --i >= 0;)
            {
                v[i] = new Complex(data[map1(i * 2, N * 2)], data[map1(i * 2 + 1, N * 2)]);
            }

            FullForwardFFT(v, data);

            var scalej = new Complex(INV_SQRT_2, -INV_SQRT_2);
            for (int i = 1, j = N - 1; i <= j; i++, j--)
            {
                var vi = v[i];
                var vj = Complex.Conjugate(v[j]);
                var t1 = (vi + vj) * Complex.FromPolarCoordinates(INV_SQRT_2, (Math.PI * 2 / 8) * i / N);
                var t2 = (vi - vj) * Complex.FromPolarCoordinates(INV_SQRT_2, (Math.PI * 10 / 8) * i / N).MulI();
                var a = (t1 - t2);
                var b = (t1 + t2) * scalej;
                data[i] = a.Real;
                data[N * 2 - i] = a.Imaginary;
                data[j] = b.Real;
                data[N * 2 - j] = -b.Imaginary;
            }

            var number = v[0];
            data[N] = number.Real - number.Imaginary;
            data[0] = number.Real + number.Imaginary;
        }

        public static void EvenDCTType3(this double[] data)
        {
            int N = data.Length;
            if ((N & 1) != 0)
            {
                throw new Exception("Inverse DCT operates only on even data sizes.");
            }
            N >>= 1;
            var v = new Complex[N];
            v[0] = new Complex(data[0] + data[N], data[0] - data[N]);
            //A.Conjug * B = (A * B.Conjug).Conjug

            var scalej = new Complex(INV_SQRT_2, INV_SQRT_2);
            for (int i = 1, j = N - 1; i <= j; i++, j--)
            {
                var a = new Complex(data[i], data[N * 2 - i]);
                var b = new Complex(data[j], -data[N * 2 - j]) * scalej;
                var t1 = (a + b) * Complex.FromPolarCoordinates(INV_SQRT_2, (-Math.PI * 2 / 8) * i / N);
                var t2 = (a - b) * Complex.FromPolarCoordinates(INV_SQRT_2, (-Math.PI * 10 / 8) * i / N).DivI();
                v[i] = (t1 - t2);
                v[j] = Complex.Conjugate(t1 + t2);
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

    public static class Fourier235DoubleUnitTest
    {
        private static void SlowDCT2(double[] data)
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
                data[i] = result[i] * scale;
            }
            data[0] = result[0];
        }

        private static void SlowFFT(Complex[] input)
        {
            int n = input.Length;
            var result = new Complex[n];
            var powers = new Complex[n];
            var array = new Complex[n];
            for (int i = n; --i >= 0;)
            {
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
                input[i] = result[i];
            }
        }

        public static bool DCTUnitTest()
        {
            Random random = new Random(1001);
            int n = 15 * 21 * 2;
            var dcti1 = new double[n];
            var dcti2 = new double[n];
            var dcti3 = new double[n];
            for (int i = n; --i >= 0;)
            {
                dcti1[i] = dcti2[i] = dcti3[i] = random.NextDouble() * 256 - 128;
            }
            dcti2.EvenDCTType2();
            SlowDCT2(dcti3);
            double maximumError1 = dcti2.GetMaximumError(dcti3) * Math.Sqrt(1.0 / n);
            dcti2.EvenDCTType3();
            double scale = 1.0 / n;
            for (int i = n; --i >= 0;)
            {
                dcti2[i] *= scale;
            }
            double maximumError2 = dcti1.GetMaximumError(dcti2);
            return maximumError1 < 1E-12 && maximumError2 < 1E-12;
        }

        public static bool FFTUnitTest()
        {
            Random random = new Random(1001);
            int n = 15 * 21 * 4;
            var ffti1 = new Complex[n];
            var ffti2 = new Complex[n];
            var ffti3 = new Complex[n];
            for (int i = n; --i >= 0;)
            {
                ffti1[i] = ffti2[i] = ffti3[i] = new Complex(
                    random.NextDouble() * 256 - 128,
                    random.NextDouble() * 256 - 128);
            }
            ffti2.FullForwardFFT();
            SlowFFT(ffti3);
            double maximumError1 = ffti2.GetMaximumError(ffti3) * Math.Sqrt(1.0 / n);
            ffti2.FullForwardIFFT();
            double scale = 1.0 / n;
            for (int i = n; --i >= 0;)
            {
                ffti2[i] *= scale;
            }
            double maximumError2 = ffti1.GetMaximumError(ffti2);
            return maximumError1 < 1E-12 && maximumError2 < 1E-12;
        }
    }
}
