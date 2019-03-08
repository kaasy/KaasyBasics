using System;
using System.Collections.Generic;
using System.Linq;
using System.Numerics;
using System.Runtime.InteropServices;
using System.Text;
using System.Threading.Tasks;

namespace Utilities
{
    public class FFTRealConstants
    {
        public int PrecisionDigits { get; private set; }
        public RealNumber M_SQRT_3_4 { get; private set; }  // = 0.86602540378443864676372317075294;   // == sqrt(3.0/4) = SIN(PI/3);
        public RealNumber M_ROOT5_C2 { get; private set; }  // = 0.55901699437494742410229341718282;   // M_SQRT5 / 4;
        public RealNumber M_ROOT5_C3 { get; private set; }  // = 0.58778525229247312916870595463907;   // = sqrt((5 - M_SQRT5) / 8);   // = Im[w^2] = SIN(PI*2/10)
        public RealNumber M_ROOT5_C4 { get; private set; }  // = 0.95105651629515357211643933337938;   // = sqrt((5 + M_SQRT5) / 8);   // = Im[w  ] = SIN(2*PI/5)
        public RealNumber PI { get; private set; }

        public FFTRealConstants(int precisionDigits)
        {
            this.PrecisionDigits = precisionDigits;
            RealNumber sqrt5 = new RealNumber(5, precisionDigits).GetSqrt();
            this.M_SQRT_3_4 = new RealNumber(3, precisionDigits).GetSqrt() >> 1;
            this.M_ROOT5_C2 = sqrt5 >> 2;
            this.M_ROOT5_C3 = ((5 - sqrt5) >> 3).GetSqrt();
            this.M_ROOT5_C4 = ((5 + sqrt5) >> 3).GetSqrt();
            this.PI = RealNumber.GetPI(precisionDigits);
        }
    }

    public static class FourierTransform235Real
    {
        private static void getIncrement(ref ComplexNumber result,
            List<int> decimalRepresentation, List<int> factorDigits,
            ComplexNumber[] exactProduct)
        {
            int index = decimalRepresentation.Count - 1;
            while (index >= 0 && decimalRepresentation[index] + 1 == factorDigits[index])
            {
                decimalRepresentation[index] = 0;
                //result *= decrementalProduct[index];
                index--;
            }
            if (index >= 0)
            {
                decimalRepresentation[index]++;
                result *= exactProduct[index];
            }
        }

        public static void FullForwardFFT(this ComplexNumber[] data, RealNumber[] temporaryArray = null)
        {
            int N = data.Length;

            int precisionDigits = 1;
            for (int i = data.Length; --i >= 0;)
            {
                precisionDigits = Math.Max(precisionDigits, data[i].GetPrecisionDigits());
            }
            FFTRealConstants precision = new FFTRealConstants(precisionDigits);

            List<int> factors = FourierTransform235.GetSlowFactorization(N).OrderBy(x => x).ToList();
            ComplexNumber[] roots = null;
            ComplexNumber[] temp = null;
            ComplexNumber[] simpleRoots = null;

            List<int> factorDigits = new List<int>();
            List<int> decimalRepresentation = new List<int>();
            List<ComplexNumber> incrementalProduct = new List<ComplexNumber>();

            int totalProduct = 1;
            int N1 = N;
            int count = 1;
            foreach (int factor in factors)
            {
                int N2 = N1;
                N1 /= factor;
                totalProduct *= factor;

                factorDigits.Add(factor);
                incrementalProduct.Clear();
                int partialProduct = totalProduct;
                for (int i = 0; i < factorDigits.Count; i++)
                {
                    incrementalProduct.Add(ComplexNumber.FromPolarAngle((precision.PI << 1) / partialProduct));
                    partialProduct /= factorDigits[i];
                }
                ComplexNumber[] exactProduct = new ComplexNumber[factorDigits.Count - 1];
                if (factorDigits.Count > 1)
                {
                    exactProduct[factorDigits.Count - 2] = incrementalProduct[factorDigits.Count - 2];
                    for (int i = factorDigits.Count - 2; --i >= 0;)
                    {
                        exactProduct[i] = exactProduct[i + 1] * (incrementalProduct[i] * incrementalProduct[i + 2].Conjugate);
                    }
                }

                decimalRepresentation.Clear();
                decimalRepresentation.AddRange(Enumerable.Repeat(0, factorDigits.Count - 1));

                ComplexNumber exactPower = ComplexNumber.One;
                for (int i = 0; i < count; i++)
                {
                    DirectFFTWithStartTwiddle(data, i * N2, N1, N1, factor, precision, ref roots, ref temp, ref simpleRoots, exactPower);
                    getIncrement(ref exactPower, decimalRepresentation, factorDigits, exactProduct);
                }

                count = totalProduct;
            }

            if (temporaryArray == null || temporaryArray.Length < N * 2)
            {
                Array.Resize(ref temporaryArray, N * 2);
            }
            for (int i = N; --i >= 0;)
            {
                var number = data[i];
                temporaryArray[i * 2] = number.Real;
                temporaryArray[i * 2 + 1] = number.Imaginary;
            }
            FourierTransform235.ReversedBaseIterate(factors, (i, reversed) => data[reversed] = new ComplexNumber(temporaryArray[i * 2], temporaryArray[i * 2 + 1]));
        }

        public static void FullForwardIFFT(this ComplexNumber[] data, RealNumber[] temporaryArray = null)
        {
            FullForwardFFT(data, temporaryArray);
            for (int i = 1, j = data.Length - 1; i < j; i++, j--)
            {
                var auxiliary = data[i];
                data[i] = data[j];
                data[j] = auxiliary;
            }
        }

        public static ComplexNumber MulI(this ComplexNumber x)
        {
            return new ComplexNumber(-x.Imaginary, x.Real);
        }

        public static ComplexNumber DivI(this ComplexNumber x)
        {
            return new ComplexNumber(x.Imaginary, -x.Real);
        }

        private static void DirectFFTWithStartTwiddle(ComplexNumber[] data, int index, int multiplier, int length, int prime,
            FFTRealConstants precision, ref ComplexNumber[] roots, ref ComplexNumber[] temp, ref ComplexNumber[] simpleRoots, ComplexNumber requiredPower)
        {
            if (roots == null || roots.Length < prime)
            {
                Array.Resize(ref roots, prime);
            }
            roots[1] = requiredPower;
            for (int i = 2; i < prime; i++)
            {
                roots[i] = roots[i >> 1] * roots[(i + 1) >> 1];
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
                            ComplexNumber
                                z0 = data[index],
                                z1 = data[index + multiplier] * roots[1],
                                z2 = data[index + 2 * multiplier] * roots[2];
                            ComplexNumber
                                t1 = z1 + z2,
                                t2 = z0 - (t1 >> 1),
                                t3 = (z1 - z2).MulI() * precision.M_SQRT_3_4;
                            data[index] = z0 + t1;
                            data[index + multiplier] = t2 + t3;
                            data[index + 2 * multiplier] = t2 - t3;
                        }
                        break;
                    case 5:
                        //http://www2.itap.physik.uni-stuttgart.de/lehre/vorlesungen/SS08/simmeth/fftalgorithms.pdf
                        {
                            ComplexNumber
                                z0 = data[index],
                                z1 = data[index + multiplier] * roots[1],
                                z2 = data[index + 2 * multiplier] * roots[2],
                                z3 = data[index + 3 * multiplier] * roots[3],
                                z4 = data[index + 4 * multiplier] * roots[4];
                            ComplexNumber
                                t1 = z1 + z4,
                                t2 = z2 + z3,
                                t3 = z1 - z4,
                                t4 = z2 - z3,
                                t5 = t1 + t2,
                                t6 = (t1 - t2) * precision.M_ROOT5_C2,
                                t7 = z0 - (t5 >> 2),
                                t8 = t7 + t6,
                                t9 = t7 - t6,
                                t10 = (t3 * precision.M_ROOT5_C4 + t4 * precision.M_ROOT5_C3).MulI(),
                                t11 = (t3 * precision.M_ROOT5_C3 - t4 * precision.M_ROOT5_C4).MulI();
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
                            simpleRoots[1] = ComplexNumber.FromPolarAngle((precision.PI << 1) / prime);
                            for (int i = 2; i < prime; i++)
                            {
                                simpleRoots[i] = simpleRoots[i >> 1] * simpleRoots[(i + 1) >> 1];
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

        public static void EvenDCTType2(this RealNumber[] data, ComplexNumber[] v = null)
        {
            int N = data.Length;
            if ((N & 1) != 0)
            {
                throw new Exception("DCT operates only on even data sizes.");
            }
            N >>= 1;    //N = number of complex numbers.
            Array.Resize(ref v, N);
            for (int i = N; --i >= 0;)
            {
                v[i] = new ComplexNumber(data[map1(i * 2, N * 2)], data[map1(i * 2 + 1, N * 2)]);
            }

            FullForwardFFT(v, data);

            int precisionDigits = 1;
            for (int i = data.Length; --i >= 0;)
            {
                precisionDigits = Math.Max(precisionDigits, data[i].GetPrecisionDigits());
            }
            RealNumber INV_SQRT_2 = new RealNumber(0.5, precisionDigits).GetSqrt();
            RealNumber PI = RealNumber.GetPI(precisionDigits);

            var scalej = new ComplexNumber(INV_SQRT_2, -INV_SQRT_2);
            ComplexNumber firstRoot = ComplexNumber.FromPolarAngle((PI >> 2) / N);
            ComplexNumber secondRoot = ComplexNumber.FromPolarAngle((PI * 5 >> 2) / N);
            ComplexNumber firstIterator = INV_SQRT_2 * firstRoot;
            ComplexNumber secondIterator = INV_SQRT_2 * secondRoot.MulI();
            for (int i = 1, j = N - 1; i <= j; i++, j--)
            {
                var vi = v[i];
                var vj = v[j].Conjugate;
                var t1 = (vi + vj) * firstIterator;
                var t2 = (vi - vj) * secondIterator;
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

        public static void EvenDCTType3(this RealNumber[] data, ComplexNumber[] v = null)
        {
            int N = data.Length;
            if ((N & 1) != 0)
            {
                throw new Exception("Inverse DCT operates only on even data sizes.");
            }
            N >>= 1;
            Array.Resize(ref v, N);
            v[0] = new ComplexNumber(data[0] + data[N], data[0] - data[N]);
            //A.Conjug * B = (A * B.Conjug).Conjug

            int precisionDigits = 1;
            for (int i = data.Length; --i >= 0;)
            {
                precisionDigits = Math.Max(precisionDigits, data[i].GetPrecisionDigits());
            }
            RealNumber INV_SQRT_2 = new RealNumber(0.5, precisionDigits).GetSqrt();
            RealNumber PI = RealNumber.GetPI(precisionDigits);

            var scalej = new ComplexNumber(INV_SQRT_2, INV_SQRT_2);
            ComplexNumber firstRoot = ComplexNumber.FromPolarAngle((-PI >> 2) / N);
            ComplexNumber secondRoot = ComplexNumber.FromPolarAngle((-PI * 5 >> 2) / N);
            ComplexNumber firstIterator = INV_SQRT_2 * firstRoot;
            ComplexNumber secondIterator = INV_SQRT_2 * secondRoot.DivI();
            for (int i = 1, j = N - 1; i <= j; i++, j--)
            {
                var a = new ComplexNumber(data[i], data[N * 2 - i]);
                var b = new ComplexNumber(data[j], -data[N * 2 - j]) * scalej;
                var t1 = (a + b) * firstIterator;
                var t2 = (a - b) * secondIterator;
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

    public static class Fourier235RealUnitTest
    {
        private static void SlowDCT2(RealNumber[] data)
        {
            int precisionDigits = 1;
            for (int i = data.Length; --i >= 0;)
            {
                precisionDigits = Math.Max(precisionDigits, data[i].GetPrecisionDigits());
            }
            RealNumber pi = RealNumber.GetPI(precisionDigits + 1);

            int n = data.Length;
            ComplexNumber[] roots = new ComplexNumber[n * 4];
            roots[0] = ComplexNumber.One;
            roots[1] = ComplexNumber.FromPolarAngle((pi << 1) / (n * 4));
            for (int i = 2; i < n * 4; i++)
            {
                roots[i] = roots[i >> 1] * roots[(i + 1) >> 1];
            }
            var result = new RealNumber[n];
            for (int i = n; --i >= 0;)
            {
                RealNumber sum = RealNumber.Zero;
                for (int j = n; --j >= 0;)
                {
                    sum += data[j] * roots[(2 * j + 1) * i % (n * 4)].Real;
                }
                result[i] = sum;
            }
            RealNumber scale = new RealNumber(2.0, precisionDigits + 1).GetSqrt();
            for (int i = n; --i > 0;)
            {
                data[i] = result[i] * scale;
            }
            data[0] = result[0];
        }

        private static void SlowFFT(ComplexNumber[] input)
        {
            int precisionDigits = 1;
            for (int i = input.Length; --i >= 0;)
            {
                precisionDigits = Math.Max(precisionDigits, input[i].GetPrecisionDigits());
            }
            RealNumber pi = RealNumber.GetPI(precisionDigits + 1);

            int n = input.Length;
            var result = new ComplexNumber[n];
            var powers = new ComplexNumber[n];
            var array = new ComplexNumber[n];
            powers[0] = ComplexNumber.One;
            powers[1] = ComplexNumber.FromPolarAngle((pi << 1) / n);
            for (int i = 2; i < n; i++)
            {
                powers[i] = powers[i >> 1] * powers[(i + 1) >> 1];
            }
            for (int i = n; --i >= 0;)
            {
                ComplexNumber sum = ComplexNumber.Zero;
                for (int j = n; --j >= 0;)
                {
                    sum += input[j] * powers[i * j % n];
                }
                result[i] = sum;
            }
            for (int i = n; --i >= 0;)
            {
                input[i] = result[i];
            }
        }

        public static bool DCTUnitTest()
        {
            int precisionDigits = 10;
            var tolerance = new RealNumber(32 - 64 * precisionDigits, IntegerNumber.One, precisionDigits);
            Random random = new Random(1001);
            int n = 3 * 5 * 7 * 2 * 2;
            var dcti1 = new RealNumber[n];
            var dcti2 = new RealNumber[n];
            var dcti3 = new RealNumber[n];
            for (int i = n; --i >= 0;)
            {
                dcti1[i] = dcti2[i] = dcti3[i] = new RealNumber(random.NextDouble() * 256 - 128, precisionDigits);
            }
            dcti2.EvenDCTType2();
            SlowDCT2(dcti3);
            var maximumError1 = dcti2.GetMaximumError(dcti3) * new RealNumber(n, precisionDigits).GetInverseSqrt();
            dcti2.EvenDCTType3();
            RealNumber scale = new RealNumber(n, precisionDigits).Inverse();
            for (int i = n; --i >= 0;)
            {
                dcti2[i] *= scale;
            }
            var maximumError2 = dcti1.GetMaximumError(dcti2);
            return maximumError1 < tolerance && maximumError2 < tolerance;
        }

        public static bool FFTUnitTest()
        {
            int precisionDigits = 10;
            var tolerance = new RealNumber(32 - 64 * precisionDigits, IntegerNumber.One, precisionDigits);
            Random random = new Random(1001);
            int n = 3 * 5 * 7 * 2;
            var ffti1 = new ComplexNumber[n];
            var ffti2 = new ComplexNumber[n];
            var ffti3 = new ComplexNumber[n];
            for (int i = n; --i >= 0;)
            {
                ffti1[i] = ffti2[i] = ffti3[i] = new ComplexNumber(
                    new RealNumber(random.NextDouble() * 256 - 128, precisionDigits),
                    new RealNumber(random.NextDouble() * 256 - 128, precisionDigits));
            }
            ffti2.FullForwardFFT();
            SlowFFT(ffti3);
            var maximumError1 = ffti2.GetMaximumError(ffti3) * new RealNumber(n, precisionDigits).GetInverseSqrt();
            ffti2.FullForwardIFFT();
            RealNumber scale = new RealNumber(n, precisionDigits).Inverse();
            for (int i = n; --i >= 0;)
            {
                ffti2[i] *= scale;
            }
            var maximumError2 = ffti1.GetMaximumError(ffti2);
            return maximumError1 < tolerance && maximumError2 < tolerance;
        }
    }
}
