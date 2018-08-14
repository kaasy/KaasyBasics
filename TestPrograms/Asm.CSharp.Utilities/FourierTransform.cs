using System;
using System.Collections.Generic;
using System.Linq;
using System.Numerics;
using System.Text;

namespace Utilities
{
    public static class FourierTransform
    {
        static Complex[] roots64 = new Complex[60];
        static Complex[] invertedRoots64 = new Complex[60];
        static ComplexFloat[] roots32 = new ComplexFloat[30];
        static ComplexFloat[] invertedRoots32 = new ComplexFloat[30];

        static FourierTransform()
        {
            for (int i = roots64.Length; --i >= 0;)
            {
                roots64[i] = -ComplexExtensions.Root(4L << i, 3L);
            }
            for (int i = invertedRoots64.Length; --i >= 0;)
            {
                invertedRoots64[i] = -ComplexExtensions.Root(4L << i, -3L);
            }

            for (int i = roots32.Length; --i >= 0;)
            {
                roots32[i] = -ComplexFloat.Root(4L << i, 3L);
            }
            for (int i = invertedRoots32.Length; --i >= 0;)
            {
                invertedRoots32[i] = -ComplexFloat.Root(4L << i, -3L);
            }
        }

        public static void FFT(this Complex[] list)
        {
            int n = list.Length;
            if ((n & (n - 1)) != 0)
            {
                throw new InvalidOperationException("Binary FFT (double) does not operate on non power of 2 lengths.");
            }
            for (int k = n; (k >>= 1) > 0;)
            {
                list.butterflyCore(0, k);
                Complex wpower = Complex.One;
                for (int j = k << 1, x = 0; j < n; x++, j += k << 1)
                {
                    int pos = 0, bit = 1;
                    while ((x & bit) != 0)
                    {
                        bit <<= 1;
                        pos++;
                    }
                    wpower *= roots64[pos];
                    list.forwardFFTCore(j, k, wpower);
                }
            }
            list.BinaryReverseSwap();
        }

        public static void IFFT(this Complex[] list)
        {
            int n = list.Length;
            if ((n & (n - 1)) != 0)
            {
                throw new InvalidOperationException("Binary IFFT (double) does not operate on non power of 2 lengths.");
            }
            list.BinaryReverseSwap();
            for (int k = 1; k < n; k <<= 1)
            {
                list.butterflyCore(0, k);
                Complex wpower = Complex.One;
                for (int j = k << 1, x = 0; j < n; x++, j += k << 1)
                {
                    int pos = 0, bit = 1;
                    while ((x & bit) != 0)
                    {
                        bit <<= 1;
                        pos++;
                    }
                    wpower *= invertedRoots64[pos];
                    list.inverseFFTCore(j, k, wpower);
                }
            }
        }

        private static int bitReverse(int x, int n)
        {
            int result = 0;
            while ((n >>= 1) > 0)
            {
                result = (result << 1) + (x & 1);
                x >>= 1;
            }
            return result;
        }

        private static void knuthFFTRecursive(this ComplexFloat[] list, ComplexFloat[] powers, int start, int length, int powersPosition)
        {
            if ((length >>= 1) <= 0)
            {
                return;
            }
            list.forwardFFTCore(start, length, powers[powersPosition]);
            knuthFFTRecursive(list, powers, start, length, powersPosition << 1);
            knuthFFTRecursive(list, powers, start + length, length, (powersPosition << 1) + 1);
        }

        private static void knuthIFFTRecursive(this ComplexFloat[] list, ComplexFloat[] powers, int start, int length, int powersPosition)
        {
            if ((length >>= 1) <= 0)
            {
                return;
            }
            knuthIFFTRecursive(list, powers, start, length, (powersPosition << 1));
            knuthIFFTRecursive(list, powers, start + length, length, (powersPosition << 1) + 1);
            list.inverseFFTCore(start, length, powers[powersPosition]);
        }

        private static void initializeFFTPowers(int n, ref ComplexFloat[] powers, bool forward)
        {
            int requiredLength = n >> 1;
            if (powers != null && powers.Length >= requiredLength)
            {
                return;
            }
            int existingLength = powers == null ? 0 : powers.Length;
            Array.Resize(ref powers, requiredLength);
            int j = bitReverse(existingLength, requiredLength);
            for (int i = existingLength; i < requiredLength; i++)
            {
                powers[i] = ComplexFloat.Root(n, forward ? j : -j);
                int bit = requiredLength;
                do
                {
                    bit >>= 1;
                } while ((j & bit) != 0);
                j += bit * 3 - requiredLength;
            }
        }

        public static void FFT(this ComplexFloat[] list)
        {
            int n = list.Length;
            if ((n & (n - 1)) != 0)
            {
                throw new InvalidOperationException("Binary FFT (float) does not operate on non power of 2 lengths.");
            }
            ComplexFloat[] powers = null;
            initializeFFTPowers(n, ref powers, true);
            knuthFFTRecursive(list, powers, 0, n, 0);
            list.BinaryReverseSwap();
        }

        public static void IFFT(this ComplexFloat[] list)
        {
            int n = list.Length;
            if ((n & (n - 1)) != 0)
            {
                throw new InvalidOperationException("Binary IFFT (float) does not operate on non power of 2 lengths.");
            }
            ComplexFloat[] powers = null;
            initializeFFTPowers(n, ref powers, false);
            list.BinaryReverseSwap();
            knuthIFFTRecursive(list, powers, 0, n, 0);
        }

        public static void BinaryReverseSwap<T>(this T[] list)
        {
            int n = list.Length;
            for (int i = 1, j = 0; i < n; i++)
            {
                int bit = n;
                do
                {
                    bit >>= 1;
                } while ((j & bit) != 0);
                j += bit * 3 - n;
                if (i < j)
                {
                    T auxiliary = list[i];
                    list[i] = list[j];
                    list[j] = auxiliary;
                }
            }
        }

        public static bool FFTUnitTest()
        {
            Random random = new Random(100111);
            int n = 128 * 1024;
            Complex[] a0 = new Complex[n];
            Complex[] a1 = new Complex[n];
            ComplexFloat[] b0 = new ComplexFloat[n];
            ComplexFloat[] b1 = new ComplexFloat[n];
            for (int i = n; --i >= 0;)
            {
                a0[i] = a1[i] = new Complex(random.NextDouble() * 256 - 128, random.NextDouble() * 256 - 128);
                b0[i] = b1[i] = new ComplexFloat((float)a0[i].Real, (float)a0[i].Imaginary);
            }
            a1.FFT();
            b1.FFT();
            a1.IFFT();
            b1.IFFT();
            float scale = 1.0f / n;
            for (int i = n; --i >= 0;)
            {
                a1[i] /= n;
                b1[i] *= scale;
            }
            double maxErrorA = a0.GetMaximumError(a1);
            double maxErrorB = b0.GetMaximumError(b1);
            return maxErrorA < 1E-9 && maxErrorB < 1E-4f;
        }

        public static double GetMaximumError(this Complex[] list1, Complex[] list2)
        {
            double maximumErrorSquared = 0;
            for (int i = Math.Max(list1.Length, list2.Length); --i >= 0;)
            {
                maximumErrorSquared = Math.Max(maximumErrorSquared, (list1[i] - list2[i]).Energy());
            }
            return Math.Sqrt(maximumErrorSquared);
        }

        public static double GetMaximumError(this ComplexFloat[] list1, ComplexFloat[] list2)
        {
            double maximumErrorSquared = 0;
            for (int i = Math.Max(list1.Length, list2.Length); --i >= 0;)
            {
                maximumErrorSquared = Math.Max(maximumErrorSquared, (list1[i] - list2[i]).EnergyDouble());
            }
            return Math.Sqrt(maximumErrorSquared);
        }

        public static double GetMaximumError(this float[] a, float[] b)
        {
            double maxerr = 0;
            for (int i = Math.Max(a.Length, b.Length); --i >= 0;)
            {
                maxerr = Math.Max(maxerr, Math.Abs(a[i] - b[i]));
            }
            return maxerr;
        }

        public static double GetMaximumError(this double[] a, double[] b)
        {
            double maxerr = 0;
            for (int i = Math.Max(a.Length, b.Length); --i >= 0;)
            {
                maxerr = Math.Max(maxerr, Math.Abs(a[i] - b[i]));
            }
            return maxerr;
        }

        public static RealNumber GetMaximumError(this RealNumber[] a, RealNumber[] b)
        {
            RealNumber maxerr = 0;
            for (int i = Math.Max(a.Length, b.Length); --i >= 0;)
            {
                maxerr = RealNumber.Max(maxerr, RealNumber.Abs(a[i] - b[i]));
            }
            return maxerr;
        }

        public static RealNumber GetMaximumError(this ComplexNumber[] a, ComplexNumber[] b)
        {
            RealNumber maxerr = 0;
            for (int i = Math.Max(a.Length, b.Length); --i >= 0;)
            {
                maxerr = RealNumber.Max(maxerr, ComplexNumber.Abs(a[i] - b[i]));
            }
            return maxerr;
        }
    }

    public static class ComplexExtensions
    {
        public static double Energy(this Complex c)
        {
            return c.Real * c.Real + c.Imaginary * c.Imaginary;
        }

        public static Complex Root(long n, long k)
        {
            double angle = (Math.PI * 2.0) * k / n;
            return new Complex(Math.Cos(angle), Math.Sin(angle));
        }

        public static void butterflyCore(this Complex[] data, int startIndex, int count)
        {
            int middleIndex = startIndex + count;
            for (int i = count; --i >= 0;)
            {
                Complex input1 = data[startIndex + i];
                Complex input2 = data[middleIndex + i];
                data[startIndex + i] = input1 + input2;
                data[middleIndex + i] = input1 - input2;
            }
        }

        public static void forwardFFTCore(this Complex[] data, int startIndex, int count, Complex wpower)
        {
            int middleIndex = startIndex + count;
            for (int i = count; --i >= 0;)
            {
                Complex input1 = data[startIndex + i];
                Complex input2 = data[middleIndex + i] * wpower;
                data[startIndex + i] = input1 + input2;
                data[middleIndex + i] = input1 - input2;
            }
        }

        public static void inverseFFTCore(this Complex[] data, int startIndex, int count, Complex wpower)
        {
            int middleIndex = startIndex + count;
            for (int i = count; --i >= 0;)
            {
                Complex input1 = data[startIndex + i];
                Complex input2 = data[middleIndex + i];
                data[startIndex + i] = input1 + input2;
                data[middleIndex + i] = (input1 - input2) * wpower;
            }
        }

        public static void butterflyCore(this ComplexFloat[] data, int startIndex, int count)
        {
            int middleIndex = startIndex + count;
            for (int i = count; --i >= 0;)
            {
                ComplexFloat input1 = data[startIndex + i];
                ComplexFloat input2 = data[middleIndex + i];
                data[startIndex + i] = input1 + input2;
                data[middleIndex + i] = input1 - input2;
            }
        }

        public static void forwardFFTCore(this ComplexFloat[] data, int startIndex, int count, ComplexFloat wpower)
        {
            int middleIndex = startIndex + count;
            for (int i = count; --i >= 0;)
            {
                ComplexFloat input1 = data[startIndex + i];
                ComplexFloat input2 = data[middleIndex + i] * wpower;
                data[startIndex + i] = input1 + input2;
                data[middleIndex + i] = input1 - input2;
            }
        }

        public static void inverseFFTCore(this ComplexFloat[] data, int startIndex, int count, ComplexFloat wpower)
        {
            int middleIndex = startIndex + count;
            for (int i = count; --i >= 0;)
            {
                ComplexFloat input1 = data[startIndex + i];
                ComplexFloat input2 = data[middleIndex + i];
                data[startIndex + i] = input1 + input2;
                data[middleIndex + i] = (input1 - input2) * wpower;
            }
        }
    }

    public struct ComplexFloat
    {
        public float Real, Imaginary;

        public static readonly ComplexFloat One = new ComplexFloat(1.0f, 0.0f);
        public static readonly ComplexFloat Zero = new ComplexFloat(0.0f, 0.0f);

        public ComplexFloat Conjugate { get { return new ComplexFloat(this.Real, -this.Imaginary); } }

        public ComplexFloat(float real, float imaginary)
        {
            this.Real = real;
            this.Imaginary = imaginary;
        }

        public static ComplexFloat operator -(ComplexFloat x)
        {
            return new ComplexFloat(-x.Real, -x.Imaginary);
        }

        public static ComplexFloat operator +(ComplexFloat x, ComplexFloat y)
        {
            return new ComplexFloat(x.Real + y.Real, x.Imaginary + y.Imaginary);
        }

        public static ComplexFloat operator -(ComplexFloat x, ComplexFloat y)
        {
            return new ComplexFloat(x.Real - y.Real, x.Imaginary - y.Imaginary);
        }

        public static ComplexFloat operator *(ComplexFloat x, float y)
        {
            return new ComplexFloat(x.Real * y, x.Imaginary * y);
        }

        public static ComplexFloat operator *(ComplexFloat x, ComplexFloat y)
        {
            return new ComplexFloat(x.Real * y.Real - x.Imaginary * y.Imaginary, x.Imaginary * y.Real + x.Real * y.Imaginary);
        }

        public static ComplexFloat Root(long n, long k)
        {
            double angle = (Math.PI * 2.0) * k / n;
            return new ComplexFloat((float)Math.Cos(angle), (float)Math.Sin(angle));
        }

        public float Energy()
        {
            return this.Real * this.Real + this.Imaginary * this.Imaginary;
        }

        public double EnergyDouble()
        {
            return (double)this.Real * this.Real + (double)this.Imaginary * this.Imaginary;
        }

        public override string ToString()
        {
            return this.Real.ToString("R") + " + i * " + this.Imaginary.ToString("R");
        }

        public ComplexFloat MulI()
        {
            return new ComplexFloat(-this.Imaginary, this.Real);
        }

        public ComplexFloat DivI()
        {
            return new ComplexFloat(this.Imaginary, -this.Real);
        }

        public static ComplexFloat Root(double angle)
        {
            return new ComplexFloat((float)Math.Cos(angle), (float)Math.Sin(angle));
        }

        public static ComplexFloat FromPolarCoordinates(double magnitude, double angle)
        {
            return new ComplexFloat((float)(magnitude * Math.Cos(angle)), (float)(magnitude * Math.Sin(angle)));
        }
    }
}
