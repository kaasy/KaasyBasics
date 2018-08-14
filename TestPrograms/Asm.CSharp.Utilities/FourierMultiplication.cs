using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Numerics;
using System.Runtime.InteropServices;
using System.Text;
using System.Windows.Forms;

namespace Utilities
{
    [StructLayout(LayoutKind.Sequential)]
    public struct FourierPoint
    {
        [DllImport("AsmX64BitBasicOperations.dll", CallingConvention = CallingConvention.Cdecl)]
        public static extern void mulmod(ref FourierPoint a, ref FourierPoint b);

        [DllImport("AsmX64BitBasicOperations.dll", CallingConvention = CallingConvention.Cdecl)]
        public static extern void butterfly(ref FourierPoint a, ref FourierPoint b);

        [DllImport("AsmX64BitBasicOperations.dll", CallingConvention = CallingConvention.Cdecl)]
        public static extern void addmod(ref FourierPoint a, ref FourierPoint b);

        [DllImport("AsmX64BitBasicOperations.dll", CallingConvention = CallingConvention.Cdecl)]
        public static extern void submod(ref FourierPoint a, ref FourierPoint b);

        [DllImport("AsmX64BitBasicOperations.dll", CallingConvention = CallingConvention.Cdecl)]
        public static extern void forwardFFTCore(IntPtr number, int startIndex, int n, ref FourierPoint wpower);

        [DllImport("AsmX64BitBasicOperations.dll", CallingConvention = CallingConvention.Cdecl)]
        public static extern void inverseFFTCore(IntPtr number, int startIndex, int n, ref FourierPoint wpower);

        [DllImport("AsmX64BitBasicOperations.dll", CallingConvention = CallingConvention.Cdecl)]
        public static extern void butterflyCore(IntPtr number, int startIndex, int n);

        [DllImport("AsmX64BitBasicOperations.dll", CallingConvention = CallingConvention.Cdecl)]
        public static extern void PropagateCarries(IntPtr number, int bitsPerDigit, int numberLength);

        [DllImport("AsmX64BitBasicOperations.dll", CallingConvention = CallingConvention.Cdecl)]
        public static extern void OpMov(IntPtr destination, IntPtr source, int log2BitsPerOp, int n);

        [DllImport("AsmX64BitBasicOperations.dll", CallingConvention = CallingConvention.Cdecl)]
        public static extern void OpNot(IntPtr sourceAndDestination, int log2BitsPerOp, int n);

        [DllImport("AsmX64BitBasicOperations.dll", CallingConvention = CallingConvention.Cdecl)]
        public static extern void OpAdd(IntPtr source1, IntPtr source2, IntPtr destination, int log2BitsPerOp, int n);

        [DllImport("AsmX64BitBasicOperations.dll", CallingConvention = CallingConvention.Cdecl)]
        public static extern void OpSub(IntPtr source1, IntPtr source2, IntPtr destination, int log2BitsPerOp, int n);

        [DllImport("AsmX64BitBasicOperations.dll", CallingConvention = CallingConvention.Cdecl)]
        public static extern void OpAnd(IntPtr source1, IntPtr source2, IntPtr destination, int log2BitsPerOp, int n);

        [DllImport("AsmX64BitBasicOperations.dll", CallingConvention = CallingConvention.Cdecl)]
        public static extern void OpOr(IntPtr source1, IntPtr source2, IntPtr destination, int log2BitsPerOp, int n);

        [DllImport("AsmX64BitBasicOperations.dll", CallingConvention = CallingConvention.Cdecl)]
        public static extern void OpXor(IntPtr source1, IntPtr source2, IntPtr destination, int log2BitsPerOp, int n);

        [DllImport("AsmX64BitBasicOperations.dll", CallingConvention = CallingConvention.Cdecl)]
        public static extern void OpMul(IntPtr source1, IntPtr source2, IntPtr destination, int log2BitsPerOp, int n);

        [DllImport("AsmX64BitBasicOperations.dll", CallingConvention = CallingConvention.Cdecl)]
        public static extern void OpCLMul(IntPtr source1, IntPtr source2, IntPtr destination, int log2BitsPerOp, int n);

        [DllImport("AsmX64BitBasicOperations.dll", CallingConvention = CallingConvention.Cdecl)]
        public static extern void OpFPAdd(IntPtr source1, IntPtr source2, IntPtr destination, int log2BitsPerOp, int n);

        [DllImport("AsmX64BitBasicOperations.dll", CallingConvention = CallingConvention.Cdecl)]
        public static extern void OpFPSub(IntPtr source1, IntPtr source2, IntPtr destination, int log2BitsPerOp, int n);

        [DllImport("AsmX64BitBasicOperations.dll", CallingConvention = CallingConvention.Cdecl)]
        public static extern void OpFPMul(IntPtr source1, IntPtr source2, IntPtr destination, int log2BitsPerOp, int n);

        [DllImport("AsmX64BitBasicOperations.dll", CallingConvention = CallingConvention.Cdecl)]
        public static extern void OpFPDiv(IntPtr source1, IntPtr source2, IntPtr destination, int log2BitsPerOp, int n);

        [DllImport("AsmX64BitBasicOperations.dll", CallingConvention = CallingConvention.Cdecl)]
        public static extern void PowMod(ref FourierPoint a, ref FourierPoint power);

        public static readonly FourierPoint PrimeModulo = new FourierPoint((1023UL << 54) + 1UL, ulong.MaxValue);
        public static readonly FourierPoint One = new FourierPoint(1UL, 0UL);
        public static readonly FourierPoint Zero = new FourierPoint(0UL, 0UL);
        public static readonly FourierPoint MinusOne = new FourierPoint(PrimeModulo.Low - 1, PrimeModulo.High);
        public static readonly FourierPoint HalfOne = new FourierPoint((2047UL << 53) + 1UL, ulong.MaxValue >> 1);  // = 1 / 2
        public static readonly FourierPoint Generator = new FourierPoint(5UL, 0UL);
        public static readonly FourierPoint Root = new FourierPoint(0x4087ea9b3d2c84a8UL, 0xfaaf0d013359bbb3UL);    //Root = PowMod(Generator, new FourierPoint(ulong.MaxValue, 1023UL));
        public static readonly FourierPoint[] Roots;
        public static readonly FourierPoint[] RootsMultipliers, RootsInverters;
        public const ulong InverseLowP = 0x0040000000000001UL;
        public const ulong NegativeInverseLowP = 0xFFBFFFFFFFFFFFFFUL;
        public static readonly ulong[] PrimeModuloDigits = new ulong[2] { PrimeModulo.Low, PrimeModulo.High };

        static FourierPoint()
        {
            FourierPoint p = PrimeModulo;
            p.Low -= 1;
            FourierPoint halfp = new FourierPoint((p.Low >> 1) | (p.High << 63), p.High >> 1);
            var fp1 = PowMod(Generator, p);
            var fp2 = PowMod(Generator, halfp);
            remodP(ref fp1);
            remodP(ref fp2);
            var fp3 = PowMod(Root, new FourierPoint(1UL << 53, 0));
            if (fp1 != One || fp2 != MinusOne || fp3 != MinusOne)
            {
                throw new InvalidOperationException("Generator does not pass unit test.");
            }
            Roots = new FourierPoint[54];
            FourierPoint[] inverters = new FourierPoint[Roots.Length];
            Roots[Roots.Length - 1] = Root;
            inverters[Roots.Length - 1] = FourierPoint.PowMod(Root, new FourierPoint(PrimeModulo.Low - 2, PrimeModulo.High));
            for (int i = Roots.Length - 1; --i >= 0;)
            {
                Roots[i] = Roots[i + 1] * Roots[i + 1];
                remodP(ref Roots[i]);
                inverters[i] = inverters[i + 1] * inverters[i + 1];
                remodP(ref inverters[i]);
            }
            RootsMultipliers = new FourierPoint[Roots.Length - 1];
            RootsInverters = new FourierPoint[Roots.Length - 1];
            for (int i = 0; i < RootsMultipliers.Length; i++)
            {
                RootsMultipliers[i] = PrimeModulo - Roots[i] * Roots[i + 1];
                RootsInverters[i] = PrimeModulo - inverters[i] * inverters[i + 1];
            }
        }

        public static void remodP(ref FourierPoint a)
        {
            if (a.High == ulong.MaxValue && a.Low >= (1023UL << 54) + 1UL)
            {
                a.Low -= (1023UL << 54) + 1UL;
                a.High = 0;
            }
        }

        public static FourierPoint PowMod(FourierPoint a, FourierPoint power)
        {
            PowMod(ref a, ref power);
            return a;
        }

        //public static FourierPoint PowMod(FourierPoint a, FourierPoint power)
        //{
        //    if (power.High == 0 && power.Low <= 1)
        //    {
        //        return power.Low == 0 ? new FourierPoint(1UL, 0UL) : a;
        //    }
        //    FourierPoint result = PowMod(a, new FourierPoint((power.Low >> 1) | (power.High << 63), power.High >> 1));
        //    mulmod(ref result, ref result);
        //    if ((power.Low & 1UL) != 0)
        //    {
        //        mulmod(ref result, ref a);
        //    }
        //    return result;
        //}

        public ulong Low;
        public ulong High;

        public static FourierPoint operator *(FourierPoint a, FourierPoint b)
        {
            mulmod(ref a, ref b);
            return a;
        }

        public static FourierPoint operator +(FourierPoint a, FourierPoint b)
        {
            addmod(ref a, ref b);
            return a;
        }

        public static FourierPoint operator -(FourierPoint a, FourierPoint b)
        {
            submod(ref a, ref b);
            return a;
        }

        public static FourierPoint operator >>(FourierPoint x, int shift)
        {
            ulong[] number = new ulong[2] { x.Low, x.High };
            while (shift >= 64)
            {
                ulong carry = AsmX64Operations.MultiplyDigitAndAdd(PrimeModuloDigits, NegativeInverseLowP * number[0], number, 2);
                number[0] = number[1];
                number[1] = carry;
                shift -= 64;
            }
            if (shift != 0)
            {
                ulong mask = (1UL << shift) - 1UL;
                ulong carry = AsmX64Operations.MultiplyDigitAndAdd(PrimeModuloDigits, NegativeInverseLowP * number[0] & mask, number, 2);
                return new FourierPoint((number[0] >> shift) | (number[1] << (64 - shift)), (number[1] >> shift) | (carry << (64 - shift)));
            }
            return new FourierPoint(number[0], number[1]);
        }

        public static void ShiftRightDirect(ref FourierPoint x, int shift)
        {
            x = new FourierPoint((x.Low >> shift) | (x.High << (64 - shift)), x.High >> shift);
        }

        public static bool operator ==(FourierPoint a, FourierPoint b)
        {
            return a.Low == b.Low && a.High == b.High;
        }

        public static bool operator !=(FourierPoint a, FourierPoint b)
        {
            return !(a == b);
        }

        public FourierPoint(ulong lo, ulong hi) : this()
        {
            this.Low = lo;
            this.High = hi;
        }

        public override string ToString()
        {
            return this.High.ToString("X16") + this.Low.ToString("X16");
        }

        public override int GetHashCode()
        {
            return this.Low.GetHashCode() + this.High.GetHashCode();
        }

        public override bool Equals(object obj)
        {
            if (!(obj is FourierPoint))
            {
                return false;
            }
            FourierPoint other = (FourierPoint)obj;
            return this == other;
        }
    }

    public static class FourierMultiplication
    {
        public static void FFT(this FourierPoint[] list, bool performBase2Swap = true)
        {
            int n = list.Length;
            if ((n & (n - 1)) != 0)
            {
                throw new InvalidOperationException("This FFT version is only supported on power of 2 lengths.");
            }
            list.fft();
            if (performBase2Swap)
            {
                list.BinaryReverseSwap();
            }
        }

        public static void IFFT(this FourierPoint[] list, bool performBase2Swap = true)
        {
            int n = list.Length;
            if ((n & (n - 1)) != 0)
            {
                throw new InvalidOperationException("This IFFT version is only supported on power of 2 lengths.");
            }
            if (performBase2Swap)
            {
                list.BinaryReverseSwap();
            }
            list.ifft();
        }

        private static void fft(this FourierPoint[] list)
        {
            int n = list.Length;
            var handle = GCHandle.Alloc(list, GCHandleType.Pinned);
            try
            {
                IntPtr listAddress = handle.AddrOfPinnedObject();
                for (int k = n; (k >>= 1) > 0;)
                {
                    FourierPoint.butterflyCore(listAddress, 0, k);
                    FourierPoint wpower = FourierPoint.One;
                    for (int j = k << 1, x = 0; j < n; x++, j += k << 1)
                    {
                        int pos = 0;
                        while ((x & (1 << pos)) != 0)
                        {
                            pos++;
                        }
                        FourierPoint.mulmod(ref wpower, ref FourierPoint.RootsMultipliers[pos]);
                        FourierPoint.forwardFFTCore(listAddress, j, k, ref wpower);
                    }
                }
            }
            finally
            {
                handle.Free();
            }
        }

        private static void ifft(this FourierPoint[] list)
        {
            int n = list.Length;
            var handle = GCHandle.Alloc(list, GCHandleType.Pinned);
            try
            {
                IntPtr listAddress = handle.AddrOfPinnedObject();
                for (int k = 1; k < n; k <<= 1)
                {
                    FourierPoint.butterflyCore(listAddress, 0, k);
                    FourierPoint wpower = FourierPoint.One;
                    for (int j = k << 1, x = 0; j < n; x++, j += k << 1)
                    {
                        int pos = 0;
                        while ((x & (1 << pos)) != 0)
                        {
                            pos++;
                        }
                        FourierPoint.mulmod(ref wpower, ref FourierPoint.RootsInverters[pos]);
                        FourierPoint.inverseFFTCore(listAddress, j, k, ref wpower);
                    }
                }
            }
            finally
            {
                handle.Free();
            }
        }

        public static bool UnitTest(int seed)
        {
            Random random = new Random(seed);
            int bits = 17;
            int n = 1 << bits;
            FourierPoint[] a0 = new FourierPoint[n];
            FourierPoint[] a1 = new FourierPoint[n];
            byte[] bytes = new byte[16];
            for (int i = n; --i >= 0;)
            {
                random.NextBytes(bytes);
                a0[i] = a1[i] = new FourierPoint(BitConverter.ToUInt64(bytes, 0), BitConverter.ToUInt64(bytes, 8));
            }
            a0.FFT(false);
            a0.IFFT(false);
            //var scale = FourierPoint.PowMod(FourierPoint.HalfOne, new FourierPoint((ulong)bits, 0));
            for (int i = n; --i >= 0;)
            {
                //a0[i] *= scale;
                a0[i] >>= bits;
            }
            for (int i = n; --i >= 0;)
            {
                if (a0[i] != a1[i])
                {
                    return false;
                }
            }
            return true;
        }

        public static bool UnitTest()
        {
            FourierPoint a, b;
            Random random = new Random(1001);
            byte[] ra = new byte[32];
            BigInteger p = (BigInteger.One << 128) - (BigInteger.One << 54) + BigInteger.One;
            for (int i = 20 * 1000; --i >= 0;)
            {
                random.NextBytes(ra);
                ra[ra.Length - 1] = 0;
                a = new FourierPoint(BitConverter.ToUInt64(ra, 0), BitConverter.ToUInt64(ra, 8));
                b = new FourierPoint(BitConverter.ToUInt64(ra, 16), BitConverter.ToUInt64(ra, 24));
                BigInteger ba = new BigInteger(ra.Take(16).Concat(Enumerable.Repeat((byte)0, 1)).ToArray());
                BigInteger bb = new BigInteger(ra.Skip(16).Take(16).Concat(Enumerable.Repeat((byte)0, 1)).ToArray());
                FourierPoint.mulmod(ref a, ref b);
                ba = (ba * bb) % p;
                byte[] rbb1 = ba.ToByteArray();
                Array.Resize(ref rbb1, 16);
                byte[] rbb2 = new byte[16];
                BitConverter.GetBytes(a.Low).CopyTo(rbb2, 0);
                BitConverter.GetBytes(a.High).CopyTo(rbb2, 8);
                for (int k = 16; --k >= 0;)
                {
                    if (rbb1[k] != rbb2[k])
                    {
                        return false;
                    }
                }
            }
            return true;
        }

        private static int CompareBaseBits(int bitsPerDigit, long n)
        {   //middle digit <= n * (base - 1) ^ 2 + [2^(128 - bitsPerDigit)] <= P - 1
            ulong[] value = new ulong[1];
            ulong[] square = new ulong[value.Length * 2];
            ulong[] multiplier = new ulong[square.Length];
            ulong[] result = new ulong[square.Length * 2];
            ulong[] adder = new ulong[result.Length];
            value[0] = (1UL << bitsPerDigit) - 1;
            multiplier[0] = (ulong)n;
            AsmX64Operations.Square(value, square, value.Length);
            AsmX64Operations.Multiply(square, multiplier, result, square.Length);
            if (128 - bitsPerDigit < 64)
            {
                adder[0] = 1UL << (128 - bitsPerDigit);
            }
            else
            {
                adder[1] = 1UL << (64 - bitsPerDigit);
            }
            AsmX64Operations.Add(result, adder, result, 0, result.Length);

            if (result[2] > 0 || result[3] > 0 || (result[0] >= FourierPoint.PrimeModulo.Low && result[1] == ulong.MaxValue))
            {
                return 1;
            }
            if (result[0] == FourierPoint.PrimeModulo.Low - 1 && result[1] == ulong.MaxValue)
            {
                return 0;
            }
            return -1;
        }

        public const double LOG2_E = 1.4426950408889634073599246810019;    // = log2(e) = 1 / ln(2)

        public static void GetFFTParameters(int n, out int bitsPerDigit, out int logFFTSize)
        {
            long inputBits = n * 64L;
            bitsPerDigit = 63;
            long fftSize = (inputBits + bitsPerDigit - 1) / bitsPerDigit;
            bitsPerDigit = Math.Min(63, (int)Math.Ceiling(64 - Math.Log(fftSize) * (LOG2_E * 0.5)));  //use bitsPerDigit <= 63
            while (true)
            {
                fftSize = (inputBits + bitsPerDigit - 1) / bitsPerDigit;
                if (CompareBaseBits(bitsPerDigit, fftSize) <= 0)
                {
                    logFFTSize = 1 + (int)Math.Ceiling(Math.Log(fftSize) * LOG2_E);
                    return;
                }
                bitsPerDigit--;
            }
        }

        private const double LOG2_3 = 1.5849625007211561814537389439478;   // == log2(3) == ln(3) / ln(2)

        public static ulong ReadBits(ulong[] number, long bitPosition, int count)
        {
            long index = bitPosition >> 6;
            int shift = (int)bitPosition & 63;
            ulong result = number[index] >> shift;
            shift = 64 - shift;
            if (shift < count && index + 1 < number.Length)
            {
                result |= number[index + 1] << shift;
            }
            result &= (1UL << count) - 1;
            return result;
        }

        public static void WriteBits(ulong[] number, long bitPosition, int count, ulong value)
        {
            long index = bitPosition >> 6;
            int shift = (int)bitPosition & 63;
            ulong mask = (1UL << count) - 1UL;
            value &= mask;
            number[index] = (number[index] & (~(mask << shift))) | (value << shift);
            shift = 64 - shift;
            if (shift < count && index + 1 < number.Length)
            {
                number[index + 1] = (number[index + 1] & (~(mask >> shift))) | (value >> shift);
            }
        }

        private static void GroupDigits(ulong[] number, int bitsPerDigit, FourierPoint[] result)
        {
            long numberOfDigits = Math.Min(result.Length, (number.LongLength * 64 + bitsPerDigit - 1) / bitsPerDigit);
            long position = numberOfDigits * bitsPerDigit;
            for (long i = numberOfDigits; --i >= 0;)
            {
                position -= bitsPerDigit;
                result[i] = new FourierPoint(ReadBits(number, position, bitsPerDigit), 0UL);
            }
        }

        private static void UngroupDigits(FourierPoint[] number, int bitsPerDigit, ulong[] result)
        {
            long numberOfDigits = Math.Min(number.Length, (result.Length * 64 + bitsPerDigit - 1) / bitsPerDigit);
            long position = numberOfDigits * bitsPerDigit;
            for (long i = numberOfDigits; --i >= 0;)
            {
                position -= bitsPerDigit;
                WriteBits(result, position, bitsPerDigit, number[i].Low);
            }
        }

        public static void Multiply(ulong[] input1, ulong[] input2, ulong[] result, int n)
        {
            // AsmX64Operations.FastestMultiplication(input1, input2, result, n);
            // return;

            int fftBitsPerDigit, fftLog2Size;
            GetFFTParameters(n, out fftBitsPerDigit, out fftLog2Size);
            long fftN = 1L << fftLog2Size;

            //double directComputations = (double)n * n;
            //double karatsubaComputations = 4.7 * Math.Pow(n, LOG2_3);
            //double fftComputations = 7.2 * (3.0 * fftLog2Size + 4.0) * fftN;

            FourierPoint[] number1 = new FourierPoint[fftN];
            GroupDigits(input1, fftBitsPerDigit, number1);
            FourierPoint[] number2 = new FourierPoint[fftN];
            GroupDigits(input2, fftBitsPerDigit, number2);

            //FourierPoint scale = FourierPoint.PowMod(FourierPoint.HalfOne, new FourierPoint((ulong)fftLog2Size, 0UL));
            number1.FFT(false);
            number2.FFT(false);
            for (int i = number1.Length; --i >= 0;)
            {
                FourierPoint point = number1[i];
                FourierPoint.mulmod(ref point, ref number2[i]);
                //FourierPoint.mulmod(ref point, ref scale);
                point >>= fftLog2Size;
                number1[i] = point;
            }
            number1.IFFT(false);

            var number1Handle = GCHandle.Alloc(number1, GCHandleType.Pinned);
            try
            {
                FourierPoint.PropagateCarries(number1Handle.AddrOfPinnedObject(), fftBitsPerDigit, number1.Length);
            }
            finally
            {
                number1Handle.Free();
            }

            UngroupDigits(number1, fftBitsPerDigit, result);
        }

        public static bool UnitTestBigMul(int seed)
        {
            Random random = new Random(seed);
            byte[] data = new byte[16];
            int length = 6 * 1024;   //1024;
            ulong[] xx = new ulong[length];
            ulong[] yy = new ulong[length];
            ulong[] zz = new ulong[length * 2];
            ulong[] tt = new ulong[length * 2];
            ulong[] uu = new ulong[length * 2];

            int tries = 1;

            ulong[] temporaryBuffer = new ulong[AsmX64Operations.GetKaratsubaMultiplicationBufferSize(xx.Length)];
            TimeSpan directTime = TimeSpan.Zero;
            TimeSpan karatsubaTime = TimeSpan.Zero;
            TimeSpan fourierTime = TimeSpan.Zero;
            for (int step = tries; --step >= 0;)
            {
                for (int i = xx.Length; --i >= 0;)
                {
                    random.NextBytes(data);
                    xx[i] = BitConverter.ToUInt64(data, 0);
                    yy[i] = BitConverter.ToUInt64(data, 8);
                }

                directTime += AsmX64Operations.MeasureTime(() =>
                {
                    AsmX64Operations.Multiply(xx, yy, zz, length);
                });
                karatsubaTime += AsmX64Operations.MeasureTime(() =>
                {
                    AsmX64Operations.Karatsuba(xx, yy, tt, length);
                });
                fourierTime += AsmX64Operations.MeasureTime(() =>
                {
                    AsmX64Operations.FastestMultiplication(xx, yy, uu, length);
                });

                for (int i = length * 2; --i >= 0;)
                {
                    if (zz[i] != tt[i] || tt[i] != uu[i])
                    {
                        return false;
                    }
                }
            }

            int fftBitsPerDigit, fftLog2Size;
            GetFFTParameters(length, out fftBitsPerDigit, out fftLog2Size);

            double directOpsSecond = (1.0 * tries * length * length / directTime.TotalSeconds);
            double karaOpsSecond = tries * Math.Pow(length, LOG2_3) / karatsubaTime.TotalSeconds;
            double fftOpsSecond = tries * Math.Pow(2.0, fftLog2Size) * (3.0 * fftLog2Size + 4.0) / fourierTime.TotalSeconds;
            //MessageBox.Show(
            //    "Direct: " + directTime.ToString() + "   " + directOpsSecond.ToString("N3") + " ops/sec\r\n" +
            //    "Karatsuba: " + karatsubaTime.ToString() + "   " + (directOpsSecond / karaOpsSecond).ToString("N3") + " x constant\r\n" +
            //    "Fourier: " + fourierTime.ToString() + "   " + (directOpsSecond / fftOpsSecond).ToString("N3") + " x constant\r\n", "Timings");

            return true;
        }
    }
}