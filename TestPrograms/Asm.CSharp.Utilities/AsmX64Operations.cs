using Bitsy.Core;
using System;
using System.Diagnostics;
using System.Linq;
using System.Numerics;
using System.Runtime.InteropServices;
using System.Windows.Forms;

namespace Utilities
{
    public static class AsmX64Operations
    {
        [DllImport("AsmX64BitBasicOperations.dll", CallingConvention = CallingConvention.Cdecl)]
        public static extern void SetKaratsubaBreakPoint(int karatsubaBreakPoint);

        [DllImport("AsmX64BitBasicOperations.dll", CallingConvention = CallingConvention.Cdecl)]
        public static extern void SetFourierBreakPoint(int fourierBreakPoint);

        [DllImport("AsmX64BitBasicOperations.dll", CallingConvention = CallingConvention.Cdecl)]
        public static extern int GetKaratsubaMultiplicationBufferSize(int n);

        [DllImport("AsmX64BitBasicOperations.dll", CallingConvention = CallingConvention.Cdecl)]
        public static extern int GetKaratsubaSquareBufferSize(int n);

        [DllImport("AsmX64BitBasicOperations.dll", CallingConvention = CallingConvention.Cdecl)]
        public static extern void Negate(ulong[] result, int n);

        [DllImport("AsmX64BitBasicOperations.dll", CallingConvention = CallingConvention.Cdecl)]
        public static extern ulong IntShr(ulong[] result, int shift, int n, bool isSigned = false);

        [DllImport("AsmX64BitBasicOperations.dll", CallingConvention = CallingConvention.Cdecl)]
        public static extern ulong IntShl(ulong[] result, int shift, int n);

        [DllImport("AsmX64BitBasicOperations.dll", CallingConvention = CallingConvention.Cdecl)]
        public static extern ulong AddShl(ulong[] inputOutput, ulong[] input2, byte shift, int n);

        [DllImport("AsmX64BitBasicOperations.dll", CallingConvention = CallingConvention.Cdecl)]
        public static extern ulong SubShl(ulong[] inputOutput, ulong[] input2, byte shift, int n);

        [DllImport("AsmX64BitBasicOperations.dll", CallingConvention = CallingConvention.Cdecl)]
        public static extern byte AddDigit(ulong[] result, ulong digit, int n);

        [DllImport("AsmX64BitBasicOperations.dll", CallingConvention = CallingConvention.Cdecl)]
        public static extern void MultiplyDigit(ulong[] input, ulong digit, ulong[] resultDigitsPlus1, int n, bool isSigned = false);

        [DllImport("AsmX64BitBasicOperations.dll", CallingConvention = CallingConvention.Cdecl)]
        public static extern ulong DivideDigit(ulong[] result, ulong digit, int n, bool isSigned = false);

        [DllImport("AsmX64BitBasicOperations.dll", CallingConvention = CallingConvention.Cdecl)]
        public static extern byte SubtractDigit(ulong[] result, ulong digit, int n);

        [DllImport("AsmX64BitBasicOperations.dll", CallingConvention = CallingConvention.Cdecl)]
        public static extern byte Add(ulong[] input1, ulong[] input2, ulong[] result, byte initialCarry, int n);

        [DllImport("AsmX64BitBasicOperations.dll", CallingConvention = CallingConvention.Cdecl)]
        public static extern byte Subtract(ulong[] input1, ulong[] input2, ulong[] result, byte initialBorrow, int n);

        [DllImport("AsmX64BitBasicOperations.dll", CallingConvention = CallingConvention.Cdecl)]
        public static extern ulong MultiplyDigitAndAdd(ulong[] input, ulong digit, ulong[] result, int n);

        [DllImport("AsmX64BitBasicOperations.dll", CallingConvention = CallingConvention.Cdecl)]
        public static extern ulong MultiplyDigitAndSubtract(ulong[] input, ulong digit, ulong[] result, int n);

        [DllImport("AsmX64BitBasicOperations.dll", CallingConvention = CallingConvention.Cdecl)]
        public static extern void Multiply(ulong[] input1, ulong[] input2, ulong[] result, int n);

        [DllImport("AsmX64BitBasicOperations.dll", CallingConvention = CallingConvention.Cdecl)]
        public static extern void Square(ulong[] input, ulong[] result, int n);

        [DllImport("AsmX64BitBasicOperations.dll", CallingConvention = CallingConvention.Cdecl)]
        public static extern void Karatsuba(ulong[] input1, ulong[] input2, ulong[] result, int n, ulong[] temporaryBuffer = null);

        [DllImport("AsmX64BitBasicOperations.dll", CallingConvention = CallingConvention.Cdecl)]
        public static extern void KaratsubaSquare(ulong[] input, ulong[] result, int n, ulong[] temporaryBuffer = null);

        /// <summary>
        /// Restriction: modulo[0] must / has to be odd otherwise the method won't work due to the Montgomery reduction requirement.
        /// </summary>
        [DllImport("AsmX64BitBasicOperations.dll", CallingConvention = CallingConvention.Cdecl)]
        private static extern void MontgomeryExpMod(ulong[] baseValue, ulong[] power, ulong[] modulo, int n, int cacheBits);

        [DllImport("AsmX64BitBasicOperations.dll", CallingConvention = CallingConvention.Cdecl)]
        private static extern void ExpMod(ulong[] baseValue, ulong[] power, ulong[] modulo, int n, int cacheBits);

        public static void ModularExponentiation(ulong[] baseValue, ulong[] power, ulong[] modulo, int n, int cacheBits)
        {
            if ((modulo[0] & 1) != 0)
            {
                MontgomeryExpMod(baseValue, power, modulo, n, cacheBits);
                return;
            }
            ExpMod(baseValue, power, modulo, n, cacheBits);
        }

        public static void ModularMultiplication(ulong[] input1AndResult, ulong[] input2, ulong[] modulo, int n)
        {
            ulong[] temporary2x = new ulong[n * 2];
            FastestMultiplication(input1AndResult, input2, temporary2x, n, false);
            GetDivMod(temporary2x, modulo, n, null, false);
            for (int i = n; --i >= 0;)
            {
                input1AndResult[i] = temporary2x[i];
            }
        }

        public static void ModularSquare(ulong[] inputAndResult, ulong[] modulo, int n)
        {
            ulong[] temporary2x = new ulong[n * 2];
            FastestSquare(inputAndResult, temporary2x, n, false);
            GetDivMod(temporary2x, modulo, n, null, false);
            for (int i = n; --i >= 0;)
            {
                inputAndResult[i] = temporary2x[i];
            }
        }

        /// <summary>
        /// input2x[0..n-1] = input2x[0..2n-1] 'modulo' modulo[0..n-1]
        /// </summary>
        [DllImport("AsmX64BitBasicOperations.dll", CallingConvention = CallingConvention.Cdecl)]
        public static extern void GetDivMod(ulong[] input2xAndRemainder, ulong[] modulo, int n, ulong[] quotient2x = null, bool isSigned = false);

        /// <summary>
        /// a = a / x 'modulo' p
        /// </summary>
        [DllImport("AsmX64BitBasicOperations.dll", CallingConvention = CallingConvention.Cdecl)]
        private static extern bool DivideModP(ulong[] a, ulong[] x, ulong[] p, int n);
        /// <summary>
        /// a = a / x 'modulo' p
        /// </summary>
        [DllImport("AsmX64BitBasicOperations.dll", CallingConvention = CallingConvention.Cdecl)]
        private static extern void DivideModPSlow(ulong[] a, ulong[] x, ulong[] p, int n);

        [DllImport("AsmX64BitBasicOperations.dll", CallingConvention = CallingConvention.Cdecl)]
        public static extern void GetModFriendlyPrime(ulong[] input2x, ulong primeLowPart, int n);

        [DllImport("AsmX64BitBasicOperations.dll", CallingConvention = CallingConvention.Cdecl)]
        public static extern void RemodFriendlyPrime(ulong[] input1x, ulong primeLowPart, int n);

        [DllImport("AsmX64BitBasicOperations.dll", CallingConvention = CallingConvention.Cdecl)]
        public static extern void GetModFriendlyPrime521(ulong[] input2x, int n);

        [DllImport("AsmX64BitBasicOperations.dll", CallingConvention = CallingConvention.Cdecl)]
        public static extern void RemodFriendlyPrime521(ulong[] input1x, int n);

        [DllImport("AsmX64BitBasicOperations.dll", CallingConvention = CallingConvention.Cdecl)]
        public static extern void AddModP(ulong[] input1, ulong[] input2, ulong[] destination, ulong LowP, int n);

        [DllImport("AsmX64BitBasicOperations.dll", CallingConvention = CallingConvention.Cdecl)]
        public static extern void SubtractModP(ulong[] input1, ulong[] input2, ulong[] destination, ulong LowP, int n);

        [DllImport("AsmX64BitBasicOperations.dll", CallingConvention = CallingConvention.Cdecl)]
        public static extern void AddScaledModP(ulong[] destination, ulong[] input, ulong multiplier, ulong LowP, int n);

        [DllImport("AsmX64BitBasicOperations.dll", CallingConvention = CallingConvention.Cdecl)]
        public static extern void SubScaledModP(ulong[] destination, ulong[] input, ulong multiplier, ulong LowP, int n);


        [DllImport("AsmX64BitBasicOperations.dll", CallingConvention = CallingConvention.Cdecl)]
        public static extern void AddModP521(ulong[] input1, ulong[] input2, ulong[] destination, int n);

        [DllImport("AsmX64BitBasicOperations.dll", CallingConvention = CallingConvention.Cdecl)]
        public static extern void SubtractModP521(ulong[] input1, ulong[] input2, ulong[] destination, int n);

        [DllImport("AsmX64BitBasicOperations.dll", CallingConvention = CallingConvention.Cdecl)]
        public static extern void AddScaledModP521(ulong[] destination, ulong[] input, ulong multiplier, int n);

        [DllImport("AsmX64BitBasicOperations.dll", CallingConvention = CallingConvention.Cdecl)]
        public static extern void SubScaledModP521(ulong[] destination, ulong[] input, ulong multiplier, int n);

        [DllImport("AsmX64BitBasicOperations.dll", CallingConvention = CallingConvention.Cdecl)]
        public static extern void FastestMultiplication(ulong[] input1, ulong[] input2, ulong[] result, int n, bool isSigned = false);

        [DllImport("AsmX64BitBasicOperations.dll", CallingConvention = CallingConvention.Cdecl)]
        public static extern ulong CarrylessMultipyAndXor(ulong[] input, ulong digit, ulong[] result, int n);

        [DllImport("AsmX64BitBasicOperations.dll", CallingConvention = CallingConvention.Cdecl)]
        public static extern void CarrylessMultiplication(ulong[] input1, ulong[] input2, ulong[] result, int n);

        [DllImport("AsmX64BitBasicOperations.dll", CallingConvention = CallingConvention.Cdecl)]
        public static extern void FastestSquare(ulong[] input, ulong[] result, int n, bool isSigned = false);

        [DllImport("AsmX64BitBasicOperations.dll", CallingConvention = CallingConvention.Cdecl)]
        public static extern void FourierMultiplication(ulong[] input1, ulong[] input2, ulong[] result, int n);

        [DllImport("AsmX64BitBasicOperations.dll", CallingConvention = CallingConvention.Cdecl)]
        public static extern void FourierSquare(ulong[] input, ulong[] result, int n);

        public static void DivideModuloPrime(ulong[] a, ulong[] x, ulong[] p, int n)
        {
            DivideModP(a, x, p, n);
        }

        public static ulong[] BytesToULong(this byte[] array)
        {
            ulong[] result = new ulong[(array.Length + 7) >> 3];
            for (int i = array.Length >> 3; --i >= 0;)
            {
                result[i] = BitConverter.ToUInt64(array, i * 8);
            }
            if ((array.Length & 7) != 0)
            {
                ulong lastDigit = 0;
                int start = array.Length & -8;
                for (int i = array.Length; --i >= start;)
                {
                    lastDigit = (lastDigit << 8) | array[i];
                }
                result[result.Length - 1] = lastDigit;
            }
            return result;
        }

        public static ulong[] ToULong(this BigInteger value)
        {
            return value.ToByteArray().BytesToULong();
        }

        public static ulong[] ToULong(this FastInteger value)
        {
            return value.ToByteArray().BytesToULong();
        }

        public static BigInteger GetExtendedEuclideanGCD(BigInteger a, BigInteger b, out BigInteger x)
        {
            BigInteger t1 = a, t2 = b;   //x1 * a + x3 * b = t1
            BigInteger x1 = 1, x2 = 0;   //x2 * a + x4 * b = t2
            while (t2 != 0)
            {
                BigInteger q = BigInteger.DivRem(t1, t2, out t1);
                x1 -= q * x2;
                BigInteger auxiliary;
                auxiliary = t1; t1 = t2; t2 = auxiliary;
                auxiliary = x1; x1 = x2; x2 = auxiliary;
            }
            x = x1;
            return t1;
        }

        public static BigInteger Inverse(BigInteger a, BigInteger p)
        {
            BigInteger x, gcd = GetExtendedEuclideanGCD(a, p, out x);
            if ((x < 0) != (p < 0))
            {
                return x + p;
            }
            return x;
        }

        public static TimeSpan MeasureTime(Action action)
        {
            Stopwatch stopwatch = new Stopwatch();
            stopwatch.Restart();
            action();
            stopwatch.Stop();
            return stopwatch.Elapsed;
        }

        public static void getRandom(ulong[] number, Random random)
        {
            byte[] randomBytes = new byte[number.Length * 8];
            random.NextBytes(randomBytes);
            for (int i = number.Length; --i >= 0;)
            {
                number[i] = BitConverter.ToUInt64(randomBytes, i * 8);
            }
        }

        public static bool UnitTest()
        {
            Random random = new Random(101);
            int n = 2048 / 64;
            ulong[] a = new ulong[n];
            ulong[] b = new ulong[n];
            ulong[] m = new ulong[n];
            ulong[] am1 = new ulong[n];
            ulong[] am3 = new ulong[n];
            ulong[] c1 = new ulong[n * 2];
            ulong[] c2 = new ulong[n * 2];
            ulong[] c3 = new ulong[n * 2];
            ulong[] s0 = new ulong[n * 2];
            ulong[] s1 = new ulong[n * 2];
            ulong[] s2 = new ulong[n * 2];

            BigInteger a0, b0, m0;

            TimeSpan
                elapsedMontgomeryExpMod = TimeSpan.Zero,
                elapsedExpMod = TimeSpan.Zero,
                elapsedBigIntegerExpMod = TimeSpan.Zero;
            bool ok = true;

            for (int iteration = 1; --iteration >= 0;)
            {
                getRandom(a, random);
                getRandom(b, random);
                getRandom(m, random); m[0] &= ulong.MaxValue - 1;
                a0 = new BigInteger(a.SelectMany(l => BitConverter.GetBytes(l)).Concat(Enumerable.Repeat((byte)0, 1)).ToArray());
                b0 = new BigInteger(b.SelectMany(l => BitConverter.GetBytes(l)).Concat(Enumerable.Repeat((byte)0, 1)).ToArray());
                m0 = new BigInteger(m.SelectMany(l => BitConverter.GetBytes(l)).Concat(Enumerable.Repeat((byte)0, 1)).ToArray());

                a.CopyTo(am1, 0);
                ExpMod(am1, b, m, n, 5);

                BigInteger am2 = BigInteger.Zero;
                elapsedBigIntegerExpMod += MeasureTime(() =>
                {
                    am2 = BigInteger.ModPow(a0, b0, m0);
                });
                var bytes1 = am1.SelectMany(l => BitConverter.GetBytes(l)).ToArray();
                var bytes2 = am2.ToByteArray();
                ok &= Enumerable.Range(0, Math.Min(bytes1.Length, bytes2.Length)).All(idx => bytes1[idx] == bytes2[idx]);
            }

            for (int iteration = 1; --iteration >= 0;)
            {
                getRandom(a, random);
                getRandom(b, random);
                getRandom(m, random); m[0] |= 1;
                a.CopyTo(am1, 0);
                a.CopyTo(am3, 0);

                elapsedMontgomeryExpMod += MeasureTime(() =>
                {
                    MontgomeryExpMod(am1, b, m, n, 6);
                });
                elapsedExpMod += MeasureTime(() =>
                {
                    ExpMod(am3, b, m, n, 5);
                });

                ok &= Enumerable.Range(0, n).All(idx => am3[idx] == am1[idx]);
            }

            TimeSpan
                elapsedMulDirect = TimeSpan.Zero,
                elapsedKaratsuba = TimeSpan.Zero;
            ulong[] temporaryKaratsubaBuffer = new ulong[GetKaratsubaMultiplicationBufferSize(n)];
            for (int iteration = 128; --iteration >= 0;)
            {
                getRandom(a, random);
                getRandom(b, random);
                AsmX64Operations.Multiply(a, a, s0, n);

                elapsedMulDirect += MeasureTime(() =>
                {
                    AsmX64Operations.Multiply(a, b, c1, n);
                    AsmX64Operations.Square(a, s1, n);
                });
                elapsedKaratsuba += MeasureTime(() =>
                {
                    AsmX64Operations.Karatsuba(a, b, c2, n, temporaryKaratsubaBuffer);
                    AsmX64Operations.KaratsubaSquare(a, s2, n, temporaryKaratsubaBuffer);
                });

                ok &= Enumerable.Range(0, n * 2).All(idx => c1[idx] == c2[idx]);
                ok &= Enumerable.Range(0, n * 2).All(idx => s0[idx] == s1[idx]);
                ok &= Enumerable.Range(0, n * 2).All(idx => s1[idx] == s2[idx]);
            }
            if (!ok)
            {
                //MessageBox.Show("not ok - error");
                return false;
            }
            //MessageBox.Show(
            //    "elapsedMontgomeryExpMod: " + elapsedMontgomeryExpMod.ToString() + "\r\n" +
            //    "elapsedExpMod: " + elapsedExpMod.ToString() + "\r\n" +
            //    "elapsedBigIntegerExpMod: " + elapsedBigIntegerExpMod.ToString() + "\r\n" +
            //    "normal: " + elapsedMulDirect.ToString() + "  karatsuba: " + elapsedKaratsuba.ToString());
            return true;
        }

        public static bool ECCUnitTest()
        {
            Random random = new Random(104);
            int n = 8;
            bool ok = true;

            TimeSpan
               elapsedInvertFast = TimeSpan.Zero,
               elapsedInvertSlow = TimeSpan.Zero;
            byte[] randomBytes = new byte[16 * n];
            for (int i = 20; --i >= 0;)
            {
                ulong[] a0 = new ulong[n * 2]; a0[0] = 1UL;
                random.NextBytes(randomBytes);
                ulong[] x = Enumerable.Range(0, n * 2).Select(idx => idx >= n ? 0 : BitConverter.ToUInt64(randomBytes, idx * 8)).ToArray();
                random.NextBytes(randomBytes);
                ulong[] p = Enumerable.Range(0, n).Select(idx => BitConverter.ToUInt64(randomBytes, idx * 8)).ToArray();
                p[0] |= 1;
                GetDivMod(a0, p, n);
                GetDivMod(x, p, n);
                ulong[] a1 = a0.ToArray();
                ulong[] a2 = a0.ToArray();

                var fx = new FastInteger(x);
                var fp = new FastInteger(p);
                var f0 = new FastInteger(a0);
                if (fx.ToString() != x.ToBigInteger().ToString() || fp.ToString() != p.ToBigInteger().ToString() || f0.ToString() != a0.ToBigInteger().ToString())
                {
                    ok = false;
                }
                FastInteger idv = f0.DivideModulo(fx, fp);
                elapsedInvertFast += MeasureTime(() =>
                {
                    AsmX64Operations.DivideModuloPrime(a1, x, p, n);
                });
                elapsedInvertSlow += MeasureTime(() =>
                {
                    AsmX64Operations.DivideModPSlow(a2, x, p, n);
                });

                if (a1.ToBigInteger().ToString() != idv.ToString())
                {
                    ok = false;
                }

                ulong[] r1 = new ulong[n * 2];
                Multiply(a1, x, r1, n);
                GetDivMod(r1, p, n);

                ulong[] r2 = new ulong[n * 2];
                Multiply(a2, x, r2, n);
                GetDivMod(r2, p, n);

                if (!Enumerable.Range(0, n).All(idx => r1[idx] == r2[idx]))
                {
                    ok = false;
                }
            }

            TimeSpan
               elapsedExponentiationSlow = TimeSpan.Zero,
               elapsedExponentiation521 = TimeSpan.Zero,
               elapsedExponentiationBigInt = TimeSpan.Zero,
               elapsedExponentiationFastInt = TimeSpan.Zero,
               elapsedExponentiationBigInt521 = TimeSpan.Zero,
               elapsedExponentiationFast = TimeSpan.Zero;
            SECP256K1Base ecc = new SECP256K1Base();
            SECP521R1Base ecc2 = new SECP521R1Base();
            ECCSecP256K1FastInteger fastEcc = new ECCSecP256K1FastInteger();
            ECCSecP521R1 ecc521 = new ECCSecP521R1();

            random = new Random(105);
            OwnECPoint point = ecc.G;
            OwnECPoint point521 = ecc2.G;
            ok = ok && ecc.Verify(point);
            ok = ok && ecc2.Verify(point521);
            ok = ok && ecc2.UnitTest();
            randomBytes = new byte[ecc.N * 8];
            ECC other = new ECC();
            int ops = 1;
            ulong[] order521 = ecc2.Order.ToArray();
            order521[0]--;
            for (int i = ops; --i >= 0;)
            {
                random.NextBytes(randomBytes);
                BigInteger orderInteger = randomBytes.ToBigInteger();
                ulong[] order = Enumerable.Range(0, ecc.N).Select(idx => BitConverter.ToUInt64(randomBytes, idx * 8)).ToArray();
                OwnECPoint point2 = ecc.MultiplyReference(point, order);
                OwnECPoint point3 = new OwnECPoint();
                OwnECPoint point521_O = new OwnECPoint();
                ECPoint point8_521 = null;
                elapsedExponentiation521 += MeasureTime(() =>
                {
                    point521_O = ecc2.MultiplyWithEndomorphism(point521, order521);
                });
                elapsedExponentiationBigInt521 += MeasureTime(() =>
                {
                    point8_521 = ecc521.ECMultiplication(order521.ToBigInteger());
                });
                ok = ok && point521_O.X.IsEqual(point521.X) && point521_O.Y.ToBigInteger() + point521.Y.ToBigInteger() == (BigInteger.One << 521) - BigInteger.One;
                ok = ok && point8_521.X == point521.X.ToBigInteger() && point8_521.Y + point521.Y.ToBigInteger() == (BigInteger.One << 521) - BigInteger.One;
                FastECPoint point6 = ecc.ECMultiplication(randomBytes.ToFastInteger());
                OwnECPoint point4 = new OwnECPoint();
                ECPoint point5 = null;
                ECPoint point7 = null;
                elapsedExponentiationSlow += MeasureTime(() =>
                {
                    point3 = ecc.Multiply(point, order);
                });
                elapsedExponentiationFast += MeasureTime(() =>
                {
                    point4 = ecc.MultiplyWithEndomorphism(point, order);
                });
                elapsedExponentiationBigInt += MeasureTime(() =>
                {
                    point5 = other.ECMultiplication(orderInteger);
                });
                elapsedExponentiationFastInt += MeasureTime(() =>
                {
                    point7 = fastEcc.ECMultiplication(orderInteger);
                });
                ok &= ecc.Verify(point2);
                ok &= ecc.Verify(point3);
                ok &= ecc.Verify(point4);
                ok &= other.Verify(point5);

                if (!ok || !point2.IsSameAs(point3) || !point2.IsSameAs(point4) || !point2.IsSameAs(point5) ||
                    !point2.IsSameAs(point6) || !point2.IsSameAs(point7))
                {
                    ok = false;
                }
            }
            //MessageBox.Show(
            //    "secP521r1 ecc multiplication: " + elapsedExponentiation521.ToString() + " ops/sec = " + (ops / elapsedExponentiation521.TotalSeconds).ToString("N3") + "\r\n" +
            //    "P521R1 BigInt ecc multiplication: " + elapsedExponentiationBigInt521.ToString() + " ops/sec = " + (ops / elapsedExponentiationBigInt521.TotalSeconds).ToString("N3") + "\r\n" +
            //    "Fast ecc multiplication: " + elapsedExponentiationFast.ToString() + " ops/sec = " + (ops / elapsedExponentiationFast.TotalSeconds).ToString("N3") + "\r\n" +
            //    "BigInteger ecc multiplication: " + elapsedExponentiationBigInt.ToString() + " ops/sec = " + (ops / elapsedExponentiationBigInt.TotalSeconds).ToString("N3") + "\r\n" +
            //    "FastInteger ecc multiplication: " + elapsedExponentiationFastInt.ToString() + " ops/sec = " + (ops / elapsedExponentiationFastInt.TotalSeconds).ToString("N3") + "\r\n" +
            //    "Slow ecc multiplication: " + elapsedExponentiationSlow.ToString() + " ops/sec = " + (ops / elapsedExponentiationSlow.TotalSeconds).ToString("N3") + "\r\n", "Time",
            //     MessageBoxButtons.OK, MessageBoxIcon.Information);

            return ok;
        }
    }
}
