using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Utilities
{
    public static class CarrylessMultiplication
    {
        private static ushort[][] digitMultiplier;

        static CarrylessMultiplication()
        {
            digitMultiplier = new ushort[256][];
            for (int i = 256; --i >= 0;)
            {
                ushort[] multiplier = new ushort[256];
                for (int j = 256; --j >= 0;)
                {
                    multiplier[j] = (ushort)slowMultiply(i, j);
                }
                digitMultiplier[i] = multiplier;
            }
        }

        private static int slowMultiply(int a, int b)
        {
            int result = 0;
            for (int i = 32; --i >= 0;)
            {
                if ((a & (1 << i)) != 0)
                {
                    result ^= b << i;
                }
            }
            return result;
        }

        private static ulong multiply(ulong a, ulong b, out ulong productHigh)
        {
            byte[] result = new byte[16];
            for (int i = 8; --i >= 0;)
            {
                ushort[] multiplier = digitMultiplier[(a >> (i * 8)) & 0xFF];
                for (int j = 8; --j >= 0;)
                {
                    ushort x = multiplier[(b >> (j * 8)) & 0xFF];
                    result[i + j] ^= (byte)x;
                    result[i + j + 1] ^= (byte)(x >> 8);
                }
            }
            productHigh = BitConverter.ToUInt64(result, 8);
            return BitConverter.ToUInt64(result, 0);
        }

        private static void carrylessMultiplication(ulong[] input1, ulong[] input2, ulong[] result, int n)
        {
            for (int i = n * 2; --i >= 0;)
            {
                result[i] = 0UL;
            }
            for (int i = n; --i >= 0;)
            {
                for (int j = n; --j >= 0;)
                {
                    ulong hi;
                    result[i + j] ^= multiply(input1[i], input2[j], out hi);
                    result[i + j + 1] ^= hi;
                }
            }
        }

        private static ulong carrylessMultiplyAndXor(ulong[] input, ulong digit, ulong[] result, int n)
        {
            ulong carry = 0UL;
            for (int i = 0; i < n; i++)
            {
                ulong hi;
                result[i] ^= carry ^ multiply(digit, input[i], out hi);
                carry = hi;
            }
            return carry;
        }

        public static bool UnitTest(int seed)
        {
            Random random = new Random(seed);
            int n = 891;
            ulong[] a = new ulong[n];
            ulong[] b = new ulong[n];
            ulong[] s1 = new ulong[n + 1];
            ulong[] s2 = new ulong[n + 1];
            ulong[] r1 = new ulong[n * 2];
            ulong[] r2 = new ulong[n * 2];
            for (int i = n; --i >= 0;)
            {
                byte[] data = new byte[16];
                random.NextBytes(data);
                a[i] = BitConverter.ToUInt64(data, 0);
                b[i] = BitConverter.ToUInt64(data, 8);
            }
            bool ok = true;

            s1[n] = AsmX64Operations.CarrylessMultipyAndXor(a, b[0], s1, n);
            s2[n] = carrylessMultiplyAndXor(a, b[0], s2, n);
            for (int i = 0; i <= n ; i++)
            {
                if (s1[i] != s2[i])
                {
                    ok = false;
                }
            }

            AsmX64Operations.CarrylessMultiplication(a, b, r1, n);
            carrylessMultiplication(a, b, r2, n);
            for (int i = 0; i < n * 2; i++)
            {
                if (r1[i] != r2[i])
                {
                    ok = false;
                }
            }

            return ok;
        }
    }
}
