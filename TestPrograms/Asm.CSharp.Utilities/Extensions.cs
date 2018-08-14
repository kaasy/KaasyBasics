using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Numerics;
using Bitsy.Core;

namespace Utilities
{
    public static class Extensions
    {
        public static byte[] ToByteArray(this ulong[] data)
        {
            return data.SelectMany(x => BitConverter.GetBytes(x)).ToArray();
        }

        public static BigInteger ToBigInteger(this IEnumerable<byte> bytes)
        {
            return new BigInteger(bytes.Concat(Enumerable.Repeat((byte)0, 1)).ToArray());
        }

        public static FastInteger ToFastInteger(this IEnumerable<byte> bytes)
        {
            return new FastInteger(bytes.Concat(Enumerable.Repeat((byte)0, 1)).ToArray());
        }

        public static FastInteger Square(this FastInteger x)
        {
            return x * x;
        }

        public static bool IsEqual(this ulong[] data1, ulong[] data2)
        {
            if (data1.Length != data2.Length)
            {
                return false;
            }
            for (int i = data1.Length; --i >= 0;)
            {
                if (data1[i] != data2[i])
                {
                    return false;
                }
            }
            return true;
        }

        public static bool IsSameAs(this OwnECPoint point1, ECPoint point2)
        {
            if (point1.IsInfinity || point2.IsInfinity)
            {
                return point1.IsInfinity && point2.IsInfinity;
            }
            point2 = point2.Normalize();
            return point1.X.ToBigInteger() == point2.X && point1.Y.ToBigInteger() == point2.Y;
        }

        public static bool IsSameAs(this IList<byte> data1, IList<byte> data2)
        {
            if (data1.Count != data2.Count)
            {
                return false;
            }
            return Enumerable.Range(0, data1.Count).All(idx => data1[idx] == data2[idx]);
        }
    }

    public static class BigIntegerExtensions
    {
        private const double LG_E = 1.4426950408889634073599246810019;  // = log2(e) = 1/ln(2).

        public static BigInteger Sqrt(this BigInteger n)
        {
            if (n.IsZero)
            {
                return BigInteger.Zero;
            }
            if (n.Sign < 0)
            {
                throw new ArithmeticException("NaN");
            }
            double lg = BigInteger.Log(n) * (LG_E * 0.5);
            int shift = (int)Math.Floor(lg) - 53;
            double exp = Math.Pow(2, lg - shift);
            BigInteger root = shift <= 0 ? new BigInteger(exp) << shift : new BigInteger(exp) >> (-shift);
            while (!isSqrt(n, root))
            {
                root = (root + n / root) >> 1;
            }
            return root;
        }

        private static Boolean isSqrt(BigInteger n, BigInteger root)
        {
            BigInteger difference = n - root * root;
            return difference.Sign >= 0 && difference <= (root << 1);
        }
    }
}
