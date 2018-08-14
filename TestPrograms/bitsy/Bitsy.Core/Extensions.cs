using System;
using System.Collections.Generic;
using System.Linq;
using System.Numerics;

namespace Bitsy.Core
{
    public static class Extensions
    {
        public static BigInteger ModInverse(this BigInteger n, BigInteger p)
        {
            BigInteger x = 1;   //x * n + u1 * p = a
            BigInteger y = 0;   //y * n + u2 * p = b
            BigInteger a = n;
            BigInteger b = p;
            while (b != 0)
            {
                BigInteger q = BigInteger.DivRem(a, b, out a);
                x -= q * y;
                BigInteger temp;
                temp = a; a = b; b = temp;
                temp = x; x = y; y = temp;
            }
            return x < 0 ? x + p : x;
        }

        public static bool TestBit(this BigInteger i, int n)
        {
            int bitLength = i.BitLength();
            return !(i >> n).IsEven;
        }

        public static int BitLength(this BigInteger i)
        {
            int bitLength = 0;
            do
            {
                bitLength++;
            } while ((i >>= 1) != 0);
            return bitLength;
        }

        public static byte[] ToByteArrayUnsigned(this BigInteger i, bool bigEndian)
        {
            byte[] bytes = i.ToByteArray();
            if (bytes[bytes.Length - 1] == 0x00)
            {
                Array.Resize(ref bytes, bytes.Length - 1);
            }
            if (bigEndian)
            {
                Array.Reverse(bytes, 0, bytes.Length);
            }

            return bytes;
        }

        public static byte[] ToByteArrayUnsigned(this FastInteger i, bool bigEndian)
        {
            byte[] bytes = i.ToByteArray();
            int length = bytes.Length;
            while (length > 1 && ((bytes[length - 1] == 0 && bytes[length - 2] < 128) || (bytes[length - 1] == 255 && bytes[length - 2] >= 128)))
            {
                length--;
            }
            Array.Resize(ref bytes, length);
            if (bigEndian)
            {
                Array.Reverse(bytes, 0, bytes.Length);
            }
            return bytes;
        }

        private static void FindSE(BigInteger p, out BigInteger s, out BigInteger e)
        {
            s = p - 1;
            e = 0;
            while (s.IsEven)
            {
                s >>= 1;
                e += 1;
            }
        }


        public static BigInteger Order(this BigInteger b, BigInteger p)
        {
            BigInteger m = 1;
            BigInteger e = 0;

            while (BigInteger.ModPow(b, m, p) != 1)
            {
                m *= 2;
                e++;
            }

            return e;
        }

        private static BigInteger TwoExp(BigInteger e)
        {
            BigInteger a = 1;

            while (e > 0)
            {
                a *= 2;
                e--;
            }

            return a;
        }

        public static string ToHex(this BigInteger b)
        {
            return Hex.BigIntegerToHex(b);
        }

        public static string ToHex(this FastInteger b)
        {
            return Hex.FastIntegerToHex(b);
        }

        public static string ToHex(this byte[] bytes)
        {
            return Hex.BytesToHex(bytes);
        }

        public static BigInteger HexToBigInteger(this string s)
        {
            return Hex.HexToBigInteger(s);
        }

        public static FastInteger HexToFastInteger(this string s)
        {
            return Hex.HexToFastInteger(s);
        }

        public static byte[] HexToBytes(this string s)
        {
            return Hex.HexToBytes(s);
        }

        public static BigInteger ToBigInteger(this IEnumerable<byte> bytes, bool bigEndian)
        {
            byte[] clone = bytes.ToArray();
            if (bigEndian)
            {
                Array.Reverse(clone);
            }
            return new BigInteger(clone);
        }

        public static BigInteger ToBigIntegerUnsigned(this IEnumerable<byte> bytes, bool bigEndian)
        {
            byte[] clone = bigEndian ?
                Enumerable.Repeat((byte)0, 1).Concat(bytes).ToArray() :
                bytes.Concat(Enumerable.Repeat((byte)0, 1)).ToArray();
            if (bigEndian)
            {
                Array.Reverse(clone);
            }
            return new BigInteger(clone);
        }

        public static FastInteger ToFastIntegerUnsigned(this IEnumerable<byte> bytes, bool bigEndian)
        {
            byte[] clone = bigEndian ?
                Enumerable.Repeat((byte)0, 1).Concat(bytes).ToArray() :
                bytes.Concat(Enumerable.Repeat((byte)0, 1)).ToArray();
            if (bigEndian)
            {
                Array.Reverse(clone);
            }
            return new FastInteger(clone);
        }

        public static BigInteger Square(this BigInteger x)
        {
            return x * x;
        }

        public static FastInteger Square(this FastInteger x)
        {
            return x * x;
        }

        public static FastInteger ToFastInteger(this IEnumerable<byte> bytes)
        {
            return new FastInteger(bytes.Concat(Enumerable.Repeat((byte)0, 1)).ToArray());
        }

        public static BigInteger ToBigInteger(this ulong[] bits)
        {
            BigInteger result = new BigInteger(bits.SelectMany(l => BitConverter.GetBytes(l)).Concat(Enumerable.Repeat((byte)0, 1)).ToArray());
            return result;
        }

        public static FastInteger ToFastInteger(this ulong[] bits)
        {
            FastInteger result = new FastInteger(bits.SelectMany(l => BitConverter.GetBytes(l)).Concat(Enumerable.Repeat((byte)0, 1)).ToArray());
            return result;
        }

        public static string GetBitcoinAddress(this ECPoint publicKey, bool compressed = true)
        {
            var pubKeyHash = Hash160.Hash(publicKey.EncodePoint(compressed));

            byte[] addressBytes = new byte[pubKeyHash.Length + 1];
            Buffer.BlockCopy(pubKeyHash, 0, addressBytes, 1, pubKeyHash.Length);
            return Base58.EncodeWithCheckSum(addressBytes);
        }

        public static BigInteger ShanksSqrt(this BigInteger a, BigInteger p)
        {
            if (BigInteger.ModPow(a, (p - 1) / 2, p) == (p - 1))
            {
                return -1;
            }   //No Sqrt Exists

            if ((int)(p & 3) == 3)
            {
                return BigInteger.ModPow(a, (p + 1) / 4, p);
            }

            //Initialize 
            BigInteger s, e;
            FindSE(p, out s, out e);

            BigInteger n, m, x, b, g, r;
            n = 2;
            while (BigInteger.ModPow(n, (p - 1) / 2, p) == 1)
            {
                n++;
            }//Finds Generator

            x = BigInteger.ModPow(a, (s + 1) / 2, p);
            b = BigInteger.ModPow(a, s, p);
            g = BigInteger.ModPow(n, s, p);
            r = e;
            m = b.Order(p);
            if (m == 0)
            {
                return x;
            }

            //For Debugging
            //Console.WriteLine("{0}, {1}, {2}, {3}, {4}",m, x, b, g, r);
            while (m > 0)
            {
                x = (x * BigInteger.ModPow(g, TwoExp(r - m - 1), p)) % p;
                b = (b * BigInteger.ModPow(g, TwoExp(r - m), p)) % p;
                g = BigInteger.ModPow(g, TwoExp(r - m), p);
                r = m;
                m = b.Order(p);
                //For Debugging
                //Cnsole.WriteLine("{0}, {1}, {2}, {3}, {4}", m, x, b, g, r);
            }

            return x;
        }
    }
}
