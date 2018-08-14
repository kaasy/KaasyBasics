using System;
using System.Collections.Generic;
using System.Linq;
using System.Numerics;
using System.Text;
using System.Threading;
using System.Threading.Tasks;
using System.Windows.Forms;

namespace Utilities
{
    public struct FastInteger
    {
        public static readonly FastInteger One = new FastInteger(1L);
        public static readonly FastInteger Zero = new FastInteger(0L);

        private int[] bits;
        //store just 31 bits per int as unsigned (rule)
        //Only this top highest integer should have negative value if the full number is negative.
        public bool IsNegative { get { return !this.IsZero && this.bits[this.bits.Length - 1] < 0; } }
        public bool IsEven { get { return this.IsZero || (this.bits[0] & 1) == 0; } }
        public bool IsZero { get { return this.bits.Length == 0; } }

        public int[] Bits { get { return this.bits; } }

        public static FastInteger Abs(FastInteger number)
        {
            return number.IsNegative ? -number : number;
        }

        private const int LowMask = int.MaxValue;
        public int Sign { get { return this.IsZero ? 0 : (this.IsNegative ? -1 : 1); } }

        public byte LowestByte { get { return this.IsZero ? (byte)0 : (byte)this.bits[0]; } }

        private int Extender { get { return this.IsNegative ? LowMask : 0; } }
        private int IntExtender { get { return this.IsNegative ? -1 : 0; } }

        public FastInteger(ulong[] data) : this()
        {
            int bitsCount = data.Length * 64;
            int[] bits = new int[(bitsCount + 30) / 31];
            for (int i = 0, k = 0; i < bitsCount; i += 31, k++)
            {
                int cell = i >> 6;
                int bpos = i & 63;
                int value = (int)(data[cell] >> bpos);
                int need = bpos - 33;
                if (need > 0 && cell + 1 < data.Length)
                {
                    value |= (int)data[cell + 1] << (31 - need);
                }
                bits[k] = value & LowMask;
            }
            initialize(bits);
        }

        public FastInteger(byte[] data) : this()
        {
            int bitsCount = data.Length * 8;
            int[] bits = new int[(bitsCount + 30) / 31];
            int destIndex = 0, low = 0;
            for (int i = 0; i < data.Length; i++)
            {
                if (low >= 31)
                {
                    destIndex++;
                    low -= 31;
                }
                int unset = low - 23;
                bits[destIndex] |= ((int)data[i] << low) & LowMask;
                if (unset > 0)
                {
                    bits[destIndex + 1] |= (int)data[i] >> (8 - unset);
                }
                low += 8;
            }
            if (data.Length > 0 && (data[data.Length - 1] & 128) != 0)
            {
                bits[destIndex] |= -1 << low;
            }
            initialize(bits);
        }

        public FastInteger(int[] bits) : this()
        {
            initialize(bits);
        }

        public override bool Equals(object obj)
        {
            return obj is FastInteger && (FastInteger)obj == this;
        }

        public override int GetHashCode()
        {
            int hashCode = 1640531527;
            for (int i = 0; i < this.bits.Length; i++)
            {
                hashCode = (hashCode + this.bits[i]) * -1640531527;
            }
            return hashCode;
        }

        private void initialize(int[] bits)
        {
            int digits = bits.Length;
            if (digits > 1 && bits[digits - 1] == -1)
            {
                digits--;
                while (digits > 1 && bits[digits - 1] == LowMask)
                {
                    digits--;
                }
                bits[digits - 1] |= 1 << 31;
            }
            while (digits > 0 && bits[digits - 1] == 0)
            {
                digits--;
            }
            if (digits != bits.Length)
            {
                Array.Resize(ref bits, digits);
            }
            this.bits = bits;
        }

        public static implicit operator FastInteger(long x)
        {
            return new FastInteger(x);
        }

        public static implicit operator FastInteger(Bitsy.Core.FastInteger fi)
        {
            return new FastInteger(fi.Bits);
        }

        public static implicit operator Bitsy.Core.FastInteger(FastInteger fi)
        {
            return new Bitsy.Core.FastInteger(fi.bits);
        }

        public FastInteger(long number) : this()
        {
            if (number == 0)
            {
                this.bits = new int[0];
                return;
            }
            int digits = 1;
            for (long temp = number; temp < int.MinValue || temp > int.MaxValue; temp >>= 31)
            {
                digits++;
            }
            var bits = new int[digits];
            for (int i = 0; i + 1 < digits; i++)
            {
                bits[i] = (int)number & LowMask;
                number >>= 31;
            }
            bits[digits - 1] = (int)number;
            this.bits = bits;
        }

        public int Digits
        {
            get
            {
                int result = this.bits.Length;
                return result;
            }
        }

        private static int bsr(ulong x)
        {
            uint hi = (uint)(x >> 32);
            return hi != 0 ? bsr(hi) + 32 : bsr((uint)x);
        }

        private static int bsr(uint x)
        {
            if (x == 0)
            {
                return -1;
            }
            int result = 0;
            for (int position = 16; position > 0; position >>= 1)
            {
                uint y = x >> position;
                if (y != 0)
                {
                    result += position;
                    x = y;
                }
            }
            return result;
        }

        public int BitsCount
        {
            get
            {
                if (this.IsZero)
                {
                    return 0;
                }
                int top = this.bits[this.Digits - 1];
                int bitsCount = (this.Digits - 1) * 31 + 1 + bsr((uint)Math.Abs(top));
                return bitsCount;
            }
        }

        public static FastInteger operator +(FastInteger number1, FastInteger number2)
        {
            int digits1 = number1.Digits;
            int digits2 = number2.Digits;
            if (digits1 < digits2)
            {
                int temp1 = digits1; digits1 = digits2; digits2 = temp1;
                FastInteger temp2 = number1; number1 = number2; number2 = temp2;
            }
            if (number2.IsZero)
            {
                return number1;
            }
            int bitcount = 1 + Math.Max(number1.BitsCount, number2.BitsCount);
            int[] result = new int[(bitcount + 30) / 31];
            int carry = 0;
            for (int i = 0; i < digits2; i++)
            {
                carry += (number1.bits[i] & LowMask) + (number2.bits[i] & LowMask);
                result[i] = carry & LowMask;
                carry = (carry >> 31) & 1;
            }
            for (int i = digits2; i < digits1; i++)
            {
                carry += (number1.bits[i] & LowMask) + number2.Extender;
                result[i] = carry & LowMask;
                carry = (carry >> 31) & 1;
            }
            carry += (number1.IsNegative ? -1 : 0) + (number2.IsNegative ? -1 : 0);
            if (carry == -1)
            {
                result[digits1 - 1] |= -1 << 31;
            }
            if (digits1 < result.Length)
            {
                result[digits1] = carry;
            }
            return new FastInteger(result);
        }

        public static FastInteger operator -(FastInteger number1, FastInteger number2)
        {
            int bitcount = 1 + Math.Max(number1.BitsCount, number2.BitsCount);
            int[] result = new int[(bitcount + 30) / 31];
            int digits1 = number1.Digits;
            int digits2 = number2.Digits;
            int min = Math.Min(digits1, digits2);
            int max = Math.Max(digits1, digits2);

            int borrow = 0;
            for (int i = 0; i < min; i++)
            {
                borrow += (number1.bits[i] & LowMask) - (number2.bits[i] & LowMask);
                result[i] = borrow & LowMask;
                borrow >>= 31;
            }
            for (int i = min; i < digits1; i++)
            {
                borrow += (number1.bits[i] & LowMask) - number2.Extender;
                result[i] = borrow & LowMask;
                borrow >>= 31;
            }
            for (int i = min; i < digits2; i++)
            {
                borrow += number1.Extender - (number2.bits[i] & LowMask);
                result[i] = borrow & LowMask;
                borrow >>= 31;
            }
            borrow += (number1.IsNegative ? -1 : 0) + (number2.IsNegative ? 1 : 0);
            if (borrow == -1)
            {
                result[max - 1] |= -1 << 31;
            }
            if (max < result.Length)
            {
                result[max] = borrow;
            }
            return new FastInteger(result);
        }

        public static FastInteger operator -(FastInteger number)
        {
            int[] result = new int[number.bits.Length];
            number.bits.CopyTo(result, 0);
            int carry = 1;
            for (int i = 0; i < result.Length; i++)
            {
                carry += result[i] ^ LowMask;
                result[i] = carry & LowMask;
                carry = (carry >> 31) & 1;
            }
            result[result.Length - 1] |= (carry ^ 1) << 31;
            return new FastInteger(result);
        }

        private static int multiplyAndAdd(int[] source, int multiplier, int[] destination, int destinationIndex)
        {
            int carry = 0;
            for (int i = 0; i < source.Length; i++)
            {
                long multiplication = Math.BigMul(source[i], multiplier) + carry + destination[destinationIndex + i];
                destination[destinationIndex + i] = (int)multiplication & LowMask;
                carry = (int)(multiplication >> 31);
            }
            return carry;
        }

        public static FastInteger operator *(FastInteger number1, FastInteger number2)
        {
            int digits1 = number1.Digits;
            int digits2 = number2.Digits;
            int[] result = new int[digits1 + digits2];

            for (int i = 0; i < digits1; i++)
            {
                result[digits2 + i] = multiplyAndAdd(number2.bits, number1.bits[i], result, i);
            }
            return new FastInteger(result);
        }

        public static FastInteger operator >>(FastInteger x, int shift)
        {
            if (shift < 0)
            {
                return x << (-shift);
            }
            int blocks = shift / 31;
            if (blocks >= x.Digits)
            {
                return new FastInteger(new int[] { x.IntExtender });
            }
            int[] result = new int[x.Digits - blocks];
            for (int i = x.Digits - blocks; --i >= 0;)
            {
                result[i] = x.bits[i + blocks];
            }
            shift %= 31;
            if (shift != 0)
            {
                int backShift = 31 - shift;
                int carry = 0;
                for (int i = result.Length; --i >= 0;)
                {
                    int current = result[i];
                    result[i] = carry | (current >> shift);
                    carry = (current << backShift) & LowMask;
                }
            }
            return new FastInteger(result);
        }

        public static FastInteger operator <<(FastInteger x, int shift)
        {
            if (shift < 0)
            {
                return x >> (-shift);
            }
            int[] result = new int[x.Digits + (shift + 30) / 31];
            int blocks = shift / 31;
            for (int i = x.bits.Length; --i >= 0;)
            {
                result[i + blocks] = x.bits[i];
            }
            shift %= 31;
            if (shift != 0)
            {
                int backShift = 31 - shift;
                int previous = 0;
                for (int i = blocks; i < result.Length; i++)
                {
                    int current = result[i];
                    result[i] = ((current << shift) & LowMask) | previous;
                    previous = current >> backShift;
                }
            }
            return new FastInteger(result);
        }

        public static FastInteger operator &(FastInteger x, FastInteger y)
        {
            if (x.Digits < y.Digits)
            {
                FastInteger temp = x; x = y; y = temp;
            }
            if (y.IsZero)
            {
                return y;
            }
            int[] result = new int[x.Digits];
            for (int i = y.Digits; --i >= 0;)
            {
                result[i] = (x.bits[i] & y.bits[i]) & LowMask;
            }
            for (int i = y.Digits; i < x.Digits; i++)
            {
                result[i] = (x.bits[i] & y.Extender) & LowMask;
            }
            result[result.Length - 1] |= (x.Extender & y.Extender) << 31;
            return new FastInteger(result);
        }

        public static FastInteger operator ^(FastInteger x, FastInteger y)
        {
            if (x.Digits < y.Digits)
            {
                FastInteger temp = x; x = y; y = temp;
            }
            if (y.IsZero)
            {
                return x;
            }
            int[] result = new int[x.Digits];
            for (int i = y.Digits; --i >= 0;)
            {
                result[i] = (x.bits[i] ^ y.bits[i]) & LowMask;
            }
            for (int i = y.Digits; i < x.Digits; i++)
            {
                result[i] = (x.bits[i] ^ y.Extender) & LowMask;
            }
            result[result.Length - 1] |= (x.Extender ^ y.Extender) << 31;
            return new FastInteger(result);
        }

        public static FastInteger operator *(FastInteger number, int scalar)
        {
            int digits = number.Digits;
            int[] result = new int[digits + 1];
            result[digits] = multiplyAndAdd(number.bits, scalar, result, 0);
            return new FastInteger(result);
        }

        public static bool operator ==(FastInteger a, FastInteger b)
        {
            int digits = a.Digits;
            if (digits != b.Digits)
            {
                return false;
            }
            for (int i = 0; i < digits; i++)
            {
                if (a.bits[i] != b.bits[i])
                {
                    return false;
                }
            }
            return true;
        }

        public static bool operator !=(FastInteger a, FastInteger b)
        {
            return !(a == b);
        }

        /// <summary>
        /// return this / n modulo p.
        /// </summary>
        /// <param name="n"></param>
        /// <param name="p"></param>
        /// <returns></returns>
        public FastInteger DivideModulo(FastInteger n, FastInteger p)
        {
            if (n.IsZero)
            {
                return n;
            }
            if (p.IsZero)
            {
                return p;
            }
            FastInteger a = this;   //x * a + u1 * b = n
            FastInteger c = 0;      //y * c + u2 * d = p
            FastInteger x = n;
            while (p.IsEven)
            {
                p >>= 1;
            }
            FastInteger y = p;
            while (x.IsEven)
            {
                x >>= 1;
                a = a.IsEven ? a >> 1 : (a + p) >> 1;
            }
            while (!x.IsZero)
            {
                while (x.IsEven)
                {
                    x >>= 1;
                    a = a.IsEven ? a >> 1 : (a + p) >> 1;
                }
                FastInteger diff = x - y;
                if (diff.IsNegative)
                {
                    y = x;
                    x = -diff;
                    FastInteger temp = a; a = c; c = temp;
                }
                else
                {
                    x = diff;
                }
                a -= c;
                if (a.IsNegative)
                {
                    a += p;
                }
            }
            return c;
        }

        public override string ToString()
        {
            FastInteger x = this;
            bool negative = x.IsNegative;
            if (negative)
            {
                x = -x;
            }
            StringBuilder sb = new StringBuilder();
            int[] bits = new int[x.bits.Length];
            x.bits.CopyTo(bits, 0);

            do
            {
                long mod = 0;
                for (int i = bits.Length; --i >= 0;)
                {
                    mod = (mod << 31) + bits[i];
                    bits[i] = (int)Math.DivRem(mod, 10L, out mod);
                }
                sb.Append((char)('0' + mod));
            } while (bits.Any(b => b != 0));
            if (negative)
            {
                sb.Append('-');
            }
            for (int i = 0, j = sb.Length - 1; i < j; i++, j--)
            {
                char c = sb[i];
                sb[i] = sb[j];
                sb[j] = c;
            }
            return sb.ToString();
        }

        public FastInteger ModInverse(FastInteger p)
        {
            return One.DivideModulo(this, p);
        }

        public byte[] ToByteArray()
        {
            int totalBits = this.Digits * 31 + 1;
            byte[] result = new byte[(totalBits + 7) >> 3];
            int destPosition = 0, low = 0;
            for (int i = 0; i < this.Digits; i++)
            {
                int digit = this.bits[i];
                result[destPosition++] |= (byte)(digit << low);
                int free = low + 23;
                digit >>= 8 - low;
                while (free >= 8)
                {
                    result[destPosition++] |= (byte)digit;
                    free -= 8;
                    digit >>= 8;
                }
                if (free > 0)
                {
                    result[destPosition] = (byte)digit;
                }
                low = (low + 31) & 7;
            }
            if (this.IsNegative)
            {
                for (int k = destPosition; k < result.Length; k++)
                {
                    result[k] |= (byte)(-1 << low);
                    low = 0;
                }
            }
            return result;
        }

        public static bool operator <=(FastInteger a, FastInteger b)
        {
            var diff = a - b;
            return diff.IsNegative || diff.IsZero;
        }
        public static bool operator <(FastInteger a, FastInteger b)
        {
            return (a - b).IsNegative;
        }
        public static bool operator >=(FastInteger a, FastInteger b)
        {
            return !(a < b);
        }
        public static bool operator >(FastInteger a, FastInteger b)
        {
            return !(a <= b);
        }

        public bool GetBitAt(int index)
        {
            int mod, div = Math.DivRem(index, 31, out mod);
            return (this.bits[div] & (1 << mod)) != 0;
        }

        public FastInteger ModPow(FastInteger power, FastInteger prime, Func<FastInteger, FastInteger> modPFunction)
        {
            FastInteger a = this;
            if (power.IsNegative)
            {
                a = a.ModInverse(prime);
                power = -power;
            }
            int highestBitPosition = power.Digits * 31 - 1;
            while (highestBitPosition >= 0 && !power.GetBitAt(highestBitPosition))
            {
                highestBitPosition--;
            }
            FastInteger result = highestBitPosition < 0 ? One : a;
            for (int i = highestBitPosition; --i >= 0;)
            {
                result = modPFunction(result.Square());
                if (power.GetBitAt(i))
                {
                    result = modPFunction(result * a);
                }
            }
            var diff = result - prime;
            return diff.IsNegative ? result : diff;
        }

        public FastInteger ShanksSqrt(FastInteger p, Func<FastInteger, FastInteger> modPFunction)
        {
            if (this.ModPow((p - 1) >> 1, p, modPFunction) == (p - 1))
            {
                return -1;
            }   //No Sqrt Exists
            if ((p.LowestByte & 3) == 3)
            {
                var result = this.ModPow((p + 1) >> 2, p, modPFunction);
                return result;
            }
            throw new NotImplementedException("Not a special prime modulo.");
        }
    }

    public static class FastIntegerUnitTest
    {
        public static bool UnitTest()
        {
            //    SPECIAL prime: 2^64 - 2^32 + 1
            //prime: 2^64 - 2^34 + 1
            //prime: 2^64 - 2^40 + 1
            //    SPECIAL prime: 2 ^ 128 - 2 ^ 54 + 1
            //prime: 2 ^ 128 - 2 ^ 108 + 1
            //prime: 2 ^ 256 - 2 ^ 168 + 1
            //prime: 2 ^ 256 - 2 ^ 174 + 1
            //prime: 2 ^ 512 - 2 ^ 32 + 1
            //prime: 2 ^ 512 - 2 ^ 288 + 1
            //    SPECIAL prime: 2 ^ 1024 - 2 ^ 142 + 1
            //    SPECIAL prime: 2 ^ 1024 - 2 ^ 226 + 1
            //prime: 2 ^ 1024 - 2 ^ 562 + 1
            //prime: 2 ^ 1024 - 2 ^ 718 + 1
            //prime: 2 ^ 1024 - 2 ^ 856 + 1
            //prime: 2 ^ 1024 - 2 ^ 880 + 1
            //prime: 2 ^ 4096 - 2 ^ 3510 + 1
            //prime: 2 ^ 4096 - 2 ^ 3708 + 1
            if ("not exec".Trim() == "exec")
            {
                StringBuilder primes = new StringBuilder();
                for (int exp = 6; exp <= 12; exp++)
                {
                    int shift = 1 << exp;
                    Parallel.For(32, shift, i =>
                    {
                        ulong[] baseValue = new ulong[shift / 64];
                        ulong[] expValue = new ulong[shift / 64];
                        ulong[] modulo = new ulong[shift / 64];
                        ulong[] subtractor = new ulong[shift / 64];
                        ulong[] one = new ulong[shift / 64];
                        ulong[] pPlusOne = new ulong[shift / 64];
                        subtractor[i / 64] |= 1UL << (i & 63);
                        modulo[0] = 1UL;
                        one[0] = 1UL;
                        AsmX64Operations.Subtract(modulo, subtractor, modulo, 0, modulo.Length);
                        modulo.CopyTo(expValue, 0);
                        modulo.CopyTo(pPlusOne, 0);
                        expValue[0]--;
                        pPlusOne[0]++;

                        int isPrime = 1;
                        for (int a = 2; a <= 128; a++)
                        {
                            one.CopyTo(baseValue, 0);
                            baseValue[0] = (ulong)a;
                            AsmX64Operations.ModularExponentiation(baseValue, expValue, modulo, modulo.Length, 5);
                            if (Enumerable.Range(0, modulo.Length).All(idx => baseValue[idx] == one[idx]) ||
                                Enumerable.Range(0, modulo.Length).All(idx => baseValue[idx] == pPlusOne[idx]))
                            {
                                continue;
                            }
                            isPrime = 0;
                            break;
                        }

                        if (isPrime != 0)
                        {
                            lock (primes)
                            {
                                primes.AppendLine("prime: 2^" + shift.ToString() + " - 2^" + i.ToString() + " + 1");
                            }
                        }
                    });
                }
                Clipboard.SetText(primes.ToString());
                System.Windows.Forms.MessageBox.Show(primes.ToString(), "message");
            }

            Random random = new Random(1001);
            bool ok = true;
            for (int i = 400; --i >= 0;)
            {
                byte[] bytes1 = new byte[31 + random.Next(2)];
                byte[] bytes2 = new byte[31 + random.Next(2)];
                random.NextBytes(bytes1);
                random.NextBytes(bytes2);

                BigInteger n1 = new BigInteger(bytes1);
                BigInteger n2 = new BigInteger(bytes2);
                FastInteger f1 = new FastInteger(bytes1);
                FastInteger f2 = new FastInteger(bytes2);
                if (n1.ToString() != f1.ToString())
                {
                    ok = false;
                }
                if (n2.ToString() != f2.ToString())
                {
                    ok = false;
                }
                BigInteger a1 = n1 + n2;
                FastInteger a2 = f1 + f2;
                if (a1.ToString() != a2.ToString())
                {
                    ok = false;
                }
                BigInteger s1 = n1 - n2;
                FastInteger s2 = f1 - f2;
                if (s1.ToString() != s2.ToString())
                {
                    ok = false;
                }
                BigInteger m1 = n1 * n2;
                FastInteger m2 = f1 * f2;
                if (m1.ToString() != m2.ToString())
                {
                    ok = false;
                }
                int shrvalue = random.Next(256) + 1;
                BigInteger sh1 = n1 >> shrvalue;
                FastInteger sh2 = f1 >> shrvalue;
                if (sh1.ToString() != sh2.ToString())
                {
                    ok = false;
                }
                if ((-f1).ToString() != (-n1).ToString() || (-f2).ToString() != (-n2).ToString())
                {
                    ok = false;
                }
                int shlvalue = random.Next(256) + 1;
                BigInteger shl1 = n1 << shlvalue;
                FastInteger shl2 = f1 << shlvalue;
                if (shl1.ToString() != shl2.ToString())
                {
                    ok = false;
                }
                BigInteger and1 = n1 & n2;
                FastInteger and2 = f1 & f2;
                if (and1.ToString() != and2.ToString())
                {
                    ok = false;
                }
                BigInteger xor1 = n1 ^ n2;
                FastInteger xor2 = f1 ^ f2;
                if (xor1.ToString() != xor2.ToString())
                {
                    ok = false;
                }

                if (!ok)
                {
                    break;
                }
            }
            return ok;
        }
    }
}
