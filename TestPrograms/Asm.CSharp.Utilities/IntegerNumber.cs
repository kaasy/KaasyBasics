using System;
using System.Collections.Generic;
using System.Linq;
using System.Numerics;
using System.Runtime.InteropServices;
using System.Text;
using System.Threading;
using System.Threading.Tasks;
using System.Windows.Forms;

namespace Utilities
{
    public struct IntegerNumber
    {
        private const ulong SignBit = 1UL << 63;
        public static readonly IntegerNumber One = new IntegerNumber(1L);
        public static readonly IntegerNumber MinusOne = new IntegerNumber(-1L);
        public static readonly IntegerNumber Zero = new IntegerNumber(0L);
        private const double LG_E = 1.4426950408889634073599246810019;  // = log2(e) = 1/ln(2).
        private static readonly double LOG10_2 = Math.Log10(2);

        private ulong[] bits;

        public IntegerNumber(ulong[] data, bool copy) : this()
        {
            if (copy)
            {
                this.bits = new ulong[data.Length];
                data.CopyTo(this.bits, 0);
                return;
            }
            this.bits = data;
        }

        public IntegerNumber(long number) : this()
        {
            this.bits = new ulong[1];
            this.bits[0] = (ulong)number;
        }

        public IntegerNumber(byte[] data) : this()
        {
            int digits = (data.Length + 7) >> 3;
            this.bits = new ulong[digits];
            for (int i = data.Length >> 3; --i >= 0;)
            {
                this.bits[i] = BitConverter.ToUInt64(data, i << 3);
            }
            int lastBytes = data.Length & 7;
            if (lastBytes == 0)
            {
                return;
            }
            ulong lastDigit = 0UL;
            for (int i = lastBytes, j = data.Length - 1; --i >= 0; j--)
            {
                lastDigit |= ((ulong)data[j]) << (i << 3);
            }
            if (data[data.Length - 1] >= 128)
            {
                lastDigit |= ulong.MaxValue << (lastBytes << 3);
            }
            this.bits[digits - 1] = lastDigit;
        }

        public long ToLong { get { return (long)this.bits[0]; } }

        public ulong Extender { get { return GetExtender(this.bits[this.bits.Length - 1]); } }
        public int Digits { get { return this.bits.Length; } }

        public int RealDigits
        {
            get
            {
                int length = this.Digits;
                while (length > 1 && (long)this.bits[length - 1] == ((long)this.bits[length - 2]) >> 63)
                {
                    length--;
                }
                return length;
            }
        }

        public int BitsCount
        {
            get
            {
                int digits = this.RealDigits;
                long highestDigit = (long)this.bits[digits - 1];
                highestDigit = highestDigit < 0 ? -highestDigit : highestDigit;
                return BSR((ulong)highestDigit) + 1 + (digits - 1) * 64;
            }
        }

        private long HighestDigit
        {
            get { return (long)this.bits[this.Digits - 1]; }
            set { this.bits[this.Digits - 1] = (ulong)value; }
        }

        public bool IsNegative { get { return this.HighestDigit < 0; } }
        public bool IsEven { get { return (this.bits[0] & 1) == 0; } }
        public bool IsZero { get { return isZero(0, this.Digits); } }

        private bool isZero(int startIndex, int count)
        {
            for (int i = startIndex, j = startIndex + count - 1; i <= j; i++, j--)
            {
                if (this.bits[i] != 0UL || this.bits[j] != 0UL)
                {
                    return false;
                }
            }
            return true;
        }

        public bool HasBit(int bitIndex)
        {
            int cell = bitIndex >> 6;
            return cell >= 0 && cell < this.bits.Length && (this.bits[cell] & (1UL << (bitIndex & 63))) != 0;
        }

        public static IntegerNumber Abs(IntegerNumber number)
        {
            return number.IsNegative ? -number : number;
        }

        public int Sign { get { return this.IsNegative ? -1 : (this.IsZero ? 0 : 1); } }
        public byte LowestByte { get { return (byte)this.bits[0]; } }

        public override bool Equals(object obj)
        {
            return obj is IntegerNumber && (IntegerNumber)obj == this;
        }

        public override int GetHashCode()
        {
            int hashCode = 1640531527;
            for (int i = 0; i < this.bits.Length; i++)
            {
                hashCode = (hashCode + this.bits[i].GetHashCode()) * -1640531527;
            }
            return hashCode;
        }

        public static implicit operator IntegerNumber(long x)
        {
            return new IntegerNumber(x);
        }

        public static int BSR(ulong x)
        {
            uint hi = (uint)(x >> 32);
            return hi != 0 ? BSR(hi) + 32 : BSR((uint)x);
        }

        public static int BSR(uint x)
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

        public IntegerNumber Square()
        {
            ulong[] result = new ulong[this.Digits * 2];
            AsmX64Operations.FastestSquare(this.bits, result, this.Digits, true);
            return new IntegerNumber(result, false);
        }

        private IntegerNumber ExtendTo(int digits)
        {
            if (digits <= this.Digits)
            {
                return this;
            }
            ulong[] value = new ulong[digits];
            this.bits.CopyTo(value, 0);
            ulong extender = this.Extender;
            for (int i = this.Digits; i < digits; i++)
            {
                value[i] = extender;
            }
            return new IntegerNumber(value, false);
        }

        public static IntegerNumber operator +(IntegerNumber number1, IntegerNumber number2)
        {
            if (number1.IsNegative != number2.IsNegative)
            {
                return number2.IsZero ? number1 : number1 - (-number2);
            }
            IntegerNumber source1 = number1.ExtendTo(number2.Digits);
            IntegerNumber source2 = number2.ExtendTo(number1.Digits);
            ulong[] result = new ulong[source1.Digits];
            byte carry = AsmX64Operations.Add(source1.bits, source2.bits, result, 0, result.Length);
            if (carry != (byte)(result[result.Length - 1] >> 63))
            {
                Array.Resize(ref result, result.Length + 1);
                result[result.Length - 1] = carry == 0 ? 0UL : ulong.MaxValue;
            }
            return new IntegerNumber(result, false);
        }

        public static IntegerNumber operator -(IntegerNumber number1, IntegerNumber number2)
        {
            if (number1.IsNegative != number2.IsNegative)
            {
                return number2.IsZero ? number1 : number1 + (-number2);
            }   //when subtracting numbers with the same sign, the absolute value of the result is less than the absoulute value of each input. (number size decreases).
            IntegerNumber source1 = number1.ExtendTo(number2.Digits);
            IntegerNumber source2 = number2.ExtendTo(number1.Digits);
            ulong[] result = new ulong[source1.Digits];
            AsmX64Operations.Subtract(source1.bits, source2.bits, result, 0, result.Length);
            return new IntegerNumber(result, false);
        }

        public static IntegerNumber operator -(IntegerNumber number)
        {
            if (number.bits[number.bits.Length - 1] == SignBit && number.isZero(0, number.bits.Length - 1))
            {   //make positive the highly compact negative numbers - as a special case.
                ulong[] bits = new ulong[number.bits.Length + 1];
                bits[number.bits.Length - 1] = SignBit;
                return new IntegerNumber(bits, false);
            }
            IntegerNumber result = new IntegerNumber(number.bits, true);
            AsmX64Operations.Negate(result.bits, result.Digits);
            return result;
        }

        public static IntegerNumber operator *(IntegerNumber number1, IntegerNumber number2)
        {
            IntegerNumber source1 = number1.ExtendTo(number2.Digits);
            IntegerNumber source2 = number2.ExtendTo(number1.Digits);
            ulong[] result = new ulong[source1.Digits * 2];
            AsmX64Operations.FastestMultiplication(source1.bits, source2.bits, result, source1.Digits, true);
            return new IntegerNumber(result, false);
        }

        public static IntegerNumber operator *(IntegerNumber number, long scalar)
        {
            ulong[] result = new ulong[number.Digits + 1];
            AsmX64Operations.MultiplyDigit(number.bits, (ulong)scalar, result, number.Digits, true);
            return new IntegerNumber(result, false);
        }

        public static IntegerNumber operator /(IntegerNumber number1, IntegerNumber number2)
        {
            IntegerNumber remainder, quotient = IntegerNumber.DivRem(number1, number2, out remainder);
            return quotient;
        }

        public static IntegerNumber operator /(IntegerNumber number1, long number2)
        {
            long remainder;
            IntegerNumber quotient = IntegerNumber.DivRem(number1, number2, out remainder);
            return quotient;
        }

        public static IntegerNumber operator %(IntegerNumber number1, IntegerNumber number2)
        {
            IntegerNumber remainder, quotient = IntegerNumber.DivRem(number1, number2, out remainder);
            return remainder;
        }

        public static IntegerNumber operator >>(IntegerNumber x, int shift)
        {
            if (shift == 0)
            {
                return x;
            }
            if (shift < 0)
            {
                return x << (-shift);
            }
            IntegerNumber result = new IntegerNumber(x.bits, true);
            AsmX64Operations.IntShr(result.bits, shift, result.Digits, true);
            return result;
        }

        public static IntegerNumber operator <<(IntegerNumber x, int shift)
        {
            if (shift == 0)
            {
                return x;
            }
            if (shift < 0)
            {
                return x >> (-shift);
            }

            int xlen = x.RealDigits;
            int resultDigits = xlen + (shift >> 6);
            int singleDigitShift = shift & 63;
            if (singleDigitShift != 0)
            {
                long hix = (long)x.bits[xlen - 1] >> (63 - singleDigitShift);
                if (hix != 0L && hix != -1L)
                {
                    resultDigits++;
                }
            }
            ulong[] result = x.bits;
            ExtendToExact(ref result, resultDigits, true);
            AsmX64Operations.IntShl(result, shift, resultDigits);
            return new IntegerNumber(result, false);
        }

        public static bool operator ==(IntegerNumber a, IntegerNumber b)
        {
            if (a.Digits < b.Digits)
            {
                return b == a;
            }
            //now a.Digits >= b.Digits
            for (int i = 0, j = b.Digits - 1; i <= j; i++, j--)
            {
                if (a.bits[i] != b.bits[i] || a.bits[j] != b.bits[j])
                {
                    return false;
                }
            }
            var extender = b.Extender;
            var startOffset = b.Digits;
            for (int i = a.Digits; --i >= startOffset;)
            {
                if (a.bits[i] != extender)
                {
                    return false;
                }
            }
            return true;
        }

        public static bool operator !=(IntegerNumber a, IntegerNumber b)
        {
            return !(a == b);
        }

        /// <summary>
        /// return this / n modulo p.
        /// </summary>
        /// <param name="n"></param>
        /// <param name="p"></param>
        /// <returns></returns>
        public IntegerNumber DivideModulo(IntegerNumber n, IntegerNumber p)
        {
            if (n.IsZero)
            {
                return n;
            }
            if (p.IsZero)
            {
                return p;
            }
            IntegerNumber a = this;   //x * a + u1 * b = n
            IntegerNumber c = 0;      //y * c + u2 * d = p
            IntegerNumber x = n;
            while (p.IsEven)
            {
                p >>= 1;
            }
            IntegerNumber y = p;
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
                IntegerNumber diff = x - y;
                if (diff.IsNegative)
                {
                    y = x;
                    x = -diff;
                    IntegerNumber temp = a; a = c; c = temp;
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

        private static void getStringRepresentation(IntegerNumber x, int base10Digits, StringBuilder result)
        {
            if (base10Digits <= 18)
            {
                string value = x.ToLong.ToString();
                int diff = value.Length - base10Digits;
                if (diff != 0)
                {
                    if (diff < 0)
                    {
                        value = "".PadLeft(-diff, '0') + value;
                    }
                    else
                    {
                        value = value.Remove(0, diff);
                    }
                }
                result.Append(value);
                return;
            }
            int secondHalfDigits = (base10Digits + 1) / 2;
            IntegerNumber half = IntegerNumber.Pow(10, secondHalfDigits);
            IntegerNumber reminder, quotient = IntegerNumber.DivRem(x, half, out reminder);
            getStringRepresentation(quotient, base10Digits - secondHalfDigits, result);
            getStringRepresentation(reminder, secondHalfDigits, result);
        }

        public override string ToString()
        {
            IntegerNumber x = this;
            bool negative = x.IsNegative;
            if (negative)
            {
                x = -x;
            }
            StringBuilder result = new StringBuilder();
            getStringRepresentation(x, Math.Max(1, (int)Math.Ceiling(x.BitsCount * LOG10_2)), result);
            int lowzeros = 0;
            while (lowzeros + 1 < result.Length && result[lowzeros] == '0')
            {
                lowzeros++;
            }
            if (lowzeros > 0)
            {
                result = result.Remove(0, lowzeros);
            }
            if (negative)
            {
                result = result.Insert(0, '-');
            }
            return result.ToString();
        }

        public IntegerNumber ModInverse(IntegerNumber p)
        {
            return One.DivideModulo(this, p);
        }

        public byte[] ToByteArray()
        {
            byte[] result = new byte[this.Digits * 8];
            GCHandle handle = GCHandle.Alloc(this.bits, GCHandleType.Pinned);
            try
            {
                Marshal.Copy(handle.AddrOfPinnedObject(), result, 0, result.Length);
            }
            finally
            {
                handle.Free();
            }
            return result;
        }

        public static bool operator <=(IntegerNumber a, IntegerNumber b)
        {
            var diff = a - b;
            return diff.IsNegative || diff.IsZero;
        }
        public static bool operator <(IntegerNumber a, IntegerNumber b)
        {
            return (a - b).IsNegative;
        }
        public static bool operator >=(IntegerNumber a, IntegerNumber b)
        {
            return !(a < b);
        }
        public static bool operator >(IntegerNumber a, IntegerNumber b)
        {
            return !(a <= b);
        }

        public bool GetBitAt(int index)
        {
            return (this.bits[index >> 6] & (1UL << (index & 63))) != 0;
        }

        public IntegerNumber ModPow(IntegerNumber power, IntegerNumber prime, Func<IntegerNumber, IntegerNumber> modPFunction)
        {
            IntegerNumber a = this;
            if (power.IsNegative)
            {
                a = a.ModInverse(prime);
                power = -power;
            }
            int highestBitPosition = power.Digits * 64 - 1;
            while (highestBitPosition >= 0 && !power.GetBitAt(highestBitPosition))
            {
                highestBitPosition--;
            }
            IntegerNumber result = highestBitPosition < 0 ? One : a;
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

        public IntegerNumber ShanksSqrt(IntegerNumber p, Func<IntegerNumber, IntegerNumber> modPFunction)
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

        public static IntegerNumber Pow(int baseNumber, int exponent)
        {
            int highestBitPosition = BSR((uint)exponent);
            if (highestBitPosition < 0)
            {
                return IntegerNumber.One;
            }
            IntegerNumber number = new IntegerNumber(baseNumber);
            IntegerNumber result = number;
            for (int i = highestBitPosition; --i >= 0;)
            {
                result = result.Square();
                shrink(ref result);
                if ((exponent & (1 << i)) != 0)
                {
                    result *= baseNumber;
                }
            }
            return result;
        }

        public static IntegerNumber DivRem(IntegerNumber number, long digit, out long remainder)
        {
            IntegerNumber result = new IntegerNumber(number.bits, true);
            remainder = (long)AsmX64Operations.DivideDigit(result.bits, (ulong)digit, result.Digits, true);
            return result;
        }

        private static ulong GetExtender(ulong highestDigit)
        {
            ulong extender = (ulong)(((long)highestDigit) >> 63);
            return extender;
        }

        private static void ExtendToExact(ref ulong[] number, int digits, bool makeNewCopyAlways)
        {
            int numberLength = number.Length;
            if (numberLength == digits && !makeNewCopyAlways)
            {
                return;
            }
            ulong[] result = new ulong[digits];
            for (int i = Math.Min(digits, numberLength); --i >= 0;)
            {
                result[i] = number[i];
            }
            if (numberLength < digits)
            {
                ulong extender = GetExtender(number[numberLength - 1]);
                for (int i = numberLength; i < digits; i++)
                {
                    result[i] = extender;
                }
            }
            number = result;
        }

        const int IntegerDigitsBreakPoint = 64;

        public static IntegerNumber DivRem(IntegerNumber dividend, IntegerNumber divisor, out IntegerNumber remainder)
        {
            shrink(ref dividend);
            shrink(ref divisor);
            int dividendDigits = dividend.Digits;
            int divisorDigits = divisor.Digits;
            IntegerNumber quotient;
            int n = Math.Max(divisorDigits, (dividendDigits + 1) >> 1);
            if (divisorDigits <= IntegerDigitsBreakPoint)
            {
                dividend = dividend.ExtendTo(n * 2);
                divisor = divisor.ExtendTo(n);
                quotient = new IntegerNumber(new ulong[n * 2], false);
                remainder = new IntegerNumber(dividend.bits, true);
                AsmX64Operations.GetDivMod(remainder.bits, divisor.bits, n, quotient.bits, true);
                shrink(ref quotient);
                shrink(ref remainder);
                return quotient;
            }
            bool sign = false;
            if (dividend.IsNegative != divisor.IsNegative)
            {
                dividend = IntegerNumber.Abs(dividend);
                divisor = IntegerNumber.Abs(divisor);
                sign = true;
            }
            else if (dividend.IsNegative && divisor.IsNegative)
            {
                dividend = -dividend;
                divisor = -divisor;
            }

            int delta = n - divisorDigits;
            IntegerNumber x = (divisor << (delta * 128)).Inverse();
            quotient = dividend * x >> (n * 128);
            shrink(ref quotient);
            remainder = dividend - quotient * divisor;
            shrink(ref remainder);
            int count = 0;
            while (remainder.IsNegative)
            {
                remainder += divisor;
                quotient -= One;
                count++;
                if (count >= 2)
                {
                    System.Diagnostics.Debugger.Break();
                }
            }
            while (remainder >= divisor)
            {
                remainder -= divisor;
                quotient += One;
                count++;
                if (count >= 2)
                {
                    System.Diagnostics.Debugger.Break();
                }
            }

            if (sign)
            {
                quotient = -quotient;
                remainder = -remainder;
            }
            return quotient;
        }

        public IntegerNumber Inverse()
        {
            shrink(ref this);
            IntegerNumber remainder, result = Inverse(this, out remainder);
            return result;
        }

        public IntegerNumber GetDigits(int start, int count, bool alwaysPositive)
        {
            ulong[] result = alwaysPositive && (this.bits[start + count - 1] & SignBit) != 0 ? new ulong[count + 1] : new ulong[count];
            int low = Math.Max(0, -start);
            for (int i = count; --i >= low;)
            {
                result[i] = this.bits[start + i];
            }
            return new IntegerNumber(result, false);
        }

        public static IntegerNumber Inverse(IntegerNumber number, out IntegerNumber remainder)
        {
            int n = number.Digits;
            if (n <= IntegerDigitsBreakPoint)
            {
                ulong[] a = new ulong[(n + 1) * 2];
                a[n * 2] = 1UL;
                IntegerNumber b = number.ExtendTo(n + 1);
                ulong[] q = new ulong[(n + 1) * 2];
                AsmX64Operations.GetDivMod(a, b.bits, n + 1, q, true);
                remainder = new IntegerNumber(a, false);
                shrink(ref remainder);
                IntegerNumber result = new IntegerNumber(q, false);
                shrink(ref result);
                return result;
            }

            bool isNegative = false;
            if (number.IsNegative)
            {
                number = -number;
                isNegative = true;
                n = number.Digits;
            }
            //Newton iteration: x <- x + x * (1 - d * x^2)
            int m = (n + 5) >> 1;
            IntegerNumber dx, dhi = number.GetDigits(n - m, m, false);
            IntegerNumber x = Inverse(dhi, out dx);

            //IntegerNumber test = x * d + dx;
            //shrink(ref test);
            //if (test.Digits != d.Digits * 2 + 1 || test.bits[d.Digits * 2] != 1UL || Enumerable.Range(0, d.Digits * 2).Any(idx => test.bits[idx] != 0UL))
            //{
            //    System.Diagnostics.Debugger.Break();
            //}

            const int CorrectionDigits = 4;
            IntegerNumber dlo = number.GetDigits(0, n - m, true);

            IntegerNumber dp = (dx << ((n - m) * 64)) - dlo * x;
            IntegerNumber delta = dp >> ((m - CorrectionDigits) * 64);  //keep one additional correction digit.
            shrink(ref delta);

            delta *= x;
            delta >>= (m + CorrectionDigits) * 64;
            shrink(ref delta);

            x <<= (n - m) * 64;
            x += delta;

            //FACT: number * x == (IntegerNumber.One << (n * 128)) - ((dp << ((n - m) * 64)) - delta * number);
            //but :     remainder = number * x;
            //then:     Array.Resize(ref remainder.bits, n * 2);
            //then:     shrink(ref remainder);
            //then:     remainder = -remainder;
            remainder = ((dp - delta * dhi) << ((n - m) * 64)) - delta * dlo;
            shrink(ref remainder);

            int count = 0;
            while (remainder.IsNegative)
            {
                remainder += number;
                x -= One;
                count++;
                if (count >= 2)
                {
                    System.Diagnostics.Debugger.Break();
                }
            }
            while (remainder >= number)
            {
                remainder -= number;
                x += One;
                count++;
                if (count >= 2)
                {
                    System.Diagnostics.Debugger.Break();
                }
            }

            if (isNegative)
            {
                x = -x;
                remainder = -remainder;
            }
            return x;
        }

        public long GetHighestPackedDigit(out int associatedShift)
        {
            shrink(ref this);
            int n = this.Digits;
            if (n <= 1)
            {
                associatedShift = 0;
                return (long)this.bits[0];
            }
            long result = (long)this.bits[n - 1];
            long highestDigit = result < 0 ? -result : result;
            int shift = 62 - BSR((ulong)highestDigit);
            associatedShift = (n - 1) * 64;
            if (shift > 0)
            {
                result = (result << shift) | (long)(this.bits[n - 2] >> (64 - shift));
                associatedShift -= shift;
            }
            return result;
        }

        private const double Power2_64 = 18446744073709551616.0;
        public double ToDouble
        {
            get
            {
                shrink(ref this);
                int n = this.Digits;
                if (n <= 1)
                {
                    return this.HighestDigit;
                }
                return (this.HighestDigit * Power2_64 + this.bits[n - 2]) * Math.Pow(Power2_64, n - 2);
            }
        }

        public static double Log2(IntegerNumber number)
        {
            if (number.IsNegative)
            {
                throw new InvalidOperationException("Real logarithm of a negative number does not exist.");
            }
            shrink(ref number);
            int n = number.Digits;
            if (n <= 1)
            {
                return number.IsZero ? 0.0 : Math.Log(number.HighestDigit) * LG_E;
            }
            return Math.Log(number.HighestDigit * Power2_64 + number.bits[n - 2]) * LG_E + (n - 2) * 64;
        }

        private static bool isSquareRoot(IntegerNumber number, IntegerNumber root)
        {
            IntegerNumber difference = number - root.Square();
            return !difference.IsNegative && difference <= (root << 1);
        }

        //f[x_] := x^2 - y
        //FullSimplify[x - f[x] / D[f[x], x]] = (x + y / x) / 2  [Newton method]
        //FullSimplify[x - 2*f[x]*D[f[x], x]/(2*D[f[x], x]^2 - f[x]*D[f[x], {x, 2}])] = x * (x^2 + 3 y) / (3 x^2 + y) [Halley]
        //In practice Newton method is faster than Halley 3'rd order convergence for SQRT.
        public IntegerNumber Sqrt()
        {
            if (this.IsNegative)
            {
                throw new ArithmeticException("NaN");
            }
            if (this.IsZero)
            {
                return IntegerNumber.Zero;
            }
            shrink(ref this);
            int n = this.Digits;
            int computedDigits = (n & 1) != 0 ? 1 : 2;
            if (computedDigits + 2 <= n)
            {
                computedDigits += 2;
            }

            ulong[] inputBits = this.bits;
            IntegerNumber t = this.GetDigits(n - computedDigits, computedDigits, true);
            double estimate = Math.Sqrt(t.ToDouble);
            double highDouble = Math.Floor(estimate / Power2_64);
            double lowDouble = estimate - highDouble * Power2_64;
            IntegerNumber root = new IntegerNumber(new ulong[] { (ulong)lowDouble, (ulong)highDouble }, false);
            if (root.IsNegative)
            {
                Array.Resize(ref root.bits, 3);
            }
            do
            {
                root = (root + t / root) >> 1;
            } while (!isSquareRoot(t, root));
            if (computedDigits == n)
            {
                return root;
            }

            do
            {
                int digits = Math.Min(computedDigits, n - computedDigits + 1) >> 1;
                root <<= digits * 64;
                computedDigits += digits * 2;
                t = computedDigits >= n ? this : this.GetDigits(n - computedDigits, computedDigits, true);
                root = (root + t / root) >> 1;
            } while (computedDigits < n || !isSquareRoot(t, root));
            return root;
        }

        //InverseSqrt[this] = Sqrt[(this << (this.Digits * 64)).Inverse()]
        //InverseSqrt[this] = Sqrt[(1 << (this.Digits * 192)) / this]
        //f[x_] := (1/x)^2 - y
        //FullSimplify[x - f[x] / D[f[x], x]] = x * (1/2) * (3 - x^2 y)
        public IntegerNumber InverseSqrt()
        {
            if (this.IsNegative || this.IsZero)
            {
                throw new ArithmeticException("NaN");
            }
            shrink(ref this);
            int n = this.Digits;
            int computedDigits = (n & 1) != 0 ? 1 : 2;
            if (computedDigits + 2 <= n)
            {
                computedDigits += 2;
            }

            ulong[] inputBits = this.bits;
            IntegerNumber t = this.GetDigits(n - computedDigits, computedDigits, true);
            double estimate = Math.Pow(Power2_64, computedDigits * 1.5) / Math.Sqrt(t.ToDouble);
            List<ulong> initialDigits = new List<ulong>();
            while (estimate != 0)
            {
                initialDigits.Add((ulong)Math.IEEERemainder(estimate, Power2_64));
                estimate = Math.Floor(estimate / Power2_64);
            }
            IntegerNumber root = new IntegerNumber(initialDigits.ToArray(), false);
            IntegerNumber three = new IntegerNumber(3L);
            //repeat iteration once to reach full precision.
            int iterations1 = (int)Math.Ceiling(Math.Log(initialDigits.Count * 64.0 / 52, 2.0));
            for (int step = iterations1; --step >= 0;)
            {
                int shift = computedDigits * 64 * 3;
                root = root * ((three << shift) - root.Square() * t) >> (1 + shift);
                shrink(ref root);
            }
            if (computedDigits == n)
            {
                return root;
            }

            int additionalDigits = 2;
            int totalDigits = n + additionalDigits;

            do
            {
                int digits = Math.Min(Math.Max(computedDigits - 1, 1), totalDigits - computedDigits);
                computedDigits += digits;
                t = this.GetDigits(n - computedDigits, computedDigits, true);

                //root <<= digits * 64;
                //int shift = computedDigits * 3 * 64;
                //root = root * ((three << shift) - root.Square() * t) >> (1 + shift);

                int shift1 = (computedDigits * 3 - digits * 2) * 64;
                int shift2 = (computedDigits * 3 - digits * 3) * 64;
                root = root * ((three << shift1) - root.Square() * t) >> (1 + shift2);

                shrink(ref root);
            } while (computedDigits < totalDigits);

            root >>= additionalDigits * 64;
            shrink(ref root);

            return root;
        }

        public static void shrink(ref IntegerNumber number)
        {
            int digits = number.RealDigits;
            if (digits != number.Digits)
            {
                Array.Resize(ref number.bits, digits);
            }
        }
    }

    public static class IntegerNumberUnitTest
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

            Random random = new Random(1002);
            BigInteger maxx = 0;
            BigInteger minn = 0;
            bool ok = true;
            for (int i = 2; --i >= 0;)
            {
                byte[] bytes1 = new byte[112 * 8 + random.Next(32 * 8)];
                byte[] bytes2 = new byte[112 * 8 + random.Next(32 * 8)];
                random.NextBytes(bytes1);
                random.NextBytes(bytes2);

                BigInteger n1 = new BigInteger(bytes1);
                BigInteger n2 = new BigInteger(bytes2);
                IntegerNumber f1 = new IntegerNumber(bytes1);
                IntegerNumber f2 = new IntegerNumber(bytes2);
                if (n1.ToString() != f1.ToString())
                {
                    ok = false;
                }
                if (n2.ToString() != f2.ToString())
                {
                    ok = false;
                }

                BigInteger a1 = n1 + n2;
                IntegerNumber a2 = f1 + f2;
                if (a1.ToString() != a2.ToString())
                {
                    ok = false;
                }

                BigInteger s1 = n1 - n2;
                IntegerNumber s2 = f1 - f2;
                if (s1.ToString() != s2.ToString())
                {
                    ok = false;
                }

                BigInteger m1 = n1 * n2;
                IntegerNumber m2 = f1 * f2;
                if (m1.ToString() != m2.ToString())
                {
                    ok = false;
                }

                int shrvalue = random.Next(256) + 1;
                BigInteger sh1 = n1 >> shrvalue;
                IntegerNumber sh2 = f1 >> shrvalue;
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
                IntegerNumber shl2 = f1 << shlvalue;
                if (shl1.ToString() != shl2.ToString())
                {
                    ok = false;
                }

                byte[] bytesINV = new byte[(192 + 32) * 8 + random.Next(64 * 8)];
                random.NextBytes(bytesINV);
                BigInteger num1 = new BigInteger(bytesINV);
                IntegerNumber num2 = new IntegerNumber(bytesINV);
                if (num1.ToString() != num2.ToString())
                {
                    ok = false;
                }
                BigInteger inv0 = (BigInteger.One << (num2.Digits * 64 * 2)) / num1;
                IntegerNumber inv1 = num2.Inverse();
                if (inv0.ToString() != inv1.ToString())
                {
                    ok = false;
                }

                byte[] bytes4 = new byte[bytes1.Length * 4];
                random.NextBytes(bytes4);
                BigInteger n4 = new BigInteger(bytes4);
                IntegerNumber f4 = new IntegerNumber(bytes4);
                BigInteger qb4 = n4 / n1;
                IntegerNumber qn4 = f4 / f1;
                if (qb4.ToString() != qn4.ToString())
                {
                    ok = false;
                }
                byte[] bytes3 = new byte[(bytes1.Length + bytes2.Length) & -8];
                random.NextBytes(bytes3);
                BigInteger square1 = BigInteger.Abs(new BigInteger(bytes3));
                IntegerNumber square2 = IntegerNumber.Abs(new IntegerNumber(bytes3));
                BigInteger root1 = square1.Sqrt();
                IntegerNumber root2 = square2.Sqrt();
                if (root1.ToString() != root2.ToString())
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
