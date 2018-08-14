using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Numerics;
using System.Diagnostics;

namespace Utilities.ReferenceRealNumber
{
    [DebuggerDisplay("{ToString()}")]
    public struct RealNumber
    {
        internal const int RealMaxBits = 2048;
        private const int ToStringDigits = 100;
        public static readonly double LN_2 = Math.Log(2);
        public static readonly double LG_E = 1.0 / Math.Log(2);
        public static readonly double LOG10_2 = Math.Log10(2);
        public static readonly double LOG5_2 = Math.Log(2, 5);

        public double ToDouble
        {
            get
            {
                int shl = 53 - this.mantissaBits;
                double mantissa = (double)(shl >= 0 ? this.mantissa << shl : this.mantissa >> (-shl));
                double result = mantissa * Math.Exp((this.shift - shl) * LN_2);
                return result;
            }
        }

        public long ToLong
        {
            get
            {
                long highestBit = this.shift + this.mantissaBits;
                if (highestBit < -1)
                {
                    return 0;
                }
                if (highestBit > 65)
                {
                    throw new ArgumentOutOfRangeException("ToLong");
                }
                BigInteger mantissa = this.shift >= 0 ? this.mantissa << (int)this.shift : this.mantissa >> (-(int)this.shift);
                if (mantissa < long.MinValue || mantissa > long.MaxValue)
                {
                    throw new ArgumentOutOfRangeException("ToLong");
                }
                return (long)mantissa;
            }
        }

        private long shift;
        private BigInteger mantissa;
        //value = 2^(shift) * mantissa

        public bool IsNegative { get { return this.mantissa.Sign < 0; } }
        public bool IsZero { get { return this.mantissa.Sign == 0; } }

        public static implicit operator RealNumber(long x)
        {
            return new RealNumber(x);
        }

        public static implicit operator RealNumber(double x)
        {
            return new RealNumber(x);
        }

        public RealNumber(long mantissa)
        {
            this.shift = 0;
            this.mantissa = new BigInteger(mantissa);
            this.simplify();
        }

        public RealNumber(long shift, BigInteger mantissa)
        {
            this.shift = shift;
            this.mantissa = mantissa;
            this.simplify();
        }

        private void simplify()
        {
            if (this.mantissa.IsZero)
            {
                return;
            }
            int overflowBits = this.mantissaBits - RealMaxBits;
            if (overflowBits > 0)
            {
                this.shift += overflowBits;
                this.mantissa = (this.mantissa + (BigInteger.One << (overflowBits - 1))) >> overflowBits;
            }
            if (this.mantissa.IsEven)
            {
                BigInteger bit = BigInteger.One;
                int bitShift = 0;
                while ((this.mantissa & bit).IsZero)
                {
                    bitShift++;
                    bit <<= 1;
                }
                this.mantissa >>= bitShift;
                this.shift += bitShift;
            }
        }

        private int mantissaBits { get { return (int)this.mantissaBitsDouble; } }
        private double mantissaBitsDouble { get { return this.mantissa.IsZero ? 0 : Math.Ceiling(BigInteger.Log(BigInteger.Abs(this.mantissa)) * LG_E); } }

        public RealNumber(double x)
            : this()
        {
            if (x == 0)
            {
                this.shift = 0;
                this.mantissa = BigInteger.Zero;
                return;
            }
            bool sign = x < 0;
            x = Math.Abs(x);
            this.shift = (long)Math.Ceiling(Math.Log(x) * LG_E - 53);
            x *= Math.Pow(0.5, this.shift);
            while (x != Math.Floor(x))
            {
                this.shift--;
                x *= 2;
            }
            this.mantissa = new BigInteger(sign ? -x : x);
            this.simplify();
        }

        public static RealNumber operator *(RealNumber x, RealNumber y)
        {
            RealNumber result = new RealNumber(x.shift + y.shift, x.mantissa * y.mantissa);
            return result;
        }

        public static RealNumber operator /(RealNumber x, RealNumber y)
        {
            int shift = RealMaxBits - x.mantissaBits + y.mantissaBits;
            BigInteger resultMantissa = BigInteger.Divide(x.mantissa << shift, y.mantissa);
            return new RealNumber(x.shift - y.shift - shift, resultMantissa);
        }

        public static readonly RealNumber Zero = new RealNumber(0, BigInteger.Zero);

        public RealNumber Abs { get { return new RealNumber(this.shift, BigInteger.Abs(this.mantissa)); } }

        public RealNumber Sqrt
        {
            get
            {
                if (this.mantissa.IsZero)
                {
                    return RealNumber.Zero;
                }
                if (this.IsNegative)
                {
                    throw new ArithmeticException("NaN");
                }
                int shift = RealMaxBits * 2 - this.mantissaBits;
                shift += ((int)this.shift ^ shift) & 1;
                RealNumber root = new RealNumber((this.shift - shift) / 2, (this.mantissa << shift).Sqrt());
                return root;
            }
        }

        public static RealNumber operator -(RealNumber x)
        {
            return new RealNumber(x.shift, -x.mantissa);
        }

        public static RealNumber operator <<(RealNumber x, int shift)
        {
            return new RealNumber(x.shift + shift, x.mantissa);
        }

        public static RealNumber operator >>(RealNumber x, int shift)
        {
            return new RealNumber(x.shift - shift, x.mantissa);
        }

        public static RealNumber operator +(RealNumber x, RealNumber y)
        {
            if (x.shift < y.shift)
            {
                return y + x;
            }
            if (x.IsZero || y.IsZero)
            {
                return x.IsZero ? y : x;
            }
            //now: x.shift >= y.shift
            double highestBitX = x.mantissaBitsDouble + x.shift;
            double highestBitY = y.mantissaBitsDouble + y.shift;
            if (Math.Abs(highestBitX - highestBitY) >= RealMaxBits + 4)
            {
                return highestBitX > highestBitY ? x : y;
            }
            int diff = (int)(x.shift - y.shift);
            BigInteger mantissa = (x.mantissa << diff) + y.mantissa;
            return new RealNumber(y.shift, mantissa);
        }

        public static RealNumber operator -(RealNumber x, RealNumber y)
        {
            return x + (-y);
        }

        public static bool operator <(RealNumber x, RealNumber y)
        {
            RealNumber diff = x - y;
            return diff.IsNegative;
        }

        public static bool operator >(RealNumber x, RealNumber y)
        {
            RealNumber diff = y - x;
            return diff.IsNegative;
        }

        public static bool operator >=(RealNumber x, RealNumber y)
        {
            RealNumber diff = y - x;
            return diff.IsNegative || diff.IsZero;
        }

        public static bool operator <=(RealNumber x, RealNumber y)
        {
            RealNumber diff = x - y;
            return diff.IsNegative || diff.IsZero;
        }

        public RealNumber Log
        {
            get
            {
                if (this.IsNegative)
                {
                    throw new ArithmeticException("NaN");
                }
                if (this.IsZero)
                {
                    return new RealNumber(int.MaxValue, BigInteger.MinusOne);   // - infinity
                }
                RealNumber x = this;
                long shiftLeft = RealMaxBits / 2 - (x.mantissaBits + x.shift);
                x = new RealNumber(x.shift + shiftLeft, x.mantissa);
                RealNumber log = (PI >> 1) * x / AGM(x, 4) - LN2 * shiftLeft;   //Error = O(1 + 1 / x^2) = O(1 + 1 / 2^RealMaxBits) [OK]
                return log;
            }
        }

        #region Exponential

        /*
        https://en.wikipedia.org/wiki/Householder%27s_method

        solve f(x) = ln(x) - y = 0; find x so f(x) = 0 ==> x = exp(y)
        iteration x(n+1) = x(n) + d * (1/f)'[d-1](x(n)) / (1/f)'[d](x(n))

        inv_f(f(x)) = x => inv_f(ln(x)-y) = x => inv_f(x-y)= exp(x) => inv_f(x) = exp(x+y) = e^x * e^y
        TEST:
        inv_f(f(x)) = exp(ln(x)-y+y) = x [ok]
        f(inv_f(x)) = ln(exp(x+y)) - y = x [ok].

        Identities:
        (1/f)'(x)  = -f'(x)/f(x)^2
        (1/f)''(x) = -f''(x)/f(x)^2 + 2 * f'(x)^2 / f(x)^3

        ==> x_next <-- x + 2 *[ -f'(x)/f(x)^2 ] / [ -f''(x)/f(x)^2 + 2 * f'(x)^2 / f(x)^3 ]
            x_next <-- x + 2 *[ -f'(x)*f(x)   ] / [ -f''(x)*f(x) + 2 * f'(x)^2 ]


        f'(x) = 1/x
        f''(x) = -1/x^2

            x_next <-- x + 2 *[ -(1/x)*(ln(x) - y)   ] / [ -(-1/x^2)*(ln(x) - y) + 2 * (1/x)^2 ]
            x_next <-- x + 2 *[ -x * (ln(x) - y)   ] / [ ln(x) - y + 2 ]
            x_next <-- x * { 1 - 2 *[ ln(x) - y ] / [ ln(x) - y + 2 ] }
            x_next <-- x * { [ 2 - (ln(x) - y) ]  / [ (ln(x) - y) + 2 ] }
        */

        public RealNumber Exp
        {
            get
            {   //2^x = exp(this) ==> exp(x * ln(2)) = exp(this) => x * ln(2) = this => x = this * log2(e) = this / ln(2)
                RealNumber x = this / LN2;
                long k = x.ToLong;
                RealNumber log = x - k;
                if (log.IsZero)
                {
                    return new RealNumber(k, BigInteger.One);
                }
                //x = k + log ==> 2^x == 2^(k+log) == 2^k * 2^log = 2^k * exp(log * ln(2))
                log *= LN2;
                RealNumber exp = new RealNumber(Math.Exp(log.ToDouble));
                for (int bits = 50; bits < RealMaxBits; bits *= 3)
                {
                    RealNumber cache = exp.Log - log;
                    exp *= (2 - cache) / (2 + cache);
                }
                RealNumber result = new RealNumber(exp.shift + k, exp.mantissa);
                return result;
            }
        }

        #endregion

        private static RealNumber GetPI()
        {
            RealNumber sqrt2 = new RealNumber(2).Sqrt;
            RealNumber y = sqrt2 - 1;
            RealNumber a = 6 - 4 * sqrt2;
            RealNumber error = 1;
            for (int n = 1; ; n++)
            {
                RealNumber y2 = y * y;
                RealNumber root4 = (RealNumber.One - y2 * y2).Sqrt.Sqrt;
                y = (RealNumber.One - root4) / (RealNumber.One + root4);
                y2 = y + 1;
                y2 *= y2;
                RealNumber save = y2;
                y2 *= y2;
                RealNumber nexta = y2 * a - (((save - y) * y) << (n * 2 + 1));
                RealNumber currentError = (nexta - a).Abs;
                if (currentError >= error)
                {
                    return RealNumber.One / a;
                }
                a = nexta;
                error = currentError;
            }
        }

        private static RealNumber GetLN2()
        {
            int shift = RealMaxBits / 2;
            RealNumber x = new RealNumber(shift, BigInteger.One);
            RealNumber ln2 = (PI >> 1) * x / AGM(x, 4) / shift;
            return ln2;
        }

        public static readonly RealNumber One = new RealNumber(0, BigInteger.One);
        public static readonly RealNumber PI = GetPI();
        public static readonly RealNumber LN2 = GetLN2();
        public static readonly RealNumber LN10 = new RealNumber(0, new BigInteger(10)).Log;

        public static RealNumber AGM(RealNumber x, RealNumber y)
        {
            RealNumber previousError = (x - y).Abs;
            while (true)
            {
                RealNumber am = (x + y) >> 1;
                RealNumber gm = (x * y).Sqrt;
                RealNumber error = (am - gm).Abs;
                if (error >= previousError)
                {
                    return am;
                }
                x = am;
                y = gm;
                previousError = error;
            }
        }

        private static void shiftRightInBase10(ref BigInteger number, ref long shift10)
        {
            if (number.IsZero)
            {
                return;
            }
            int keepDigits = (int)Math.Floor(LOG10_2 * RealMaxBits);
            int shift = (int)Math.Floor(BigInteger.Log10(BigInteger.Abs(number))) - keepDigits;
            number = shift >= 0 ? number / BigInteger.Pow(10, shift) : number * BigInteger.Pow(10, -shift);
            shift10 += shift;
        }

        private static BigInteger mostSignificantExponent(BigInteger number, long power, out long shift10)
        {
            shift10 = 0;
            if (power == 0)
            {
                return BigInteger.One;
            }
            int bit = 63;
            while ((power & (1L << bit)) == 0)
            {
                bit--;
            }
            BigInteger result = number;
            int keepDigits = (int)Math.Floor(LOG10_2 * RealMaxBits);
            while (--bit >= 0)
            {
                result *= result;
                shift10 *= 2;
                shiftRightInBase10(ref result, ref shift10);
                if ((power & (1L << bit)) != 0)
                {
                    result *= number;
                    shiftRightInBase10(ref result, ref shift10);
                }
            }
            return result;
        }

        private string getStringRepresentation()
        {
            long shift = this.shift, shift10;
            BigInteger scale;
            if (shift >= 0)
            {
                int bucketBits = RealMaxBits;
                BigInteger start = BigInteger.One << bucketBits;
                long exponent = shift / bucketBits;
                scale = mostSignificantExponent(start, exponent, out shift10);
                scale <<= (int)(shift % bucketBits);
                shiftRightInBase10(ref scale, ref shift10);
            }
            else
            {
                shift = -shift;
                int bucketFives = (int)Math.Ceiling(LOG5_2 * RealMaxBits);
                BigInteger start = BigInteger.Pow(5, bucketFives);
                long exponent = shift / bucketFives;
                scale = mostSignificantExponent(start, exponent, out shift10);
                scale *= BigInteger.Pow(5, (int)(shift % bucketFives));
                shiftRightInBase10(ref scale, ref shift10);
                shift10 -= shift;
            }
            BigInteger mantissa = scale * this.mantissa;
            int digits = (int)Math.Floor(BigInteger.Log10(BigInteger.Abs(mantissa))) - (ToStringDigits + 1);
            shift10 += digits;
            mantissa = digits >= 0 ? mantissa / BigInteger.Pow(10, digits) : mantissa * BigInteger.Pow(10, -digits);

            BigInteger rem, div = BigInteger.DivRem(mantissa, 10, out rem);
            mantissa = div + (rem <= -5 ? BigInteger.MinusOne : (rem >= 5 ? BigInteger.One : BigInteger.Zero));   //Rounding in base 10.
            shift10++;

            while (!mantissa.IsZero)
            {
                div = BigInteger.DivRem(mantissa, 10, out rem);
                if (!rem.IsZero)
                {
                    break;
                }
                mantissa = div;
                shift10++;
            }
            string result = mantissa.ToString();
            int dotIndex = (mantissa.Sign < 0 ? 1 : 0) + 1;
            result = result.Insert(dotIndex, ".");
            shift10 += result.Length - (dotIndex + 1);
            if (result[result.Length - 1] == '.')
            {
                result = result.Remove(result.Length - 1);
            }
            if (shift10 != 0)
            {
                result += "e" + (shift10 > 0 ? "+" : "") + shift10.ToString();
            }
            return result;
        }

        public override string ToString()
        {
            if (this.IsZero)
            {
                return "0";
            }
            string firstDigits = this.getStringRepresentation();
            return firstDigits;
        }

        public static RealNumber Max(RealNumber x, RealNumber y)
        {
            RealNumber diff = x - y;
            return diff.IsNegative ? y : x;
        }
    }

    [DebuggerDisplay("{ToString()}")]
    public struct ComplexNumber
    {
        public RealNumber Real { get; set; }
        public RealNumber Imaginary { get; set; }

        public static readonly ComplexNumber ImaginaryOne = new ComplexNumber(RealNumber.Zero, RealNumber.One);
        public static readonly ComplexNumber One = new ComplexNumber(RealNumber.One, RealNumber.Zero);

        public static implicit operator ComplexNumber(RealNumber real)
        {
            return new ComplexNumber(real, RealNumber.Zero);
        }

        public ComplexNumber(RealNumber real, RealNumber imaginary)
            : this()
        {
            this.Real = real;
            this.Imaginary = imaginary;
        }

        public static ComplexNumber operator *(ComplexNumber x, ComplexNumber y)
        {
            return new ComplexNumber(
                x.Real * y.Real - x.Imaginary * y.Imaginary,
                x.Real * y.Imaginary + x.Imaginary * y.Real);
        }

        public static ComplexNumber operator *(ComplexNumber x, RealNumber y)
        {
            return new ComplexNumber(x.Real * y, x.Imaginary * y);
        }

        public static ComplexNumber operator +(ComplexNumber x, ComplexNumber y)
        {
            return new ComplexNumber(x.Real + y.Real, x.Imaginary + y.Imaginary);
        }

        public static ComplexNumber operator -(ComplexNumber x, ComplexNumber y)
        {
            return new ComplexNumber(x.Real - y.Real, x.Imaginary - y.Imaginary);
        }

        public RealNumber Energy { get { return this.Real * this.Real + this.Imaginary * this.Imaginary; } }
        public RealNumber Magnitude { get { return this.Energy.Sqrt; } }

        public override string ToString()
        {
            return "{" + this.Real.ToString() + " ; " + this.Imaginary.ToString() + "}";
        }
    }

    public static class RealNumbersUnitTest
    {
        public static bool UnitTest(int seed)
        {
            bool ok = true;
            RealNumber log0 = new RealNumber(1L << 62);
            RealNumber exp1 = log0.Exp;
            RealNumber log1 = exp1.Log;
            ok &= (log1 - log0).IsZero;

            RealNumber n_2 = RealNumber.LN2.Exp;
            RealNumber n_1a = new RealNumber(1.0).Exp.Log;
            RealNumber n_1b = new RealNumber(10.0).Log.Exp;
            var s1 = n_1b.ToString();
            RealNumber error = RealNumber.Max(RealNumber.Max((n_2 - 2).Abs, (n_1a - 1).Abs), (n_1b - 10).Abs);
            var s2 = error.ToString();
            RealNumber tolerance = RealNumber.One >> 2000;
            ok &= error < tolerance;
            var s3 = new RealNumber(10000, BigInteger.Pow(3, 10000)).ToString();
            ok &= s3.StartsWith("3.25464658549366205883184429132") && s3.EndsWith("+7781");
            var s4 = new RealNumber(-25000, BigInteger.Pow(1, 1)).ToString();
            ok &= s4.StartsWith("1.778723326301851865925552009580") && s4.EndsWith("-7526");
            return ok;
        }
    }
}
