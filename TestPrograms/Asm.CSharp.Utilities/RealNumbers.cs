using System;
using System.Collections.Concurrent;
using System.Collections.Generic;
using System.Diagnostics;
using System.Numerics;

namespace Utilities
{
    [DebuggerDisplay("{ToString()}")]
    public struct RealNumber
    {
        //internal const int RealMaxBits = 2048;
        //private const int ToStringDigits = 100 * 100;
        private const double LN_2 = 0.69314718055994530941723212145818;  // = ln(2).
        private const double LG_E = 1.4426950408889634073599246810019;  // = log2(e) = 1/ln(2).
        private static readonly double LOG10_2 = Math.Log10(2);
        private static readonly double LOG5_2 = Math.Log(2, 5);

        private int RealMaxBits { get { return this.totalDigits * 64; } }
        private int ToStringDigits { get { return (int)Math.Ceiling(RealMaxBits * LOG10_2); } }

        public static readonly RealNumber Zero = new RealNumber(0L, 2);

        public static implicit operator RealNumber(long x)
        {
            return new RealNumber(x, 2);
        }

        public static implicit operator RealNumber(double x)
        {
            return new RealNumber(x, 2);
        }

        //value = 2^(shift) * mantissa
        private long shift;
        private IntegerNumber mantissa;
        internal int totalDigits;

        public RealNumber(int numberOf64BitDigits) : this()
        {   //value of the number is 0.
            this.totalDigits = numberOf64BitDigits;
        }

        public RealNumber(long mantissa, int numberOf64BitDigits) : this(numberOf64BitDigits)
        {
            this.mantissa = new IntegerNumber(mantissa);
        }

        public RealNumber(long shift, IntegerNumber mantissa, int numberOf64BitDigits) : this(numberOf64BitDigits)
        {
            this.shift = shift;
            this.mantissa = mantissa;
            this.simplify();
        }

        public RealNumber(double x, int numberOf64BitDigits) : this(numberOf64BitDigits)
        {
            if (x == 0)
            {
                this.shift = 0;
                this.mantissa = IntegerNumber.Zero;
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
            long mantissa = (long)x;
            this.mantissa = new IntegerNumber(sign ? -mantissa : mantissa);
            this.simplify();
        }

        private void simplify()
        {
            IntegerNumber.shrink(ref this.mantissa);
            if (this.mantissa.IsZero)
            {
                return;
            }
            int overflowBits = this.mantissaBits + 2 - RealMaxBits;
            if (overflowBits > 0)
            {
                this.shift += overflowBits;
                this.mantissa = (this.mantissa + (IntegerNumber.One << (overflowBits - 1))) >> overflowBits;
            }
            if (this.mantissa.IsEven)
            {
                int bitShift = 0;
                while (!this.mantissa.HasBit(bitShift))
                {
                    bitShift++;
                }
                this.mantissa >>= bitShift;
                this.shift += bitShift;
            }
        }

        public double ToDouble
        {
            get
            {
                int shl = 53 - this.mantissaBits;
                long mantissa = (shl >= 0 ? this.mantissa << shl : this.mantissa >> (-shl)).ToLong;
                double result = mantissa * Math.Exp((this.shift - shl) * LN_2);
                return result;
            }
        }

        public double EstimateLog2 { get { return IntegerNumber.Log2(this.mantissa) + this.shift; } }

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
                IntegerNumber mantissa = this.shift >= 0 ? this.mantissa << (int)this.shift : this.mantissa >> (-(int)this.shift);
                if (mantissa < long.MinValue || mantissa > long.MaxValue)
                {
                    throw new ArgumentOutOfRangeException("ToLong");
                }
                return mantissa.ToLong;
            }
        }

        public bool IsNegative { get { return this.mantissa.IsNegative; } }
        public bool IsZero { get { return this.mantissa.IsZero; } }
        private int mantissaBits { get { return this.mantissa.BitsCount; } }

        public static RealNumber operator *(RealNumber x, RealNumber y)
        {
            RealNumber result = new RealNumber(x.shift + y.shift, x.mantissa * y.mantissa, Math.Max(x.totalDigits, y.totalDigits));
            return result;
        }

        public static RealNumber operator *(RealNumber x, IntegerNumber y)
        {
            RealNumber result = new RealNumber(x.shift, x.mantissa * y, x.totalDigits);
            return result;
        }

        public static RealNumber operator *(RealNumber x, long y)
        {
            RealNumber result = new RealNumber(x.shift, x.mantissa * y, x.totalDigits);
            return result;
        }

        public static RealNumber operator *(long x, RealNumber y)
        {
            RealNumber result = new RealNumber(y.shift, y.mantissa * x, y.totalDigits);
            return result;
        }

        public RealNumber Square()
        {
            RealNumber result = new RealNumber(this.shift * 2, this.mantissa.Square(), this.totalDigits);
            return result;
        }

        public static RealNumber operator /(RealNumber x, RealNumber y)
        {
            int resultDigits = Math.Max(x.totalDigits, y.totalDigits);
            int shift = resultDigits * 64 - x.mantissaBits + y.mantissaBits;
            IntegerNumber remainder;
            IntegerNumber quotient = IntegerNumber.DivRem(x.mantissa << shift, y.mantissa, out remainder);
            return new RealNumber(x.shift - y.shift - shift, quotient, resultDigits);
        }

        public static RealNumber operator /(RealNumber x, long y)
        {
            RealNumber result = new RealNumber(x.shift - 64, (x.mantissa << 64) / y, x.totalDigits);
            return result;
        }

        public static RealNumber operator /(RealNumber x, IntegerNumber y)
        {
            int shift = x.totalDigits * 64 - x.mantissaBits + y.BitsCount;
            IntegerNumber remainder;
            IntegerNumber quotient = IntegerNumber.DivRem(x.mantissa << shift, y, out remainder);
            return new RealNumber(x.shift - shift, quotient, x.totalDigits);
        }

        public RealNumber Inverse()
        {
            IntegerNumber mantissa = this.mantissa;
            int shl = this.RealMaxBits - mantissa.BitsCount - 1;
            if (shl > 0)
            {
                mantissa <<= shl;
            }
            IntegerNumber quotient = mantissa.Inverse();    //quotient = (2^64 << this.mantissa.Digits * 2) / this.mantissa
            return new RealNumber(shl - this.shift - mantissa.Digits * 128L, quotient, this.totalDigits);
        }

        public static RealNumber Abs(RealNumber x)
        {
            return new RealNumber(x.shift, IntegerNumber.Abs(x.mantissa), x.totalDigits);
        }

        public RealNumber GetSlowerSqrt()
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
            RealNumber root = new RealNumber((this.shift - shift) / 2, (this.mantissa << shift).Sqrt(), this.totalDigits);
            return root;
        }

        public RealNumber GetSqrt()
        {
            if (this.IsZero)
            {
                return this;
            }
            return this.GetInverseSqrt().Inverse();
        }

        //f[x_] := (1/x)^2 - y
        //FullSimplify[x - f[x] / D[f[x], x]] = x * (1/2) * (3 - x^2 y)
        //FullSimplify[x - 2 * f[x] * D[f[x], x] / (2 * D[f[x], x] ^ 2 - f[x] * D[f[x], {x, 2}])] = x * (3 + x^2 y) / (1 + 3 x^2 y)
        public RealNumber GetInverseSqrt()
        {
            if (this.mantissa.IsZero)
            {
                throw new ArithmeticException("1 / sqrt(0) error.");
            }
            if (this.IsNegative)
            {
                throw new ArithmeticException("NaN");
            }
            //this = 2^(shift+mantissaBits) * [mantissa >> mantissaBits]
            //sqrt(this) = 2^([shift+mantissaBits]/2) * sqrt[mantissa >> mantissaBits]
            //1/sqrt(this) = 2^[-(shift+mantissaBits)/2] / { sqrt[mantissa >> mantissaBits] }
            int mantissaBits = this.mantissaBits;
            long resultShift = this.shift + mantissaBits;
            if ((resultShift & 1L) != 0)
            {
                mantissaBits++;
                resultShift++;
            }
            resultShift >>= 1;
            RealNumber source = new RealNumber(-mantissaBits, this.mantissa, this.totalDigits + 2);
            int totalBits = source.totalDigits * 64;
            RealNumber result = new RealNumber(1.0 / Math.Sqrt(source.ToDouble), 2);
            RealNumber three = new RealNumber(3L, 1);
            for (int bits = 52; bits < totalBits; bits *= 2)
            {
                int precision = Math.Min(source.totalDigits, (bits * 2 + 63) >> 6);
                result = result.ChangePrecision(precision);
                RealNumber x = source.ChangePrecision(precision);
                result *= (three - x * result.Square()) >> 1;
            }
            RealNumber inverseSqrt = new RealNumber(result.shift - resultShift, result.mantissa, this.totalDigits);
            return inverseSqrt;
        }

        public static RealNumber operator -(RealNumber x)
        {
            return new RealNumber(x.shift, -x.mantissa, x.totalDigits);
        }

        public static RealNumber operator <<(RealNumber x, int shift)
        {
            return new RealNumber(x.shift + shift, x.mantissa, x.totalDigits);
        }

        public static RealNumber operator >>(RealNumber x, int shift)
        {
            return new RealNumber(x.shift - shift, x.mantissa, x.totalDigits);
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
            int realMaxBits = Math.Max(x.RealMaxBits, y.RealMaxBits);
            long highestBitX = x.mantissaBits + x.shift;
            long highestBitY = y.mantissaBits + y.shift;
            if (Math.Abs(highestBitX - highestBitY) >= realMaxBits + 4)
            {
                return highestBitX > highestBitY ? x : y;
            }
            int diff = (int)(x.shift - y.shift);
            IntegerNumber mantissa = (x.mantissa << diff) + y.mantissa;
            return new RealNumber(y.shift, mantissa, Math.Max(x.totalDigits, y.totalDigits));
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

        private static RealNumber getLN2Value(int numberOf64BitDigits)
        {
            long shift = ((numberOf64BitDigits + 2) * 64L + 35) / 36;

            RealNumber q = new RealNumber(-shift, IntegerNumber.One, numberOf64BitDigits + 2);
            RealNumber q4 = q.Square().Square();
            RealNumber q8 = q4.Square();
            RealNumber q16 = q8.Square();
            RealNumber q24 = q16 * q8;
            RealNumber rootFirstArgument = (RealNumber.One + q8 + q24) * q << 1;
            RealNumber rootSecondArgument = RealNumber.One + ((q4 + q16) << 1);
            RealNumber arg1 = rootFirstArgument * rootSecondArgument << 1;
            RealNumber arg2 = (rootFirstArgument + rootSecondArgument).Square() - arg1;

            RealNumber pi = GetPI(numberOf64BitDigits + 2);
            RealNumber ln2 = (pi >> 1) / AGMReal(arg1, arg2) / shift;
            RealNumber result = new RealNumber(ln2.shift, ln2.mantissa, numberOf64BitDigits);
            return result;
        }

        public RealNumber GetLog()
        {
            if (this.IsNegative)
            {
                throw new ArithmeticException("NaN");
            }
            if (this.IsZero)
            {
                return new RealNumber(int.MaxValue, IntegerNumber.MinusOne, this.totalDigits);   // - infinity
            }
            long shiftLeft = ((this.totalDigits + 2) * 64L + 35) / 36 - (this.mantissaBits + this.shift);

            RealNumber x = new RealNumber(this.shift + shiftLeft, this.mantissa, this.totalDigits + 2);
            RealNumber q = x.Inverse();
            RealNumber q4 = q.Square().Square();
            RealNumber q8 = q4.Square();
            RealNumber q16 = q8.Square();
            RealNumber q24 = q16 * q8;
            RealNumber rootFirstArgument = (RealNumber.One + q8 + q24) * q << 1;
            RealNumber rootSecondArgument = RealNumber.One + ((q4 + q16) << 1);
            RealNumber arg1 = rootFirstArgument * rootSecondArgument << 1;
            RealNumber arg2 = (rootFirstArgument + rootSecondArgument).Square() - arg1;

            RealNumber pi = GetPI(this.totalDigits + 2);
            RealNumber ln2 = GetLN2(this.totalDigits + 2);
            RealNumber log = (pi >> 1) / AGMReal(arg1, arg2) - ln2 * shiftLeft;   //need: x >= 2^(d/36) so O(q^36) are negligible.
            RealNumber result = log.ChangePrecision(this.totalDigits);
            return result;
        }

        //Sasaki and Kanada's formula: (for log computation)
        //ln(1/q) = PI / [AGM(t2(q)^2, t3(q)^2] ; 1/q = x
        //ln(x) = PI / 4 / [AGM(t2(q^4)^2, t3(q^4)^2]
        //t2(q^4) = 2 [ q + q^9 + q^25 + O(q^49) ]
        //t3(q^4) = 1 + 2 [ q^4 + q^16 + O(q^36) ] ; so we need x >= 2^[d/36]
        //AGM(t2^2, t3^2) = AGM(2*t2*t3, t2^2 + t3^2) / 2   ;   --> saves one square root.    
        public static ComplexNumber GetComplexLog(ComplexNumber argument)
        {
            int totalDigits = Math.Max(argument.Real.totalDigits, argument.Imaginary.totalDigits);
            long shiftLeft = Math.Max(
                ((totalDigits + 2) * 64L + 35) / 36 - (argument.Real.mantissaBits + argument.Real.shift),
                ((totalDigits + 2) * 64L + 35) / 36 - (argument.Imaginary.mantissaBits + argument.Imaginary.shift));

            ComplexNumber x = new ComplexNumber(
                new RealNumber(argument.Real.shift + shiftLeft, argument.Real.mantissa, totalDigits + 2),
                new RealNumber(argument.Imaginary.shift + shiftLeft, argument.Imaginary.mantissa, totalDigits + 2));
            ComplexNumber q = x.Inverse();
            ComplexNumber q4 = q.Square().Square();
            ComplexNumber q8 = q4.Square();
            ComplexNumber q16 = q8.Square();
            ComplexNumber q24 = q16 * q8;
            ComplexNumber rootFirstArgument = (ComplexNumber.One + q8 + q24) * q << 1;
            ComplexNumber rootSecondArgument = ComplexNumber.One + ((q4 + q16) << 1);
            ComplexNumber arg1 = rootFirstArgument * rootSecondArgument << 1;
            ComplexNumber arg2 = (rootFirstArgument + rootSecondArgument).Square() - arg1;

            RealNumber pi = GetPI(totalDigits + 2);
            RealNumber ln2 = GetLN2(totalDigits + 2);
            ComplexNumber log = (pi >> 1) / AGMComplex(arg1, arg2);     //need: x >= 2^(d/36) so O(q^36) are negligible.
            ComplexNumber result = new ComplexNumber(
                (log.Real - ln2 * shiftLeft).ChangePrecision(totalDigits),
                log.Imaginary.ChangePrecision(totalDigits));
            return result;
        }

        #region Taylor expansion exponentiation

        private static double getApproximativeLogFactorial(int n)
        {
            const double ln_sqrt_2pi = 0.91893853320467274178032973640562;  // = Math.Log(Math.PI * 2) * 0.5
            double lnN = Math.Log(n);
            return ln_sqrt_2pi + lnN * 0.5 + n * (lnN - 1);
        }

        private static int getNumberOfExpFactors(int powerOf2, int totalDigits)
        {   //find n so: [2^powerOf2]^n / n! <= 2^(-64 * totalDigits) or:
            //<==> n * powerOf2 - Ln[n!] / Ln[2] <= -64 * totalDigits 
            double rightValue = -64.0 * totalDigits;
            double lgExp = powerOf2;
            int n = 1;
            while (n * lgExp - getApproximativeLogFactorial(n) * LG_E > rightValue)
            {
                n <<= 1;
            }
            int left = n >> 1;
            int right = n;
            while (left < right)
            {   // [st m0 M1 en]
                n = left + ((right - left + 1) >> 1);   //rightmost strictly less.
                if (n * lgExp - getApproximativeLogFactorial(n) * LG_E > rightValue)
                {
                    left = n;
                }
                else
                {
                    right = n - 1;
                }
            }
            return left;
        }

        private static IntegerNumber getIntegerProduct(int start, int end)
        {
            if (start >= end)
            {
                return start == end ? new IntegerNumber((long)start) : IntegerNumber.One;
            }
            int middle = start + ((end - start) >> 1);
            IntegerNumber left = getIntegerProduct(start, middle);
            IntegerNumber right = getIntegerProduct(middle + 1, end);
            IntegerNumber result = left * right;
            IntegerNumber.shrink(ref result);
            return result;
        }

        //result = sum[i = start to end: of startPower * exp^[i-start] * product[i + 1 to end]
        private static RealNumber getLocalExponentialSummation(int start, int end, ref RealNumber startPower, RealNumber exp, out IntegerNumber product)
        {
            if (end == start)
            {
                product = new IntegerNumber((long)start);
                return startPower;
            }
            RealNumber save = startPower;
            startPower *= exp;
            RealNumber result = getLocalExponentialSummation(start + 1, end, ref startPower, exp, out product);
            result += save * product;
            product *= start;
            IntegerNumber.shrink(ref product);
            return result;
        }

        /// <summary>
        /// returns A = sum(i = start to end of: exp ^ i / i!) 
        /// returns B = [end! / (start - 1)!] = product from [start to end] inclusive
        /// A * B = Exp[exp]
        /// </summary>
        private static RealNumber taylorSeriesExponent(int start, int end,
            RealNumber exp,
            ref RealNumber startPower,          // startPower = exp ^ start [is input]
            out IntegerNumber factorial)	    // = end! / (start - 1)! = product from [start to end] inclusive
        {
            //if (start == end)
            //{
            //    factorial = start == 0 ? IntegerNumber.One : new IntegerNumber((long)start);
            //    return startPower;
            //}
            if (start + 24 >= end)
            {
                RealNumber result = getLocalExponentialSummation(start, end, ref startPower, exp, out factorial);
                return result;
                ///result = sum[i = start to end: of startPower * exp^(i-start) * product[i + 1 to end]
                //factorial = getIntegerProduct(start, end);
                //IntegerNumber product = factorial / start;
                //RealNumber result = startPower * product;
                //for (int i = start + 1; i <= end; i++)
                //{
                //    startPower *= exp;
                //    product /= i;
                //    result += startPower * product;
                //}
                //return result;
            }

            int middle = start + ((end - start) >> 1);
            IntegerNumber secondFactorial;
            RealNumber left = taylorSeriesExponent(start, middle, exp, ref startPower, out factorial);
            startPower *= exp;      //we have startPower computed at "middle", need to multiply once here to compute it for the "middle + 1" exponent.
            RealNumber right = taylorSeriesExponent(middle + 1, end, exp, ref startPower, out secondFactorial);
            factorial *= secondFactorial;
            IntegerNumber.shrink(ref factorial);
            //left / A! + right / B! = [left * (B!/A!) + right] / B!
            return left * secondFactorial + right;
        }

        public IntegerNumber getBitBurst(int index, int count)
        {   //value = this.mantissa << this.shift
            IntegerNumber result = this.mantissa.GetBitsAt(index - (int)this.shift, count);
            return result;
        }

        /// <summary>
        /// Using the Taylor expansion with the bit burst algorithm from: https://members.loria.fr/PZimmermann/talks/arctan.pdf
        /// </summary>
        /// <returns></returns>
        public RealNumber GetTaylorExp()
        {
            if (this.IsZero)
            {
                return new RealNumber(1L, this.totalDigits);
            }

            int totalDigits = this.totalDigits + 2;
            bool isNegative = this.IsNegative;
            RealNumber log = RealNumber.Abs(this).ChangePrecision(totalDigits);

            double estimateLog2 = log.EstimateLog2;
            //Exp[n] = Exp[n / 2] ^ 2 ; -> Repeat until argument n is in (0 .. 0.25]
            int squarings = Math.Max(0, (int)Math.Floor(estimateLog2) + 2);
            log >>= squarings;

            int totalInputBits = Math.Min(totalDigits * 64, 2 + this.mantissa.BitsCount);
            RealNumber result = new RealNumber(1L, totalDigits);
            for (int count = 1; count < totalInputBits; count <<= 1)
            {
                int start = count << 1;
                IntegerNumber y = log.getBitBurst(-start, count);
                if (y.IsZero)
                {
                    continue;
                }
                RealNumber x = new RealNumber(-start, y, totalDigits);
                int factors = getNumberOfExpFactors(-count, totalDigits);
                //if (-1.0 * factors * count - getApproximativeLogFactorial(factors) * LG_E < -totalDigits * 64.0 ||
                //    -1.0 * (factors + 1) * count - getApproximativeLogFactorial(factors + 1) * LG_E > -totalDigits * 64.0)
                //{
                //    Debugger.Break();
                //}

                x = isNegative ? -x : x;
                IntegerNumber factorial;
                RealNumber partial = taylorSeriesExponent(1, factors, x, ref x, out factorial);
                result *= RealNumber.One + partial / factorial;
            }

            for (int i = squarings; --i >= 0;)
            {
                result = result.Square();
            }
            result = result.ChangePrecision(this.totalDigits);
            return result;
        }

        #endregion

        #region Exponential

        //https://en.wikipedia.org/wiki/Householder%27s_method
        //f[x_] := Log[x] - y
        //D[f[x], x]
        //D[f[x], {x, 2}]
        //FullSimplify[x - f[x]*D[f[x], x]/(D[f[x], x]^2 - f[x]*D[f[x], {x, 2}]/2)] = x (2 + y - Log[x])/(2 - y + Log[x]) = x*(2 - f[x])*(2 + f[x])
        public RealNumber GetExp()
        {   //2^x = exp(this) ==> exp(x * ln(2)) = exp(this) => x * ln(2) = this => x = this * log2(e) = this / ln(2)            
            RealNumber ln2 = GetLN2(this.totalDigits + 2);
            long k = (this / ln2).ToLong;
            RealNumber log = this - k * ln2;
            if (log.IsZero)
            {
                return new RealNumber(k, IntegerNumber.One, this.totalDigits);
            }
            //x = k + log ==> 2^x == 2^(k+log) == 2^k * 2^log = 2^k * exp(log * ln(2))
            RealNumber exp = new RealNumber(Math.Exp(log.ToDouble), 1 + 2);
            for (int bits = 52; bits < this.RealMaxBits; bits *= 3)
            {
                int stepPrecision = Math.Min(this.totalDigits + 2, (bits * 3 + 63) >> 6);
                exp.totalDigits = stepPrecision;
                RealNumber cache = exp.GetLog() - log.ChangePrecision(stepPrecision);
                exp *= (2 - cache) / (2 + cache);
            }
            RealNumber result = new RealNumber(exp.shift + k, exp.mantissa, this.totalDigits);
            return result;
        }

        private static void reduce(List<RealNumber> list, List<Func<ComplexNumber, ComplexNumber>> translations,
            ref RealNumber argument, RealNumber factor, Func<ComplexNumber, ComplexNumber> translation)
        {
            if ((factor.IsNegative && argument <= factor) || (!factor.IsNegative && argument >= factor))
            {
                argument -= factor;
                list.Add(factor);
                translations.Add(translation);
            }
        }

        //f[x_] := Log[x] - y
        //D[f[x], x] = 1/x
        //D[f[x], {x, 2}] = -1/x^2
        //FullSimplify[x - 2*f[x]*D[f[x], x]/(2*D[f[x], x]^2 - f[x]*D[f[x], {x, 2}])] = x (-1 + 4/(2 - y + Log[x])) = x*(2 - f[x])/(2 + f[x])
        public static ComplexNumber GetComplexExp(ComplexNumber argument)
        {
            int totalDigits = Math.Max(argument.Real.totalDigits, argument.Imaginary.totalDigits);
            RealNumber pi = GetPI(totalDigits + 2);
            RealNumber argDiv2PI = argument.Imaginary / (pi << 1);
            RealNumber imaginaryPart = (argDiv2PI - argDiv2PI.Floor()) * (pi << 1);
            List<RealNumber> list = new List<RealNumber>();
            List<Func<ComplexNumber, ComplexNumber>> functionsList = new List<Func<ComplexNumber, ComplexNumber>>();
            reduce(list, functionsList, ref imaginaryPart, -pi, x => -x);
            reduce(list, functionsList, ref imaginaryPart, -pi >> 1, x => x.DivI());
            reduce(list, functionsList, ref imaginaryPart, pi, x => -x);
            reduce(list, functionsList, ref imaginaryPart, pi >> 1, x => x.MulI());
            RealNumber ln2 = GetLN2(totalDigits + 2);
            long k = (argument.Real / ln2).ToLong;
            ComplexNumber log = new ComplexNumber(argument.Real - ln2 * k, imaginaryPart);
            double factor = Math.Exp(log.Real.ToDouble);
            double imaginary = imaginaryPart.ToDouble;
            ComplexNumber exp = new ComplexNumber(
                new RealNumber(factor * Math.Cos(imaginary), 2 + 1),
                new RealNumber(factor * Math.Sin(imaginary), 2 + 1));
            for (int bits = 52; bits < (totalDigits + 2) * 64; bits *= 3)
            {
                int stepPrecision = Math.Min(totalDigits + 2, (bits * 3 + 63) >> 6);
                exp.SetTotalDigits(stepPrecision);
                ComplexNumber cache = GetComplexLog(exp) - log.ChangePrecision(stepPrecision);
                exp *= (2 - cache) / (2 + cache);
            }
            ComplexNumber result = new ComplexNumber(
                new RealNumber(exp.Real.shift + k, exp.Real.mantissa, totalDigits),
                new RealNumber(exp.Imaginary.shift + k, exp.Imaginary.mantissa, totalDigits));
            for (int i = list.Count; --i >= 0;)
            {
                result = functionsList[i](result);
            }
            return result;
        }

        #endregion

        public RealNumber Floor()
        {
            if (this.shift >= 0)
            {
                return this;
            }
            if (this.shift < -this.mantissaBits - 1)
            {
                return RealNumber.Zero;
            }
            IntegerNumber mantissa = this.mantissa >> (-(int)this.shift);
            RealNumber result = new RealNumber(0L, mantissa, this.totalDigits);
            return result;
        }

        public RealNumber ChangePrecision(int newNumberOf64BitDigits)
        {
            RealNumber result = new RealNumber(this.shift, this.mantissa, newNumberOf64BitDigits);
            return result;
        }

        private static ConcurrentDictionary<int, RealNumber> piDictionary = new ConcurrentDictionary<int, RealNumber>();
        public static RealNumber GetPI(int numberOf64BitDigits)
        {
            RealNumber pi = piDictionary.GetOrAdd(numberOf64BitDigits, getPIValue);
            return pi;
        }

        private static RealNumber getPIValue(int numberOf64BitDigits)
        {
            RealNumber sqrt2 = new RealNumber(2L, numberOf64BitDigits + 2).GetSqrt();
            RealNumber one = new RealNumber(1L, numberOf64BitDigits + 2);
            RealNumber y = sqrt2 - 1;
            RealNumber a = 6 - (sqrt2 << 2);
            RealNumber error = one;
            for (int n = 1; ; n++)
            {
                RealNumber y2 = y.Square();
                RealNumber root4 = (one - y2.Square()).GetSqrt().GetSqrt();
                y = (one - root4) / (one + root4);
                y2 = (y + 1).Square();
                RealNumber save = y2;
                y2 = y2.Square();
                RealNumber nexta = y2 * a - (((save - y) * y) << (n * 2 + 1));
                RealNumber currentError = RealNumber.Abs(nexta - a);
                if (currentError >= error)
                {
                    RealNumber inverse = a.Inverse();
                    RealNumber result = new RealNumber(inverse.shift, inverse.mantissa, numberOf64BitDigits);
                    return result;
                }
                a = nexta;
                error = currentError;
            }
        }

        private static ConcurrentDictionary<int, RealNumber> ln2Dictionary = new ConcurrentDictionary<int, RealNumber>();
        public static RealNumber GetLN2(int numberOf64BitDigits)
        {
            RealNumber ln2 = ln2Dictionary.GetOrAdd(numberOf64BitDigits, getLN2Value);
            return ln2;
        }

        public static readonly RealNumber One = new RealNumber(0, IntegerNumber.One, 1);
        public static readonly RealNumber Infinity = new RealNumber(long.MaxValue, IntegerNumber.One, 1);
        //public static readonly RealNumber PI = GetPI();
        //public static readonly RealNumber LN2 = GetLN2();
        //public static readonly RealNumber LN10 = new RealNumber(0, new IntegerNumber(10), ).Log;

        public static RealNumber AGMReal(RealNumber x, RealNumber y)
        {
            RealNumber previousError = RealNumber.Abs(x - y);
            while (true)
            {
                RealNumber am = (x + y) >> 1;
                RealNumber gm = (x * y).GetSqrt();
                RealNumber error = RealNumber.Abs(am - gm);
                if (error >= previousError)
                {
                    return am;
                }
                x = am;
                y = gm;
                previousError = error;
            }
        }

        public static ComplexNumber AGMComplex(ComplexNumber x, ComplexNumber y)
        {
            ComplexNumber am, gm;
            RealNumber error, previousError = (x - y).Energy;
            while (true)
            {
                am = (x + y) >> 1;
                gm = (x * y).GetSqrtWithPositiveRealPart();
                //https://maths-people.anu.edu.au/~brent/pd/RNC7t.pdf
                //https://books.google.ro/books?id=-8wuH5AwbwMC&pg=PA163&lpg=PA163&dq=complex+agm+sin+cos&source=bl&ots=tY80dAM7T3&sig=Il_GHgnl_rbFLIz4gLTFUK7XCJM&hl=en&sa=X&ved=2ahUKEwiPiNXtldncAhUJyqQKHT2FC40Q6AEwBnoECAkQAQ#v=onepage&q=complex%20agm%20sin%20cos&f=false
                error = (am - gm).Energy;
                if (error >= previousError)
                {
                    break;
                }
                x = am;
                y = gm;
                previousError = error;
            }
            return am;
        }

        private static void shiftRightInBase10(int keepDigits, ref IntegerNumber number, ref long shift10)
        {
            if (number.IsZero)
            {
                return;
            }
            int shift = (int)Math.Floor(number.BitsCount * LOG10_2) - keepDigits;
            number = shift >= 0 ? number / IntegerNumber.Pow(10, shift) : number * IntegerNumber.Pow(10, -shift);
            shift10 += shift;
            IntegerNumber.shrink(ref number);
        }

        private static IntegerNumber mostSignificantExponent(IntegerNumber number, int keepDigits, long power, out long shift10)
        {
            shift10 = 0;
            if (power == 0)
            {
                return IntegerNumber.One;
            }
            int bit = IntegerNumber.BSR((ulong)power);
            IntegerNumber result = number;
            while (--bit >= 0)
            {
                result *= result;
                shift10 *= 2;
                shiftRightInBase10(keepDigits, ref result, ref shift10);
                if ((power & (1L << bit)) != 0)
                {
                    result *= number;
                    shiftRightInBase10(keepDigits, ref result, ref shift10);
                }
            }
            return result;
        }

        private string getStringRepresentation()
        {
            long shift = this.shift, shift10;
            IntegerNumber scale;
            int keepDigits = (int)Math.Ceiling(this.RealMaxBits * LOG10_2);
            if (shift >= 0)
            {
                int bucketBits = RealMaxBits;
                IntegerNumber start = IntegerNumber.One << bucketBits;
                long exponent = shift / bucketBits;
                scale = mostSignificantExponent(start, keepDigits, exponent, out shift10);
                scale <<= (int)(shift % bucketBits);
                shiftRightInBase10(keepDigits, ref scale, ref shift10);
            }
            else
            {
                shift = -shift;
                int bucketFives = (int)Math.Ceiling(LOG5_2 * RealMaxBits);
                IntegerNumber start = IntegerNumber.Pow(5, bucketFives);
                long exponent = shift / bucketFives;
                scale = mostSignificantExponent(start, keepDigits, exponent, out shift10);
                scale *= IntegerNumber.Pow(5, (int)(shift % bucketFives));
                shiftRightInBase10(keepDigits, ref scale, ref shift10);
                shift10 -= shift;
            }
            IntegerNumber mantissa = scale * this.mantissa;
            int digits = (int)Math.Floor(mantissa.BitsCount * LOG10_2) - (this.ToStringDigits + 1);
            shift10 += digits;
            mantissa = digits >= 0 ? mantissa / IntegerNumber.Pow(10, digits) : mantissa * IntegerNumber.Pow(10, -digits);

            long rem;
            IntegerNumber div = IntegerNumber.DivRem(mantissa, 1000000000L, out rem);
            mantissa = div + (rem <= -500000000L ? IntegerNumber.MinusOne : (rem >= 500000000L ? IntegerNumber.One : IntegerNumber.Zero));   //Rounding in base 10.
            shift10 += 9;

            string result = mantissa.ToString();
            int hizeroes = 0;
            while (hizeroes + 1 < result.Length && result[result.Length - 1 - hizeroes] == '0')
            {
                hizeroes++;
            }
            if (hizeroes > 0)
            {
                result = result.Remove(result.Length - hizeroes);
                shift10 += hizeroes;
            }

            int dotIndex = mantissa.IsNegative ? 2 : 1;
            shift10 += result.Length - dotIndex;
            if (shift10 >= 0 && shift10 <= 20)
            {
                dotIndex += (int)shift10;
                shift10 = 0;
            }
            else if (shift10 >= -20 && shift10 < 0)
            {
                result = result.Insert(dotIndex - 1, "".PadRight(-(int)shift10, '0'));
                shift10 = 0;
            }
            if (dotIndex >= result.Length)
            {
                result = result + "".PadRight(dotIndex - result.Length, '0');
            }
            else
            {
                result = result.Insert(dotIndex, ".");
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

        public int GetPrecisionDigits()
        {
            return this.totalDigits;
        }
    }

    [DebuggerDisplay("{ToString()}")]
    public struct ComplexNumber
    {
        public static readonly ComplexNumber Zero = new ComplexNumber(RealNumber.Zero, RealNumber.Zero);
        public static readonly ComplexNumber One = new ComplexNumber(RealNumber.One, RealNumber.Zero);

        public static implicit operator ComplexNumber(RealNumber real)
        {
            return new ComplexNumber(real, new RealNumber(0L, real.totalDigits));
        }

        public static implicit operator ComplexNumber(long real)
        {
            return new ComplexNumber(new RealNumber(real, 2), RealNumber.Zero);
        }

        public RealNumber Real { get; set; }
        public RealNumber Imaginary { get; set; }

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

        public static ComplexNumber operator *(RealNumber x, ComplexNumber y)
        {
            return new ComplexNumber(x * y.Real, x * y.Imaginary);
        }

        public ComplexNumber Square()
        {
            return new ComplexNumber((this.Real + this.Imaginary) * (this.Real - this.Imaginary), this.Real * this.Imaginary << 1);
        }

        public static ComplexNumber operator *(ComplexNumber x, RealNumber y)
        {
            return new ComplexNumber(x.Real * y, x.Imaginary * y);
        }

        public static ComplexNumber operator /(ComplexNumber x, RealNumber y)
        {
            return new ComplexNumber(x.Real / y, x.Imaginary / y);
        }

        public static ComplexNumber operator /(RealNumber x, ComplexNumber y)
        {
            return y.Conjugate * (x * y.Energy.Inverse());
        }

        public static ComplexNumber operator /(ComplexNumber x, ComplexNumber y)
        {
            return x * y.Conjugate * y.Energy.Inverse();
        }

        public ComplexNumber Conjugate { get { return new ComplexNumber(this.Real, -this.Imaginary); } }

        public static ComplexNumber operator +(ComplexNumber x, ComplexNumber y)
        {
            return new ComplexNumber(x.Real + y.Real, x.Imaginary + y.Imaginary);
        }

        public static ComplexNumber operator >>(ComplexNumber x, int shift)
        {
            return new ComplexNumber(x.Real >> shift, x.Imaginary >> shift);
        }

        public static ComplexNumber operator <<(ComplexNumber x, int shift)
        {
            return new ComplexNumber(x.Real << shift, x.Imaginary << shift);
        }

        public static ComplexNumber operator -(ComplexNumber x, ComplexNumber y)
        {
            return new ComplexNumber(x.Real - y.Real, x.Imaginary - y.Imaginary);
        }

        public static ComplexNumber operator -(ComplexNumber x)
        {
            return new ComplexNumber(-x.Real, -x.Imaginary);
        }

        //* Calculates the square root of this object.
        //Adapted from Numerical Recipes in C - The Art of Scientific Computing
        public ComplexNumber GetSqrtWithPositiveRealPart()
        {
            if (this.Real.IsZero && this.Imaginary.IsZero)
            {
                return this;
            }
            RealNumber absRe = RealNumber.Abs(this.Real);
            RealNumber absIm = RealNumber.Abs(this.Imaginary);
            RealNumber r, w;
            if (absRe >= absIm)
            {
                r = absIm / absRe;
                w = absRe.GetSqrt() * ((RealNumber.One + (RealNumber.One + r.Square()).GetSqrt()) >> 1).GetSqrt();
            }
            else
            {
                r = absRe / absIm;
                w = absIm.GetSqrt() * ((r + (RealNumber.One + r.Square()).GetSqrt()) >> 1).GetSqrt();
            }

            ComplexNumber result;
            if (!this.Real.IsNegative)
            {
                result = new ComplexNumber(w, this.Imaginary / w >> 1);
            }
            else
            {
                result = new ComplexNumber(absIm / w >> 1, this.Imaginary.IsNegative ? -w : w);
            }
            if (result.Real.IsNegative)
            {
                result = -result;
            }
            return result;
        }

        public RealNumber Energy { get { return this.Real.Square() + this.Imaginary.Square(); } }
        public RealNumber Magnitude { get { return this.Energy.GetSqrt(); } }

        public static RealNumber Abs(ComplexNumber argument)
        {
            return argument.Magnitude;
        }

        public ComplexNumber ComplexLog()
        {
            ComplexNumber result = RealNumber.GetComplexLog(this);
            return result;
        }

        public ComplexNumber ComplexExp()
        {
            ComplexNumber result = RealNumber.GetComplexExp(this);
            return result;
        }

        public override string ToString()
        {
            return "{" + this.Real.ToString() + " ; " + this.Imaginary.ToString() + "}";
        }

        public void SetTotalDigits(int totalDigits)
        {
            this.Real = this.Real.ChangePrecision(totalDigits);
            this.Imaginary = this.Imaginary.ChangePrecision(totalDigits);
        }

        public ComplexNumber ChangePrecision(int newNumberOf64BitDigits)
        {
            ComplexNumber result = new ComplexNumber(
                this.Real.ChangePrecision(newNumberOf64BitDigits),
                this.Imaginary.ChangePrecision(newNumberOf64BitDigits));
            return result;
        }

        public static ComplexNumber FromPolarAngle(RealNumber angle)
        {
            return RealNumber.GetComplexExp(new ComplexNumber(RealNumber.Zero, angle));
        }

        public int GetPrecisionDigits()
        {
            return Math.Max(this.Real.GetPrecisionDigits(), this.Imaginary.GetPrecisionDigits());
        }

        public ComplexNumber Inverse()
        {
            return this.Conjugate * this.Energy.Inverse();
        }

        public ComplexNumber ToPower(long power)
        {
            if (power == 0 || power == 1)
            {
                return power == 0 ? ComplexNumber.One : this;
            }
            ComplexNumber result = this.ToPower(power >> 1).Square();
            if ((power & 1L) != 0)
            {
                result *= this;
            }
            return result;
        }
    }

    public static class RealNumbersUnitTest
    {
        public static bool UnitTest(int seed)
        {
            bool ok = true;
            int totalDigits = 2048 / 64;
            RealNumber log0 = new RealNumber(1L << 62, totalDigits);
            RealNumber exp1 = log0.GetExp();
            RealNumber log1 = exp1.GetLog();
            ok &= (log1 - log0).IsZero;

            RealNumber n_2 = RealNumber.GetLN2(totalDigits).GetExp();
            RealNumber n_1a = new RealNumber(1.0, totalDigits).GetExp().GetLog();
            RealNumber n_1b = new RealNumber(10.0, totalDigits).GetLog().GetExp();
            var s1 = n_1b.ToString();
            RealNumber error = RealNumber.Max(RealNumber.Max(RealNumber.Abs(n_2 - 2), RealNumber.Abs(n_1a - 1)), RealNumber.Abs(n_1b - 10));
            var s2 = error.ToString();
            RealNumber tolerance = RealNumber.One >> 2000;
            ok &= error < tolerance;
            var s3 = new RealNumber(10000, IntegerNumber.Pow(3, 10000), totalDigits).ToString();
            ok &= s3.StartsWith("3.25464658549366205883184429132") && s3.EndsWith("+7781");
            var s4 = new RealNumber(-25000, IntegerNumber.Pow(1, 1), totalDigits).ToString();
            ok &= s4.StartsWith("1.778723326301851865925552009580") && s4.EndsWith("-7526");
            return ok;
        }

        public static bool NumericIntegerOperationsUnitTest()
        {
            Random random = new Random(1001);
            byte[] data = new byte[257];
            bool ok = true;
            for (int step = 2000; --step >= 0;)
            {
                random.NextBytes(data);
                data[data.Length - 1] = 0;
                IntegerNumber int10 = new IntegerNumber(data);
                IntegerNumber.shrink(ref int10);
                IntegerNumber test1 = int10.Sqrt();
                IntegerNumber test2 = test1.Square();
                ok &= 0 <= int10 - test2 && int10 - test2 <= test1 * 2;
                IntegerNumber t1 = int10.Inverse();
                IntegerNumber t2 = (IntegerNumber.One << (int10.Digits * 128)) / int10;
                ok &= t1 == t2;
                IntegerNumber t3 = int10.InverseSqrt();
                IntegerNumber t4 = ((IntegerNumber.One << (int10.Digits * 64 * 3)) / int10).Sqrt();
                IntegerNumber t5 = ((int10 << (int10.Digits * 64)).Inverse()).Sqrt();
                ok &= t3 == t4 && t4 == t5;
            }
            return ok;
        }

        public static bool SqrtUnitTest()
        {
            RealNumber ln10 = new RealNumber(10L, 4).GetLog();
            List<string> representation = new List<string>();
            for (int power = -25; power <= 25; power++)
            {
                RealNumber rn1 = (ln10 * power).GetExp();
                RealNumber rn2 = (ln10 * power).GetTaylorExp();
                if (!(rn1 - rn2).IsZero)
                {
                    return false;
                }
                double dn = Math.Pow(10, power);
                for (int i = 0; i <= 100; i++)
                {
                    RealNumber number1 = rn1 * i / 100L;
                    RealNumber number2 = rn2 * i / 100L;
                    double num1 = dn * i / 100.0;
                    double num2 = double.Parse(number1.ToString());
                    while (num1.ToString("R").Contains("9999"))
                    {
                        num1 = BitConverter.Int64BitsToDouble(BitConverter.DoubleToInt64Bits(num1) + 1L);
                    }
                    while (num1.ToString("R").Contains("0000"))
                    {
                        num1 = BitConverter.Int64BitsToDouble(BitConverter.DoubleToInt64Bits(num1) - 1L);
                    }
                    while (num2.ToString("R").Contains("9999"))
                    {
                        num2 = BitConverter.Int64BitsToDouble(BitConverter.DoubleToInt64Bits(num2) + 1L);
                    }
                    while (num2.ToString("R").Contains("0000"))
                    {
                        num2 = BitConverter.Int64BitsToDouble(BitConverter.DoubleToInt64Bits(num2) - 1L);
                    }
                    if (num1 != num2)
                    {
                        return false;
                    }
                }
            }

            IntegerNumber zero1 = new IntegerNumber(0L);
            RealNumber zero2 = new RealNumber(0L, 100);

            bool ok = zero1.ToString() == "0" && zero2.ToString() == "0";
            if (!ok)
            {
                return false;
            }
            int qwords = 1 * 256;
            RealNumber maxError = new RealNumber(-(qwords * 64 - 16), IntegerNumber.One, 1);
            for (int i = 2; i < 100; i++)
            {
                RealNumber x = new RealNumber((long)i, qwords);
                RealNumber sqrtX = x.GetSqrt();
                RealNumber invSqrtX = x.GetInverseSqrt();
                RealNumber one = sqrtX * invSqrtX;
                RealNumber zeroA = sqrtX.Square() - x;
                RealNumber zeroB = invSqrtX.Square().Inverse() - x;
                RealNumber error = RealNumber.Max(RealNumber.Abs(one - 1), RealNumber.Max(RealNumber.Abs(zeroA), RealNumber.Abs(zeroB)));
                if (error > maxError)
                {
                    return false;
                }
            }
            return true;
        }

        public static bool ComplexArcTanSinCosUnitTest()
        {
            RealNumber epsilon = 1E-34;
            ComplexNumber test1 = new ComplexNumber(-100, -1E-10);
            ComplexNumber logt1 = test1.ComplexLog();
            ComplexNumber tres1 = logt1.ComplexExp();
            if (RealNumber.Abs(tres1.Real - test1.Real) > epsilon || RealNumber.Abs(tres1.Imaginary - test1.Imaginary) > epsilon || logt1.Imaginary > 0)
            {
                return false;
            }

            ComplexNumber test2 = new ComplexNumber(-100, 1E-10);
            ComplexNumber logt2 = test2.ComplexLog();
            ComplexNumber tres2 = logt2.ComplexExp();
            if (RealNumber.Abs(tres2.Real - test2.Real) > epsilon || RealNumber.Abs(tres2.Imaginary - test2.Imaginary) > epsilon || logt2.Imaginary < 0)
            {
                return false;
            }

            Random random = new Random(1001);
            int bytes = 100;
            int precisionQWords = (bytes + 7) / 8 + 10;
            RealNumber maxError = new RealNumber(-precisionQWords * 64 + 64, IntegerNumber.One, 2);

            for (int step = 30; --step >= 0;)
            {
                byte[] data1 = new byte[bytes];
                byte[] data2 = new byte[bytes];
                random.NextBytes(data1);
                random.NextBytes(data2);
                IntegerNumber integ1 = new IntegerNumber(data1);
                IntegerNumber integ2 = new IntegerNumber(data2);
                RealNumber arg1 = new RealNumber(-bytes * 8 + random.Next(32), integ1, precisionQWords);
                RealNumber arg2 = new RealNumber(-bytes * 8 + random.Next(32), integ2, precisionQWords);
                ComplexNumber comp1 = new ComplexNumber(arg1, arg2);

                ComplexNumber log = comp1.ComplexLog();
                ComplexNumber exp = log.ComplexExp();

                RealNumber error = RealNumber.Max(
                    RealNumber.Abs(RealNumber.One - comp1.Real / exp.Real),
                    RealNumber.Abs(RealNumber.One - comp1.Imaginary / exp.Imaginary));
                if (error > maxError)
                {
                    return false;
                }
            }
            return true;
        }
    }
}
