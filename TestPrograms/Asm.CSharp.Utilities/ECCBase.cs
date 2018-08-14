using Bitsy.Core;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Numerics;
using System.Text;

namespace Utilities
{
    //http://www.secg.org/sec2-v2.pdf
    public abstract class ECCBase
    {
        //EC-function: Y^2 = X^3 + a * X + b 'modulo' P
        public abstract int N { get; }
        public abstract int BitsCount { get; }
        public abstract int BytesCount { get; }
        public abstract int KaratsubaBufferN { get; }

        public abstract void GetModFriendlyPrime(ulong[] inputOutput2x);
        public abstract void RemodFriendlyPrime(ulong[] inputOutput1x);

        public abstract ulong[] P { get; }
        public abstract ulong[] A { get; }
        public abstract ulong[] B { get; }
        public abstract ulong[] One { get; }

        public abstract ulong[] W1ModOrder { get; }
        public abstract ulong[] W2ModOrder { get; }
        public abstract ulong[] W1ModP { get; }
        public abstract ulong[] W2ModP { get; }

        public abstract OwnECPoint G { get; }
        public abstract ulong[] Order { get; }
        public abstract int H { get; }
        public abstract bool SupportsEndomorphism { get; }

        public abstract JacobianECPoint InfinityJacobian { get; }
        public abstract OwnECPoint Infinity { get; }

        public abstract bool AIsZero { get; }
        public abstract bool AIsMinus3 { get; }

        private ulong[] temporary2x;
        private ulong[] temporaryKaratsubaBuffer;
        private ulong[] x1;
        private ulong[] y1;
        private ulong[] S1;
        private ulong[] S2;
        private ulong[] Z2;
        private ulong[] Z3;
        private ulong[] R1;
        protected ulong[] specialComparer;

        public Action<ulong[], ulong[], ulong[]> AdditionModP;
        public Action<ulong[], ulong[], ulong> AddScaledModP;
        public Action<ulong[], ulong[], ulong[]> SubtractModP;
        public Action<ulong[], ulong[], ulong> SubScaledModP;

        public void AddModP(ulong[] destination, ulong[] source)
        {
            this.AdditionModP(destination, source, destination);
        }
        public void SubModP(ulong[] destination, ulong[] source)
        {
            this.SubtractModP(destination, source, destination);
        }

        public ECCBase()
        {
            this.temporary2x = new ulong[N * 2];
            this.temporaryKaratsubaBuffer = new ulong[KaratsubaBufferN];
            this.x1 = new ulong[N];
            this.y1 = new ulong[N];
            this.S1 = new ulong[N];
            this.S2 = new ulong[N];
            this.Z2 = new ulong[N];
            this.Z3 = new ulong[N];
            this.R1 = new ulong[N];
            this.specialComparer = new ulong[N];
        }

        public static int bsr(ulong x)
        {
            int result = 0;
            for (int shift = 32; shift > 0; shift >>= 1)
            {
                if (x >= (1UL << shift))
                {
                    x >>= shift;
                    result += shift;
                }
            }
            return result;
        }

        public static int HighBit(ulong[] number)
        {
            for (int i = number.Length; --i >= 0;)
            {
                if (number[i] == 0)
                {
                    continue;
                }
                return i * 64 + bsr(number[i]);
            }
            return -1;
        }

        public void SquareModP(ulong[] inputOutput)
        {
            AsmX64Operations.KaratsubaSquare(inputOutput, temporary2x, N, temporaryKaratsubaBuffer);
            this.GetModFriendlyPrime(temporary2x);
            for (int i = N; --i >= 0;)
            {
                inputOutput[i] = temporary2x[i];
            }
        }

        public void MultiplyModP(ulong[] inputOutput, ulong[] input2)
        {
            AsmX64Operations.Karatsuba(inputOutput, input2, temporary2x, N, temporaryKaratsubaBuffer);
            this.GetModFriendlyPrime(temporary2x);
            for (int i = N; --i >= 0;)
            {
                inputOutput[i] = temporary2x[i];
            }
        }

        public static void Let(ulong[] destination, ulong[] source)
        {
            for (int i = destination.Length; --i >= 0;)
            {
                destination[i] = source[i];
            }
        }

        public OwnECPoint Negative(OwnECPoint point)
        {
            OwnECPoint result = point.Copy();
            SubtractModP(this.P, result.Y, result.Y);
            return result;
        }

        public bool IsZero(ulong[] number)
        {
            this.RemodFriendlyPrime(number);
            for (int i = N; --i >= 0;)
            {
                if (number[i] != 0)
                {
                    return false;
                }
            }
            return true;
        }

        public bool IsEqual(ulong[] number1, ulong[] number2)
        {
            SubtractModP(number1, number2, this.specialComparer);
            return IsZero(this.specialComparer);
        }

        public void Subtract(ref OwnECPoint p, OwnECPoint q)
        {
            Add(ref p, q, true);
        }

        public void Subtract(ref JacobianECPoint p, OwnECPoint q)
        {
            Add(ref p, q, true);
        }

        public bool Verify(OwnECPoint point)
        {   //Y^2 = X^3 + a*X + b 'modulo' P
            if (point.IsInfinity)
            {
                return true;
            }
            Let(x1, point.X);
            Let(y1, point.Y);
            SquareModP(y1);
            SquareModP(x1);
            AddModP(x1, A);
            MultiplyModP(x1, point.X);
            AddModP(x1, B);
            SubModP(x1, y1);
            return IsZero(x1);
        }

        //https://en.wikibooks.org/wiki/Cryptography/Prime_Curve/Jacobian_Coordinates#Mixed_Addition_(with_affine_point)_(8M_+_3S)
        //U1 = X1  [*Z2^2]
        //U2 = X2*Z1^2
        //S1 = Y1  [*Z2^3]
        //S2 = Y2*Z1^3
        //if (U1 == U2)
        //  if (S1 != S2)
        //    return POINT_AT_INFINITY
        //  else 
        //    return POINT_DOUBLE(X1, Y1, Z1)
        //H = U2 - U1
        //R = S2 - S1
        //X3 = R^2 - H^3 - 2*U1*H^2
        //Y3 = R*(U1*H^2 - X3) - S1*H^3
        //Z3 = H*Z1*Z2
        //return (X3, Y3, Z3)

        //http://www.hyperelliptic.org/EFD/oldefd/jacobian.html#mADD
        //Z2 = [p.Z]^2
        //S1 = [q.X] * Z2 - [p.X]
        //R1 = [+/-q.Y] * [p.Z] * Z2 - [p.Y]
        //S2 = S1^2
        //Z3 = ([p.Z]+S1)^2 - S2 - Z2
        //S2 *= 4
        //S1 *= S2
        //R1 = 2*R1
        //Y3 = [p.X]*S2
        //X3 = R1^2-S1-2*Y3
        //Y3 = R1*(Y3-X3)-2*[p.y]*S1
        public void Add(ref JacobianECPoint p, OwnECPoint q, bool subtract = false)
        {   //Mixed addition returns p + q in "7M + 4S"
            if (p.IsInfinity || q.IsInfinity)
            {
                if (!q.IsInfinity)
                {
                    p.SetFrom(q.X, q.Y, One, q.IsInfinity);
                    if (subtract)
                    {
                        SubtractModP(P, p.Y, p.Y);
                    }
                }
                return;
            }

            Let(Z2, p.Z); SquareModP(Z2);           //Z2 = [p.Z]^2
            Let(S1, q.X); MultiplyModP(S1, Z2); SubModP(S1, p.X);   //S1 = [q.X] * Z2 - [p.X]
            Let(R1, q.Y);
            if (subtract)
            {
                SubtractModP(P, R1, R1);
            }
            MultiplyModP(R1, p.Z); MultiplyModP(R1, Z2); SubModP(R1, p.Y);   //R1 = [+/-q.Y] * [p.Z] * Z2 - [p.Y]

            if (IsZero(S1))
            {
                if (IsZero(R1))
                {   //now: p == q
                    GetDouble(ref p);
                    return;
                }
                p.SetInfinity();
                return;
            }

            Let(S2, S1); SquareModP(S2);        //S2 = S1^2
            Let(Z3, p.Z); AddModP(Z3, S1); SquareModP(Z3); SubModP(Z3, S2); SubModP(Z3, Z2);    //Z3 = ([p.Z]+S1)^2 - S2 - Z2

            AddScaledModP(S2, S2, 3);     //S2 *= 4
            MultiplyModP(S1, S2);   //S1 *= S2
            AddModP(R1, R1);        //R1 = 2*R1
            Let(y1, p.X); MultiplyModP(y1, S2);        //Y3 = [p.X]*S2
            Let(x1, R1); SquareModP(x1); SubModP(x1, S1); SubScaledModP(x1, y1, 2);        //X3 = R1^2-S1-2*Y3
            MultiplyModP(S1, p.Y); SubModP(y1, x1); MultiplyModP(y1, R1); SubScaledModP(y1, S1, 2);   //Y3 = R1*(Y3-X3)-2*[p.y]*S1

            p.SetFrom(x1, y1, Z3, false);
        }

        public void Add(ref OwnECPoint p, OwnECPoint q, bool subtract = false)
        {   //returns: P + Q
            if (p.IsInfinity || q.IsInfinity)
            {
                if (!q.IsInfinity)
                {
                    p.SetFrom(q.X, q.Y, q.IsInfinity);
                    if (subtract)
                    {
                        SubtractModP(P, p.Y, p.Y);
                    }
                }
                return;
            }

            SubtractModP(p.X, q.X, S2); //s2 = Px - Qx
            if (subtract)
            {   //s1 = Py - Qy; in case of subtract real Q.Y is -Q.Y so we add the value modulo P.
                AdditionModP(p.Y, q.Y, S1); //s1 = Py - Qy
            }
            else
            {
                SubtractModP(p.Y, q.Y, S1);
            }

            if (IsZero(S2))
            {   //Px == Qx
                if (IsZero(S1))
                {   // P == Q
                    GetDouble(ref p);
                    return;
                }
                p.SetInfinity();
                return;
            }

            this.RemodFriendlyPrime(S1);
            this.RemodFriendlyPrime(S2);
            AsmX64Operations.DivideModuloPrime(S1, S2, P, N);        //S = s1 = s1 / s2 'modulo' P.

            Let(R1, S1);
            SquareModP(S1);
            SubModP(S1, p.X);
            SubModP(S1, q.X);      //s1 = Result.x = S * S - Px - Qx

            SubtractModP(p.X, S1, S2);      //s2 = Px - Rx
            MultiplyModP(S2, R1);
            SubModP(S2, p.Y);      //s2 = S(Px - Rx) - Py

            p.SetFrom(S1, S2, false);
        }

        //      ec-function: Y^2 = X^3 + a*X + b 'modulo' P
        //P = {Xp, Yp}; Q = {Xq, Yq}
        //ADD:    S = (Yp - Yq) / (Xp - Xq)
        //        RESULT = {Xr = S * S - Xp - Xq; Yr = S(Xp - Xr) - Yp}
        //DUB:    S = (3 * Xp * Xp + a) / (2 * Yp)
        //        RESULT = {Xr = S * S - 2 * Xp ; Yr = S(Xp - Xr) - Yp}
        public void GetDouble(ref OwnECPoint p)
        {   //return: 2 * p
            if (p.IsInfinity)
            {
                return;
            }
            if (IsZero(p.Y))
            {
                p.SetInfinity();
                return;
            }

            Let(S1, p.X);
            SquareModP(S1); AddScaledModP(S1, S1, 2);
            AddModP(S1, A);
            Let(S2, p.Y);
            AddModP(S2, S2);

            this.RemodFriendlyPrime(S1);
            this.RemodFriendlyPrime(S2);
            AsmX64Operations.DivideModuloPrime(S1, S2, P, N);

            Let(R1, S1);
            SquareModP(S1);
            SubScaledModP(S1, p.X, 2);    //S1 = Rx = S * S - 2 * Px

            SubtractModP(p.X, S1, S2);      //S2 = Px - Rx
            MultiplyModP(S2, R1);
            SubModP(S2, p.Y);       //S2 = S(Px - Rx) - Py

            p.SetFrom(S1, S2, false);
        }

        public void GetDouble(ref JacobianECPoint p)
        {   //return: 2 * p in "3M + 4S" for a = 0 and "3M + 5S" for a = -3.
            if (p.IsInfinity)
            {
                return;
            }
            if (IsZero(p.Y))
            {
                p.SetInfinity();
                return;
            }

            if (AIsZero)
            {
                //https://en.wikibooks.org/wiki/Cryptography/Prime_Curve/Jacobian_Coordinates
                //if (Y == 0)
                //  return POINT_AT_INFINITY
                //S = 4 * X * Y^2
                //M = 3*X^2 [+ a* Z^4]
                //X' = M^2 - 2*S
                //Y' = M*(S - X') - 8*Y^4
                //Z' = 2*Y*Z = (Y+Z)^2 - Y^2 - Z^2
                //return (X', Y', Z')
                Let(R1, p.Y); SquareModP(R1); AddModP(R1, R1);          //R1 = 2 * Y^2
                Let(Z3, R1); SquareModP(Z3);                            //Z3 = 4 * Y^4
                Let(S1, p.X); MultiplyModP(S1, R1); AddModP(S1, S1);    //S1 = S = 4 * X * Y^2

                Let(S2, p.X); SquareModP(S2); AddScaledModP(S2, S2, 2);       //S2 = M = 3*X^2 [+ a* Z^4]   because a == 0.

                Let(x1, S2); SquareModP(x1); SubScaledModP(x1, S1, 2);                            //X' = M^2 - 2*S
                Let(y1, S1); SubModP(y1, x1); MultiplyModP(y1, S2); SubScaledModP(y1, Z3, 2);     //Y' = M*(S - X') - 8*Y^4
                Let(Z3, p.Y); AddModP(Z3, Z3); MultiplyModP(Z3, p.Z);                       //Z' = 2*Y*Z

                p.SetFrom(x1, y1, Z3, false);
            }
            else if (AIsMinus3)
            {
                //http://www.hyperelliptic.org/EFD/oldefd/jacobian.html#a3DBL
                //delta:=Z1^2;
                //gamma:=Y1^2;
                //beta:=4*X1*gamma;
                //alpha:=3*(X1-delta)*(X1+delta);
                //X3:=alpha^2-2*beta;
                //Z3:=(Y1+Z1)^2-gamma-delta;
                //Y3:=alpha*(beta-X3)-8*gamma^2;
                Let(Z2, p.Z); SquareModP(Z2);           //Z2:= delta:=Z1^2;
                Let(y1, p.Y); SquareModP(y1);           //y1:= gamma:=Y1^2;
                Let(R1, p.X); MultiplyModP(R1, y1); AddScaledModP(R1, R1, 3);     //R1:= beta:=4*X1*gamma;

                SubtractModP(p.X, Z2, S1);
                AdditionModP(p.X, Z2, S2);
                MultiplyModP(S1, S2);
                AddScaledModP(S1, S1, 2);     //S1:= alpha:=3*(X1-delta)*(X1+delta)       because a == -3

                AddModP(Z2, y1);    //delta = delta + gamma
                Let(x1, S1); SquareModP(x1); SubScaledModP(x1, R1, 2);    //X3:= alpha^2-2*beta;
                Let(Z3, p.Y); AddModP(Z3, p.Z); SquareModP(Z3); SubModP(Z3, Z2);   //Z3 = (Y1+Z1)^2-delta
                SubModP(R1, x1); MultiplyModP(R1, S1); SquareModP(y1); SubScaledModP(R1, y1, 8);     //Y3 = alpha*(beta-X3)-8*gamma^2

                p.SetFrom(x1, R1, Z3, false);
            }
            else
            {
                throw new NotImplementedException();
            }
        }

        //return: P * factor on the elliptic curves domain.
        public OwnECPoint MultiplyReference(OwnECPoint p, ulong[] factor)
        {
            int high = HighBit(factor);
            if (high < 0)
            {
                return new OwnECPoint(p.X, p.Y, true);
            }
            OwnECPoint result = p.Copy();
            for (int i = high; --i >= 0;)
            {
                GetDouble(ref result);
                if ((factor[i >> 6] & (1UL << (i & 63))) != 0)
                {
                    Add(ref result, p);
                }
            }
            RemodFriendlyPrime(result.X);
            RemodFriendlyPrime(result.Y);
            return result;
        }

        //return: P * factor on the elliptic curves domain.
        public OwnECPoint Multiply(OwnECPoint p, ulong[] factor)
        {
            ulong[] exp = new ulong[factor.Length + 1];
            ulong[] exp3 = new ulong[factor.Length + 1];
            factor.CopyTo(exp, 0);
            exp3[factor.Length] = AsmX64Operations.MultiplyDigitAndAdd(exp, 3, exp3, factor.Length);
            for (int i = exp.Length; --i >= 0;)
            {
                exp3[i] ^= exp[i];  //exp3= (exp*3) ^ exp
            }
            int high = HighBit(exp3);
            if (high < 0)
            {
                return new OwnECPoint(p.X, p.Y, true);
            }
            JacobianECPoint x = new JacobianECPoint(p.X, p.Y, true);
            ulong mask = 1UL << (high & 63);
            for (int i = high; --i >= 0;)
            {
                GetDouble(ref x);
                if ((exp3[(i + 1) >> 6] & mask) != 0)
                {
                    if ((exp[(i + 1) >> 6] & mask) != 0)
                    {
                        Subtract(ref x, p);
                    }
                    else
                    {
                        Add(ref x, p);
                    }
                }
                mask = (mask >> 1) | (mask << 63);
            }
            OwnECPoint result = x.ToECPoint(this);
            return result;
        }

        public static List<Equation> SplitNumberByExtendedEuclid(BigInteger x, BigInteger p, BigInteger quitBound)
        {
            //rules:
            //a*original(x) + b*original(p) = x
            //c*original(x) + d*original(p) = p
            BigInteger a = BigInteger.One;
            BigInteger b = BigInteger.Zero;
            BigInteger c = BigInteger.Zero;
            BigInteger d = BigInteger.One;
            List<Equation> result = new List<Equation>();
            result.Add(new Equation(a, b, x));
            while (x > quitBound)
            {
                BigInteger quotient = BigInteger.DivRem(p, x, out p);   // p -= [p/x] * x
                c -= quotient * a;
                d -= quotient * b;
                BigInteger auxiliary;
                auxiliary = x; x = p; p = auxiliary;
                auxiliary = a; a = c; c = auxiliary;
                auxiliary = b; b = d; d = auxiliary;
                result.Add(new Equation(a, b, x));
                if (result.Count > 3)
                {
                    result.RemoveAt(0);
                }
            }
            return result;
        }

        //find k1, k2 so k = k1 + k2 * lambda mod n; k1 and k2 have b bits and k, n and lambda have 2*b bits each.
        public void SplitKToHalfBits(BigInteger k, BigInteger lambda, out BigInteger k1, out BigInteger k2)
        {
            BigInteger n = Order.ToBigInteger();
            BigInteger quitBound = BigInteger.One << (N * 64 / 2);
            List<Equation> list = SplitNumberByExtendedEuclid(n, lambda, quitBound);
            BigInteger a1 = list[1].R;
            BigInteger b1 = -list[1].B;
            int listIndex = list[0].NormRBSquared <= list[2].NormRBSquared ? 0 : 2;
            BigInteger a2 = list[listIndex].R, b2 = -list[listIndex].B;

            //BigInteger test1 = (a1 + b1 * lambda) % n;    //should be == 0.
            //BigInteger test2 = (a2 + b2 * lambda) % n;    //should be == 0.

            BigInteger halfn = n >> 1;
            BigInteger c1 = b2 * k; c1 = c1 >= 0 ? (c1 + halfn) / n : (c1 - halfn) / n;
            BigInteger c2 = -b1 * k; c2 = c2 >= 0 ? (c2 + halfn) / n : (c2 - halfn) / n;
            k1 = k + c1 * a1 + c2 * a2;
            k2 = c1 * b1 + c2 * b2;

            BigInteger test = (k - (k1 + k2 * lambda)) % n;
            if (test != 0)
            {
                throw new InvalidOperationException("Cannot split k to half size k1 and k2.");
            }
        }

        public ECPoint ECMultiplication(BigInteger factor)
        {
            var result = MultiplyWithEndomorphism(G, factor.ToULong());
            return new ECPoint(result.X.ToBigInteger(), result.Y.ToBigInteger(), result.IsInfinity);
        }

        public FastECPoint ECMultiplication(FastInteger factor)
        {
            var result = MultiplyWithEndomorphism(G, factor.ToULong());
            return new FastECPoint(result.X.ToFastInteger(), result.Y.ToFastInteger(), result.IsInfinity);
        }

        //return: P * factor on the elliptic curves domain.
        public OwnECPoint MultiplyWithEndomorphism(OwnECPoint p, ulong[] factor)
        {
            OwnECPoint result;
            if (!this.SupportsEndomorphism)
            {
                result = this.Multiply(p, factor);
                return result;
            }

            BigInteger lambda = W1ModOrder.ToBigInteger();
            BigInteger k = factor.ToBigInteger(), k1, k2;
            SplitKToHalfBits(k, lambda, out k1, out k2);

            OwnECPoint q = p.Copy();
            MultiplyModP(q.X, W1ModP);

            //OwnECPoint test1 = Multiply(p, W1ModOrder);
            //bool testok = test1.ToString() == q.ToString();
            //if (!testok)
            //{
            //    Debugger.Break();
            //}

            result = DualMultiply(p, k1, q, k2);
            return result;
        }

        //return: P * factor on the elliptic curves domain.
        public OwnECPoint DualMultiply(OwnECPoint p1, BigInteger k1, OwnECPoint p2, BigInteger k2)
        {
            if (k1 < 0)
            {
                p1 = Negative(p1);
                k1 = -k1;
            }
            if (k2 < 0)
            {
                p2 = Negative(p2);
                k2 = -k2;
            }

            BigInteger exp1 = k1 ^ (k1 * 3);
            BigInteger exp2 = k2 ^ (k2 * 3);
            BigInteger max = BigInteger.Max(exp1, exp2);
            int high = max <= 0 ? 0 : 1 + (int)Math.Ceiling(BigInteger.Log(max, 2));
            JacobianECPoint x = new JacobianECPoint(p1.X, p1.Y, true);
            OwnECPoint cacheAdd = p1.Copy();
            OwnECPoint cacheSub = p1.Copy();
            Add(ref cacheAdd, p2);
            Subtract(ref cacheSub, p2);

            BigInteger bit = BigInteger.One << high;
            for (int i = high; --i >= 0; bit >>= 1)
            {
                GetDouble(ref x);
                if ((exp1 & bit) != 0 && (exp2 & bit) != 0)
                {
                    int branch = ((k1 & bit) != 0 ? 0 : 1) + ((k2 & bit) != 0 ? 0 : 2);
                    switch (branch)
                    {
                        case 0: Subtract(ref x, cacheAdd); break;
                        case 1: Add(ref x, cacheSub); break;
                        case 2: Subtract(ref x, cacheSub); break;
                        case 3: Add(ref x, cacheAdd); break;
                    }
                    continue;
                }
                if ((exp1 & bit) != 0)
                {
                    if ((k1 & bit) != 0)
                    {
                        Subtract(ref x, p1);
                    }
                    else
                    {
                        Add(ref x, p1);
                    }
                    continue;
                }
                if ((exp2 & bit) != 0)
                {
                    if ((k2 & bit) != 0)
                    {
                        Subtract(ref x, p2);
                    }
                    else
                    {
                        Add(ref x, p2);
                    }
                    continue;
                }
            }

            OwnECPoint result = x.ToECPoint(this);
            return result;
        }
    }
}
