using Bitsy.Core;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Numerics;
using System.Text;

namespace Bitsy.Core
{
    public class ECC521Point : ICloneable
    {
        private const int bits = 521;
        public static readonly BigInteger P = (BigInteger.One << bits) - BigInteger.One;

        private static readonly BigInteger B = new ulong[] {
            0xef451fd46b503f00UL, 0x3573df883d2c34f1UL, 0x1652c0bd3bb1bf07UL, 0x56193951ec7e937bUL,
            0xb8b489918ef109e1UL, 0xa2da725b99b315f3UL, 0x929a21a0b68540eeUL, 0x953eb9618e1c9a1fUL, 0x0051UL }.ToBigInteger();
        public static readonly ECC521Point G = new ECC521Point(new ulong[] {
            0xf97e7e31c2e5bd66UL, 0x3348b3c1856a429bUL, 0xfe1dc127a2ffa8deUL, 0xa14b5e77efe75928UL,
            0xf828af606b4d3dbaUL, 0x9c648139053fb521UL, 0x9e3ecb662395b442UL, 0x858e06b70404e9cdUL, 0x00c6UL }.ToBigInteger(), new ulong[] {
            0x88be94769fd16650UL, 0x353c7086a272c240UL, 0xc550b9013fad0761UL, 0x97ee72995ef42640UL,
            0x17afbd17273e662cUL, 0x98f54449579b4468UL, 0x5c8a5fb42c7d1bd9UL, 0x39296a789a3bc004UL, 0x0118UL }.ToBigInteger(), false);
        public static readonly BigInteger OrderN = new ulong[] {
            0xBB6FB71E91386409UL, 0x3BB5C9B8899C47AEUL, 0x7FCC0148F709A5D0UL, 0x51868783BF2F966BUL,
            0xFFFFFFFFFFFFFFFAUL, 0xFFFFFFFFFFFFFFFFUL, 0xFFFFFFFFFFFFFFFFUL, 0xFFFFFFFFFFFFFFFFUL, 0x01FFUL }.ToBigInteger();

        private BigInteger _x;
        private BigInteger _y;
        private BigInteger _z = BigInteger.One;

        public BigInteger X
        {
            get
            {
                return _x;
            }
        }

        public BigInteger Y
        {
            get
            {
                return _y;
            }
        }

        public BigInteger Z
        {
            get
            {
                return _z;
            }
        }

        public static ECC521Point Infinity
        {
            get
            {
                return new ECC521Point(true);
            }
        }

        private bool _isInfinity = false;
        public bool IsInfinity
        {
            get
            {
                return _isInfinity;
            }
        }

        private ECC521Point(bool infinity)
        {
            if (!infinity)
            {
                throw new ArgumentException("This constructor is only for creating the point Infinity");
            }

            _isInfinity = true;
        }

        public ECC521Point(BigInteger x, BigInteger y)
        {
            _x = x;
            _y = y;
        }

        public ECC521Point(BigInteger x, BigInteger y, BigInteger z, bool isInfinity) : this(x, y, isInfinity)
        {
            _z = z;
        }

        public ECC521Point Negate()
        {
            ECC521Point r = (ECC521Point)Clone();
            r._y = P - r._y;
            return r;
        }

        public ECC521Point Subtract(ECC521Point b)
        {
            return Add(b.Negate());
        }

        private static BigInteger modP(BigInteger x)
        {
            while (true)
            {
                var high = x >> bits;
                if (high.IsZero)
                {
                    return x;
                }
                x = (x & P) + high;
            }
        }

        private static void remodP(ref BigInteger x)
        {
            while (x >= P)
            {
                x -= P;
            }
            while (x < 0)
            {
                x += P;
            }
        }

        //U1 = X1*Z2^2
        //U2 = X2*Z1^2
        //S1 = Y1*Z2^3
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
        public ECC521Point Add(ECC521Point q)
        {   // (8M + 3S)
            if (this.IsInfinity)
            {
                return q;
            }
            if (q.IsInfinity)
            {
                return this;
            }

            var z2 = modP(this.Z.Square());
            var z3 = modP(z2 * this.Z);
            var U2 = modP(q.X * z2);
            var S2 = modP(q.Y * z3);
            var H = U2 - this.X;
            var R = S2 - this.Y;
            if (H.IsZero)
            {
                if (!R.IsZero)
                {
                    return ECC521Point.Infinity;
                }
                return this.Twice();
            }

            var H2 = modP(H.Square());
            var H3 = modP(H2 * H);
            var U1 = modP(this.X * H2);

            var X3 = modP(R.Square() - (H3 + (U1 << 1)));
            var Y3 = modP(R * (U1 - X3) - this.Y * H3);
            var Z3 = modP(H * this.Z);
            var result = new ECC521Point(X3, Y3, Z3, false);
            return result;
        }

        //if (Y == 0)
        //  return POINT_AT_INFINITY
        //S = 4*X*Y^2
        //M = 3*X^2 + a*Z^4 = 3(X - Z^2) * (X + Z^2)
        //X' = M^2 - 2*S
        //Y' = M*(S - X') - 8*Y^4
        //Z' = 2*Y*Z
        //return (X', Y', Z')
        public ECC521Point Twice()
        {   // (4M + 4S)
            if (this.Y.IsZero || this.IsInfinity)
            {
                return ECC521Point.Infinity;
            }
            var y2 = modP(this.Y.Square());
            var z2 = modP(this.Z.Square());
            var S = modP(this.X * y2 << 2);
            var M = (this.X - z2) * (this.X + z2);
            M = modP(M + (M << 1));
            var x = modP(M.Square() - (S << 1));
            var y = modP(M * (S - x) - (y2.Square() << 3));
            var z = modP(this.Z * (this.Y << 1));
            var result = new ECC521Point(x, y, z, false);
            return result;
        }

        public ECC521Point Multiply(BigInteger b)
        {
            if (b.Sign == -1)
            {
                throw new FormatException("The multiplicator cannot be negative");
            }
            if (b >= OrderN)
            {
                b %= OrderN;
            }
            BigInteger exp = (b * 3) ^ b;
            ECC521Point result = ECC521Point.Infinity;
            ECC521Point affine = this.Normalize();
            ECC521Point negative = affine.Negate();

            int high = exp <= 0 ? 0 : (int)Math.Floor(BigInteger.Log(exp, 2));
            BigInteger bit = BigInteger.One << high;
            for (int i = high; --i >= 0; bit >>= 1)
            {
                result = result.Twice();
                if (!(exp & bit).IsZero)
                {
                    result = result.Add(!(b & bit).IsZero ? negative : affine);
                }
            }

            result = result.Normalize();
            return result;
        }

        public ECC521Point Normalize()
        {
            if (this.Z == BigInteger.One || this.IsInfinity)
            {
                return this;
            }
            var inverse = this.Z.ModInverse(P);
            var squareInverse = modP(inverse.Square());
            var power3Inverse = modP(inverse * squareInverse);
            var x = modP(this.X * squareInverse);
            var y = modP(this.Y * power3Inverse);
            remodP(ref x);
            remodP(ref y);
            var result = new ECC521Point(x, y, this.IsInfinity);
            return result;
        }

        public ECC521Point(BigInteger x, BigInteger y, bool isInfinity)
        {
            _x = x;
            _y = y;
            _isInfinity = isInfinity;
        }

        public object Clone()
        {
            return new ECC521Point(_x, _y, _z, _isInfinity);
        }

        public override string ToString()
        {
            return this.IsInfinity ? "Infinity" : "X = " + this.X.ToHex().ToUpperInvariant() + " ; Y = " + this.Y.ToHex().ToUpperInvariant() + " ; Z = " + this.Z.ToHex().ToUpperInvariant();
        }

        public byte[] CompressedBytes
        {
            get
            {
                byte[] data = this.X.ToByteArray();
                Array.Resize(ref data, (bits + 1 + 7) / 8);    //since P = odd then the parity bit determines the correct Y coordinate.
                data[bits >> 3] |= (byte)((int)(this.Y & 0xFF) << (bits & 7));
                return data;
            }
        }

        public byte[] UncompressedBytes
        {
            get
            {
                byte[] dataX = this.X.ToByteArray();
                byte[] dataY = this.Y.ToByteArray();
                Array.Resize(ref dataX, (bits + 7) / 8);
                Array.Resize(ref dataY, (bits + 7) / 8);
                return dataX.Concat(dataY).ToArray();
            }
        }

        public ECC521Point(byte[] bytes, bool isCompressed)
        {
            byte[] bx = bytes.Take((bits + 7) / 8).ToArray();
            if ((bits & 7) != 0)
            {
                bx[bx.Length - 1] &= (byte)((1 << (bits & 7)) - 1);
            }
            this._x = bx.ToBigIntegerUnsigned(false);
            if (isCompressed)
            {
                int checkPosition = bits >> 3;
                byte incoming = bytes[checkPosition];
                var x2 = modP(this.X.Square());
                this._y = modP(x2 * this.X - 3 * this.X + B).ShanksSqrt(P);
                byte lowerY = (byte)(this._y & 0xFF);
                byte mask = 0xFF;
                lowerY <<= bits & 7;
                mask <<= bits & 7;
                if (lowerY != (incoming & mask))
                {
                    this._y = P - this._y;
                    lowerY = (byte)(this._y & 0xFF);
                    lowerY <<= bits & 7;
                    if (lowerY != (incoming & mask))
                    {
                        throw new InvalidOperationException("Compressed bytes verification byte failed.");
                    }
                }
            }
            else
            {
                this._y = bytes.Skip(bx.Length).Take(bx.Length).ToBigIntegerUnsigned(false);
            }
        }
    }
}
