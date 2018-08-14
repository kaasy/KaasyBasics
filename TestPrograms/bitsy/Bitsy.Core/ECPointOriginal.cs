using System;
using System.Collections.Generic;
using System.Linq;
using System.Numerics;

namespace Bitsy.Core
{
    public class ECPoint : ICloneable
    {
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

        public static ECPoint Infinity
        {
            get
            {
                return new ECPoint(true);
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

        private ECPoint(bool infinity)
        {
            if (!infinity)
            {
                throw new ArgumentException("This constructor is only for creating the point Infinity");
            }

            _isInfinity = true;
        }

        public ECPoint(BigInteger x, BigInteger y)
        {
            _x = x;
            _y = y;
        }

        public ECPoint(BigInteger x, BigInteger y, BigInteger z, bool isInfinity) : this(x, y, isInfinity)
        {
            _z = z;
        }

        public byte[] EncodePoint(bool compressed)
        {
            if (IsInfinity)
            {
                return new byte[1];
            }

            byte[] x = X.ToByteArrayUnsigned(true);
            byte[] encoded;
            if (!compressed)
            {
                byte[] y = Y.ToByteArrayUnsigned(true);
                encoded = new byte[65];
                encoded[0] = 0x04;
                Buffer.BlockCopy(y, 0, encoded, 33 + (32 - y.Length), y.Length);
            }
            else
            {
                encoded = new byte[33];
                encoded[0] = (byte)(Y.TestBit(0) ? 0x03 : 0x02);
            }

            Buffer.BlockCopy(x, 0, encoded, 1 + (32 - x.Length), x.Length);
            return encoded;
        }

        public static ECPoint DecodePoint(byte[] encoded)
        {
            if (encoded == null || (encoded.Length != 33 && encoded.Length != 65))
            {
                throw new FormatException("Invalid encoded point");
            }

            if (encoded[0] == 0x04)
            {
                // uncompressed
                byte[] unsigned = new byte[32];

                Buffer.BlockCopy(encoded, 1, unsigned, 0, 32);

                var x = unsigned.ToBigIntegerUnsigned(true);

                Buffer.BlockCopy(encoded, 33, unsigned, 0, 32);
                var y = unsigned.ToBigIntegerUnsigned(true);

                return new ECPoint(x, y);
            }
            else if (encoded[0] == 0x02 || encoded[0] == 0x03)
            {
                // compressed
                byte[] unsigned = new byte[32];

                Buffer.BlockCopy(encoded, 1, unsigned, 0, 32);
                var x = unsigned.ToBigIntegerUnsigned(true);

                // solve y
                var y = ((x * x * x + 7) % Secp256k1.P).ShanksSqrt(Secp256k1.P);

                bool negate = false;
                if (y.TestBit(0))
                {
                    if (encoded[0] == 0x02)
                    {
                        negate = true;
                    }
                }
                else
                {
                    if (encoded[0] == 0x03)
                    {
                        negate = true;
                    }
                }

                if (negate)
                {
                    // negate
                    y = -y + Secp256k1.P;
                }

                return new ECPoint(x, y);
            }
            else
            {
                throw new FormatException("Invalid encoded point");
            }
        }

        public ECPoint Negate()
        {
            ECPoint r = (ECPoint)Clone();
            r._y = Secp256k1.P - r._y;
            return r;
        }

        public ECPoint Subtract(ECPoint b)
        {
            return Add(b.Negate());
        }

        private static readonly BigInteger LowMask = (BigInteger.One << 256) - BigInteger.One;

        public static BigInteger modP(BigInteger x)
        {
            while (true)
            {
                var high = x >> 256;
                if (high.IsZero)
                {
                    return x;
                }
                x = (x & LowMask) + high * 0x1000003D1;
            }
        }

        public static void remodP(ref BigInteger x)
        {
            while (x >= Secp256k1.P)
            {
                x -= Secp256k1.P;
            }
            while (x < 0)
            {
                x += Secp256k1.P;
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
        public ECPoint Add(ECPoint q)
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
                    return ECPoint.Infinity;
                }
                return this.Twice();
            }

            var H2 = modP(H.Square());
            var H3 = modP(H2 * H);
            var U1 = modP(this.X * H2);

            var X3 = modP(R.Square() - (H3 + (U1 << 1)));
            var Y3 = modP(R * (U1 - X3) - this.Y * H3);
            var Z3 = modP(H * this.Z);
            var result = new ECPoint(X3, Y3, Z3, false);
            return result;
        }

        //if (Y == 0)
        //  return POINT_AT_INFINITY
        //S = 4*X*Y^2
        //M = 3*X^2 + a*Z^4
        //X' = M^2 - 2*S
        //Y' = M*(S - X') - 8*Y^4
        //Z' = 2*Y*Z
        //return (X', Y', Z')
        public ECPoint Twice()
        {   // (3M + 4S)
            if (this.Y.IsZero || this.IsInfinity)
            {
                return ECPoint.Infinity;
            }
            var y2 = modP(this.Y.Square());
            var S = modP(this.X * (y2 << 2));
            var M = this.X.Square();
            M = modP(M + (M << 1));
            var x = modP(M.Square() - (S << 1));
            var y = modP(M * (S - x) - (y2.Square() << 3));
            var z = modP(this.Z * (this.Y << 1));
            var result = new ECPoint(x, y, z, false);
            return result;
        }

        public ECPoint Multiply(BigInteger b)
        {
            if (b.Sign == -1)
            {
                throw new FormatException("The multiplicator cannot be negative");
            }

            b = b % Secp256k1.N;
            BigInteger exp = (b * 3) ^ b;
            ECPoint result = ECPoint.Infinity;
            ECPoint affine = this.Normalize();
            ECPoint negative = affine.Negate();

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

        public ECPoint Normalize()
        {
            if (this.Z == BigInteger.One || this.IsInfinity)
            {
                return this;
            }
            var inverse = this.Z.ModInverse(Secp256k1.P);
            var squareInverse = modP(inverse.Square());
            var power3Inverse = modP(inverse * squareInverse);
            var x = modP(this.X * squareInverse);
            var y = modP(this.Y * power3Inverse);
            remodP(ref x);
            remodP(ref y);
            var result = new ECPoint(x, y, this.IsInfinity);
            return result;
        }

        public ECPoint(BigInteger x, BigInteger y, bool isInfinity)
        {
            _x = x;
            _y = y;
            _isInfinity = isInfinity;
        }

        public object Clone()
        {
            return new ECPoint(_x, _y, _z, _isInfinity);
        }

        public override string ToString()
        {
            var point = this.Normalize();
            return point.IsInfinity ? "Infinity" : "X = " + point.X.ToHex().ToUpperInvariant() + " ; Y = " + point.Y.ToHex().ToUpperInvariant();
        }

        public byte[] CompressedBytes
        {
            get
            {
                byte[] data = this.X.ToByteArrayUnsigned(false);
                Array.Resize(ref data, 256 / 8 + 1);    //since P = odd then the parity bit determines the correct Y coordinate.
                data[data.Length - 1] = (byte)(this.Y & 0xFF);
                return data;
            }
        }

        public byte[] UncompressedBytes
        {
            get
            {
                byte[] dataX = this.X.ToByteArrayUnsigned(false);
                byte[] dataY = this.Y.ToByteArrayUnsigned(false);
                Array.Resize(ref dataX, 256 / 8);
                Array.Resize(ref dataY, 256 / 8);
                return dataX.Concat(dataY).ToArray();
            }
        }

        public ECPoint(byte[] bytes, bool isCompressed)
        {
            this._x = bytes.Take(256 / 8).ToBigIntegerUnsigned(false);
            if (isCompressed)
            {
                var x2 = modP(this.X.Square());
                var x3 = modP(x2 * this.X);
                this._y = (x3 + 7).ShanksSqrt(Secp256k1.P);
                byte lowerY = (byte)(this._y & 255);
                if (lowerY != bytes[256 / 8])
                {
                    this._y = Secp256k1.P - this._y;
                    lowerY = (byte)(this._y & 255);
                    if (lowerY != bytes[256 / 8])
                    {
                        throw new InvalidOperationException("Compressed bytes verification byte failed.");
                    }
                }
            }
            else
            {
                this._y = bytes.Skip(256 / 8).Take(256 / 8).ToBigIntegerUnsigned(false);
            }
        }
    }
}
