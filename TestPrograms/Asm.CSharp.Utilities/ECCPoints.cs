using Bitsy.Core;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Numerics;
using System.Text;

namespace Utilities
{
    [DebuggerDisplay("{ToString()}")]
    public struct OwnECPoint
    {
        public ulong[] X { get; set; }
        public ulong[] Y { get; set; }
        public bool IsInfinity { get; set; }

        public OwnECPoint(ulong[] x, ulong[] y, bool isInfinity) : this()
        {
            this.X = new ulong[x.Length];
            this.Y = new ulong[y.Length];
            x.CopyTo(this.X, 0);
            y.CopyTo(this.Y, 0);
            this.IsInfinity = isInfinity;
        }

        const string hexa = "0123456789ABCDEF";

        public static string ToHexa(ulong[] number)
        {
            string hs = string.Join("", Enumerable.Range(0, number.Length * 16).Reverse().
                Select(idx => hexa[(int)(number[idx >> 4] >> ((idx & 15) << 2)) & 15]));
            return hs;
        }

        public override string ToString()
        {
            return this.IsInfinity ? "Infinity" : "X = " + ToHexa(this.X) + " ; Y = " + ToHexa(this.Y);
        }

        public JacobianECPoint ToJacobian
        {
            get
            {
                return new JacobianECPoint(this.X, this.Y, this.IsInfinity);
            }
        }

        public void SetInfinity()
        {
            this.IsInfinity = true;
        }

        public void SetFrom(ulong[] x, ulong[] y, bool isInfinity)
        {
            ECCBase.Let(this.X, x);
            ECCBase.Let(this.Y, y);
            this.IsInfinity = isInfinity;
        }

        public OwnECPoint Copy()
        {
            return new OwnECPoint(this.X, this.Y, this.IsInfinity);
        }

        public bool IsSameAs(OwnECPoint other)
        {
            if (this.IsInfinity != other.IsInfinity)
            {
                return false;
            }
            if (this.IsInfinity)
            {
                return true;
            }
            return this.X.IsEqual(other.X) && this.Y.IsEqual(other.Y);
        }

        public bool IsSameAs(FastECPoint other)
        {
            if (this.IsInfinity != other.IsInfinity)
            {
                return false;
            }
            if (this.IsInfinity)
            {
                return true;
            }
            other = other.Normalize();
            return new Bitsy.Core.FastInteger(this.X) == other.X && new Bitsy.Core.FastInteger(this.Y) == other.Y;
        }
    }

    [DebuggerDisplay("{ToString()}")]
    public struct JacobianECPoint
    {
        public ulong[] X { get; set; }
        public ulong[] Y { get; set; }
        public ulong[] Z { get; set; }
        public bool IsInfinity { get; set; }

        public JacobianECPoint(ulong[] x, ulong[] y, bool isInfinity) : this()
        {
            this.X = new ulong[x.Length];
            this.Y = new ulong[y.Length];
            this.Z = new ulong[x.Length];
            x.CopyTo(this.X, 0);
            y.CopyTo(this.Y, 0);
            this.Z[0] = 1;
            this.IsInfinity = isInfinity;
        }

        public override string ToString()
        {
            return this.IsInfinity ? "Infinity" : "X = " + OwnECPoint.ToHexa(this.X) + " ; Y = " + OwnECPoint.ToHexa(this.Y) + " ; Z = " + OwnECPoint.ToHexa(this.Z);
        }

        public OwnECPoint ToECPoint(ECCBase ecc)
        {
            ulong[] inverse = new ulong[ecc.N];
            inverse[0] = 1;
            ecc.RemodFriendlyPrime(this.Z);
            AsmX64Operations.DivideModuloPrime(inverse, this.Z, ecc.P, ecc.N);
            ulong[] inverseSquare = inverse.ToArray();
            ecc.SquareModP(inverseSquare);
            ecc.MultiplyModP(inverse, inverseSquare);
            ulong[] nx = this.X.ToArray(); ecc.MultiplyModP(nx, inverseSquare);
            ulong[] ny = this.Y.ToArray(); ecc.MultiplyModP(ny, inverse);
            ecc.RemodFriendlyPrime(nx);
            ecc.RemodFriendlyPrime(ny);
            return new OwnECPoint(nx, ny, this.IsInfinity);
        }

        public void SetInfinity()
        {
            this.IsInfinity = true;
        }

        public void SetFrom(ulong[] x, ulong[] y, ulong[] z, bool isInfinity)
        {
            ECCBase.Let(this.X, x);
            ECCBase.Let(this.Y, y);
            ECCBase.Let(this.Z, z);
            this.IsInfinity = false;
        }
    }

    public struct Equation
    {
        public BigInteger A { get; set; }
        public BigInteger B { get; set; }
        public BigInteger R { get; set; }

        public BigInteger NormRBSquared { get { return this.R * this.R + this.B * this.B; } }

        public Equation(BigInteger a, BigInteger b, BigInteger r) : this()
        {
            this.A = a;
            this.B = b;
            this.R = r;
        }
    }
}
