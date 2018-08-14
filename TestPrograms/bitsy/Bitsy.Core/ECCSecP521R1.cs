using System;
using System.Collections.Generic;
using System.Linq;
using System.Numerics;
using System.Text;

namespace Bitsy.Core
{
    public class ECCSecP521R1 : ECCBaseClass
    {
        private static readonly ECPoint _G = new ECPoint(ECC521Point.G.X, ECC521Point.G.Y, ECC521Point.G.IsInfinity);
        public override int BitsCount => 521;

        public override int BytesCount => 66;

        public override BigInteger P => ECC521Point.P;

        public override BigInteger N => ECC521Point.OrderN;

        public override ECPoint G => _G;

        public override ECPoint ECMultiplication(BigInteger factor)
        {
            ECC521Point point = ECC521Point.G.Multiply(factor);
            return new ECPoint(point.X, point.Y, point.IsInfinity);
        }
    }
}
