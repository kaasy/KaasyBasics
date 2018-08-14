using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Numerics;
using Bitsy.Core;

namespace Utilities
{
    public class ECC
    {
        public BigInteger OrderN { get; private set; }

        public ECC()
        {
            this.OrderN = Secp256k1.N;
        }

        public ECPoint ECMultiplication(BigInteger factor)
        {
            ECPoint result = Secp256k1.G.Multiply(factor);
            return result;
        }

        public bool Verify(ECPoint point)
        {
            if (point.IsInfinity)
            {
                return true;
            }
            var x2 = ECPoint.modP(point.X.Square());
            var x3 = ECPoint.modP(x2 * point.X);
            var y2 = ECPoint.modP(point.Y.Square());
            var difference = x3 + 7 - y2;
            ECPoint.remodP(ref difference);
            return difference.IsZero;
        }
    }
}
