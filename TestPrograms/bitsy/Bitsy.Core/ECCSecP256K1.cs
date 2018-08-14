using System;
using System.Collections.Generic;
using System.Linq;
using System.Numerics;
using System.Text;

namespace Bitsy.Core
{
    public class ECCSecP256K1 : ECCBaseClass
    {
        public override int BitsCount => 256;

        public override int BytesCount => 32;

        public override BigInteger P => Secp256k1.P;

        public override BigInteger N => Secp256k1.N;

        public override ECPoint G => Secp256k1.G;

        public override ECPoint ECMultiplication(BigInteger factor)
        {
            ECPoint result = Secp256k1.G.Multiply(factor);
            return result;
        }
    }

    public class ECCSecP256K1FastInteger : ECCBaseClass
    {
        public override int BitsCount => 256;

        public override int BytesCount => 32;

        public override BigInteger P => Secp256k1.P;

        public override BigInteger N => Secp256k1.N;

        public override ECPoint G => Secp256k1.G;

        public override ECPoint ECMultiplication(BigInteger factor)
        {
            FastECPoint result = Secp256k1.FG.Multiply(new FastInteger(factor.ToByteArray()));
            return new ECPoint(new BigInteger(result.X.ToByteArray()), new BigInteger(result.Y.ToByteArray()));
        }
    }
}
