using System;
using System.Collections.Generic;
using System.Linq;
using System.Numerics;
using System.Text;

namespace Bitsy.Core
{
    public abstract class ECCBaseClass
    {
        public abstract int BitsCount { get; }
        public abstract int BytesCount { get; }
        public int CompressedBytesCount { get { return (this.BitsCount + 1 + 7) >> 3; } }
        public abstract BigInteger P { get; }
        public abstract BigInteger N { get; }
        public abstract ECPoint G { get; }

        public FastInteger GetRandomFactorModN_FastInteger()
        {
            BigInteger random = this.GetRandomFactorModN();
            return new FastInteger(random.ToByteArray());
        }

        public BigInteger GetRandomFactorModN()
        {
            byte[] bytes = new byte[(this.BitsCount + 7) / 8 + 1];
            using (var rng = System.Security.Cryptography.RandomNumberGenerator.Create())
            {
                while (true)
                {
                    rng.GetBytes(bytes);
                    bytes[bytes.Length - 1] = 0;
                    if ((this.BitsCount & 7) != 0)
                    {
                        bytes[bytes.Length - 2] &= (byte)((1 << (this.BitsCount & 7)) - 1);
                    }
                    var number = new BigInteger(bytes);
                    if (number < this.P)
                    {
                        return number;
                    }
                }
            }
        }

        public abstract ECPoint ECMultiplication(BigInteger factor);
    }
}
