using System;
using System.Collections.Generic;
using System.Linq;
using System.Numerics;
using System.Text;
using Bitsy.Core;

namespace Utilities
{
    //Y^2 = X^3 + a*X + b 'modulo' P
    public class SECP256K1Base : ECCBase
    {
        private void addModP(ulong[] input1, ulong[] input2, ulong[] result)
        {
            AsmX64Operations.AddModP(input1, input2, result, LowP, N);
        }
        private void subModP(ulong[] input1, ulong[] input2, ulong[] result)
        {
            AsmX64Operations.SubtractModP(input1, input2, result, LowP, N);
        }
        private void addScaledModP(ulong[] result, ulong[] input, ulong multiplier)
        {
            AsmX64Operations.AddScaledModP(result, input, multiplier, LowP, N);
        }
        private void subScaledModP(ulong[] result, ulong[] input, ulong multiplier)
        {
            AsmX64Operations.SubScaledModP(result, input, multiplier, LowP, N);
        }

        public SECP256K1Base() : base()
        {
            this.AdditionModP = this.addModP;
            this.SubtractModP = this.subModP;
            this.AddScaledModP = this.addScaledModP;
            this.SubScaledModP = this.subScaledModP;
        }

        private const int numberOfQwords = 256 / 64;
        public override int N { get { return numberOfQwords; } }
        public override int BitsCount { get { return 256; } }
        public override int BytesCount { get { return 32; } }
        private static int karatsubaBufferN = AsmX64Operations.GetKaratsubaMultiplicationBufferSize(numberOfQwords);
        public override int KaratsubaBufferN { get { return karatsubaBufferN; } }
        private const ulong LowP = 0x1000003D1UL;

        public static readonly ulong[] _P = new ulong[] { 0xFFFFFFFEFFFFFC2F, ulong.MaxValue, ulong.MaxValue, ulong.MaxValue };  //P = 2^256 - (2^32 + 2^9 + 2^8 + 2^7 + 2^6 + 2^4 + 1) = 2^256 - 0x1'000003D1 = 2^256 - 4294968273
        public static readonly ulong[] _A = new ulong[] { 0, 0, 0, 0 };
        public static readonly ulong[] _B = new ulong[] { 7, 0, 0, 0 };
        public static readonly ulong[] _One = new ulong[] { 1, 0, 0, 0 };
        public static readonly ulong[] _W1ModP = new ulong[] { 0x3ec693d68e6afa40, 0x630fb68aed0a766a, 0x919bb86153cbcb16, 0x851695d49a83f8ef }; //W2ModP ^ 3 == 1 mod P ; W2ModP ^ 2 == W1ModP mod P
        public static readonly ulong[] _W2ModP = new ulong[] { 0xc1396c28719501ee, 0x9cf0497512f58995, 0x6e64479eac3434e9, 0x7ae96a2b657c0710 }; //W1ModP ^ 3 == 1 mod P ; W1ModP ^ 2 == W2ModP mod P

        public static readonly ulong[] _W1ModOrder = new ulong[] { 0xe0cfc810b51283ce, 0xa880b9fc8ec739c2, 0x5ad9e3fd77ed9ba4, 0xac9c52b33fa3cf1f }; //W1ModOrder ^ 3 == 1 mod Order
        public static readonly ulong[] _W2ModOrder = new ulong[] { 0xdf02967c1b23bd72, 0x122e22ea20816678, 0xa5261c028812645a, 0x5363ad4cc05c30e0 }; //W2ModOrder ^ 3 == 1 mod Order

        public static readonly OwnECPoint _G = new OwnECPoint(
            new ulong[] { 0x59F2815B16F81798, 0x029BFCDB2DCE28D9, 0x55A06295CE870B07, 0x79BE667EF9DCBBAC },
            new ulong[] { 0x9C47D08FFB10D4B8, 0xFD17B448A6855419, 0x5DA4FBFC0E1108A8, 0x483ADA7726A3C465 }, false);
        public static readonly ulong[] _Order = new ulong[] { 0xBFD25E8CD0364141, 0xBAAEDCE6AF48A03B, ulong.MaxValue - 1, ulong.MaxValue };

        public override bool AIsZero { get { return true; } }
        public override bool AIsMinus3 { get { return false; } }

        public override void GetModFriendlyPrime(ulong[] inputOutput2x)
        {
            AsmX64Operations.GetModFriendlyPrime(inputOutput2x, LowP, N);
        }

        public override void RemodFriendlyPrime(ulong[] inputOutput1x)
        {
            AsmX64Operations.RemodFriendlyPrime(inputOutput1x, LowP, N);
        }

        public override bool SupportsEndomorphism { get { return true; } }

        public override ulong[] P { get { return _P; } }
        public override ulong[] A { get { return _A; } }
        public override ulong[] B { get { return _B; } }
        public override ulong[] One { get { return _One; } }

        public override ulong[] W1ModP { get { return _W1ModP; } }
        public override ulong[] W2ModP { get { return _W2ModP; } }
        public override ulong[] W1ModOrder { get { return _W1ModOrder; } }
        public override ulong[] W2ModOrder { get { return _W2ModOrder; } }

        public override OwnECPoint G { get { return _G; } }
        public override ulong[] Order { get { return _Order; } }
        public override int H { get { return 1; } }

        private static readonly JacobianECPoint _InfinityJacobian = new JacobianECPoint(_One, _One, true);
        private static readonly OwnECPoint _Infinity = new OwnECPoint(_One, _One, true);

        public override JacobianECPoint InfinityJacobian { get { return _InfinityJacobian; } }
        public override OwnECPoint Infinity { get { return _Infinity; } }
    }

    public class SECP256K1 : ECCBaseClass
    {
        private static readonly BigInteger PrimeModulo = SECP256K1Base._P.ToBigInteger();
        private static readonly BigInteger OrderN = SECP256K1Base._Order.ToBigInteger();
        private static readonly ECPoint Generator = new ECPoint(SECP256K1Base._G.X.ToBigInteger(), SECP256K1Base._G.Y.ToBigInteger(), SECP256K1Base._G.IsInfinity);

        public override int BitsCount => 256;

        public override int BytesCount => 32;

        public override BigInteger P => PrimeModulo;

        public override BigInteger N => OrderN;

        public override ECPoint G => Generator;

        private SECP256K1Base ecc;

        public SECP256K1() : base()
        {
            this.ecc = new SECP256K1Base();
        }

        public override ECPoint ECMultiplication(BigInteger factor)
        {
            var result = this.ecc.MultiplyWithEndomorphism(this.ecc.G, factor.ToULong());
            return new ECPoint(result.X.ToBigInteger(), result.Y.ToBigInteger(), result.IsInfinity);
        }
    }
}
