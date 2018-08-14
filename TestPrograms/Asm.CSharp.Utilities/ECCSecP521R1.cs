using System;
using System.Collections.Generic;
using System.Linq;
using System.Numerics;
using System.Text;
using Bitsy.Core;

namespace Utilities
{
    //Y^2 = X^3 + a*X + b 'modulo' P
    public class SECP521R1Base : ECCBase
    {
        private void addModP(ulong[] input1, ulong[] input2, ulong[] result)
        {
            AsmX64Operations.AddModP521(input1, input2, result, N);
        }
        private void subModP(ulong[] input1, ulong[] input2, ulong[] result)
        {
            AsmX64Operations.SubtractModP521(input1, input2, result, N);
        }
        private void addScaledModP(ulong[] result, ulong[] input, ulong multiplier)
        {
            AsmX64Operations.AddScaledModP521(result, input, multiplier, N);
        }
        private void subScaledModP(ulong[] result, ulong[] input, ulong multiplier)
        {
            AsmX64Operations.SubScaledModP521(result, input, multiplier, N);
        }

        public SECP521R1Base() : base()
        {
            this.AdditionModP = this.addModP;
            this.SubtractModP = this.subModP;
            this.AddScaledModP = this.addScaledModP;
            this.SubScaledModP = this.subScaledModP;
        }

        private const int numberOfQwords = (521 + 63) / 64;
        public override int N { get { return numberOfQwords; } }
        public override int BitsCount { get { return 521; } }
        public override int BytesCount { get { return 66; } }
        private static int karatsubaBufferN = AsmX64Operations.GetKaratsubaMultiplicationBufferSize(numberOfQwords);
        public override int KaratsubaBufferN { get { return karatsubaBufferN; } }

        public static readonly ulong[] _P = new ulong[] { ulong.MaxValue, ulong.MaxValue, ulong.MaxValue, ulong.MaxValue, ulong.MaxValue, ulong.MaxValue, ulong.MaxValue, ulong.MaxValue, 0x01FFUL };
        //P = 2^521 - 1
        public static readonly ulong[] _A = new ulong[] { ulong.MaxValue - 3, ulong.MaxValue, ulong.MaxValue, ulong.MaxValue, ulong.MaxValue, ulong.MaxValue, ulong.MaxValue, ulong.MaxValue, 0x01FFUL };
        public static readonly ulong[] _B = new ulong[] {
            0xef451fd46b503f00UL, 0x3573df883d2c34f1UL, 0x1652c0bd3bb1bf07UL, 0x56193951ec7e937bUL,
            0xb8b489918ef109e1UL, 0xa2da725b99b315f3UL, 0x929a21a0b68540eeUL, 0x953eb9618e1c9a1fUL, 0x0051UL };
        public static readonly ulong[] _One = new ulong[] { 1, 0, 0, 0, 0, 0, 0, 0, 0 };

        public static readonly OwnECPoint _G = new OwnECPoint(new ulong[] {
            0xf97e7e31c2e5bd66UL, 0x3348b3c1856a429bUL, 0xfe1dc127a2ffa8deUL, 0xa14b5e77efe75928UL,
            0xf828af606b4d3dbaUL, 0x9c648139053fb521UL, 0x9e3ecb662395b442UL, 0x858e06b70404e9cdUL, 0x00c6UL }, new ulong[] {
            0x88be94769fd16650UL, 0x353c7086a272c240UL, 0xc550b9013fad0761UL, 0x97ee72995ef42640UL,
            0x17afbd17273e662cUL, 0x98f54449579b4468UL, 0x5c8a5fb42c7d1bd9UL, 0x39296a789a3bc004UL, 0x0118UL }, false);
        public static readonly ulong[] _Order = new ulong[] {
            0xBB6FB71E91386409UL, 0x3BB5C9B8899C47AEUL, 0x7FCC0148F709A5D0UL, 0x51868783BF2F966BUL,
            0xFFFFFFFFFFFFFFFAUL, 0xFFFFFFFFFFFFFFFFUL, 0xFFFFFFFFFFFFFFFFUL, 0xFFFFFFFFFFFFFFFFUL, 0x01FFUL };

        public override bool AIsZero { get { return false; } }
        public override bool AIsMinus3 { get { return true; } }

        public override void GetModFriendlyPrime(ulong[] inputOutput2x)
        {
            AsmX64Operations.GetModFriendlyPrime521(inputOutput2x, N);
        }

        public override void RemodFriendlyPrime(ulong[] inputOutput1x)
        {
            AsmX64Operations.RemodFriendlyPrime521(inputOutput1x, N);
        }

        public override bool SupportsEndomorphism { get { return false; } }

        public override ulong[] P { get { return _P; } }
        public override ulong[] A { get { return _A; } }
        public override ulong[] B { get { return _B; } }
        public override ulong[] One { get { return _One; } }

        public override ulong[] W1ModP { get { throw new NotSupportedException(); } }
        public override ulong[] W2ModP { get { throw new NotSupportedException(); } }
        public override ulong[] W1ModOrder { get { throw new NotSupportedException(); } }
        public override ulong[] W2ModOrder { get { throw new NotSupportedException(); } }

        public override OwnECPoint G { get { return _G; } }
        public override ulong[] Order { get { return _Order; } }
        public override int H { get { return 1; } }

        private static readonly JacobianECPoint _InfinityJacobian = new JacobianECPoint(_One, _One, true);
        private static readonly OwnECPoint _Infinity = new OwnECPoint(_One, _One, true);

        public override JacobianECPoint InfinityJacobian { get { return _InfinityJacobian; } }
        public override OwnECPoint Infinity { get { return _Infinity; } }

        private static ulong[] FromHex(string hex)
        {
            hex = hex.Replace(" ", "");
            int bits = hex.Length * 4;
            int ulongs = (bits + 63) / 64;
            ulong[] result = new ulong[ulongs];
            for (int i = hex.Length; --i >= 0;)
            {
                char code = char.ToUpperInvariant(hex[hex.Length - i - 1]);
                ulong value = code >= 'A' ? (ulong)(code - 'A' + 10) : (ulong)(code - '0');
                result[i >> 4] |= value << ((i << 2) & 63);
            }
            return result;
        }

        public bool UnitTest()
        {
            OwnECPoint S = new OwnECPoint(
                FromHex(@"000001d5 c693f66c 08ed03ad 0f031f93 7443458f 601fd098 d3d0227b 4bf62873 af50740b 0bb84aa1 57fc847b cf8dc16a 8b2b8bfd 8e2d0a7d 39af04b0 89930ef6 dad5c1b4"),
                FromHex(@"00000144 b7770963 c63a3924 8865ff36 b074151e ac33549b 224af5c8 664c5401 2b818ed0 37b2b7c1 a63ac89e baa11e07 db89fcee 5b556e49 764ee3fa 66ea7ae6 1ac01823"), false);
            OwnECPoint T = new OwnECPoint(
                FromHex(@"000000f4 11f2ac2e b971a267 b80297ba 67c322db a4bb21ce c8b70073 bf88fc1c a5fde3ba 09e5df6d 39acb2c0 762c03d7 bc224a3e 197feaf7 60d63240 06fe3be9 a548c7d5"),
                FromHex(@"000001fd f842769c 707c93c6 30df6d02 eff399a0 6f1b36fb 9684f0b3 73ed0648 89629abb 92b1ae32 8fdb4553 42683849 43f0e922 2afe0325 9b32274d 35d1b958 4c65e305"), false);
            OwnECPoint ST = new OwnECPoint(
                FromHex(@"00000126 4ae115ba 9cbc2ee5 6e6f0059 e24b52c8 04632160 2c59a339 cfb757c8 9a59c358 a9a8e1f8 6d384b3f 3b255ea3 f73670c6 dc9f45d4 6b6a196d c37bbe0f 6b2dd9e9"),
                FromHex(@"00000062 a9c72b8f 9f88a271 690bfa01 7a6466c3 1b9cadc2 fc544744 aeb81707 2349cfdd c5ad0e81 b03f1897 bd9c8c6e fbdf6823 7dc3bb00 445979fb 373b20c9 a967ac55"), false);
            OwnECPoint STSub = new OwnECPoint(
                FromHex(@"00000129 2cb58b17 95ba4770 63fef7cd 22e42c20 f57ae94c eaad86e0 d21ff229 18b0dd3b 076d63be 253de24b c20c6da2 90fa54d8 3771a225 deecf914 9f79a8e6 14c3c4cd"),
                FromHex(@"00000169 5e3821e7 2c7cacaa dcf62909 cd83463a 21c6d033 93c527c6 43b36239 c46af117 ab7c7ad1 9a4c8cf0 ae95ed51 72988546 1aa2ce27 00a6365b ca3733d2 920b2267"), false);
            OwnECPoint S2X = new OwnECPoint(
                FromHex(@"00000128 79442f24 50c119e7 119a5f73 8be1f1eb a9e9d7c6 cf41b325 d9ce6d64 3106e9d6 1124a91a 96bcf201 305a9dee 55fa7913 6dc70083 1e54c3ca 4ff2646b d3c36bc6"),
                FromHex(@"00000198 64a8b885 5c2479cb efe375ae 553e2393 271ed36f adfc4494 fc0583f6 bd035988 96f39854 abeae5f9 a6515a02 1e2c0eef 139e71de 610143f5 3382f410 4dccb543"), false);
            JacobianECPoint ST1 = S.ToJacobian;
            JacobianECPoint ST2 = S.ToJacobian;
            JacobianECPoint S2 = S.ToJacobian;
            Add(ref ST1, T, false);
            Add(ref ST2, T, true);
            GetDouble(ref S2);
            OwnECPoint ST1a = ST1.ToECPoint(this);
            OwnECPoint ST2a = ST2.ToECPoint(this);
            OwnECPoint S2a = S2.ToECPoint(this);
            if (!IsEqual(ST1a.X, ST.X) || !IsEqual(ST1a.Y, ST.Y))
            {
                return false;
            }
            if (!IsEqual(ST2a.X, STSub.X) || !IsEqual(ST2a.Y, STSub.Y))
            {
                return false;
            }
            if (!IsEqual(S2a.X, S2X.X) || !IsEqual(S2a.Y, S2X.Y))
            {
                return false;
            }
            ulong[] factorD = FromHex(@"000001eb 7f81785c 9629f136 a7e8f8c6 74957109 73555411 1a2a866f a5a16669 9419bfa9 936c78b6 2653964d f0d6da94 0a695c72 94d41b2d 6600de6d fcf0edcf c89fdcb1");
            OwnECPoint dS = new OwnECPoint(
                FromHex(@"00000091 b15d09d0 ca0353f8 f96b93cd b13497b0 a4bb582a e9ebefa3 5eee61bf 7b7d041b 8ec34c6c 00c0c067 1c4ae063 318fb75b e87af4fe 859608c9 5f0ab477 4f8c95bb"),
                FromHex(@"00000130 f8f8b5e1 abb4dd94 f6baaf65 4a2d5810 411e77b7 423965e0 c7fd79ec 1ae563c2 07bd255e e9828eb7 a03fed56 5240d2cc 80ddd2ce cbb2eb50 f0951f75 ad87977f"), false);
            OwnECPoint rsd = Multiply(S, factorD);
            if (!IsEqual(rsd.X, dS.X) || !IsEqual(rsd.Y, dS.Y))
            {
                return false;
            }
            ulong[] factorE = FromHex(@"00000137 e6b73d38 f153c3a7 57561581 2608f2ba b3229c92 e21c0d1c 83cfad92 61dbb17b b77a6368 2000031b 9122c2f0 cdab2af7 2314be95 254de429 1a8f85f7 c70412e3");
            OwnECPoint dSeT = new OwnECPoint(
                FromHex(@"0000009d 3802642b 3bea152b eb9e05fb a247790f 7fc16807 2d363340 133402f2 585588dc 1385d40e bcb8552f 8db02b23 d687cae4 6185b275 28adb1bf 9729716e 4eba653d"),
                FromHex(@"0000000f e44344e7 9da6f49d 87c10637 44e5957d 9ac0a505 bafa8281 c9ce9ff2 5ad53f8d a084a2de b0923e46 501de579 7850c61b 229023dd 9cf7fc7f 04cd35eb b026d89d"), false);
            OwnECPoint ret = Multiply(T, factorE);
            Add(ref ret, rsd, false);
            if (!IsEqual(ret.X, dSeT.X) || !IsEqual(ret.Y, dSeT.Y))
            {
                return false;
            }
            return true;
        }
    }

    public class SECP521R1 : ECCBaseClass
    {
        private static readonly BigInteger PrimeModulo = SECP521R1Base._P.ToBigInteger();
        private static readonly BigInteger OrderN = SECP521R1Base._Order.ToBigInteger();
        private static readonly ECPoint Generator = new ECPoint(SECP521R1Base._G.X.ToBigInteger(), SECP521R1Base._G.Y.ToBigInteger(), SECP521R1Base._G.IsInfinity);

        public override int BitsCount => 521;

        public override int BytesCount => 66;

        public override BigInteger P => PrimeModulo;

        public override BigInteger N => OrderN;

        public override ECPoint G => Generator;

        private SECP521R1Base ecc;

        public SECP521R1() : base()
        {
            this.ecc = new SECP521R1Base();
        }

        public override ECPoint ECMultiplication(BigInteger factor)
        {
            var result = this.ecc.MultiplyWithEndomorphism(this.ecc.G, factor.ToULong());
            return new ECPoint(result.X.ToBigInteger(), result.Y.ToBigInteger(), result.IsInfinity);
        }
    }
}
