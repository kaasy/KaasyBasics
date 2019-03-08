using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using Utilities;

namespace X64Operations.UnitTest
{
    [TestClass]
    public class MultiplicationClass
    {
        private static ulong getUlong(Random random)
        {
            ulong i1 = (ulong)random.Next(1 << 22);
            ulong i2 = (ulong)random.Next(1 << 22);
            ulong i3 = (ulong)random.Next(1 << 20);
            return i1 | (i2 << 22) | (i3 << 44);
        }

        [TestMethod]
        public void LongMultiplicationNumberTest()
        {
            int n = 23 * 1024 - 731;
            ulong[] input1 = new ulong[n];
            ulong[] input2 = new ulong[n];
            ulong[] result1 = new ulong[n * 2];
            ulong[] result2 = new ulong[n * 2];
            ulong[] result3 = new ulong[n * 2];
            Random random = new Random(1001);
            for (int i = n; --i >= 0;)
            {
                input1[i] = getUlong(random);
                input2[i] = getUlong(random);
            }
            AsmX64Operations.SetKaratsubaBreakPoint(16);
            AsmX64Operations.SetFourierBreakPoint(12 * 1024);
            AsmX64Operations.FastestMultiplication(input1, input2, result1, n, false);
            AsmX64Operations.Karatsuba(input1, input2, result2, n);
            AsmX64Operations.Multiply(input1, input2, result3, n);

            for (int i = n * 2; --i >= 0;)
            {
                if (result1[i] != result2[i] || result2[i] != result3[i])
                {
                    Assert.IsTrue(false);
                    return;
                }
            }
            Assert.IsTrue(true);
        }

        [TestMethod]
        public void ProgressiveMultiplicationNumberTest()
        {
            AsmX64Operations.SetKaratsubaBreakPoint(16);
            AsmX64Operations.SetFourierBreakPoint(512);
            for (int n = 1; n < 2048; n++)
            {
                ulong[] input1 = new ulong[n];
                ulong[] input2 = new ulong[n];
                ulong[] result1 = new ulong[n * 2];
                ulong[] result2 = new ulong[n * 2];
                ulong[] result3 = new ulong[n * 2];
                Random random = new Random(1001);
                for (int i = n; --i >= 0;)
                {
                    input1[i] = getUlong(random);
                    input2[i] = getUlong(random);
                }
                AsmX64Operations.FastestMultiplication(input1, input2, result1, n, false);
                AsmX64Operations.Karatsuba(input1, input2, result2, n);
                AsmX64Operations.Multiply(input1, input2, result3, n);

                for (int i = n * 2; --i >= 0;)
                {
                    if (result1[i] != result2[i] || result2[i] != result3[i])
                    {
                        Assert.IsTrue(false);
                        return;
                    }
                }
            }
            Assert.IsTrue(true);
        }

        [TestMethod]
        public void LongSquareNumberTest()
        {
            int n = 23 * 1024 - 731;
            ulong[] input1 = new ulong[n];
            ulong[] result1 = new ulong[n * 2];
            ulong[] result2 = new ulong[n * 2];
            ulong[] result3 = new ulong[n * 2];
            Random random = new Random(1001);
            for (int i = n; --i >= 0;)
            {
                input1[i] = getUlong(random);
            }
            AsmX64Operations.SetKaratsubaBreakPoint(16);
            AsmX64Operations.SetFourierBreakPoint(12 * 1024);
            AsmX64Operations.FastestSquare(input1, result1, n, false);
            AsmX64Operations.KaratsubaSquare(input1, result2, n);
            AsmX64Operations.Square(input1, result3, n);

            for (int i = n * 2; --i >= 0;)
            {
                if (result1[i] != result2[i] || result2[i] != result3[i])
                {
                    Assert.IsTrue(false);
                    return;
                }
            }
            Assert.IsTrue(true);
        }

        [TestMethod]
        public void ProgressiveSquaringNumberTest()
        {
            AsmX64Operations.SetKaratsubaBreakPoint(16);
            AsmX64Operations.SetFourierBreakPoint(512);
            for (int n = 1; n < 2048; n++)
            {
                ulong[] input1 = new ulong[n];
                ulong[] result1 = new ulong[n * 2];
                ulong[] result2 = new ulong[n * 2];
                ulong[] result3 = new ulong[n * 2];
                Random random = new Random(1001);
                for (int i = n; --i >= 0;)
                {
                    input1[i] = getUlong(random);
                }
                AsmX64Operations.FastestSquare(input1, result1, n, false);
                AsmX64Operations.KaratsubaSquare(input1, result2, n);
                AsmX64Operations.Square(input1, result3, n);

                for (int i = n * 2; --i >= 0;)
                {
                    if (result1[i] != result2[i] || result2[i] != result3[i])
                    {
                        Assert.IsTrue(false);
                        return;
                    }
                }
            }
            Assert.IsTrue(true);
        }
    }
}
