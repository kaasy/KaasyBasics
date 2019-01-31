using Utilities;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Windows.Forms;

namespace Demo.App.Windows.Forms
{
    static class Program
    {
        /// <summary>
        /// The main entry point for the application.
        /// </summary>
        [STAThread]
        static void Main()
        {
            //IntegerNumber zzz = new IntegerNumber(0L);
            //string ssttrr = zzz.ToString();
            bool ok = true;
            //testMultiplicationSpeed();

            RealNumber ten = new RealNumber(10, 1000);
            RealNumber logTen = ten.GetLog();
            RealNumber tenExp1 = logTen.GetExp();
            RealNumber tenExp2 = logTen.GetTaylorExp();
            IntegerNumber i1 = IntegerNumber.One << 100000;
            IntegerNumber ia = i1.Sqrt();
            IntegerNumber ib = i1.InverseSqrt();

            RealNumber ssqrt10 = ten.GetSlowerSqrt();
            RealNumber sqrt10 = ten.GetSqrt();
            RealNumber invSqrt10 = ten.GetInverseSqrt();
            RealNumber one1 = sqrt10 * invSqrt10;
            ok &= (one1 - 1).IsZero;

            IntegerNumber int10 = IntegerNumber.Pow(10, 19 * 5000);
            IntegerNumber t1 = int10.Inverse();
            IntegerNumber t2 = (IntegerNumber.One << (int10.Digits * 128)) / int10;
            ok &= t1 == t2;
            IntegerNumber t3 = int10.InverseSqrt();
            IntegerNumber t4 = ((IntegerNumber.One << (int10.Digits * 96 * 2)) / int10).Sqrt();
            ok &= t3 == t4;

            ok &= Fourier235RealUnitTest.FFTUnitTest();
            ok &= Fourier235RealUnitTest.DCTUnitTest();
            ok &= Fourier235DoubleUnitTest.FFTUnitTest();
            ok &= Fourier235DoubleUnitTest.DCTUnitTest();
            ok &= Fourier235UnitTest.FFTUnitTest();
            ok &= Fourier235UnitTest.DCTUnitTest();

            ok &= RealNumbersUnitTest.ComplexArcTanSinCosUnitTest();
            ok &= RealNumbersUnitTest.NumericIntegerOperationsUnitTest();

            ok &= RealNumbersUnitTest.SqrtUnitTest();
            ok &= MemoryUnitTest.UnitTest(1001);
            ok &= CarrylessMultiplication.UnitTest(1001);

            int decimalDigits = 10100;
            int qwords = (int)Math.Ceiling(Math.Log(10, 2) * decimalDigits / 64);
            RealNumber pi = RealNumber.Zero;
            RealNumber e = RealNumber.Zero;
            RealNumber one = new RealNumber(1L, qwords);
            var computePITime = AsmX64Operations.MeasureTime(() =>
            {
                pi = RealNumber.GetPI(qwords);
            });
            var computeETime = AsmX64Operations.MeasureTime(() =>
            {
                e = one.GetExp();
            });
            string pivalue = "";
            string evalue = "";
            var baseConvertTimePI = AsmX64Operations.MeasureTime(() =>
            {
                pivalue = pi.ToString();
            });
            var baseConvertTimeE = AsmX64Operations.MeasureTime(() =>
            {
                evalue = e.ToString();
            });
            System.IO.File.WriteAllText(@"..\pi10k.txt", pivalue);
            System.IO.File.WriteAllText(@"..\e_10k.txt", evalue);

            MessageBox.Show(
                "Compute PI Time: " + computePITime.ToString() + "\r\n" +
                "Compute E  Time: " + computeETime.ToString() + "\r\n" +
                "base convert time PI: " + baseConvertTimePI.ToString() + "\r\n" +
                "base convert time E : " + baseConvertTimeE.ToString() + "\r\n",
                "Info", MessageBoxButtons.OK, MessageBoxIcon.Information);

            try
            {
                var no0 = new RealNumber(1L << 62, 2048 / 64);
                var no1 = no0.GetExp();
                var no2 = no1.GetLog();
                ok &= (no2 - no0).IsZero;

                ulong[] xx = Enumerable.Repeat(ulong.MaxValue, 1024 * 100).ToArray();
                ulong[] yy = Enumerable.Repeat(ulong.MaxValue, 1024 * 100).ToArray();
                ulong[] zz = Enumerable.Repeat((ulong)0, 1024 * 200).ToArray();
                AsmX64Operations.FastestMultiplication(xx, yy, zz, xx.Length);
                if (zz[0] != 1 || zz[zz.Length / 2 - 1] != 0 || zz[zz.Length / 2] + 2 != 0 || zz[zz.Length - 1] + 1 != 0)
                {
                    ok = false;
                }

                ok &= RealNumbersUnitTest.UnitTest(1001);
                ok &= FourierReal.UnitTest(10012);
                ok &= FourierMultiplication.UnitTest();
                ok &= IntegerNumberUnitTest.UnitTest();
                ok &= FourierMultiplication.UnitTestBigMul(1001);
                ok &= FastIntegerUnitTest.UnitTest();
                ok &= FourierMultiplication.UnitTest(1001);
                ok &= AsmX64Operations.ECCUnitTest();
                ok &= AsmX64Operations.UnitTest();
                ok &= HeapUnitTest.UnitTest();
                ok &= FourierTransform.FFTUnitTest();
                ok &= AccurateSummation.UnitTest();
                ok &= BinarySearchUnitTest.UnitTest();
                ok &= AVLTree<int>.UnitTest();
                ok &= AVLTreeSorted<int>.UnitTest();
            }
            catch (Exception ex)
            {
                MessageBox.Show(ex.ToString());
            }
            if (!ok)
            {
                MessageBox.Show("Unit tests failed.");
            }
            Application.EnableVisualStyles();
            Application.SetCompatibleTextRenderingDefault(false);
            Application.Run(new MainForm());
        }

        static void testMultiplicationSpeed()
        {
            Random random = new Random(1001);
            byte[] number = new byte[16];
            int n = 4;
            int iterations = 100000;
            while (true)
            {
                AsmX64Operations.SetKaratsubaBreakPoint(n / 2);
                ulong[] a = new ulong[n];
                ulong[] b = new ulong[n];
                ulong[] c = new ulong[n * 2];
                ulong[] tk = new ulong[AsmX64Operations.GetKaratsubaMultiplicationBufferSize(n)];
                for (int i = n; --i >= 0;)
                {
                    random.NextBytes(number);
                    a[i] = BitConverter.ToUInt64(number, 0);
                    b[i] = BitConverter.ToUInt64(number, 8);
                }
                var directTime = AsmX64Operations.MeasureTime(() =>
                {
                    for (int i = iterations; --i >= 0;)
                    {
                        AsmX64Operations.Multiply(a, b, c, n);
                    }
                });
                var karaTime = AsmX64Operations.MeasureTime(() =>
                {
                    for (int i = iterations; --i >= 0;)
                    {
                        AsmX64Operations.Karatsuba(a, b, c, n, tk);
                    }
                });
                if (karaTime.TotalSeconds < directTime.TotalSeconds)
                {
                    break;
                }
                n += 2;
            }
            MessageBox.Show("Kara break point = " + (n / 2).ToString());

            iterations = 300;
            AsmX64Operations.SetKaratsubaBreakPoint(14);
            n = 1 * 1024;
            while (true)
            {
                AsmX64Operations.SetFourierBreakPoint(n);
                ulong[] a = new ulong[n];
                ulong[] b = new ulong[n];
                ulong[] c = new ulong[n * 2];
                ulong[] tk = new ulong[AsmX64Operations.GetKaratsubaMultiplicationBufferSize(n)];
                for (int i = n; --i >= 0;)
                {
                    random.NextBytes(number);
                    a[i] = BitConverter.ToUInt64(number, 0);
                    b[i] = BitConverter.ToUInt64(number, 8);
                }
                var karaTime = AsmX64Operations.MeasureTime(() =>
                {
                    for (int i = iterations; --i >= 0;)
                    {
                        AsmX64Operations.Karatsuba(a, b, c, n, tk);
                    }
                });
                var fourierTime = AsmX64Operations.MeasureTime(() =>
                {
                    for (int i = iterations; --i >= 0;)
                    {
                        AsmX64Operations.FourierMultiplication(a, b, c, n);
                    }
                });
                if (fourierTime.TotalSeconds < karaTime.TotalSeconds)
                {
                    break;
                }
                n += 1024;
            }

            MessageBox.Show("Fourier break point = " + n.ToString());
        }
    }
}
