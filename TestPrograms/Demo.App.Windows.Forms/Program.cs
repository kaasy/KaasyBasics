using Utilities;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Windows.Forms;
using System.IO;
using System.Diagnostics;
using System.Threading.Tasks;

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
            //testMultiplicationSpeed();
            Application.EnableVisualStyles();
            Application.SetCompatibleTextRenderingDefault(false);
            Application.Run(new MainForm());
        }

        static void testMultiplicationSpeed()
        {
            Random random = new Random(1001);
            byte[] number = new byte[16];
            int n = 4;
            int iterations = 1000000;
            while (true)
            {
                AsmX64Operations.SetKaratsubaBreakPoint(n);
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
                n++;
            }
            MessageBox.Show("Kara break point = " + n.ToString());

            iterations = 100;
            //AsmX64Operations.SetKaratsubaBreakPoint(14);
            n = 4 * 1024;
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
                n += 512;
            }

            MessageBox.Show("Fourier break point = " + n.ToString());
        }
    }
}
