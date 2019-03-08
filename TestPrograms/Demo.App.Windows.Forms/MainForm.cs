using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Data;
using System.Drawing;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows.Forms;

namespace Utilities
{
    public partial class MainForm : Form
    {
        public MainForm()
        {
            InitializeComponent();
            this.pbarOPS.Enabled = false;
        }

        private static int differences(string s1, string s2)
        {
            int diffs = 0;
            for (int i = Math.Min(s1.Length, s2.Length) - 1; --i >= 0;)
            {
                diffs += s1[i] != s2[i] ? 1 : 0;
            }
            return diffs;
        }

        private void invoke(Action action)
        {
            if (action == null)
            {
                return;
            }
            if (this.InvokeRequired)
            {
                this.Invoke(action);
            }
            else
            {
                action();
            }
        }

        private void testAll()
        {
            int maxDigits = 1000;
            long totalOps = (long)maxDigits * (maxDigits + 1) / 2;
            long ops = 0;
            for (int digits = 1; digits < maxDigits; digits++)
            {
                RealNumber b_one = new RealNumber(2L, digits);
                RealNumber log1 = RealNumber.Zero, log2 = RealNumber.Zero;
                Parallel.Invoke(
                () =>
                {
                    log1 = b_one.GetLog();
                },
                () =>
                {
                    log2 = b_one.GetTaylorLog();
                });
                string log1String = log1.ToString();
                string log2String = log2.ToString();
                int diffs = differences(log2String, log1String);
                if (diffs > 0)
                {
                    MessageBox.Show("Error: Differences count = " + diffs.ToString("N3"));
                    break;
                }
                ops += digits;
                this.invoke(() =>
                {
                    this.pbarOPS.Value = Math.Min(Math.Max((int)(100.0 * ops / totalOps), 0), 100);
                });
            }
            this.invoke(() =>
            {
                this.pbarOPS.Value = 0;
                this.pbarOPS.Enabled = false;
                this.btnTestArithmeticOPS.Enabled = true;
            });
        }

        private void btnTestArithmeticOPS_Click(object sender, EventArgs e)
        {
            Task.Factory.StartNew(() =>
            {
                this.testAll();
            });
            this.btnTestArithmeticOPS.Enabled = false;
        }
    }
}
