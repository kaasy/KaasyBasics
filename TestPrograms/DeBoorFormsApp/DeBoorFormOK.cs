using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Data;
using System.Drawing;
using System.Linq;
using System.Text;
using System.Windows.Forms;
using Utilities;

namespace DeBoorFormsApp
{
    public partial class DeBoorFormOK : Form
    {
        public DeBoorFormOK()
        {
            InitializeComponent();
        }

        private void chkPeriodic_CheckedChanged(object sender, EventArgs e)
        {
            this.ctrlPaint.Periodic = this.chkPeriodic.Checked;
            this.ctrlPaint.CompleteRedraw();
        }

        private void DeBoorFormOK_Load(object sender, EventArgs e)
        {
            chkPeriodic_CheckedChanged(sender, e);
        }
    }


}
