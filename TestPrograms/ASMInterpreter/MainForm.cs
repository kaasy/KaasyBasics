using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Data;
using System.Drawing;
using System.IO;
using System.Linq;
using System.Text;
using System.Windows.Forms;
using Utilities;

namespace ASMInterpreter
{
    public partial class MainForm : Form
    {
        public MainForm()
        {
            InitializeComponent();
        }

        private void MainForm_Load(object sender, EventArgs e)
        {
            string inputFileName = Environment.GetCommandLineArgs()[1];
            string asmText = File.ReadAllText(inputFileName);
            ASMDecomposer asm = new ASMDecomposer(asmText);
            asm.Execute();
        }
    }
}
