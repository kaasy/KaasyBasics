namespace Utilities
{
    partial class MainForm
    {
        /// <summary>
        /// Required designer variable.
        /// </summary>
        private System.ComponentModel.IContainer components = null;

        /// <summary>
        /// Clean up any resources being used.
        /// </summary>
        /// <param name="disposing">true if managed resources should be disposed; otherwise, false.</param>
        protected override void Dispose(bool disposing)
        {
            if (disposing && (components != null))
            {
                components.Dispose();
            }
            base.Dispose(disposing);
        }

        #region Windows Form Designer generated code

        /// <summary>
        /// Required method for Designer support - do not modify
        /// the contents of this method with the code editor.
        /// </summary>
        private void InitializeComponent()
        {
            this.btnTestArithmeticOPS = new System.Windows.Forms.Button();
            this.pbarOPS = new System.Windows.Forms.ProgressBar();
            this.SuspendLayout();
            // 
            // btnTestArithmeticOPS
            // 
            this.btnTestArithmeticOPS.Location = new System.Drawing.Point(13, 13);
            this.btnTestArithmeticOPS.Name = "btnTestArithmeticOPS";
            this.btnTestArithmeticOPS.Size = new System.Drawing.Size(122, 23);
            this.btnTestArithmeticOPS.TabIndex = 0;
            this.btnTestArithmeticOPS.Text = "Test Arithmetic OPS";
            this.btnTestArithmeticOPS.UseVisualStyleBackColor = true;
            this.btnTestArithmeticOPS.Click += new System.EventHandler(this.btnTestArithmeticOPS_Click);
            // 
            // pbarOPS
            // 
            this.pbarOPS.Anchor = ((System.Windows.Forms.AnchorStyles)(((System.Windows.Forms.AnchorStyles.Top | System.Windows.Forms.AnchorStyles.Left) 
            | System.Windows.Forms.AnchorStyles.Right)));
            this.pbarOPS.Location = new System.Drawing.Point(13, 42);
            this.pbarOPS.Name = "pbarOPS";
            this.pbarOPS.Size = new System.Drawing.Size(775, 23);
            this.pbarOPS.TabIndex = 1;
            // 
            // MainForm
            // 
            this.AutoScaleDimensions = new System.Drawing.SizeF(6F, 13F);
            this.AutoScaleMode = System.Windows.Forms.AutoScaleMode.Font;
            this.ClientSize = new System.Drawing.Size(800, 450);
            this.Controls.Add(this.pbarOPS);
            this.Controls.Add(this.btnTestArithmeticOPS);
            this.Name = "MainForm";
            this.Text = "Main Form";
            this.ResumeLayout(false);

        }

        #endregion

        private System.Windows.Forms.Button btnTestArithmeticOPS;
        private System.Windows.Forms.ProgressBar pbarOPS;
    }
}

