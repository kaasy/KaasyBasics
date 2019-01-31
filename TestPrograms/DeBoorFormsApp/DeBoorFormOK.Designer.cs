namespace DeBoorFormsApp
{
    partial class DeBoorFormOK
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
            this.ctrlPaint = new Utilities.ctrlPaint();
            this.chkPeriodic = new System.Windows.Forms.CheckBox();
            this.SuspendLayout();
            // 
            // ctrlPaint
            // 
            this.ctrlPaint.Anchor = ((System.Windows.Forms.AnchorStyles)((((System.Windows.Forms.AnchorStyles.Top | System.Windows.Forms.AnchorStyles.Bottom) 
            | System.Windows.Forms.AnchorStyles.Left) 
            | System.Windows.Forms.AnchorStyles.Right)));
            this.ctrlPaint.Degree = 3;
            this.ctrlPaint.Periodic = false;
            this.ctrlPaint.Location = new System.Drawing.Point(0, 32);
            this.ctrlPaint.Name = "ctrlPaint";
            this.ctrlPaint.Size = new System.Drawing.Size(676, 440);
            this.ctrlPaint.TabIndex = 0;
            // 
            // chkPeriodic
            // 
            this.chkPeriodic.AutoSize = true;
            this.chkPeriodic.Location = new System.Drawing.Point(12, 9);
            this.chkPeriodic.Name = "chkPeriodic";
            this.chkPeriodic.Size = new System.Drawing.Size(94, 17);
            this.chkPeriodic.TabIndex = 1;
            this.chkPeriodic.Text = "Periodic spline";
            this.chkPeriodic.UseVisualStyleBackColor = true;
            this.chkPeriodic.CheckedChanged += new System.EventHandler(this.chkPeriodic_CheckedChanged);
            // 
            // DeBoorFormOK
            // 
            this.AutoScaleDimensions = new System.Drawing.SizeF(6F, 13F);
            this.AutoScaleMode = System.Windows.Forms.AutoScaleMode.Font;
            this.ClientSize = new System.Drawing.Size(676, 472);
            this.Controls.Add(this.chkPeriodic);
            this.Controls.Add(this.ctrlPaint);
            this.Name = "DeBoorFormOK";
            this.Text = "De Boor Form OK - BSpline - Click to add points";
            this.Load += new System.EventHandler(this.DeBoorFormOK_Load);
            this.ResumeLayout(false);
            this.PerformLayout();

        }

        #endregion

        private Utilities.ctrlPaint ctrlPaint;
        private System.Windows.Forms.CheckBox chkPeriodic;
    }
}

