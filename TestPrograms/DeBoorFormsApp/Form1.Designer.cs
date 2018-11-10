namespace DeBoorFormsApp
{
    partial class Form1
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
            this.SuspendLayout();
            // 
            // ctrlPaint
            // 
            this.ctrlPaint.Dock = System.Windows.Forms.DockStyle.Fill;
            this.ctrlPaint.Location = new System.Drawing.Point(0, 0);
            this.ctrlPaint.Name = "ctrlPaint";
            this.ctrlPaint.Size = new System.Drawing.Size(936, 518);
            this.ctrlPaint.TabIndex = 0;
            // 
            // Form1
            // 
            this.AutoScaleDimensions = new System.Drawing.SizeF(6F, 13F);
            this.AutoScaleMode = System.Windows.Forms.AutoScaleMode.Font;
            this.ClientSize = new System.Drawing.Size(936, 518);
            this.Controls.Add(this.ctrlPaint);
            this.Name = "Form1";
            this.Text = "Form1";
            this.ResumeLayout(false);

        }

        #endregion
        private Utilities.ctrlPaint ctrlPaint;
    }
}

