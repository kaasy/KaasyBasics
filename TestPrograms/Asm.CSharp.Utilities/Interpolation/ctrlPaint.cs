using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Drawing;
using System.Data;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows.Forms;
using System.Numerics;

namespace Utilities
{
    public partial class ctrlPaint : UserControl
    {
        public ctrlPaint()
        {
            InitializeComponent();
            this.Degree = 3;
            this.LocalSuppot = true;
        }

        public int Degree { get; set; }
        public bool LocalSuppot { get; set; }

        private void ctrlPaint_Paint(object sender, PaintEventArgs e)
        {
            try
            {
                Redraw(e);
            }
            catch (Exception ex)
            {
                this.toolStripStatusLabel1.Text = ex.ToString();
            }
        }

        private List<Complex> pointsXY = new List<Complex>();

        private static double interpolate(double x1, double x, double x2, double y1, double y2)
        {
            return (x - x1) / (x2 - x1) * (y2 - y1) + y1;
        }

        private List<double> getDistances()
        {
            List<double> distances = new List<double>();
            double firstDistance = 0;
            distances.Add(firstDistance);
            for (int i = 1; i < this.pointsXY.Count; i++)
            {
                firstDistance += this.pointsXY[i - 1].Distance(this.pointsXY[i]);
                distances.Add(firstDistance);
            }
            if (this.pointsXY.Count > 1)
            {
                firstDistance += this.pointsXY[0].Distance(this.pointsXY[this.pointsXY.Count - 1]);
                distances.Add(firstDistance);
            }
            return distances;
        }

        private void Redraw(PaintEventArgs e)
        {
            int width = this.pictureBox1.Width;
            int height = this.pictureBox1.Height;
            e.Graphics.Clear(Color.Black);
            if (width <= 1 || height <= 1)
            {
                return;
            }
            Graphics G = e.Graphics;
            Pen thePen = new Pen(Color.White, 2.0f);

            var distances = this.getDistances();
            List<PointF> list = new List<PointF>();
            GeneralSpline spline;
            if (this.LocalSuppot)
            {
                spline = new DeBoorAlgorithm(distances, this.pointsXY, this.Degree);
            }
            else
            {
                spline = new PeriodicCubicSpline(distances, this.pointsXY);
                double err = (spline as PeriodicCubicSpline).UnitTestError();
                if (err > 1E-6)
                {
                    System.Diagnostics.Debugger.Break();
                }
            }

            int oldDistance = 0;
            for (int i = 1; i < distances.Count; i++)
            {
                double currentDistance = distances[i];
                while (oldDistance < currentDistance)
                {
                    Complex ev = spline.Evaluate(oldDistance, i - 1);
                    list.Add(new PointF((float)ev.Real, (float)ev.Imaginary));
                    oldDistance++;
                }
            }
            if (list.Count >= 2)
            {
                e.Graphics.DrawLines(thePen, list.ToArray());
            }
        }

        private void pictureBox1_MouseClick(object sender, MouseEventArgs e)
        {
            HashSet<Complex> toRemove = new HashSet<Complex>();
            foreach (Complex point in this.pointsXY)
            {
                double distanceSquared = (point.Real - e.X).Square() + (point.Imaginary - e.Y).Square();
                if (distanceSquared < 200)
                {
                    toRemove.Add(point);
                }
            }
            this.pointsXY.RemoveAll(p => toRemove.Contains(p));
            this.pointsXY.Add(new Complex(e.X, e.Y));
            this.pictureBox1.Invalidate();
        }
    }
}