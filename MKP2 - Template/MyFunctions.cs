using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using OpenTK;
using System.Windows;

namespace MKP2___Template
{
    static class MyFunctions
    {
        private static int factorial(int n)
        {
            int res = 1;
            for(int i = 2; i<=n; ++i)
                res *= i;
            return res;
        }

        // combination number for triangle patch
        private static int comb2(int n, int i, int j, int k)
        {
            return factorial(n) / (factorial(i) * factorial(j) * factorial(k));
        }
        // Bernstein polynomial for triangle patch
        private static double Bernstein2(int n, int i, int j, int k, double s, double t, double u)
        {
            return comb2(n, i, j, k) * Math.Pow(s, i) * Math.Pow(t, j) * Math.Pow(u, k);
        }


        // GA : a, b, c, d, e, f
        // B200 B110 B020 B101 B011 B002 (a)   p200      f200
        // B200 B110 B020 B101 B011 B002 (b)   p110      f110
        // B200 B110 B020 B101 B011 B002 (c) * p020   =  f020
        // B200 B110 B020 B101 B011 B002 (d)   p101      f101
        // B200 B110 B020 B101 B011 B002 (e)   p011      f011
        // B200 B110 B020 B101 B011 B002 (f)   p002      f002

        // matrix to invert for triangle patch
        public static double[,] CreateBezierMatrix(int Degree, List<Vector3> GA)
        {
            // 300 210 120 030    201 111 021    102 012    003
            double[,] matrix = new double[GA.Count, GA.Count];
            for (int l = 0; l<GA.Count; ++l)
            {
                int counter = 0;
                for (int i = 0; i <= Degree; ++i)
                {
                    for (int j = 0; j <= Degree - i; ++j)
                    {
                        matrix[l, counter] = Bernstein2(Degree, Degree - i - j, j, i, GA[l].X, GA[l].Y, GA[l].Z);
                        counter++;
                    }
                }
            }
            return matrix;
            
        }

        // combination number for isocurves
        private static int comb(int n, int k)
        {
            return factorial(n) / (factorial(k) * factorial(n-k));
        }
        // Bernstein polynomial for isocurves
        private static double Bernstein(int n, int k, double a)
        {
            return comb(n, k) * Math.Pow(a, k) * Math.Pow(1-a, n-k);
        }

        // matrix to invert for isocurves
        public static double[,] CreateCurveMatrix(int Degree)
        {
            double[,] matrix = new double[Degree+1, Degree+1];
            //List<float> GA = new List<float>();
            for (int i = 0; i < Degree + 1; ++i)
            {
                for (int j = 0; j < Degree + 1; ++j)
                {
                    matrix[i, j] = Bernstein(Degree, j, (float)i / Degree);
                }
            }
            return matrix;
        }

        // barycentric coordinates od point P
        public static Tuple<double, double, double> GetBC(Vector3 P, Vector3 A, Vector3 B, Vector3 C)
        {
            Point _P = new Point(P.X, P.Y);
            Point _A = new Point(A.X, A.Y);
            Point _B = new Point(B.X, B.Y);
            Point _C = new Point(C.X, C.Y);
            Vector AB = _B - _A;
            Vector AC = _C - _A;
            Vector PA = _A - _P;
            Vector PB = _B - _P;
            Vector PC = _C - _P;

            double areaABC = Vector.CrossProduct(AB, AC);
            double areaPBC = Vector.CrossProduct(PB, PC);
            double areaPCA = Vector.CrossProduct(PC, PA);

            double s = areaPBC / areaABC;
            double t = areaPCA / areaABC;
            double u = 1 - s - t;

            return new Tuple<double,double,double>(s, t, u);
        }
        // Multiply 2 matrices
        public static double[] Multiply(ref double[,] M, ref double[] v, int S)
        {
            double[] res = new double[S];
            for (int i = 0; i < S; i++)
            {
                res[i] = 0;
                for (int j = 0; j < S; j++)
                {
                    res[i] += M[i, j] * v[j];
                }
            }
            return res;
        }
    }

}
