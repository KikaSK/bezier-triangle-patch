using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using OpenTK;

namespace MKP2___Template
{
    // enum type  determines the type of a patch
    public enum type
    {
        TENSOR,  // Bezier tensor product patch
        TRIANGLE, // Bezier triangle patch
        SPHERE, // one eighth of a sphere
        CYLINDER, // one fourth of a cylinder
        CONE // one fourth of a cone
    };

    // enum placement determines the position of a patch
    public enum placement 
    {
        LEFT,
        MIDDLE,
        RIGHT
    }

    class Patch
    {    
        // Class CNet describes the control vertices of a patch
        public class CNet
        {
            public List<int> Indices;
            public List<Vector3> Coordinates;
            public CNet()
            {
                Indices = new List<int>(); // Indices of a control net
                Coordinates = new List<Vector3>(); // Coordinates of vertices of the control net
            }
        }

        // Class IPoints describes the interpolated points of a patch
        public class IPoints
        {
            public List<Vector3> Coordinates;
            public IPoints()
            {
                Coordinates = new List<Vector3>(); // Coordinates of interpolated points
            }
        }

        // sampled patch (points computed using the de Casteljau algorithm)
        public class Sampl
        {
            public List<int> Indices;
            public List<Vector3> Coordinates;
            public List<Vector3> Normals;
            public Sampl()
            {
                Indices = new List<int>();
                Coordinates = new List<Vector3>();
                Normals = new List<Vector3>();
            }
        }

        public type TypeOfPatch; 
        public placement Place;
        public int NumberOfSamples, Degree;

        public CNet ControlNet;
        public IPoints InterpolatedPoints;
        public Sampl Sampling;
        public float[] Color;

        public bool? DisplayIsoCurves;

        public double hlS, hlT, hlU; // parameter S and T, respectively, for the isocurves

        public List<Vector3> IsoS_ControlPolygon, IsoT_ControlPolygon, IsoU_ControlPolygon;
        public List<Vector3> IsoS_Sampling, IsoT_Sampling, IsoU_Sampling;

        // Initialization of a patch
        public Patch(type _TypeOfPatch, int _Degree, int _NumberOfSamples, bool? _DisplayIsoCurves, double _hlS, double _hlT, double _hlU, float[] _Color, placement _Place)
        {
            TypeOfPatch = _TypeOfPatch;
            NumberOfSamples = _NumberOfSamples;
            Degree = _Degree;
            Color = _Color;
            Place = _Place;

            DisplayIsoCurves = _DisplayIsoCurves;
            hlS = _hlS;
            hlT = _hlT;
            hlU = _hlU;

            ControlNet = new CNet();
            Sampling = new Sampl();
            InterpolatedPoints = new IPoints();

            IsoS_ControlPolygon = new List<Vector3>();
            IsoT_ControlPolygon = new List<Vector3>();
            IsoU_ControlPolygon = new List<Vector3>();
            IsoS_Sampling = new List<Vector3>();
            IsoT_Sampling = new List<Vector3>();
            IsoU_Sampling = new List<Vector3>();

            // Initial sampling of a patch
            Sampling.Coordinates = Sample(TypeOfPatch, NumberOfSamples);
            Sampling.Indices = GetIndices(TypeOfPatch, NumberOfSamples, true);
            Sampling.Normals = GetNormals(TypeOfPatch, NumberOfSamples);

            // Initial control net
            ControlNet.Coordinates = Sample(TypeOfPatch, Degree);
            ControlNet.Indices = GetIndices(TypeOfPatch, Degree, true);            

            // Initial interpolated vertices
            InterpolatedPoints.Coordinates = Sample(TypeOfPatch, Degree);            
        }

        // the vertices of the triangle grid of the degree 3 are inserted in the following order
        //         9
        //        / \
        //       7---8
        //      / \ / \
        //     4---5---6
        //    / \ / \ / \
        //   0---1---2---3
        //

        // sampling of the initial patch with the given number of samples eU
        private List<Vector3> Sample(type _TypeOfPatch, int eU)
        {
            List<Vector3> SampleList = new List<Vector3>();

            // With respect to the order of insertion described above, sample the initial triangle
            // with vertices (-1, -1, 0)   (1, -1, 0)    (0, 1, 0) into the list "SampleList", i.e.
            // determine the coordanites of the sampled points

            Vector3 A = new Vector3(-1, -1, 0);
            Vector3 B = new Vector3(1, -1, 0);
            Vector3 C = new Vector3(0, 1, 0);

            for (int i = 0; i<=eU; ++i)
            {
                // Points on the edges in the ith storey of the triangle pyramid
                Vector3 P = A + 1.0f * i / eU * (C - A);
                Vector3 Q = B + 1.0f * i / eU * (C - B);

                for (int j = 0; j < (eU - i); ++j)
                {
                    Vector3 R = P + 1.0f * j / (eU - i) * (Q - P);
                    SampleList.Add(R);
                }
                SampleList.Add(Q);
            }

            return SampleList;
        }

        //                  0 0 1
        //
        //          1/3 0 2/3     0 1/3 2/3
        //                   
        //     2/3 0 1/3   1/3 1/3 1/3   0 2/3 1/3
        //   
        // 1 0 0    2/3 1/3 0     1/3 2/3 0     0 1 0
        //
        // get Grevill Abscissae for triangle patch
        private List<Vector3> GrevillAbscissae(int eU)
        {
            List<Vector3> Grevill = new List<Vector3>();

            float c;
            for(int i = 0; i<=eU; ++i)
            {
                c = 1.0f * i / eU;
                for (int j = 0; j <= (eU - i); ++j)
                {
                    float b = 1.0f * j / eU;
                    float a = 1.0f - b - c;

                    Grevill.Add(new Vector3(a, b, c));
                }
            }
            return Grevill;
        }

        // getting indicies for a patch with given number of samples eU, eV
        // in the direction u, v, respectively 
        private List<int> GetIndices(type _TypeOfPatch, int eU, bool DrawAll)
        {
            List<int> IndList = new List<int>();

            // With respect to the order of insertion described above,
            // create a sequence of indices, where each consecutive triplet
            // determines the indices of a triangle. The triangles may be
            // drawn in random order, however you should be consistent with
            // the orientation (the order of vertices within a triangle).
            // One way to do it is to create a sequence (for the degree 3)
            //    0 1 4 1 2 5 2 3 6 4 5 7 5 6 8 7 8 9 (1 5 4 2 6 5 5 8 7)
            // 

            // -----------------------------------------------
            //         9
            //        / \
            //       7---8
            //      / \ / \               22
            //     4---5---6            20  21
            //    / \ / \ / \         10  11  12
            //   0---1---2---3      00  01  02  03
            // -----------------------------------------------

            int pos = 0;
            for(int i = 0; i<eU; ++i) // rows
            {
                for(int j = 0; j<eU-i; ++j) // cols
                {
                    IndList.Add(pos);
                    IndList.Add(pos + 1);
                    IndList.Add(pos + eU - i + 1);
                    pos++;
                }
                pos++;
            }

            if(DrawAll)
            {
                pos = 0;
                for(int i=0; i<eU-1; ++i)
                {
                    for(int j=0; j<eU-i-1; ++j)
                    {
                        IndList.Add(pos + 1);
                        IndList.Add(pos + eU + 2 - i);
                        IndList.Add(pos + eU + 1 - i);
                        pos++;
                    }
                    pos+=2;
                }
            }
            return IndList;            
        }
         
        // initialize normals
        private List<Vector3> GetNormals(type _TypeOfPatch, int eU)
        {
            List<Vector3> NormalList = new List<Vector3>();
            // 1 + 2 + 3 + ... + eU ... (eU+1)*(eU+2)/2
            for(int i = 0; i<(eU+1)*(eU+2)/2; ++i)
            {
                NormalList.Add(new Vector3(0, 0, 1));
            }
            // For each sample (sampled point) initialize the normal vector to the value (0, 0, 1).
            
            return NormalList;
        }

        // Evaluate point by Casteljau for triangle patch
        public Tuple<Vector3, Vector3> EvaluateCasteljau(float s, float t, float u)
        {
            int n = Degree;
            List<Vector3> CP = new List<Vector3>(ControlNet.Coordinates);
            Vector3 normal = new Vector3();
            Vector3 point = new Vector3();
            while (n >= 1)
            {
                List<Vector3> newCP = new List<Vector3>();
                List<int> indices = GetIndices(TypeOfPatch, n, false);
                for (int j = 0; j < indices.Count; j += 3)
                {
                    Vector3 A = CP[indices[j]];
                    Vector3 B = CP[indices[j + 1]];
                    Vector3 C = CP[indices[j + 2]];
                    newCP.Add(s * A + t * B + u * C);
                }
                CP = new List<Vector3>(newCP);
                newCP.Clear();
                n--;
                if (n == 1)
                {
                    Vector3 a = CP[1] - CP[0];
                    Vector3 b = CP[2] - CP[0];
                    Vector3 my_normal = (Vector3.Cross(a, b)).Normalized();
                    normal = my_normal;
                }
            }
            point = CP[0];
            Tuple<Vector3, Vector3> res = new Tuple<Vector3, Vector3>(point, normal);
            return res;
        }

        // Computation of points of the patch 
        public void RecomputePatch()
        {

            if (TypeOfPatch == type.TRIANGLE)
            {
                // Get the coordinates of the control net of the INTERPOLATING BEZIER PATCH,
                // which are stored in the vector "ControlNet.Coordinates"

                // getting control net by multipication by inverted regular matrix of bezier polynoms
                // evaluated in grevill abscissae

                ControlNet.Coordinates.Clear();
                List<Vector3> GA = GrevillAbscissae(Degree);
                double[,] BezierMatrix = MyFunctions.CreateBezierMatrix(Degree, GA);
                alglib.rmatrixinverse(ref BezierMatrix, out int Info, out alglib.matinvreport Report);
                if (Info == -1)
                {
                    throw new Exception("ERROR: Computing inverse of singular matrix!");
                }
                else
                {
                    double[] heights = new double[GA.Count];
                    for (int i = 0; i < heights.Count(); ++i)
                        heights[i] = InterpolatedPoints.Coordinates[i].Z;

                    double[] res = MyFunctions.Multiply(ref BezierMatrix, ref heights, GA.Count());
                    ControlNet.Coordinates = new List<Vector3>(InterpolatedPoints.Coordinates);
                    for (int i = 0; i < res.Count(); ++i)
                        ControlNet.Coordinates[i] = new Vector3(ControlNet.Coordinates[i].X, ControlNet.Coordinates[i].Y, (float)res[i]);
                }
            }

            if (TypeOfPatch == type.SPHERE)
            {
                // Get the coordinates of the control net of the SPHERE
             
                //control points of sphere patch

                ControlNet.Coordinates.Clear();

                ControlNet.Coordinates.Add(new Vector3(1, -1, -1));
                ControlNet.Coordinates.Add(new Vector3(1, 1, -1));
                ControlNet.Coordinates.Add(new Vector3(-1, 1, -1));

                ControlNet.Coordinates.Add(new Vector3(1, -1, 1));
                ControlNet.Coordinates.Add(new Vector3(-1, 1, 1));

                ControlNet.Coordinates.Add(new Vector3(-1, -1, 1));          
            }           

            if (TypeOfPatch == type.CONE)
            {
                // Get the coordinates of the control net of the CONE
             
                //control points of cone patch

                ControlNet.Coordinates.Clear();

                ControlNet.Coordinates.Add(new Vector3(1, -1, -1));
                ControlNet.Coordinates.Add(new Vector3(1, 1, -1));
                ControlNet.Coordinates.Add(new Vector3(-1, 1, -1));

                ControlNet.Coordinates.Add(new Vector3(0, -1, 0));
                ControlNet.Coordinates.Add(new Vector3(-1, 0, 0));

                ControlNet.Coordinates.Add(new Vector3(-1, -1, 1));
            }

            // Using the control net compute the samples of the patch using the DE CASTELJAU ALGORITHM,
            // also for each sample compute the respective normal vector
            
            Sampling.Coordinates = Sample(TypeOfPatch, NumberOfSamples);
            for (int i = 0; i<Sampling.Coordinates.Count; ++i)
            {
                float s, t, u;

                // barycentric coordinates
                Tuple<double, double, double> BC = MyFunctions.GetBC(Sampling.Coordinates[i], new Vector3(-1, -1, 0),new Vector3(1, -1, 0),new Vector3(0, 1, 0));
                s = (float)BC.Item1;
                t = (float)BC.Item2;
                u = (float)BC.Item3;

                // evaluate using Casteljau
                Tuple<Vector3, Vector3> eval = EvaluateCasteljau(s, t, u);
                Sampling.Normals[i] = eval.Item2;
                Sampling.Coordinates[i] = eval.Item1;
            }

            if (DisplayIsoCurves == true)
            {
                // Compute the control polygon and sampled points of the S-isocurve, T-isocurve and U-isocurve
            
                IsoS_ControlPolygon.Clear();
                IsoT_ControlPolygon.Clear();
                IsoU_ControlPolygon.Clear();

                IsoS_Sampling.Clear();
                IsoT_Sampling.Clear();
                IsoU_Sampling.Clear();

                float s, t, u;

                //sampling of isocurves

                for (int i = 0; i < NumberOfSamples + 1; ++i)
                {
                    s = (float)hlS;
                    t = (float)i * (1.0f - s) / NumberOfSamples;
                    u = 1 - s - t;
                    Vector3 sampled_point = EvaluateCasteljau(s, t, u).Item1;
                    IsoS_Sampling.Add(sampled_point);

                    t = (float)hlT;
                    s = (float)i * (1.0f - t) / NumberOfSamples;
                    u = 1 - s - t;
                    sampled_point = EvaluateCasteljau(s, t, u).Item1;
                    IsoT_Sampling.Add(sampled_point);


                    u = (float)hlU;
                    t = (float)i * (1.0f - u) / NumberOfSamples;
                    s = 1 - u - t;
                    sampled_point = EvaluateCasteljau(s, t, u).Item1;
                    IsoU_Sampling.Add(sampled_point);
                }

                // control vertices of isocurves, again by inverting matrix

                double[] toInterpolateSx = new double[Degree + 1];
                double[] toInterpolateSy = new double[Degree + 1];
                double[] toInterpolateSz = new double[Degree + 1];
                double[] toInterpolateTx = new double[Degree + 1];
                double[] toInterpolateTy = new double[Degree + 1];
                double[] toInterpolateTz = new double[Degree + 1];
                double[] toInterpolateUx = new double[Degree + 1];
                double[] toInterpolateUy = new double[Degree + 1];
                double[] toInterpolateUz = new double[Degree + 1];

                Vector3 eval = new Vector3();

                // obtaining Degree+1 points on isocurve to interpolate
                for (int i = 0; i < Degree + 1; ++i)
                {
                    s = (float)hlS;
                    t = (float)i * (1.0f - s) / Degree;
                    u = 1 - s - t;
                    eval = EvaluateCasteljau(s, t, u).Item1;
                    toInterpolateSx[i] = eval.X;
                    toInterpolateSy[i] = eval.Y;
                    toInterpolateSz[i] = eval.Z;

                    t = (float)hlT;
                    s = (float)i * (1.0f - t) / Degree;
                    u = 1 - s - t;
                    eval = EvaluateCasteljau(s, t, u).Item1;
                    toInterpolateTx[i] = eval.X;
                    toInterpolateTy[i] = eval.Y;
                    toInterpolateTz[i] = eval.Z;


                    u = (float)hlU;
                    t = (float)i * (1.0f - u) / Degree;
                    s = 1 - u - t;
                    eval = EvaluateCasteljau(s, t, u).Item1;
                    toInterpolateUx[i] = eval.X;
                    toInterpolateUy[i] = eval.Y;
                    toInterpolateUz[i] = eval.Z;
                }

                // obtaining matrix
                double[,] matrix = MyFunctions.CreateCurveMatrix(Degree);
                alglib.rmatrixinverse(ref matrix, out int Info, out alglib.matinvreport Report);
                if (Info == -1)
                {
                    throw new Exception("ERROR: Computing inverse of singular matrix!");
                }

                // calculating control vertices
                double[] controlSx = MyFunctions.Multiply(ref matrix, ref toInterpolateSx, Degree + 1);
                double[] controlSy = MyFunctions.Multiply(ref matrix, ref toInterpolateSy, Degree + 1);
                double[] controlSz = MyFunctions.Multiply(ref matrix, ref toInterpolateSz, Degree + 1);
                double[] controlTx = MyFunctions.Multiply(ref matrix, ref toInterpolateTx, Degree + 1);
                double[] controlTy = MyFunctions.Multiply(ref matrix, ref toInterpolateTy, Degree + 1);
                double[] controlTz = MyFunctions.Multiply(ref matrix, ref toInterpolateTz, Degree + 1);
                double[] controlUx = MyFunctions.Multiply(ref matrix, ref toInterpolateUx, Degree + 1);
                double[] controlUy = MyFunctions.Multiply(ref matrix, ref toInterpolateUy, Degree + 1);
                double[] controlUz = MyFunctions.Multiply(ref matrix, ref toInterpolateUz, Degree + 1);

                for (int i = 0; i < Degree + 1; ++i)
                {
                    IsoS_ControlPolygon.Add(new Vector3((float)controlSx[i], (float)controlSy[i], (float)controlSz[i]));
                    IsoT_ControlPolygon.Add(new Vector3((float)controlTx[i], (float)controlTy[i], (float)controlTz[i]));
                    IsoU_ControlPolygon.Add(new Vector3((float)controlUx[i], (float)controlUy[i], (float)controlUz[i]));
                }

            }
        }

        // Affine trnasformaions of the control net
        public void Translate(float eX, float eY, float eZ)
        {
            for(int i = 0; i < ControlNet.Coordinates.Count; i++)
            {
                ControlNet.Coordinates[i] = new Vector3(ControlNet.Coordinates[i].X + eX, ControlNet.Coordinates[i].Y + eY, ControlNet.Coordinates[i].Z + eZ);
            }
        }

        public void Scale(float eX, float eY, float eZ)
        {
            for (int i = 0; i < ControlNet.Coordinates.Count; i++)
            {
                ControlNet.Coordinates[i] = new Vector3(ControlNet.Coordinates[i].X * eX, ControlNet.Coordinates[i].Y * eY, ControlNet.Coordinates[i].Z * eZ);
            }
        }
    }
}
