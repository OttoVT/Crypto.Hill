using System;
using System.Collections.Generic;
using System.Text;

namespace Crypto.Hill.Core.Matricies
{
    public class Matrix//<T> where T : System.Double
    {
        private readonly double[,] _matrix;

        public double this[int x, int y]
        {
            get { return _matrix[x, y]; }
        }

        public int DimensionX
        {
            get { return _matrix.GetLength(0); }
        }

        public int DimensionY
        {
            get { return _matrix.GetLength(1); }
        }

        public Matrix(int dimensionX, int dimensionY)
        {
            _matrix = MatrixCreate(dimensionX, dimensionY);
            var random = new Random();

            for (int i = 0; i < dimensionX; i++)
            {
                for (int j = 0; j < DimensionY; j++)
                {
                    _matrix[i, j] = random.Next(0, int.MaxValue - 1);
                }
            }
        }

        public Matrix(double[,] matrix)
        {
            _matrix = matrix;
        }

        public double Determinant()
        {
            var result = MatrixDeterminant(this._matrix);

            return result;
        }

        public Matrix Mod(double n)
        {
            for (int i = 0; i < _matrix.GetLength(0); i++)
            {
                for (int j = 0; j < _matrix.GetLength(1); j++)
                {
                    var calculated = _matrix[i, j] % n;
                    _matrix[i, j] = (calculated < 0 ? calculated + n : calculated) % n;
                }
            }

            return this;
        }

        public Matrix MatrixInverse()
        {
            double[,] matrix = _matrix;
            // assumes determinant is not 0
            // that is, the matrix does have an inverse
            int n = matrix.GetLength(0);
            double[,] result = MatrixCreate(n, n); // make a copy of matrix
            for (int i = 0; i < n; ++i)
            for (int j = 0; j < n; ++j)
                result[i, j] = matrix[i, j];

            double[,] lum; // combined lower & upper
            int[] perm;
            int toggle;
            toggle = MatrixDecompose(matrix, out lum, out perm);

            double[] b = new double[n];
            for (int i = 0; i < n; ++i)
            {
                for (int j = 0; j < n; ++j)
                    if (i == perm[j])
                        b[j] = 1.0;
                    else
                        b[j] = 0.0;

                double[] x = Helper(lum, b); // 
                for (int j = 0; j < n; ++j)
                    result[j, i] = x[j];
            }
            return new Matrix(result);
        } // MatrixInverse

        public static Matrix operator *(Matrix b, Matrix c)
        {
            var newMatrix = new Matrix(c.DimensionX, b.DimensionX);

            for (int i = 0; i < newMatrix.DimensionX; i++)
            {
                double summ = 0;
                for (int j = 0; j < newMatrix.DimensionY; j++)
                {
                    newMatrix._matrix[i, j] = GetMatrixMultiplication(i, j, c, b);
                }
            }

            return newMatrix;
        }

        static int MatrixDecompose(double[,] m, out double[,] lum, out int[] perm)
        {
            // Crout's LU decomposition for matrix determinant and inverse
            // stores combined lower & upper in lum[][]
            // stores row permuations into perm[]
            // returns +1 or -1 according to even or odd number of row permutations
            // lower gets dummy 1.0s on diagonal (0.0s above)
            // upper gets lum values on diagonal (0.0s below)

            int toggle = +1; // even (+1) or odd (-1) row permutatuions
            int n = m.GetLength(0);

            // make a copy of m[][] into result lu[][]
            lum = MatrixCreate(n, n);
            for (int i = 0; i < n; ++i)
            for (int j = 0; j < n; ++j)
                lum[i,j] = m[i,j];


            // make perm[]
            perm = new int[n];
            for (int i = 0; i < n; ++i)
                perm[i] = i;

            for (int j = 0; j < n - 1; ++j) // process by column. note n-1 
            {
                double max = Math.Abs(lum[j,j]);
                int piv = j;

                for (int i = j + 1; i < n; ++i) // find pivot index
                {
                    double xij = Math.Abs(lum[i,j]);
                    if (xij > max)
                    {
                        max = xij;
                        piv = i;
                    }
                } // i

                if (piv != j)
                {
                    double[] tmp = new double[lum.GetLength(1)];

                    for (int i = 0; i < lum.GetLength(1); i++)
                    {
                        tmp[i] = lum[piv, i];
                    }

                    for (int i = 0; i < lum.GetLength(1); i++)
                    {
                        lum[piv, i] = lum[j, i];
                    }

                    for (int i = 0; i < lum.GetLength(1); i++)
                    {
                        lum[j, i] = tmp[i];
                    }

                    int t = perm[piv]; // swap perm elements
                    perm[piv] = perm[j];
                    perm[j] = t;

                    toggle = -toggle;
                }

                double xjj = lum[j,j];
                if (xjj != 0.0)
                {
                    for (int i = j + 1; i < n; ++i)
                    {
                        double xij = lum[i,j] / xjj;
                        lum[i,j] = xij;
                        for (int k = j + 1; k < n; ++k)
                            lum[i,k] -= xij * lum[j,k];
                    }
                }

            } // j

            return toggle;
        } // MatrixDecompose

        static double[] Helper(double[,] luMatrix, double[] b) // helper
        {
            int n = luMatrix.GetLength(0);
            double[] x = new double[n];
            b.CopyTo(x, 0);

            for (int i = 1; i < n; ++i)
            {
                double sum = x[i];
                for (int j = 0; j < i; ++j)
                    sum -= luMatrix[i,j] * x[j];
                x[i] = sum;
            }

            x[n - 1] /= luMatrix[n - 1,n - 1];
            for (int i = n - 2; i >= 0; --i)
            {
                double sum = x[i];
                for (int j = i + 1; j < n; ++j)
                    sum -= luMatrix[i,j] * x[j];
                x[i] = sum / luMatrix[i,i];
            }

            return x;
        } // Helper

        static double MatrixDeterminant(double[,] matrix)
        {
            double[,] lum;
            int[] perm;
            int toggle = MatrixDecompose(matrix, out lum, out perm);
            double result = toggle;
            for (int i = 0; i < lum.GetLength(0); ++i)
                result *= lum[i,i];
            return result;
        }

        // ----------------------------------------------------------------

        static double[,] MatrixCreate(int rows, int cols)
        {
            double[,] result = new double[rows, cols];

            return result;
        }

        static double[,] MatrixProduct(double[,] matrixA,
          double[,] matrixB)
        {
            int aRows = matrixA.GetLength(1);
            int aCols = matrixA.GetLength(0);
            int bRows = matrixB.GetLength(1);
            int bCols = matrixB.GetLength(0);
            if (aCols != bRows)
                throw new Exception("Non-conformable matrices");

            double[,] result = MatrixCreate(aRows, bCols);

            for (int i = 0; i < aRows; ++i) // each row of A
                for (int j = 0; j < bCols; ++j) // each col of B
                    for (int k = 0; k < aCols; ++k) // could use k < bRows
                        result[i,j] += matrixA[i,k] * matrixB[k,j];

            return result;
        }

        public override string ToString()
        {
            string s = "";
            for (int i = 0; i < _matrix.GetLength(0); ++i)
            {
                for (int j = 0; j < _matrix.GetLength(1); ++j)
                    s += _matrix[i,j].ToString("F3").PadLeft(8) + " ";
                s += Environment.NewLine;
            }
            return s;
        }

        static double[,] ExtractLower(double[,] lum)
        {
            // lower part of an LU Doolittle decomposition (dummy 1.0s on diagonal, 0.0s above)
            int n = lum.GetLength(0);
            double[,] result = MatrixCreate(n, n);
            for (int i = 0; i < n; ++i)
            {
                for (int j = 0; j < n; ++j)
                {
                    if (i == j)
                        result[i,j] = 1.0;
                    else if (i > j)
                        result[i,j] = lum[i,j];
                }
            }
            return result;
        }

        static double[,] ExtractUpper(double[,] lum)
        {
            // upper part of an LU (lu values on diagional and above, 0.0s below)
            int n = lum.GetLength(0);
            double[,] result = MatrixCreate(n, n);
            for (int i = 0; i < n; ++i)
            {
                for (int j = 0; j < n; ++j)
                {
                    if (i <= j)
                        result[i,j] = lum[i,j];
                }
            }
            return result;
        }


        private static double GetMatrixMultiplication(int column, int row, Matrix b, Matrix c)
        {
            double summ = 0;
            for (int i = 0; i < b.DimensionY; i++)
            {
                summ += b._matrix[column, i] * c._matrix[row, i];
            }

            return summ;
        }


    }
}
