using System;
using System.Collections.Generic;
using System.Text;
using Crypto.Hill.Core.Matricies;

namespace Crypto.Hill
{
    class Program
    {
        static void Main(string[] args)
        {
            Console.WriteLine((int)'a');
            Console.WriteLine((int)'A');
            Matrix keyMatrix = new Matrix(new double[2,2] {{1, 3}, {2, 7}});
            var determinant = keyMatrix.Determinant();
            Matrix mInverse = keyMatrix.MatrixInverse();
            var alphabetLength = 27;
            
            Matrix toEncrypt = new Matrix(new double[1,2]{{1, 14}}); //BO
            Matrix encrypted = (keyMatrix * toEncrypt).Mod(alphabetLength);
            Console.WriteLine(encrypted.ToString());
            Matrix dencrypted = (encrypted * mInverse).Mod(alphabetLength);
            Console.WriteLine(dencrypted.ToString());

            Console.Write("Input string to encrypt: ");
            string toEncr = Console.ReadLine();

            var groupSize =  keyMatrix.DimensionX;
            toEncr = toEncr.Length % 2 != 0 ? toEncr+'@' : toEncr;
            List<Matrix> groupedInput = new List<Matrix> ();
            var nextMatrix = new double[1, groupSize];
            var groupCounter = 0;
            for (int i = 0; i < toEncr.Length; i++, groupCounter++)
            {
                if (groupCounter == groupSize)
                {
                    groupedInput.Add(new Matrix(nextMatrix));
                    nextMatrix = new double[1, groupSize];

                    groupCounter = 0;
                }

                nextMatrix[0, groupCounter] = toEncr[i] - 'A';
            }

            groupedInput.Add(new Matrix(nextMatrix));
            StringBuilder sb = new StringBuilder();

            List<Matrix> encryptedMatrices = new List<Matrix>();
            foreach (var groupMatrix in groupedInput)
            {
                var encryptedMatrix = (keyMatrix * groupMatrix).Mod(alphabetLength);
                encryptedMatrices.Add(encryptedMatrix);
                for (int i = 0; i < groupSize; i++)
                {
                    sb.Append((char)(encryptedMatrix[0, i] + 'A'));
                }
            }

            Console.WriteLine(sb.ToString());

            sb = new StringBuilder();

            foreach (var matrix in encryptedMatrices)
            {
                var decryptedMatrix = (matrix * mInverse).Mod(alphabetLength);

                for (int i = 0; i < groupSize; i++)
                {
                    sb.Append((char)(decryptedMatrix[i, 0] + 'A'));
                }

                if (sb[sb.Length - 1] == '@')
                    sb.Remove(sb.Length - 1, 1);
            }

            Console.WriteLine(sb.ToString());

            Console.ReadKey();
        }
    }
}
