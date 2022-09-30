namespace MathCore;

/// <summary>
/// Class <c>Matrix</c> representing a rectangular matrix.
/// <remarks>
/// T are the types that implement the INumber interface.
/// This list of types includes the unmanaged.
/// </remarks>
/// </summary>
public class Matrix<T> where T : INumber<T>
{
    private readonly T[][] _storage;
    public int Rows { get; }
    public int Columns { get; }

    public T this[int i, int j]
    {
        get => _storage[i][j];
        set => _storage[i][j] = value;
    }

    public Matrix(int rows, int columns)
    {
        Rows = rows;
        Columns = columns;
        _storage = new T[rows].Select(_ => new T[columns]).ToArray();
    }

    /// <summary>
    /// Method <c>Copy</c> with a new matrix object returned.
    /// </summary>
    /// <param name="otherMatrix">A matrix from which the values are copied.</param>
    /// <returns>A copied matrix.</returns>
    public static Matrix<T> Copy(Matrix<T> otherMatrix)
    {
        Matrix<T> newMatrix = new(otherMatrix.Rows, otherMatrix.Columns);

        for (int i = 0; i < otherMatrix.Rows; i++)
        {
            for (int j = 0; j < otherMatrix.Columns; j++)
            {
                newMatrix[i, j] = otherMatrix[i, j];
            }
        }

        return newMatrix;
    }
}

/// <summary>
/// Class <c>SquareMatrix</c> represents a square matrix.
/// <remarks>
/// T are the types that implement the INumber interface.
/// This list of types includes the unmanaged.
/// </remarks>
/// </summary>
public class SquareMatrix<T> where T : INumber<T>
{
    private readonly T[,] _storage;
    public int Size { get; }

    public T this[int i, int j]
    {
        get => _storage[i, j];
        set => _storage[i, j] = value;
    }

    public SquareMatrix(int size)
        => (Size, _storage) = (size, new T[size, size]);

    public void Clear()
        => Array.Clear(_storage, 0, _storage.Length);

    /// <summary>
    /// Method <c>Copy</c> with a new matrix object returned.
    /// </summary>
    /// <param name="otherMatrix">A matrix from which the values are copied.</param>
    /// <returns>A copied matrix.</returns>
    public static SquareMatrix<T> Copy(SquareMatrix<T> otherMatrix)
    {
        SquareMatrix<T> newMatrix = new(otherMatrix.Size);

        for (int i = 0; i < otherMatrix.Size; i++)
        {
            for (int j = 0; j < otherMatrix.Size; j++)
            {
                newMatrix[i, j] = otherMatrix[i, j];
            }
        }

        return newMatrix;
    }
}

/// <summary>
/// Class <c>SparseMatrix</c> representing a sparse matrix.
/// </summary>
public class SparseMatrix
{
    /// <summary>
    /// Property <c>Ig</c> contains the indexes from which
    /// the elements of the k-th row (column) begin in the arrays Jg, Ggl and Ggu.
    /// </summary>
    public int[] Ig { get; set; }

    /// <summary>
    /// Property <c>Jg</c> contains the numbers of columns (rows) of the stored off-diagonal
    /// elements of the lower (upper) triangle matrix.
    /// </summary>
    public int[] Jg { get; set; }

    /// <summary>
    /// Property <c>Di</c> contains
    /// the diagonal elements of the matrix.
    /// </summary>
    public double[] Di { get; set; }

    /// <summary>
    /// Property <c>Ggl</c> to store off-diagonal
    /// elements of the lower triangle matrix by rows.
    /// </summary>
    public double[] Ggl { get; set; }

    /// <summary>
    /// Property <c>Ggu</c> to store off-diagonal
    /// elements of the upper triangle matrix by columns.
    /// </summary>
    public double[] Ggu { get; set; }

    public int Size { get; }

    /// <summary>
    /// Constructor of sparse matrix.
    /// </summary>
    /// <param name="size">Matrix size.</param>
    /// <param name="sizeOffDiag">Number of off-diagonal elements.</param>
    public SparseMatrix(int size, int sizeOffDiag)
    {
        Size = size;
        Ig = new int[size + 1];
        Jg = new int[sizeOffDiag];
        Ggl = new double[sizeOffDiag];
        Ggu = new double[sizeOffDiag];
        Di = new double[size];
    }

    public static Vector<double> operator *(SparseMatrix matrix, Vector<double> vector)
    {
        Vector<double> product = new(vector.Length);

        for (int i = 0; i < vector.Length; i++)
        {
            product[i] = matrix.Di[i] * vector[i];

            for (int j = matrix.Ig[i]; j < matrix.Ig[i + 1]; j++)
            {
                product[i] += matrix.Ggl[j] * vector[matrix.Jg[j]];
                product[matrix.Jg[j]] += matrix.Ggu[j] * vector[i];
            }
        }

        return product;
    }

    /// <summary>
    /// Output the sparse matrix in dense format.
    /// </summary>
    /// <param name="path">File path.</param>
    public void SaveAsDense(string path)
    {
        double[,] a = new double[Size, Size];

        for (int i = 0; i < Size; i++)
        {
            a[i, i] = Di[i];

            for (int j = Ig[i]; j < Ig[i + 1]; j++)
            {
                a[i, Jg[j]] = Ggl[j];
                a[Jg[j], i] = Ggu[j];
            }
        }

        using var sw = new StreamWriter(path);
        for (int i = 0; i < Size; i++)
        {
            for (int j = 0; j < Size; j++)
            {
                sw.Write(a[i, j].ToString("0.0000") + "\t");
            }

            sw.WriteLine();
        }
    }

    /// <summary>
    /// Conversion of the sparse matrix to the profile format.
    /// Needed for direct solvers.
    /// </summary>
    public void ToProfileMatrix()
    {
        int[] ignew = Ig.ToArray();

        for (int i = 0; i < Size; i++)
        {
            int i0 = Ig[i];
            int i1 = Ig[i + 1];

            int profile = i1 - i0;

            if (profile > 0)
            {
                int count = i - Jg[i0];
                ignew[i + 1] = ignew[i] + count;
            }
            else
            {
                ignew[i + 1] = ignew[i];
            }
        }

        double[] gglnew = new double[ignew[^1]];
        double[] ggunew = new double[ignew[^1]];

        for (int i = 0; i < Size; i++)
        {
            int i0 = ignew[i];
            int i1 = ignew[i + 1];

            int j = i - (i1 - i0);

            int i0Old = Ig[i];

            for (int ik = i0; ik < i1; ik++, j++)
            {
                if (j == Jg[i0Old])
                {
                    gglnew[ik] = Ggl[i0Old];
                    ggunew[ik] = Ggu[i0Old];
                    i0Old++;
                }
                else
                {
                    gglnew[ik] = 0.0;
                    ggunew[ik] = 0.0;
                }
            }
        }

        Ig = ignew;
        Ggl = gglnew;
        Ggu = ggunew;
    }
}