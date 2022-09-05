namespace MathCore;

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

    public Matrix(int size)
    {
        Rows = size;
        Columns = size;
        _storage = new T[size].Select(_ => new T[size]).ToArray();
    }

    public Matrix(int rows, int columns)
    {
        Rows = rows;
        Columns = columns;
        _storage = new T[rows].Select(_ => new T[columns]).ToArray();
    }

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
}

public class SparseMatrix
{
    // public fields - its bad, but the readability is better
    public int[] ig;
    public int[] jg;
    public double[] di;
    public double[] ggl;
    public double[] ggu;
    public int Size { get; }

    public SparseMatrix(int size, int sizeOffDiag)
    {
        Size = size;
        ig = new int[size + 1];
        jg = new int[sizeOffDiag];
        ggl = new double[sizeOffDiag];
        ggu = new double[sizeOffDiag];
        di = new double[size];
    }

    public static Vector<double> operator *(SparseMatrix matrix, Vector<double> vector)
    {
        Vector<double> product = new(vector.Length);

        for (int i = 0; i < vector.Length; i++)
        {
            product[i] = matrix.di[i] * vector[i];

            for (int j = matrix.ig[i]; j < matrix.ig[i + 1]; j++)
            {
                product[i] += matrix.ggl[j] * vector[matrix.jg[j]];
                product[matrix.jg[j]] += matrix.ggu[j] * vector[i];
            }
        }

        return product;
    }

    public void PrintDense(string path)
    {
        double[,] a = new double[Size, Size];

        for (int i = 0; i < Size; i++)
        {
            a[i, i] = di[i];

            for (int j = ig[i]; j < ig[i + 1]; j++)
            {
                a[i, jg[j]] = ggl[j];
                a[jg[j], i] = ggu[j];
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

    public void AsProfileMatrix()
    {
        int[] ignew = ig.ToArray();

        for (int i = 0; i < Size; i++)
        {
            int i0 = ig[i];
            int i1 = ig[i + 1];

            int profile = i1 - i0;

            if (profile > 0)
            {
                int count = i - jg[i0];
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

            int i0Old = ig[i];

            for (int ik = i0; ik < i1; ik++, j++)
            {
                if (j == jg[i0Old])
                {
                    gglnew[ik] = ggl[i0Old];
                    ggunew[ik] = ggu[i0Old];
                    i0Old++;
                }
                else
                {
                    gglnew[ik] = 0.0;
                    ggunew[ik] = 0.0;
                }
            }
        }

        ig = ignew;
        ggl = gglnew;
        ggu = ggunew;
    }
}