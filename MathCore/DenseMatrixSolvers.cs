namespace MathCore;

public abstract class DenseMatrixSolver
{
    protected Vector<double>? _solution;
    protected Vector<double> _vector = default!;
    protected Matrix<double> _matrix = default!;
    public ImmutableArray<double>? Solution => _solution?.ToImmutableArray();

    public void SetVector(Vector<double> vector)
        => _vector = Vector<double>.Copy(vector);

    public void SetMatrix(Matrix<double> matrix)
        => _matrix = Matrix<double>.Copy(matrix);

    protected DenseMatrixSolver(Matrix<double> matrix, Vector<double> vector)
        => (_matrix, _vector) = (matrix, vector);

    protected DenseMatrixSolver()
    {
    }

    public abstract void Compute();
}

public class Gauss : DenseMatrixSolver
{
    public Gauss(Matrix<double> matrix, Vector<double> vector) : base(matrix, vector)
    {
    }

    public Gauss()
    {
    }

    public override void Compute()
    {
        try
        {
            ArgumentNullException.ThrowIfNull(_matrix, $"{nameof(_matrix)} cannot be null, set the Matrix");
            ArgumentNullException.ThrowIfNull(_vector, $"{nameof(_vector)} cannot be null, set the Vector");

            if (_matrix.Rows != _matrix.Columns)
            {
                throw new NotSupportedException("The Gaussian method will not be able to solve this system");
            }

            const double eps = 1E-15;

            for (int k = 0; k < _matrix.Rows; k++)
            {
                var max = Math.Abs(_matrix[k, k]);
                int index = k;

                for (int i = k + 1; i < _matrix.Rows; i++)
                {
                    if (Math.Abs(_matrix[i, k]) > max)
                    {
                        max = Math.Abs(_matrix[i, k]);
                        index = i;
                    }
                }

                for (int j = 0; j < _matrix.Rows; j++)
                {
                    (_matrix[k, j], _matrix[index, j]) =
                        (_matrix[index, j], _matrix[k, j]);
                }

                (_vector[k], _vector[index]) = (_vector[index], _vector[k]);

                for (int i = k; i < _matrix.Rows; i++)
                {
                    double temp = _matrix[i, k];

                    if (Math.Abs(temp) < eps)
                    {
                        throw new Exception("Zero element of the column");
                    }

                    for (int j = 0; j < _matrix.Rows; j++)
                    {
                        _matrix[i, j] /= temp;
                    }

                    _vector[i] /= temp;

                    if (i != k)
                    {
                        for (int j = 0; j < _matrix.Rows; j++)
                        {
                            _matrix[i, j] -= _matrix[k, j];
                        }

                        _vector[i] -= _vector[k];
                    }
                }
            }

            _solution = new(_vector.Length);

            for (int k = _matrix.Rows - 1; k >= 0; k--)
            {
                _solution![k] = _vector[k];

                for (int i = 0; i < k; i++)
                {
                    _vector[i] -= _matrix[i, k] * _solution[k];
                }
            }
        }
        catch (Exception ex)
        {
            Console.WriteLine(ex.Message);
        }
    }
}