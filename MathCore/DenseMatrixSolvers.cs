namespace MathCore;

/// <summary>
/// Abstract base class <c>DenseMatrixSolver</c> to solve the SLAE for a dense matrix.
/// </summary>
public abstract class DenseMatrixSolver
{
    protected Vector<double>? _solution;
    protected Vector<double> _vector = default!;
    protected SquareMatrix<double> _matrix = default!;

    /// <summary>
    /// Property <c>Solution</c>represents the solution of a given SLAE.
    /// </summary>
    public ImmutableArray<double>? Solution => _solution?.ToImmutableArray();

    /// <summary>
    /// Set a dense matrix for the solver.
    /// </summary>
    /// <param name="matrix">Dense matrix.</param>
    public void SetMatrix(SquareMatrix<double> matrix)
        => _matrix = SquareMatrix<double>.Copy(matrix);

    /// <summary>
    /// Set a vector for the solver.
    /// </summary>
    /// <param name="vector">Vector of the right part.</param>
    public void SetVector(Vector<double> vector)
        => _vector = Vector<double>.Copy(vector);

    protected DenseMatrixSolver(SquareMatrix<double> matrix, Vector<double> vector)
        => (_matrix, _vector) = (SquareMatrix<double>.Copy(matrix), Vector<double>.Copy(vector));

    protected DenseMatrixSolver()
    {
    }

    /// <summary>
    /// Start solving the SLAE.
    /// </summary>
    public abstract void Compute();
}

/// <summary>
/// Class <c>Gauss</c> to solve the SLAE by the Gaussian method.
/// </summary>
public class Gauss : DenseMatrixSolver
{
    public Gauss(SquareMatrix<double> matrix, Vector<double> vector) : base(matrix, vector)
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

            const double eps = 1E-15;

            for (int k = 0; k < _matrix.Size; k++)
            {
                var max = Math.Abs(_matrix[k, k]);
                int index = k;

                for (int i = k + 1; i < _matrix.Size; i++)
                {
                    if (Math.Abs(_matrix[i, k]) > max)
                    {
                        max = Math.Abs(_matrix[i, k]);
                        index = i;
                    }
                }

                for (int j = 0; j < _matrix.Size; j++)
                {
                    (_matrix[k, j], _matrix[index, j]) =
                        (_matrix[index, j], _matrix[k, j]);
                }

                (_vector[k], _vector[index]) = (_vector[index], _vector[k]);

                for (int i = k; i < _matrix.Size; i++)
                {
                    double temp = _matrix[i, k];

                    if (Math.Abs(temp) < eps)
                    {
                        throw new Exception("Zero element of the column");
                    }

                    for (int j = 0; j < _matrix.Size; j++)
                    {
                        _matrix[i, j] /= temp;
                    }

                    _vector[i] /= temp;

                    if (i != k)
                    {
                        for (int j = 0; j < _matrix.Size; j++)
                        {
                            _matrix[i, j] -= _matrix[k, j];
                        }

                        _vector[i] -= _vector[k];
                    }
                }
            }

            _solution = new(_vector.Length);

            for (int k = _matrix.Size - 1; k >= 0; k--)
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