namespace MathCore;

/// <summary>
/// Abstract base class <c>DirectSolver</c> for direct methods of solving SLAEs.
/// </summary>
public abstract class DirectSolver
{
    protected TimeSpan? _runningTime;
    protected SparseMatrix _matrix = default!;
    protected Vector<double> _vector = default!;
    protected Vector<double>? _solution;

    /// <value>Property <c>RunningTime</c>to determine the time in which the SLAE was solved. </value>
    public TimeSpan? RunningTime => _runningTime;

    /// <value>Property <c>Solution</c>represents the solution of a given SLAE. </value>
    public ImmutableArray<double>? Solution => _solution?.ToImmutableArray();

    /// <summary>
    /// Set a sparse matrix for the solver.
    /// </summary>
    /// <param name="matrix">Sparse matrix converted to profile format.</param>
    public void SetMatrix(SparseMatrix matrix)
        => _matrix = matrix;

    /// <summary>
    /// Set a vector for the solver.
    /// </summary>
    /// <param name="vector">Vector of the right part.</param>
    public void SetVector(Vector<double> vector)
        => _vector = vector;

    /// <summary>
    /// Start solving the SLAE.
    /// </summary>
    public abstract void Compute();

    protected void LDU()
    {
        for (int i = 0; i < _matrix.Size; i++)
        {
            double sumdi = 0.0;

            int i0 = _matrix.Ig[i];
            int i1 = _matrix.Ig[i + 1];

            int j = i - (i1 - i0);

            for (int ij = i0; ij < i1; ij++, j++)
            {
                double suml = 0.0;
                double sumu = 0.0;

                int j0 = _matrix.Ig[j];
                int j1 = _matrix.Ig[j + 1];

                int ik = i0;
                int jk = j0;

                int k = i - (i1 - i0);

                int ci = ij - i0;
                int cj = j1 - j0;

                int cij = ci - cj;

                if (cij > 0)
                {
                    ik += cij;
                    k += cij;
                }
                else
                {
                    jk -= cij;
                }

                for (; ik < ij; ik++, jk++, k++)
                {
                    suml += _matrix.Ggl[ik] * _matrix.Di[k] * _matrix.Ggu[jk];
                    sumu += _matrix.Ggu[ik] * _matrix.Di[k] * _matrix.Ggl[jk];
                }

                _matrix.Ggl[ij] = (_matrix.Ggl[ij] - suml) / _matrix.Di[j];
                _matrix.Ggu[ij] = (_matrix.Ggu[ij] - sumu) / _matrix.Di[j];

                sumdi += _matrix.Ggl[ij] * _matrix.Ggu[ij] * _matrix.Di[k];
            }

            _matrix.Di[i] -= sumdi;
        }
    }

    protected void DiagonalStroke()
    {
        for (int i = 0; i < _matrix.Size; i++)
        {
            _solution![i] /= _matrix.Di[i];
        }
    }

    protected void ForwardElimination()
    {
        for (int i = 0; i < _matrix.Size; i++)
        {
            double sum = 0.0;

            int i0 = _matrix.Ig[i];
            int i1 = _matrix.Ig[i + 1];

            int j = i - (i1 - i0);

            for (int ij = i0; ij < i1; ij++, j++)
            {
                sum += _matrix.Ggl[ij] * _solution![j];
            }

            _solution![i] -= sum;
        }
    }

    protected void BackSubstitution()
    {
        for (int i = _matrix.Size - 1; i >= 0; i--)
        {
            int i0 = _matrix.Ig[i];
            int i1 = _matrix.Ig[i + 1];

            for (int j = i - 1, ij = i1 - 1; ij >= i0; ij--, j--)
            {
                _solution![j] -= _matrix.Ggu[ij] * _solution[i];
            }
        }
    }
}

/// <summary>
/// Class <c>DecomposerLDU</c> to solve the SLAE using the LDU decomposition..
/// </summary>
public class DecomposerLDU : DirectSolver
{
    public override void Compute()
    {
        _solution = new(_matrix.Size);
        Vector<double>.Copy(_vector, _solution);

        Stopwatch sw = Stopwatch.StartNew();

        LDU();
        ForwardElimination();
        DiagonalStroke();
        BackSubstitution();

        sw.Stop();

        _runningTime = sw.Elapsed;
    }
}