namespace MathCore;

public abstract class DirectSolver
{
    protected TimeSpan? _runningTime;
    protected SparseMatrix _matrix = default!;
    protected Vector<double> _vector = default!;
    protected Vector<double>? _solution;
    public TimeSpan? RunningTime => _runningTime;
    public ImmutableArray<double>? Solution => _solution?.ToImmutableArray();

    public void SetMatrix(SparseMatrix matrix)
        => _matrix = matrix;

    public void SetVector(Vector<double> vector)
        => _vector = vector;

    public abstract void Compute();

    protected void LDU()
    {
        for (int i = 0; i < _matrix.Size; i++)
        {
            double sumdi = 0.0;

            int i0 = _matrix.ig[i];
            int i1 = _matrix.ig[i + 1];

            int j = i - (i1 - i0);

            for (int ij = i0; ij < i1; ij++, j++)
            {
                double suml = 0.0;
                double sumu = 0.0;

                int j0 = _matrix.ig[j];
                int j1 = _matrix.ig[j + 1];

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
                    suml += _matrix.ggl[ik] * _matrix.di[k] * _matrix.ggu[jk];
                    sumu += _matrix.ggu[ik] * _matrix.di[k] * _matrix.ggl[jk];
                }

                _matrix.ggl[ij] = (_matrix.ggl[ij] - suml) / _matrix.di[j];
                _matrix.ggu[ij] = (_matrix.ggu[ij] - sumu) / _matrix.di[j];

                sumdi += _matrix.ggl[ij] * _matrix.ggu[ij] * _matrix.di[k];
            }

            _matrix.di[i] -= sumdi;
        }
    }

    protected void DiagonalStroke()
    {
        for (int i = 0; i < _matrix.Size; i++)
        {
            _solution![i] /= _matrix.di[i];
        }
    }

    protected void ForwardElimination()
    {
        for (int i = 0; i < _matrix.Size; i++)
        {
            double sum = 0.0;

            int i0 = _matrix.ig[i];
            int i1 = _matrix.ig[i + 1];

            int j = i - (i1 - i0);

            for (int ij = i0; ij < i1; ij++, j++)
            {
                sum += _matrix.ggl[ij] * _solution![j];
            }

            _solution![i] -= sum;
        }
    }

    protected void BackSubstitution()
    {
        for (int i = _matrix.Size - 1; i >= 0; i--)
        {
            int i0 = _matrix.ig[i];
            int i1 = _matrix.ig[i + 1];

            for (int j = i - 1, ij = i1 - 1; ij >= i0; ij--, j--)
            {
                _solution![j] -= _matrix.ggu[ij] * _solution[i];
            }
        }
    }
}

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