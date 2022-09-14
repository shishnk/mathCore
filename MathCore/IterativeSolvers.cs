namespace MathCore;

/// <summary>
/// Abstract base class <c>IterativeSolver</c> for iterative methods of solving SLAEs.
/// </summary>
public abstract class IterativeSolver
{
    protected TimeSpan? _runningTime;
    protected SparseMatrix _matrix = default!;
    protected Vector<double> _vector = default!;
    protected Vector<double>? _solution;

    /// <value>Property <c>MaxIters</c> to set maximum number of iterations the solver can make. </value>
    public int MaxIters { get; }

    /// <value>Property <c>Eps</c> for specifying the accuracy of the solution. </value>
    public double Eps { get; }

    /// <value>Property <c>RunningTime</c>to determine the time in which the SLAE was solved. </value>
    public TimeSpan? RunningTime => _runningTime;

    /// <value>Property <c>Solution</c>represents the solution of a given SLAE. </value>
    public ImmutableArray<double>? Solution => _solution?.ToImmutableArray();

    protected IterativeSolver(int maxIters, double eps)
        => (MaxIters, Eps) = (maxIters, eps);


    /// <summary>
    /// Set a sparse matrix for the solver.
    /// </summary>
    /// <param name="matrix">SparseMatrix</param>
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

    protected void Cholesky(double[] ggnew, double[] dinew)
    {
        double suml = 0.0;
        double sumdi = 0.0;

        for (int i = 0; i < _matrix.Size; i++)
        {
            int i0 = _matrix.Ig[i];
            int i1 = _matrix.Ig[i + 1];

            for (int k = i0; k < i1; k++)
            {
                int j = _matrix.Jg[k];
                int j0 = _matrix.Ig[j];
                int j1 = _matrix.Ig[j + 1];
                int ik = i0;
                int kj = j0;

                while (ik < k && kj < j1)
                {
                    if (_matrix.Jg[ik] == _matrix.Jg[kj])
                    {
                        suml += ggnew[ik] * ggnew[kj];
                        ik++;
                        kj++;
                    }
                    else
                    {
                        if (_matrix.Jg[ik] > _matrix.Jg[kj])
                            kj++;
                        else
                            ik++;
                    }
                }

                ggnew[k] = (ggnew[k] - suml) / dinew[j];
                sumdi += ggnew[k] * ggnew[k];
                suml = 0.0;
            }

            dinew[i] = Math.Sqrt(dinew[i] - sumdi);
            sumdi = 0.0;
        }
    }

    protected Vector<double> MoveForCholesky(Vector<double> vector, double[] ggnew, double[] dinew)
    {
        Vector<double> y = new(vector.Length);
        Vector<double> x = new(vector.Length);
        Vector<double>.Copy(vector, y);

        double sum = 0.0;

        for (int i = 0; i < _matrix.Size; i++) // Прямой ход
        {
            int i0 = _matrix.Ig[i];
            int i1 = _matrix.Ig[i + 1];

            for (int k = i0; k < i1; k++)
                sum += ggnew[k] * y[_matrix.Jg[k]];

            y[i] = (y[i] - sum) / dinew[i];
            sum = 0.0;
        }

        Vector<double>.Copy(y, x);

        for (int i = _matrix.Size - 1; i >= 0; i--) // Обратный ход
        {
            int i0 = _matrix.Ig[i];
            int i1 = _matrix.Ig[i + 1];
            x[i] = y[i] / dinew[i];

            for (int k = i0; k < i1; k++)
                y[_matrix.Jg[k]] -= ggnew[k] * x[i];
        }

        return x;
    }

    protected Vector<double> Direct(Vector<double> vector, double[] gglnew, double[] dinew)
    {
        Vector<double> y = new(vector.Length);
        Vector<double>.Copy(vector, y);

        double sum = 0.0;

        for (int i = 0; i < _matrix.Size; i++)
        {
            int i0 = _matrix.Ig[i];
            int i1 = _matrix.Ig[i + 1];

            for (int k = i0; k < i1; k++)
                sum += gglnew[k] * y[_matrix.Jg[k]];

            y[i] = (y[i] - sum) / dinew[i];
            sum = 0.0;
        }

        return y;
    }

    protected Vector<double> Reverse(Vector<double> vector, double[] ggunew)
    {
        Vector<double> result = new(vector.Length);
        Vector<double>.Copy(vector, result);

        for (int i = _matrix.Size - 1; i >= 0; i--)
        {
            int i0 = _matrix.Ig[i];
            int i1 = _matrix.Ig[i + 1];

            for (int k = i0; k < i1; k++)
                result[_matrix.Jg[k]] -= ggunew[k] * result[i];
        }

        return result;
    }

    protected void LU(double[] gglnew, double[] ggunew, double[] dinew)
    {
        double suml = 0.0;
        double sumu = 0.0;
        double sumdi = 0.0;

        for (int i = 0; i < _matrix.Size; i++)
        {
            int i0 = _matrix.Ig[i];
            int i1 = _matrix.Ig[i + 1];

            for (int k = i0; k < i1; k++)
            {
                int j = _matrix.Jg[k];
                int j0 = _matrix.Ig[j];
                int j1 = _matrix.Ig[j + 1];
                int ik = i0;
                int kj = j0;

                while (ik < k && kj < j1)
                {
                    if (_matrix.Jg[ik] == _matrix.Jg[kj])
                    {
                        suml += gglnew[ik] * ggunew[kj];
                        sumu += ggunew[ik] * gglnew[kj];
                        ik++;
                        kj++;
                    }
                    else if (_matrix.Jg[ik] > _matrix.Jg[kj])
                    {
                        kj++;
                    }
                    else
                    {
                        ik++;
                    }
                }

                gglnew[k] -= suml;
                ggunew[k] = (ggunew[k] - sumu) / dinew[j];
                sumdi += gglnew[k] * ggunew[k];
                suml = 0.0;
                sumu = 0.0;
            }

            dinew[i] -= sumdi;
            sumdi = 0.0;
        }
    }
}

/// <summary>
/// Local-optimal scheme
/// </summary>
public class LOS : IterativeSolver
{
    public LOS(int maxIters, double eps) : base(maxIters, eps)
    {
    }

    public override void Compute()
    {
        try
        {
            ArgumentNullException.ThrowIfNull(_matrix, $"{nameof(_matrix)} cannot be null, set the matrix");
            ArgumentNullException.ThrowIfNull(_vector, $"{nameof(_vector)} cannot be null, set the vector");

            _solution = new(_vector.Length);

            Vector<double> z = new(_vector.Length);

            Stopwatch sw = Stopwatch.StartNew();

            var r = _vector - (_matrix * _solution);

            Vector<double>.Copy(r, z);

            var p = _matrix * z;

            var squareNorm = r * r;

            for (int index = 0; index < MaxIters && squareNorm > Eps; index++)
            {
                var alpha = p * r / (p * p);
                _solution += alpha * z;
                squareNorm = (r * r) - (alpha * alpha * (p * p));
                r -= alpha * p;

                var tmp = _matrix * r;

                var beta = -(p * tmp) / (p * p);
                z = r + (beta * z);
                p = tmp + (beta * p);
            }

            sw.Stop();

            _runningTime = sw.Elapsed;
        }
        catch (Exception ex)
        {
            Console.WriteLine($"We had problem: {ex.Message}");
        }
    }
}

/// <summary>
/// Local-optimal scheme with LU decomposition preconditioning.
/// </summary>
public class LOSLU : IterativeSolver
{
    public LOSLU(int maxIters, double eps) : base(maxIters, eps)
    {
    }

    public override void Compute()
    {
        try
        {
            ArgumentNullException.ThrowIfNull(_matrix, $"{nameof(_matrix)} cannot be null, set the matrix");
            ArgumentNullException.ThrowIfNull(_vector, $"{nameof(_vector)} cannot be null, set the vector");

            _solution = new(_vector.Length);

            double[] gglnew = new double[_matrix.Ggl.Length];
            double[] ggunew = new double[_matrix.Ggu.Length];
            double[] dinew = new double[_matrix.Di.Length];

            _matrix.Ggl.Copy(gglnew);
            _matrix.Ggu.Copy(ggunew);
            _matrix.Di.Copy(dinew);

            Stopwatch sw = Stopwatch.StartNew();

            LU(gglnew, ggunew, dinew);

            var r = Direct(_vector - (_matrix * _solution), gglnew, dinew);
            var z = Reverse(r, ggunew);
            var p = Direct(_matrix * z, gglnew, dinew);

            var squareNorm = r * r;

            for (int iter = 0; iter < MaxIters && squareNorm > Eps; iter++)
            {
                var alpha = p * r / (p * p);
                squareNorm = (r * r) - (alpha * alpha * (p * p));
                _solution += alpha * z;
                r -= alpha * p;

                var tmp = Direct(_matrix * Reverse(r, ggunew), gglnew, dinew);

                var beta = -(p * tmp) / (p * p);
                z = Reverse(r, ggunew) + (beta * z);
                p = tmp + (beta * p);
            }

            sw.Stop();

            _runningTime = sw.Elapsed;
        }
        catch (Exception ex)
        {
            Console.WriteLine($"We had problem: {ex.Message}");
        }
    }
}

/// <summary>
/// Biconjugate gradient method with LU decomposition preconditioning.
/// </summary>
public class BCGSTABLU : IterativeSolver
{
    public BCGSTABLU(int maxIters, double eps) : base(maxIters, eps)
    {
    }

    public override void Compute()
    {
        try
        {
            ArgumentNullException.ThrowIfNull(_matrix, $"{nameof(_matrix)} cannot be null, set the matrix");
            ArgumentNullException.ThrowIfNull(_vector, $"{nameof(_vector)} cannot be null, set the vector");

            double alpha = 1.0;
            double omega = 1.0;
            double rho = 1.0;

            double vectorNorm = _vector.Norm();

            _solution = new(_vector.Length);

            double[] gglnew = new double[_matrix.Ggl.Length];
            double[] ggunew = new double[_matrix.Ggu.Length];
            double[] dinew = new double[_matrix.Di.Length];

            _matrix.Ggl.Copy(gglnew);
            _matrix.Ggu.Copy(ggunew);
            _matrix.Di.Copy(dinew);

            Vector<double> r0 = new(_vector.Length);
            Vector<double> p = new(_vector.Length);
            Vector<double> v = new(_vector.Length);

            Stopwatch sw = Stopwatch.StartNew();

            LU(gglnew, ggunew, dinew);

            var r = Direct(_vector - (_matrix * _solution), gglnew, dinew);

            Vector<double>.Copy(r, r0);

            for (int iter = 0; iter < MaxIters && r.Norm() / vectorNorm >= Eps; iter++)
            {
                var temp = rho;
                rho = r0 * r;
                var beta = rho / temp * (alpha / omega);
                p = r + (beta * (p - (omega * v)));
                v = Direct(_matrix * Reverse(p, ggunew), gglnew, dinew);
                alpha = rho / (r0 * v);
                var s = r - (alpha * v);
                var t = Direct(_matrix * Reverse(s, ggunew), gglnew, dinew);
                omega = t * s / (t * t);
                _solution += (omega * s) + (alpha * p);
                r = s - (omega * t);
            }

            _solution = Reverse(_solution, ggunew);

            sw.Stop();

            _runningTime = sw.Elapsed;
        }
        catch (Exception ex)
        {
            Console.WriteLine($"We had problem: {ex.Message}");
        }
    }

    /// <summary>
    /// Nonlinear conjugate gradient method.
    /// </summary>
    public class CGM : IterativeSolver
    {
        public CGM(int maxIters, double eps) : base(maxIters, eps)
        {
        }

        public override void Compute()
        {
            try
            {
                ArgumentNullException.ThrowIfNull(_matrix, $"{nameof(_matrix)} cannot be null, set the matrix");
                ArgumentNullException.ThrowIfNull(_vector, $"{nameof(_vector)} cannot be null, set the vector");

                double vectorNorm = _vector.Norm();

                _solution = new(_vector.Length);

                Vector<double> z = new(_vector.Length);

                Stopwatch sw = Stopwatch.StartNew();

                var r = _vector - (_matrix * _solution);

                Vector<double>.Copy(r, z);

                for (int iter = 0; iter < MaxIters && (r.Norm() / vectorNorm) >= Eps; iter++)
                {
                    var tmp = _matrix * z;
                    var alpha = r * r / (tmp * z);
                    _solution += alpha * z;
                    var squareNorm = r * r;
                    r -= alpha * tmp;
                    var beta = r * r / squareNorm;
                    z = r + (beta * z);
                }

                sw.Stop();

                _runningTime = sw.Elapsed;
            }
            catch (Exception ex)
            {
                Console.WriteLine($"We had problem: {ex.Message}");
            }
        }
    }

    /// <summary>
    /// Nonlinear conjugate gradient method with Cholesky(LL^T) decomposition preconditioning.
    /// </summary>
    public class CGMCholesky : IterativeSolver
    {
        public CGMCholesky(int maxIters, double eps) : base(maxIters, eps)
        {
        }

        public override void Compute()
        {
            try
            {
                ArgumentNullException.ThrowIfNull(_matrix, $"{nameof(_matrix)} cannot be null, set the matrix");
                ArgumentNullException.ThrowIfNull(_vector, $"{nameof(_vector)} cannot be null, set the vector");

                double vectorNorm = _vector.Norm();

                _solution = new(_vector.Length);

                double[] ggnew = new double[_matrix.Ggu.Length];
                double[] dinew = new double[_matrix.Di.Length];

                _matrix.Ggu.Copy(ggnew);
                _matrix.Di.Copy(dinew);

                Stopwatch sw = Stopwatch.StartNew();

                Cholesky(ggnew, dinew);

                var r = _vector - (_matrix * _solution);
                var z = MoveForCholesky(r, ggnew, dinew);

                for (int iter = 0; iter < MaxIters && r.Norm() / vectorNorm >= Eps; iter++)
                {
                    var tmp = MoveForCholesky(r, ggnew, dinew) * r;
                    var sndTemp = _matrix * z;
                    var alpha = tmp / (sndTemp * z);
                    _solution += alpha * z;
                    r -= alpha * sndTemp;
                    var fstTemp = MoveForCholesky(r, ggnew, dinew);
                    var beta = fstTemp * r / tmp;
                    z = fstTemp + (beta * z);
                }

                sw.Stop();

                _runningTime = sw.Elapsed;
            }
            catch (Exception ex)
            {
                Console.WriteLine($"We had problem: {ex.Message}");
            }
        }
    }
}