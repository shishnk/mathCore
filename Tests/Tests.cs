namespace Tests;

public class Testing
{
    private readonly ITestOutputHelper _output;

    public Testing(ITestOutputHelper output) => _output = output;

    [Fact]
    public void TestWithIterativeAndDirectSolvers()
    {
        SparseMatrix matrix = new(5, 8);
        Vector<double> vector = new(5); // 70, 85, 118, 125, 86

        matrix.Ig[0] = 0;
        matrix.Ig[1] = 0;
        matrix.Ig[2] = 1;
        matrix.Ig[3] = 2;
        matrix.Ig[4] = 4;
        matrix.Ig[5] = 8;

        matrix.Jg[0] = 0;
        matrix.Jg[1] = 1;
        matrix.Jg[2] = 0;
        matrix.Jg[3] = 2;
        matrix.Jg[4] = 0;
        matrix.Jg[5] = 1;
        matrix.Jg[6] = 2;
        matrix.Jg[7] = 3;

        matrix.Di[0] = 6;
        matrix.Di[1] = 7;
        matrix.Di[2] = 9;
        matrix.Di[3] = 10;
        matrix.Di[4] = -5;

        matrix.Ggl[0] = 2;
        matrix.Ggl[1] = 1;
        matrix.Ggl[2] = 4;
        matrix.Ggl[3] = 2;
        matrix.Ggl[4] = 7;
        matrix.Ggl[5] = 9;
        matrix.Ggl[6] = 10;
        matrix.Ggl[7] = 14;

        matrix.Ggu[0] = 2;
        matrix.Ggu[1] = 3;
        matrix.Ggu[2] = 5;
        matrix.Ggu[3] = 6;
        matrix.Ggu[4] = 8;
        matrix.Ggu[5] = 12;
        matrix.Ggu[6] = 13;
        matrix.Ggu[7] = 15;

        vector[0] = 70;
        vector[1] = 85;
        vector[2] = 118;
        vector[3] = 125;
        vector[4] = 86;

        // Массивы ig, jg, di, ggl, ggu матрицы должны быть заполнены значениями (код заполнения опущен).

        /* Матрица в плотном виде (возможно использовать метод PrintDense(string path) для проверки)
           6   2   0   5   8
           2   7   3   0   12
           0   1   9   6   13
           4   0   2   10  15
           0   1   9   6   13
        */

        // Создание итерационного решателя (возможны LOS, LOSLU, BCGSTABLU, CGM, CGMCholesky)
        IterativeSolver iterativeSolver = new BCGSTABLU(1000, 1E-14);

        // Создание прямого решателя
        DirectSolver directSolver = new DecomposerLDU();

        iterativeSolver.SetMatrix(matrix);
        // iterativeSolver.SetVector(vector);
        iterativeSolver.Compute();

        // Переводим матрицу в профильный формат
        matrix.ToProfileMatrix();

        directSolver.SetMatrix(matrix);
        directSolver.SetVector(vector);
        directSolver.Compute();

        _output.WriteLine(iterativeSolver.RunningTime.ToString());

        foreach (var (value, idx) in iterativeSolver.Solution!.Value.Select((value, idx) => (value, idx)))
        {
            _output.WriteLine($"b{idx} = {value}");
        }

        _output.WriteLine(directSolver.RunningTime.ToString());

        foreach (var (value, idx) in directSolver.Solution!.Value.Select((value, idx) => (value, idx)))
        {
            _output.WriteLine($"b{idx} = {value}");
        }

        /* Полученное решение (итерационный метод)
        b0 = 0.9999999999999978
        b1 = 2.000000000000005
        b2 = 2.999999999999996
        b3 = 4.0000000000000036
        b4 = 4.999999999999999
        */

        /* Полученное решение (прямой метод)
         b0 = 0.9999999999999996
         b1 = 2.000000000000001
         b2 = 2.999999999999999
         b3 = 4
         b4 = 5
         */

        /* Аналитическое решение
        b0 = 1
        b1 = 2
        b2 = 3
        b3 = 4
        b4 = 5
        */
    }

    [Fact]
    public void TestWithDenseMatrixSolvers()
    {
        SquareMatrix<double> matrix = new(4); // представляет собой двумерный массив
        Vector<double> vector = new(4); // 40, 30, 38, 68 

        /* Матрица в плотном виде
        1   4   1   2
        1   2   2   1
        -4   1   3   3
        1   1   5   4
        */

        // Создаем решатель СЛАУ для плотных матриц (в нашем случае для решения используется метод Гаусса)
        DenseMatrixSolver solver = new Gauss();

        vector[0] = 40;
        vector[1] = 30;
        vector[2] = 38;
        vector[3] = 68;

        matrix[0, 0] = 1;
        matrix[0, 1] = 4;
        matrix[0, 2] = 1;
        matrix[0, 3] = 2;

        matrix[1, 0] = 1;
        matrix[1, 1] = 2;
        matrix[1, 2] = 2;
        matrix[1, 3] = 1;

        matrix[2, 0] = -4;
        matrix[2, 1] = 1;
        matrix[2, 2] = 3;
        matrix[2, 3] = 3;

        matrix[3, 0] = 1;
        matrix[3, 1] = 1;
        matrix[3, 2] = 5;
        matrix[3, 3] = 4;

        solver.SetMatrix(matrix);
        solver.SetVector(vector);
        solver.Compute();

        foreach (var (value, idx) in solver.Solution!.Value.Select((value, idx) => (value, idx)))
        {
            _output.WriteLine($"b{idx} = {value}");
        }
        
        /* Полученное решение
         b0 = 2.0000000000000013
         b1 = 3.9999999999999987
         b2 = 6.000000000000002
         b3 = 8
         */

        /* Аналитическое решение
        b0 = 2
        b1 = 4
        b2 = 6
        b3 = 8
        */
    }
}