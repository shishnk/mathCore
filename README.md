# MathCore
.NET библиотека для численных методов.
 
## Модули библиотеки
* Классы, представляющий матрицы в плотном формате, разреженном формате
* Прямые решатели СЛАУ
* Итерационные решатели СЛАУ
* Класс вектора
 
## Примеры
 
### Использование прямых и итерационных решателей СЛАУ
```csharp
SparseMatrix matrix = new(5, 8);
Vector<double> vector = new(5); // 70, 85, 118, 125, 86
 
// Массивы ig, jg, di, ggl, ggu матрицы должны быть заполнены значениями (код заполнения опущен).
 
/* Матрица в плотном виде (возможно использовать метод PrintDense(string path) для проверки)
*   6   2   0   5   8
*   2   7   3   0   12
*   0   1   9   6   13
*   4   0   2   10  15
*   0   1   9   6   13
*/
 
// Создание итерационного решателя (возможны LOS, LOSLU, BCGSTABLU, CGM, CGMCholesky)
IterativeSolver iterativeSolver = new BCGSTABLU(1000, 1E-14);
 
// Создание прямого решателя
DirectSolver directSolver = new DecomposerLDU();
 
iterativeSolver.SetMatrix(matrix);
iterativeSolver.SetVector(vector);
iterativeSolver.Compute();
 
// Переводим матрицу в профильный формат
matrix.ToProfileMatrix();
 
directSolver.SetMatrix(matrix);
directSolver.SetVector(vector);
directSolver.Compute();
 
Console.WriteLine(iterativeSolver.RunningTime);
 
foreach (var (value, idx) in iterativeSolver.Solution!.Value.Select((value, idx) => (value, idx)))
{
    Console.WriteLine($"b{idx} = {value}");
}
 
Console.WriteLine(directSolver.RunningTime);
 
foreach (var (value, idx) in directSolver.Solution!.Value.Select((value, idx) => (value, idx)))
{
    Console.WriteLine($"b{idx} = {value}");
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
```
 
### Использование методов решения СЛАУ для плотной матрицы (example by Gross)
```csharp
SquareMatrix<double> matrix = new(4);
Vector<double> vector = new(4); // 40, 28, 38, 68
 
/* Матрица в плотном виде
   1   4   1   2
   1   2   2   1
  -4   1   3   3
   1   1   5   4
*/
 
// Создаем решатель СЛАУ для плотных матриц (в нашем случае для решения используется метод Гаусса)
DenseMatrixSolver solver = new Gauss();
 
solver.SetMatrix(matrix);
solver.SetVector(vector);
solver.Compute();
 
foreach (var (value, idx) in solver.Solution!.Value.Select((value, idx) => (value, idx)))
{
    Console.WriteLine($"b{idx} = {value}");
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
```
 
## NuGet
```
Install-Package Hukutka.MathCore
```
 