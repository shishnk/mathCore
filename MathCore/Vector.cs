namespace MathCore;

public class Vector<T> where T : INumber<T>
{
    private readonly T[] _storage;
    public int Length { get; }

    public T this[int idx]
    {
        get => _storage[idx];
        set => _storage[idx] = value;
    }

    public Vector(int length)
        => (Length, _storage) = (length, new T[length]);

    public static T operator *(Vector<T> a, Vector<T> b)
    {
        T result = T.Zero;

        for (int i = 0; i < a.Length; i++)
        {
            result += a[i] * b[i];
        }

        return result;
    }

    public static Vector<T> operator *(double constant, Vector<T> vector)
    {
        Vector<T> result = new(vector.Length);

        for (int i = 0; i < vector.Length; i++)
        {
            result[i] = vector[i] * T.Create(constant);
        }

        return result;
    }

    public static Vector<T> operator +(Vector<T> a, Vector<T> b)
    {
        Vector<T> result = new(a.Length);

        for (int i = 0; i < a.Length; i++)
        {
            result[i] = a[i] + b[i];
        }

        return result;
    }

    public static Vector<T> operator -(Vector<T> a, Vector<T> b)
    {
        Vector<T> result = new(a.Length);

        for (int i = 0; i < a.Length; i++)
        {
            result[i] = a[i] - b[i];
        }

        return result;
    }

    public static void Copy(Vector<T> source, Vector<T> destination)
    {
        for (int i = 0; i < source.Length; i++)
        {
            destination[i] = source[i];
        }
    }

    public static Vector<T> Copy(Vector<T> otherVector)
    {
        Vector<T> newVector = new(otherVector.Length);

        Array.Copy(otherVector._storage, newVector._storage, otherVector.Length);

        return newVector;
    }

    public void Fill(double value)
    {
        for (int i = 0; i < Length; i++)
        {
            _storage[i] = T.Create(value);
        }
    }

    public double Norm()
    {
        T result = T.Zero;

        for (int i = 0; i < Length; i++)
        {
            result += _storage[i] * _storage[i];
        }

        return Math.Sqrt(Convert.ToDouble(result));
    }

    public ImmutableArray<T> ToImmutableArray()
        => ImmutableArray.Create(_storage);
}