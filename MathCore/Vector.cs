﻿namespace MathCore;

/// <summary>
/// Class <c>Vector</c> is a collection with overloads and methods.
/// <remarks>
/// T are the types that implement the INumber interface.
/// This list of types includes the unmanaged.
/// </remarks>
/// </summary>
public class Vector<T> : IEnumerable<T> where T : INumber<T>
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
        var result = T.Zero;

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
            result[i] = vector[i] * T.CreateChecked(constant);
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

    /// <summary>
    /// Method <c>Copy</c> copies the values to another vector object.
    /// </summary>
    /// <param name="source">A source vector.</param>
    /// <param name="destination">A destination vector.</param>
    public static void Copy(Vector<T> source, Vector<T> destination)
    {
        for (int i = 0; i < source.Length; i++)
        {
            destination[i] = source[i];
        }
    }

    /// <summary>
    /// Method <c>Copy</c> with a new vector object returned.
    /// </summary>
    /// <param name="otherVector">A vector from which the values are copied.</param>
    /// <returns>A vector with copied values.</returns>
    public static Vector<T> Copy(Vector<T> otherVector)
    {
        Vector<T> newVector = new(otherVector.Length);

        Array.Copy(otherVector._storage, newVector._storage, otherVector.Length);

        return newVector;
    }

    /// <summary>
    /// Method <c>Fill</c> for filling a vector with one value.
    /// </summary>
    /// <param name="value">Fill-in value.</param>
    public void Fill(double value)
    {
        for (int i = 0; i < Length; i++)
        {
            _storage[i] = T.CreateChecked(value);
        }
    }

    /// <summary>
    /// Method <c>Norm</c> that returns the Euclidean norm of a vector.
    /// </summary>
    /// <returns>Euclidean norm.</returns>
    public double Norm()
    {
        var result = T.Zero;

        for (int i = 0; i < Length; i++)
        {
            result += _storage[i] * _storage[i];
        }

        return Math.Sqrt(Convert.ToDouble(result));
    }

    public ImmutableArray<T> ToImmutableArray()
        => ImmutableArray.Create(_storage);

    public IEnumerator<T> GetEnumerator()
    {
        foreach (var value in _storage)
        {
            yield return value;
        }
    }

    IEnumerator IEnumerable.GetEnumerator() => GetEnumerator();

    /// <summary>
    ///  Method for adding values from a collection.
    /// </summary>
    /// <param name="collection"> A collection whose values will be added to the vector</param>
    /// <exception cref="ArgumentOutOfRangeException"></exception>
    public void Add(IEnumerable<T> collection)
    {
        var enumerable = collection as T[] ?? collection.ToArray();

        if (Length != enumerable.Length)
        {
            throw new ArgumentOutOfRangeException(nameof(collection), "Sizes of vector and collection not equal");
        }

        for (int i = 0; i < Length; i++)
        {
            _storage[i] = enumerable[i];
        }
    }
}