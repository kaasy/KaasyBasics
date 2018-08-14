using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace Utilities
{
    public static class BinarySearch
    {
        /// <summary>
        /// Performs a binary search on a sorted list.
        /// </summary>
        /// <typeparam name="T">Type of list items.</typeparam>
        /// <typeparam name="V">Type of the sought value.</typeparam>
        /// <param name="sortedList">Input sorted list to perform search on</param>
        /// <param name="value">The sought value.</param>
        /// <param name="convertor">Convertor function from list type to sought value type.</param>
        /// <param name="start">Start index of the search.</param>
        /// <param name="count">The number of elements to include in the search. Use count == -1 to search until 'list.Count - 1' inclusive.</param>
        /// <returns>The rightmost index in the list that contains a value lesser or equal than the sought value. 
        /// If count == -1 then the return index is in [start - 1.. list.Count - 1] 
        /// Else the return index is in [start - 1 .. start + count - 1] 
        /// </returns>
        public static int BinarySearchRightmostLesserOrEqual<T, V>(this IList<T> sortedList, V value,
            Func<T, V> convertor, int start = 0, int count = -1) where V : IComparable<V>
        {
            int bit, end;
            if (count == -1)
            {
                end = sortedList.Count;
                count = end - start;
            }
            else
            {
                end = start + count;
            }
            do
            {
                bit = count;
            } while ((count &= count - 1) != 0);
            int index = start - 1;
            for (; bit != 0; bit >>= 1)
            {
                int jump = index + bit;
                index = jump < end && convertor(sortedList[jump]).CompareTo(value) <= 0 ? jump : index;
            }
            return index;
        }

        /// <summary>
        /// Performs a binary search on a sorted list.
        /// </summary>
        /// <typeparam name="T">Type of list items.</typeparam>
        /// <typeparam name="V">Type of the sought value.</typeparam>
        /// <param name="sortedList">Input sorted list to perform search on</param>
        /// <param name="value">The sought value.</param>
        /// <param name="convertor">Convertor function from list type to sought value type.</param>
        /// <param name="start">Start index of the search.</param>
        /// <param name="count">The number of elements to include in the search. Use count == -1 to search until 'list.Count - 1' inclusive.</param>
        /// <returns>The rightmost index in the list that contains a value strictly lesser than the sought value. 
        /// If count == -1 then the return index is in [start - 1.. list.Count - 1] 
        /// Else the return index is in [start - 1 .. start + count - 1] 
        /// </returns>
        public static int BinarySearchRightmostLesser<T, V>(this IList<T> sortedList, V value,
            Func<T, V> convertor, int start = 0, int count = -1) where V : IComparable<V>
        {
            int bit, end;
            if (count == -1)
            {
                end = sortedList.Count;
                count = end - start;
            }
            else
            {
                end = start + count;
            }
            do
            {
                bit = count;
            } while ((count &= count - 1) != 0);
            int index = start - 1;
            for (; bit != 0; bit >>= 1)
            {
                int jump = index + bit;
                index = jump < end && convertor(sortedList[jump]).CompareTo(value) < 0 ? jump : index;
            }
            return index;
        }

        /// <summary>
        /// Performs a binary search on a sorted list.
        /// </summary>
        /// <typeparam name="T">Type of list items.</typeparam>
        /// <typeparam name="V">Type of the sought value.</typeparam>
        /// <param name="sortedList">Input sorted list to perform search on</param>
        /// <param name="value">The sought value.</param>
        /// <param name="convertor">Convertor function from list type to sought value type.</param>
        /// <param name="start">Start index of the search.</param>
        /// <param name="count">The number of elements to include in the search. Use count == -1 to search until 'list.Count - 1' inclusive.</param>
        /// <returns>The leftmost index in the list that contains a value greater or equal than the sought value. 
        /// If count == -1 then the return index is in [start .. list.Count] 
        /// Else the return index is in [start .. start + count] 
        /// </returns>
        public static int BinarySearchLeftmostGreaterOrEqual<T, V>(this IList<T> sortedList, V value,
            Func<T, V> convertor, int start = 0, int count = -1) where V : IComparable<V>
        {
            return 1 + sortedList.BinarySearchRightmostLesser(value, convertor, start, count);
        }

        /// <summary>
        /// Performs a binary search on a sorted list.
        /// </summary>
        /// <typeparam name="T">Type of list items.</typeparam>
        /// <typeparam name="V">Type of the sought value.</typeparam>
        /// <param name="sortedList">Input sorted list to perform search on</param>
        /// <param name="value">The sought value.</param>
        /// <param name="convertor">Convertor function from list type to sought value type.</param>
        /// <param name="start">Start index of the search.</param>
        /// <param name="count">The number of elements to include in the search. Use count == -1 to search until 'list.Count - 1' inclusive.</param>
        /// <returns>The leftmost index in the list that contains a value strictly greater than the sought value. 
        /// If count == -1 then the return index is in [start .. list.Count] 
        /// Else the return index is in [start .. start + count] 
        /// </returns>
        public static int BinarySearchLeftmostGreater<T, V>(this IList<T> sortedList, V value,
            Func<T, V> convertor, int start = 0, int count = -1) where V : IComparable<V>
        {
            return 1 + sortedList.BinarySearchRightmostLesserOrEqual(value, convertor, start, count);
        }

        /// <summary>
        /// Performs a binary search on a sorted list.
        /// </summary>
        /// <typeparam name="T">Type of list items.</typeparam>
        /// <param name="sortedList">Input sorted list to perform search on</param>
        /// <param name="value">The sought value.</param>
        /// <param name="convertor">Convertor function from list type to the 'double' value type.</param>
        /// <param name="start">Start index of the search.</param>
        /// <param name="count">The number of elements to include in the search. Use count == -1 to search until 'list.Count - 1' inclusive.</param>
        /// <returns>The index of the closest element in the sorted list to the sought value. 
        /// If count == -1 then the return index is in [start .. list.Count - 1] 
        /// Else the return index is in [start .. start + count - 1].
        /// If count == 0 then the returned value is 'start'.
        /// </returns>
        public static int BinarySearchLeftmostGreater<T>(this IList<T> sortedList, double value,
            Func<T, double> convertor, int start = 0, int count = -1)
        {
            int end = count == -1 ? sortedList.Count : start + count;
            int index = sortedList.BinarySearchRightmostLesserOrEqual(value, convertor, start, count);
            index = index < start || (index + 1 < end &&
                convertor(sortedList[index + 1]) - value < value - convertor(sortedList[index])) ? index + 1 : index;
            return index;
        }
    }

    public static class BinarySearchUnitTest
    {
        public static bool UnitTest()
        {
            Random random = new Random(10011);
            int n = 10 * 1000;
            List<int> list = new List<int>();
            for (int i = n; --i >= 0;)
            {
                list.Add(random.Next(i + 1));
            }
            list.Sort();
            for (int i = 5000; --i >= 0;)
            {
                int index1 = list.BinarySearchLeftmostGreaterOrEqual(i, x => x);
                int index2 = list.BinarySearchRightmostLesserOrEqual(i, x => x);
                int idx1 = list.IndexOf(i);
                int idx2 = list.LastIndexOf(i);
                if (idx1 == -1 || idx2 == -1)
                {
                    if (index2 - index1 != -1)
                    {
                        return false;
                    }
                    continue;
                }
                if (index1 != idx1 || index2 != idx2)
                {
                    return false;
                }
            }
            return true;
        }
    }
}
