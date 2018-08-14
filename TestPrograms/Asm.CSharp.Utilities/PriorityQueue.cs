using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;


namespace Utilities
{
    public class PriorityQueue<T> where T : IComparable<T>
    {
        private List<T> cells = new List<T>();

        public PriorityQueue()
        {
        }

        public PriorityQueue(IEnumerable<T> initialSet) : this()
        {
            this.cells.AddRange(initialSet);
            for (int i = this.cells.Count / 2; --i >= 0;)
            {
                this.sift(i);
            }
        }

        public T PeekMin()
        {
            int n = this.cells.Count;
            if (n <= 0)
            {
                throw new InvalidOperationException("Cannot peek minimum because the priority queue is empty.");
            }
            T result = this.cells[0];
            return result;
        }

        public T ExtractMin()
        {
            int n = this.cells.Count;
            if (n <= 0)
            {
                throw new InvalidOperationException("Cannot extract minimum because the priority queue is empty.");
            }
            T result = this.cells[0];
            this.cells[0] = this.cells[n - 1];
            this.cells.RemoveAt(n - 1);
            if (n > 1)
            {
                this.sift(0);
            }
            return result;
        }

        public void Add(T item)
        {
            this.cells.Add(item);
            this.climb(this.cells.Count - 1, 0);
        }

        private void sift(int inputIndex)
        {
            int son, n = this.cells.Count, index = inputIndex;
            T cache = this.cells[index];
            while ((son = index * 2 + 1) < n)
            {
                son = son + 1 < n && this.cells[son + 1].CompareTo(this.cells[son]) < 0 ? son + 1 : son;
                this.cells[index] = this.cells[son];
                index = son;
            }
            this.cells[index] = cache;
            climb(index, inputIndex);
        }

        private void climb(int index, int stopIndex)
        {
            T cache = this.cells[index];
            while (index > stopIndex)
            {
                int father = (index - 1) >> 1;
                if (cache.CompareTo(this.cells[father]) >= 0)
                {
                    break;
                }
                this.cells[index] = this.cells[father];
                index = father;
            }
            this.cells[index] = cache;
        }
    }

    public static class HeapUnitTest
    {
        public static bool UnitTest()
        {
            Random random = new Random(1001);
            List<double> initialList = Enumerable.Range(0, 100 * 1000).Select(idx => random.NextDouble() * random.NextDouble()).ToList();
            PriorityQueue<double> heap = new PriorityQueue<double>(initialList);
            for (int i = initialList.Count; --i >= 0;)
            {
                double value = random.NextDouble() * random.NextDouble();
                initialList.Add(value);
                heap.Add(value);
            }
            initialList.Sort();
            for (int i = 0; i < initialList.Count; i++)
            {
                if (initialList[i] != heap.ExtractMin())
                {
                    return false;
                }
            }
            return true;
        }
    }
}
