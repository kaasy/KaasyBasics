using System;
using System.Collections;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Text;

namespace Utilities
{
    public class AVLTreeSorted<T> : IList<T> where T : IComparable<T>
    {
        class TreeNode : IAVLTreeNode<T>
        {
            public TreeNode Left { get; set; }
            public TreeNode Right { get; set; }
            public T Key { get; set; }
            public int Height { get; set; }
            public int LeftCount { get; set; }

            public IAVLTreeNode<T> LeftChild { get { return this.Left; } }
            public IAVLTreeNode<T> RightChild { get { return this.Right; } }

            public void RecomputeHeight()
            {
                this.Height = 1 + Math.Max(this.Left == null ? 0 : this.Left.Height, this.Right == null ? 0 : this.Right.Height);
            }
        }

        #region Tree operations

        private static TreeNode rotateLeft(TreeNode node)
        {
            TreeNode right = node.Right;
            node.Right = right.Left;
            right.Left = node;
            node.RecomputeHeight();
            right.RecomputeHeight();
            right.LeftCount += node.LeftCount + 1;
            return right;
        }

        private static TreeNode rotateRight(TreeNode node)
        {
            TreeNode left = node.Left;
            node.Left = left.Right;
            left.Right = node;
            node.RecomputeHeight();
            left.RecomputeHeight();
            node.LeftCount -= left.LeftCount + 1;
            return left;
        }

        private static int getHeight(TreeNode node)
        {
            return node == null ? 0 : node.Height;
        }

        private static TreeNode balance(TreeNode node)
        {
            TreeNode left = node.Left;
            TreeNode right = node.Right;
            int leftHeight = getHeight(left);
            int rightHeight = getHeight(right);
            node.Height = 1 + Math.Max(leftHeight, rightHeight);
            if (leftHeight - rightHeight > 1)
            {
                int leftLeftHeight = getHeight(left.Left);
                int leftRightHeight = getHeight(left.Right);
                if (leftRightHeight > leftLeftHeight)
                {
                    node.Left = rotateLeft(left);
                }
                node = rotateRight(node);
            }
            else if (rightHeight - leftHeight > 1)
            {
                int rightRightHeight = getHeight(right.Right);
                int rightLeftHeight = getHeight(right.Left);
                if (rightLeftHeight > rightRightHeight)
                {
                    node.Right = rotateRight(right);
                }
                node = rotateLeft(node);
            }
            //{   //verification.
            //    leftHeight = getHeight(node.Left);
            //    rightHeight = getHeight(node.Right);
            //    if (Math.Abs(leftHeight - rightHeight) > 1)
            //    {
            //        throw new InvalidOperationException();
            //    }
            //}
            return node;
        }

        private static TreeNode insertKey(TreeNode node, T key)
        {
            if (node == null)
            {
                return new TreeNode
                {
                    Height = 1,
                    Key = key
                };
            }
            if (key.CompareTo(node.Key) < 0)
            {
                node.Left = insertKey(node.Left, key);
                node.LeftCount++;
            }
            else
            {   //case key >= node.Key (insert last).
                node.Right = insertKey(node.Right, key);
            }
            node = balance(node);
            return node;
        }

        private static TreeNode getLeftmostNode(TreeNode node)
        {
            while (node.Left != null)
            {
                node = node.Left;
            }
            return node;
        }

        private static TreeNode removeAt(TreeNode node, int index)
        {
            if (node == null)
            {
                throw new IndexOutOfRangeException();
            }
            if (index < node.LeftCount)
            {
                node.Left = removeAt(node.Left, index);
                node.LeftCount--;
            }
            else if (index > node.LeftCount)
            {
                node.Right = removeAt(node.Right, index - node.LeftCount - 1);
            }
            else
            {
                if (node.Left == null || node.Right == null)
                {
                    node.Key = default(T);
                    return node.Left != null ? node.Left : node.Right;
                }
                TreeNode next = getLeftmostNode(node.Right);
                node.Key = next.Key;
                //next.Key = default(T);    //not necessary because it is set by removeAt at a deeper recursion.
                node.Right = removeAt(node.Right, 0);
            }
            node = balance(node);
            return node;
        }

        private static TreeNode getNodeAt(TreeNode node, int index)
        {
            if (node == null)
            {
                throw new IndexOutOfRangeException();
            }
            TreeNode result;
            if (index < node.LeftCount)
            {
                result = getNodeAt(node.Left, index);
            }
            else if (index == node.LeftCount)
            {
                result = node;
            }
            else
            {
                result = getNodeAt(node.Right, index - node.LeftCount - 1);
            }
            return result;
        }

        private static int getLeftmostGreaterOrEqualIndex(TreeNode node, T key, bool onlyExactMatch, out TreeNode foundNode)
        {
            if (node == null)
            {
                foundNode = null;
                return 0;
            }
            int index;
            int comparison = key.CompareTo(node.Key);
            if (comparison <= 0)
            {   //delete first.
                index = getLeftmostGreaterOrEqualIndex(node.Left, key, onlyExactMatch, out foundNode);
                foundNode = foundNode == null && (!onlyExactMatch || comparison == 0) ? node : foundNode;
            }
            else
            {
                index = node.LeftCount + 1 + getLeftmostGreaterOrEqualIndex(node.Right, key, onlyExactMatch, out foundNode);
            }
            return index;
        }

        private static TreeNode removeFirstOccurence(TreeNode node, T key, out bool wasRemoved)
        {
            if (node == null)
            {
                wasRemoved = false;
                return null;
            }
            int comparison = key.CompareTo(node.Key);
            if (comparison <= 0)
            {   //delete first.
                node.Left = removeFirstOccurence(node.Left, key, out wasRemoved);
                if (wasRemoved)
                {
                    node.LeftCount--;
                }
                else if (comparison == 0)
                {
                    wasRemoved = true;
                    if (node.Left == null || node.Right == null)
                    {
                        node.Key = default(T);
                        return node.Left != null ? node.Left : node.Right;
                    }
                    TreeNode next = getLeftmostNode(node.Right);
                    node.Key = next.Key;
                    //next.Key = default(T);    //not necessary because it is set by removeAt.
                    node.Right = removeAt(node.Right, 0);
                    return balance(node);
                }
            }
            else
            {
                node.Right = removeFirstOccurence(node.Right, key, out wasRemoved);
            }
            return wasRemoved ? balance(node) : node;
        }

        private static void clear(TreeNode node)
        {
            if (node == null)
            {
                return;
            }
            node.Key = default(T);
            clear(node.Left);
            clear(node.Right);
            node.Left = node.Right = null;
        }

        private static bool verify(TreeNode node, out int count)
        {
            if (node == null)
            {
                count = 0;
                return true;
            }
            int leftCount, rightCount;
            bool okLeft = verify(node.Left, out leftCount);
            bool okRight = verify(node.Right, out rightCount);
            count = leftCount + 1 + rightCount;
            int leftHeight = getHeight(node.Left);
            int rightHeight = getHeight(node.Right);
            if (okLeft && okRight && leftCount == node.LeftCount && Math.Abs(leftHeight - rightHeight) <= 1)
            {
                return true;
            }
            return false;
        }

        #endregion

        private TreeNode root;

        public int Count { get; private set; }

        public bool IsReadOnly { get { return false; } }

        public T this[int index]
        {
            get { return getNodeAt(this.root, index).Key; }
            set { throw new NotSupportedException(); }
        }

        public AVLTreeSorted()
        {
        }

        public AVLTreeSorted(IEnumerable<T> initialSet) : this()
        {
            foreach (T item in initialSet)
            {
                this.Add(item);
            }
        }

        public int IndexOf(T item)
        {
            TreeNode foundNode;
            int index = getLeftmostGreaterOrEqualIndex(this.root, item, true, out foundNode);
            if (foundNode != null)
            {
                return index;
            }
            return -1;
        }

        public int GetLeftmostGreaterOrEqualIndex(T item, out T foundItem)
        {
            TreeNode foundNode;
            int index = getLeftmostGreaterOrEqualIndex(this.root, item, false, out foundNode);
            foundItem = foundNode == null ? default(T) : foundNode.Key;
            return index;
        }

        public void Insert(int index, T item)
        {
            throw new NotSupportedException();
        }

        public void RemoveAt(int index)
        {
            this.root = removeAt(this.root, index);
            this.Count--;
        }

        public void Add(T item)
        {
            this.root = insertKey(this.root, item);
            this.Count++;
        }

        public void Clear()
        {
            clear(this.root);
            this.root = null;
            this.Count = 0;
        }

        public bool Contains(T item)
        {
            return this.IndexOf(item) >= 0;
        }

        private static void copyTo(TreeNode node, T[] array, ref int arrayIndex)
        {
            if (node == null)
            {
                return;
            }
            copyTo(node.Left, array, ref arrayIndex);
            array[arrayIndex++] = node.Key;
            copyTo(node.Right, array, ref arrayIndex);
        }

        public void CopyTo(T[] array, int arrayIndex)
        {
            copyTo(this.root, array, ref arrayIndex);
        }

        public bool Remove(T item)
        {
            bool wasRemoved;
            this.root = removeFirstOccurence(this.root, item, out wasRemoved);
            if (wasRemoved)
            {
                this.Count--;
                return true;
            }
            return false;
            //slower - double iteration code below :
            //int index = this.IndexOf(item);
            //if (index >= 0)
            //{
            //    this.RemoveAt(index);
            //    return true;
            //}
            //return false;
        }

        public IEnumerator<T> GetEnumerator()
        {
            return new AVLTreeEnumerator<T>(this.root);
        }

        IEnumerator IEnumerable.GetEnumerator()
        {
            return this.GetEnumerator();
        }

        private static TimeSpan getAddTime(IList<int> list, int count, Random random, Action postAction)
        {
            Stopwatch stopwatch = new Stopwatch();
            stopwatch.Restart();
            for (int i = count; --i >= 0;)
            {
                int value = random.Next();
                list.Add(value);
            }
            postAction();
            stopwatch.Stop();
            return stopwatch.Elapsed;
        }

        private static TimeSpan getDeleteTime(IList<int> list, int count, Random random)
        {
            Stopwatch stopwatch = new Stopwatch();
            stopwatch.Restart();
            for (int i = count; --i >= 0;)
            {
                int index = random.Next(list.Count);
                list.Remove(list[index]);
            }
            stopwatch.Stop();
            return stopwatch.Elapsed;
        }

        private static bool areEqual(IList<int> list1, IList<int> list2)
        {
            if (list1.Count != list2.Count)
            {
                return false;
            }
            for (int i = list1.Count; --i >= 0;)
            {
                if (list1[i] != list2[i])
                {
                    return false;
                }
            }
            return true;
        }

        private static bool compareAndVerify(IList<int> list, AVLTreeSorted<int> tree)
        {
            if (!areEqual(list, tree))
            {
                return false;
            }
            int treeNodesCount;
            bool ok = AVLTreeSorted<int>.verify(tree.root, out treeNodesCount);
            if (ok && treeNodesCount == tree.Count)
            {
                return true;
            }
            return false;
        }

        private static TimeSpan getIndexOfTime(IList<int> list, int count, Random random)
        {
            Stopwatch stopwatch = new Stopwatch();
            stopwatch.Restart();

            int[] startCount = new int[count + 1];
            int total = 0;
            for (int i = 0; i < count; i++)
            {
                int times = random.Next(3);
                for (int k = times; --k >= 0;)
                {
                    list.Add(i);
                }
                startCount[i] = total;
                total += times;
            }
            startCount[count] = total;
            for (int i = 0; i < count; i++)
            {
                int index = list.IndexOf(i);
                if ((index >= 0 && index != startCount[i]) || (index == -1 && startCount[i] != startCount[i + 1]))
                {
                    throw new InvalidOperationException();
                }
            }

            stopwatch.Stop();
            return stopwatch.Elapsed;
        }

        public static bool UnitTest()
        {
            int n = 10 * 1000;
            int randomSeed = 123455;
            Random randomList = new Random(randomSeed);
            Random randomTree = new Random(randomSeed);

            List<int> list = new List<int>();
            AVLTreeSorted<int> tree = new AVLTreeSorted<int>();
            bool ok;

            TimeSpan timeIndexOfList = getIndexOfTime(list, n, randomList);
            TimeSpan timeIndexOfTree = getIndexOfTime(tree, n, randomTree);
            ok = compareAndVerify(list, tree);
            if (!ok)
            {
                return false;
            }

            TimeSpan timeDeleteList = getDeleteTime(list, list.Count, randomList);
            TimeSpan timeDeleteTree = getDeleteTime(tree, tree.Count, randomTree);
            ok = compareAndVerify(list, tree);
            if (!ok)
            {
                return false;
            }

            TimeSpan timeInsertList = getAddTime(list, n, randomList, () => list.Sort());
            TimeSpan timeInsertTree = getAddTime(tree, n, randomTree, () => { });
            ok = compareAndVerify(list, tree);
            if (!ok)
            {
                return false;
            }

            timeDeleteList += getDeleteTime(list, n / 2, randomList);
            timeDeleteTree += getDeleteTime(tree, n / 2, randomTree);
            ok = compareAndVerify(list, tree);
            if (!ok)
            {
                return false;
            }

            int[] array = new int[list.Count + 100];
            tree.CopyTo(array, 100);
            ok = areEqual(array.Skip(100).ToList(), list);
            if (!ok)
            {
                return false;
            }

            List<int> enumeratorList = new List<int>();
            foreach (int item in tree)
            {
                enumeratorList.Add(item);
            }
            ok = areEqual(enumeratorList, list);
            if (!ok)
            {
                return false;
            }

            timeDeleteList += getDeleteTime(list, n / 2, randomList);
            timeDeleteTree += getDeleteTime(tree, n / 2, randomTree);
            ok = compareAndVerify(list, tree);
            if (!ok)
            {
                return false;
            }

            return true;
        }
    }
}
