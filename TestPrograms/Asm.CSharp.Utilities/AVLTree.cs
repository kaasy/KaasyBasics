using System;
using System.Collections;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Text;

namespace Utilities
{
    public interface IAVLTreeNode<T>
    {
        IAVLTreeNode<T> LeftChild { get; }
        IAVLTreeNode<T> RightChild { get; }
        T Key { get; }
    }

    public class AVLTree<T> : IList<T>
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

        private static TreeNode insertAt(TreeNode node, int index, T key)
        {
            if (node == null)
            {
                return new TreeNode
                {
                    Height = 1,
                    Key = key
                };
            }
            if (index <= node.LeftCount)
            {
                node.Left = insertAt(node.Left, index, key);
                node.LeftCount++;
            }
            else
            {
                node.Right = insertAt(node.Right, index - node.LeftCount - 1, key);
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
            set { getNodeAt(this.root, index).Key = value; }
        }

        public AVLTree()
        {
        }

        public AVLTree(IEnumerable<T> initialSet) : this()
        {
            foreach (T item in initialSet)
            {
                this.Add(item);
            }
        }

        public int IndexOf(T item)
        {
            throw new NotSupportedException();
        }

        public void Insert(int index, T item)
        {
            this.root = insertAt(this.root, index, item);
            this.Count++;
        }

        public void RemoveAt(int index)
        {
            this.root = removeAt(this.root, index);
            this.Count--;
        }

        public void Add(T item)
        {
            this.Insert(this.Count, item);
        }

        public void Clear()
        {
            clear(this.root);
            this.root = null;
            this.Count = 0;
        }

        public bool Contains(T item)
        {
            throw new NotSupportedException();
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
            throw new NotSupportedException();
        }

        public IEnumerator<T> GetEnumerator()
        {
            return new AVLTreeEnumerator<T>(this.root);
        }

        IEnumerator IEnumerable.GetEnumerator()
        {
            return this.GetEnumerator();
        }

        private static TimeSpan getInsertTime(IList<int> list, int count, Random random)
        {
            Stopwatch stopwatch = new Stopwatch();
            stopwatch.Restart();
            for (int i = count; --i >= 0;)
            {
                int index = random.Next(list.Count + 1);
                int value = random.Next();
                list.Insert(index, value);
            }
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
                list.RemoveAt(index);
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

        private static bool compareAndVerify(IList<int> list, AVLTree<int> tree)
        {
            if (!areEqual(list, tree))
            {
                return false;
            }
            int treeNodesCount;
            bool ok = AVLTree<int>.verify(tree.root, out treeNodesCount);
            if (ok && treeNodesCount == tree.Count)
            {
                return true;
            }
            return false;
        }

        public static bool UnitTest()
        {
            int n = 10 * 1000;
            int randomSeed = 123455;
            Random randomList = new Random(randomSeed);
            Random randomTree = new Random(randomSeed);

            List<int> list = new List<int>();
            AVLTree<int> tree = new AVLTree<int>();
            bool ok;

            TimeSpan timeInsertList = getInsertTime(list, n, randomList);
            TimeSpan timeInsertTree = getInsertTime(tree, n, randomTree);
            ok = compareAndVerify(list, tree);
            if (!ok)
            {
                return false;
            }
            TimeSpan timeDeleteList = getDeleteTime(list, n / 2, randomList);
            TimeSpan timeDeleteTree = getDeleteTime(tree, n / 2, randomTree);
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

    public class AVLTreeEnumerator<T> : IEnumerator<T>
    {
        private IAVLTreeNode<T> rootNode, pendingNode;
        private Stack<IAVLTreeNode<T>> stack = new Stack<IAVLTreeNode<T>>();

        public AVLTreeEnumerator(IAVLTreeNode<T> node)
        {
            this.rootNode = node;
            this.Reset();
        }

        public T Current { get { return this.stack.Count > 0 ? this.stack.Peek().Key : default(T); } }

        object IEnumerator.Current { get { return this.Current; } }

        public void Dispose()
        {
            this.stack = null;
            this.rootNode = this.pendingNode = null;
        }

        private void addToStack(IAVLTreeNode<T> node)
        {
            while (node != null)
            {
                this.stack.Push(node);
                node = node.LeftChild;
            }
        }

        public bool MoveNext()
        {
            if (this.pendingNode != null)
            {
                this.addToStack(this.pendingNode);
                this.pendingNode = null;
                return this.stack.Count > 0;
            }
            if (this.stack.Count > 0)
            {
                IAVLTreeNode<T> node = this.stack.Pop();
                this.addToStack(node.RightChild);
            }
            return this.stack.Count > 0;
        }

        public void Reset()
        {
            this.pendingNode = this.rootNode;
            this.stack.Clear();
        }
    }
}
