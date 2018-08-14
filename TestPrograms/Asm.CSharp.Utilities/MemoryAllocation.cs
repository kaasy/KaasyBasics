using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Runtime.InteropServices;
using System.Text;
using System.Threading.Tasks;

namespace Utilities
{
    public struct MemoryBlockBySize : IComparable<MemoryBlockBySize>
    {
        public long Address;
        public long Size;

        public long End { get { return this.Address + this.Size; } }

        public MemoryBlockBySize(long address, long size) : this()
        {
            this.Address = address;
            this.Size = size;
        }

        public MemoryBlockByAddress ToAddress { get { return new MemoryBlockByAddress(this.Address, this.Size); } }

        public int CompareTo(MemoryBlockBySize other)
        {
            int compare = this.Size.CompareTo(other.Size);
            if (compare != 0)
            {
                return compare;
            }
            return this.Address.CompareTo(other.Address);
        }

        public long GetOverlap(MemoryBlockByAddress address)
        {
            return Math.Min(this.End, address.End) - Math.Max(this.Address, address.Address);
        }

        public override string ToString()
        {
            return this.Address.ToString() + " (" + this.Size.ToString() + ")";
        }
    }

    public class MemoryBlockByAddress : IComparable<MemoryBlockByAddress>
    {
        public long Address;
        public long Size;
        public long End { get { return this.Address + this.Size; } }

        public MemoryBlockByAddress()
        {
        }

        public MemoryBlockByAddress(long address, long size)
        {
            this.Address = address;
            this.Size = size;
        }

        public MemoryBlockBySize ToSize { get { return new MemoryBlockBySize(this.Address, this.Size); } }

        public int CompareTo(MemoryBlockByAddress other)
        {
            return this.Address.CompareTo(other.Address);
        }

        public override string ToString()
        {
            return this.Address.ToString() + " (" + this.Size.ToString() + ")";
        }

        public long GetOverlap(MemoryBlockByAddress address)
        {
            return Math.Min(this.End, address.End) - Math.Max(this.Address, address.Address);
        }
    }

    public class MemoryAllocation
    {
        private byte[] memory;
        private AVLTreeSorted<MemoryBlockBySize> structureBySize = new AVLTreeSorted<MemoryBlockBySize>();
        private AVLTreeSorted<MemoryBlockByAddress> structureByAddress = new AVLTreeSorted<MemoryBlockByAddress>();
        private AVLTreeSorted<MemoryBlockByAddress> allocationsByAddress = new AVLTreeSorted<MemoryBlockByAddress>();
        private long free;

        public int AllocatedBlocksCount { get { return this.allocationsByAddress.Count; } }
        public int FreeBlocksCount { get { return this.structureByAddress.Count; } }

        public long FreeBytes { get { return this.free; } }

        public MemoryAllocation(long bytes)
        {
            this.memory = new byte[bytes];
            this.structureByAddress.Add(new MemoryBlockByAddress(0L, bytes));
            this.structureBySize.Add(new MemoryBlockBySize(0L, bytes));
            this.free = bytes;
        }

        public MemoryBlockByAddress Allocate(long bytes)
        {
            if (bytes <= 0L)
            {
                return new MemoryBlockByAddress(0L, 0L);
            }
            MemoryBlockBySize existing, sizeRequest = new MemoryBlockBySize(0L, bytes);
            this.structureBySize.GetLeftmostGreaterOrEqualIndex(sizeRequest, out existing);
            while (existing.Size < bytes)
            {
                this.assureFreeMemory(bytes);
                this.structureBySize.GetLeftmostGreaterOrEqualIndex(sizeRequest, out existing);
            }
            MemoryBlockByAddress result = new MemoryBlockByAddress(existing.Address, bytes);
            MemoryBlockBySize replacer = new MemoryBlockBySize(existing.Address + bytes, existing.Size - bytes);
            this.structureBySize.Remove(existing);
            this.structureByAddress.Remove(existing.ToAddress);
            if (replacer.Size > 0L)
            {
                this.structureBySize.Add(replacer);
                this.structureByAddress.Add(replacer.ToAddress);
            }
            this.free -= bytes;
            //this.memFill(result.Address, result.Size, (byte)result.Size);
            this.allocationsByAddress.Add(result);
            return result;
        }

        private struct MemoryRange
        {
            public long TotalFreeBytes;
            public long TotalCost;

            public MemoryRange(long totalFreeBytes, long totalCost)
            {
                this.TotalFreeBytes = totalFreeBytes;
                this.TotalCost = totalCost;
            }
        }

        private void assureFreeMemory(long bytes)
        {
            if (this.free < bytes)
            {
                throw new OutOfMemoryException();
            }
            this.assureFreeMemoryByMinimumMemoryCopy(bytes);
            //this.assureFreeMemoryByStartAccumulation(bytes);
        }

        private void assureFreeMemoryByMinimumMemoryCopy(long bytes)
        {
            int n = this.structureByAddress.Count;
            MemoryRange[] memory = new MemoryRange[n];

            int storeIndex = 0;
            long lastPointer = 0L, totalFreeSoFar = 0L, totalCostSoFar = 0L;
            foreach (MemoryBlockByAddress block in this.structureByAddress)
            {
                long occupiedBytes = block.Address - lastPointer;
                totalFreeSoFar += block.Size;
                totalCostSoFar += occupiedBytes;
                memory[storeIndex++] = new MemoryRange(totalFreeSoFar, totalCostSoFar);
                lastPointer = block.End;
            }

            if (totalFreeSoFar != this.free)
            {
                throw new InvalidOperationException("Free memory mismatch.");
            }

            int startIndex = -1, endIndex = -1;
            long maxFreeBlockCost = long.MaxValue;
            int j = memory.BinarySearchLeftmostGreaterOrEqual(bytes, mr => mr.TotalFreeBytes);
            for (int i = 0; i < n && j < n; i++)
            {
                long cost = memory[j].TotalCost - memory[i].TotalCost;
                long freeBytes = memory[j].TotalFreeBytes - (i == 0 ? 0L : memory[i - 1].TotalFreeBytes);
                if (freeBytes >= bytes && cost < maxFreeBlockCost)
                {
                    startIndex = i;
                    endIndex = j;
                    maxFreeBlockCost = cost;
                }
                while (j < n && memory[j].TotalFreeBytes - memory[i].TotalFreeBytes < bytes)
                {
                    j++;
                }
            }

            long actualCost = 0L;
            MemoryBlockByAddress resultBlock = this.structureByAddress[startIndex];

            for (int i = startIndex + 1; i <= endIndex; i++)
            {
                MemoryBlockByAddress foundNode, block = this.structureByAddress[i];
                long count = block.Address - resultBlock.End;

                this.memoryMove(resultBlock.End, resultBlock.Address, count);
                int sti = this.allocationsByAddress.GetLeftmostGreaterOrEqualIndex(new MemoryBlockByAddress(resultBlock.End, 0L), out foundNode);
                int eni = this.allocationsByAddress.GetLeftmostGreaterOrEqualIndex(new MemoryBlockByAddress(block.Address, 0L), out foundNode);
                actualCost += count;
                long delta = resultBlock.Size;
                for (int k = eni; --k >= sti;)
                {
                    this.allocationsByAddress[k].Address -= delta;
                }
                resultBlock = new MemoryBlockByAddress(block.Address - resultBlock.Size, resultBlock.Size + block.Size);
            }

            if (actualCost != maxFreeBlockCost)
            {
                Debugger.Break();
            }

            for (int i = endIndex + 1; --i >= startIndex;)
            {
                MemoryBlockBySize memoryBlockSize = this.structureByAddress[i].ToSize;
                this.structureBySize.Remove(memoryBlockSize);
                this.structureByAddress.RemoveAt(i);
            }
            this.structureByAddress.Add(resultBlock);
            this.structureBySize.Add(resultBlock.ToSize);
        }

        private void assureFreeMemoryByStartAccumulation(long bytes)
        {
            MemoryBlockByAddress freeBlock = new MemoryBlockByAddress(-1L, 0L);
            int removedBlocks = 0;
            foreach (MemoryBlockByAddress block in this.structureByAddress)
            {
                if (freeBlock.Address < 0L)
                {
                    freeBlock = block;
                    removedBlocks++;
                    continue;
                }
                long occupiedBytes = block.Address - freeBlock.End;
                if (occupiedBytes < 0)
                {
                    throw new InvalidOperationException("Allocated memory is overlapping.");
                }
                if (occupiedBytes == 0)
                {
                    freeBlock.Size += block.Size;
                }
                else
                {
                    this.memoryMove(freeBlock.End, freeBlock.Address, occupiedBytes);
                    MemoryBlockByAddress foundNode;
                    int startIndex = this.allocationsByAddress.GetLeftmostGreaterOrEqualIndex(new MemoryBlockByAddress(freeBlock.End, 0L), out foundNode);
                    int endIndexExclusive = this.allocationsByAddress.GetLeftmostGreaterOrEqualIndex(new MemoryBlockByAddress(block.Address, 0L), out foundNode);
                    long delta = freeBlock.Size;
                    for (int i = endIndexExclusive; --i >= startIndex;)
                    {
                        this.allocationsByAddress[i].Address -= delta;
                    }
                    freeBlock = new MemoryBlockByAddress(block.Address - freeBlock.Size, block.Size + freeBlock.Size);
                }
                removedBlocks++;
                if (freeBlock.Size >= bytes)
                {
                    break;
                }
            }
            for (int i = removedBlocks; --i >= 0;)
            {
                MemoryBlockBySize blockBySize = this.structureByAddress[i].ToSize;
                this.structureBySize.Remove(blockBySize);
                this.structureByAddress.RemoveAt(i);
            }
            this.structureByAddress.Add(freeBlock);
            this.structureBySize.Add(freeBlock.ToSize);
        }

        private void check()
        {
            var l1 = this.structureByAddress.ToList();
            var l2 = this.structureBySize.OrderBy(x => x.Address).ToList();
            for (int i = Math.Max(l1.Count, l2.Count); --i >= 0;)
            {
                if (l1[i].Address != l2[i].Address || l1[i].Size != l2[i].Size)
                {
                    Debugger.Break();
                }
                if (i > 0 && l1[i - 1].End > l1[i].Address)
                {
                    Debugger.Break();
                }
            }
        }

        private void memFill(long destinationIndex, long count, byte value)
        {
            for (long i = count; --i >= 0;)
            {
                this.memory[destinationIndex + i] = value;
            }
        }

        private bool memCheck(long destinationIndex, long count, byte value)
        {
            for (long i = count; --i >= 0;)
            {
                if (this.memory[destinationIndex + i] != value)
                {
                    return false;
                }
            }
            return true;
        }

        private void memoryMove(long sourceIndex, long destinationIndex, long count)
        {
            if (sourceIndex == destinationIndex || count <= 0)
            {
                return;
            }
            if (destinationIndex > sourceIndex)
            {
                for (long i = count; --i >= 0;)
                {
                    this.memory[i + destinationIndex] = this.memory[i + sourceIndex];
                }
            }
            else
            {
                for (long i = 0; i < count; i++)
                {
                    this.memory[i + destinationIndex] = this.memory[i + sourceIndex];
                }
            }
        }

        public void Free(MemoryBlockByAddress memoryBlock)
        {
            if (memoryBlock.Size <= 0L)
            {
                return;
            }
            //if (!memCheck(memoryBlock.Address, memoryBlock.Size, (byte)memoryBlock.Size))
            //{
            //    Debugger.Break();
            //}
            bool ok = this.allocationsByAddress.Remove(memoryBlock);
            if (!ok)
            {
                throw new InvalidOperationException("Free called to not allocated memory pointer.");
            }
            MemoryBlockByAddress existing;
            int index = this.structureByAddress.GetLeftmostGreaterOrEqualIndex(memoryBlock, out existing);
            MemoryBlockByAddress previous = index <= 0 ? new MemoryBlockByAddress(-1L, 0L) : this.structureByAddress[index - 1];
            long overlap = memoryBlock.GetOverlap(previous);
            if (previous.End == memoryBlock.Address)
            {
                memoryBlock = new MemoryBlockByAddress(previous.Address, previous.Size + memoryBlock.Size);
                this.structureByAddress.Remove(previous);
                this.structureBySize.Remove(previous.ToSize);
                this.free -= previous.Size;
            }
            if (existing != null)
            {
                overlap = Math.Max(overlap, memoryBlock.GetOverlap(existing));
                if (memoryBlock.End == existing.Address)
                {
                    memoryBlock = new MemoryBlockByAddress(memoryBlock.Address, memoryBlock.Size + existing.Size);
                    this.structureByAddress.Remove(existing);
                    this.structureBySize.Remove(existing.ToSize);
                    this.free -= existing.Size;
                }
            }
            this.structureBySize.Add(memoryBlock.ToSize);
            this.structureByAddress.Add(memoryBlock);
            this.free += memoryBlock.Size;
            if (overlap > 0)
            {
                Debugger.Break();
                throw new InvalidOperationException("Overlap between dispose and already free memory.");
            }
        }
    }

    public static class MemoryUnitTest
    {
        public static bool UnitTest(int seed)
        {
            Random random = new Random(seed);
            int totalBytes = 256 * 1024;
            MemoryAllocation memory = new MemoryAllocation(totalBytes);
            List<MemoryBlockByAddress> allocations = new List<MemoryBlockByAddress>();
            for (int i = 20 * 1000; --i >= 0;)
            {
                bool allocate = allocations.Count == 0 || (memory.FreeBytes > 0 && random.Next(2) == 0);
                if (allocate)
                {
                    allocations.Add(memory.Allocate(1 + random.Next((int)memory.FreeBytes)));
                }
                else
                {
                    int index = random.Next(allocations.Count);
                    memory.Free(allocations[index]);
                    allocations.RemoveAt(index);
                }
            }
            foreach (var item in allocations)
            {
                memory.Free(item);
            }
            if (memory.AllocatedBlocksCount != 0 || memory.FreeBlocksCount != 1 || memory.FreeBytes != totalBytes)
            {
                return false;
            }
            return true;
        }
    }
}
