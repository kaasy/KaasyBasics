using System;
using System.Collections.Generic;
using System.Linq;
using System.Runtime.InteropServices;
using System.Text;
using System.Threading.Tasks;

namespace Utilities
{
    public enum RegisterLocation
    {
        Input,
        Work,
        Output,
        Memory,
        Jump
    }

    public enum RegisterType
    {
        Integer,
        Real,
        Complex
    }

    public enum ComplexRegisterPart
    {
        Full,
        Real,
        Imaginary,
        Magnitude,
        Phase
    }

    public struct RegisterReference
    {
        public RegisterLocation RegisterLocation;
        public RegisterType RegisterType;
        public ComplexRegisterPart ComplexRegisterPart;
        public int RegisterIndex;
        public string JumpLabel;

        public RegisterReference(RegisterLocation location, RegisterType type, int index, ComplexRegisterPart part)
        {
            this.RegisterLocation = location;
            this.RegisterType = type;
            this.RegisterIndex = index;
            this.ComplexRegisterPart = part;
            this.JumpLabel = null;
        }

        private static List<Pair> registersComplexParts = new List<Pair>() {
            new Pair(".real", "Real"),
            new Pair(".imag", "Imaginary"),
            new Pair(".phase", "Phase"),
            new Pair(".abs", "Magnitude")
        };

        public RegisterReference(string compact)
        {
            string name = compact.ToLowerInvariant().Trim(' ', '\t');
            this.ComplexRegisterPart = ComplexRegisterPart.Full;

            foreach (Pair pair in registersComplexParts)
            {
                if (name.EndsWith(pair.Key))
                {
                    ComplexRegisterPart currentPart;
                    if (Enum.TryParse(pair.Member, out currentPart))
                    {
                        this.ComplexRegisterPart = currentPart;
                    }
                    name = name.Remove(name.Length - pair.Key.Length);
                    break;
                }
            }
            while (name.Length > 1 && char.IsDigit(name[name.Length - 1]))
            {
                name = name.Remove(name.Length - 1);
            }
            //ii, ir, ic, i, r, c, oi, or, oc
            switch (name)
            {
                case "ii": this.RegisterLocation = RegisterLocation.Input; this.RegisterType = RegisterType.Integer; break;
                case "ir": this.RegisterLocation = RegisterLocation.Input; this.RegisterType = RegisterType.Real; break;
                case "ic": this.RegisterLocation = RegisterLocation.Input; this.RegisterType = RegisterType.Complex; break;
                case "i": this.RegisterLocation = RegisterLocation.Work; this.RegisterType = RegisterType.Integer; break;
                case "r": this.RegisterLocation = RegisterLocation.Work; this.RegisterType = RegisterType.Real; break;
                case "c": this.RegisterLocation = RegisterLocation.Work; this.RegisterType = RegisterType.Complex; break;
                case "oi": this.RegisterLocation = RegisterLocation.Output; this.RegisterType = RegisterType.Integer; break;
                case "or": this.RegisterLocation = RegisterLocation.Output; this.RegisterType = RegisterType.Real; break;
                case "oc": this.RegisterLocation = RegisterLocation.Output; this.RegisterType = RegisterType.Complex; break;
                default: throw new InvalidOperationException();
            }
            this.RegisterIndex = int.Parse(compact.Substring(name.Length));
            this.JumpLabel = null;
        }

        public RegisterReference(int jumpLine, string label)
        {
            this.ComplexRegisterPart = ComplexRegisterPart.Full;
            this.RegisterLocation = RegisterLocation.Input;
            this.RegisterType = RegisterType.Integer;
            this.RegisterIndex = jumpLine;
            this.JumpLabel = label;
        }
    }

    public class ASMRegisters
    {
        public IntegerNumber[][] IntegerRegisters = null;
        public RealNumber[][] RealRegisters = null;
        public ComplexNumber[][] ComplexRegisters = null;
        public ASMRegisters()
        {
        }

        public void Apply(bool isFloatingPoint, int subregistersCount, int realCount, int complexCount)
        {
            if (!isFloatingPoint)
            {
                Array.Resize(ref this.RealRegisters, 0);
                Array.Resize(ref this.ComplexRegisters, 0);
                Array.Resize(ref this.IntegerRegisters, realCount);
                for (int i = realCount; --i >= 0;)
                {
                    Array.Resize(ref this.IntegerRegisters[i], subregistersCount);
                }
            }
            else
            {
                Array.Resize(ref this.IntegerRegisters, 0);
                Array.Resize(ref this.RealRegisters, realCount);
                for (int i = realCount; --i >= 0;)
                {
                    Array.Resize(ref this.RealRegisters[i], subregistersCount);
                }
                Array.Resize(ref this.ComplexRegisters, complexCount);
                for (int i = complexCount; --i >= 0;)
                {
                    Array.Resize(ref this.ComplexRegisters[i], subregistersCount);
                }
            }
        }
    }

    public class ASMConfig
    {
        public int SubregisterIsFloatingPoint { get; set; }
        public int SubregisterWidthBits { get; set; }
        public int SubregistersPerRegister { get; set; }
        public int NumberOfRealRegisters { get; set; }
        public int NumberOfComplexRegisters { get; set; }
        public int NumberOfInputRealRegisters { get; set; }
        public int NumberOfInputComplexRegisters { get; set; }
        public int NumberOfOutputRealRegisters { get; set; }
        public int NumberOfOutputComplexRegisters { get; set; }
        public int MemorySizeKB { get; set; }
        public int NumberOfExecutionThreads { get; set; }

        public ASMConfig()
        {
            this.SubregisterIsFloatingPoint = 1;
            this.SubregisterWidthBits = 32;
            this.SubregistersPerRegister = 32;
            this.NumberOfRealRegisters = 16;
            this.NumberOfComplexRegisters = 16;
            this.NumberOfInputRealRegisters = 0;
            this.NumberOfInputComplexRegisters = 1;
            this.NumberOfOutputRealRegisters = 1;
            this.NumberOfOutputComplexRegisters = 0;
            this.MemorySizeKB = 1024;
            this.NumberOfExecutionThreads = 32;
        }

        public ASMRegisters Input = new ASMRegisters();
        public ASMRegisters Output = new ASMRegisters();
        public ASMRegisters Work = new ASMRegisters();
        private MemoryAllocation memory;

        public void ApplyConfiguration()
        {
            this.memory = new MemoryAllocation(this.MemorySizeKB << 10);
            this.Input.Apply(this.SubregisterIsFloatingPoint != 0, this.SubregistersPerRegister, this.NumberOfInputRealRegisters, this.NumberOfInputComplexRegisters);
            this.Output.Apply(this.SubregisterIsFloatingPoint != 0, this.SubregistersPerRegister, this.NumberOfOutputRealRegisters, this.NumberOfOutputComplexRegisters);
            this.Work.Apply(this.SubregisterIsFloatingPoint != 0, this.SubregistersPerRegister, this.NumberOfRealRegisters, this.NumberOfComplexRegisters);
        }
    }

    public struct Pair
    {
        public string Key { get; set; }
        public string Member { get; set; }
        public Pair(string key, string member) : this()
        {
            this.Key = key;
            this.Member = member;
        }
    }

    public enum ASMInstructionType
    {
        None,
        Configuration,
        Execution
    }

    public class ASMDecomposer
    {
        private ASMConfig config = new ASMConfig();
        private string[] lines;
        private int instructionPointer = 0;
        private ASMInstructionType asmMode = ASMInstructionType.None;

        private Dictionary<string, int> labels = new Dictionary<string, int>();

        public ASMDecomposer(string asmText)
        {
            this.lines = asmText.Split(new char[] { '\r', '\n' }, StringSplitOptions.RemoveEmptyEntries);
            this.findLabels();
        }

        private void findLabels()
        {
            this.labels.Clear();
            for (int i = this.lines.Length; --i >= 0;)
            {
                string line = this.lines[i].TrimStart(' ');
                int index = line.IndexOf(':');
                if (index < 0)
                {
                    continue;
                }
                this.labels[line.Substring(0, index)] = i;
            }
        }

        private bool executionEnded { get { return this.instructionPointer < 0 || this.instructionPointer >= this.lines.Length; } }

        private bool tryChangeRegion()
        {
            string line = this.lines[this.instructionPointer];
            line = line.TrimStart(' ').ToLowerInvariant();
            if (line.TrimEnd(' ') == "#endregion")
            {
                if (this.asmMode == ASMInstructionType.Configuration)
                {
                    this.config.ApplyConfiguration();
                }
                this.asmMode = ASMInstructionType.None;
                return true;
            }
            string regionToken = "#region ";
            if (line.StartsWith(regionToken))
            {
                switch (line.Substring(regionToken.Length).Trim(' ').ToLowerInvariant())
                {
                    case "config":
                        this.asmMode = ASMInstructionType.Configuration;
                        return true;
                    case "execution":
                        this.asmMode = ASMInstructionType.Execution;
                        return true;
                    default: throw new NotSupportedException();
                }
            }
            return false;
        }

        private static List<Pair> globalSettings = new List<Pair>()
        {
            new Pair("subregister_is_floating_point", "SubregisterIsFloatingPoint"),
            new Pair("subregister_width_bits", "SubregisterWidthBits"),
            new Pair("subregisters_per_register", "SubregistersPerRegister"),
            new Pair("number_of_real_registers", "NumberOfRealRegisters"),
            new Pair("number_of_complex_registers", "NumberOfComplexRegisters"),
            new Pair("number_of_input_real_registers", "NumberOfInputRealRegisters"),
            new Pair("number_of_input_complex_registers", "NumberOfInputComplexRegisters"),
            new Pair("number_of_output_real_registers", "NumberOfOutputRealRegisters"),
            new Pair("number_of_output_complex_registers", "NumberOfOutputComplexRegisters"),
            new Pair("memory_size_kb", "MemorySizeKB"),
            new Pair("number_of_execution_threads", "NumberOfExecutionThreads")
        };

        private static Dictionary<string, System.Reflection.PropertyInfo> hashSettingToMemberInformation;
        private static HashSet<string> opcodes;

        static ASMDecomposer()
        {
            hashSettingToMemberInformation = new Dictionary<string, System.Reflection.PropertyInfo>();
            foreach (Pair pair in globalSettings)
            {
                hashSettingToMemberInformation.Add(pair.Key, typeof(ASMConfig).GetProperty(pair.Member));
            }
            opcodes = new HashSet<string>();
            foreach (string opcode in globalOperations)
            {
                opcodes.Add(opcode);
            }
        }

        private bool tryChangeSetting()
        {
            string line = this.lines[this.instructionPointer];
            line = line.TrimStart(' ').ToLowerInvariant();
            string[] values = line.Split(new char[] { ' ', '\t' }, StringSplitOptions.RemoveEmptyEntries);
            System.Reflection.PropertyInfo propertyInfo;
            if (hashSettingToMemberInformation.TryGetValue(values[0], out propertyInfo))
            {
                propertyInfo.SetValue(this.config, int.Parse(values[1]), null);
                return true;
            }
            return false;
        }

        private static List<string> globalOperations = new List<string>()
        {
            "mov",
            "and",
            "or",
            "xor",
            "not",
            "clmul",

            "add",
            "adc",
            "neg",
            "sub",
            "sbb",
            "mul",
            "div",
            "lldiv",
            "modinv",
            "sqrt",
            "shl",
            "shr",
            "expmod",

            "invsqrt",
            "inverse",
            "exp",
            "log",

            "fft",
            "dct2",
            "dct3",
            "new",

            "jmp_equal",
            "jmp_not_equal",
            "jmp_less",
            "jmp_less_equal",
            "jmp_greater",
            "jmp_greater_equal",
            "jmp_not_less",
            "jmp_not_less_equal",
            "jmp_not_greater",
            "jmp_not_greater_equal"
        };

        private bool tryExecute()
        {
            string line = this.lines[this.instructionPointer];
            line = line.TrimStart(' ').ToLowerInvariant();
            string[] ops = line.Split(new char[] { '\t', ' ' }, StringSplitOptions.RemoveEmptyEntries);
            string mainop = ops[0];
            if (!opcodes.Contains(mainop))
            {
                return false;
            }
            line = string.Join(" ", ops.Skip(1));
            ops = line.Split(new char[] { ',' }, StringSplitOptions.RemoveEmptyEntries);
            for (int i = ops.Length; --i >= 0;)
            {
                ops[i] = ops[i].Trim(' ', '\t');
            }

            List<RegisterReference> registers = new List<RegisterReference>();
            for (int i = 0; i < ops.Length; i++)
            {
                int jumpLocation;
                if (labels.TryGetValue(ops[i], out jumpLocation))
                {
                    registers.Add(new RegisterReference(jumpLocation, ops[i]));
                }
                else
                {
                    registers.Add(new RegisterReference(ops[i]));
                }
            }

            bool ok = doExecute(mainop, registers, this.config);
            return ok;
        }

        public delegate void GenericOpFunction0(IntPtr sourceAndDestination, int log2BitsPerOp, int n);
        public delegate void GenericOpFunction1(IntPtr destination, IntPtr source, int log2BitsPerOp, int n);
        public delegate void GenericOpFunction2(IntPtr source1, IntPtr source2, IntPtr destination, int log2BitsPerOp, int n);

        private static bool doExecute(string mainop, List<RegisterReference> registers, ASMConfig config)
        {
            RegisterReference destination = registers.Count > 0 ? registers[0] : new RegisterReference();
            RegisterReference source = registers.Count > 1 ? registers[1] : new RegisterReference();

            switch (mainop)
            {
                //case "mov":
                //case "and": fn2 = FourierPoint.OpAnd; break;
                //case "or": fn2 = FourierPoint.OpOr; break;
                //case "xor": fn2 = FourierPoint.OpXor; break;
                //case "not": fn0 = FourierPoint.OpNot; break;
                //case "clmul": fn2 = FourierPoint.OpCLMul; break;

                //case "add": fn2 = FourierPoint.OpAdd; fn3 = FourierPoint.OpFPAdd; break;
                //case "adc": throw new NotImplementedException(); break;
                //case "neg": throw new NotImplementedException(); break;
                //case "sub": fn2 = FourierPoint.OpSub; fn3 = FourierPoint.OpFPSub; break;
                //case "sbb": throw new NotImplementedException(); break;
                //case "mul": fn2 = FourierPoint.OpMul; fn3 = FourierPoint.OpFPMul; break;
                //case "div": fn3 = FourierPoint.OpFPDiv; break;
                case "lldiv": throw new NotImplementedException();
                case "modinv": throw new NotImplementedException();
                case "sqrt": throw new NotImplementedException();
                case "shl": throw new NotImplementedException();
                case "shr": throw new NotImplementedException();
                case "expmod": throw new NotImplementedException();

                case "invsqrt": throw new NotImplementedException();
                case "inverse": throw new NotImplementedException();
                case "exp": throw new NotImplementedException();
                case "log": throw new NotImplementedException();

                case "fft": throw new NotImplementedException();
                case "dct2": throw new NotImplementedException();
                case "dct3": throw new NotImplementedException();
                case "new": throw new NotImplementedException();

                case "jmp_equal": break;
                case "jmp_not_equal": break;
                case "jmp_less": break;
                case "jmp_less_equal": break;
                case "jmp_greater": break;
                case "jmp_greater_equal": break;
                case "jmp_not_less": break;
                case "jmp_not_less_equal": break;
                case "jmp_not_greater": break;
                case "jmp_not_greater_equal": break;
            }
            return false;
        }

        private void nextInstruction()
        {
            this.instructionPointer++;
        }

        public void Execute()
        {
            while (!executionEnded)
            {
                if (this.tryChangeRegion())
                {
                    nextInstruction();
                    continue;
                }
                if (this.asmMode == ASMInstructionType.Configuration)
                {
                    if (tryChangeSetting())
                    {
                        nextInstruction();
                        continue;
                    }
                    else
                    {
                        throw new NotSupportedException();
                    }
                }
                else if (this.asmMode == ASMInstructionType.Execution)
                {
                    if (tryExecute())
                    {
                        nextInstruction();
                        continue;
                    }
                    else
                    {
                        throw new NotSupportedException();
                    }
                }
                else if (this.asmMode == ASMInstructionType.None)
                {
                    nextInstruction();
                }
                else
                {
                    throw new NotImplementedException();
                }
            }
        }
    }
}
