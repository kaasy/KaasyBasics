using System;
using System.Collections.Generic;
using System.Linq;
using System.Net.Sockets;
using System.Text;
using System.Numerics;
using Bitsy.Core;
using System.Collections.Concurrent;
using System.IO;
using System.Xml.Serialization;

namespace Utilities
{
    public enum NetworkCommand : int
    {
        Ping = 0,
        PingBack = 1,
        RegisterClient = 2,
        ComputeECPoint = 3,
        Response_ComputeECPoint = 4,
        SolveServerChallenge = 5,
        FollowsServerChallengeSolution = 6,
        IncomingMessageReplyToAll = 7,
        IncomingMessageNoReply = 8,
        RegisterNewUserName = 9
    }

    public class NetworkProtocol : IDisposable
    {
        private ECCBaseClass ecc;
        private Socket socket;
        private System.Security.Cryptography.AesManaged aes = new System.Security.Cryptography.AesManaged();
        private System.Security.Cryptography.ICryptoTransform encryptor;
        private System.Security.Cryptography.ICryptoTransform decryptor;
        private BigInteger serverChallenge = BigInteger.One;
        private BigInteger clientSecretKey = BigInteger.One;
        public bool EncryptionAcknowledged { get; private set; }
        public Action<NetworkProtocol, byte[]> ProcessBytesCommand = (networkProtocol, bytes) => { };
        public Action<NetworkProtocol, string> ProcessStringCommand = (networkProtocol, stringArgument) => { };
        public Action<NetworkProtocol, string> LogString = (networkProtocol, stringArgument) => { };
        public string RegisteredUserName { get; private set; }

        public NetworkCommand ReceiveCommand()
        {
            NetworkCommand command = (NetworkCommand)this.ReceiveInt();
            return command;
        }

        public long ReceiveLong()
        {
            byte[] data = this.ReceiveBytes(8);
            long result = BitConverter.ToInt64(data, 0);
            return result;
        }

        public int ReceiveInt()
        {
            byte[] data = this.ReceiveBytes(4);
            int result = BitConverter.ToInt32(data, 0);
            return result;
        }

        public byte[] ReceiveBytes(int bytesSize)
        {
            byte[] data = new byte[bytesSize];
            int ofs = 0;
            while (ofs < bytesSize)
            {
                SocketError error;
                int size = this.socket.Receive(data, ofs, bytesSize - ofs, SocketFlags.None, out error);
                if (size <= 0)
                {
                    throw new SocketException(-2);
                }
                ofs += Math.Max(size, 0);
                if (error != SocketError.Success)
                {
                    throw new SocketException((int)error);
                }
            }
            return data;
        }

        public NetworkProtocol(Socket client, ECCBaseClass ecc)
        {
            this.ecc = ecc;
            this.socket = client;
            this.aes.KeySize = 256;
            this.aes.BlockSize = 128;
            this.aes.FeedbackSize = 128;
            this.aes.Mode = System.Security.Cryptography.CipherMode.CBC;
            this.aes.Padding = System.Security.Cryptography.PaddingMode.ISO10126;
            this.aes.Key = new byte[256 / 8];
            this.aes.IV = new byte[128 / 8];
            this.encryptor = aes.CreateEncryptor();
            this.decryptor = aes.CreateDecryptor();
        }

        private void SendBytes(byte[] data)
        {
            int ofs = 0;
            while (ofs < data.Length)
            {
                SocketError error;
                int size = this.socket.Send(data, ofs, data.Length - ofs, SocketFlags.None, out error);
                ofs += Math.Max(size, 0);
                if (error != SocketError.Success)
                {
                    throw new Exception("Send SocketError: " + error.ToString());
                }
            }
        }

        public void SendCommand(NetworkCommand command, IEnumerable<byte> bytes)
        {
            byte[] flux = BitConverter.GetBytes((int)command).Concat(bytes).ToArray();
            this.SendBytes(flux);
        }

        public void SendCommand(NetworkCommand command)
        {
            this.SendBytes(BitConverter.GetBytes((int)command));
        }

        public void SendCommand(NetworkCommand command, int parameter1)
        {
            byte[] flux = new byte[8];
            BitConverter.GetBytes((int)command).CopyTo(flux, 0);
            BitConverter.GetBytes(parameter1).CopyTo(flux, 4);
            this.SendBytes(flux);
        }

        public void SendCommand(NetworkCommand command, int parameter1, int parameter2)
        {
            byte[] flux = new byte[12];
            BitConverter.GetBytes((int)command).CopyTo(flux, 0);
            BitConverter.GetBytes(parameter1).CopyTo(flux, 4);
            BitConverter.GetBytes(parameter2).CopyTo(flux, 8);
            this.SendBytes(flux);
        }

        public void SendCommand(NetworkCommand command, long parameter1)
        {
            byte[] flux = new byte[12];
            BitConverter.GetBytes((int)command).CopyTo(flux, 0);
            BitConverter.GetBytes(parameter1).CopyTo(flux, 4);
            this.SendBytes(flux);
        }

        public void Dispose()
        {
            if (this.socket != null)
            {
                try
                {
                    this.socket.Disconnect(true);
                }
                catch { }
                this.socket.Close();
                this.socket.Dispose();
                this.socket = null;
            }
        }

        public NetworkCommand LastReceivedCommand { get; private set; }
        public NetworkCommand ProcessCommand()
        {
            NetworkCommand command = this.ReceiveCommand();
            this.LastReceivedCommand = command;
            switch (command)
            {
                case NetworkCommand.Ping:
                    {
                        long time = this.ReceiveLong();
                        this.SendCommand(NetworkCommand.PingBack, time);
                    }
                    break;
                case NetworkCommand.PingBack:
                    {
                        long time = this.ReceiveLong();
                        TimeSpan elapsed = TimeSpan.FromTicks(DateTime.Now.Ticks - time);
                        this.LogString(this, "Returning ping time: " + elapsed.ToString() + "\r\n");
                    }
                    break;
                case NetworkCommand.ComputeECPoint:
                    {
                        var factor = ecc.GetRandomFactorModN();
                        var multiplication = ecc.ECMultiplication(factor);
                        var factorBytes = factor.ToByteArray();
                        Array.Resize(ref factorBytes, ecc.BytesCount);
                        var eccBytes = multiplication.UncompressedBytes;
                        this.SendCommand(NetworkCommand.Response_ComputeECPoint, factorBytes.Concat(eccBytes));
                    }
                    break;
                case NetworkCommand.Response_ComputeECPoint:
                    {
                        byte[] data = this.ReceiveBytes(ecc.BytesCount * 3);
                        var factor = data.Take(ecc.BytesCount).ToBigInteger();
                        var verification = ecc.ECMultiplication(factor).UncompressedBytes;
                        bool ok = Enumerable.Range(0, verification.Length).All(idx => data[ecc.BytesCount + idx] == verification[idx]);
                        if (!ok)
                        {
                            throw new InvalidOperationException("ECC verification failed!!");
                        }
                    }
                    break;
                case NetworkCommand.RegisterClient:
                    {
                        var point = new ECC521Point(this.ReceiveBytes(ecc.CompressedBytesCount), true);
                        string clientIdentifiers = this.ReadString();
                        this.serverChallenge = ecc.GetRandomFactorModN();
                        var clientChallengePoint = point.Multiply(this.serverChallenge);
                        string userName = ConfigurationFile.LoadFromFile().GetUserName(clientIdentifiers);
                        var currentTime = DateTime.Now.Ticks;
                        this.SendCommand(NetworkCommand.SolveServerChallenge,
                            clientChallengePoint.UncompressedBytes.
                            Concat(BitConverter.GetBytes(currentTime)).
                            Concat(userName.StringToVariableLengthBytes()));
                    }
                    break;
                case NetworkCommand.SolveServerChallenge:
                    {
                        var point = new ECC521Point(this.ReceiveBytes(ecc.BytesCount * 2), false);
                        var serverTimeBytes = this.ReceiveBytes(8);
                        this.RegisteredUserName = this.ReadString();
                        var secretPoint = point.Multiply(this.clientSecretKey.ModInverse(ecc.N));
                        var solution = RIPEMD160.Hash(secretPoint.CompressedBytes);
                        this.SendCommand(NetworkCommand.FollowsServerChallengeSolution, solution.Concat(serverTimeBytes));
                        this.SetSharedKey(SHA256.Hash(secretPoint.UncompressedBytes), solution.Take(128 / 8).ToArray());
                    }
                    break;
                case NetworkCommand.FollowsServerChallengeSolution:
                    {
                        byte[] ripemd160Solution = this.ReceiveBytes(160 / 8);
                        var challengeDuration = TimeSpan.FromTicks(DateTime.Now.Ticks - BitConverter.ToInt64(this.ReceiveBytes(8), 0));
                        var secretPoint = ecc.G.Multiply(serverChallenge);
                        byte[] verification = RIPEMD160.Hash(secretPoint.CompressedBytes);
                        if (!ripemd160Solution.IsSameAs(verification))
                        {
                            throw new InvalidOperationException("Authentication failed.");
                        }
                        this.SetSharedKey(SHA256.Hash(secretPoint.UncompressedBytes), verification.Take(128 / 8).ToArray());
                        this.LogString(this, "Client challenge duration: " + challengeDuration.ToString() + "\r\n");
                    }
                    break;
                case NetworkCommand.IncomingMessageReplyToAll:
                case NetworkCommand.IncomingMessageNoReply:
                    {
                        string incomingMessage = this.DecodeIncomingMessage();
                        this.ProcessStringCommand(this, incomingMessage);
                    }
                    break;
                case NetworkCommand.RegisterNewUserName:
                    {
                        string clientIdentifiers = this.ReadString();
                        string userName = this.ReadString();
                        UserInformation userInformation = new UserInformation()
                        {
                            ClientIdentifiers = clientIdentifiers,
                            UserName = userName
                        };
                        var configurationFile = ConfigurationFile.LoadFromFile();
                        configurationFile.AddUser(userInformation);
                        configurationFile.SaveToFile();
                    }
                    break;
            }
            return command;
        }

        public void RegisterNewUserName(string mPhoneFullIdentifiers, string userName)
        {
            this.SendCommand(NetworkCommand.RegisterNewUserName,
                mPhoneFullIdentifiers.StringToVariableLengthBytes().
                Concat(userName.StringToVariableLengthBytes()));
            this.RegisteredUserName = userName;
        }

        private void SetSharedKey(byte[] sharedKey, byte[] iv)
        {
            this.aes.Key = sharedKey;
            this.aes.IV = iv;
            this.encryptor = aes.CreateEncryptor();
            this.decryptor = aes.CreateDecryptor();
            this.EncryptionAcknowledged = true;
        }

        public void RegisterNewClient(string clientIdentifiers)
        {
            this.clientSecretKey = ecc.GetRandomFactorModN();
            var publicKey = ecc.ECMultiplication(this.clientSecretKey);
            this.SendCommand(NetworkCommand.RegisterClient,
                publicKey.CompressedBytes.Concat(clientIdentifiers.StringToVariableLengthBytes()));
        }

        public string DecodeIncomingMessage()
        {
            int length = this.ReceiveInt();
            byte[] encryptedMessage = this.ReceiveBytes(length);
            byte[] decodedMessage = ApplyCryptoTransform(this.decryptor, encryptedMessage);
            byte[] incomingChecksum = Enumerable.Range(decodedMessage.Length - 256 / 8, 256 / 8).Select(idx => decodedMessage[idx]).ToArray();
            Array.Resize(ref decodedMessage, decodedMessage.Length - 256 / 8);
            byte[] currentChecksum = SHA256.Hash(decodedMessage);
            if (!currentChecksum.IsSameAs(incomingChecksum))
            {
                throw new InvalidOperationException("Error in communication - checksum error.");
            }
            string stringArgument = Encoding.UTF8.GetString(decodedMessage);
            return stringArgument;
        }

        public void SendEncryptedMessage(string message, NetworkCommand command)
        {
            message = message == null ? "" : message;
            byte[] messageBytes = Encoding.UTF8.GetBytes(message);
            byte[] checksum = SHA256.Hash(messageBytes);
            byte[] communicationToEncrypt = messageBytes.Concat(checksum).ToArray();
            byte[] encryptedCommunication = ApplyCryptoTransform(this.encryptor, communicationToEncrypt);
            this.SendCommand(command,
                BitConverter.GetBytes(encryptedCommunication.Length).Concat(encryptedCommunication).ToArray());
        }

        public static byte[] ApplyCryptoTransform(System.Security.Cryptography.ICryptoTransform transform, byte[] inputData)
        {
            byte[] result = transform.TransformFinalBlock(inputData, 0, inputData.Length);
            return result;
        }
    }

    public class UserInformation
    {
        public string ClientIdentifiers { get; set; }
        public string UserName { get; set; }
    }

    public class ConfigurationFile
    {
        private const string ConfigurationFileName = "RoundChat.config.xml";
        private static XmlSerializer xmlSerializer = new XmlSerializer(typeof(ConfigurationFile));

        public UserInformation[] Users = new UserInformation[0];

        public string GetUserName(string clientIdentifiers)
        {
            string[] identifiers = clientIdentifiers.Split(new char[] { ' ' }, StringSplitOptions.RemoveEmptyEntries);
            HashSet<string> hashIdentifiers = new HashSet<string>(identifiers);
            foreach (UserInformation userInformation in this.Users)
            {
                identifiers = userInformation.ClientIdentifiers.Split(new char[] { ' ' }, StringSplitOptions.RemoveEmptyEntries);
                if (identifiers.Any(hashIdentifiers.Contains))
                {
                    return userInformation.UserName;
                }
            }
            return "";
        }

        public void AddUser(UserInformation userInformation)
        {
            Array.Resize(ref this.Users, this.Users.Length + 1);
            this.Users[this.Users.Length - 1] = userInformation;
        }

        public void SaveToFile()
        {
            StringBuilder contents = new StringBuilder();
            using (var writer = new StringWriter(contents))
            {
                xmlSerializer.Serialize(writer, this);
            }
            File.WriteAllText(ConfigurationFileName, contents.ToString());
        }

        public static ConfigurationFile LoadFromFile()
        {
            if (!File.Exists(ConfigurationFileName))
            {
                return new ConfigurationFile();
            }
            try
            {
                string contents = File.ReadAllText(ConfigurationFileName);
                using (var reader = new StringReader(contents))
                {
                    ConfigurationFile configurationFile = xmlSerializer.Deserialize(reader) as ConfigurationFile;
                    return configurationFile;
                }
            }
            catch
            {
                return new ConfigurationFile();
            }
        }
    }

    public static class NetworkProtocolExtensions
    {
        public static byte[] StringToVariableLengthBytes(this string text)
        {
            var result = BitConverter.GetBytes(text.Length).Concat(Encoding.UTF8.GetBytes(text)).ToArray();
            return result;
        }

        public static string ReadString(this NetworkProtocol protocol)
        {
            int length = protocol.ReceiveInt();
            byte[] bytes = protocol.ReceiveBytes(length);
            string result = Encoding.UTF8.GetString(bytes);
            return result;
        }

        public static string VariableLengthBytesToString(this byte[] bytes)
        {
            int length = BitConverter.ToInt32(bytes, 0);
            var result = Encoding.UTF8.GetString(bytes, 0, length);
            return result;
        }
    }
}
