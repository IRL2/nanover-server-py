using System;
using System.Collections.Generic;
using System.Net;
using System.Net.Sockets;
using System.Text;
using System.Threading.Tasks;
using System.Timers;

namespace Essd
{
    /// <summary>
    /// Implementation of an Extremely Simple Service Discovery (ESSD) client. 
    /// </summary>
    public class Client
    {
        /// <summary>
        /// Default port upon which ESSD clients listen for services.
        /// </summary>
        public const int DefaultListenPort = 54545;
        
        private UdpClient udpClient;


        /// <summary>
        /// Initialises an ESSD client.
        /// </summary>
        /// <param name="listenPort">Port at which to listen for services.</param>
        public Client(int listenPort=DefaultListenPort)
        {
            udpClient = new UdpClient();
    
            udpClient.Client.Bind(new IPEndPoint(IPAddress.Any, listenPort));
        }

        /// <summary>
        /// Searches for services for the given duration, blocking.
        /// </summary>
        /// <param name="duration">Duration to search for, in milliseconds.</param>
        /// <returns>Collection of unique services found during search.</returns>
        /// <remarks>
        /// There is no guarantee that services found during the search will still
        /// be up after the search ends.
        /// </remarks>
        public ICollection<ServiceHub> SearchForServices(int duration=3000)
        {
            HashSet<ServiceHub> servicesFound = new HashSet<ServiceHub>();

            var from = new IPEndPoint(0, 0);
            
            var timer = new Timer(duration);
            udpClient.Client.ReceiveTimeout = duration;
            timer.Elapsed += OnTimerElapsed;
            bool timerElapsed = false;
            timer.Start();
            
            while(!timerElapsed)
            {
                byte[] messageBytes;
                try
                {
                    messageBytes = udpClient.Receive(ref from);
                }
                catch (SocketException e)
                {
                    if(e.SocketErrorCode == SocketError.TimedOut)
                        break;
                    throw;
                }
                
                var message = DecodeMessage(messageBytes);
                if (message is null)
                {
                    Console.WriteLine($"ESSD: Received invalid message, not a valid UTF8 string.");
                    continue;
                }

                ServiceHub service;
                try
                {
                    service = new ServiceHub(message);
                }
                catch (ArgumentException)
                {
                    Console.WriteLine($"ESSD: Received invalid service definition.");
                    continue;
                }

                servicesFound.Add(service);
            }

            return servicesFound;

            void OnTimerElapsed(object sender, ElapsedEventArgs e)
            {
                timerElapsed = true;
            }
            
        }



        private string DecodeMessage(byte[] messageBytes)
        {
            try
            {
                var message = Encoding.UTF8.GetString(messageBytes);
                return message;
            }
            catch (ArgumentException)
            {
                return null;
            }
            
        }


    }

    /// <summary>
    ///     A definition of a ServiceHub that can be discovered or broadcast.
    /// 
    ///     A service hub consists of properties that must at least consist of a name and ip address.
    ///     The payload can optionally include additional information on the services provided.
    /// </summary>
    public class ServiceHub : IEquatable<ServiceHub>
    {

        /// <summary>
        /// The key used for the service name field, which is always required.
        /// </summary>
        public const string NameKey = "name";
        
        /// <summary>
        /// The key used for the service address field, which is always required.
        /// </summary>
        public const string AddressKey = "address";
        
        /// <summary>
        /// The name of the service.
        /// </summary>
        public string Name => (string) Properties[NameKey];
        
        /// <summary>
        /// The address of the service.
        /// </summary>
        public string Address => (string) Properties[AddressKey];
        
        private Dictionary<string, object> Properties;

        /// <summary>
        /// Initialises a service hub from a JSON string describing the service hub.
        /// </summary>
        /// <param name="serviceHubJson">JSON string describing the service.</param>
        public ServiceHub(string serviceHubJson)
        {
            Properties = Newtonsoft.Json.JsonConvert.DeserializeObject<Dictionary<string, object>>(serviceHubJson);
            ValidateProperties(Properties);
        }

        private void ValidateProperties(Dictionary<string, object> properties)
        {
            try
            {
                var name = properties[NameKey];
                var address = properties[AddressKey];
            }
            catch (KeyNotFoundException e)
            {
                throw new ArgumentException(e.Message);
            }
        }

        private void ValidateField(Dictionary<string, object> properties, string key)
        {
            try
            {
                var field = properties[key];
            }
            catch (KeyNotFoundException)
            {
                throw new ArgumentException($"Field {key} not found in service hub definition.");
            }
        }


        /// <inheritdoc />
        /// <remarks>
        /// For the purposes of discovery, two service hubs are considered equal if
        /// they have the same name and are serving at the same address.
        /// </remarks>
        /// TODO remove this constraint? It doesn't really make sense.
        public bool Equals(ServiceHub other)
        {
            if (ReferenceEquals(null, other)) return false;
            if (ReferenceEquals(this, other)) return true;
            return Name.Equals(other.Name) || Address.Equals(other.Address);
        }

        /// <inheritdoc />
        public override bool Equals(object obj)
        {
            if (ReferenceEquals(null, obj)) return false;
            if (ReferenceEquals(this, obj)) return true;
            if (obj.GetType() != this.GetType()) return false;
            return Equals((ServiceHub) obj);
        }

        /// <inheritdoc />
        public override int GetHashCode()
        {
            return (Name + Address).GetHashCode();
        }
    }
}