using System;
using System.Collections.Generic;
using Newtonsoft.Json;

namespace Essd
{
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
        
        public SortedDictionary<string, object> Properties;

        /// <summary>
        /// Initialises a service hub from a JSON string describing the service hub.
        /// </summary>
        /// <param name="serviceHubJson">JSON string describing the service.</param>
        public ServiceHub(string serviceHubJson)
        {
            Properties = JsonConvert.DeserializeObject<SortedDictionary<string, object>>(serviceHubJson);
            ValidateProperties(Properties);
        }

        public ServiceHub(string name, string address)
        {
            Properties = new SortedDictionary<string, object>();
            Properties[NameKey] = name;
            Properties[AddressKey] = address;
            ValidateProperties(Properties);
        }
        
        private void ValidateProperties(IDictionary<string, object> properties)
        {
            try
            {
                ValidateField(properties, NameKey);
                ValidateField(properties, AddressKey);
            }
            catch (ArgumentException e)
            {
                throw new ArgumentException(e.Message);
            }
        }

        private void ValidateField(IDictionary<string, object> properties, string key)
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
        public bool Equals(ServiceHub other)
        {
            if (ReferenceEquals(null, other)) return false;
            if (ReferenceEquals(this, other)) return true;
            return (Name == other.Name && Address == other.Address);
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

        public string ToJson()
        {
            return JsonConvert.SerializeObject(Properties);
        }

        public override string ToString()
        {
            return $"{Name}:{Address}";
        }
    }
}