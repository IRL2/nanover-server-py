using System;
using System.Collections.Generic;
using Newtonsoft.Json;

namespace Essd
{
    /// <summary>
    ///     A definition of a ServiceHub that can be discovered or broadcast.
    ///     A service hub consists of properties that must at least consist of a name and IP address.
    ///     The payload can optionally include additional information on the services provided.
    /// </summary>
    public class ServiceHub : IEquatable<ServiceHub>
    {
        /// <summary>
        ///     The key used for the service name field, which is always required.
        /// </summary>
        public const string NameKey = "name";

        /// <summary>
        ///     The key used for the service address field, which is always required.
        /// </summary>
        public const string AddressKey = "address";

        /// <summary>
        ///     The key used for the ESSD version field, which will be generated if not provided.
        /// </summary>
        public const string VersionKey = "essd_version";

        /// <summary>
        ///     The key used for the services field.
        /// </summary>
        public const string ServicesKey = "services";

        /// <summary>
        ///     The key used for the ID field.
        /// </summary>
        public const string IdKey = "id";

        /// <summary>
        /// The raw properties of this service hub.
        /// </summary>
        public Dictionary<string, object> Properties;

        /// <summary>
        ///     The name of the service hub.
        /// </summary>
        public string Name => (string) Properties[NameKey];

        /// <summary>
        ///     The address of the service hub.
        /// </summary>
        public string Address => (string) Properties[AddressKey];

        /// <summary>
        ///     The ESSD version of this service hub.
        /// </summary>
        public string Version => (string) Properties[VersionKey];

        /// <summary>
        ///     The ID of the service hub.
        /// </summary>
        public string Id => (string) Properties[IdKey];
        
        /// <summary>
        ///     Initialises a service hub from a JSON string describing the service hub.
        /// </summary>
        /// <param name="serviceHubJson">JSON string describing the service.</param>
        public ServiceHub(string serviceHubJson)
        {
            Properties = JsonConvert.DeserializeObject<Dictionary<string, object>>(serviceHubJson);
            ValidateProperties(Properties);

            if (!Properties.ContainsKey(VersionKey))
                Properties[VersionKey] = GetVersion();
            if (!Properties.ContainsKey(IdKey))
                Properties[IdKey] = GenerateUUID();
            if (!Properties.ContainsKey(ServicesKey))
                Properties[ServicesKey] = new Dictionary<string, object>();
        }

        /// <summary>
        /// Initialises a Service Hub.
        /// </summary>
        /// <param name="name"> The name of the service hub. </param>
        /// <param name="address"> The IP address of the service hub. </param>
        /// <param name="id"> The ID of the service hub, by default this is generated. </param>
        /// <param name="version"> The version of the service hub, by default this is automatically
        /// determined.
        /// </param>
        public ServiceHub(string name, string address, string id = null, string version = null)
        {
            Properties = new Dictionary<string, object>();
            Properties[NameKey] = name;
            Properties[AddressKey] = address;

            Properties[VersionKey] = version == null ? GetVersion() : version;
            Properties[IdKey] = id == null ? GenerateUUID() : id;
        }

        /// <summary>
        /// Determines whether this service hub definition is equivalent to another service hub.
        /// Equality is determined based on the name and address of the hub. 
        /// </summary>
        /// <param name="other"></param>
        /// <returns><c>true</c>, if the service hubs are the same, <c>false</c> otherwise.</returns>
        /// <remarks>
        /// When receiving the same service hub on multiple interfaces, they will have different addresses.
        /// The ID field can be used in this case to determine equality in origin.
        /// </remarks>
        public bool Equals(ServiceHub other)
        {
            if (ReferenceEquals(null, other)) return false;
            if (ReferenceEquals(this, other)) return true;
            return Name == other.Name && Address == other.Address && Id == other.Id;
        }

        private string GetVersion()
        {
            return GetType().Assembly.GetName().Version.ToString();
        }

        private string GenerateUUID()
        {
            return Guid.NewGuid().ToString();
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
        public override bool Equals(object obj)
        {
            if (ReferenceEquals(null, obj)) return false;
            if (ReferenceEquals(this, obj)) return true;
            if (obj.GetType() != GetType()) return false;
            return Equals((ServiceHub) obj);
        }

        /// <inheritdoc />
        public override int GetHashCode()
        {
            return Id.GetHashCode();
        }

        /// <summary>
        /// Converts this <see cref="ServiceHub"/> to its JSON representation, for transmission.
        /// </summary>
        /// <returns> JSON string representation of this <see cref="ServiceHub"/>.</returns>
        public string ToJson()
        {
            return JsonConvert.SerializeObject(Properties);
        }

        public override string ToString()
        {
            return $"{Name}:{Address}:{Id}";
        }
    }
}