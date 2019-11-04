using System.Net.Sockets;
using System.Text;
using System.Threading.Tasks;
using Essd;

namespace Essd.Test
{
    internal class SimpleServer
    {
        private const string broadcastIPAddress = "255.255.255.255";
        private int port;
        private UdpClient client = new UdpClient(); 
        
        public SimpleServer(int port = Client.DefaultListenPort)
        {
            this.port = port;
        }

        public async Task<int> BroadcastAsync(ServiceHub hub)
        {
            string hubJson = hub.ToJson();
            var message = Encoding.UTF8.GetBytes(hubJson);
            return await client.SendAsync(message, message.Length, broadcastIPAddress , port);
        }
        
    }
}