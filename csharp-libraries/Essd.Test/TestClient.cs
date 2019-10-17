using System;
using System.Linq;
using System.Net.Sockets;
using System.Text;
using System.Threading;
using System.Threading.Tasks;
using NUnit.Framework;

namespace Essd.Test
{
    public class Tests
    {

        private ServiceHub testService = new ServiceHub("test","3.4.5.6");
        
        [SetUp]
        public void Setup()
        {
            
        }

        [Test]
        public async Task TestClientBlockingSearch()
        {
            var client = new Client(54544);
            var server = new SimpleServer(54544);
            var blockingSearch = Task.Run(() => client.SearchForServices(duration: 500));
            await server.BroadcastAsync(testService);
            var services = await blockingSearch;
            Assert.AreEqual(1, services.Count);
            Assert.AreEqual(services.First(), testService);
        }
        

        
        [Test]
        public async Task TestClientBlockingSearchReceiveSameService()
        {
            var client = new Client(54546);
            var server = new SimpleServer(54546);
            var blockingSearch = Task.Run(() => client.SearchForServices(duration: 500));
            for (int i = 0; i < 5; i++)
            {
                await server.BroadcastAsync(testService);
                await Task.Delay(100);
            }
                
            var services = await blockingSearch;
            Assert.AreEqual(1, services.Count);
        }
        
        [Test]
        public async Task TestClientStopNoSearch()
        {
            var client = new Client(54547);
            try
            {
                await client.StopSearch();
            }
            catch (InvalidOperationException)
            {
                Assert.Pass();
            }
            Assert.Fail();
        }

        [Test]
        public async Task TestClientStartStopSearch()
        {
            var client = new Client(54548);
            client.StartSearch();
            Assert.True(client.Searching);
            await client.StopSearch();
            Assert.False(client.Searching);
        }

        [Test]
        public async Task TestClientSearchAsync()
        {
            int port = 54549;
            var client = new Client(port);
            var server = new SimpleServer(port);
            int serviceFound =0;
            client.ServiceFound += (sender, hub) => serviceFound++;
            client.StartSearch();
            for (int i = 0; i < 5; i++)
            {
                await Task.Delay(100);
                await server.BroadcastAsync(testService);
            }
            await Task.Delay(100);
            Assert.AreEqual(5, serviceFound);
            await client.StopSearch();
        }
        
        [Test]
        public async Task TestClientInvalidString()
        {
            var client = new Client(54550);
            var server = new UdpClient();
            var blockingSearch = Task.Run(() => client.SearchForServices(duration: 500));
            byte[] invalidString = Encoding.UTF32.GetBytes("somestring");
            await server.SendAsync(invalidString, invalidString.Length, "255.255.255.255", 54550);
            var services = await blockingSearch;
            Assert.AreEqual(0, services.Count);
        }
        
        [Test]
        public async Task TestClientInvalidJson()
        {
            var client = new Client(54551);
            var server = new UdpClient();
            var blockingSearch = Task.Run(() => client.SearchForServices(duration: 500));
            byte[] invalidString = Encoding.UTF8.GetBytes("somestring");
            await server.SendAsync(invalidString, invalidString.Length, "255.255.255.255", 54550);
            var services = await blockingSearch;
            Assert.AreEqual(0, services.Count);
        }

        [Test]
        public async Task TestClientDoubleStartSearch()
        {
            var client = new Client(54552);
            client.StartSearch();
            try
            {
                client.StartSearch();
            }
            catch (InvalidOperationException)
            {
                Assert.Pass();
            }
            Assert.Fail();
        }

        [Test]
        public async Task TestClientStartBlockingWhileSearching()
        {
            var client = new Client(54553);
            client.StartSearch();
            try
            {
                client.SearchForServices();
            }
            catch (InvalidOperationException)
            {
                Assert.Pass();
            }
            Assert.Fail();
        }
    }
}