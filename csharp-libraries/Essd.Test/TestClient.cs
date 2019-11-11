using System;
using System.Collections.Generic;
using System.Linq;
using System.Net.Sockets;
using System.Text;
using System.Threading;
using System.Threading.Tasks;
using NUnit.Framework;

namespace Essd.Test
{
    internal class Tests
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
        public async Task TestMultipleClientSocketReuse()
        {
            var client = new Client(54544);
            var server = new SimpleServer(54544);
            var anotherClient = new Client(54544);

            var blockingSearch = Task.Run(() => client.SearchForServices(duration: 500));
            var anotherSearch = Task.Run((() => anotherClient.SearchForServices(500)));
            await server.BroadcastAsync(testService);
            var services = await blockingSearch;
            var secondServices = await anotherSearch;
            Assert.AreEqual(1, services.Count);
            Assert.AreEqual(services.First(), testService);
            Assert.AreEqual(services, secondServices);
            
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

        [Test]
        public void TestFilterServicesById()
        {
            var serviceB = new ServiceHub(testService.Name, "123.434.345.234", testService.Id);
            var services = new List<ServiceHub>();
            services.Add(testService);
            services.Add(serviceB);
            var filtered = Client.FilterServicesById(services);
            Assert.AreEqual(1, filtered.Count);
        }
        
        [Test]
        public void TestFilterServicesByIdLocalHost()
        {
            var serviceB = new ServiceHub(testService.Name, "127.0.0.1", testService.Id);
            var services = new List<ServiceHub>();
            services.Add(testService);
            services.Add(serviceB);
            var filtered = Client.FilterServicesById(services);
            Assert.AreEqual(1, filtered.Count);
            Assert.True(filtered.First().Address == "127.0.0.1");
        }
        
        [Test]
        public void TestFilterServicesByIdNotLocalHost()
        {
            var serviceB = new ServiceHub(testService.Name, "127.0.0.1", testService.Id);
            var services = new List<ServiceHub>();
            services.Add(testService);
            services.Add(serviceB);
            var filtered = Client.FilterServicesById(services, prioritiseLocalHost:false);
            Assert.AreEqual(1, filtered.Count);
            Assert.False(filtered.First().Address == "127.0.0.1");
        }
        
        [Test]
        public void TestFilterServicesByIdMultipleServices()
        {
            var serviceB = new ServiceHub(testService.Name, "123.434.345.234", testService.Id);
            var serviceC = new ServiceHub(testService.Name, "1.3.4.5", "5483");
            var serviceD = new ServiceHub(testService.Name, "24.4.5.8", "5483");
            var services = new List<ServiceHub>();
            services.Add(testService);
            services.Add(serviceB);
            services.Add(serviceC);
            services.Add(serviceD);
            
            var filtered = Client.FilterServicesById(services);
            Assert.AreEqual(2, filtered.Count);
        }
        
        [Test]
        public void TestFilterServicesByIdMultipleServicesLocalhost()
        {
            var serviceB = new ServiceHub(testService.Name, "123.434.345.234", testService.Id);
            var serviceC = new ServiceHub(testService.Name, "1.3.4.5", "5483");
            var serviceD = new ServiceHub(testService.Name, "127.0.0.1", "5483");
            var services = new List<ServiceHub>();
            services.Add(testService);
            services.Add(serviceB);
            services.Add(serviceC);
            services.Add(serviceD);
            
            var filtered = Client.FilterServicesById(services);
            Assert.AreEqual(2, filtered.Count);
            var localHostCount = filtered.Where(s => s.Address == "127.0.0.1").Count();
            Assert.AreEqual(1, localHostCount);
        }

        [Test]
        public async Task TestSendUTF8()
        {
            testService.Properties["name"] = testService.Name + "한국어😀";
            var client = new Client(54554);
            var server = new SimpleServer(54554);
            var blockingSearch = Task.Run(() => client.SearchForServices(duration: 500));
            await server.BroadcastAsync(testService);
            var services = await blockingSearch;
            Assert.AreEqual(1, services.Count);
            Assert.AreEqual(services.First(), testService);
        }
    }
}