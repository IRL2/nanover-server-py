using System;
using System.Collections.Generic;
using System.Linq;
using System.Net;
using System.Net.Sockets;
using System.Text;
using System.Threading;
using System.Threading.Tasks;
using System.Timers;
using Newtonsoft.Json;
using Timer = System.Timers.Timer;

namespace Essd
{
    // https://stackoverflow.com/questions/19404199/how-to-to-make-udpclient-receiveasync-cancelable
    public static class AsyncExtensions
    {
        /// <summary>
        /// Enables a task to be cancellable.
        /// </summary>
        /// <param name="task">Task that usually cannot be canceled.</param>
        /// <param name="cancellationToken">Cancellation token to use for cancellation.</param>
        /// <typeparam name="T">Task result type.</typeparam>
        /// <returns>Task</returns>
        /// <exception cref="OperationCanceledException">If the task is cancelled before completion.</exception>
        public static async Task<T> WithCancellation<T>(this Task<T> task, CancellationToken cancellationToken)
        {
            var tcs = new TaskCompletionSource<bool>();
            using (cancellationToken.Register(s => ((TaskCompletionSource<bool>) s).TrySetResult(true), tcs))
            {
                if (task != await Task.WhenAny(task, tcs.Task)) throw new OperationCanceledException(cancellationToken);
            }

            return task.Result;
        }
    }

    /// <summary>
    /// Implementation of an Extremely Simple Service Discovery (ESSD) client.
    /// Provides methods for searching for discoverable services.
    /// </summary>
    public class Client
    {
        /// <summary>
        /// Default port upon which ESSD clients listen for services.
        /// </summary>
        public const int DefaultListenPort = 54545;

        private readonly CancellationTokenSource cancellationTokenSource = new CancellationTokenSource();
        private Task searchTask;

        private readonly UdpClient udpClient;

        /// <summary>
        /// Initialises an ESSD client.
        /// </summary>
        /// <param name="listenPort">Port at which to listen for services.</param>
        public Client(int listenPort = DefaultListenPort)
        {
            udpClient = new UdpClient();
            udpClient.Client.SetSocketOption(SocketOptionLevel.Socket, SocketOptionName.ReuseAddress, true);

            udpClient.Client.Bind(new IPEndPoint(IPAddress.Any, listenPort));
        }

        /// <summary>
        /// Whether this client is currently searching for services.
        /// </summary>
        public bool Searching { get; private set; }

        /// <summary>
        /// Event triggered when a service is found.
        /// </summary>
        /// <remarks>
        /// Note that this event is only triggered when using the background search, and it will be
        /// triggered every time a service is found, even if that service has been previously encountered.
        /// </remarks>
        public event EventHandler<ServiceHub> ServiceFound;

        /// <summary>
        /// Start searching for discoverable services in the background.
        /// </summary>
        /// <returns> Task representing the search. </returns>
        /// <exception cref="InvalidOperationException"> If a search is already happening. </exception>
        public Task StartSearch()
        {
            if (Searching)
                throw new InvalidOperationException("Already searching!");
            searchTask = SearchForServicesAsync(cancellationTokenSource.Token);
            return searchTask;
        }


        /// <summary>
        /// Stops searching for discoverable services, if a search was underway.
        /// </summary>
        /// <returns> The task representing the search, which will have been instructed to terminate. </returns>
        /// <exception
        /// cref="InvalidOperationException">
        /// If an attempt to stop a non-existent search is made.
        /// </exception>
        public async Task StopSearch()
        {
            if (!Searching)
                throw new InvalidOperationException("Attempted to stop a non-existent search for services");
            cancellationTokenSource.Cancel();
            await searchTask;
        }

        /// <summary>
        /// Method to filter a collection of service hubs to select only one service hub per ID. Useful when
        /// a client receives the same service hub from multiple interfaces with different addresses. 
        /// </summary>
        /// <param name="services"> The collection of service hubs to filter.</param>
        /// <param name="prioritiseLocalHost">
        /// Whether to prioritise the localhost IP address, if available,
        /// in selecting from services with the same ID.
        /// </param>
        /// <returns>
        /// A collection of service hubs such that there is only one service hub for each ID
        /// found in the original collection.
        /// </returns>
        public static ICollection<ServiceHub> FilterServicesById(ICollection<ServiceHub> services, bool prioritiseLocalHost=true)
        {
            var idToServiceMap = GetIdToServiceMap(services);
            var filteredServices = new HashSet<ServiceHub>();
            foreach (var id in idToServiceMap)
            {
                ServiceHub selectedService = null;
                if (prioritiseLocalHost)
                    selectedService = GetLocalHostService(id.Value);
                if(selectedService == null)
                    selectedService = id.Value.First();
                filteredServices.Add(selectedService);
            }

            return filteredServices;
        }

        private static ServiceHub GetLocalHostService(IEnumerable<ServiceHub> services)
        {
            var selectedServices = services.Where(s => s.Address == "127.0.0.1");
            if (!selectedServices.Any())
                return null;
            return selectedServices.First();
        }

        private static IDictionary<string, HashSet<ServiceHub>> GetIdToServiceMap(ICollection<ServiceHub> services)
        {
            var idToServiceMap = new Dictionary<string, HashSet<ServiceHub>>();
            foreach (var service in services)
            {
                if(idToServiceMap.ContainsKey(service.Id) == false)
                    idToServiceMap[service.Id] = new HashSet<ServiceHub>();
                idToServiceMap[service.Id].Add(service);
            }

            return idToServiceMap;
        }

        /// <summary>
        /// Searches for services for the given duration, blocking while searching.
        /// </summary>
        /// <param name="duration">Duration to search for, in milliseconds.</param>
        /// <returns>Collection of unique services found during search.</returns>
        /// <remarks>
        /// There is no guarantee that services found during the search will still
        /// be up after the search ends.
        /// </remarks>
        public ICollection<ServiceHub> SearchForServices(int duration = 3000)
        {
            if (Searching)
                throw new InvalidOperationException(
                    "Cannot start a blocking search while running another search in the background");

            var servicesFound = new HashSet<ServiceHub>();

            var from = new IPEndPoint(0, 0);

            var timer = new Timer(duration);
            udpClient.Client.ReceiveTimeout = duration;
            timer.Elapsed += OnTimerElapsed;
            var timerElapsed = false;
            timer.Start();

            while (!timerElapsed)
            {
                var (timedOut, messageBytes) = WaitForMessage(ref from);
                if (timedOut)
                    break;

                ServiceHub service;
                try
                {
                    service = DecodeServiceHub(messageBytes);
                }
                catch (ArgumentException e)
                {
                    Console.WriteLine($"ESSD: Exception passing service definition: {e.Message}");
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


        private async Task SearchForServicesAsync(CancellationToken token)
        {
            Searching = true;
            while (!token.IsCancellationRequested)
            {
                UdpReceiveResult message;
                try
                {
                    message = await udpClient.ReceiveAsync().WithCancellation(token);
                }
                catch (OperationCanceledException)
                {
                    break;
                }

                ServiceHub service;
                try
                {
                    service = DecodeServiceHub(message.Buffer);
                }
                catch (ArgumentException e)
                {
                    Console.WriteLine($"ESSD: Exception passing service definition: {e.Message}");
                    continue;
                }

                ServiceFound?.Invoke(this, service);
            }

            Searching = false;
        }

        private ServiceHub DecodeServiceHub(byte[] messageBytes)
        {
            var message = DecodeMessage(messageBytes);
            ServiceHub service;
            try
            {
                service = new ServiceHub(message);
            }
            catch (JsonException e)
            {
                throw new ArgumentException("Invalid JSON string encountered.");
            }

            return service;
        }

        private (bool, byte[]) WaitForMessage(ref IPEndPoint from)
        {
            byte[] messageBytes;
            try
            {
                messageBytes = udpClient.Receive(ref from);
            }
            catch (SocketException e)
            {
                if (e.SocketErrorCode == SocketError.TimedOut)
                    return (true, null);
                throw;
            }

            return (false, messageBytes);
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
                throw new ArgumentException("ESSD: Received invalid message, not a valid UTF8 string.");
            }
        }
    }
}