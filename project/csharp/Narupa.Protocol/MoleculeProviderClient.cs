using System;
using System.Collections.Generic;
using System.Threading;
using System.Threading.Tasks;
using Grpc.Core;
using Narupa.Protocol.Instance;
using Narupa.Protocol.Trajectory;

namespace Narupa.Protocol
{
    public class MoleculeReceiver
    {
        private MoleculeProvider.MoleculeProviderClient client;

        public MoleculeReceiver(string channel)
        {
            client = new MoleculeProvider.MoleculeProviderClient(new Channel(channel, ChannelCredentials.Insecure));
            Task.Run((Action) SubscribeFrames);
        }

        private List<(uint, List<FrameData>)> frame_queue = new List<(uint, List<FrameData>)>();

        public event Action<uint, List<FrameData>> OnFrameData;

        public void Update()
        {
            foreach (var frame in frame_queue)
                OnFrameData?.Invoke(frame.Item1, frame.Item2);
            frame_queue.Clear();
        }

        private async void SubscribeFrames()
        {
            var subscription = client.SubscribeFrame(new GetFrameRequest());
            var stream = subscription.ResponseStream;
            uint frameIndex = 0;
            List<FrameData> frames = null;
            while (await stream.MoveNext(CancellationToken.None))
            {
                var response = stream.Current;
                if (response.Delimiter == Delimiter.Start)
                {
                    frames = new List<FrameData>();
                    frameIndex = response.FrameIndex;
                }
                else if (response.Delimiter == Delimiter.End)
                {
                    frame_queue.Add((frameIndex, frames));
                    frames = null;
                }
                else
                {
                    frames?.Add(response.Frame);
                }
            }
        }
    }
}