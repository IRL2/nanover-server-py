using System.Collections.Generic;
using System.Threading.Tasks;
using Grpc.Core;
using Narupa.Protocol.Instance;
using Narupa.Protocol.Topology;

namespace Narupa.Protocol
{
    public class Test : InstanceService.InstanceServiceBase
    {
        public override async Task StreamTopology(StreamTopologyRequest request, IServerStreamWriter<StreamTopologyResponse> responseStream, ServerCallContext context)
        {
            while (true)
            {
                await Task.Delay(1);
            }
        }

        private class Stream
        {
            public IServerStreamWriter<StreamTopologyResponse> Stream;
        }

        private List<StreamTopologyResponse> streams;
        
        public void StreamTopologyPublish(uint frameIndex, List<TopologyInfo> info)
        {
            StreamTopologySendToAll(new StreamTopologyResponse() {FrameIndex = frameIndex, Delimiter = Delimiter.Start});
            foreach (var response in responses)
            {
                
            }
            StreamTopologySendToAll(new StreamTopologyResponse() {FrameIndex = frameIndex, Delimiter = Delimiter.Start});
        }

        private void StreamTopologySendToAll(StreamTopologyResponse streamTopologyResponse)
        {
            foreach (var stream in streams)
            {
                Stream.writeasync();
            }
        }
    }
}