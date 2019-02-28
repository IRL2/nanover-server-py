using System;
using System.Collections.Generic;
using System.Threading.Tasks;
using Grpc.Core;
using Narupa.Protocol.Renderer;

namespace Narupa.Protocol
{
    public class RendererProviderService : RendererProvider.RendererProviderBase
    {
        private readonly List<Action> actions = new List<Action>();

        public Action<RunGraphicsBenchmarkRequest, Action<RunGraphicsBenchmarkResponse>> GraphicsBenchmark { get; set; }

        public override Task<RunGraphicsBenchmarkResponse> RunGraphicsBenchmark(RunGraphicsBenchmarkRequest request,
            ServerCallContext context)
        {
            var task = new TaskCompletionSource<RunGraphicsBenchmarkResponse>();
            actions.Add(() => GraphicsBenchmark?.Invoke(request, response => task.SetResult(response)));
            return task.Task;
        }

        public void Update(float dt)
        {
            foreach (var action in actions)
                action.Invoke();
            actions.Clear();
        }
    }
}