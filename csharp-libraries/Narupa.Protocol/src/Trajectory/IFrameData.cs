using System.Collections.Generic;

namespace Narupa.Protocol.Trajectory
{
    /// <summary>
    /// General representation of a frame consisting of key-value and key-array pairs.
    /// </summary>
    public interface IFrameData
    {
        bool TryGetFloatArray(string id, out IReadOnlyList<float> value);

        void AddFloatArray(string id, IEnumerable<float> value);

        bool TryGetIndexArray(string id, out IReadOnlyList<uint> value);

        void AddIndexArray(string id, IEnumerable<uint> value);

        bool TryGetStringArray(string id, out IReadOnlyList<string> value);

        void AddStringArray(string id, IEnumerable<string> value);

        bool TryGetNumericValue(string id, out double value);

        void AddNumericValue(string id, double value);
    }
}