using System.Collections.Generic;

namespace Narupa.Protocol.Trajectory
{
    /// <summary>
    /// General representation of a frame consisting of key-value and key-array pairs.
    /// </summary>
    public interface IFrameData
    {
        /// <summary>
        /// Try to get a float array from the frame data, keyed by the given ID, producing the value.
        /// </summary>
        /// <param name="id">The ID of the float array. </param>
        /// <param name="value">The resulting array. </param>
        /// <returns><c> True </c>, if successful, <c> False </c> otherwise. </returns>
        bool TryGetFloatArray(string id, out IReadOnlyList<float> value);

        /// <summary>
        /// Add an array of floats to the frame data, keyed by the given ID.
        /// </summary>
        /// <param name="id">The ID of the array. </param>
        /// <param name="value">The array to add to the frame data. </param>
        /// <returns><c> True </c>, if successful, <c> False </c> otherwise. </returns>
        void AddFloatArray(string id, IEnumerable<float> value);

        /// <summary>
        /// Try to get an array of indices from the frame data, keyed by the given ID, producing the value.
        /// </summary>
        /// <param name="id">The ID of the array. </param>
        /// <param name="value">The resulting array. </param>
        /// <returns><c> True </c>, if successful, <c> false </c> otherwise. </returns>
        bool TryGetIndexArray(string id, out IReadOnlyList<uint> value);
        
        /// <summary>
        /// Add an array of indices to the frame data, keyed by the given ID.
        /// </summary>
        /// <param name="id">The ID of the array. </param>
        /// <param name="value">The array to add to the frame data. </param>
        /// <returns><c> True </c>, if successful, <c> False </c> otherwise. </returns>
        void AddIndexArray(string id, IEnumerable<uint> value);

        /// <summary>
        /// Try to get an array of strings from the frame data, keyed by the given ID, producing the value.
        /// </summary>
        /// <param name="id">The ID of the array. </param>
        /// <param name="value">The resulting array. </param>
        /// <returns><c> True </c>, if successful, <c> false </c> otherwise. </returns>
        bool TryGetStringArray(string id, out IReadOnlyList<string> value);

        /// <summary>
        /// Add an array of strings to the frame data, keyed by the given ID.
        /// </summary>
        /// <param name="id">The ID of the array. </param>
        /// <param name="value">The array to add to the frame data. </param>
        /// <returns><c> True </c>, if successful, <c> False </c> otherwise. </returns>
        void AddStringArray(string id, IEnumerable<string> value);

        /// <summary>
        /// Try to get a of numeric value from the frame data, keyed by the given ID, producing the value.
        /// </summary>
        /// <param name="id">The ID of the array. </param>
        /// <param name="value">The resulting numeric value. </param>
        /// <returns><c> True </c>, if successful, <c> false </c> otherwise. </returns>
        bool TryGetNumericValue(string id, out double value);

        /// <summary>
        /// Add a of numeric value to the frame data, keyed by the given ID.
        /// </summary>
        /// <param name="id">The ID of the array. </param>
        /// <param name="value">The value to add to the frame data. </param>
        /// <returns><c> True </c>, if successful, <c> False </c> otherwise. </returns>
        void AddNumericValue(string id, double value);
    }
}