// Copyright (c) Intangible Realities Laboratory. All rights reserved.
// Licensed under the GPL. See License.txt in the project root for license information.

using System;
using System.Collections;
using System.Collections.Generic;
using Google.Protobuf.Collections;

namespace Narupa.Protocol.Trajectory
{
    public partial class FrameData : IEnumerable
    {
        public const string BondArrayKey = "bond";
        public const string ParticleElementArrayKey = "particle.element";
        public const string ParticleTypeArrayKey = "particle.type";
        public const string ParticlePositionArrayKey = "particle.position";

        /// <summary>
        ///     Float array of particle positions. Nominally grouped in sets of three to form 3D vectors
        /// </summary>
        public IReadOnlyList<float> ParticlePositions
        {
            get => GetFloatArray(ParticlePositionArrayKey);
            set => AddFloatArray(ParticlePositionArrayKey, value);
        }

        /// <summary>
        ///     Array of indices representing bonds between particles
        /// </summary>
        public IReadOnlyList<uint> Bonds
        {
            get => GetIndexArray(BondArrayKey);
            set => AddIndexArray(BondArrayKey, value);
        }

        /// <summary>
        ///     Array of indices representing particles that are atomic elements. Use 0 to indicate an particle that is not
        ///     an atom
        /// </summary>
        public IReadOnlyList<uint> ParticleElements
        {
            get => GetIndexArray(ParticleElementArrayKey);
            set => AddIndexArray(ParticleElementArrayKey, value);
        }

        /// <summary>
        ///     Array of string representing atom names
        /// </summary>
        public IReadOnlyList<string> ParticleTypes
        {
            get => GetStringArray(ParticleTypeArrayKey);
            set => AddStringArray(ParticleTypeArrayKey, value);
        }

        public IEnumerator GetEnumerator()
        {
            throw new NotImplementedException();
        }

        /// <summary>
        ///     Get an array of uint[] values stored with the key 'id' and assigns it to value, returning true if the
        ///     key existed and false otherwise
        /// </summary>
        /// <param name="id"></param>
        /// <param name="value"></param>
        /// <returns></returns>
        public bool TryGetIndexArray(string id, out RepeatedField<uint> value)
        {
            if (Arrays.TryGetValue(id, out var array) && array.ValuesCase == ValueArray.ValuesOneofCase.IndexValues)
            {
                value = array.IndexValues.Values;
                return true;
            }

            value = default(RepeatedField<uint>);
            return false;
        }

        /// <summary>
        ///     Get an array of float[] values stored with the key 'id' and assigns it to value, returning true if the
        ///     key existed and false otherwise
        /// </summary>
        /// <param name="id"></param>
        /// <param name="value"></param>
        /// <returns></returns>
        public bool TryGetFloatArray(string id, out RepeatedField<float> value)
        {
            if (Arrays.TryGetValue(id, out var array) && array.ValuesCase == ValueArray.ValuesOneofCase.FloatValues)
            {
                value = array.FloatValues.Values;
                return true;
            }

            value = default(RepeatedField<float>);
            return false;
        }

        /// <summary>
        ///     Get an array of string[] values stored with the key 'id' and assigns it to value, returning true if the
        ///     key existed and false otherwise
        /// </summary>
        /// <param name="id"></param>
        /// <param name="value"></param>
        /// <returns></returns>
        public bool TryGetStringArray(string id, out RepeatedField<string> value)
        {
            if (Arrays.TryGetValue(id, out var array) && array.ValuesCase == ValueArray.ValuesOneofCase.StringValues)
            {
                value = array.StringValues.Values;
                return true;
            }

            value = default(RepeatedField<string>);
            return false;
        }

        /// <summary>
        ///     Get an array of uint[] values stored with the key 'id', returning null if it does not exist
        /// </summary>
        /// <param name="id"></param>
        /// <returns></returns>
        public RepeatedField<uint> GetIndexArray(string id)
        {
            return TryGetIndexArray(id, out var value) ? value : null;
        }

        /// <summary>
        ///     Get an array of float[] values stored with the key 'id', returning null if it does not exist
        /// </summary>
        /// <param name="id"></param>
        /// <returns></returns>
        public RepeatedField<float> GetFloatArray(string id)
        {
            return TryGetFloatArray(id, out var value) ? value : null;
        }

        /// <summary>
        ///     Get an array of string[] values stored with the key 'id', returning null if it does not exist
        /// </summary>
        /// <param name="id"></param>
        /// <returns></returns>
        public RepeatedField<string> GetStringArray(string id)
        {
            return TryGetStringArray(id, out var value) ? value : null;
        }

        /// <summary>
        ///     Add a general object 'item' to FrameData. If 'item' cannot be stored, throws an InvalidOperationException
        /// </summary>
        /// <param name="id"></param>
        /// <param name="item"></param>
        /// <exception cref="InvalidOperationException"></exception>
        public void Add(string id, object item)
        {
            switch (item)
            {
                case float[] floatArray:
                    AddFloatArray(id, floatArray);
                    break;
                case uint[] uintArray:
                    AddIndexArray(id, uintArray);
                    break;
                case string[] stringArray:
                    AddStringArray(id, stringArray);
                    break;
                default:
                    throw new ArgumentException($"Invalid FrameData Item with key {id}");
            }
        }

        /// <summary>
        ///     Add an array of float[] values, storing with the key 'id'
        /// </summary>
        /// <param name="id"></param>
        /// <param name="array"></param>
        public void AddFloatArray(string id, IEnumerable<float> array)
        {
            var valueArray = new ValueArray {FloatValues = new FloatArray()};
            var values = valueArray.FloatValues.Values;
            values.AddRange(array);
            Arrays.Add(id, valueArray);
        }

        /// <summary>
        ///     Add an array of uint[] values, storing with the key 'id'
        /// </summary>
        /// <param name="id"></param>
        /// <param name="array"></param>
        public void AddIndexArray(string id, IEnumerable<uint> array)
        {
            var valueArray = new ValueArray {IndexValues = new IndexArray()};
            var values = valueArray.IndexValues.Values;
            values.AddRange(array);
            Arrays.Add(id, valueArray);
        }

        /// <summary>
        ///     Add an array of string[] values, storing with the key 'id'
        /// </summary>
        /// <param name="id"></param>
        /// <param name="array"></param>
        public void AddStringArray(string id, IEnumerable<string> array)
        {
            var valueArray = new ValueArray {StringValues = new StringArray()};
            var values = valueArray.StringValues.Values;
            values.AddRange(array);
            Arrays.Add(id, valueArray);
        }
    }
}