// Copyright (c) Intangible Realities Laboratory. All rights reserved.
// Licensed under the GPL. See License.txt in the project root for license information.

using System;
using System.Collections;
using System.Collections.Generic;
using Google.Protobuf.Collections;
using Google.Protobuf.WellKnownTypes;

namespace Narupa.Protocol.Trajectory
{
    public partial class FrameData : IEnumerable
    {
        public const string BondArrayKey = "bond.pairs";
        public const string BondOrderArrayKey = "bond.orders";

        public const string ParticlePositionArrayKey = "particle.positions";
        public const string ParticleElementArrayKey = "particle.elements";
        public const string ParticleTypeArrayKey = "particle.types";
        public const string ParticleNameArrayKey = "particle.names";
        public const string ParticleResidueArrayKey = "particle.residues";
        public const string ParticleCountValueKey = "particle.count";

        public const string ResidueNameArrayKey = "residue.names";
        public const string ResidueIdArrayKey = "residue.ids";
        public const string ResidueChainArrayKey = "residue.chains";
        public const string ResidueCountValueKey = "residue.count";
        
        public const string ChainNameArrayKey = "chain.names";
        public const string ChainCountValueKey = "chain.count";
        
        public const string KineticEnergyValueKey = "energy.kinetic";
        public const string PotentialEnergyValueKey = "energy.potential";

        public IEnumerator GetEnumerator()
        {
            throw new NotImplementedException();
        }

        /// <inheritdoc cref="IFrameData.TryGetIndexArray"/>
        public bool TryGetIndexArray(string id, out IReadOnlyList<uint> value)
        {
            if (Arrays.TryGetValue(id, out var array) &&
                array.ValuesCase == ValueArray.ValuesOneofCase.IndexValues)
            {
                value = array.IndexValues.Values;
                return true;
            }

            value = default(RepeatedField<uint>);
            return false;
        }

        /// <inheritdoc cref="IFrameData.TryGetFloatArray"/>
        public bool TryGetFloatArray(string id, out IReadOnlyList<float> value)
        {
            if (Arrays.TryGetValue(id, out var array) &&
                array.ValuesCase == ValueArray.ValuesOneofCase.FloatValues)
            {
                value = array.FloatValues.Values;
                return true;
            }

            value = default(RepeatedField<float>);
            return false;
        }

        /// <inheritdoc cref="IFrameData.TryGetStringArray"/>
        public bool TryGetStringArray(string id, out IReadOnlyList<string> value)
        {
            if (Arrays.TryGetValue(id, out var array) &&
                array.ValuesCase == ValueArray.ValuesOneofCase.StringValues)
            {
                value = array.StringValues.Values;
                return true;
            }

            value = default(RepeatedField<string>);
            return false;
        }


        /// <summary>
        /// Add a general object 'item' to FrameData. If 'item' cannot be stored, throws an
        /// InvalidOperationException
        /// </summary>
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
                case double value:
                    AddNumericValue(id, value);
                    break;
                default:
                    throw new ArgumentException($"Invalid FrameData Item with key {id}");
            }
        }

        /// <inheritdoc cref="IFrameData.AddFloatArray"/>
        public void AddFloatArray(string id, IEnumerable<float> array)
        {
            var valueArray = new ValueArray { FloatValues = new FloatArray() };
            var values = valueArray.FloatValues.Values;
            values.AddRange(array);
            Arrays.Add(id, valueArray);
        }

        /// <inheritdoc cref="IFrameData.AddIndexArray"/>
        public void AddIndexArray(string id, IEnumerable<uint> array)
        {
            var valueArray = new ValueArray { IndexValues = new IndexArray() };
            var values = valueArray.IndexValues.Values;
            values.AddRange(array);
            Arrays.Add(id, valueArray);
        }

        /// <inheritdoc cref="IFrameData.AddStringArray"/>
        public void AddStringArray(string id, IEnumerable<string> array)
        {
            var valueArray = new ValueArray { StringValues = new StringArray() };
            var values = valueArray.StringValues.Values;
            values.AddRange(array);
            Arrays.Add(id, valueArray);
        }

        /// <inheritdoc cref="IFrameData.AddNumericValue"/>
        public void AddNumericValue(string id, double value)
        {
            var valuePb = Value.ForNumber(value);
            Values.Add(id, valuePb);
        }

        /// <inheritdoc cref="IFrameData.TryGetNumericValue"/>
        public bool TryGetNumericValue(string id, out double value)
        {
            if (Values.TryGetValue(id, out var array) &&
                array.KindCase == Value.KindOneofCase.NumberValue)
            {
                value = array.NumberValue;
                return true;
            }

            value = default(double);
            return false;
        }
    }
}