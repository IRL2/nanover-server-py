// Copyright (c) Intangible Realities Laboratory. All rights reserved.
// Licensed under the GPL. See License.txt in the project root for license information.

using Google.Protobuf.WellKnownTypes;
using JetBrains.Annotations;

namespace Narupa.Protocol.Protobuf.Extensions
{
    public static class StructExtensions
    {

        /// <summary>
        /// Gets a float value from the given <see cref="Struct" /> and key,
        /// returning null if no such key or value type exists.
        /// </summary>
        /// <param name="structure">Protobuf struct type from which to extract a value.</param>
        /// <param name="key">Field name in the struct to fetch value for.</param>
        /// <returns>Float value corresponding to field, null otherwise.</returns>
        public static float? GetFloatValue(this Struct structure, string key)
        {
            var value = GetValue(structure, key);
            return value?.KindCase == Value.KindOneofCase.NumberValue ? (float?) value.NumberValue : null;
        }

        /// <summary>
        /// Gets a double value from the given <see cref="Struct" /> and key,
        /// returning null if no such key or value type exists.
        /// </summary>
        /// <param name="structure">Protobuf struct type from which to extract a value.</param>
        /// <param name="key">Field name in the struct to fetch value for.</param>
        /// <returns>double value corresponding to field, null otherwise.</returns>
        public static double? GetDoubleValue(this Struct structure, string key)
        {
            var value = GetValue(structure, key);
            return value?.KindCase == Value.KindOneofCase.NumberValue ? (double?) value.NumberValue : null;
        }

        /// <summary>
        /// Gets an integer value from the given <see cref="Struct" /> and key,
        /// returning null if no such key or value type exists.
        /// </summary>
        /// <param name="structure">Protobuf struct type from which to extract a value.</param>
        /// <param name="key">Field name in the struct to fetch value for.</param>
        /// <returns>integer value corresponding to field, null otherwise.</returns>
        public static int? GetIntValue(this Struct structure, string key)
        {
            var value = GetValue(structure, key);
            return value?.KindCase == Value.KindOneofCase.NumberValue ? (int?) value.NumberValue : null;
        }

        /// <summary>
        /// Gets an unsigned integer value from the given <see cref="Struct" /> and key,
        /// returning null if no such key or value type exists.
        /// </summary>
        /// <param name="structure">Protobuf struct type from which to extract a value.</param>
        /// <param name="key">Field name in the struct to fetch value for.</param>
        /// <returns>unsigned integer value corresponding to field, null otherwise.</returns>
        public static uint? GetUIntValue(this Struct structure, string key)
        {
            var value = GetValue(structure, key);
            return value?.KindCase == Value.KindOneofCase.NumberValue ? (uint?) value.NumberValue : null;
        }

        /// <summary>
        /// Gets a string value from the given <see cref="Struct" /> and key,
        /// returning null if no such key or value type exists.
        /// </summary>
        /// <param name="structure">Protobuf struct type from which to extract a value.</param>
        /// <param name="key">Field name in the struct to fetch value for.</param>
        /// <returns>string value corresponding to field, null otherwise.</returns>
        public static string GetStringValue(this Struct structure, string key)
        {
            var value = GetValue(structure, key);
            return value?.KindCase == Value.KindOneofCase.StringValue ? value.StringValue : null;
        }

        /// <summary>
        /// Gets a bool value from the given <see cref="Struct" /> and key,
        /// returning null if no such key or value type exists.
        /// </summary>
        /// <param name="structure">Protobuf struct type from which to extract a value.</param>
        /// <param name="key">Field name in the struct to fetch value for.</param>
        /// <returns>bool value corresponding to field, null otherwise.</returns>
        public static bool? GetBoolValue(this Struct structure, string key)
        {
            var value = GetValue(structure, key);
            return value?.KindCase == Value.KindOneofCase.BoolValue ? (bool?) value.BoolValue : null;
        }

        /// <summary>
        /// Gets a <see cref="Struct" /> value from the given <see cref="Struct" /> and key,
        /// returning null if no such key or value type exists.
        /// </summary>
        /// <param name="structure">Protobuf struct type from which to extract a value.</param>
        /// <param name="key">Field name in the struct to fetch value for.</param>
        /// <returns><see cref="Struct" /> value corresponding to field, null otherwise.</returns>
        public static Struct GetStructValue(this Struct structure, string key)
        {
            var value = GetValue(structure, key);
            return value?.KindCase == Value.KindOneofCase.StructValue ? value.StructValue : null;
        }

        /// <summary>
        /// Gets a <see cref="ListValue" /> value from the given <see cref="Struct" /> and key,
        /// returning null if no such key or value type exists.
        /// </summary>
        /// <param name="structure">Protobuf struct type from which to extract a value.</param>
        /// <param name="key">Field name in the struct to fetch value for.</param>
        /// <returns><see cref="ListValue" /> value corresponding to field, null otherwise.</returns>
        public static ListValue GetListValue(this Struct structure, string key)
        {
            var value = GetValue(structure, key);
            return value?.KindCase == Value.KindOneofCase.ListValue ? value.ListValue : null;
        }

        /// <summary>
        /// Gets a <see cref="Value" /> from the given <see cref="Struct" /> and key,
        /// returning null if no such key or the structure itself exists.
        /// </summary>
        /// <param name="structure">Protobuf struct type from which to extract a value.</param>
        /// <param name="key">Field name in the struct to fetch value for.</param>
        /// <returns><see cref="Value" /> corresponding to field name, null otherwise.</returns>
        public static Value GetValue([NotNull] this Struct structure, string key)
        {
            var success = structure.Fields.TryGetValue(key, out var value);
            return success ? value : null;
        }

        /// <summary>
        /// Sets a boolean value in the given <see cref="Struct" /> and field, replacing any existing value.
        /// </summary>
        /// <param name="structure">Protobuf struct type in which to set a value.</param>
        /// <param name="key">Field name in the struct to set value for.</param>
        /// <param name="value">Value to be set.</param>
        public static void SetBoolValue(this Struct structure, string key, bool value)
        {
            structure.Fields[key] = Value.ForBool(value);
        }

        /// <summary>
        /// Sets a string value in the given <see cref="Struct" /> and field, replacing any existing value.
        /// </summary>
        /// <param name="structure">Protobuf struct type in which to set a value.</param>
        /// <param name="key">Field name in the struct to set value for.</param>
        /// <param name="value">Value to be set.</param>
        public static void SetStringValue(this Struct structure, string key, string value)
        {
            structure.Fields[key] = Value.ForString(value);
        }

        /// <summary>
        /// Sets a <see cref="Struct" /> value in the given <see cref="Struct" /> and field, replacing any existing value.
        /// </summary>
        /// <param name="structure">Protobuf struct type in which to set a value.</param>
        /// <param name="key">Field name in the struct to set value for.</param>
        /// <param name="value">Value to be set.</param>
        public static void SetStructValue(this Struct structure, string key, Struct value)
        {
            structure.Fields[key] = Value.ForStruct(value);
        }

        /// <summary>
        /// Sets a float value in the given <see cref="Struct" /> and field, replacing any existing value.
        /// </summary>
        /// <param name="structure">Protobuf struct type in which to set a value.</param>
        /// <param name="key">Field name in the struct to set value for.</param>
        /// <param name="value">Value to be set.</param>
        public static void SetFloatValue(this Struct structure, string key, float value)
        {
            SetDoubleValue(structure, key, value);
        }

        /// <summary>
        /// Sets an integer value in the given <see cref="Struct" /> and field, replacing any existing value.
        /// </summary>
        /// <param name="structure">Protobuf struct type in which to set a value.</param>
        /// <param name="key">Field name in the struct to set value for.</param>
        /// <param name="value">Value to be set.</param>
        public static void SetIntValue(this Struct structure, string key, int value)
        {
            SetDoubleValue(structure, key, value);
        }

        /// <summary>
        /// Sets an unsigned integer value in the given <see cref="Struct" /> and field, replacing any existing value.
        /// </summary>
        /// <param name="structure">Protobuf struct type in which to set a value.</param>
        /// <param name="key">Field name in the struct to set value for.</param>
        /// <param name="value">Value to be set.</param>
        public static void SetUIntvalue(this Struct structure, string key, uint value)
        {
            SetDoubleValue(structure, key, value);
        }

        /// <summary>
        /// Sets a double value in the given <see cref="Struct" /> and field, replacing any existing value.
        /// </summary>
        /// <param name="structure">Protobuf struct type in which to set a value.</param>
        /// <param name="key">Field name in the struct to set value for.</param>
        /// <param name="value">Value to be set.</param>
        public static void SetDoubleValue(this Struct structure, string key, double value)
        {
            structure.Fields[key] = new Value(Value.ForNumber(value));
        }

        /// <summary>
        /// Sets a <see cref="ListValue" /> in the given <see cref="Struct" /> and field, replacing any existing value.
        /// </summary>
        /// <param name="structure">Protobuf struct type in which to set a value.</param>
        /// <param name="key">Field name in the struct to set value for.</param>
        /// <param name="value">Value to be set.</param>
        public static void SetListValue(this Struct structure, string key, Value[] value)
        {
            structure.Fields[key] = Value.ForList(value);
        }

        /// <summary>
        /// Sets a <see cref="ListValue" /> in the given <see cref="Struct" /> and field, replacing any existing value.
        /// </summary>
        /// <param name="structure">Protobuf struct type in which to set a value.</param>
        /// <param name="key">Field name in the struct to set value for.</param>
        /// <param name="value">Value to be set.</param>
        public static void SetListValue(this Struct structure, string key, ListValue value)
        {
            structure.Fields[key] = new Value
            {
                ListValue = value
            };
        }
    }
}