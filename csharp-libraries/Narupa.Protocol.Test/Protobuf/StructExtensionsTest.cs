// Copyright (c) Intangible Realities Laboratory. All rights reserved.
// Licensed under the GPL. See License.txt in the project root for license information.

using System;
using System.Collections.Generic;
using Google.Protobuf.WellKnownTypes;
using Narupa.Protocol.Protobuf.Extensions;
using NUnit.Framework;

namespace Narupa.Protocol.Test.Protobuf
{
    public class StructExtensionsTest
    {
        /// <summary>
        /// Set of simple structs for testing getting and setting values.
        /// </summary>
        private static IEnumerable<SetValueParameter>
            GetSetValueParameters()
        {
            yield return new SetValueParameter()
            {
                Value = 5.0f,
                SetValue = (structure, value) => structure.SetFloatValue("test", (float) value),
                
                GetValue = (structure) => structure.GetFloatValue("test"),
                GetWrongValueType = (structure) => structure.GetStringValue("test"),
            };
            yield return new SetValueParameter()
            {
                Value = 1,
                SetValue = (structure, value) => structure.SetIntValue("test", (int) value),
                GetValue = (structure) => structure.GetIntValue("test"),
                GetWrongValueType = (structure) => structure.GetStructValue("test"),
            };
            yield return new SetValueParameter()
            {
                Value = 5.0,
                SetValue = (structure, value) => structure.SetDoubleValue("test", (double) value),
                GetValue = (structure) => structure.GetDoubleValue("test"),
                GetWrongValueType = (structure) => structure.GetListValue("test"),

            };
            yield return new SetValueParameter()
            {
                Value = (uint) 5,
                SetValue = (structure, value) => structure.SetUIntvalue("test", (uint) value),
                GetValue = (structure) => structure.GetUIntValue("test"),
                GetWrongValueType = (structure) => structure.GetListValue("test"),
            };
            yield return new SetValueParameter()
            {
                Value = true,
                SetValue = (structure, value) => structure.SetBoolValue("test", (bool) value),
                GetValue = (structure) => structure.GetBoolValue("test"),
                GetWrongValueType = (structure) => structure.GetIntValue("test"),

            };
            yield return new SetValueParameter()
            {
                Value = false,
                SetValue = (structure, value) => structure.SetBoolValue("test", (bool) value),
                GetValue = (structure) => structure.GetBoolValue("test"),
                GetWrongValueType = (structure) => structure.GetStructValue("test"),

            };
            yield return new SetValueParameter()
            {
                Value = "testing",
                SetValue = (structure, value) => structure.SetStringValue("test", (string) value),
                GetValue = (structure) => structure.GetStringValue("test"),
                GetWrongValueType = (structure) => structure.GetDoubleValue("test"),

            };
            yield return new SetValueParameter()
            {
                Value = new Struct(),
                SetValue = (structure, value) => structure.SetStructValue("test", (Struct) value),
                GetValue = (structure) => structure.GetStructValue("test"),
                GetWrongValueType = (structure) => structure.GetUIntValue("test"),
            };
            
            var value1 = Value.ForNumber(5);
            var value2 = Value.ForString("test");
            yield return new SetValueParameter()
            {
                Value = new ListValue()
                {
                    Values = { value1, value2, value2}
                },
                SetValue = (structure, value) => structure.SetListValue("test", (ListValue) value),
                GetValue = (structure) => structure.GetListValue("test"),
                GetWrongValueType = (structure) => structure.GetBoolValue("test"),

            };
        }
        
        
        /// <summary>
        /// Tests that setting and getting values from a struct works.
        /// </summary>
        /// <param name="testParameters"></param>
        [Test]
        public void TestSetValue(
            [ValueSource(nameof(GetSetValueParameters))]
            SetValueParameter testParameters)
        {
            var structure = new Struct();
            testParameters.SetValue(structure, testParameters.Value);

            var existing = testParameters.GetValue(structure);

            Assert.IsNotNull(existing);
            Assert.AreEqual(testParameters.Value, existing);
        }
        
        
        /// <summary>
        /// Tests that getting the wrong value type from a struct returns null.
        /// </summary>
        /// <param name="testParameters"></param>
        [Test]
        public void TestGetWrongValueType(
            [ValueSource(nameof(GetSetValueParameters))]
            SetValueParameter testParameters)
        {
            var structure = new Struct();
            testParameters.SetValue(structure, testParameters.Value);

            var existing = testParameters.GetWrongValueType(structure);

            Assert.IsNull(existing);
        }
        
        public struct SetValueParameter
        {
            public object Value;
            public Action<Struct, object> SetValue;
            public Func<Struct, object> GetValue;
            /// <summary>
            /// A test function for getting the wrong value type, should return null.
            /// </summary>
            public Func<Struct, object> GetWrongValueType;
        }
    }
}