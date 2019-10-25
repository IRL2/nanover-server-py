using System;
using System.Collections.Generic;
using Narupa.Protocol.Trajectory;
using NUnit.Framework;

namespace Narupa.Protocol.Test
{
    public class FrameDataTest
    {
        private static IEnumerable<AddArrayParameter>
            GetAddArrayParameters()
        {
            yield return new AddArrayParameter
            {
                Array = new[] { 0u, 1u, 2u },
                AddArray = (data, array) => data.AddIndexArray("test", (IEnumerable<uint>) array),
                GetArray = data => data.Arrays["test"].IndexValues.Values
            };
            yield return new AddArrayParameter
            {
                Array = new[] { 0f, 1f, 2f },
                AddArray = (data, array) => data.AddFloatArray("test", (IEnumerable<float>) array),
                GetArray = data => data.Arrays["test"].FloatValues.Values
            };
            yield return new AddArrayParameter
            {
                Array = new[] { "a", "b", "c" },
                AddArray =
                    (data, array) => data.AddStringArray("test", (IEnumerable<string>) array),
                GetArray = data => data.Arrays["test"].StringValues.Values
            };
        }

        private static IEnumerable<ShortcutParameter>
            ShortcutParameters()
        {
            yield return new ShortcutParameter
            {
                Array = new[] { 0f, 1f, 2f },
                Get = data => data.GetParticlePositions(),
                Set = (data, value) => data.SetParticlePositions((IEnumerable<float>) value)
            };
            yield return new ShortcutParameter
            {
                Array = new[] { "a", "b", "c" },
                Get = data => data.GetParticleTypes(),
                Set = (data, value) => data.SetParticleTypes((IEnumerable<string>) value)
            };
            yield return new ShortcutParameter
            {
                Array = new[] { "a", "b", "c" },
                Get = data => data.GetParticleNames(),
                Set = (data, value) => data.SetParticleNames((IEnumerable<string>) value)
            };
            yield return new ShortcutParameter
            {
                Array = new[] { "a", "b", "c" },
                Get = data => data.GetResidueNames(),
                Set = (data, value) => data.SetResidueNames((IEnumerable<string>) value)
            };
            yield return new ShortcutParameter
            {
                Array = new[] { "1", "2", "3"},
                Get = data => data.GetResidueIds(),
                Set = (data, value) => data.SetResidueIds((IEnumerable<string>) value)
            };
            yield return new ShortcutParameter
            {
                Array = new[] { "a", "b", "c" },
                Get = data => data.GetChainNames(),
                Set = (data, value) => data.SetChainNames((IEnumerable<string>) value)
            };
            yield return new ShortcutParameter
            {
                Array = new[] { 0u, 1u, 2u },
                Get = data => data.GetBondPairs(),
                Set = (data, value) => data.SetBondPairs((IReadOnlyList<uint>) value)
            };
            yield return new ShortcutParameter
            {
                Array = new[] { 0u, 1u, 2u },
                Get = data => data.GetParticleElements(),
                Set = (data, value) => data.SetParticleElements((IReadOnlyList<uint>) value)
            };
        }

        private static IEnumerable<TryGetArrayParameter>
            GetTryGetArrayParameters()
        {
            yield return new TryGetArrayParameter
            {
                Array = new[] { 0u, 1u, 2u },
                AddArray = (data, array) => data.AddIndexArray("test", (IEnumerable<uint>) array),
                AddArrayWrongName = (data, array)
                    => data.AddIndexArray("test2", (IEnumerable<uint>) array),
                AddArrayWrongType = data => data.AddStringArray("test", new[] { "a", "b", "c" }),
                TryGetArray = data
                    => data.TryGetIndexArray("test", out var value) ? (true, value) : (false, null),
                GetArray = data => data.GetIndexArray("test")
            };
            yield return new TryGetArrayParameter
            {
                Array = new[] { 0f, 1f, 2f },
                AddArray = (data, array) => data.AddFloatArray("test", (IEnumerable<float>) array),
                AddArrayWrongName = (data, array)
                    => data.AddFloatArray("test2", (IEnumerable<float>) array),
                AddArrayWrongType = data => data.AddStringArray("test", new[] { "a", "b", "c" }),
                TryGetArray = data
                    => data.TryGetFloatArray("test", out var value) ? (true, value) : (false, null),
                GetArray = data => data.GetFloatArray("test")
            };
            yield return new TryGetArrayParameter
            {
                Array = new[] { "a", "b", "c" },
                AddArray =
                    (data, array) => data.AddStringArray("test", (IEnumerable<string>) array),
                AddArrayWrongName = (data, array)
                    => data.AddStringArray("test2", (IEnumerable<string>) array),
                AddArrayWrongType = data => data.AddFloatArray("test", new[] { 0f, 1f, 2f }),
                TryGetArray = data
                    => data.TryGetStringArray("test", out var value)
                           ? (true, value)
                           : (false, null),
                GetArray = data => data.GetStringArray("test")
            };
        }

        [Test]
        public void AddArray(
            [ValueSource(nameof(GetAddArrayParameters))] AddArrayParameter testParameters)
        {
            var data = new FrameData();
            testParameters.AddArray(data, testParameters.Array);

            var existing = testParameters.GetArray(data);

            Assert.IsNotNull(existing);
            Assert.AreEqual(testParameters.Array, existing);
        }

        [Test]
        public void TryGetArray(
            [ValueSource(nameof(GetTryGetArrayParameters))]
            TryGetArrayParameter testParameters)
        {
            var data = new FrameData();
            testParameters.AddArray(data, testParameters.Array);

            var (success, field) = testParameters.TryGetArray(data);

            Assert.AreEqual(true, success);

            Assert.IsNotNull(field);
            Assert.AreEqual(testParameters.Array, field);
        }

        [Test]
        public void TryGetIndexValues_Missing(
            [ValueSource(nameof(GetTryGetArrayParameters))]
            TryGetArrayParameter testParameters)
        {
            var data = new FrameData();

            var (success, _) = testParameters.TryGetArray(data);

            Assert.AreEqual(false, success);
        }

        [Test]
        public void TryGetIndexValues_WrongType(
            [ValueSource(nameof(GetTryGetArrayParameters))]
            TryGetArrayParameter testParameters)
        {
            var data = new FrameData();
            testParameters.AddArrayWrongType(data);

            var (success, _) = testParameters.TryGetArray(data);

            Assert.AreEqual(false, success);
        }

        [Test]
        public void TryGetIndexValues_WrongName(
            [ValueSource(nameof(GetTryGetArrayParameters))]
            TryGetArrayParameter testParameters)
        {
            var data = new FrameData();
            testParameters.AddArrayWrongName(data, testParameters.Array);

            var (success, _) = testParameters.TryGetArray(data);

            Assert.AreEqual(false, success);
        }

        [Test]
        public void GetArray(
            [ValueSource(nameof(GetTryGetArrayParameters))]
            TryGetArrayParameter testParameters)
        {
            var data = new FrameData();
            testParameters.AddArray(data, testParameters.Array);

            var field = testParameters.GetArray(data);

            Assert.IsNotNull(field);
            Assert.AreEqual(testParameters.Array, field);
        }

        [Test]
        public void GetIndexValues_Missing(
            [ValueSource(nameof(GetTryGetArrayParameters))]
            TryGetArrayParameter testParameters)
        {
            var data = new FrameData();

            var field = testParameters.GetArray(data);

            Assert.IsNull(field);
        }

        [Test]
        public void GetIndexValues_WrongType(
            [ValueSource(nameof(GetTryGetArrayParameters))]
            TryGetArrayParameter testParameters)
        {
            var data = new FrameData();

            var field = testParameters.GetArray(data);

            Assert.IsNull(field);
        }

        [Test]
        public void GetIndexValues_WrongName(
            [ValueSource(nameof(GetTryGetArrayParameters))]
            TryGetArrayParameter testParameters)
        {
            var data = new FrameData();
            testParameters.AddArrayWrongName(data, testParameters.Array);

            var field = testParameters.GetArray(data);

            Assert.IsNull(field);
        }

        [Test]
        public void SetAndGetShortcut(
            [ValueSource(nameof(ShortcutParameters))] ShortcutParameter testParameters)
        {
            var data = new FrameData();
            testParameters.Set(data, testParameters.Array);

            var field = testParameters.Get(data);

            Assert.AreEqual(testParameters.Array, field);
        }

        [Test]
        public void GetShortcut_Missing(
            [ValueSource(nameof(ShortcutParameters))] ShortcutParameter testParameters)
        {
            var data = new FrameData();

            var field = testParameters.Get(data);

            Assert.IsNull(field);
        }

        public struct AddArrayParameter
        {
            public object Array;
            public Action<FrameData, object> AddArray;
            public Func<FrameData, object> GetArray;
        }

        public struct TryGetArrayParameter
        {
            public object Array;
            public Action<FrameData, object> AddArray;
            public Action<FrameData> AddArrayWrongType;
            public Action<FrameData, object> AddArrayWrongName;
            public Func<FrameData, (bool, object)> TryGetArray;
            public Func<FrameData, object> GetArray;
        }

        public struct ShortcutParameter
        {
            public object Array;
            public Action<FrameData, object> Set;
            public Func<FrameData, object> Get;
        }
    }
}