using System.Collections.Generic;
using Narupa.Protocol.Imd;
using NUnit.Framework;

namespace Narupa.Protocol.Test.Imd
{
    public class InteractionTests
    {
        /// <summary>
        ///     Generates an interaction with fields set to default values, using both constructors.
        /// </summary>
        private static IEnumerable<TestCaseData>
            DefaultInteraction()
        {
            yield return new TestCaseData(new Interaction("1"));
            yield return new TestCaseData(new Interaction());
        }

        /// <summary>
        ///     Generates an interaction with some different values et.
        /// </summary>
        private static IEnumerable<TestCaseData>
            TestInteraction()
        {
            yield return new TestCaseData(new Interaction("2", "1", "harmonic",
                1000f, false));
        }

        private static IEnumerable<TestCaseData>
            DefaultConstructor()
        {
            yield return new TestCaseData(new Interaction());
        }

        [Test]
        [TestCaseSource(nameof(DefaultInteraction))]
        public void TestDefaultType(Interaction interaction)
        {
            Assert.AreEqual("gaussian", interaction.Type);
        }

        [Test]
        [TestCaseSource(nameof(DefaultInteraction))]
        public void TestDefaultScale(Interaction interaction)
        {
            Assert.AreEqual(1f, interaction.Scale);
        }

        [Test]
        [TestCaseSource(nameof(DefaultInteraction))]
        public void TestDefaultMassWeighted(Interaction interaction)
        {
            Assert.AreEqual(true, interaction.MassWeighted);
        }

        [Test]
        [TestCaseSource(nameof(TestInteraction))]
        public void TestConstructor_PlayerId(Interaction testInteraction)
        {
            Assert.AreEqual("2", testInteraction.PlayerId);
        }

        [Test]
        [TestCaseSource(nameof(TestInteraction))]
        public void TestConstructor_Type(Interaction testInteraction)
        {
            Assert.AreEqual("harmonic", testInteraction.Type);
        }

        [Test]
        [TestCaseSource(nameof(TestInteraction))]
        public void TestConstructor_Interact(Interaction testInteraction)
        {
            Assert.AreEqual("1", testInteraction.InteractionId);
        }

        [Test]
        [TestCaseSource(nameof(TestInteraction))]
        public void TestConstructor_Scale(Interaction testInteraction)
        {
            Assert.That(testInteraction.Scale, Is.EqualTo(1000f).Within(1e-8f));
        }

        [Test]
        [TestCaseSource(nameof(TestInteraction))]
        public void TestConstructor_MassWeighted(Interaction testInteraction)
        {
            Assert.AreEqual(false, testInteraction.MassWeighted);
        }


        [Test]
        [TestCaseSource(nameof(DefaultInteraction))]
        public void TestSetInteractionType(Interaction interaction)
        {
            var targetType = "harmonic";
            interaction.Type = targetType;

            Assert.AreEqual(targetType, interaction.Type);
            Assert.AreEqual(targetType, interaction.Properties.Fields[Interaction.TypeKey].StringValue);
        }

        [Test]
        [TestCaseSource(nameof(DefaultInteraction))]
        public void TestSetScale(Interaction interaction)
        {
            var scale = 1000f;
            interaction.Scale = scale;

            Assert.AreEqual(scale, interaction.Scale);
            Assert.That(scale, Is.EqualTo(interaction.Scale).Within(1e-8f));
            Assert.That((float) interaction.Properties.Fields[Interaction.ScaleKey].NumberValue,
                Is.EqualTo(scale).Within(1e-8f));
        }

        [Test]
        [TestCaseSource(nameof(DefaultInteraction))]
        public void TestSetMassWeighted(Interaction interaction)
        {
            var massWeighted = false;
            interaction.MassWeighted = massWeighted;

            Assert.AreEqual(massWeighted, interaction.MassWeighted);
            Assert.AreEqual(massWeighted, interaction.Properties.Fields[Interaction.MassWeightedKey].BoolValue);
        }
    }
}