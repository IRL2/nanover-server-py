using System.Collections.Generic;
using Narupa.Protocol.Imd;
using NUnit.Framework;

namespace Narupa.Protocol.Test.Imd
{
    public class ParticleInteractionTests
    {
        /// <summary>
        /// Generates an interaction with fields set to default values, using both constructors.
        /// </summary>
        private static IEnumerable<TestCaseData>
            DefaultInteraction()
        {
            yield return new TestCaseData(new ParticleInteraction("1"));
            yield return new TestCaseData(new ParticleInteraction());
        }

        /// <summary>
        /// Generates an interaction with some different values et.
        /// </summary>
        private static IEnumerable<TestCaseData>
            TestInteraction()
        {
            yield return new TestCaseData(new ParticleInteraction("2", "1", "harmonic",
                1000f, false, true));
        }

        private static IEnumerable<TestCaseData>
            DefaultConstructor()
        {
            yield return new TestCaseData(new ParticleInteraction());
        }

        [Test]
        [TestCaseSource(nameof(DefaultInteraction))]
        public void TestDefaultType(ParticleInteraction interaction)
        {
            Assert.AreEqual("gaussian", interaction.Type);
        }

        [Test]
        [TestCaseSource(nameof(DefaultInteraction))]
        public void TestDefaultScale(ParticleInteraction interaction)
        {
            Assert.AreEqual(1f, interaction.Scale);
        }

        [Test]
        [TestCaseSource(nameof(DefaultInteraction))]
        public void TestDefaultMassWeighted(ParticleInteraction interaction)
        {
            Assert.AreEqual(true, interaction.MassWeighted);
        }

        [Test]
        [TestCaseSource(nameof(DefaultInteraction))]
        public void TestDefaultResetVels(ParticleInteraction interaction)
        {
            Assert.AreEqual(false, interaction.ResetVelocities);
        }
        
        [Test]
        [TestCaseSource(nameof(TestInteraction))]
        public void TestConstructor_PlayerId(ParticleInteraction testInteraction)
        {
            Assert.AreEqual("2", testInteraction.PlayerId);
        }

        [Test]
        [TestCaseSource(nameof(TestInteraction))]
        public void TestConstructor_Type(ParticleInteraction testInteraction)
        {
            Assert.AreEqual("harmonic", testInteraction.Type);
        }

        [Test]
        [TestCaseSource(nameof(TestInteraction))]
        public void TestConstructor_Interact(ParticleInteraction testInteraction)
        {
            Assert.AreEqual("1", testInteraction.InteractionId);
        }

        [Test]
        [TestCaseSource(nameof(TestInteraction))]
        public void TestConstructor_Scale(ParticleInteraction testInteraction)
        {
            Assert.That(testInteraction.Scale, Is.EqualTo(1000f).Within(1e-8f));
        }

        [Test]
        [TestCaseSource(nameof(TestInteraction))]
        public void TestConstructor_MassWeighted(ParticleInteraction testInteraction)
        {
            Assert.AreEqual(false, testInteraction.MassWeighted);
        }

        [Test]
        [TestCaseSource(nameof(TestInteraction))]
        public void TestConstructor_ResetVelocities(ParticleInteraction testInteraction)
        {
            Assert.AreEqual(true, testInteraction.ResetVelocities);
        }

        [Test]
        [TestCaseSource(nameof(DefaultInteraction))]
        public void TestSetInteractionType(ParticleInteraction interaction)
        {
            var targetType = "harmonic";
            interaction.Type = targetType;

            Assert.AreEqual(targetType, interaction.Type);
            Assert.AreEqual(targetType, interaction.Properties.Fields[ParticleInteraction.TypeKey].StringValue);
        }

        [Test]
        [TestCaseSource(nameof(DefaultInteraction))]
        public void TestSetScale(ParticleInteraction interaction)
        {
            var scale = 1000f;
            interaction.Scale = scale;

            Assert.AreEqual(scale, interaction.Scale);
            Assert.That(scale, Is.EqualTo(interaction.Scale).Within(1e-8f));
            Assert.That((float) interaction.Properties.Fields[ParticleInteraction.ScaleKey].NumberValue,
                Is.EqualTo(scale).Within(1e-8f));
        }

        [Test]
        [TestCaseSource(nameof(DefaultInteraction))]
        public void TestSetMassWeighted(ParticleInteraction interaction)
        {
            var massWeighted = false;
            interaction.MassWeighted = massWeighted;

            Assert.AreEqual(massWeighted, interaction.MassWeighted);
            Assert.AreEqual(massWeighted, interaction.Properties.Fields[ParticleInteraction.MassWeightedKey].BoolValue);
        }
        
        [Test]
        [TestCaseSource(nameof(DefaultInteraction))]
        public void TestSetResetVels(ParticleInteraction interaction)
        {
            var value = true;
            interaction.ResetVelocities = value;

            Assert.AreEqual(value, interaction.ResetVelocities);
            Assert.AreEqual(value, interaction.Properties.Fields[ParticleInteraction.ResetVelocitiesKey].BoolValue);
        }
    }
}