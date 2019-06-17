// Copyright (c) Intangible Realities Laboratory. All rights reserved.
// Licensed under the GPL. See License.txt in the project root for license information.

using Google.Protobuf.WellKnownTypes;
using Narupa.Protocol.Protobuf.Extensions;

namespace Narupa.Protocol.Imd
{
    public partial class ParticleInteraction
    {
        /// <summary>
        /// Field name in the underlying properties structure that corresponds to interaction type.
        /// </summary>
        public const string TypeKey = "type";
        /// <summary>
        /// Field name in the underlying properties structure that corresponds to scale.
        /// </summary>
        public const string ScaleKey = "scale";
        /// <summary>
        /// Field name in the underlying properties structure that corresponds to whether or not to use
        /// mass weighting.
        /// </summary>
        public const string MassWeightedKey = "mass_weighted";

        /// <summary>
        ///     Constructor for an interactive molecular dynamics interaction that provides shortcuts
        ///     for setting commonly used parameters.
        /// </summary>
        /// <param name="playerId">The player ID performing the interaction. </param>
        /// <param name="interactionId">An identifier of the interaction.</param>
        /// <param name="interactionType">The type of interaction to apply.</param>
        /// <param name="scale">The scale to be applied.</param>
        /// <param name="massWeighted">Whether mass weighting should be used.</param>
        public ParticleInteraction(string playerId, string interactionId = "0", string interactionType = "gaussian",
            float scale = 1.0f, bool massWeighted = true)
        {
            Properties = new Struct();
            PlayerId = playerId;
            InteractionId = interactionId;
            Type = interactionType;
            Scale = scale;
            MassWeighted = massWeighted;
        }

        /// <summary>
        ///     The type of interaction potential to be used with this interaction.
        /// </summary>
        /// <remarks>
        ///     Typically set to "gaussian" or "harmonic".
        /// </remarks>
        public string Type
        {
            get
            {
                EnsurePropertiesExists();
                return Properties.GetStringValue(TypeKey) ?? "gaussian";
            }
            set
            {
                EnsurePropertiesExists();
                Properties.SetStringValue(TypeKey, value);
            }
        }

        /// <summary>
        ///     The scale factor to apply to the interaction, adjusting the strength.
        /// </summary>
        public float Scale
        {
            get
            {
                EnsurePropertiesExists();
                return Properties.GetFloatValue(ScaleKey) ?? 1.0f;
            }
            set
            {
                EnsurePropertiesExists();
                Properties.SetFloatValue(ScaleKey, value);
            }
        }

        /// <summary>
        ///     Whether the interaction should be mass weighted according to the mass of the particles it is applied to.
        /// </summary>
        /// <remarks>
        ///     For classical molecular mechanics simulations, mass weighting generally provides stability.
        ///     For reactive simulations not using mass weighting can make it easier to break/form bonds.
        /// </remarks>
        public bool MassWeighted
        {
            get
            {
                EnsurePropertiesExists();
                return Properties.GetBoolValue(MassWeightedKey) ?? true;
            }
            set
            {
                EnsurePropertiesExists();
                Properties.SetBoolValue(MassWeightedKey, value);
            }
        }
        
        private void EnsurePropertiesExists()
        {
            if (Properties == null) Properties = new Struct();
        }
    }
}