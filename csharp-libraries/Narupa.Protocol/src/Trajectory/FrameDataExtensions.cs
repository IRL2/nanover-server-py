using System.Collections.Generic;

namespace Narupa.Protocol.Trajectory
{
    public static class FrameDataExtensions
    {
        #region Bonds

        /// <summary>
        /// Does this frame have bonds?
        /// </summary>
        public static bool HasBonds(this IFrameData data)
        {
            return data.TryGetIndexArray(FrameData.BondArrayKey, out _);
        }

        /// <summary>
        /// Set the bonds of this frame.
        /// </summary>
        public static void SetBonds(this IFrameData data,
                                    IEnumerable<uint> values)
        {
            data.AddIndexArray(FrameData.BondArrayKey, values);
        }

        /// <summary>
        /// Try to get the bonds array, returning true and setting the out
        /// variable values to the array if found.
        /// </summary>
        public static bool TryGetBonds(this IFrameData data,
                                       out IReadOnlyList<uint> values)
        {
            return data.TryGetIndexArray(FrameData.BondArrayKey, out values);
        }

        /// <summary>
        /// Get the bonds if present, else returning null.
        /// </summary>
        public static IReadOnlyList<uint> GetBonds(this IFrameData data)
        {
            return data.GetIndexArray(FrameData.BondArrayKey);
        }

        #endregion

        #region Bond Orders

        /// <summary>
        /// Does this frame have bond orders?
        /// </summary>
        public static bool HasBondOrders(this IFrameData data)
        {
            return data.TryGetFloatArray(FrameData.BondOrderArrayKey, out _);
        }

        /// <summary>
        /// Set the bond orders of this frame.
        /// </summary>
        public static void SetBondOrders(this IFrameData data,
                                         IEnumerable<float> values)
        {
            data.AddFloatArray(FrameData.BondOrderArrayKey, values);
        }

        /// <summary>
        /// Try to get the bond orders array, returning true and setting the out
        /// variable values to the array if found.
        /// </summary>
        public static bool TryGetBondOrders(this IFrameData data,
                                            out IReadOnlyList<float> values)
        {
            return data.TryGetFloatArray(FrameData.BondOrderArrayKey, out values);
        }

        /// <summary>
        /// Get the bond orders if present, else returning null.
        /// </summary>
        public static IReadOnlyList<float> GetBondOrders(this IFrameData data)
        {
            return data.GetFloatArray(FrameData.BondOrderArrayKey);
        }

        #endregion

        #region ParticlePositions

        /// <summary>
        /// Does this frame have particle positions?
        /// </summary>
        public static bool HasParticlePositions(this IFrameData data)
        {
            return data.TryGetFloatArray(FrameData.ParticlePositionArrayKey, out _);
        }

        /// <summary>
        /// Set the particle positions of this frame.
        /// </summary>
        public static void SetParticlePositions(this IFrameData data,
                                                IEnumerable<float> values)
        {
            data.AddFloatArray(FrameData.ParticlePositionArrayKey, values);
        }

        /// <summary>
        /// Try to get the particle positions array, returning true and setting the out
        /// variable values to the array if found.
        /// </summary>
        public static bool TryGetParticlePositions(this IFrameData data,
                                                   out IReadOnlyList<float> values)
        {
            return data.TryGetFloatArray(FrameData.ParticlePositionArrayKey, out values);
        }

        /// <summary>
        /// Get the particle positions if present, else returning null.
        /// </summary>
        public static IReadOnlyList<float> GetParticlePositions(this IFrameData data)
        {
            return data.GetFloatArray(FrameData.ParticlePositionArrayKey);
        }

        #endregion

        #region ParticleElements

        /// <summary>
        /// Does this frame have particle elements?
        /// </summary>
        public static bool HasParticleElements(this IFrameData data)
        {
            return data.TryGetIndexArray(FrameData.ParticleElementArrayKey, out _);
        }

        /// <summary>
        /// Set the particle elements of this frame.
        /// </summary>
        public static void SetParticleElements(this IFrameData data,
                                               IEnumerable<uint> values)
        {
            data.AddIndexArray(FrameData.ParticleElementArrayKey, values);
        }

        /// <summary>
        /// Try to get the particle elements array, returning true and setting the out
        /// variable values to the array if found.
        /// </summary>
        public static bool TryGetParticleElements(this IFrameData data,
                                                  out IReadOnlyList<uint> values)
        {
            return data.TryGetIndexArray(FrameData.ParticleElementArrayKey, out values);
        }

        /// <summary>
        /// Get the particle elements if present, else returning null.
        /// </summary>
        public static IReadOnlyList<uint> GetParticleElements(this IFrameData data)
        {
            return data.GetIndexArray(FrameData.ParticleElementArrayKey);
        }

        #endregion

        #region Particle Types

        /// <summary>
        /// Does this frame have particle types?
        /// </summary>
        public static bool HasParticleTypes(this IFrameData data)
        {
            return data.TryGetStringArray(FrameData.ParticleTypeArrayKey, out _);
        }

        /// <summary>
        /// Set the particle types of this frame.
        /// </summary>
        public static void SetParticleTypes(this IFrameData data,
                                            IEnumerable<string> values)
        {
            data.AddStringArray(FrameData.ParticleTypeArrayKey, values);
        }

        /// <summary>
        /// Try to get the particle types array, returning true and setting the out
        /// variable values to the array if found.
        /// </summary>
        public static bool TryGetParticleTypes(this IFrameData data,
                                               out IReadOnlyList<string> values)
        {
            return data.TryGetStringArray(FrameData.ParticleTypeArrayKey, out values);
        }

        /// <summary>
        /// Get the particle types if present, else returning null.
        /// </summary>
        public static IReadOnlyList<string> GetParticleTypes(this IFrameData data)
        {
            return data.GetStringArray(FrameData.ParticleTypeArrayKey);
        }

        #endregion

        #region ParticleNames

        /// <summary>
        /// Does this frame have particle names?
        /// </summary>
        public static bool HasParticleNames(this IFrameData data)
        {
            return data.TryGetStringArray(FrameData.ParticleNameArrayKey, out _);
        }

        /// <summary>
        /// Set the particle names of this frame.
        /// </summary>
        public static void SetParticleNames(this IFrameData data,
                                            IEnumerable<string> values)
        {
            data.AddStringArray(FrameData.ParticleNameArrayKey, values);
        }

        /// <summary>
        /// Try to get the particle names array, returning true and setting the out
        /// variable values to the array if found.
        /// </summary>
        public static bool TryGetParticleNames(this IFrameData data,
                                               out IReadOnlyList<string> values)
        {
            return data.TryGetStringArray(FrameData.ParticleNameArrayKey, out values);
        }

        /// <summary>
        /// Get the particle names if present, else returning null.
        /// </summary>
        public static IReadOnlyList<string> GetParticleNames(this IFrameData data)
        {
            return data.GetStringArray(FrameData.ParticleNameArrayKey);
        }

        #endregion

        #region ParticleResidues

        /// <summary>
        /// Does this frame have particle residues?
        /// </summary>
        public static bool HasParticleResidues(this IFrameData data)
        {
            return data.TryGetIndexArray(FrameData.ParticleResidueArrayKey, out _);
        }

        /// <summary>
        /// Set the particle residues of this frame.
        /// </summary>
        public static void SetParticleResidues(this IFrameData data,
                                               IEnumerable<uint> values)
        {
            data.AddIndexArray(FrameData.ParticleResidueArrayKey, values);
        }

        /// <summary>
        /// Try to get the particle residues array, returning true and setting the out
        /// variable values to the array if found.
        /// </summary>
        public static bool TryGetParticleResidues(this IFrameData data,
                                                  out IReadOnlyList<uint> values)
        {
            return data.TryGetIndexArray(FrameData.ParticleResidueArrayKey, out values);
        }

        /// <summary>
        /// Get the particle residues if present, else returning null.
        /// </summary>
        public static IReadOnlyList<uint> GetParticleResidues(this IFrameData data)
        {
            return data.GetIndexArray(FrameData.ParticleResidueArrayKey);
        }

        #endregion

        #region ResidueNames

        /// <summary>
        /// Does this frame have residue names?
        /// </summary>
        public static bool HasResidueNames(this IFrameData data)
        {
            return data.TryGetStringArray(FrameData.ResidueNameArrayKey, out _);
        }

        /// <summary>
        /// Set the residue names of this frame.
        /// </summary>
        public static void SetResidueNames(this IFrameData data,
                                           IEnumerable<string> values)
        {
            data.AddStringArray(FrameData.ResidueNameArrayKey, values);
        }

        /// <summary>
        /// Try to get the residue names array, returning true and setting the out
        /// variable values to the array if found.
        /// </summary>
        public static bool TryGetResidueNames(this IFrameData data,
                                              out IReadOnlyList<string> values)
        {
            return data.TryGetStringArray(FrameData.ResidueNameArrayKey, out values);
        }

        /// <summary>
        /// Get the residue names if present, else returning null.
        /// </summary>
        public static IReadOnlyList<string> GetResidueNames(this IFrameData data)
        {
            return data.GetStringArray(FrameData.ResidueNameArrayKey);
        }

        #endregion

        #region ResidueChains

        /// <summary>
        /// Does this frame have residue chains?
        /// </summary>
        public static bool HasResidueChains(this IFrameData data)
        {
            return data.TryGetIndexArray(FrameData.ResidueChainArrayKey, out _);
        }

        /// <summary>
        /// Set the residue chains of this frame.
        /// </summary>
        public static void SetResidueChains(this IFrameData data,
                                            IEnumerable<uint> values)
        {
            data.AddIndexArray(FrameData.ResidueChainArrayKey, values);
        }

        /// <summary>
        /// Try to get the residue chains array, returning true and setting the out
        /// variable values to the array if found.
        /// </summary>
        public static bool TryGetResidueChains(this IFrameData data,
                                               out IReadOnlyList<uint> values)
        {
            return data.TryGetIndexArray(FrameData.ResidueChainArrayKey, out values);
        }

        /// <summary>
        /// Get the residue chains if present, else returning null.
        /// </summary>
        public static IReadOnlyList<uint> GetResidueChains(this IFrameData data)
        {
            return data.GetIndexArray(FrameData.ResidueChainArrayKey);
        }

        #endregion

        #region ChainNames

        /// <summary>
        /// Does this frame have chain names?
        /// </summary>
        public static bool HasChainNames(this IFrameData data)
        {
            return data.TryGetStringArray(FrameData.ChainNameArrayKey, out _);
        }

        /// <summary>
        /// Set the chain names of this frame.
        /// </summary>
        public static void SetChainNames(this IFrameData data,
                                         IEnumerable<string> values)
        {
            data.AddStringArray(FrameData.ChainNameArrayKey, values);
        }

        /// <summary>
        /// Try to get the chain names array, returning true and setting the out
        /// variable values to the array if found.
        /// </summary>
        public static bool TryGetChainNames(this IFrameData data,
                                            out IReadOnlyList<string> values)
        {
            return data.TryGetStringArray(FrameData.ChainNameArrayKey, out values);
        }

        /// <summary>
        /// Get the chain names if present, else returning null.
        /// </summary>
        public static IReadOnlyList<string> GetChainNames(this IFrameData data)
        {
            return data.GetStringArray(FrameData.ChainNameArrayKey);
        }

        #endregion

        #region Particle Count

        /// <summary>
        /// Does this frame have a particle count?
        /// </summary>
        public static bool HasParticleCount(this IFrameData data)
        {
            return data.TryGetNumericValue(FrameData.ParticleCountValueKey, out _);
        }

        /// <summary>
        /// Set the particle count of this frame.
        /// </summary>
        public static void SetParticleCount(this IFrameData data,
                                            int value)
        {
            data.AddNumericValue(FrameData.ParticleCountValueKey, value);
        }

        /// <summary>
        /// Try to get the particle count, returning true and setting the out
        /// variable values to the value if found.
        /// </summary>
        public static bool TryGetParticleCount(this IFrameData data,
                                               out int value)
        {
            return data.TryGetIntegerValue(FrameData.ParticleCountValueKey, out value);
        }

        /// <summary>
        /// Get the particle count if present, else returning null.
        /// </summary>
        public static int? GetParticleCount(this IFrameData data)
        {
            return (int?) data.GetNumericValue(FrameData.ParticleCountValueKey);
        }

        #endregion

        #region Kinetic Energy

        /// <summary>
        /// Does this frame have a kinetic energy?
        /// </summary>
        public static bool HasKineticEnergy(this IFrameData data)
        {
            return data.TryGetNumericValue(FrameData.KineticEnergyValueKey, out _);
        }

        /// <summary>
        /// Set the kinetic energy of this frame.
        /// </summary>
        public static void SetKineticEnergy(this IFrameData data,
                                            double value)
        {
            data.AddNumericValue(FrameData.KineticEnergyValueKey, value);
        }

        /// <summary>
        /// Try to get the kinetic energy, returning true and setting the out
        /// variable values to the value if found.
        /// </summary>
        public static bool TryGetKineticEnergy(this IFrameData data,
                                               out int value)
        {
            return data.TryGetIntegerValue(FrameData.KineticEnergyValueKey, out value);
        }

        /// <summary>
        /// Get the kinetic energy if present, else returning null.
        /// </summary>
        public static float? GetKineticEnergy(this IFrameData data)
        {
            return (float?) data.GetNumericValue(FrameData.KineticEnergyValueKey);
        }

        #endregion

        #region Potential Energy

        /// <summary>
        /// Does this frame have a potential energy?
        /// </summary>
        public static bool HasPotentialEnergy(this IFrameData data)
        {
            return data.TryGetNumericValue(FrameData.PotentialEnergyValueKey, out _);
        }

        /// <summary>
        /// Set the potential energy of this frame.
        /// </summary>
        public static void SetPotentialEnergy(this IFrameData data,
                                              double value)
        {
            data.AddNumericValue(FrameData.PotentialEnergyValueKey, value);
        }

        /// <summary>
        /// Try to get the potential energy, returning true and setting the out
        /// variable values to the value if found.
        /// </summary>
        public static bool TryGetPotentialEnergy(this IFrameData data,
                                                 out int value)
        {
            return data.TryGetIntegerValue(FrameData.PotentialEnergyValueKey, out value);
        }

        /// <summary>
        /// Get the potential energy if present, else returning null.
        /// </summary>
        public static float? GetPotentialEnergy(this IFrameData data)
        {
            return (float?) data.GetNumericValue(FrameData.PotentialEnergyValueKey);
        }

        #endregion

        public static IReadOnlyList<string> GetStringArray(this IFrameData data,
                                                           string id)
        {
            return data.TryGetStringArray(id, out var array) ? array : null;
        }

        public static IReadOnlyList<uint> GetIndexArray(this IFrameData data,
                                                        string id)
        {
            return data.TryGetIndexArray(id, out var array) ? array : null;
        }

        public static IReadOnlyList<float> GetFloatArray(this IFrameData data,
                                                         string id)
        {
            return data.TryGetFloatArray(id, out var array) ? array : null;
        }

        public static double? GetNumericValue(this IFrameData data,
                                              string id)
        {
            return data.TryGetNumericValue(id, out var array) ? (double?) array : null;
        }

        public static bool TryGetIntegerValue(this IFrameData data,
                                              string id,
                                              out int output)
        {
            if (data.TryGetNumericValue(id, out var value))
            {
                output = (int) value;
                return true;
            }

            output = 0;
            return false;
        }
    }
}