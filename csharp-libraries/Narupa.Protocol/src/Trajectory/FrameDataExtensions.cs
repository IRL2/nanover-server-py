using System.Collections.Generic;

namespace Narupa.Protocol.Trajectory
{
    public static class FrameDataExtensions
    {
        #region Bonds

        /// <summary>
        /// Does this frame have bonds?
        /// </summary>
        public static bool HasBonds(this FrameData data)
        {
            return data.TryGetIndexArray(FrameData.BondArrayKey, out _);
        }

        /// <summary>
        /// Set the bonds of this frame.
        /// </summary>
        public static void SetBondPairs(this FrameData data,
                                    IEnumerable<uint> values)
        {
            data.AddIndexArray(FrameData.BondArrayKey, values);
        }

        /// <summary>
        /// Try to get the bond indices array, returning true and setting the out
        /// variable values to the array if found.
        /// </summary>
        public static bool TryGetBondPairs(this FrameData data,
                                       out IReadOnlyList<uint> values)
        {
            return data.TryGetIndexArray(FrameData.BondArrayKey, out values);
        }

        /// <summary>
        /// Get the bond indices if present, else returning null.
        /// </summary>
        public static IReadOnlyList<uint> GetBondPairs(this FrameData data)
        {
            return data.GetIndexArray(FrameData.BondArrayKey);
        }

        #endregion

        #region Bond Orders

        /// <summary>
        /// Does this frame have bond orders?
        /// </summary>
        public static bool HasBondOrders(this FrameData data)
        {
            return data.TryGetFloatArray(FrameData.BondOrderArrayKey, out _);
        }

        /// <summary>
        /// Set the bond orders of this frame.
        /// </summary>
        public static void SetBondOrders(this FrameData data,
                                         IEnumerable<float> values)
        {
            data.AddFloatArray(FrameData.BondOrderArrayKey, values);
        }

        /// <summary>
        /// Try to get the bond orders array, returning true and setting the out
        /// variable values to the array if found.
        /// </summary>
        public static bool TryGetBondOrders(this FrameData data,
                                            out IReadOnlyList<float> values)
        {
            return data.TryGetFloatArray(FrameData.BondOrderArrayKey, out values);
        }

        /// <summary>
        /// Get the bond orders if present, else returning null.
        /// </summary>
        public static IReadOnlyList<float> GetBondOrders(this FrameData data)
        {
            return data.GetFloatArray(FrameData.BondOrderArrayKey);
        }

        #endregion

        #region ParticlePositions

        /// <summary>
        /// Does this frame have particle positions?
        /// </summary>
        public static bool HasParticlePositions(this FrameData data)
        {
            return data.TryGetFloatArray(FrameData.ParticlePositionArrayKey, out _);
        }

        /// <summary>
        /// Set the particle positions of this frame.
        /// </summary>
        public static void SetParticlePositions(this FrameData data,
                                                IEnumerable<float> values)
        {
            data.AddFloatArray(FrameData.ParticlePositionArrayKey, values);
        }

        /// <summary>
        /// Try to get the particle positions array, returning true and setting the out
        /// variable values to the array if found.
        /// </summary>
        public static bool TryGetParticlePositions(this FrameData data,
                                                   out IReadOnlyList<float> values)
        {
            return data.TryGetFloatArray(FrameData.ParticlePositionArrayKey, out values);
        }

        /// <summary>
        /// Get the particle positions if present, else returning null.
        /// </summary>
        public static IReadOnlyList<float> GetParticlePositions(this FrameData data)
        {
            return data.GetFloatArray(FrameData.ParticlePositionArrayKey);
        }

        #endregion

        #region ParticleElements

        /// <summary>
        /// Does this frame have particle elements?
        /// </summary>
        public static bool HasParticleElements(this FrameData data)
        {
            return data.TryGetIndexArray(FrameData.ParticleElementArrayKey, out _);
        }

        /// <summary>
        /// Set the particle elements of this frame.
        /// </summary>
        public static void SetParticleElements(this FrameData data,
                                               IEnumerable<uint> values)
        {
            data.AddIndexArray(FrameData.ParticleElementArrayKey, values);
        }

        /// <summary>
        /// Try to get the particle elements array, returning true and setting the out
        /// variable values to the array if found.
        /// </summary>
        public static bool TryGetParticleElements(this FrameData data,
                                                  out IReadOnlyList<uint> values)
        {
            return data.TryGetIndexArray(FrameData.ParticleElementArrayKey, out values);
        }

        /// <summary>
        /// Get the particle elements if present, else returning null.
        /// </summary>
        public static IReadOnlyList<uint> GetParticleElements(this FrameData data)
        {
            return data.GetIndexArray(FrameData.ParticleElementArrayKey);
        }

        #endregion

        #region Particle Types

        /// <summary>
        /// Does this frame have particle types?
        /// </summary>
        public static bool HasParticleTypes(this FrameData data)
        {
            return data.TryGetStringArray(FrameData.ParticleTypeArrayKey, out _);
        }

        /// <summary>
        /// Set the particle types of this frame.
        /// </summary>
        public static void SetParticleTypes(this FrameData data,
                                            IEnumerable<string> values)
        {
            data.AddStringArray(FrameData.ParticleTypeArrayKey, values);
        }

        /// <summary>
        /// Try to get the particle types array, returning true and setting the out
        /// variable values to the array if found.
        /// </summary>
        public static bool TryGetParticleTypes(this FrameData data,
                                               out IReadOnlyList<string> values)
        {
            return data.TryGetStringArray(FrameData.ParticleTypeArrayKey, out values);
        }

        /// <summary>
        /// Get the particle types if present, else returning null.
        /// </summary>
        public static IReadOnlyList<string> GetParticleTypes(this FrameData data)
        {
            return data.GetStringArray(FrameData.ParticleTypeArrayKey);
        }

        #endregion

        #region ParticleNames

        /// <summary>
        /// Does this frame have particle names?
        /// </summary>
        public static bool HasParticleNames(this FrameData data)
        {
            return data.TryGetStringArray(FrameData.ParticleNameArrayKey, out _);
        }

        /// <summary>
        /// Set the particle names of this frame.
        /// </summary>
        public static void SetParticleNames(this FrameData data,
                                            IEnumerable<string> values)
        {
            data.AddStringArray(FrameData.ParticleNameArrayKey, values);
        }

        /// <summary>
        /// Try to get the particle names array, returning true and setting the out
        /// variable values to the array if found.
        /// </summary>
        public static bool TryGetParticleNames(this FrameData data,
                                               out IReadOnlyList<string> values)
        {
            return data.TryGetStringArray(FrameData.ParticleNameArrayKey, out values);
        }

        /// <summary>
        /// Get the particle names if present, else returning null.
        /// </summary>
        public static IReadOnlyList<string> GetParticleNames(this FrameData data)
        {
            return data.GetStringArray(FrameData.ParticleNameArrayKey);
        }

        #endregion

        #region ParticleResidues

        /// <summary>
        /// Does this frame have particle residues?
        /// </summary>
        public static bool HasParticleResidues(this FrameData data)
        {
            return data.TryGetIndexArray(FrameData.ParticleResidueArrayKey, out _);
        }

        /// <summary>
        /// Set the particle residues of this frame.
        /// </summary>
        public static void SetParticleResidues(this FrameData data,
                                               IEnumerable<uint> values)
        {
            data.AddIndexArray(FrameData.ParticleResidueArrayKey, values);
        }

        /// <summary>
        /// Try to get the particle residues array, returning true and setting the out
        /// variable values to the array if found.
        /// </summary>
        public static bool TryGetParticleResidues(this FrameData data,
                                                  out IReadOnlyList<uint> values)
        {
            return data.TryGetIndexArray(FrameData.ParticleResidueArrayKey, out values);
        }

        /// <summary>
        /// Get the particle residues if present, else returning null.
        /// </summary>
        public static IReadOnlyList<uint> GetParticleResidues(this FrameData data)
        {
            return data.GetIndexArray(FrameData.ParticleResidueArrayKey);
        }

        #endregion

        #region ResidueNames

        /// <summary>
        /// Does this frame have residue names?
        /// </summary>
        public static bool HasResidueNames(this FrameData data)
        {
            return data.TryGetStringArray(FrameData.ResidueNameArrayKey, out _);
        }

        /// <summary>
        /// Set the residue names of this frame.
        /// </summary>
        public static void SetResidueNames(this FrameData data,
                                           IEnumerable<string> values)
        {
            data.AddStringArray(FrameData.ResidueNameArrayKey, values);
        }

        /// <summary>
        /// Try to get the residue names array, returning true and setting the out
        /// variable values to the array if found.
        /// </summary>
        public static bool TryGetResidueNames(this FrameData data,
                                              out IReadOnlyList<string> values)
        {
            return data.TryGetStringArray(FrameData.ResidueNameArrayKey, out values);
        }

        /// <summary>
        /// Get the residue names if present, else returning null.
        /// </summary>
        public static IReadOnlyList<string> GetResidueNames(this FrameData data)
        {
            return data.GetStringArray(FrameData.ResidueNameArrayKey);
        }

        #endregion

        #region ResidueIds

        /// <summary>
        /// Does this frame have residue ids?
        /// </summary>
        public static bool HasResidueIds(this FrameData data)
        {
            return data.TryGetStringArray(FrameData.ResidueIdArrayKey, out _);
        }

        /// <summary>
        /// Set the residue ids of this frame.
        /// </summary>
        public static void SetResidueIds(this FrameData data,
            IEnumerable<string> values)
        {
            data.AddStringArray(FrameData.ResidueIdArrayKey, values);
        }

        /// <summary>
        /// Try to get the residue ids array, returning true and setting the out
        /// variable values to the array if found.
        /// </summary>
        public static bool TryGetResidueIds(this FrameData data,
            out IReadOnlyList<string> values)
        {
            return data.TryGetStringArray(FrameData.ResidueIdArrayKey, out values);
        }

        /// <summary>
        /// Get the residue ids if present, else returning null.
        /// </summary>
        public static IReadOnlyList<string> GetResidueIds(this FrameData data)
        {
            return data.GetStringArray(FrameData.ResidueIdArrayKey);
        }

        #endregion
        
        #region ResidueChains

        /// <summary>
        /// Does this frame have residue chains?
        /// </summary>
        public static bool HasResidueChains(this FrameData data)
        {
            return data.TryGetIndexArray(FrameData.ResidueChainArrayKey, out _);
        }

        /// <summary>
        /// Set the residue chains of this frame.
        /// </summary>
        public static void SetResidueChains(this FrameData data,
                                            IEnumerable<uint> values)
        {
            data.AddIndexArray(FrameData.ResidueChainArrayKey, values);
        }

        /// <summary>
        /// Try to get the residue chains array, returning true and setting the out
        /// variable values to the array if found.
        /// </summary>
        public static bool TryGetResidueChains(this FrameData data,
                                               out IReadOnlyList<uint> values)
        {
            return data.TryGetIndexArray(FrameData.ResidueChainArrayKey, out values);
        }

        /// <summary>
        /// Get the residue chains if present, else returning null.
        /// </summary>
        public static IReadOnlyList<uint> GetResidueChains(this FrameData data)
        {
            return data.GetIndexArray(FrameData.ResidueChainArrayKey);
        }

        #endregion

        #region Residue Count

        /// <summary>
        /// Does this frame have a residue count?
        /// </summary>
        public static bool HasResidueCount(this FrameData data)
        {
            return data.TryGetNumericValue(FrameData.ResidueCountValueKey, out _);
        }

        /// <summary>
        /// Set the residue count of this frame.
        /// </summary>
        public static void SetResidueCount(this FrameData data,
            int value)
        {
            data.AddNumericValue(FrameData.ResidueCountValueKey, value);
        }

        /// <summary>
        /// Try to get the residue count, returning true and setting the out
        /// variable values to the value if found.
        /// </summary>
        public static bool TryGetResidueCount(this FrameData data,
            out int value)
        {
            return data.TryGetIntegerValue(FrameData.ResidueCountValueKey, out value);
        }

        /// <summary>
        /// Get the residue count if present, else returning null.
        /// </summary>
        public static int? GetResidueCount(this FrameData data)
        {
            return (int?) data.GetNumericValue(FrameData.ResidueCountValueKey);
        }

        #endregion
        
        #region ChainNames

        /// <summary>
        /// Does this frame have chain names?
        /// </summary>
        public static bool HasChainNames(this FrameData data)
        {
            return data.TryGetStringArray(FrameData.ChainNameArrayKey, out _);
        }

        /// <summary>
        /// Set the chain names of this frame.
        /// </summary>
        public static void SetChainNames(this FrameData data,
                                         IEnumerable<string> values)
        {
            data.AddStringArray(FrameData.ChainNameArrayKey, values);
        }

        /// <summary>
        /// Try to get the chain names array, returning true and setting the out
        /// variable values to the array if found.
        /// </summary>
        public static bool TryGetChainNames(this FrameData data,
                                            out IReadOnlyList<string> values)
        {
            return data.TryGetStringArray(FrameData.ChainNameArrayKey, out values);
        }

        /// <summary>
        /// Get the chain names if present, else returning null.
        /// </summary>
        public static IReadOnlyList<string> GetChainNames(this FrameData data)
        {
            return data.GetStringArray(FrameData.ChainNameArrayKey);
        }

        #endregion

        #region Chain Count

        /// <summary>
        /// Does this frame have a chain count?
        /// </summary>
        public static bool HasChainCount(this FrameData data)
        {
            return data.TryGetNumericValue(FrameData.ChainCountValueKey, out _);
        }

        /// <summary>
        /// Set the chain count of this frame.
        /// </summary>
        public static void SetChainCount(this FrameData data,
            int value)
        {
            data.AddNumericValue(FrameData.ChainCountValueKey, value);
        }

        /// <summary>
        /// Try to get the chain count, returning true and setting the out
        /// variable values to the value if found.
        /// </summary>
        public static bool TryGetChainCount(this FrameData data,
            out int value)
        {
            return data.TryGetIntegerValue(FrameData.ChainCountValueKey, out value);
        }

        /// <summary>
        /// Get the chain count if present, else returning null.
        /// </summary>
        public static int? GetChainCount(this FrameData data)
        {
            return (int?) data.GetNumericValue(FrameData.ChainCountValueKey);
        }

        #endregion
        
        #region Particle Count

        /// <summary>
        /// Does this frame have a particle count?
        /// </summary>
        public static bool HasParticleCount(this FrameData data)
        {
            return data.TryGetNumericValue(FrameData.ParticleCountValueKey, out _);
        }

        /// <summary>
        /// Set the particle count of this frame.
        /// </summary>
        public static void SetParticleCount(this FrameData data,
                                            int value)
        {
            data.AddNumericValue(FrameData.ParticleCountValueKey, value);
        }

        /// <summary>
        /// Try to get the particle count, returning true and setting the out
        /// variable values to the value if found.
        /// </summary>
        public static bool TryGetParticleCount(this FrameData data,
                                               out int value)
        {
            return data.TryGetIntegerValue(FrameData.ParticleCountValueKey, out value);
        }

        /// <summary>
        /// Get the particle count if present, else returning null.
        /// </summary>
        public static int? GetParticleCount(this FrameData data)
        {
            return (int?) data.GetNumericValue(FrameData.ParticleCountValueKey);
        }

        #endregion

        #region Kinetic Energy

        /// <summary>
        /// Does this frame have a kinetic energy?
        /// </summary>
        public static bool HasKineticEnergy(this FrameData data)
        {
            return data.TryGetNumericValue(FrameData.KineticEnergyValueKey, out _);
        }

        /// <summary>
        /// Set the kinetic energy of this frame.
        /// </summary>
        public static void SetKineticEnergy(this FrameData data,
                                            double value)
        {
            data.AddNumericValue(FrameData.KineticEnergyValueKey, value);
        }

        /// <summary>
        /// Try to get the kinetic energy, returning true and setting the out
        /// variable values to the value if found.
        /// </summary>
        public static bool TryGetKineticEnergy(this FrameData data,
                                               out int value)
        {
            return data.TryGetIntegerValue(FrameData.KineticEnergyValueKey, out value);
        }

        /// <summary>
        /// Get the kinetic energy if present, else returning null.
        /// </summary>
        public static float? GetKineticEnergy(this FrameData data)
        {
            return (float?) data.GetNumericValue(FrameData.KineticEnergyValueKey);
        }

        #endregion

        #region Potential Energy

        /// <summary>
        /// Does this frame have a potential energy?
        /// </summary>
        public static bool HasPotentialEnergy(this FrameData data)
        {
            return data.TryGetNumericValue(FrameData.PotentialEnergyValueKey, out _);
        }

        /// <summary>
        /// Set the potential energy of this frame.
        /// </summary>
        public static void SetPotentialEnergy(this FrameData data,
                                              double value)
        {
            data.AddNumericValue(FrameData.PotentialEnergyValueKey, value);
        }

        /// <summary>
        /// Try to get the potential energy, returning true and setting the out
        /// variable values to the value if found.
        /// </summary>
        public static bool TryGetPotentialEnergy(this FrameData data,
                                                 out int value)
        {
            return data.TryGetIntegerValue(FrameData.PotentialEnergyValueKey, out value);
        }

        /// <summary>
        /// Get the potential energy if present, else returning null.
        /// </summary>
        public static float? GetPotentialEnergy(this FrameData data)
        {
            return (float?) data.GetNumericValue(FrameData.PotentialEnergyValueKey);
        }

        #endregion

        public static IReadOnlyList<string> GetStringArray(this FrameData data,
                                                           string id)
        {
            return data.TryGetStringArray(id, out var array) ? array : null;
        }

        public static IReadOnlyList<uint> GetIndexArray(this FrameData data,
                                                        string id)
        {
            return data.TryGetIndexArray(id, out var array) ? array : null;
        }

        public static IReadOnlyList<float> GetFloatArray(this FrameData data,
                                                         string id)
        {
            return data.TryGetFloatArray(id, out var array) ? array : null;
        }

        public static double? GetNumericValue(this FrameData data,
                                              string id)
        {
            return data.TryGetNumericValue(id, out var array) ? (double?) array : null;
        }

        public static bool TryGetIntegerValue(this FrameData data,
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