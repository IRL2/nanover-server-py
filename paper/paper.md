---
title: 'NanoVer Server: A Python Package for Serving Real-Time Multi-User Interactive
 Molecular Dynamics in Virtual Reality'
tags:
  - Python
  - virtual reality
  - interactive molecular dynamics
  - iMD-VR
  - iMD
authors:
  - name: Harry J. Stroud
    affiliation: '1'
  - name: Mark D. Wonnacott
    affiliation: '1'
  - name: Jonathan Barnoud
    affiliation: '1'
  - name: Rhoslyn Roebuck Williams
    affiliation: '1'
  - name: Mohamed Dhouioui
    affiliation: '1'
  - name: Adam McSloy
    affiliation: '2'
  - name: Ludovica Aisa
    affiliation: '1'
  - name: Luis E. Toledo
    affiliation: '1'
  - name: Phil Bates
    affiliation: '2'
  - name: Adrian Mulholland
    affiliation: '2'
  - name: David R. Glowacki
    affiliation: '1'
affiliations:
  - index: 1
    name: Intangible Realities Laboratory, Centro Singular de Investigación en 
          Technoloxías Intelixentes (CiTIUS), Universidade de Santiago de Compostela, Spain
  - index: 2
    name: Centre for Computational Chemistry, School of Chemistry, 
          University of Bristol, United Kingdom
date: SUBMISSION DATE, TDB
bibliography: paper.bib

---

# Summary

(New version)
NanoVer Server is a Python package facilitating real-time multi-user interactive molecular dynamics (iMD) simulations.

<!---
- Facilitating?
- Is the sentence too wordy?
--->

It is the primary software of NanoVer, interfacing with standard molecular dynamics packages to provide a familiar Python and Jupyter based workflow for running iMD simulations and serving them to local and remote clients over the network. This forms the base of NanoVer's primary application: the exploration of molecular systems using interactive molecular dynamics in virtual reality (iMD-VR), which is achieved by connecting a dedicated VR client from the NanoVer iMD-VR package.

<!---
- Express succinctly difference between Python and VR clients
--->

(Original version)
NanoVer Server is a Python package that facilitates multi-user interactive molecular dynamics (iMD) simulations in real-time. It is part of the NanoVer software ecosystem, comprising the server component of a server-client architecture that allows multiple clients to connect and interact with a real-time molecular simulation. NanoVer Server provides an interface to stream data between a molecular simulation engine and connected clients, and includes a client that can connect to the server via a Python script or Jupyter notebook. Furthermore, NanoVer Server interfaces with the NanoVer iMD-VR package to facilitate its current primary application: the exploration of molecular systems using interactive molecular dynamics in virtual reality (iMD-VR).

(Suggested version)
NanoVer Server is a Python package facilitating real-time multi-user interactive molecular dynamics (iMD) simulations. It is the primary software of NanoVer, interfacing with standard molecular dynamics packages to provide a familiar Python and Jupyter based workflow for running iMD simulations and serving them to local and remote clients over the network. This forms the base of NanoVer's primary application: the exploration of molecular systems using interactive molecular dynamics in virtual reality (iMD-VR), which can be achieved by connecting a dedicated VR client from the NanoVer iMD-VR package.

# Statement of need

For decades, the family of simulation methods encompassed by the umbrella term "molecular dynamics" (MD) have been indispensable for exploring the temporal evolution and properties of atomic and molecular systems[@alder_phase_1957; @alder_molecular_1958; @alder_studies_1959; @rahman_correlations_1964; @verlet_computer_1967; @mccammon_dynamics_1977]. MD has been used to study a plethora of chemical and biological systems[@van_gunsteren_validation_2018]; applications include the prediction of protein structures [@geng_applications_2019], simulation of drug docking in protein-ligand systems [@de_vivo_recent_2017], and the characterisation and nano-engineering of materials [@lau_nano_engineering_2018].

<!---
- Should has been be changed to make it more contemporary
One of the major challenges when using MD simulations to explore the dynamics of chemical and biological systems is the sampling of rare events.
- Reword sentence to remove repetition of chemical and biological 
- change sampling?
--->

Many significant molecular processes occur over timescales substantially longer than are computationally feasible to model for most researchers, even when harnessing the power of high-performance computing [@yang_enhanced_2019; @kamenik_enhanced_2022]. This is particularly true for biological systems [@hollingsworth_molecular_2018], which often comprise hundreds of thousands of atoms [@brooks_biomolecular_2024]. A wide variety of enhanced sampling techniques have been developed to tackle the issue of insufficient sampling by brute force MD simulations [@henin_enhanced_2022; @yang_enhanced_2019; @kamenik_enhanced_2022]. Families of such techniques include (though are by no means limited to) umbrella sampling[@torrie_nonphysical_1977; @kastner_umbrella_2011], metadynamics [@laio_escaping_2002; @barducci_well_tempered_2008; @valsson_enhancing_2016], steered molecular dynamics [@izrailev_steered_1999; @park_free_2003], replica exchange approaches [@swendsen_replica_1986; @geyer_markov_1991; @sugita_replica_exchange_1999]and adaptive biasing force methods [@darve_calculating_2001; @comer_adaptive_2015]. Though effective, many of these techniques have certain limitations, such as (a) requiring a priori definition of a reaction coordinate, collective variable(s) and/or constraints or restraints on the system, (b) needing a large number of simulation steps and/or multiple parallel simulations, or (c)employing adaptive strategies that are not guaranteed to sample the desired behaviour.

<!---
- Think about "have certain limitations", e.g. specific, suffer from, etc.
--->

Interactive molecular dynamics (iMD) is a nonequilibrium approach that allows researchers to bias atomic/molecular simulations on-the-fly, without the need to pre-define the behaviour of the system[@stone_system_2001; @rapaport_interactive_1997; @oconnor_sampling_2018; @oconnor_interactive_2019]. 

<!---
- Nonequilibrium (do we need to explain further)
--->

iMD enables researchers to interact with molecular systems in real-time, to bias their dynamics and sample rare events of interest [@stone_system_2001].

<!---
- Can we merge the first two sentences?
--->

A number of programs have implementations of iMD, including for example TeraChem [@luehr_ab_2015], SCINE [@weymuth_immersive_2021], and NAMD/VMD [@stone_system_2001]. These iMD implementations enable researchers to interact with inherently 3-D systems via haptic interfaces visualized on 2-D screens. Fewer environments are available that enable researchers to interact with molecular dynamics simulations in their native 3-D environment, which is especially important for exploring complex molecular structural transformations [@oconnor_sampling_2018]. This problem can be overcome by performing iMD in virtual reality (iMD-VR), which provides a natural interface for visualising molecular systems by mapping the 3-D simulation space of the molecular system to a 3-D virtual environment that the researcher can inhabit. While several programs already enable molecular visualisation in VR [@doutreligne_unitymol_2014; @bennie_virtual_2023; @pettersen_ucsf_2021; @ozvoldik_yasara_2023; @cassidy_proteinvr_2020; @cortes_rodriguez_molecularwebxr_2025], fewer support the combination of VR-enabled 3-D visualisation in VR with real-time iMD. The iMD-VR combination enables researchers to interact directly with molecular simulations in an environment which is natively 3-D, where researchers can reach out and ‘touch’ molecular simulations as if they were tangible objects. 

<!---
- Native?
- Do we need "as if they were tangible objects"? 
--->

Many studies have demonstrated the utility of iMD-VR for research applications,  in areas spanning protein-drug binding, [@deeks_interactive_2020; @deeks_interactive_docking_2020; @henry_chan_discovery_2021; @walters_emerging_2022], protein conformational dynamics [@juarez_jimenez_combining_2020], machine-learning potential energy surfaces [@amabilino_training_2019; @amabilino_training_2020], discovering reaction networks[@shannon_exploring_2021], and chemistry education [@bennie_teaching_2019].

<!---
- Subheading?
--->

NanoVer is a software ecosystem for performing quantitative real-time multi-user iMD simulations, with an emphasis on iMD-VR. In order to facilitate this, NanoVer operates using a server-client architecture. NanoVer Server facilitates iMD by providing an interface between the simulation engine used for the molecular system and the researcher(s) connecting to the simulation, equipping the researcher(s) with tools to interact with the system in real time. NanoVer Server already interfaces with multiple programs used to perform molecular simulations [@eastman_openmm_2024; @larsen_atomic_2017], with particular support for OpenMM. NanoVer Server is written in Python, making it easy to integrate with many of the existing tools of the computational chemistry community and beyond. Additionally, being free and open-source, NanoVer Server can be customised and extended to interface with other molecular simulation engines.

NanoVer Server performs quantitative iMD simulations, delivering rigorous on-the-fly metrics about both the molecular simulation e.g. energies, particle forces and particle velocities, and the perturbations induced by the user e.g. the collective user forces, associated potential energy and cumulative work done.

<!---
- return to this
--->

A user can tune the relationship between simulation time and real time during iMD simulations by specifying the number of simulation steps to perform between the publishing of each frame (a data structure containing information about the molecular system in its current state) to the clients connecting to it. This feature enables accurate integration of the equations of motion between each frame without significantly slowing the simulation from the user perspective.

<!---
- retain accuracy and usability (responsiveness)
- the point is not to watch simulation in slow motion
- significantly slowing simulation, do we have any evidence
--->

To achieve quantitative iMD, the server adopts the following blueprint between the publishing of each frame:

  1) Perform _n_ simulation steps (as specified by the user) using existing iMD forces and energies

  2) Calculate all iMD forces and energies applied to the system using its current configuration, passing this information to the simulation engine&mdash;these forces will be applied for the next _n_ simulation steps

  3) Compile the data about the current state of the system (including the iMD forces and energies calculated in step 2) into a frame and publish to the connected client(s)

<!---
- Point 1: existing? How to explain, current
- Replace/support with figure?
--->

Using the blueprint for quantitative iMD described above, all the information about the iMD interactions applied to the molecular system during an iMD simulation is delivered in the frames published by the server.

The server-client architecture of NanoVer enables multiple users to connect simultaneously to a server. NanoVer Server interfaces with NanoVer iMD-VR, a client that enables the researcher(s) to connect to simulations on the server via a VR interface, which enables real-time collaborative iMD-VR simulations. The server-client structure enables multi-user iMD-VR to be performed using colocated (where users interact in the same physical space) and distributed (where users can connect to the same virtual space from different physical spaces) setups. This flexible structure has great potential for collaborative research: previous iMD-VR frameworks have demonstrated that cloud computing can be used to facilitate real-time multi-user iMD-VR across large physical distances[@jamieson_binnie_narupa_2020; @deeks_interactive_2020; @deeks_free_2023; @glowacki_group_2022]. Furthermore, the protocol on which NanoVer Server is built comprises general tools for constructing multi-user VR experiences, the scope of which extend beyond application in computational chemistry [@glowacki_group_2022].

<!---
- multiple users/collaboration
- Interface with iMD-VR program
- colocated and distributed setups possible
- Possibility for cloud, collaboration more extensively
- general tools for building VR experiences beyond comp chem
- cite isness-D (second to last/last paragraph)
--->

# Availability

NanoVer Server can be installed as a [conda package](https://anaconda.org/irl/nanover-server), the documentation is
available on the project's [readthedocs](https://irl2.github.io/nanover-docs/index.html) and the source code is 
available on [GitHub](https://github.com/IRL2/nanover-server-py). More information about 
[how to install NanoVer Server](https://irl2.github.io/nanover-docs/installation.html) can be found in the 
documentation.

NanoVer Server also has a suite of 
[Jupyter notebook tutorials](https://github.com/IRL2/nanover-server-py/tree/main/examples), many of which demonstrate 
how to use the server in conjunction with the 
[NanoVer iMD-VR client](https://github.com/IRL2/nanover-imd-vr). For new users, the "basics" notebooks are a good 
place to start, which provide hands-on tutorials that introduce many of the key features of NanoVer.

# Acknowledgements

This work is supported by the European Research Council under the European Union’s Horizon 2020 research and 
innovation programme through consolidator Grant NANOVR 866559, the Axencia Galega de Innovación for funding as an 
Investigador Distinguido through the Oportunius Program, and the Xunta de Galicia (Centro de investigación de 
Galicia accreditation 2019–2022, ED431G-2019/04).

INCLUDE BRISTOL FUNDING 

LIST OF CONTRIBUTIONS MADE BY AUTHORS

We would also like to thank all previous contributors to NanoVer: a full list of contributors can be found 
[here](https://github.com/IRL2/nanover-server-py/blob/main/CONTRIBUTORS.md).

We would like to thank the developers of the Python packages used by NanoVer Server, and the Python community as a 
whole for providing a great ecosystem that facilitates development of free and open source scientific software.

# References