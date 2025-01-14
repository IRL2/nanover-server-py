---
title: 'NanoVer Server: A Package for Serving Real-Time Multi-User/Collaborative Interactive
 Molecular Dynamics in Virtual Reality'
tags:
  - Python
  - virtual reality
  - interactive molecular dynamics
  - iMD-VR
  - iMD
authors:
  - name: Jonathan Barnoud
    affiliation: '1'
  - name: Mark Wonnacott
    affiliation: '1'
  - name: Harry J. Stroud
    affiliation: '1'
  - name: Rhoslyn Roebuck Williams
    affiliation: '1'
  - name: Mohamed Dhouioui
    affiliation: '1'
  - name: Ludovica Aisa
    affiliation: '1'
  - name: Luis Toledo
    affiliation: '1'
  - name: David Glowacki
    affiliation: '1'
  - name: Adam McSloy
    affiliation: '2'
  - name: Phil Bates
    affiliation: '2'
  - name: Adrian Mulholland
    affiliation: '2'
affiliations:
  - index: 1
    name: Intangible Realities Laboratory, Centro Singular de Investigación en 
          Technoloxías Intelixentes (CiTIUS), Universidade de Santiago de Compostela, España
  - index: 2
    name: Centre for Computational Chemistry, School of Chemistry, 
          University of Bristol, United Kingdom
date: SUBMISSION DATE, TDB
bibliography: paper.bib

---

# Summary

NanoVer Server is a Python package that facilitates real-time multi-user/collaborative interactive 
molecular dynamics (iMD) simulations. It is part of the NanoVer software ecosystem, comprising the
server component of a server-client architecture that allows multiple clients to connect and 
interact with a real-time molecular simulation. NanoVer Server provides an interface to stream
data between a molecular simulation engine and connected clients, and includes a client that can 
connect to the server via a Python script or Jupyter notebook. Furthermore, NanoVer Server 
interfaces with the NanoVer iMD-VR package to facilitate its current primary application: the 
exploration of molecular systems using interactive molecular dynamics in virtual reality (iMD-VR).

# Statement of need

For decades, the family of simulation methods encompassed by the umbrella term "molecular dynamics"
(MD) have been indispensable for exploring the temporal evolution and properties of atomic and molecular systems 
[@alder_phase_1957; @alder_molecular_1958; @alder_studies_1959; @rahman_correlations_1964; @verlet_computer_1967; 
@mccammon_dynamics_1977]. MD has been used to study a plethora of chemical and biological systems 
[@van_gunsteren_validation_2018]; applications include the prediction of protein structures [@geng_applications_2019], 
simulation of drug docking in protein-ligand systems [@de_vivo_recent_2017], and the characterisation and 
nano-engineering of materials [@lau_nano_engineering_2018]. One of the major challenges when using MD simulations to 
explore the dynamics of chemical and biological systems is the sampling of rare events. Many interesting molecular 
processes occur over timescales substantially longer than are computationally feasible to model for most researchers, 
even when harnessing the power of high-performance computing [@yang_enhanced_2019; @kamenik_enhanced_2022]. This is 
particularly true for biological systems [@hollingsworth_molecular_2018], which often contain hundreds of thousands of 
atoms [@brooks_biomolecular_2024]. A wide variety of enhanced sampling techniques have been developed to tackle the 
issue of insufficient sampling by brute force MD simulations [@henin_enhanced_2022; @yang_enhanced_2019; 
@kamenik_enhanced_2022]. Families of such techniques include (though are by no means limited to) umbrella sampling
[@torrie_nonphysical_1977; @kastner_umbrella_2011], metadynamics [@laio_escaping_2002; @barducci_well_tempered_2008; 
@valsson_enhancing_2016], steered molecular dynamics [@izrailev_steered_1999; @park_free_2003],
replica exchange approaches [@swendsen_replica_1986; @geyer_markov_1991; @sugita_replica_exchange_1999] 
and adaptive biasing force methods [@darve_calculating_2001; @comer_adaptive_2015].
Though effective, many of these techniques have certain limitations, such as (a) requiring a priori definition of
a reaction coordinate, collective variable(s) and/or constraints or restraints on the 
system, (b) needing a large number of simulation steps and/or multiple parallel simulations, or (c) 
employing adaptive strategies that are not guaranteed to sample the desired behaviour.

Interactive molecular dynamics (iMD) is an out-of-equilibrium approach that allows researchers to bias molecular 
simulations on-the-fly, without the need to pre-define the behaviour of the system of interest [Citations]. It 
allows researchers to interact with molecular systems in real-time, to bias their dynamics and sample rare events of
interest [Citations]. A number of programs have implementations of iMD, including (include examples of softwares that
offer iMD, probably with emphasis on those using 2-D interfaces) [Citations]. However, a limitation of many of the 
aforementioned iMD implementations is that they force researchers to interact with inherently 3-D systems via 2-D 
interfaces, which means they have limited utility for exploring complex 3-D reaction coordinates [Citations?]. This 
problem can be overcome by performing iMD in virtual reality (iMD-VR). Virtual reality (VR) provides a natural 
interface for visualising molecular systems by mapping the 3-D simulation space of the molecular system 
to a 3-D virtual environment that the researcher can access, with several programs already supporting molecular 
visualisation in VR [Citations]. By combining 3-D visualisation in VR with iMD, iMD-VR enables researchers 
to interact directly with molecular simulations in the natural set of spatial dimensions. Many studies have 
demonstrated the utility of iMD-VR for research applications, including (cite Helen, Robin, Rhos, 
Becca, etc. as well as other iMD-VR programs (which? which ones count as iMD-VR?)).

NanoVer is a software ecosystem for performing quantitative real-time multi-user iMD simulations, with an emphasis 
on iMD-VR. In order to facilitate this, NanoVer operates using a server-client architecture, of which NanoVer Server
comprises the server component. NanoVer Server facilitates iMD simulations by providing an interface between the 
physics engine used to simulate the molecular system and the researcher(s) connecting to the simulation, equipping
the researcher(s) with tools to interact with the system in real time. NanoVer Server already interfaces with many 
multiple programs used to perform molecular simulations [@eastman_openmm_2024; @larsen_atomic_2017], with a particular 
emphasis on OpenMM. The NanoVer Server package is written in Python, making it easy to integrate with many of the 
existing tools of the computational chemistry community and beyond. Additionally, being free and open-source, 
the NanoVer Server package can easily be customised and extended to interface with other physics engines.

NanoVer Server performs quantitative iMD simulations, delivering rigorous on-the-fly analytics about both the 
molecular simulation (including energies, particle forces and particle velocities) and the perturbations induced by 
the user (including the collective user forces, associated potential energy, and cumulative work done). 
The server allows a user to tune the relationship between simulation time and real time during
iMD simulations by specifying the number of simulation steps to perform between the publishing of each frame
(a data structure containing information about the molecular system in its current state) to the clients connecting
to it. This feature enables accurate integration of the equations of motion between each frame without significantly 
slowing the simulation from the user perspective. To achieve quantitative iMD within this framework, the server
adopts the following blueprint between the publishing of each frame:

  1) Perform n simulation steps (as specified by the user) using any existing iMD forces and energies

  2) Calculate all iMD forces and energies applied to the system using its current configuration, 
    passing this information to the physics engine&mdash;these forces will be applied for the next n simulation steps

  3) Gather all the data about the current state of the system (including the iMD forces and energies calculated in
    step 2) and construct a frame to publish to the client(s) connecting to the simulation.

Using the blueprint for quantitative iMD described above, all the information about the iMD interactions applied
to the molecular system during an iMD simulation is delivered in the frames published by the server.

The server-client architecture of NanoVer enables multiple users to connect to a single simulation running on a 
NanoVer Server simultaneously, facilitating real-time collaborative iMD and/or molecular visualisation. 
NanoVer Server interfaces with NanoVer iMD-VR, a client that enables the researcher(s) to connect to simulations 
on the server via a VR interface, which enables real-time collaborative iMD-VR simulations. The server-client 
structure enables multi-user iMD-VR to be performed using colocated (where users interact in the same physical space)
and distributed (where users can connect to the same virtual space from different physical spaces) setups. This 
flexible structure has great potential for collaborative research: previous iMD-VR frameworks have demonstrated that 
cloud computing can be used to facilitate real-time multi-user iMD-VR across large physical distances 
[Cite the appropriate Narupa paper, and any other appropriate ones?]. Furthermore, the protocol on which NanoVer 
Server is built comprises general tools for constructing multi-user VR experiences, the scope of which extend beyond 
application in computational chemistry.

# Notes for "Statement of need" section

Aim of section: "A Statement of need section that clearly illustrates the research purpose of the
software and places it in the context of related work."

Ideas for specific topics:
  - MD is a key method used to explore the dynamics of molecular systems
  - Sampling rare events is difficult as many processes take place on timescales substantially longer than computationally feasible, particularly for large systems e.g. proteins
  - Interactive MD provides a method for a user to dynamically guide molecular systems to sample rare events
  - NanoVer enables researchers to harness their chemical intuition and spatial reasoning to perform real-time iMD simulations that can be interacted with via VR
  - Interfaces with popular molecular simulation programs (?) to facilitate real-time iMD-VR simulations
    and visualisation/analysis packages
  - Multi-user virtual environment (either colocated or distributed) facilitates collaboration between researchers
  - NanoVer can be used both in research and teaching settings as demonstrated by its precursor, Narupa (citations?) from which it is forked - Don't talk about Narupa, talk about iMD-VR (refer to iMD-VR frameworks rather than Narupa itself)
  - Server-client structure with primary server written in Python so it can be easily integrated with other tools in workflows in the comp chem community (e.g. Jupyter notebooks) and beyond
  - Can be used directly for rigorous quantitative simulations enabling direct use in research applications - rigorous on-the-fly analytics
  - Not only an iMD software, but also a molecular visualisation tool
  - The protocol on which NanoVer is built facilitates multi-person VR experiences, and it's scope extends beyond application in computational chemistry

Outline of paragraphs:

  - First paragraph: introduce MD simulations, give examples of their applications and use, discuss limitations regarding sampling rare events
  - Second paragraph: Introduce the concept of interactive molecular dynamics as a biasing method to explore the dynamics of rare events, and give examples of other programs that can do iMD. This might also be the right time to introduce iMD-VR as a method of performing iMD in a way that harnesses researchers' spatial intuition.
  - Third paragraph: Introduce the NanoVer server as a method of performing real-time multi-user iMD
  - Fourth paragraph: Describe the protocol for quantitative iMD employed by NanoVer
  - Fifth paragraph: Discuss the client-server architecture for collaborative research

First paragraph:
 
  - Talk about any issues to do with whatever enhanced sampling methods discussed, possibly including needing to 
    specify a pre-determined reaction coordinate, placing constraints on the system, requiring a large number of 
    samples/running multiple simulations to obtain desired information which still does not guarantee sampling of
    the desired transitions/transformations etc.
  - 


Ideas from the paper meeting (#1):

  - Focus on NanoVer for iMD (with emphasis on VR for interaction and visualisation)
  - Client/server structure facilitates multiplayer
  - View to crowd-sourced/gameified science approach in the future

Ideas from the paper meeting (#2): see whiteboard! echoed much of the above.


# Acknowledgements

From the JOSS website: "Acknowledgement of any financial support."

This is probably the place to acknowledge the funding that supports NanoVer.

This also seems to be the place to acknowledge contributions from other people that have helped with the project?

- ERC
- Xunta

# References

From the JOSS website: "A list of key references, including to other software addressing related needs. 
Note that the references should include full names of venues, e.g., journals and conferences, not 
abbreviations only understood in the context of a specific discipline."

# Notes for "Summary" section

Aim of section: "A summary describing the high-level functionality and purpose of the software for a diverse, non-specialist audience."

Q: What software are we presenting in this paper?
A: The python implementation of the NanoVer server (nanover-server-py)

The NanoVer server (nanover-server-py) facilitates real-time multi-user interactive molecular dynamics (iMD) simulations, for which the current primary application is performing collaborative iMD simulations in virtual reality (iMD-VR).

NanoVer is a software ecosystem for performing interactive molecular dynamics simulations in virtual reality (iMD-VR), consisting of:

- the NanoVer server (nanover-server-py)
- the NanoVer VR client (nanover-imd-vr)

The NanoVer server facilitates real-time iMD, interfaces with nanover-imd-vr to facilitate our primary use case of iMD-VR. The NanoVer server comprises a set of libraries that facilitate client-server collaborative interactive molecular dynamics simulations in virtual reality... It is part of the NanoVer software ecosystem for collaborative interactive molecular dynamics simulations in virtual reality (iMD-VR).


