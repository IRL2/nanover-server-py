---
title: 'NanoVer Server: A Package for Serving Real-Time Multi-User/Collaborative Interactive
 Molecular Dynamics in Virtual Reality'
tags:
  - Python
  - Rust
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

NanoVer Server is a python package that facilitates real-time multi-user/collaborative interactive 
molecular dynamics (iMD) simulations. It is part of the NanoVer software ecosystem, comprising the
server component of a server-client architecture that allows multiple clients to connect and 
interact with a real-time molecular simulation. The NanoVer Server provides an interface to stream
data between a molecular simulation engine and connected clients, and interfaces with the NanoVer 
iMD-VR package to facilitate its current primary application: the exploration of molecular systems
using interactive molecular dynamics in virtual reality (iMD-VR).

# Statement of need

For decades, the family of simulation methods encompassed by the umbrella term "molecular dynamics"
(MD) have been indispensable for exploring the temporal evolution of molecular systems [@Alder1958].
MD has been applied to a wide range of chemical and biological systems, including (include some 
applications...) [Citations]. One of the major challenges when using MD simulations to explore the dynamics
of chemical and biological systems is the sampling of rare events [Citations]. Many interesting molecular
processes occur over timescales substantially longer than are computationally feasible 
to model, even when harnessing the power of high-performance computing [Citations]. This is particularly true
for biological systems, which may contain hundreds of thousands of atoms that must be explicitly
simulated in order to accurately determine the dynamics of the system [Citations]. A number of techniques have
been developed to tackle this problem, including (examples of enhanced sampling approaches...) [Citations].
(Talk about any issues to do with whatever enhanced sampling methods discussed, possibly including
needing to specify a pre-determined reaction coordinate, placing constraints on the system, etc.)

Interactive molecular dynamics (iMD) is an approach that allows researchers to bias molecular simulations
on-the-fly, without the need to pre-define the behaviour of the system of interest [Citations]. It allows researchers
to interact with molecular systems in real-time, to bias their dynamics and sample rare events of interest [Citations].
A number of programs have implementations of iMD, including (include examples of softwares that offer iMD, probably
with emphasis on those using 2-D interfaces) [Citations]. However, a limitation of many of the aforementioned iMD 
implementations is that they force researchers to interact with inherently 3-D systems via 2-D interfaces, which
means they have limited utility for exploring complex 3-D reaction coordinates [Citations?]. This problem can
be overcome by performing iMD in virtual reality (iMD-VR). Virtual reality (VR) provides a natural interface for 
visualising molecular systems by mapping the 3-D simulation space in which the molecular system is simulated to 
a 3-D virtual space that the researcher can access, and a number of programs offer such 3-D molecular visualisation. 
When combined with the tools to interact with molecular simulations, such a mapping facilitates iMD in 3-D, 
enabling researchers to explore the behaviour of molecular systems in the natural set of spatial dimensions. 
Many studies have demonstrated the utility of iMD-VR for research applications, including (cite Helen, Robin, Rhos, 
Becca, etc. as well as other iMD-VR programs (which? which ones count as iMD-VR?)).

NanoVer is a software ecosystem for performing quantitative real-time multi-user iMD simulations, with an emphasis 
on iMD-VR. In order to facilitate this, NanoVer operates using a server-client architecture, of which NanoVer Server
comprises the server component. NanoVer Server facilitates iMD simulations by providing an interface between the 
physics engine used to simulate the molecular system and the researcher(s) connecting to the simulation, equipping
the researcher(s) with tools to interact with the system in real time. NanoVer Server already interfaces with many 
multiple programs used to perform molecular simulations [Cite OpenMM and ASE], with a particular emphasis on OpenMM. 
The NanoVer Server package is written in Python, making it easy to integrate with many of the existing tools of the 
computational chemistry community and beyond. Additionally, being free and open-source (distributed under the MIT 
license [Citation?]), the NanoVer Server package can easily be customised and extended to interface with other 
physics engines.

NanoVer Server performs quantitative iMD simulations, delivering rigorous on-the-fly analytics about both the 
molecular simulation (including energies, particle forces and particle velocities) and the perturbations induced by 
the user (including the user forces, associated potential energy, and cumulative work done by the users on the
molecular system). The server allows a user to specify the interval of simulation steps at which to
publish a frame, containing information about the molecular system in its current state, to the clients connecting
to it.

The server-client architecture of NanoVer enables multiple users to connect to a single simulation running on a 
NanoVer Server simultaneously, facilitating real-time collaborative iMD and/or molecular visualisation. 
NanoVer Server interfaces with NanoVer iMD-VR, a client that enables the researcher(s) to connect to simulations 
on the server via a VR interface, which enables real-time collaborative iMD-VR simulations. The server-client 
structure enables multi-user iMD-VR to be performed using colocated (where users interact in the same physical space)
and distributed (where users can connect to the same virtual space from different physical spaces) setups. This 
flexible structure has great potential for collaborative research: previous iMD-VR frameworks have demonstrated that 
cloud computing can be used to facilitate real-time multi-user iMD-VR across large physical distances 
[Cite the appropriate Narupa paper, and any other appropriate ones?].

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
  - Introduce the concept of interactive molecular dynamics as a biasing method to explore the dynamics of rare events, and give examples of other programs that can do iMD. This might also be the right time to introduce iMD-VR as a method of performing iMD in a way that harnesses researchers' spatial intuition.
  - Introduce the NanoVer server as a method of performing real-time multi-user iMD
  - Describe the protocol for quantitative iMD employed by NanoVer


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


