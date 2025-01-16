# Notes for "Summary" section

Aim of section: "A summary describing the high-level functionality and purpose of the software for a diverse, non-specialist audience."

Q: What software are we presenting in this paper?
A: The python implementation of the NanoVer server (nanover-server-py)

The NanoVer server (nanover-server-py) facilitates real-time multi-user interactive molecular dynamics (iMD) simulations, for which the current primary application is performing collaborative iMD simulations in virtual reality (iMD-VR).

NanoVer is a software ecosystem for performing interactive molecular dynamics simulations in virtual reality (iMD-VR), consisting of:

- the NanoVer server (nanover-server-py)
- the NanoVer VR client (nanover-imd-vr)

The NanoVer server facilitates real-time iMD, interfaces with nanover-imd-vr to facilitate our primary use case of iMD-VR. The NanoVer server comprises a set of libraries that facilitate client-server collaborative interactive molecular dynamics simulations in virtual reality... It is part of the NanoVer software ecosystem for collaborative interactive molecular dynamics simulations in virtual reality (iMD-VR).

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

# Notes for Acknowledgements section

From the JOSS website: "Acknowledgement of any financial support."

This is probably the place to acknowledge the funding that supports NanoVer.

# References

From the JOSS website: "A list of key references, including to other software addressing related needs. 
Note that the references should include full names of venues, e.g., journals and conferences, not 
abbreviations only understood in the context of a specific discipline."



