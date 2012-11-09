CoCoMANS Materials Mini-Apps
============================

Description
-----------

CoCoMANS aims at developing new applications-based, self-consistent, two-way, scale-bridging methods that have broad applicability to the targeted science, will map well to emerging heterogeneous computing models (while concurrently guiding hardware and software to maturity), and provide the algorithmic acceleration necessary to probe new scientific challenges at unprecedented scales. In the domain of Materials Science, we want to implement a multi scale method bridging both time and length scales and apply our co-design methodology on the initial mini-app. The challenge for the physics-based problem will be the transfer of information between the two chosen scales, and the set-up of simulations for the higher-order problem (i.e., the atomistic molecular dynamics simulations). Potential applications for such simulations encompass examples where the length scale, time scale, or both scales of a molecular-dynamics (MD) simulation is too large for current computers or when continuum methods require more than behavior laws to generate accurate results. In this project, we focus on the problem of shocks in materials.


Installation
------------

There are currently two mini-apps in this directory: the older `MD`, which addresses only the higher-order problem, and the newer `matminiapp`, which attempts a coupled higher-order/lower-order solution.  Each of these is implemented as a single C++ source file and associated header file.  There's no `Makefile`; simply compile each mini-app on the command line:

    g++ -O2 -g -o MD MD.cpp
    g++ -O2 -g -o matminiapp matminiapp.cpp


Usage
-----

The `MD` mini-app takes either zero or two command line arguments.  If no arguments are provided, it runs a single MD simulation.  Otherwise, the first argument is the total number of MD simulations to run and the second argument is the number of MD simulations to run concurrently.  This serves as a performance mock-up for a coupled higher-order/lower-order simulation.

The `matminiapp` mini-app takes no arguments.

Input file names are hard-wired into the code.  See the documentation for descriptions of what files are used and how to interpret and modify their contents.


License
-------

Los Alamos National Security, LLC (LANS) owns the copyright to the CoCoMANS Materials Mini-Apps, which it identifies internally as LA-CC-12-045.  The license is BSD 3.  See [LICENSE.md](https://github.com/cocomans/materials/blob/master/LICENSE.md) for the full text.
