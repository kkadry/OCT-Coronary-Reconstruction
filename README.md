# OCT-Coronary-Reconstruction
This repository contains the code utilized in the paper: "A Platform for High-Fidelity Patient-SpecificStructural Modeling of Atherosclerotic Arteries:From Intravascular Imaging to Three-DimensionalStress Distributions" 

The repository consists of three sections (Image Processing, 2D Geometry, and 3D Geometry) and are to be utilized serially in that order.
The Image Processing folder consists of multiple matlab codes which convert a point cloud representation of an optical coherence tomography image stack into a set of 2D boundary curves detailing the vessel borders and the plaque component borders.

The 2D Geometry folder consists of multiple python codes which convert the boundary curve outputs of Image Processing into 2D surface step files which are to be manually lofted in any CAD software into 3D solid geometries.

The 3D Geometry folder consists of multiple python codes which takes the 3D geometric step files associated with the vessel wall and plaque components and applies a series of boolean operations to ensure no overlapping volumes or surfaces occur. These overlapping volumes are then discretized to produce a computational mesh which will be imported into abaqus for structural modeling and simulation. The capability for stress adaptive remeshing also exists if mesh and stress data from a previosly simulated model is available.
