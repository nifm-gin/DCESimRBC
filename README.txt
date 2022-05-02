DCESIM - MR Signal simulation in the scope of DCE-MRI like experiment
----

DCESIM is a tool that simulate a DCE-MRI like experiment.
Detailed description of the algorithm and applications can be found in [1]
It runs on Matlab.

[S, C, t] = VoxelSim('Param1.txt')

produces the MR Signal S, the concentration profil C at the corresponding acquisition time points t.

The algorithm takes into account various phenomena that impact the MR signal:
- The mesoscopic magnetic field perturbations induced by the vessel and/or cells.
- The diffusion of the water
- The flow of the CA within the vascular pool
- The exchange of the CA between the vascular pool and the extravascular extracellular compartment.
- The diffusion of the CA within the extravascular extracellular compartment.

Parameters that controled for the above phenomena can be modified. 
They are stored in a text file (Param.txt) which is the input of VoxelSim.

Two Param.txt files are provided as examples.

Note also that:
The Arterial Input Function can be modified by implementing your own modification in AIFGen.m.
Different MR sequences than the one provided can be designed by replacing the parameter Seq.Name by your own @MRSequence in the parameter text file.

History:

2/7/2013   - 1.1   - RF pulse implementation updated to deal properly with RF spoiling
                   - Display updated

9/1/2013   - 1.0.1 - Minor bugs fixed
				   - Add support for a GE RF Gradient spoiled sequence (Seq_GEGradRFSpoiled)
30/11/2012 - 1.0.0 - Initial Release

References:
[1] Pannetier N., Debacker C., Mauconduit F., Christen T., Barbier EL, Plos One, 2013

Contact:
nicolas D0T pannetier AT gmail D0T com