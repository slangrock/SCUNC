#!/usr/bin/env python
import ROOT
import rat
import math

class FibreHandling:
    """Class to assign fibre specific variables

     To be used after the database is loaded but before loop over the root file.
     Fills and calculates all fibre specific variables needed for the analysis 
     before the analysis loop is started.

     Args:
       fibre (string) : name of the fibre used in the run
       wavelength (float) : wavlength the run was produced with

     Attributes:
       f (string) : fibre name
       wl (float) : wavelength
       sourcepos (TVector3) : position of the fibre in the detector
       sourcedir (TVector3) : direction the fibre is pointed at
       AV1_cross (TVector3) : AV crossing coordinates of the direct beam light
                              near side of the fibre
       AV2_cross (TVector3) : AV crossing coordinates of the direct beam light
                              far side of the fibre
       PSUP_cross (TVector3) : Coordinates at which the direct beam light hits the PSUP
       n_scint (float) : refractive index of the material inside the AV
       n_water (float) : refractive index of the material outside the AV
       timecuts (list of floats) : list of all cut values applied to the beam light to 
                                   separate light undergone different optical interactions
       spatialcuts (list of floats) : list of spatial cut values applied to direct beam 
                                      light and light reflected off the near side of the
                                      AV
    """



    def __init__(self,fibre):
        self.f = fibre
        self.wl = 0
        self.sourcepos = ROOT.TVector3()
        self.sourcedir = ROOT.TVector3()
        self.AV1_cross = ROOT.TVector3()
        self.AV2_cross = ROOT.TVector3()
        self.PSUP_cross = ROOT.TVector3()
        self.n_scint = 0
        self.n_water = 0
        self.timecuts = []
        self.spatialcuts = []

    def get_fibre_position(self):
        """ Function to obtain fibre positions and directions

         Args:
           f (string) : fibre name

         Returns:
           sourcepos (TVector3) : position of the fibre in the detector  
           sourcedir (TVector3) : direction the fibre is pointed at

         Fills self.sourcepos and self.sourcedir with values listed in file
        """
        #read file containing position and direction of all fibres
        input_file = open("/users/langrock/plotting_macros/SMELLIE/SMELLIE_analysis_framework/datatables/fibre_positions_new.txt",'r')
        for line in input_file:
            #find the fibre
            if self.f in line:
                words = line.split()
                #fill the position and direction vectors
                self.sourcepos.SetXYZ(float(words[1]),float(words[2]),float(words[3])) 
                self.sourcedir.SetXYZ(float(words[4]),float(words[5]),float(words[6])) 
            else:
                continue
        else:
            print "Reached end of file!"

        return self.sourcepos, self.sourcedir

    def get_crossing_points(self,wl):
        """ Function to obtain the crossing coordinates of the direct 
            beam light with the detector geometry

         Args:
          wl (float) : wavelength

         Returns:
           AV1_cross (TVector3) : AV crossing coordinates of the direct beam light
                              near side of the fibre
           AV2_cross (TVector3) : AV crossing coordinates of the direct beam light
                              far side of the fibre
           PSUP_cross (TVector3) : Coordinates at which the direct beam light hits the PSUP
           n_scint (float) : refractive index of the material inside the AV
           n_water (float) : refractive index of the material outside the AV

         Calculates the crossing points of the direct beam light with the near side of 
         the AV, the far side of the AV and the PSUP using the RAT LightPath class. Fills
         self.AV1_cross, self.AV2_cross, self.PSUP_cross and the refractive indices 
         self.n_scint and self.n_water.
        """

        #Define radius of AV and PSUP
        rav = 6000.
        rpsup = 8400.
        self.wl = wl

        #Get LightPath class and calculate refractive indices
        LightPath = rat.utility().GetLightPathCalculator()
        energy = LightPath.WavelengthToEnergy(self.wl*math.pow(10,-6))
        self.n_scint = LightPath.GetScintRI(energy)
        self.n_water = LightPath.GetWaterRI(energy)

        #calculate first crossing point with AV
        self.AV1_cross = LightPath.VectorToSphereEdge(self.sourcepos, self.sourcedir, rav, True)
	#calculate refracted source direction
        AV1_norm = -self.AV1_cross.Unit()
        newsourcedir = LightPath.PathRefraction(self.sourcedir, AV1_norm, self.n_water, self.n_scint)
	#calculate second crossing point with AV
        self.AV2_cross = LightPath.VectorToSphereEdge(self.AV1_cross, newsourcedir, rav, False)
	#calculate refracted source direction
        AV2_norm = -self.AV2_cross.Unit()
        newnewsourcedir = LightPath.PathRefraction(newsourcedir, AV2_norm, self.n_scint, self.n_water)
	#calculate crossing point with PSUP
        self.PSUP_cross = LightPath.VectorToSphereEdge(self.AV2_cross, newnewsourcedir, rpsup, False)

        return self.AV1_cross, self.AV2_cross, self.PSUP_cross, self.n_scint, self.n_water

    def cut_values(self):
        """Function to obtain cut values from file

         Args:
           f (string) : fibre name

         Fills self.spatialcuts and self.timecuts with values from file
        """
        #read file containg the cut values for all files
        input_file = open("/users/langrock/plotting_macros/SMELLIE/SMELLIE_analysis_framework/datatables/cut_values.txt",'r')
        for line in input_file:
            #find the fibre
            if self.f in line:
                words = line.split()
                for i in range(1,9):
                    #fill the list with the spatial cut values
                    self.spatialcuts.append(float(words[i]))
                for i in range(9,16):
                    #fill the list with the timing cut values
                    self.timecuts.append(float(words[i]))

            else:
                continue

        else:
            print "Reached end of file!"

