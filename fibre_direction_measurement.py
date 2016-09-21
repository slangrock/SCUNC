#!/usr/bin/env python
import ROOT
import rat
import math
import fibre_handling
import optparse
import os
import random

class DirectionMeasurement:
    """Class to measure the direction of the SMELLIE fibres based on the beam center and 
       the fibre position

     Args: 
       fibre (string) : label of the fibre to be measured

     Attributes:
       f (string) : fibre label
       theta (float) : theta coordinate of the beam center, determined using Gaussian fit
       theta_sigma (float) : standard deviation of the Gaussian fit on the theta coordinate
       phi (float) : phi coordinate of the beam center, determined using Gaussian fit
       phi_sigma (float) : standard deviation of the Gaussian fit on the phi coordinate
       x (float) : x coordinate of the beam center
       y (float) : y coordinate of the beam center
       z (float) : z coordinate of the beam center
       x_error (float) : beam center x coordinate uncertainty
       y_error (float) : beam center y coordinate uncertainty
       z_error (float) : beam center z coordinate uncertainty
       positive_x (boolean) : variable to aid application of direction uncertainties
       positive_y (boolean) : variable to aid application of direction uncertainties
       positive_z (boolean) : variable to aid application of direction uncertainties
    """

    def __init__(self,fibre):
        self.f = fibre
        self.theta = 0
        self.theta_sigma = 0
        self.phi = 0
        self.phi_sigma = 0
        self.x = 0
        self.y = 0
        self.z = 0
        self.x_error = 0
        self.y_error = 0
        self.z_error = 0
        self.positive_x = True
        self.positive_y = True
        self.positive_z = True


    def get_theta_phi(self):
        """Function to determine the theta and phi coordinates of the beam center of a SMELLIE fibre

         It opens the root files returned by apply_cuts.generate_data_output for all four single laser 
         wavelengths for the specified fibre and loads the theta and phi direct beam distributions. It 
         applies a Gaussian fit to all distributions and averages the mean and the standard deviation for 
        """

        #load 495nm and get theta and phi direct beam distributions
        file_495 = ROOT.TFile("/data/langrock/rat-5.0-SMELLIE_analysis/" + self.f + "/root/495_secretinput_data_water.root")
        theta_495 = file_495.Get("theta_beam")
        phi_495 = file_495.Get("phi_beam")

        #load 446nm and get theta and phi direct beam distributions
        file_446 = ROOT.TFile("/data/langrock/rat-5.0-SMELLIE_analysis/" + self.f + "/root/446_secretinput_data_water.root")
        theta_446 = file_446.Get("theta_beam")
        phi_446 = file_446.Get("phi_beam")

        #load 407nm and get theta and phi direct beam distributions
        file_407 = ROOT.TFile("/data/langrock/rat-5.0-SMELLIE_analysis/" + self.f + "/root/407_secretinput_data_water.root")
        theta_407 = file_407.Get("theta_beam")
        phi_407 = file_407.Get("phi_beam")

        #load 375nm and get theta and phi direct beam distributions
        file_375 = ROOT.TFile("/data/langrock/rat-5.0-SMELLIE_analysis/" + self.f + "/root/375_secretinput_data_water.root")
        theta_375 = file_375.Get("theta_beam")
        phi_375 = file_375.Get("phi_beam")

        #fit 495nm theta
        c1a = ROOT.TCanvas("c1a","",1)
        theta_495.Draw()
        f1a = ROOT.TF1("f1a","gaus",0,3.5)
        f1a.SetLineColor(1)
        f1a.SetLineWidth(2)
        theta_495.Fit("f1a","R")
        
        #add fit parameters to theta and theta_sigma
        self.theta += f1a.GetParameter(1)
        self.theta_sigma += math.pow(f1a.GetParameter(2),2)
        
        #fit 495nm phi
        c1b = ROOT.TCanvas("c1b","",1)
        phi_495.Draw()
        f1b = ROOT.TF1("f1b","gaus",-3.5,3.5)
        f1b.SetLineColor(1)
        f1b.SetLineWidth(2)
        phi_495.Fit("f1b","R")
        
        #add fit parameters to phi and phi_sigma
        self.phi += f1b.GetParameter(1)
        self.phi_sigma += math.pow(f1b.GetParameter(2),2)
        
        #fit 446nm theta
        c2a = ROOT.TCanvas("c2a","",1)
        theta_446.Draw()
        f2a = ROOT.TF1("f2a","gaus",0,3.5)
        f2a.SetLineColor(1)
        f2a.SetLineWidth(2)
        theta_446.Fit("f2a","R")

        #add fit parameters to theta and theta_sigma
        self.theta += f2a.GetParameter(1)
        self.theta_sigma += math.pow(f2a.GetParameter(2),2)

        #fit 446nm phi
        c2b = ROOT.TCanvas("c2b","",1)
        phi_446.Draw()
        f2b = ROOT.TF1("f2b","gaus",-3.5,3.5)
        f2b.SetLineColor(1)
        f2b.SetLineWidth(2)
        phi_446.Fit("f2b","R")

        #add fit parameters to phi and phi_sigma
        self.phi += f2b.GetParameter(1)
        self.phi_sigma += math.pow(f2b.GetParameter(2),2)

        #fit 407nm theta
        c3a = ROOT.TCanvas("c3a","",1)
        theta_407.Draw()
        f3a = ROOT.TF1("f3a","gaus",0,3.5)
        f3a.SetLineColor(1)
        f3a.SetLineWidth(2)
        theta_407.Fit("f3a","R")

        #add fit parameters to theta and theta_sigma
        self.theta += f3a.GetParameter(1)
        self.theta_sigma += math.pow(f3a.GetParameter(2),2)
        
        #fit 407nm phi
        c3b = ROOT.TCanvas("c3b","",1)
        phi_407.Draw()
        f3b = ROOT.TF1("f3b","gaus",-3.5,3.5)
        f3b.SetLineColor(1)
        f3b.SetLineWidth(2)
        phi_407.Fit("f3b","R")

        #add fit parameters to phi and phi_sigma
        self.phi += f3b.GetParameter(1)
        self.phi_sigma += math.pow(f3b.GetParameter(2),2)

        #fit 375nm theta
        c4a = ROOT.TCanvas("c4a","",1)
        theta_375.Draw()
        f4a = ROOT.TF1("f4a","gaus",0,3.5)
        f4a.SetLineColor(1)
        f4a.SetLineWidth(2)
        theta_375.Fit("f4a","R")

        #add fit parameters to theta and theta_sigma
        self.theta += f4a.GetParameter(1)
        self.theta_sigma += math.pow(f4a.GetParameter(2),2)

        #fit 375nm phi
        c4b = ROOT.TCanvas("c4b","",1)
        phi_375.Draw()
        f4b = ROOT.TF1("f4b","gaus",-3.5,3.5)
        f4b.SetLineColor(1)
        f4b.SetLineWidth(2)
        phi_375.Fit("f4b","R")

        #add fit parameters to phi and phi_sigma
        self.phi += f4b.GetParameter(1)
        self.phi_sigma += math.pow(f4b.GetParameter(2),2)

        #calculate average theta, phi, theta_sigma and phi_sigma
        self.theta = self.theta/4
        self.phi = self.phi/4
        
        self.theta_sigma = math.sqrt(self.theta_sigma)/4
        self.phi_sigma = math.sqrt(self.phi_sigma)/4

    def convert_polar(self):
        """Function to convert the theta and phi coordinates determined by get_theta_phi 
           into Cartesian coordinates x, y, z, including the uncertainties on all coordinates
        """

        self.x =  8400*math.sin(self.theta)*math.cos(self.phi)
        self.x_error = math.sqrt(math.pow(8400*math.cos(self.theta)*math.cos(self.phi)*self.theta_sigma,2) + math.pow(8400*math.sin(self.theta)*math.sin(self.phi)*self.phi_sigma,2))
        self.y =  8400*math.sin(self.theta)*math.sin(self.phi)
        self.y_error = math.sqrt(math.pow(8400*math.cos(self.theta)*math.sin(self.phi)*self.theta_sigma,2) + math.pow(8400*math.sin(self.theta)*math.cos(self.phi)*self.phi_sigma,2))
        self.z =  8400*math.cos(self.theta)
        self.z_error = math.sqrt(math.pow(8400*math.sin(self.theta)*self.theta_sigma,2))

    def get_direction(self):
        """Function to calculate the fibre direction based on the fibre position and the 
           determined beam center from the data

         Returns:
          dir_x (float) : x coordinate of fibre direction vector
          dir_y (float) : x coordinate of fibre direction vector
          dir_z (float) : x coordinate of fibre direction vector
          dir_error_x (float) : x coordinate uncertainty of fibre direction vector
          dir_error_x (float) : x coordinate uncertainty of fibre direction vector
          dir_error_x (float) : x coordinate uncertainty of fibre direction vector
        """
        
        #get fibre position of the to be measured fibre
        val = fibre_handling.FibreHandling(self.f)
        sourcepos, sourcedir = val.get_fibre_position()
        
        #fill beam center vector with the x, y, z coordinates determined from the data
        beam_center = ROOT.TVector3()
        beam_center.SetXYZ(self.x,self.y,self.z)

        #determine fibre direction from beam center vector and fibre position
        direction = beam_center - sourcepos
        dir_x = direction.X()
        dir_y = direction.Y()
        dir_z = direction.Z()

        #calculate uncertainty on fibre direction from uncertainty on the beam center coordinates
        dir_error_x = math.sqrt(math.pow((math.sqrt(math.pow(dir_x,2) + math.pow(dir_y,2) + math.pow(dir_z,2)) - math.pow(dir_x,2)*math.pow(math.pow(dir_x,2) + math.pow(dir_y,2) + math.pow(dir_z,2),-0.5))*self.x_error/(math.pow(dir_x,2) + math.pow(dir_y,2) + math.pow(dir_z,2)),2) + math.pow(dir_x*dir_y*math.pow(math.pow(dir_x,2) + math.pow(dir_y,2) + math.pow(dir_z,2),-1.5)*self.y_error,2) + math.pow(dir_x*dir_z*math.pow(math.pow(dir_x,2) + math.pow(dir_y,2) + math.pow(dir_z,2),-1.5)*self.z_error,2))
        dir_error_y = math.sqrt(math.pow((math.sqrt(math.pow(dir_x,2) + math.pow(dir_y,2) + math.pow(dir_z,2)) - math.pow(dir_y,2)*math.pow(math.pow(dir_x,2) + math.pow(dir_y,2) + math.pow(dir_z,2),-0.5))*self.y_error/(math.pow(dir_x,2) + math.pow(dir_y,2) + math.pow(dir_z,2)),2) + math.pow(dir_y*dir_x*math.pow(math.pow(dir_x,2) + math.pow(dir_y,2) + math.pow(dir_z,2),-1.5)*self.x_error,2) + math.pow(dir_y*dir_z*math.pow(math.pow(dir_x,2) + math.pow(dir_y,2) + math.pow(dir_z,2),-1.5)*self.z_error,2))
        dir_error_z = math.sqrt(math.pow((math.sqrt(math.pow(dir_x,2) + math.pow(dir_y,2) + math.pow(dir_z,2)) - math.pow(dir_z,2)*math.pow(math.pow(dir_x,2) + math.pow(dir_y,2) + math.pow(dir_z,2),-0.5))*self.z_error/(math.pow(dir_x,2) + math.pow(dir_y,2) + math.pow(dir_z,2)),2) + math.pow(dir_z*dir_x*math.pow(math.pow(dir_x,2) + math.pow(dir_y,2) + math.pow(dir_z,2),-1.5)*self.x_error,2) + math.pow(dir_z*dir_y*math.pow(math.pow(dir_x,2) + math.pow(dir_y,2) + math.pow(dir_z,2),-1.5)*self.y_error,2))

        #normalise direction vector
        direction.SetMag(1.)
        dir_x = direction.X()
        dir_y = direction.Y()
        dir_z = direction.Z()

        #return direction coordinates and uncertainties
        return dir_x, dir_y, dir_z, dir_error_x, dir_error_y, dir_error_z

    def get_random_x_y(self):
        """Function to return two random booleans

         Returns:
          positive_x (boolean) : sign for x coordinate of the direction vector
          positive_y (boolean) : sign for y coordinate of the direction vector
        """

        #assign two variables randomly with value 0 or 1
        sign_x = random.randint(0,1)
        sign_y = random.randint(0,1)

        #define in which case booleans are true or false
        if sign_x == 1:
            self.positive_x = True
        if sign_x == 0:
            self.positive_x = False

        if sign_y == 1:
            self.positive_y = True
        if sign_y == 0:
            self.positive_y = False

        return self.positive_x, self.positive_y

    def get_random_z(self):
        """Function to return one random booleans

         Returns:
          positive_z (boolean) : sign for z coordinate of the direction vector
        """

        #assign a variable randomly with value 0 or 1
        sign_z = random.randint(0,1)

        #define in which case boolean is true or false
       
        if sign_z == 1:
            self.positive_z = True
        if sign_z == 0:
            self.positive_z = False
                
        return self.positive_z


if __name__ == '__main__':
    """Script carries out fibre direction measurements for the fibre defined by the -f flag. For each fibre
       either the direction measurement can be carried out or the order at which the systematic uncertainties 
       are applied can be defined.
    """

    parser = optparse.OptionParser( usage = "python %prog [flags]")
    parser.add_option("-f", dest="fibre", help="Fibre to be measured")
    parser.add_option("-o", dest="measurement", help="Direction measurement or Random uncertainties")

    (options, args) = parser.parse_args()

    if not options.fibre:
        print 'Need fibre!'
        raise Exception
    if not options.measurement:
        print 'Need Measurement option!'
        raise Exception

    dir_meas = DirectionMeasurement(options.fibre)

    #if direction measurement should be carried out, get all direction coordinates using the DirectionMeasurement class and write the resulting coordinates to file
    if options.measurement == "direction":
        dir_meas.get_theta_phi()
        dir_meas.convert_polar()
        x, y, z, x_err, y_err, z_err =  dir_meas.get_direction()

        output = open("/users/langrock/plotting_macros/SMELLIE/SMELLIE_analysis_framework/results/direction_measurements.txt",'a')

        output.write("fibre \t x \t y \t z \t delta x \t delta y \t delta z \n")

        output.write(options.fibre)
        output.write('\t')
        output.write(str(x))
        output.write(' \t ')
        output.write(str(y))
        output.write('\t')
        output.write(str(z))
        output.write('\t')
        output.write(str(x_err))
        output.write(' \t ')
        output.write(str(y_err))
        output.write('\t')
        output.write(str(z_err))
        output.write('\n')
       
        output.close()

    #if the application of the direction uncertainties should be determined, use the DirectionMeasurement class to assign the signs of the uncertainties and write them to file
    elif options.measurement == "systematics":
        x, y = dir_meas.get_random_x_y()
        z = dir_meas.get_random_z()

        output = open("/users/langrock/plotting_macros/SMELLIE/SMELLIE_analysis_framework/results/direction_systematics.txt",'a')

        #output.write("fibre \t x \t y \t z \n")

        output.write(options.fibre)
        output.write('\t')
        output.write(str(x))
        output.write(' \t ')
        output.write(str(y))
        output.write('\t')
        output.write(str(z))
        output.write('\n')
       
        output.close()

    else:
        print "Invalid measurement option!"
