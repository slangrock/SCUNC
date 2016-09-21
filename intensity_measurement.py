#!/usr/bin/env python
import ROOT
import rat
import os
import fnmatch
import optparse

class IntensityMeasurement:
    """Class to measure the intensity of a SMELLIE run

     Attributes:
      mean (float) : mean of the nhits distribution
      rms (float) : standard deviation of the nhits distribution
      mean_nppb (TH1D) : histogram of the nhits mean versus the intensity the SMELLIE
                         root files were simulated with
      mean_nppb_up (TH1D) : histogram of the nhits mean + rms value versus the intensity
                            the SMELLIE root files were simulated with
      mean_nppb_down (TH1D) : histogram of the nhits mean - rms value versus the intensity
                              the SMELLIE root files were simulated with
    """

    def __init__(self):
        self.mean = 0
        self.rms = 0      

        self.mean_nppb = ROOT.TH1D("intensity","",10000,0,10000)
        self.mean_nppb_up = ROOT.TH1D("intensity up","",10000,0,10000)
        self.mean_nppb_down = ROOT.TH1D("intensity down","",10000,0,10000) 

    def get_nhits_hist(self, filename):
        """Function to get the mean and standard deviation of a nhits distribution

         Args:
          filename (string) : name of input SMELLIE root file

         Returns:
          mean (float) : mean of the nhits distribution
          rms (float) : standard deviation of the nhits distribution
        """

        #define nhits histogram
        h_nhits = ROOT.TH1D("number of hits","",200,0.0,1000.0)
        #loop over file
        for ds, run in rat.dsreader(filename):
            #loop through events
            for iev in range(ds.GetEVCount()):
                ev = ds.GetEV(iev)
                #extract nhits of event
                h_nhits.Fill(ev.GetNhits())

        #get mean and standard deviation from distribution
        self.mean = h_nhits.GetMean(1)
        self.rms = h_nhits.GetRMS(1)

        return self.mean, self.rms

    def fill_graphs(self, intensity):
        """Function to fill nhits versus intensiy histograms

         Args:
          intensity (float) : intensity the SMELLIE root file was simulated with
        """
        self.mean_nppb.Fill(intensity, self.mean)
        self.mean_nppb_up.Fill(intensity, self.mean + self.rms)
        self.mean_nppb_down.Fill(intensity, self.mean - self.rms)

    def fit_simulation(self, hist):       
        """Function to fit nhits versus intensity histograms using a logarithmic function.
           Fit range currently set for water-fill analysis.

         Args:
          hist (TH1D) : histogram to be fitted

         Returns:
          resulting fit parameters as floats
        """
        hist.Draw()
        f1 = ROOT.TF1("f1","[0]*log([1]+[2]*x)",0,1000)
        f1.SetParameter(0,100)
        f1.SetParameter(1,1)
        f1.SetParameter(2,0.0002)
        hist.Fit("f1","R")

        return f1.GetParameter(0), f1.GetParameter(1), f1.GetParameter(2)

    def fit_simulation_linear(self, hist):       
        """Function to fit nhits versus intensity histograms using a linear function.
           Fit range currently set for water-fill analysis.

         Args:
          hist (TH1D) : histogram to be fitted

         Returns:
          resulting fit parameters as floats
        """
        hist.Draw()
        f1 = ROOT.TF1("f1","[0]*x+[1]",0,1000)
        f1.SetParameter(0,1)
        f1.SetParameter(1,2)
        hist.Fit("f1","R")

        return f1.GetParameter(0), f1.GetParameter(1), 0

if __name__ == '__main__':
    """Script runs through all files in the directory given by the -d flag. It only selects the 
       files for the given fibre-wavelength combination (-f and -w flags). The root file names 
       have to have the structure fibre_wavelength_intensity.root
    """

    parser = optparse.OptionParser( usage = "python %prog [flags]")
    parser.add_option("-d", dest="dir_name", help="Name of directory for input files")
    parser.add_option("-f", dest="fibre", help="Fibre to be measured")
    parser.add_option("-w", dest="wavelength", help="Wavelength of fibre")
    parser.add_option("-t", dest="type", help="Data type, needs to be 'sim' or 'data'")

    (options, args) = parser.parse_args()

    if not options.dir_name:
        print 'Need directory name!'
        raise Exception
    if not options.type:
        print 'Need data type!'
        raise Exception
    if not options.fibre:
        print 'Need fibre!'
        raise Exception
    if not options.wavelength:
        print 'Need wavelength!'
        raise Exception

    #open file to save fit parameters to
    output_file = open("/data/langrock/rat-5.0-SMELLIE_analysis/" + options.fibre + "/" + options.wavelength + "_" + options.type + "_intensity_fake.txt","w")

    int_meas = IntensityMeasurement()

    #loop through given directory
    for root, dirs, filenames in os.walk(options.dir_name):
        for x in filenames:
            if fnmatch.fnmatch(x, "*.root"):
                infile = options.dir_name + x
                    
                #split the file name to get fibre, wavelength, intensity
                components = x.split("_")
                if len(components)!=0:
                    #run through all files for the chosen fibre-wavelength combination
                    if str(options.fibre) == str(components[0]) and str(options.wavelength) == str(components[1]):
                    
                        #for simulations fill nhits versus intensity histograms
                        if options.type == 'sim':
                            intensity = int(components[2].split('.')[0])
                            int_meas.get_nhits_hist(infile)
                            int_meas.fill_graphs(float(intensity))
                       
                        #for data write nhits mean and standard deviation to file
                        elif options.type == 'data':   
                            intensity = components[2].split('.')[0]
                            mean, rms = int_meas.get_nhits_hist(infile)
                            output_file.write(str(options.fibre))
                            output_file.write("\t")
                            output_file.write(str(options.wavelength))
                            output_file.write("\t")
                            output_file.write(str(intensity))
                            output_file.write("\t")
                            output_file.write(str(mean))
                            output_file.write("\t")
                            output_file.write(str(rms))
                            output_file.write("\t")
                            output_file.write(str(x))
                            output_file.write("\n")                       
                        else:
                            print "Invalid data type!"
                            raise Exception
                
                    else:
                        continue

    #for simulation after filling the nhits versus intensity histograms, fit them and write the fit parameters to file
    if options.type == 'sim':
        parameter1, parameter2, parameter3 = int_meas.fit_simulation(int_meas.mean_nppb)
        parameter1_up, parameter2_up, parameter3_up = int_meas.fit_simulation(int_meas.mean_nppb_up)
        parameter1_down, parameter2_down, parameter3_down = int_meas.fit_simulation(int_meas.mean_nppb_down)
        output_file.write(options.fibre)
        output_file.write("\t")
        output_file.write(options.wavelength)
        output_file.write("\t")
        output_file.write(str(parameter1))
        output_file.write("\t")
        output_file.write(str(parameter2))
        output_file.write("\t")
        output_file.write(str(parameter3))
        output_file.write("\t")
        output_file.write(str(parameter1_up))
        output_file.write("\t")
        output_file.write(str(parameter2_up))
        output_file.write("\t")
        output_file.write(str(parameter3_up))
        output_file.write("\t")
        output_file.write(str(parameter1_down))
        output_file.write("\t")
        output_file.write(str(parameter2_down))
        output_file.write("\t")
        output_file.write(str(parameter3_down))
        output_file.write("\n")

    output_file.close()
