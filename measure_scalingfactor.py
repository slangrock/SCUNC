#!/usr/bin/env python
import ROOT
import optparse
import fnmatch
import os
import plot_style
import math

class ScalingFactorMeasurement:
    """Class to measure the scattering length scaling factor of SMELLIE data. Needs the text files
       returned by the functions on apply_cuts.py.

     Attributes:
      gr_scatt (TGraphErrors) : graph of the ratio of scattered hits over total hits over a simulated
                                run versus the scattering length scaling factor the run was simulated 
                                with 
      gr_beam (TGraphErrors) : graph of the ratio of direct beam hits over total hits over a simulated
                               run versus the scattering length scaling factor the run was simulated 
                               with                                 
    """

    def __init__(self):

        self.gr_scatt = ROOT.TGraphErrors(5)
        self.gr_beam = ROOT.TGraphErrors(5)

    def get_sim(self, dir_name, wavelength):
        """Function to fill the scattering/direct beam event versus scaling factor graphs

         Args:
          dir_name (string) : directory that points to the text files created by the functions
                              in apply_cuts.py
          wavelength (string) : wavelength of the runs to be measured
        """

        graph_point = 0

        #run through given directory
        for root, dirs, filenames in os.walk(dir_name):
            #sort the filenames to ensure graph points are set in the right order
            filenames.sort()
            for x in filenames:
                if fnmatch.fnmatch(x, wavelength + "*_water.txt"):
                    infile = dir_name + x

                    print x
                    
                    #split filename
                    components = x.split("_")
                    if len(components)!=0:
                        if components[1] != 'secretinput':
                            meas_file = open(infile,'r')
                            #read first line of text files and fill gr_scatt and gr_beam
                            firstline = meas_file.readline()
                            values = firstline.split()
                            self.gr_scatt.SetPoint(graph_point,float(values[1]),float(values[4]))
                            self.gr_scatt.SetPointError(graph_point,0.0,float(values[6]))
                            self.gr_beam.SetPoint(graph_point,float(values[1]),float(values[9]))
                            self.gr_beam.SetPointError(graph_point,0.0,float(values[11]))
                            graph_point += 1

    def get_data(self, dir_name, wavelength):
        """Function to get the scattered events over total events and direct beam events over total
           events and accompanying uncertainty for the data files.

         Args:
          dir_name (string) : directory that points to the text files created by the functions
                              in apply_cuts.py
          wavelength (string) : wavelength of the runs to be measured

         Returns:
          hit ratios for the scattering region and direct beam region and their uncertainties as floats
        """
        for root, dirs, filenames in os.walk(dir_name):
            for x in filenames:
                if fnmatch.fnmatch(x, wavelength + "_secretinput_water.txt"):
                    infile = dir_name + x

                    components = x.split("_")
                    if len(components)!=0:
                        meas_file = open(infile,'r')
                        firstline = meas_file.readline()
                        values = firstline.split()

                        return float(values[4]), float(values[6]), float(values[9]), float(values[11])

    def fit_scaling_factor(self, fibre, wl):
        """Function to fit the scattering and direct beam region versus scaling factor graphs using
           a linear fit.

         Args:
          fibre (string) : fibre label of the measured fibre
          wl (string) : wavelength the fibre was run at

         Returns:
          fit parameters and their uncertainties for both fits as floats
        """
        style = plot_style.PlotStyle()
        style.set_style()
        ROOT.gROOT.SetStyle("clearRetro")

        #draw scattering region graph
        c1 = ROOT.TCanvas("c1","c1",1)
        c1.Draw()
        self.gr_scatt.GetXaxis().SetTitle("s")
        self.gr_scatt.GetYaxis().SetTitle("Ratio of scattered vs total hits")
        self.gr_scatt.GetXaxis().SetLabelSize(0.04)
        self.gr_scatt.GetYaxis().SetLabelSize(0.04)
        self.gr_scatt.GetXaxis().SetTitleSize(0.05)
        self.gr_scatt.GetYaxis().SetTitleSize(0.045)
        self.gr_scatt.GetYaxis().SetTitleOffset(1.0)
        self.gr_scatt.GetXaxis().SetTitleOffset(0.8)
        self.gr_scatt.SetMarkerColor(30)
        self.gr_scatt.SetLineColor(30)
        self.gr_scatt.SetLineWidth(2)
        self.gr_scatt.SetMarkerStyle(8)
        self.gr_scatt.SetMarkerSize(0.5)
        self.gr_scatt.Draw("acp")

        #apply linear fit
        f1 = ROOT.TF1("f1","pol1",0.2, 2.0)
        f1.SetLineColor(1)
        f1.SetLineWidth(2)
        self.gr_scatt.Fit("f1","R")

        l1 = ROOT.TLegend(0.2,0.59,0.55,0.8,"","brNDC")
        l1.AddEntry(self.gr_scatt,"Simulation", "ple")
        l1.AddEntry(f1,"Fit on Simulation", "l")
        l1.SetFillColor(0)
        l1.SetTextSize(0.05)
        l1.Draw("same")

        c1.Update()

        pad = ROOT.TPad("pad","pad",0.35,0.2,0.85,0.6)

        label = ROOT.TPaveText(0.1,0.0,0.89,0.3,"br")

        pad.SetFillStyle(4000)
        pad.SetFillColor(0)
        pad.Draw()
        pad.cd()

        label.SetFillColor(0)
        label.SetLineColor(0)
        label.SetShadowColor(0)
        label.SetTextSize(0.1)
        #label.AddText("SNO+ Preliminary")
        label.AddText("MC simulations, water-filled detector")
        label.Draw()

        c1.Update()

        #save canvas
        c1.Print("/data/langrock/rat-5.0-SMELLIE_analysis/" + fibre + "/plots/" + wl + "_scatt_fit_noisedown.pdf","pdf")

        #draw direct beam region graph
        c2 = ROOT.TCanvas("c2","c2",1)
        c2.Draw()
        self.gr_beam.GetXaxis().SetTitle("s")
        self.gr_beam.GetYaxis().SetTitle("Ratio of in-beam vs total hits")
        self.gr_beam.GetXaxis().SetLabelSize(0.04)
        self.gr_beam.GetYaxis().SetLabelSize(0.04)
        self.gr_beam.GetXaxis().SetTitleSize(0.05)
        self.gr_beam.GetYaxis().SetTitleSize(0.045)
        self.gr_beam.GetYaxis().SetTitleOffset(1.0)
        self.gr_beam.GetXaxis().SetTitleOffset(0.8)
        self.gr_beam.SetMarkerColor(50)
        self.gr_beam.SetLineColor(50)
        self.gr_beam.SetLineWidth(2)
        self.gr_beam.SetMarkerStyle(8)
        self.gr_beam.SetMarkerSize(0.5)
        self.gr_beam.Draw("acp")

        #apply linear fit
        f2 = ROOT.TF1("f2","x++1",0.2, 2.0)
        f2.SetLineColor(1)
        f2.SetLineWidth(2)
        self.gr_beam.Fit("f2","R")

        l2 = ROOT.TLegend(0.5,0.59,0.85,0.8,"","brNDC")
        l2.AddEntry(self.gr_beam,"Simulation", "ple")
        l2.AddEntry(f2,"Fit on Simulation", "l")
        l2.SetFillColor(0)
        l2.SetTextSize(0.05)
        l2.Draw("same")

        c2.Update()

        pad2 = ROOT.TPad("pad2","pad2",0.15,0.2,0.65,0.6)
        label_1 = ROOT.TPaveText(0.1,0.0,0.89,0.3,"br")

        pad2.SetFillStyle(4000)
        pad2.SetFillColor(0)
        pad2.Draw()
        pad2.cd()

        label_1.SetFillColor(0)
        label_1.SetLineColor(0)
        label_1.SetShadowColor(0)
        label_1.SetTextSize(0.1)
        #label_1.AddText("SNO+ Preliminary")
        label_1.AddText("MC simulations, water-filled detector")
        label_1.Draw()

        c2.Update()
        
        #save canvas
        c2.Print("/data/langrock/rat-5.0-SMELLIE_analysis/" + fibre + "/plots/" + wl + "_beam_fit_noisedown.pdf","pdf")

        #return fit parameters
        return f1.GetParameter(1), f1.GetParError(1), f1.GetParameter(0), f1.GetParError(0), f2.GetParameter(0), f2.GetParError(0), f2.GetParameter(1), f2.GetParError(1)

if __name__ == '__main__':
    """Script to measure the scattering length scaling factor of a given data file for each user defined 
       fibre-wavelength combination. To run, the text files produced by apply_cuts.py have to be saved 
       in the directory accessed via the -d flag.
    """

    parser = optparse.OptionParser( usage = "python %prog [flags]")
    parser.add_option("-d", dest="dir_name", help="Name of directory for input files")
    parser.add_option("-f", dest="fibre", help="Fibre to be measured")
    parser.add_option("-w", dest="wavelength", help="Wavelength of fibre")

    (options, args) = parser.parse_args()

    if not options.dir_name:
        print 'Need directory name!'
        raise Exception
    if not options.fibre:
        print 'Need fibre!'
        raise Exception
    if not options.wavelength:
        print 'Need wavelength!'
        raise Exception

    #each fibre has its own folder in the directory given y -d
    dirname = options.dir_name + "/" + options.fibre + "/"

    scale_factor = ScalingFactorMeasurement()
    #fill scattering region and direct beam region graph from simulated files
    scale_factor.get_sim(dirname, options.wavelength)
    #get ratios for data run
    scatt_c, scatt_c_error, beam_c, beam_c_error = scale_factor.get_data(dirname, options.wavelength)

    #fit the graphs from the simulation
    scatt_a, scatt_a_error, scatt_b, scatt_b_error, beam_a, beam_a_error, beam_b, beam_b_error = scale_factor.fit_scaling_factor(options.fibre, options.wavelength)
 
    #calculate the scaling factor from the data ratios and the fit parameters
    scaling_factor_scatt = (scatt_c - scatt_b)/scatt_a
    scaling_factor_scatt_error = math.sqrt(math.pow(scatt_b*scatt_c_error/scatt_a,2) + math.pow(scatt_c*scatt_b_error/scatt_a,2) + math.pow((scatt_c - scatt_b)*scatt_a_error/math.pow(scatt_a,2),2))
    scaling_factor_beam = (beam_c - beam_b)/beam_a
    scaling_factor_beam_error = math.sqrt(math.pow(beam_b*beam_c_error/beam_a,2) + math.pow(beam_c*beam_b_error/beam_a,2) + math.pow((beam_c - beam_b)*beam_a_error/math.pow(beam_a,2),2))

    #save all determined fit parameters to file for each fibre-wavelength combination
    output_fit = open("/users/langrock/plotting_macros/SMELLIE/SMELLIE_analysis_framework/results/fit_parameters_noisedown.txt",'a')

    output_fit.write(options.fibre)
    output_fit.write('\t')
    output_fit.write(options.wavelength)
    output_fit.write('\t')
    output_fit.write('scatt: ')
    output_fit.write(str(scatt_a))
    output_fit.write(' +/- ')
    output_fit.write(str(scatt_a_error))
    output_fit.write('\t')
    output_fit.write(str(scatt_b))
    output_fit.write(' +/- ')
    output_fit.write(str(scatt_b_error))
    output_fit.write('\t')
    output_fit.write('beam: ')
    output_fit.write(str(beam_a))
    output_fit.write(' +/- ')
    output_fit.write(str(beam_a_error))
    output_fit.write('\t')
    output_fit.write('beam: ')
    output_fit.write(str(beam_b))
    output_fit.write(' +/- ')
    output_fit.write(str(beam_b_error))
    output_fit.write('\n')
        
    output_fit.close()

    #save all calculated scaling factors to file, including their uncertainties, for all fibre-wavelength combinations
    output = open("/users/langrock/plotting_macros/SMELLIE/SMELLIE_analysis_framework/results/scaling_factors_noisedown.txt",'a')

    output.write(options.fibre)
    output.write('\t')
    output.write(options.wavelength)
    output.write('\t')
    output.write('scatt: ')
    output.write(str(scaling_factor_scatt))
    output.write(' +/- ')
    output.write(str(scaling_factor_scatt_error))
    output.write('\t')
    output.write('beam: ')
    output.write(str(scaling_factor_beam))
    output.write(' +/- ')
    output.write(str(scaling_factor_beam_error))
    output.write('\n')
        
    output.close()
