#!/usr/bin/env python
import ROOT
import rat
import optparse
import os
import fnmatch
import plot_style

def get_wavelength_bias(h_bias,infile,wl):
    """Function to return a bias histogram for all fibres at a given wavelength. It calculates the difference 
       between the assumed scattering length scaling factor (here 1.0) and the measured scaling factor of a
       range of SMELLIE simulations.

     Args:
      h_bias (TH1D) : histogram to be filled with the calculated bias
      infile (string) : input file containing the measured scaling factors of the SMELLIE simulations
      wl (string) : wavelength the bias measurement is carried out for

     Returns:
      h_bias (TH1D) : filled histogram
    """

    for line in infile:
        words = line.split()
        if len(words)!=0:

            if words[1] == wl:
                scatt = float(words[3])
                beam = float(words[7])

                h_bias.Fill(1-scatt)
                h_bias.Fill(1-beam)

    return h_bias

def get_fibre_bias(h_bias,infile,fibre):
    """Function to return a bias histogram for a specific fibre at all wavelengths. It calculates the difference 
       between the assumed scattering length scaling factor (here 1.0) and the measured scaling factor of a
       range of SMELLIE simulations.

     Args:
      h_bias (TH1D) : histogram to be filled with the calculated bias
      infile (string) : input file containing the measured scaling factors of the SMELLIE simulations
      fibre (string) : fibre the bias measurement is carried out for

     Returns:
      h_bias (TH1D) : filled histogram
    """

    for line in infile:
        words = line.split()
        if len(words)!=0:

            if words[0] == fibre:
                
                scatt = float(words[3])
                beam = float(words[7])

                h_bias.Fill(1-scatt)
                h_bias.Fill(1-beam)

    return h_bias

def get_all_bias(h_bias,infile):
    """Function to return a bias histogram for all fibre-wavelength combinations. It calculates the difference 
       between the assumed scattering length scaling factor (here 1.0) and the measured scaling factor of a
       range of SMELLIE simulations.

     Args:
      h_bias (TH1D) : histogram to be filled with the calculated bias
      infile (string) : input file containing the measured scaling factors of the SMELLIE simulations

     Returns:
      h_bias (TH1D) : filled histogram
    """

    for line in infile:
        words = line.split()
        if len(words)!=0:
            scatt = float(words[3])
            beam = float(words[7])

            h_bias.Fill(1-scatt)
            h_bias.Fill(1-beam)

    return h_bias


if __name__ == '__main__':
    """Script to determine the bias of the SMELLIE scattering length measurement. it can be carried out for
       all fibre-wavelength combinations as well as for a specific fibre or wavelength only. Fills a
       histogram and applies a Gaussian fit. The resulting plot is saved. The bias of the scaling factor 
       measurement is equivalent to the mean of hte Gaussian fit. It is designed to loop through a directory
       containing a set of files containing the measured scaling factor values of SMELLIE simulations.
    """

    parser = optparse.OptionParser( usage = "python %prog [flags]")
    parser.add_option("-t", dest="type", help="Bias type, needs to be 'full', 'fibre' or 'wavelength'")
    parser.add_option("-f", dest="fibre", help="Fibre to be measured, only necessary if type == 'fibre'")
    parser.add_option("-w", dest="wavelength", help="Wavelength of fibre, only necessary of type == 'wavelength'")

    (options, args) = parser.parse_args()

    if not options.type:
        print 'Need bias type!'
        raise Exception

    style = plot_style.PlotStyle()
    style.set_style()
    ROOT.gROOT.SetStyle("clearRetro")

    #define histogram and initialise string the canvas is saved at
    h_bias = ROOT.TH1D("h_bias","", 50,-0.1,0.1)
    canvas_string = ""

    #loop through given directory
    for root, dirs, filenames in os.walk("/users/langrock/plotting_macros/SMELLIE/SMELLIE_analysis_framework/results/"):
        for x in filenames:
            if fnmatch.fnmatch(x, "scaling_factors_bias_*.txt"):
                #open file
                infile = open("/users/langrock/plotting_macros/SMELLIE/SMELLIE_analysis_framework/results/" + x,"r")             
                #full bias measurement
                if options.type == "full":
                    h_bias.Add(get_all_bias(h_bias,infile))
                    canvas_string = "/data/langrock/rat-5.0-SMELLIE_analysis/scaling_factor_bias.pdf"

                #fibre specific bias measurement
                elif options.type == "fibre":
                    h_bias.Add(get_fibre_bias(h_bias,infile,options.fibre))
                    canvas_string = "/data/langrock/rat-5.0-SMELLIE_analysis/scaling_factor_bias" + options.fibre + ".pdf"
                #wavelength specific bias measurement
                elif options.type == "wavelength":
                    h_bias.Add(get_fibre_bias(h_bias,infile,options.wavelength))
                    canvas_string = "/data/langrock/rat-5.0-SMELLIE_analysis/scaling_factor_bias" + options.wavelength + ".pdf"
                else:
                    print "Invalid bias type!"
                    raise Exception
                    

    #draw histogram and fit with Gaussian
    c1 = ROOT.TCanvas("c1","",1)
    c1.Draw()
    h_bias.SetXTitle("s_{sim}-s_{meas}")
    h_bias.SetYTitle("Entries per bin")
    h_bias.Draw()
    f1 = ROOT.TF1("f1","gaus",-0.1,0.1)
    f1.SetLineColor(50)
    h_bias.Fit(f1,"R")

    pad = ROOT.TPad("pada","pada",0.49,0.75,0.91,0.95)
    label = ROOT.TPaveText(0.05,0.0,0.95,0.4,"br")

    pad.SetFillStyle(4000)
    pad.SetFillColor(0)
    pad.Draw()
    pad.cd()

    label.SetFillColor(0)
    label.SetLineColor(0)
    label.SetShadowColor(0)
    label.SetTextSize(0.175)
    #label.AddText("SNO+ Preliminary")
    label.AddText("MC simulations, water-filled detector")
    label.Draw()

    #save canvas
    c1.Print(canvas_string,"pdf")


