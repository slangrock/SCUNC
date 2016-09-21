#!/usr/bin/env python
import ROOT
import plot_style
import optparse
import os

def norm_hist(hist):
    """Function to normalise a histogram to its maximum

     Args:
      hist (TH1D) : one-dimensional histogram

     Returns:
      normalised histogram
    """
    norm = hist.GetMaximum()
    hist.Scale(1/norm)
    return hist

def fit_profile():
    """Function to fit the measured angular profile. Applies a Gaussian fit to 
       the bulk of the profile and an exponential fit to the tail of the profile.
       Both fits are combined to determine a function which describes the full 
       angular profile. 

     Returns:
      the function containing the combined fits 
    """

    #retrieve angular profile from file 
    file_average = ROOT.TFile("/data/langrock/rat-5.0-SMELLIE_analysis/water_5_angle_profile_average.root")
    angle = file_average.Get("average_smooth")

    #define the Gaussian and exponential fit functions and the combined function
    f1 = ROOT.TF1("f1","gaus", 0, 1.5)
    f2 = ROOT.TF1("f2","expo", 1.0, 10)
    total = ROOT.TF1("total","gaus(0)+expo(3)", 0, 10)

    #draw profile and fit with Gaussian and exponential
    c1 = ROOT.TCanvas("c1","",1)
    c1.Draw()
    angle.SetXTitle("#alpha (#circ)")
    angle.SetAxisRange(0,10,"X")
    angle.Draw()
    f1.SetLineColor(30)
    angle.Fit(f1,"R")
    f2.SetLineColor(50)
    angle.Fit(f2,"R+")

    pad = ROOT.TPad("pad","pad",0.3,0.45,0.89,0.65)
    label = ROOT.TPaveText(0.05,0.0,0.95,0.4,"br")
    pad.SetFillStyle(4000)
    pad.SetFillColor(0)
    pad.Draw()
    pad.cd()

    label.SetFillColor(0)
    label.SetLineColor(0)
    label.SetShadowColor(0)
    label.SetTextSize(0.175)
    label.AddText("Generic angular profile, water-filled detector")
    label.Draw()

    #save canvas
    c1.Print("/data/langrock/rat-5.0-SMELLIE_analysis/angular_profile_separate_fits.pdf","pdf")

    #get fit parameters from the Gaussian and exponential and set them as parameters of the combined function
    par0 = f1.GetParameter(0)
    par1 = f1.GetParameter(1)
    par2 = f1.GetParameter(2)
    par3 = f2.GetParameter(0)
    par4 = f2.GetParameter(1)

    total.SetParameter(0,par0)
    total.SetParameter(1,par1)
    total.SetParameter(2,par2)
    total.SetParameter(3,par3)
    total.SetParameter(4,par4)

    #draw profile and combined function
    c2 = ROOT.TCanvas("c2","",1)
    c2.Draw()
    angle.Draw()
    total.SetLineColor(60)
    angle.Fit(total,"R")

    pad.Draw()
    pad.cd()
    label.Draw()

    #save canvas
    c2.Print("/data/langrock/rat-5.0-SMELLIE_analysis/angular_profile_combined_fits.pdf","pdf")

    return total

def fluctuate_fit(total, fluctuate):
    """Function to fluctuate the width and the tail of the combined function by a 
       factor fluctuate. The width and tail contributions are both increased and 
       decreased by defining new separate functions, which are ultimately 
       converted back to histograms.

     Args:
      total (TF1) : combined function returned by fit_profile()
      fluctuate (float) : factor the profile is altered by

     Returns:
      two altered angular profiles in histogram format
    """

    #define function to increase width and tail contributions of the profile
    total_up = ROOT.TF1("total_up","gaus(0)+expo(3)", 0, 10)
    total_up.SetParameter(0,total.GetParameter(0))
    total_up.SetParameter(1,total.GetParameter(1))
    total_up.SetParameter(2,total.GetParameter(2)*fluctuate)
    total_up.SetParameter(3,total.GetParameter(3))
    total_up.SetParameter(4,total.GetParameter(4)*1/fluctuate)

    #define function to decrease width and tail contributions of the profile
    total_down = ROOT.TF1("total_down","gaus(0)+expo(3)", 0, 10)
    total_down.SetParameter(0,total.GetParameter(0))
    total_down.SetParameter(1,total.GetParameter(1))
    total_down.SetParameter(2,total.GetParameter(2)*1/fluctuate)
    total_down.SetParameter(3,total.GetParameter(3))
    total_down.SetParameter(4,total.GetParameter(4)*fluctuate)

    #define histograms and fill them with the combined and the two altered functions
    angle = ROOT.TH1D()
    angle.SetBins(75,0,15)
    angle.FillRandom("total",5000)

    angle_up = ROOT.TH1D()
    angle_up.SetBins(75,0,15)
    angle_up.FillRandom("total_up",5000)

    angle_down = ROOT.TH1D()
    angle_down.SetBins(75,0,15)
    angle_down.FillRandom("total_down",5000)

    #normalise histograms
    angle = norm_hist(angle)
    angle_up = norm_hist(angle_up)
    angle_down = norm_hist(angle_down)

    #draw histograms
    c1 = ROOT.TCanvas("c1","",1)
    c1.Draw()
    angle.SetXTitle("#alpha (#circ)")
    angle.SetYTitle("Intensity per 0.2#circ")
    angle.Draw()
    angle_up.SetLineColor(30)
    angle_up.Draw("same")
    angle_down.SetLineColor(50)
    angle_down.Draw("same")

    l1 = ROOT.TLegend(0.40,0.6,0.85,0.8,"","brNDC")
    l1.AddEntry(angle,"Standard profile", "l")
    l1.AddEntry(angle_up,"Widened profile", "l")
    l1.AddEntry(angle_down,"Narrower profile", "l")
    l1.SetFillColor(0)
    l1.SetTextSize(0.03)
    l1.Draw("same")

    pad = ROOT.TPad("pad","pad",0.3,0.45,0.89,0.65)
    label = ROOT.TPaveText(0.05,0.0,0.95,0.4,"br")
    pad.SetFillStyle(4000)
    pad.SetFillColor(0)
    pad.Draw()
    pad.cd()

    label.SetFillColor(0)
    label.SetLineColor(0)
    label.SetShadowColor(0)
    label.SetTextSize(0.175)
    label.AddText("Generic angular profile, water-filled detector")
    label.Draw()

    #save canvas
    c1.Print("/data/langrock/rat-5.0-SMELLIE_analysis/angular_profile_systematics_" + str(fluctuate) + ".pdf","pdf")

    return angle_up, angle_down

def fluctuate_width(total, fluctuate):
    """Function to fluctuate the width of the combined function by a factor 
       fluctuate. The width contributions are both increased and decreased by 
       defining new separate functions, which are ultimately converted back to 
       histograms.

     Args:
      total (TF1) : combined function returned by fit_profile()
      fluctuate (float) : factor the profile is altered by

     Returns:
      two altered angular profiles in histogram format
    """

    #define function to increase the width contributions of the profile
    total_up = ROOT.TF1("total_up","gaus(0)+expo(3)", 0, 10)
    total_up.SetParameter(0,total.GetParameter(0))
    total_up.SetParameter(1,total.GetParameter(1))
    total_up.SetParameter(2,total.GetParameter(2)*fluctuate)
    total_up.SetParameter(3,total.GetParameter(3))
    total_up.SetParameter(4,total.GetParameter(4))

    #define function to decrease the width contributions of the profile
    total_down = ROOT.TF1("total_down","gaus(0)+expo(3)", 0, 10)
    total_down.SetParameter(0,total.GetParameter(0))
    total_down.SetParameter(1,total.GetParameter(1))
    total_down.SetParameter(2,total.GetParameter(2)*1/fluctuate)
    total_down.SetParameter(3,total.GetParameter(3))
    total_down.SetParameter(4,total.GetParameter(4))

    #define histograms and fill them with the combined and the two altered functions
    angle = ROOT.TH1D()
    angle.SetBins(75,0,15)
    angle.FillRandom("total",5000)

    angle_up = ROOT.TH1D()
    angle_up.SetBins(75,0,15)
    angle_up.FillRandom("total_up",5000)

    angle_down = ROOT.TH1D()
    angle_down.SetBins(75,0,15)
    angle_down.FillRandom("total_down",5000)

    #normalise histograms
    angle = norm_hist(angle)
    angle_up = norm_hist(angle_up)
    angle_down = norm_hist(angle_down)

    #draw histograms
    c1 = ROOT.TCanvas("c1","",1)
    c1.Draw()
    angle.SetXTitle("#alpha (#circ)")
    angle.SetYTitle("Intensity per 0.2#circ")
    angle.SetAxisRange(0,5,"x")
    angle.Draw()
    angle_up.SetLineColor(30)
    angle_up.Draw("same")
    angle_down.SetLineColor(50)
    angle_down.Draw("same")

    l1 = ROOT.TLegend(0.40,0.6,0.85,0.8,"","brNDC")
    l1.AddEntry(angle,"Standard profile", "l")
    l1.AddEntry(angle_up,"Widened profile", "l")
    l1.AddEntry(angle_down,"Narrower profile", "l")
    l1.SetFillColor(0)
    l1.SetTextSize(0.03)
    l1.Draw("same")

    pad = ROOT.TPad("pad","pad",0.3,0.45,0.89,0.65)
    label = ROOT.TPaveText(0.05,0.0,0.95,0.4,"br")
    pad.SetFillStyle(4000)
    pad.SetFillColor(0)
    pad.Draw()
    pad.cd()

    label.SetFillColor(0)
    label.SetLineColor(0)
    label.SetShadowColor(0)
    label.SetTextSize(0.175)
    label.AddText("Generic angular profile, water-filled detector")
    label.Draw()

    #save canvas
    c1.Print("/data/langrock/rat-5.0-SMELLIE_analysis/angular_profile_systematics_" + str(fluctuate) + "_width.pdf","pdf")

    return angle_up, angle_down

def fluctuate_tail(total, fluctuate):
    """Function to fluctuate the tail of the combined function by a factor 
       fluctuate. The tail contributions are both increased and decreased by 
       defining new separate functions, which are ultimately converted back to 
       histograms.

     Args:
      total (TF1) : combined function returned by fit_profile()
      fluctuate (float) : factor the profile is altered by

     Returns:
      two altered angular profiles in histogram format
    """

    #define function to increase the tail contributions of the profile
    total_up = ROOT.TF1("total_up","gaus(0)+expo(3)", 0, 10)
    total_up.SetParameter(0,total.GetParameter(0))
    total_up.SetParameter(1,total.GetParameter(1))
    total_up.SetParameter(2,total.GetParameter(2))
    total_up.SetParameter(3,total.GetParameter(3))
    total_up.SetParameter(4,total.GetParameter(4)*1/fluctuate)

    #define function to decrease the tail contributions of the profile
    total_down = ROOT.TF1("total_down","gaus(0)+expo(3)", 0, 10)
    total_down.SetParameter(0,total.GetParameter(0))
    total_down.SetParameter(1,total.GetParameter(1))
    total_down.SetParameter(2,total.GetParameter(2))
    total_down.SetParameter(3,total.GetParameter(3))
    total_down.SetParameter(4,total.GetParameter(4)*fluctuate)

    #define histograms and fill them with the combined and the two altered functions
    angle = ROOT.TH1D()
    angle.SetBins(75,0,15)
    angle.FillRandom("total",5000)

    angle_up = ROOT.TH1D()
    angle_up.SetBins(75,0,15)
    angle_up.FillRandom("total_up",5000)

    angle_down = ROOT.TH1D()
    angle_down.SetBins(75,0,15)
    angle_down.FillRandom("total_down",5000)

    #normalise histograms
    angle = norm_hist(angle)
    angle_up = norm_hist(angle_up)
    angle_down = norm_hist(angle_down)

    #draw histograms
    c1 = ROOT.TCanvas("c1","",1)
    c1.Draw()
    angle.SetXTitle("#alpha (#circ)")
    angle.SetYTitle("Intensity per 0.2#circ")
    angle.SetAxisRange(0,10,"x")
    angle.SetAxisRange(0,0.2,"y")
    angle.SetTitleOffset(0.95,"y")
    angle.Draw()
    angle_up.SetLineColor(30)
    angle_up.Draw("same")
    angle_down.SetLineColor(50)
    angle_down.Draw("same")

    l1 = ROOT.TLegend(0.40,0.6,0.85,0.8,"","brNDC")
    l1.AddEntry(angle,"Standard profile", "l")
    l1.AddEntry(angle_up,"Higher tail contributions", "l")
    l1.AddEntry(angle_down,"Smaller tail contributions", "l")
    l1.SetFillColor(0)
    l1.SetTextSize(0.03)
    l1.Draw("same")

    pad = ROOT.TPad("pad","pad",0.3,0.45,0.89,0.65)
    label = ROOT.TPaveText(0.05,0.0,0.95,0.4,"br")
    pad.SetFillStyle(4000)
    pad.SetFillColor(0)
    pad.Draw()
    pad.cd()

    label.SetFillColor(0)
    label.SetLineColor(0)
    label.SetShadowColor(0)
    label.SetTextSize(0.175)
    label.AddText("Generic angular profile, water-filled detector")
    label.Draw()

    #save canvas
    c1.Print("/data/langrock/rat-5.0-SMELLIE_analysis/angular_profile_systematics_" + str(fluctuate) + "_tail.pdf","pdf")

    return angle_up, angle_down

def extract_fibre_profiles(angle, fluct_type):
    """Function to save the altered angular profiles into probability arrays as used in SMELLIE.ratdb.

     Args:
      angle (TH1D) : histogram of the profile to be extracted
      fluct_type (string) : type of profile alteration applied to the profile
    """

    #define txt file the extracted profile is saved to
    outputfile = open("/data/langrock/rat-5.0-SMELLIE_analysis/angular_profile_systematics_" + fluct_type + ".txt","w")

    #write intensities per bin into array for all angles below 15 degrees
    for i in range(0,angle.GetNbinsX()):
        if angle.GetBinLowEdge(i)>=15:
            continue

        if angle.GetBinContent(i) > 0:
            outputfile.write(str(angle.GetBinContent(i)))
            outputfile.write(", ")

    outputfile.write("\n \n")

    #write accompanying angle for each intensity into array for all angles below 15 degrees
    for i in range(0,angle.GetNbinsX()):
        if angle.GetBinLowEdge(i)>=15:
            continue

        if angle.GetBinContent(i) > 0:
            if angle.GetBinCenter(i) < 0:
                outputfile.write("0.0, ")

            else:
                outputfile.write(str(angle.GetBinCenter(i)))
                outputfile.write(", ")

    outputfile.close()

if __name__ == '__main__':
    """Script to fluctute the measured angluar profile for systematic considerations.
       Possible alterations include the profile width, the profile tail or both the 
       tail and the width.
    """
    

    parser = optparse.OptionParser( usage = "python %prog [flags]")
    parser.add_option("-f", dest="factor", help="factor to fluctate profile by")
    parser.add_option("-t", dest="type", help="fluctuation type, can be 'width', 'tail' or 'both'")

    
    (options, args) = parser.parse_args()

    style = plot_style.PlotStyle()
    style.set_style()
    ROOT.gROOT.SetStyle("clearRetro")

    #get combined function to describe angular profile
    total = fit_profile()

    #alter width only
    if options.type == "width":
        width_up, width_down = fluctuate_width(total, float(options.factor))
        extract_fibre_profiles(width_up,"width_up")
        extract_fibre_profiles(width_down,"width_down")
    
    #alter tail only
    elif options.type == "tail":
        tail_up, tail_down = fluctuate_tail(total, float(options.factor))
        extract_fibre_profiles(tail_up,"tail_up")
        extract_fibre_profiles(tail_down,"tail_down")
    
    #alter both tail and width
    elif options.type == "both":
        up, down = fluctuate_fit(total, float(options.factor))
        extract_fibre_profiles(up,"up")
        extract_fibre_profiles(down,"down")

    else:
        print "Invalid systematic type!"
        raise Exception
                
