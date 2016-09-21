#!/usr/bin/env python
import ROOT
import optparse
import os
import plot_style

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

def plot_average():
    """Function to plot and save averaged angular profiles (smoothed and unsmoothed in same canvas)
    """

    #get histograms
    file_average = ROOT.TFile("/data/langrock/rat-5.0-SMELLIE_analysis/water_5_angle_profile_average.root")
    average = file_average.Get("average")
    average_smooth = file_average.Get("average_smooth")

    c1 = ROOT.TCanvas("c1","c1", 1)
    c1.Draw()

    pad = ROOT.TPad("pad","pad",0.3,0.45,0.89,0.65)
    label = ROOT.TPaveText(0.05,0.0,0.95,0.4,"br")

    average = norm_hist(average)
    average_smooth = norm_hist(average_smooth)

    average.SetXTitle("#alpha (#circ)")
    average.SetYTitle("Intensity per 0.2#circ")
    average.SetLineWidth(2.2)
    average.SetLineColor(50)
    average.SetAxisRange(0,15,"X")
    average.Draw()

    average_smooth.SetLineColor(30)
    average_smooth.SetLineWidth(2.2)
    average_smooth.Draw("same")
    
    l = ROOT.TLegend(0.4,0.7,0.88,0.87,"","brNDC")
    l.AddEntry(average,"Averaged profile", "l")
    l.AddEntry(average_smooth,"Smooth averaged profile", "l")
    l.SetFillColor(0)
    l.SetTextSize(0.04)
    l.Draw("same")

    pad.SetFillStyle(4000)
    pad.SetFillColor(0)
    pad.Draw()
    pad.cd()

    label.SetFillColor(0)
    label.SetLineColor(0)
    label.SetShadowColor(0)
    label.SetTextSize(0.175)
    #label.AddText("Data from 20/03/2014, air-filled detector")
    label.AddText("Data from 20/03/2014 converted to water-filled detector")
    label.Draw()

    #save canvas
    c1.Print("/data/langrock/rat-5.0-SMELLIE_analysis/average_angle_profile.pdf","pdf")

def plot_water_air_comp(fibre, wl):
    """Function to plot the air-fill and water-fill profiles for a given fibre-wavelength combination

     Args:
      fibre (string) : fibre label
      wl (string) : wavelength fibre was run at
    """
    
    #retreive air and water profiles
    file_air = ROOT.TFile("/data/langrock/rat-5.0-SMELLIE_analysis/" + fibre + "/root/" + str(wl) + "_air_5_angle_profile.root")
    angle_air = file_air.Get("h_angle")
    file_water = ROOT.TFile("/data/langrock/rat-5.0-SMELLIE_analysis/" + fibre + "/root/" + str(wl) + "_water_5_angle_profile.root")
    angle_water = file_water.Get("h_angle")

    c1 = ROOT.TCanvas("c1","c1", 1)
    c1.Draw()

    pad = ROOT.TPad("pad","pad",0.3,0.45,0.89,0.65)
    label = ROOT.TPaveText(0.05,0.0,0.95,0.4,"br")

    angle_air = norm_hist(angle_air)
    angle_water = norm_hist(angle_water)

    angle_air.SetXTitle("#alpha (#circ)")
    angle_air.SetYTitle("Intensity per 5 PMTs")
    angle_air.SetLineWidth(2.2)
    angle_air.SetLineColor(50)
    angle_air.SetAxisRange(0,15,"X")
    angle_air.Draw()

    angle_water.SetLineColor(30)
    angle_water.SetLineWidth(2.2)
    angle_water.Draw("same")
    
    l = ROOT.TLegend(0.4,0.7,0.88,0.87,"","brNDC")
    l.AddEntry(angle_air,"Profile extracted from data", "l")
    l.AddEntry(angle_water,"Profile converted to water-fill", "l")
    l.SetFillColor(0)
    l.SetTextSize(0.04)
    l.Draw("same")

    pad.SetFillStyle(4000)
    pad.SetFillColor(0)
    pad.Draw()
    pad.cd()

    label.SetFillColor(0)
    label.SetLineColor(0)
    label.SetShadowColor(0)
    label.SetTextSize(0.175)
    label.AddText("Data from 20/03/2014, air-filled detector")
    label.AddText("Data from 20/03/2014 converted to water-filled detector")
    label.Draw()

    #save canvas
    c1.Print("/data/langrock/rat-5.0-SMELLIE_analysis/" + fibre + "/plots/" + str(wl) + "_air_water_comp_angle_profile.pdf","pdf")

def plot_diff_runs(fibre):
    """Function to plot angular profiles extracted from different runs for a given fibre

     Args:
      fibre (string) : fibre label
    """

    #get profiles from files
    file_407 = ROOT.TFile("/data/langrock/rat-5.0-SMELLIE_analysis/" + fibre + "/root/407_water_5_angle_profile_7112.root")
    angle_407 = file_407.Get("h_angle")
    file_495 = ROOT.TFile("/data/langrock/rat-5.0-SMELLIE_analysis/" + fibre + "/root/495_water_5_angle_profile_7128.root")
    angle_495 = file_495.Get("h_angle")
    file_446 = ROOT.TFile("/data/langrock/rat-5.0-SMELLIE_analysis/" + fibre + "/root/446_water_5_angle_profile_7141.root")
    angle_446 = file_446.Get("h_angle")

    c1 = ROOT.TCanvas("c1","c1", 1)
    c1.Draw()

    pad = ROOT.TPad("pad","pad",0.3,0.45,0.89,0.65)
    label = ROOT.TPaveText(0.05,0.0,0.95,0.4,"br")

    angle_407 = norm_hist(angle_407)
    angle_495 = norm_hist(angle_495)
    angle_446 = norm_hist(angle_446)

    angle_407.SetXTitle("#alpha (#circ)")
    angle_407.SetYTitle("Intensity per 5 PMTs")
    angle_407.SetLineWidth(2.2)
    angle_407.SetLineColor(50)
    angle_407.SetAxisRange(0,15,"X")
    angle_407.Draw()

    angle_495.SetLineColor(30)
    angle_495.SetLineWidth(2.2)
    angle_495.Draw("same")

    angle_446.SetLineColor(60)
    angle_446.SetLineWidth(2.2)
    angle_446.Draw("same")
    
    l = ROOT.TLegend(0.35,0.7,0.88,0.87,"","brNDC")
    l.AddEntry(angle_407,"407nm at 100% laser intensity", "l")
    l.AddEntry(angle_495,"495nm at 85% laser intensity", "l")
    l.AddEntry(angle_446,"446nm at 100% laser intensity", "l")
    l.SetFillColor(0)
    l.SetTextSize(0.04)
    l.Draw("same")

    pad.SetFillStyle(4000)
    pad.SetFillColor(0)
    pad.Draw()
    pad.cd()

    label.SetFillColor(0)
    label.SetLineColor(0)
    label.SetShadowColor(0)
    label.SetTextSize(0.175)
    #label.AddText("Data from 20/03/2014, air-filled detector")
    label.AddText("Data from 20/03/2014 converted to water-filled detector")
    label.Draw()

    #save canvas
    c1.Print("/data/langrock/rat-5.0-SMELLIE_analysis/" + fibre + "/plots/water_comp_angle_profile_diff_runs.pdf","pdf")

def plot_diff_pmt(fibre, wl, pmt1, pmt2):
    """Function to plot profiles with different bin widths for a given fibre-wavelength combination, dependent
       on the number of PMTs covered by each bin

     Args:
      fibre (string) : fibre label
      wl (string) : wavelength fibre was run at
      pmt1 (int) : number of PMTs covered per angular bin for first histogram
      pmt2 (int) : number of PMTs covered per angular bin for second histogram     
    """

    #get histograms from file
    file_5 = ROOT.TFile("/data/langrock/rat-5.0-SMELLIE_analysis/" + fibre + "/root/" + str(wl) + "_air_" + str(pmt1) + "_angle_profile.root")
    angle_5 = file_5.Get("h_angle")
    file_10 = ROOT.TFile("/data/langrock/rat-5.0-SMELLIE_analysis/" + fibre + "/root/" + str(wl) + "_air_" + str(pmt2) + "_angle_profile.root")
    angle_10 = file_10.Get("h_angle")

    c1 = ROOT.TCanvas("c1","c1", 1)
    c1.Draw()

    pad = ROOT.TPad("pad","pad",0.3,0.45,0.89,0.65)
    label = ROOT.TPaveText(0.05,0.0,0.95,0.4,"br")

    angle_5 = norm_hist(angle_5)
    angle_10 = norm_hist(angle_10)

    angle_5.SetXTitle("#alpha (#circ)")
    angle_5.SetYTitle("Intensity per n PMTs")
    angle_5.SetLineWidth(2.2)
    angle_5.SetLineColor(50)
    angle_5.SetAxisRange(0,15,"X")
    angle_5.Draw()

    angle_10.SetLineColor(30)
    angle_10.SetLineWidth(2.2)
    angle_10.Draw("same")
    
    l = ROOT.TLegend(0.4,0.7,0.87,0.87,"","brNDC")
    l.AddEntry(angle_5,str(pmt1) + " PMTs per angular bin", "l")
    l.AddEntry(angle_10,str(pmt2) + " PMTs per angular bin", "l")
    l.SetFillColor(0)
    l.SetTextSize(0.04)
    l.Draw("same")

    pad.SetFillStyle(4000)
    pad.SetFillColor(0)
    pad.Draw()
    pad.cd()

    label.SetFillColor(0)
    label.SetLineColor(0)
    label.SetShadowColor(0)
    label.SetTextSize(0.175)
    label.AddText("Data from 20/03/2014, air-filled detector")
    label.Draw()

    #save canvas
    c1.Print("/data/langrock/rat-5.0-SMELLIE_analysis/" + fibre + "/plots/" + str(wl) + "_air_" + str(pmt1) + "_" + str(pmt2) + "_angle_profile.pdf","pdf")

def extract_fibre_profiles(fibre, wl, pmt1):
    """Function to save the measured angular profiles into probability arrays as used in SMELLIE.ratdb
       for a given fibre-wavelength combination and a certain bin width.

     Args:
      fibre (string) : fibre label
      wl (string) : wavelength fibre was run at 
      pmt1 (int) : number of PMTs covered per angular bin for the histogram
    """

    #get histogram from file
    file_5 = ROOT.TFile("/data/langrock/rat-5.0-SMELLIE_analysis/" + fibre + "/root/" + str(wl) + "_air_" + str(pmt1) + "_angle_profile.root")
    angle = file_5.Get("h_angle")

    angle = norm_hist(angle)

    outputfile = open("/data/langrock/rat-5.0-SMELLIE_analysis/" + fibre + "/angular_profile.txt","w")

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

def extract_average_profiles():
    """Function to save averaged smooth profile into probability arrays as used in SMELLIE.ratdb
    """

    #get histogram from file
    file_average = ROOT.TFile("/data/langrock/rat-5.0-SMELLIE_analysis/water_5_angle_profile_average.root")
    angle = file_average.Get("average_smooth")

    angle = norm_hist(angle)

    outputfile = open("/data/langrock/rat-5.0-SMELLIE_analysis/water_5_angular_profile_average.txt","w")

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
    """Script to plot the angular profiles measured using the get_angle_profiles.py. All plots
       are automatically saved.
    """

    parser = optparse.OptionParser( usage = "python %prog [flags]")
    parser.add_option("-f", dest="fibre", help="Fibre to be measured")
    parser.add_option("-w", dest="wavelength", help="Wavelength of fibre")
    
    (options, args) = parser.parse_args()

    style = plot_style.PlotStyle()
    style.set_style()
    ROOT.gROOT.SetStyle("clearRetro")

    #plot comparison between different binned histograms
    plot_diff_pmt(options.fibre, options.wavelength, 5, 10)
    plot_diff_pmt(options.fibre, options.wavelength, 5, 7)
    plot_diff_pmt(options.fibre, options.wavelength, 5, 4)

    #plot comparison between air and water profiles
    plot_water_air_comp(options.fibre, options.wavelength)

    #plot averaged profiles (smoothed and unsmoothed)
    plot_average()

    #plot comparison between different runs
    plot_diff_runs(options.fibre)

    #extract the probability arrays from the different measured profiles as used by SMELLIE.ratdb
    extract_fibre_profiles(options.fibre, options.wavelength, 5)
    extract_average_profiles()
