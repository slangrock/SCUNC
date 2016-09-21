#!/usr/bin/env python
import ROOT
import optparse
import os
import plot_style

def plot_tres(fibre, wl, ratio):
    """Function to plot and save the time residual histograms for all cut reagions

     Args:
      fibre (string) : fibre label
      wl (string) : wavelength fibre was run at
      ratio (string) : scattering length scaling factor applied to the run  
    """

    #open file
    infile = ROOT.TFile("/data/langrock/rat-5.0-SMELLIE_analysis/" + fibre + "/root/" + str(wl) + "_" + str(ratio) + "_data_water.root")

    #get time residual histograms for each cut region
    beam = infile.Get("time_residual_beam")
    double = infile.Get("time_residual_double")
    scatt = infile.Get("time_residual_scatt")
    avin = infile.Get("time_residual_avin")
    avout = infile.Get("time_residual_avout")
    psup = infile.Get("time_residual_psup")
    multi = infile.Get("time_residual_multi")

    c1 = ROOT.TCanvas("c1","",1)
    c1.SetLogy(1)
    c1.Draw()

    beam.SetXTitle("t_{res} (ns)")
    beam.SetYTitle("Intensity per 0.64ns")
    beam.SetLineColor(50)
    beam.Draw()

    double.SetLineColor(882)
    double.Draw("same")

    scatt.SetLineColor(1)
    scatt.Draw("same")
    
    avout.SetLineColor(20)
    avout.Draw("same")
    
    avin.SetLineColor(40)
    avin.Draw("same")

    psup.SetLineColor(30)
    psup.Draw("same")

    multi.SetLineColor(60)
    multi.Draw("same")

    l1 = ROOT.TLegend(0.40,0.6,0.85,0.89,"","brNDC")
    l1.AddEntry(beam,"In-beam hits", "l")
    l1.AddEntry(double,"Late pulses", "l")
    l1.AddEntry(avout,"Reflections off the outside of the AV", "l")
    l1.AddEntry(scatt,"Scattered events", "l")
    l1.AddEntry(avin,"Reflections off the inside of the AV", "l")
    l1.AddEntry(psup,"Reflections off the PSUP", "l")
    l1.AddEntry(multi,"Multiple effects", "l")
    l1.SetFillColor(0)
    l1.SetTextSize(0.03)
    l1.Draw("same")

    pad = ROOT.TPad("pada","pada",0.48,0.5,0.9,0.7)
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
    c1.Print("/data/langrock/rat-5.0-SMELLIE_analysis/" + fibre + "/plots/" + str(wl) + "_" + str(ratio) + "_tres_data.pdf","pdf")

def plot_angle_tres(fibre, wl, ratio):
    """Function to plot and save two-dimensional histigrams showing the beam angle       versus the time residual for each cut region.

     Args:
      fibre (string) : fibre label
      wl (string) : wavelength fibre was run at
      ratio (string) : scattering length scaling factor applied to the run  
    """

    #open file
    infile = ROOT.TFile("/data/langrock/rat-5.0-SMELLIE_analysis/" + fibre + "/root/" + str(wl) + "_" + str(ratio) + "_data_water.root")

    #get beam angle versus time residual plots for each cut region
    beam = infile.Get("angle_time_beam")
    double = infile.Get("angle_time_double")
    scatt = infile.Get("angle_time_scatt")
    avin = infile.Get("angle_time_avin")
    avout = infile.Get("angle_time_avout")
    psup = infile.Get("angle_time_psup")
    multi = infile.Get("angle_time_multi")

    c1 = ROOT.TCanvas("c1","",1)
    c1.Draw()

    beam.SetXTitle("t_{res} (ns)")
    beam.SetYTitle("#alpha (#circ)")
    beam.SetAxisRange(0,95,"Y")
    beam.SetLineColor(50)
    beam.SetMarkerColor(50)
    beam.Draw()

    double.SetLineColor(882)
    double.SetMarkerColor(882)
    double.Draw("same")

    scatt.SetLineColor(1)
    scatt.SetMarkerColor(1)
    scatt.Draw("same")
    
    avout.SetLineColor(20)
    avout.SetMarkerColor(20)
    avout.Draw("same")
    
    avin.SetLineColor(40)
    avin.SetMarkerColor(40)
    avin.Draw("same")

    psup.SetLineColor(30)
    psup.SetMarkerColor(30)
    psup.Draw("same")

    multi.SetLineColor(60)
    multi.SetMarkerColor(60)
    multi.Draw("same")

    #l1 = ROOT.TLegend(0.45,0.11,0.89,0.4,"","brNDC")
    #l1.AddEntry(beam,"In-beam hits", "l")
    #l1.AddEntry(double,"Late pulses", "l")
    #l1.AddEntry(avout,"Reflections off the outside of the AV", "l")
    #l1.AddEntry(scatt,"Scattered events", "l")
    #l1.AddEntry(avin,"Reflections off the inside of the AV", "l")
    #l1.AddEntry(psup,"Reflections off the PSUP", "l")
    #l1.AddEntry(multi,"Multiple effects", "l")
    #l1.SetFillColor(0)
    #l1.SetTextSize(0.03)
    #l1.Draw("same")

    pad = ROOT.TPad("pada","pada",0.48,0.2,0.9,0.4)
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
    c1.Print("/data/langrock/rat-5.0-SMELLIE_analysis/" + fibre + "/plots/" + str(wl) + "_" + str(ratio) + "_alpha_tres_data.png","png")

def plot_z_tres(fibre, wl, ratio):
    """Function to plot and save two-dimensional histigrams showing the z 
       coordinate versus the time residual for each cut region.

     Args:
      fibre (string) : fibre label
      wl (string) : wavelength fibre was run at
      ratio (string) : scattering length scaling factor applied to the run  
    """

    #open file
    infile = ROOT.TFile("/data/langrock/rat-5.0-SMELLIE_analysis/" + fibre + "/root/" + str(wl) + "_" + str(ratio) + "_data_water.root")

    #get z versus time residual histograms for each cut region
    beam = infile.Get("z_time_beam")
    double = infile.Get("z_time_double")
    scatt = infile.Get("z_time_scatt")
    avin = infile.Get("z_time_avin")
    avout = infile.Get("z_time_avout")
    psup = infile.Get("z_time_psup")
    multi = infile.Get("z_time_multi")

    c1 = ROOT.TCanvas("c1","",1)
    c1.Draw()

    beam.SetXTitle("t_{res} (ns)")
    beam.SetYTitle("z (mm)")
    beam.SetTitleOffset(1.1,"y")
    beam.SetLabelSize(0.03,"y")
    beam.SetLineColor(50)
    beam.SetMarkerColor(50)
    beam.Draw()

    double.SetLineColor(882)
    double.SetMarkerColor(882)
    double.Draw("same")

    scatt.SetLineColor(1)
    scatt.SetMarkerColor(1)
    scatt.Draw("same")
    
    avout.SetLineColor(20)
    avout.SetMarkerColor(20)
    avout.Draw("same")
    
    avin.SetLineColor(40)
    avin.SetMarkerColor(40)
    avin.Draw("same")

    psup.SetLineColor(30)
    psup.SetMarkerColor(30)
    psup.Draw("same")

    multi.SetLineColor(60)
    multi.SetMarkerColor(60)
    multi.Draw("same")

    l1 = ROOT.TLegend(0.40,0.6,0.85,0.89,"","brNDC")
    l1.AddEntry(beam,"In-beam hits", "l")
    l1.AddEntry(double,"Late pulses", "l")
    l1.AddEntry(avout,"Reflections off the outside of the AV", "l")
    l1.AddEntry(scatt,"Scattered events", "l")
    l1.AddEntry(avin,"Reflections off the inside of the AV", "l")
    l1.AddEntry(psup,"Reflections off the PSUP", "l")
    l1.AddEntry(multi,"Multiple effects", "l")
    l1.SetFillColor(0)
    l1.SetTextSize(0.03)
    l1.Draw("same")

    #save canvas
    c1.Print("/data/langrock/rat-5.0-SMELLIE_analysis/" + fibre + "/plots/" + str(wl) + "_" + str(ratio) + "_z_tres_data.png","png")

if __name__ == '__main__':
    """Script to plot and save different histograms for each cut region in one 
       canvas. The histograms were saved to root files by the functions in 
       apply_cuts.py. This script retreives the histograms and plots them in a
       uniform style. 
    """

    parser = optparse.OptionParser( usage = "python %prog [flags]")
    parser.add_option("-f", dest="fibre", help="Fibre to be measured")
    parser.add_option("-w", dest="wavelength", help="Wavelength of fibre")
    parser.add_option("-r", dest="ratio", help="Ratio")
    
    (options, args) = parser.parse_args()

    style = plot_style.PlotStyle()
    style.set_style()
    ROOT.gROOT.SetStyle("clearRetro")

    plot_angle_tres(options.fibre, options.wavelength, options.ratio)
    plot_z_tres(options.fibre, options.wavelength, options.ratio)
    plot_tres(options.fibre, options.wavelength, options.ratio)
