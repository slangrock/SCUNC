#!/usr/bin/env python
import ROOT
import rat
import fibre_handling
import angle_measurement
import data_corrections
import optparse
import math
import os
import fnmatch
from array import array

def get_profiles(file_name, fibre, n):
    """Function which runs over the EV branch of a RAT root file to extract the angular
       beam profile. It fills angular histograms with variable bin width determined by the 
       number of PMTs n covered by the angular bin.
     
     Args:
      file_name (string) : name of the input rat file
      fibre (string) : fibre label
      n (int) : number of PMTs to be covered by each angular bin

     Returns:
      a histogram normalised by the number of working PMTs per bin
    """
    reader = ROOT.RAT.DU.DSReader(file_name,False)

    #get fibre specific variables
    val = fibre_handling.FibreHandling(fibre)
    sourcepos, sourcedir = val.get_fibre_position()

    pmt_prop = rat.utility().GetPMTInfo()  

    #get bins of the angular histogram to be filled dependent on the number of PMTs n to be covered by each bin
    angle = angle_measurement.AngleMeasurement()
    angle.get_number_of_pmts(n)
    angle.create_pmt_angle_map(sourcepos, sourcedir)
    angle.get_bin_edges()

    #get tools for data correction
    corrections = data_corrections.DataCorrections(reader)
    pmtid = corrections.mh_corr()
    nevents = reader.GetEntryCount()          

    #define angular histograms
    h_angle = ROOT.TH1D("h_angle", "", angle.nbins, angle.bin_edges)
    h_angle_pmt = ROOT.TH2D("h_angle_pmt", "", 1000, 0, 10000, angle.nbins, angle.bin_edges)

    #start looping through file
    for ievent in range(0,nevents):
        ds = reader.GetEntry(ievent)

        #loop through events
        for iev in range(ds.GetEVCount()):
            ev = ds.GetEV(iev)

            #run over PMTs
            pmts = ev.GetCalPMTs()
            for ipmt in range(0,pmts.GetCount()):
                pmt_cal = pmts.GetPMT(ipmt)
                pmt_id = pmt_cal.GetID()

                #calculate multiple hits corrections
                P_occ = pmtid.GetBinContent(pmt_id)/nevents
                mu = 0
                if P_occ == 1:
                    mu = 0
                else:
                    mu = -math.log(1-P_occ)
                c_mh = 1 + (mu - (1 - math.exp(-mu)))

                #get pmt position and direction with respect to fibre position
                pmtpos = pmt_prop.GetPosition(pmt_id)
                pmtdir = (pmtpos - sourcepos)

                #calculate angle with respect to the fibre direction
                alpha_mc_rad = math.acos((sourcedir * pmtdir)/(sourcedir.Mag() * pmtdir.Mag()))
                alpha_mc = math.degrees(alpha_mc_rad)        

                #fill histograms
                h_angle.Fill(alpha_mc,c_mh)
                h_angle_pmt.Fill(pmt_id,alpha_mc)

    #normalise angular histogram by number of working PMTs per bin and return
    return angle.normalise_bins(h_angle_pmt, h_angle)

def get_water_profiles(file_name, fibre, n):
    """Function which runs over the EV branch of a RAT root file to extract the angular
       beam profile. It fills angular histograms with variable bin width determined by the 
       number of PMTs n covered by the angular bin. The histogram to be extracted from air-fill
       data are converted to a water-fill profile for analysis purposes.
     
     Args:
      file_name (string) : name of the input rat file
      fibre (string) : fibre label
      n (int) : number of PMTs to be covered by each angular bin

     Returns:
      a histogram normalised by the number of working PMTs per bin
    """

    reader = ROOT.RAT.DU.DSReader(file_name,False)

    #get fibre specific variables
    val = fibre_handling.FibreHandling(fibre)
    sourcepos, sourcedir = val.get_fibre_position()

    pmt_prop = rat.utility().GetPMTInfo()  

    #get bins of the angular histogram to be filled dependent on the number of PMTs n to be covered by each bin
    angle = angle_measurement.AngleMeasurement()
    angle.get_number_of_pmts(n)
    angle.create_pmt_angle_map(sourcepos, sourcedir)
    angle.get_bin_edges()

    #get tools for data correction
    corrections = data_corrections.DataCorrections(reader)
    pmtid = corrections.mh_corr()
    nevents = reader.GetEntryCount()

    #define angular histograms
    h_angle = ROOT.TH1D("h_angle", "", angle.nbins, angle.bin_edges)
    h_angle_pmt = ROOT.TH2D("h_angle_pmt", "", 1000, 0, 10000, angle.nbins, angle.bin_edges)

    #start looping through file
    for ievent in range(0,nevents):
        ds = reader.GetEntry(ievent)
        
        #loop through events
        for iev in range(ds.GetEVCount()):
            ev = ds.GetEV(iev)

            #run over PMTs
            pmts = ev.GetCalPMTs()
            for ipmt in range(0,pmts.GetCount()):
                pmt_cal = pmts.GetPMT(ipmt)
                pmt_id = pmt_cal.GetID()

                #calculate multiple hits corrections
                P_occ = pmtid.GetBinContent(pmt_id)/nevents
                mu = 0
                if P_occ == 1:
                    mu = 0
                else:
                    mu = -math.log(1-P_occ)
                c_mh = 1 + (mu - (1 - math.exp(-mu)))

                #get pmt position and direction with respect to fibre position
                pmtpos = pmt_prop.GetPosition(pmt_id)
                pmtdir = (pmtpos - sourcepos)

                #calculate angle with respect to the fibre direction
                alpha_mc_rad = math.acos((sourcedir * pmtdir)/(sourcedir.Mag() * pmtdir.Mag()))
                alpha_mc = math.degrees(alpha_mc_rad)        

                #refractive indices of all involved materials
                n_quartz = 1.46
                n_air = 1.0003
                n_water = 1.34

                #snell's law - calculate incident angle alpha_quartz from air-fill angle and then use alpha_quartz to calculate the water-fill angle alpha_water
                alpha_quartz = math.asin(math.sin(alpha_mc_rad) * n_air/n_quartz  )
                alpha_water = math.degrees(math.asin(math.sin(alpha_quartz) * n_quartz/n_water ))

                #fill histograms
                h_angle.Fill(alpha_water,c_mh)
                h_angle_pmt.Fill(pmt_id,alpha_water)

    #normalise angular histogram by number of working PMTs per bin and return
    return angle.normalise_bins(h_angle_pmt, h_angle)


def average_hist(hist):
    """Function to average and smooth multiple angular profiles
     
     Args:
      hist (list) : list of histograms to be averaged

     Returns:
      average (TH1D) : averaged histogram
      average_smooth (TH1D) : smoothed averaged histogram
    """

    #define average histograms
    average = ROOT.TH1D("average", "", 600,0,120)
    average_smooth = ROOT.TH1D("average_smooth", "", 600,0,120)

    #loop through list of histograms to be averaged
    for i in range(0,len(hist)):     
        #re-bin list content to enable adding of histograms
        hist[i].SetBins(600,0,120)
        average.Add(hist[i])
        average_smooth.Add(hist[i])

    #divide the bin content of each bin of average by the number of histograms added
    for i in range(average.GetNbinsX()):
        average.SetBinContent(i,average.GetBinContent(i)/len(hist))

    #divide the bin content of each bin of average_smooth by the number of histograms added
    for i in range(average_smooth.GetNbinsX()):
        average_smooth.SetBinContent(i,average_smooth.GetBinContent(i)/len(hist))

    #normalise histograms by their maximum
    average.Scale(1/average.GetMaximum())
    average_smooth.Scale(1/average.GetMaximum())
    average_smooth.Smooth(10)
    
    return average, average_smooth

if __name__ == '__main__':
    """Script runs through all the given files to measure the angular profile of the SMELLIE fibres.
       Dependent on the medium chosen by the -m the profile is measured straight from the files or 
       converted to a water-fill profile assuming air-fill data. The -p flag determines how many PMTs
       are covered by each bin for the measurement. The measured histograms are saved to a root file.
    """

    parser = optparse.OptionParser( usage = "python %prog [flags]")
    parser.add_option("-m", dest="medium", help="Name of the inner detector fill medium")
    parser.add_option("-p", dest="n_pmt", help="Number of PMTs covered by each angle")

    (options, args) = parser.parse_args()

    if not options.medium:
        print 'Need detector medium!'
        raise Exception
    if not options.n_pmt:
        print 'Need number of PMTs!'
        raise Exception

    #take profile from root files without conversions. For wach fibre-wavelength combination the histograms are saved to a root file and the averaged profiles are saved to a separate root file
    if options.medium == "air":
        hist_list = []

        #profile from FS155, 407nm
        hist = get_profiles("/data/snoplus/slangrock/root_files/SMELLIE/data/FS155_407_100_7112.root", "FS155", int(options.n_pmt))
        hist_list.append(hist)
        outputroot = ROOT.TFile("/data/langrock/rat-5.0-SMELLIE_analysis/FS155/root/407_air_" + str(options.n_pmt) + "_angle_profile.root","recreate")
        hist.Write()
        outputroot.Close()

        #profile from FS037, 407nm
        hist = get_profiles("/data/snoplus/slangrock/root_files/SMELLIE/data/FS037_407_100_7116.root", "FS037", int(options.n_pmt))
        hist_list.append(hist)
        outputroot = ROOT.TFile("/data/langrock/rat-5.0-SMELLIE_analysis/FS037/root/407_air_" + str(options.n_pmt) + "_angle_profile.root","recreate")
        hist.Write()
        outputroot.Close()

        #profile from FS237, 446nm
        hist = get_profiles("/data/snoplus/slangrock/root_files/SMELLIE/data/FS237_446_28_7179.root", "FS237", int(options.n_pmt))
        hist_list.append(hist)
        outputroot = ROOT.TFile("/data/langrock/rat-5.0-SMELLIE_analysis/FS237/root/446_air_" + str(options.n_pmt) + "_angle_profile.root","recreate")
        hist.Write()
        outputroot.Close()

        #averaged profiles
        outputroot_average = ROOT.TFile("/data/langrock/rat-5.0-SMELLIE_analysis/" + str(options.medium) + "_" + str(options.n_pmt) + "_angle_profile_average.root","recreate")            
        average, average_smooth = average_hist(hist_list)
        average.Write()
        average_smooth.Write()
        outputroot_average.Close()

    #take profile from root files and convert them to water-fill profiles. For wach fibre-wavelength combination the histograms are saved to a root file and the averaged profiles are saved to a separate root file
    if options.medium == "water":
        hist_list = []

        #profile from FS155, 407nm
        hist = get_water_profiles("/data/snoplus/slangrock/root_files/SMELLIE/data/FS155_407_100_7112.root", "FS155", int(options.n_pmt))
        hist_list.append(hist)
        outputroot = ROOT.TFile("/data/langrock/rat-5.0-SMELLIE_analysis/FS155/root/407_water_" + str(options.n_pmt) + "_angle_profile_7112.root","recreate")
        hist.Write()
        outputroot.Close()

        #profile from FS155, 495nm
        hist = get_water_profiles("/data/snoplus/slangrock/root_files/SMELLIE/data/FS155_495_85_7128.root", "FS155", int(options.n_pmt))
        hist_list.append(hist)
        outputroot = ROOT.TFile("/data/langrock/rat-5.0-SMELLIE_analysis/FS155/root/495_water_" + str(options.n_pmt) + "_angle_profile_7128.root","recreate")
        hist.Write()
        outputroot.Close()

        #profile from FS155, 446nm
        hist = get_water_profiles("/data/snoplus/slangrock/root_files/SMELLIE/data/FS155_446_100_7141.root", "FS155", int(options.n_pmt))
        hist_list.append(hist)
        outputroot = ROOT.TFile("/data/langrock/rat-5.0-SMELLIE_analysis/FS155/root/446_water_" + str(options.n_pmt) + "_angle_profile_7141.root","recreate")
        hist.Write()
        outputroot.Close()             

        #averaged profiles
        outputroot_average = ROOT.TFile("/data/langrock/rat-5.0-SMELLIE_analysis/" + str(options.medium) + "_" + str(options.n_pmt) + "_angle_profile_average_fs155.root","recreate")            
        average, average_smooth = average_hist(hist_list)
        average.Write()
        average_smooth.Write()
        outputroot_average.Close()

