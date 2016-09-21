#!/usr/bin/env python
import ROOT
import rat
import math
import fibre_handling
from array import array

class AngleMeasurement:
    """Class to define the bins for the histograms used for angular profile
       measurements.

     Creates the number of bins and upper bin edges for histograms used to measure
     the angular profiles of the fibres. The amount of PMTs covered by each bin
     is specified by the user. The bin width is dependent on the fibre which is 
     measured.

     Attributes:
      variables shared by all instances:
       pmt_prop (RAT PMT properties) : data base of the protperties of the SNO+ PMTs

      variables unique to each instance:
       pmt (int) : number of PMTs
       nbins (int) : number of bins
       pmt_angle_map (list) : list of two element lists containing [angle, pmt], with angle being
                              the angle of the PMT with respect to the fibre direction and pmt
                              being the PMT number
       n (int) : variable to define how many PMTs should be covered by each angular bin
       bin_edges (array) : array to define the upper angular bin edges
    """

    pmt_prop = rat.utility().GetPMTInfo()

    def __init__(self):
        self.pmt = 0
        self.nbins = 0
        self.pmt_angle_map = []
        self.n = 0
        self.bin_edges = array('d')   


    def get_number_of_pmts(self,n):
        """Function to get number of PMTs and calculate the number of bins of 
           the histogram.

         Args:
          n (int) : number of PMTs which should be covered by each bin

         Returns:
          pmt (int) : number of PMTs
          nbins (int) : number of bins
        """
        self.n = n
        for ipmt in range(0,self.pmt_prop.GetCount()):
            pmt_type = self.pmt_prop.GetType(ipmt)
            if pmt_type == 1:
                self.pmt += 1

        self.nbins = self.pmt/n
        return self.pmt, self.nbins

    def create_pmt_angle_map(self, sourcepos, sourcedir):
        """Function to calculate the angle of each PMT

         Args: 
           sourcepos (TVector3) : vector of the position of the fibre in the detector
           sourcedir (TVector3) : direction vector the fibre is pointed at

         Returns:
           pmt_angle_map (list) : list assigning a PMT to each angle
        """
        #loop through all pmts
        for ipmt in range(0,self.pmt+1):
            pmtpos = self.pmt_prop.GetPosition(ipmt)
            pmtdir = pmtpos - sourcepos

            #get angle of PMT
            alpha_mc_rad = math.acos((sourcedir * pmtdir)/(sourcedir.Mag() * pmtdir.Mag()))
            alpha_mc = math.degrees(alpha_mc_rad) 

            #create a two object list of angle and pmt
            lcn_angle = [alpha_mc, ipmt]           
            self.pmt_angle_map.append(lcn_angle)
    
        #sort list of two object lists by increasing angle
        self.pmt_angle_map = sorted(self.pmt_angle_map,key=lambda x: (x[0]))
        return self.pmt_angle_map

    def get_bin_edges(self):
        """Function to define the variable angular bin edges of the histograms
        
         Returns:
           bin_edges (array) : array of bin edges

         Creates self.bin_edges dependend on the number of PMTs self.n covered by 
         each angle from the list self.pmt_angle_map
        """
        #set first entry to 0
        self.bin_edges.append(0.0)
        #set PMT counter
        num_pmt = 1
        #loop through all entries in the list of angles and pmts
        for i in range(0,len(self.pmt_angle_map)+1):
            #if PMT counter smaller than required number of PMTS, add 1
            if num_pmt < self.n:
                num_pmt +=1
            #if PMT counter equal to required number of PMTs, reset counter
            elif num_pmt==self.n:
                num_pmt = 1
                #add equivalent angle to the bin_edges array
                bin_edge = self.pmt_angle_map[i][0]
                self.bin_edges.append(bin_edge)

        return self.bin_edges

    def normalise_bins(self, angle_pmt_hist, angle_hist):
        """Function to normalise the angle histogram by number of working PMTs in bin

         Args:
           angle_pmt_hist (TH2D) : 2D histogram angle vs PMT
           angle_hist (TH1D) : 1D angle histogram

         Returns: 
           angle_hist (TH1D) : normalised 1D angle histogram
        """

        #loop through all angles in 2D histogram
        for i in range(0,angle_pmt_hist.GetNbinsY()):
            #set PMT counter
            pmtcontent = 0
            #for each angle scan through all PMTs
            for j in range(0,angle_pmt_hist.GetNbinsX()):
                #if entry is found add 1 to counter
                if angle_pmt_hist.GetBinContent(j,i) != 0:
                    pmtcontent += 1
            #if counter is not 0, normalise hist content at this angle by number found by PMT counter
            if pmtcontent != 0:
                angle_hist.SetBinContent(i,angle_hist.GetBinContent(i)/pmtcontent)

        return angle_hist

