#!/usr/bin/env python
import ROOT
import rat
import math

class DataCorrections:
    """Class to set all correction variables to be applied to data

     To be used after the database is loaded but before loop over the root file
     The fibre position and direction have to be known (see FibreHandling class).
     Calculates the delay due to the trigger time and fills PMT ID histogram for
     multiple hits correction before analysis loop is started.

     Args:
       reader (RAT DSReader) : Iterator over the RAT::DS objects stored in a file passed 
                               to the reader

     Attributes:
       variables shared by all instances:
        pmt_prop (RAT PMT properties) : data base of the protperties of the SNO+ PMTs
        LightPath (RAT LightPath class) : calculates light path distances between points
                                          within the SNO+ detector
        groupVelTime (RAT GroupVelocity class) : calculates the total transit time of 
                                                 particles using group velocities

       variables unique to each instance:                                          
        r (RAT DSReader) : data structure reader
        beam_time_mean (float) : time the direct beam light arrives in PMTs 
        pmtid (TH1D) : histogram of hits vs the PMT ID the hit was registered in 
    """

    pmt_prop = rat.utility().GetPMTInfo()
    LightPath = rat.utility().GetLightPathCalculator()
    groupVelTime = rat.utility().GetGroupVelocity()

    def __init__(self, reader):
        self.r = reader
        self.beam_time_mean = 0.0
        self.pmtid = ROOT.TH1D("pmtid","pmtid",10000,0.0,10000.0)

    def fit_beam_time(self, sourcepos, sourcedir):
        """ Function to calculate the time the direct light is detected

         Args:
           sourcepos (TVector3) : vector of the position of the fibre in the detector
           sourcedir (TVector3) : direction vector the fibre is pointed at

         Returns:
           beam_time_mean (float) : time the direct light arrives in PMT

         Loops over all PMT hits and gets the time the PMT was hit and calculates the
         angle of the PMT with respect to the fibre direction. For all PMT hits within
         a certain angle (considered as direct beam light) a time residual histogram is 
         filled (time residual = PMT hit time - tansit time). The mean of this distribution 
         is returned as self.beam_time_mean
        """

        tres = ROOT.TH1D("tres","",500,0,500)
        #loop through file
        for ievent in range(0,self.r.GetEntryCount()):
            ds = self.r.GetEntry(ievent)
            for iev in range(ds.GetEVCount()):
                ev = ds.GetEV(iev)

                #get all PMTs in event
                pmts = ev.GetCalPMTs()
                for ipmt in range(0,pmts.GetCount()):
                    pmt_cal = pmts.GetPMT(ipmt)
                    pmt_id = pmt_cal.GetID()
                    pmt_time = pmt_cal.GetTime()
                    pmtpos = self.pmt_prop.GetPosition(pmt_id)

                    #only consider PMTs with valid timing
                    if pmt_time > 0:

                        #calculate transit time of photon to PMT
                        self.LightPath.CalcByPosition(sourcepos,pmtpos) 
                        PathTime = self.groupVelTime.CalcByDistance(self.LightPath.GetDistInScint(),self.LightPath.GetDistInAV(),self.LightPath.GetDistInWater())
                    
                        #calculate angle of PMT with respect to fibre position
                        pmtdir = (pmtpos - sourcepos)
                        alpha_mc_rad = math.acos((sourcedir * pmtdir)/(sourcedir.Mag() * pmtdir.Mag()))
                
                        if alpha_mc_rad<=(13.5/180.)*math.pi:                            
                            tres.Fill(pmt_time - PathTime)

        c1 = ROOT.TCanvas("c1","c1",1)
        c1.SetLogy(1)
        tres.Draw()
        max_bin = tres.GetBinCenter(tres.GetMaximumBin())
        f1 = ROOT.TF1("f1","gaus",max_bin-20,max_bin+20)
        tres.Fit(f1,"R")
        self.beam_time_mean = f1.GetParameter(1)
        c1.Close()
        return self.beam_time_mean

    def mh_corr(self):
        """ Function to fill PMT ID histogram to carry out multiple hits correction
        
         Returns:
           pmtid (TH1D) : histogram of hits vs the PMT ID the hit was registered in 
        """
        #loop through file
        for ievent in range(0,self.r.GetEntryCount()):
            ds = self.r.GetEntry(ievent)

            #loop through events
            for iev in range(ds.GetEVCount()):
                ev = ds.GetEV(iev)

                #run over pmts
                pmts = ev.GetCalPMTs()
                for ipmt in range(pmts.GetCount()):
                    pmt_cal = pmts.GetPMT(ipmt)
                    pmt_id = pmt_cal.GetID()
                    self.pmtid.Fill(pmt_id)
        return self.pmtid    
