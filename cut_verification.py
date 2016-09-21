#!/usr/bin/env python
import ROOT
import rat
import fibre_handling
import define_histograms
import math
import optparse
import os
import fnmatch

def find_track(mc, track_id):
    """Function to find track with a given track_id

     Args:
      mc (RAT::DS::MC) : RAT mc instance
      track_id (int) : track ID of investigated photon
    """
    for track_num in range(0,mc.GetMCTrackCount()):
        if mc.GetMCTrack(track_num).GetTrackID() == track_id:
            return mc.GetMCTrack(track_num)
    return None

def find_scattered(file_name, fibre, wl, ratio):
    """Function which runs over the MC branch of a RAT root file and applies
       the SMELLIE cut selection to all photons which have a Rayleigh scattering
       flag in their track history. Writes the number of photons selected by 
       each cut to an output file and saves histograms for each cut region to 
       a root file. Can only be used on MC simulations.
     
     Args:
      file_name (string) : name of the input rat file
      fibre (string) : fibre label
      wl (string) : wavelength
      ratio (string) : string component to define which scattering length 
                       scaling factor the root file was produced with
    """

    reader = ROOT.RAT.DU.DSReader(file_name,True)   

    #get fibre specific variables
    val = fibre_handling.FibreHandling(fibre)
    val.cut_values()

    sourcepos, sourcedir = val.get_fibre_position()
    AV1_cross, AV2_cross, PSUP_cross, n_scint, n_water = val.get_crossing_points(float(wl)) 

    #path lengths for direct beam
    scint_path = (AV2_cross - AV1_cross).Mag()
    water_path = (AV1_cross - sourcepos).Mag() + (PSUP_cross - AV2_cross).Mag()

    #get cut values
    maxBeam, z_beam_min, z_beam_max, alpha_min, alpha_max, z_avout_min, z_avout_max, alpha_avin = val.spatialcuts[0], val.spatialcuts[1], val.spatialcuts[2], val.spatialcuts[3], val.spatialcuts[4], val.spatialcuts[5], val.spatialcuts[6], val.spatialcuts[7]

    tbeam, beam_tres, tAV1, t, tAV, tpsup, tmulti = val.timecuts[0], val.timecuts[1], val.timecuts[2], val.timecuts[3], val.timecuts[4], val.timecuts[5], val.timecuts[6]

    #define output root file
    outputroot = ROOT.TFile("/data/langrock/rat-5.0-SMELLIE_analysis/" + str(fibre) + "/root/" + str(wl) + "_" + ratio + "_tracks.root","recreate")

    #define output text file
    outputfile = open("/data/langrock/rat-5.0-SMELLIE_analysis/" + str(fibre) + "/" + str(wl) + "_" + ratio + "_tracks.txt","w")

    #define histograms
    hist = define_histograms.DefineHistograms()

    #speed of light
    c = 300

    #variables used to count photons in cut region
    beam = 0
    avin = 0
    avout = 0
    scatt = 0 
    psup = 0
    multi = 0
    total = 0
    double_refl = 0 

    #rayleigh scattering track flag
    flag = ROOT.RAT.DS.MCTrack.OpRayleigh

    pmt_prop = rat.utility().GetPMTInfo()  
    LightPath = rat.utility().GetLightPathCalculator()
    groupVelTime = rat.utility().GetGroupVelocity()

    #start looping through file
    for ievent in range(0,reader.GetEntryCount()):
        ds, run = reader.GetEntry(ievent), reader.GetRun()
        mc = ds.GetMC()
    
        #run over pmts
        for ipmt in range(mc.GetMCPMTCount()): 
            pmt_id = mc.GetMCPMT(ipmt).GetID()
            #get pmt position and direction with respect to fibre position
            pmtpos = pmt_prop.GetPosition(pmt_id)
            pmtdir = (pmtpos - sourcepos)

            #define spatial variables to cut on
            z = pmtpos.Z()
            theta = pmtpos.Theta()
            phi = pmtpos.Phi()
            alpha_mc_rad = math.acos((sourcedir * pmtdir)/(sourcedir.Mag() * pmtdir.Mag()))
            alpha_mc = math.degrees(alpha_mc_rad)        

            #calculate time it takes the photon in respective pmt to get there
            LightPath.CalcByPosition(sourcepos,pmtpos) 
            PathTime = groupVelTime.CalcByDistance(LightPath.GetDistInScint(),LightPath.GetDistInAV(),LightPath.GetDistInWater())

            #time for direct light to cross detector
            Beam_time = (scint_path*n_scint + water_path*n_water)/c
	    #AV1 reflection time off the outside of the AV
            AV_ref1_time = ((pmtpos - AV1_cross).Mag() + (AV1_cross - sourcepos).Mag()) * n_water /c        
            #AV2 reflection time off the inside of the AV after crossing the detector
            AV_ref2_time = (((pmtpos - AV2_cross).Mag() + (AV2_cross - sourcepos).Mag() - water_path)*n_scint + water_path*n_water) /c 
            #PSUP reflection time
            PSUP_ref_time = (((pmtpos - PSUP_cross).Mag() + scint_path - water_path)*n_scint + 2*water_path*n_water) /c

            #loop through photons in PMT
            mc_pmt = mc.GetMCPMT(ipmt)
            for photon in range(mc_pmt.GetMCPECount()):
                mc_photon = mc_pmt.GetMCPE(photon)
                pmt_time = mc_photon.GetCreationTime()
                time = pmt_time - PathTime

                #find photon track for selected photon (works only for non-noise photons)
                if not mc_photon.GetNoise():
                    track_id = mc_photon.GetPhotonTrackID()
                    track = find_track(mc,track_id)
                    
                    #if track contains Rayleigh scattering flag, apply cuts, count photons and fill histograms for each each cut
                    if track.GetSummaryFlag(flag):
            
                        #coutn total number of photons detected and fill histograms
                        total += 1
                        hist.t_res.Fill(time)
                        hist.angle_time.Fill(time,alpha_mc)
                        hist.z_time.Fill(time,z)
                        hist.theta_phi.Fill(phi,theta)
                        hist.h_theta.Fill(theta)
                        hist.h_phi.Fill(phi)

                        #apply direct beam cuts
                        if alpha_mc_rad<=(maxBeam/180.)*math.pi and z < z_beam_max and z > z_beam_min and time < Beam_time+tbeam and (pmt_time - PathTime) < beam_tres:
                            beam += 1
                    
                            hist.t_res_beam.Fill(time)
                            hist.angle_time_beam.Fill(time,alpha_mc)
                            hist.z_time_beam.Fill(time,z)    
                            hist.theta_phi_beam.Fill(phi,theta)
                            hist.h_theta_beam.Fill(theta)
                            hist.h_phi_beam.Fill(phi)

                        #apply late pulse cuts
                        elif alpha_mc_rad<=(maxBeam/180.)*math.pi and z < z_beam_max and z > z_beam_min and time < Beam_time+tbeam and (pmt_time - PathTime) > beam_tres and (pmt_time - PathTime) < 50:
                            double_refl += 1
                    
                            hist.t_res_double.Fill(time)
                            hist.angle_time_double.Fill(time,alpha_mc)
                            hist.z_time_double.Fill(time,z)    
                            hist.theta_phi_double.Fill(phi,theta)
                            hist.h_theta_double.Fill(theta)
                            hist.h_phi_double.Fill(phi)


                        else:
                            #apply cuts on outer (1st) AV reflections
                            if time < AV_ref1_time+tAV1 and alpha_mc_rad > (alpha_min/180.)*math.pi and alpha_mc_rad < (alpha_max/180.)*math.pi and (pmt_time -  PathTime) < t and z < z_avout_max and z > z_avout_min: 
                                avout += 1

                                hist.t_res_avout.Fill(time)
                                hist.angle_time_avout.Fill(time,alpha_mc)
                                hist.z_time_avout.Fill(time,z)
                                hist.theta_phi_avout.Fill(phi,theta)
                                hist.h_theta_avout.Fill(theta)
                                hist.h_phi_avout.Fill(phi)
                            
                            #apply cuts on scattered events
                            elif time < AV_ref2_time-tAV:
                                scatt += 1
                        
                                hist.t_res_scatt.Fill(time)
                                hist.angle_time_scatt.Fill(time,alpha_mc)
                                hist.z_time_scatt.Fill(time,z)
                                hist.theta_phi_scatt.Fill(phi,theta)
                                hist.h_theta_scatt.Fill(theta)
                                hist.h_phi_scatt.Fill(phi)

                            #apply cuts on inner (2nd) AV reflections
                            elif  time > AV_ref2_time-tAV and ((time < PSUP_ref_time-tpsup and alpha_mc_rad > (alpha_avin/180.)*math.pi and alpha_mc_rad < ((alpha_avin+15)/180.)*math.pi) or (time < PSUP_ref_time-tpsup+10 and alpha_mc_rad > ((alpha_avin+15)/180.)*math.pi and alpha_mc_rad < ((alpha_avin+20)/180.)*math.pi) or (time < PSUP_ref_time-tpsup+20 and alpha_mc_rad > ((alpha_avin+20)/180.)*math.pi and alpha_mc_rad < ((alpha_avin+30)/180.)*math.pi) or (time < PSUP_ref_time-tpsup+25 and alpha_mc_rad > ((alpha_avin+30)/180.)*math.pi and alpha_mc_rad < ((alpha_avin+40)/180.)*math.pi) or (time < PSUP_ref_time-tpsup+35 and alpha_mc_rad > ((alpha_avin+40)/180.)*math.pi and alpha_mc_rad < ((alpha_avin+50)/180.)*math.pi) or (time < PSUP_ref_time-tpsup+40 and alpha_mc_rad > ((alpha_avin+50)/180.)*math.pi and alpha_mc_rad < ((alpha_avin+60)/180.)*math.pi) or (time < PSUP_ref_time-tpsup+45 and alpha_mc_rad > ((alpha_avin+60)/180.)*math.pi)):
                                avin += 1
                        
                                hist.t_res_avin.Fill(time)
                                hist.angle_time_avin.Fill(time,alpha_mc)
                                hist.z_time_avin.Fill(time,z)
                                hist.theta_phi_avin.Fill(phi,theta)
                                hist.h_theta_avin.Fill(theta)
                                hist.h_phi_avin.Fill(phi)

                            #apply cuts on PSUP reflections
                            elif time >  AV_ref2_time-tAV and time < PSUP_ref_time+tmulti: 
                                psup += 1

                                hist.t_res_psup.Fill(time)
                                hist.angle_time_psup.Fill(time,alpha_mc)
                                hist.z_time_psup.Fill(time,z)
                                hist.theta_phi_psup.Fill(phi,theta)
                                hist.h_theta_psup.Fill(theta)
                                hist.h_phi_psup.Fill(phi)

                            #apply cuts on multiple effects
                            elif time > PSUP_ref_time+tmulti:
                                multi += 1

                                hist.t_res_multi.Fill(time)
                                hist.angle_time_multi.Fill(time,alpha_mc)
                                hist.z_time_multi.Fill(time,z)
                                hist.theta_phi_multi.Fill(phi,theta)
                                hist.h_theta_multi.Fill(theta)
                                hist.h_phi_multi.Fill(phi)

    #save histograms to root file
    outputroot.Write()
    outputroot.Close()

    #save all values to a txt file
    outputfile.write("total: " + str(total) + "\n")
    outputfile.write("beam: " + str(beam) + "\n")
    outputfile.write("double_refl: " + str(double_refl) + "\n")
    outputfile.write("avin: " + str(avin) + "\n")
    outputfile.write("avout: " + str(avout) + "\n")
    outputfile.write("scatt: " + str(scatt) + "\n")
    outputfile.write("psup: " + str(psup) + "\n")
    outputfile.write("multi: " + str(multi) + "\n")

    outputfile.close()

def find_noise(file_name, fibre, wl, ratio):
    """Function which runs over the MC branch of a RAT root file and applies
       the SMELLIE cut selection to all photons which originate from PMT noise. 
       Writes the number of photons selected by each cut to an output file and 
       saves histograms for each cut region to a root file. Can only be used on 
       MC simulations.
     
     Args:
      file_name (string) : name of the input rat file
      fibre (string) : fibre label
      wl (string) : wavelength
      ratio (string) : string component to define which scattering length 
                       scaling factor the root file was produced with
    """

    reader = ROOT.RAT.DU.DSReader(file_name,True)   

    #get fibre specific variables
    val = fibre_handling.FibreHandling(fibre)
    val.cut_values()

    sourcepos, sourcedir = val.get_fibre_position()
    AV1_cross, AV2_cross, PSUP_cross, n_scint, n_water = val.get_crossing_points(float(wl)) 

    #path lengths for direct beam
    scint_path = (AV2_cross - AV1_cross).Mag()
    water_path = (AV1_cross - sourcepos).Mag() + (PSUP_cross - AV2_cross).Mag()

    #get cut values
    maxBeam, z_beam_min, z_beam_max, alpha_min, alpha_max, z_avout_min, z_avout_max, alpha_avin = val.spatialcuts[0], val.spatialcuts[1], val.spatialcuts[2], val.spatialcuts[3], val.spatialcuts[4], val.spatialcuts[5], val.spatialcuts[6], val.spatialcuts[7]

    tbeam, beam_tres, tAV1, t, tAV, tpsup, tmulti = val.timecuts[0], val.timecuts[1], val.timecuts[2], val.timecuts[3], val.timecuts[4], val.timecuts[5], val.timecuts[6]

    #define output root file
    outputroot = ROOT.TFile("/data/langrock/rat-5.0-SMELLIE_analysis/" + str(fibre) + "/root/" + str(wl) + "_" + ratio + "_noise.root","recreate")

    #define output text file
    outputfile = open("/data/langrock/rat-5.0-SMELLIE_analysis/" + str(fibre) + "/" + str(wl) + "_" + ratio + "_noise.txt","w")

    #define histograms
    hist = define_histograms.DefineHistograms()

    #speed of light
    c = 300

    #variables used to count photons in cut region
    beam = 0
    avin = 0
    avout = 0
    scatt = 0 
    psup = 0
    multi = 0
    total = 0
    double_refl = 0

    pmt_prop = rat.utility().GetPMTInfo()  
    LightPath = rat.utility().GetLightPathCalculator()
    groupVelTime = rat.utility().GetGroupVelocity()

    #start looping through file
    for ievent in range(0,reader.GetEntryCount()):
        ds, run = reader.GetEntry(ievent), reader.GetRun()
        mc = ds.GetMC()
    
        #run over pmts
        for ipmt in range(mc.GetMCPMTCount()): 
            pmt_id = mc.GetMCPMT(ipmt).GetID()
            #get pmt position and direction with respect to fibre position
            pmtpos = pmt_prop.GetPosition(pmt_id)
            pmtdir = (pmtpos - sourcepos)

            #define spatial variables to cut on
            z = pmtpos.Z()
            theta = pmtpos.Theta()
            phi = pmtpos.Phi()
            alpha_mc_rad = math.acos((sourcedir * pmtdir)/(sourcedir.Mag() * pmtdir.Mag()))
            alpha_mc = math.degrees(alpha_mc_rad)        

            #calculate time it takes the photon in respective pmt to get there
            LightPath.CalcByPosition(sourcepos,pmtpos) 
            PathTime = groupVelTime.CalcByDistance(LightPath.GetDistInScint(),LightPath.GetDistInAV(),LightPath.GetDistInWater())

            #time for direct light to cross detector
            Beam_time = (scint_path*n_scint + water_path*n_water)/c
	    #AV1 reflection time off the outside of the AV
            AV_ref1_time = ((pmtpos - AV1_cross).Mag() + (AV1_cross - sourcepos).Mag()) * n_water /c        
            #AV2 reflection time off the inside of the AV after crossing the detector
            AV_ref2_time = (((pmtpos - AV2_cross).Mag() + (AV2_cross - sourcepos).Mag() - water_path)*n_scint + water_path*n_water) /c 
            #PSUP reflection time
            PSUP_ref_time = (((pmtpos - PSUP_cross).Mag() + scint_path - water_path)*n_scint + 2*water_path*n_water) /c

            #loop through photons in PMT
            mc_pmt = mc.GetMCPMT(ipmt)
            for photon in range(mc_pmt.GetMCPECount()):
                mc_photon = mc_pmt.GetMCPE(photon)
                pmt_time = mc_photon.GetCreationTime()
                time = pmt_time - PathTime

                #if photon is a noise hit, apply cuts, count photons and fill histograms for each each cut
                if mc_photon.GetNoise():

                    #count total number of photons detected and fill histograms
                    total += 1                    
                    hist.t_res.Fill(time)
                    hist.angle_time.Fill(time,alpha_mc)
                    hist.z_time.Fill(time,z)
                    hist.theta_phi.Fill(phi,theta)
                    hist.h_theta.Fill(theta)
                    hist.h_phi.Fill(phi)

                    #apply direct beam cuts
                    if alpha_mc_rad<=(maxBeam/180.)*math.pi and z < z_beam_max and z > z_beam_min and time < Beam_time+tbeam and (pmt_time - PathTime) < beam_tres:
                        beam += 1
                    
                        hist.t_res_beam.Fill(time)
                        hist.angle_time_beam.Fill(time,alpha_mc)
                        hist.z_time_beam.Fill(time,z)    
                        hist.theta_phi_beam.Fill(phi,theta)
                        hist.h_theta_beam.Fill(theta)
                        hist.h_phi_beam.Fill(phi)

                    #apply late pulse cuts
                    elif alpha_mc_rad<=(maxBeam/180.)*math.pi and z < z_beam_max and z > z_beam_min and time < Beam_time+tbeam and (pmt_time - PathTime) > beam_tres and (pmt_time - PathTime) < 50:
                        double_refl += 1
                    
                        hist.t_res_double.Fill(time)
                        hist.angle_time_double.Fill(time,alpha_mc)
                        hist.z_time_double.Fill(time,z)    
                        hist.theta_phi_double.Fill(phi,theta)
                        hist.h_theta_double.Fill(theta)
                        hist.h_phi_double.Fill(phi)

                    else:
                        #apply cuts on outer (1st) AV reflections
                        if time < AV_ref1_time+tAV1 and alpha_mc_rad > (alpha_min/180.)*math.pi and alpha_mc_rad < (alpha_max/180.)*math.pi and (pmt_time -  PathTime) < t and z < z_avout_max and z > z_avout_min: 
                            avout += 1

                            hist.t_res_avout.Fill(time)
                            hist.angle_time_avout.Fill(time,alpha_mc)
                            hist.z_time_avout.Fill(time,z)
                            hist.theta_phi_avout.Fill(phi,theta)
                            hist.h_theta_avout.Fill(theta)
                            hist.h_phi_avout.Fill(phi)

                        #apply cuts on scattered events
                        elif time < AV_ref2_time-tAV:
                            scatt += 1
                            
                            hist.t_res_scatt.Fill(time)
                            hist.angle_time_scatt.Fill(time,alpha_mc)
                            hist.z_time_scatt.Fill(time,z)
                            hist.theta_phi_scatt.Fill(phi,theta)
                            hist.h_theta_scatt.Fill(theta)
                            hist.h_phi_scatt.Fill(phi)
                            
                        #apply cuts on inner (2nd) AV reflections
                        elif time > AV_ref2_time-tAV and ((time < PSUP_ref_time-tpsup and alpha_mc_rad > (alpha_avin/180.)*math.pi and alpha_mc_rad < ((alpha_avin+15)/180.)*math.pi) or (time < PSUP_ref_time-tpsup+10 and alpha_mc_rad > ((alpha_avin+15)/180.)*math.pi and alpha_mc_rad < ((alpha_avin+20)/180.)*math.pi) or (time < PSUP_ref_time-tpsup+20 and alpha_mc_rad > ((alpha_avin+20)/180.)*math.pi and alpha_mc_rad < ((alpha_avin+30)/180.)*math.pi) or (time < PSUP_ref_time-tpsup+25 and alpha_mc_rad > ((alpha_avin+30)/180.)*math.pi and alpha_mc_rad < ((alpha_avin+40)/180.)*math.pi) or (time < PSUP_ref_time-tpsup+35 and alpha_mc_rad > ((alpha_avin+40)/180.)*math.pi and alpha_mc_rad < ((alpha_avin+50)/180.)*math.pi) or (time < PSUP_ref_time-tpsup+40 and alpha_mc_rad > ((alpha_avin+50)/180.)*math.pi and alpha_mc_rad < ((alpha_avin+60)/180.)*math.pi) or (time < PSUP_ref_time-tpsup+45 and alpha_mc_rad > ((alpha_avin+60)/180.)*math.pi)):
                            avin += 1
                        
                            hist.t_res_avin.Fill(time)
                            hist.angle_time_avin.Fill(time,alpha_mc)
                            hist.z_time_avin.Fill(time,z)
                            hist.theta_phi_avin.Fill(phi,theta)
                            hist.h_theta_avin.Fill(theta)
                            hist.h_phi_avin.Fill(phi)

                        #apply cuts on PSUP reflections
                        elif time >  AV_ref2_time-tAV and time < PSUP_ref_time+tmulti: 
                            psup += 1

                            hist.t_res_psup.Fill(time)
                            hist.angle_time_psup.Fill(time,alpha_mc)
                            hist.z_time_psup.Fill(time,z)
                            hist.theta_phi_psup.Fill(phi,theta)
                            hist.h_theta_psup.Fill(theta)
                            hist.h_phi_psup.Fill(phi)

                        #apply cuts on multiple effects
                        elif time > PSUP_ref_time+tmulti:
                            multi += 1

                            hist.t_res_multi.Fill(time)
                            hist.angle_time_multi.Fill(time,alpha_mc)
                            hist.z_time_multi.Fill(time,z)
                            hist.theta_phi_multi.Fill(phi,theta)
                            hist.h_theta_multi.Fill(theta)
                            hist.h_phi_multi.Fill(phi)

    #save histograms to root file
    outputroot.Write()
    outputroot.Close()

    #save all values to a text file
    outputfile.write("total: " + str(total) + "\n")
    outputfile.write("beam: " + str(beam) + "\n")
    outputfile.write("double_refl: " + str(double_refl) + "\n")
    outputfile.write("avin: " + str(avin) + "\n")
    outputfile.write("avout: " + str(avout) + "\n")
    outputfile.write("scatt: " + str(scatt) + "\n")
    outputfile.write("psup: " + str(psup) + "\n")
    outputfile.write("multi: " + str(multi) + "\n")

    outputfile.close()

if __name__ == '__main__':
    """Script runs through all files in the directory given by the -d flag. It only selects the 
       files for the given fibre-wavelength combination (-f and -w flags). It runs either over 
       the MC branch and checks for either scattered photons or noise hits. The root file names 
       have to have the structure fibre_wavelength_ratio.root
    """

    parser = optparse.OptionParser( usage = "python %prog [flags]")
    parser.add_option("-d", dest="dir_name", help="Name of directory for input files")
    parser.add_option("-f", dest="fibre", help="Fibre to be measured")
    parser.add_option("-w", dest="wavelength", help="Wavelength of fibre")
    parser.add_option("-t", dest="type", help="Check type, needs to be 'noise' or 'tracks'")
    
    (options, args) = parser.parse_args()

    if not options.dir_name:
        print 'Need directory name!'
        raise Exception
    if not options.type:
        print 'Need check type!'
        raise Exception
    if not options.fibre:
        print 'Need fibre!'
        raise Exception
    if not options.wavelength:
        print 'Need wavelength!'
        raise Exception

    #loop through given directory
    for root, dirs, filenames in os.walk(options.dir_name):
        for x in filenames:
            if fnmatch.fnmatch(x, "*.root"):
                infile = options.dir_name + x
                    
                #split the file name to get fibre, wavelength and ratio
                components = x.split("_")
                if len(components)!=0:
                    ratio = ""

                    for i in range(len(components)):
                        if "." in components[i]:
                            wl = components[i].split(".")
                            components[i] = wl[0]

                    #rum through all files for the chosen fibre-wavelength combination
                    if str(options.fibre) == str(components[0]) and str(options.wavelength) == str(components[1]):
                        if len(components) > 2:
                            ratio = components[2] 

                        #count noise hits
                        if options.type == "noise":
                            find_noise(infile,options.fibre,options.wavelength, ratio)

                        #count scattered photons
                        elif options.type == "tracks":
                            find_scattered(infile,options.fibre,options.wavelength, ratio)
                        else:
                            print "Invalid check type!"
                            raise Exception

                    else:
                        continue

