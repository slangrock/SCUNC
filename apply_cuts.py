#!/usr/bin/env python
import ROOT
import rat
import fibre_handling
import data_corrections
import define_histograms
import math
import optparse
import os
import fnmatch

def generate_mc_output(file_name, fibre, wl, ratio):
    """Function which runs over the MC branch of a RAT root file and applies
       the SMELLIE cut selection. Writes the number of photons selected by 
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
    outputroot = ROOT.TFile("/data/langrock/rat-5.0-SMELLIE_analysis/" + str(fibre) + "/root/" + str(wl) +  "_" + str(ratio) + "_mc.root","recreate")

    #define output text file
    outputfile = open("/data/langrock/rat-5.0-SMELLIE_analysis/" + str(fibre) + "/" + str(wl) + "_" + ratio + ".txt","w")

    #define histograms
    hist = define_histograms.DefineHistograms()
    
    #speed of light
    c = 300

    #variables used to count photons in cut region
    beam = 0
    double_refl = 0
    avin = 0
    avout = 0
    scatt = 0 
    psup = 0
    multi = 0
    total = 0


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

                #count total number of photons detected and fill histograms
                total += 1
                hist.t_res.Fill(time)
                hist.angle_time.Fill(time,alpha_mc)
                hist.z_time.Fill(time,z)
                hist.theta_phi.Fill(phi,theta)
                hist.h_theta.Fill(theta)
                hist.h_phi.Fill(phi)

                ##### apply cuts, count photons and fill histograms for each each cut ####

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

    #calculate poissonian uncertainty on all cut region counts
    total_error = math.sqrt(float(total))
    double_refl_error = math.sqrt(float(double_refl))
    scatt_error = math.sqrt(float(scatt))
    beam_error = math.sqrt(float(beam))
    avin_error = math.sqrt(float(avin))
    avout_error = math.sqrt(float(avout))
    psup_error = math.sqrt(float(psup))
    multi_error = math.sqrt(float(multi))

    #calculate the ratio of scattered/in-beam events to the total number of events and the uncertainty on the ratio
    ratio_scatt = float(scatt)/float(total)
    ratio_scatt_error = math.sqrt(math.pow(scatt_error/float(total),2) + math.pow(float(scatt)*total_error/math.pow(float(total),2),2))
    ratio_beam = float(beam)/float(total)
    ratio_beam_error = math.sqrt(math.pow(beam_error/float(total),2) + math.pow(float(beam)*total_error/math.pow(float(total),2),2))

    #save all values to a txt file
    outputfile.write(wl)
    outputfile.write('\t')
    outputfile.write(ratio.replace('p','.'))
    outputfile.write('\t scattering ratio: ')
    outputfile.write(str(ratio_scatt))
    outputfile.write('\t +/- ')
    outputfile.write(str(ratio_scatt_error))
    outputfile.write('\t beam ratio: ')
    outputfile.write(str(ratio_beam))
    outputfile.write('\t +/- ')
    outputfile.write(str(ratio_beam_error))
    outputfile.write('\n \n \n \n')

    outputfile.write("total: " + str(total) + " +/-" + str(total_error) + "\n")
    outputfile.write("beam: " + str(beam) + " +/-" + str(beam_error) + "\n")
    outputfile.write("double_refl: " + str(double_refl) + " +/-" + str(double_refl_error) + "\n")
    outputfile.write("avin: " + str(avin) + " +/-" + str(avin_error) + "\n")
    outputfile.write("avout: " + str(avout) + " +/-" + str(avout_error) + "\n")
    outputfile.write("scatt: " + str(scatt) + " +/-" + str(scatt_error) + "\n")
    outputfile.write("psup: " + str(psup) + " +/-" + str(psup_error) + "\n")
    outputfile.write("multi: " + str(multi) + " +/-" + str(multi_error) + "\n")

    outputfile.close()


def generate_data_output(file_name, fibre, wl, ratio):
    """Function which runs over the EV branch of a RAT root file and applies
       the SMELLIE cut selection. Writes the number of hitd selection by 
       each cut to an output file and saves histograms for each cut region to 
       a root file. Can be used on MC simulations and SNO+ data.
     
     Args:
      file_name (string) : name of the input rat file
      fibre (string) : fibre label
      wl (string) : wavelength
      ratio (string) : string component to define which scattering length 
                       scaling factor the root file was produced with (in case of 
                       the simulation)
    """

    reader = ROOT.RAT.DU.DSReader(file_name,True)

    #get fibre specific variabes
    val = fibre_handling.FibreHandling(fibre)
    val.cut_values()

    sourcepos, sourcedir = val.get_fibre_position()
    AV1_cross, AV2_cross, PSUP_cross, n_scint, n_water = val.get_crossing_points(float(wl)) 

    #path lengths for direct beam
    scint_path = (AV2_cross - AV1_cross).Mag()
    water_path = (AV1_cross - sourcepos).Mag() + (PSUP_cross - AV2_cross).Mag()

    maxBeam, z_beam_min, z_beam_max, alpha_min, alpha_max, z_avout_min, z_avout_max, alpha_avin = val.spatialcuts[0], val.spatialcuts[1], val.spatialcuts[2], val.spatialcuts[3], val.spatialcuts[4], val.spatialcuts[5], val.spatialcuts[6], val.spatialcuts[7]
 
    tbeam, beam_tres, tAV1, t, tAV, tpsup, tmulti = val.timecuts[0], val.timecuts[1], val.timecuts[2], val.timecuts[3], val.timecuts[4], val.timecuts[5], val.timecuts[6]

    #define output root file
    outputroot = ROOT.TFile("/data/langrock/rat-5.0-SMELLIE_analysis/" + str(fibre) + "/root/" + str(wl) + "_" + str(ratio) + "_data_water_c_mh.root","recreate")

    #define output text file
    outputfile = open("/data/langrock/rat-5.0-SMELLIE_analysis/" + str(fibre) + "/" + str(wl) + "_" + ratio + "_water.txt","w")

    #get tools for data correction
    corrections = data_corrections.DataCorrections(reader)
    beam_center = corrections.fit_beam_time(sourcepos,sourcedir)
    pmtid = corrections.mh_corr()
    nevents = reader.GetEntryCount()

    #define histograms
    hist = define_histograms.DefineHistograms()

    #speed of light
    c = 300

    #variables used to count photons in cut region
    beam = 0
    double_refl = 0
    avin = 0
    avout = 0
    scatt = 0 
    psup = 0
    multi = 0
    total = 0


    pmt_prop = rat.utility().GetPMTInfo()  
    LightPath = rat.utility().GetLightPathCalculator()
    groupVelTime = rat.utility().GetGroupVelocity()

    #start looping through file
    for ievent in range(0,nevents):
        ds, run = reader.GetEntry(ievent), reader.GetRun()

        #loop through events
        for iev in range(ds.GetEVCount()):
            ev = ds.GetEV(iev)
            hist.h_nhits.Fill(ev.GetNhits())

            pmts = ev.GetCalPMTs()
            #run over pmts
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

                #define spatial variables to cut on
                z = pmtpos.Z()
                theta = pmtpos.Theta()
                phi = pmtpos.Phi()
                alpha_mc_rad = math.acos((sourcedir * pmtdir)/(sourcedir.Mag() * pmtdir.Mag()))
                alpha_mc = math.degrees(alpha_mc_rad)        

                pmt_time = pmt_cal.GetTime()

                #calculate time it takes the photon in respective pmt to get there
                LightPath.CalcByPosition(sourcepos,pmtpos) 
                PathTime = groupVelTime.CalcByDistance(LightPath.GetDistInScint(),LightPath.GetDistInAV(),LightPath.GetDistInWater())

                time = pmt_time - PathTime

                #time for direct light to cross detector
                Beam_time = beam_center + (scint_path*n_scint + water_path*n_water)/c
		#AV1 reflection time off the outside of the AV
                AV_ref1_time = beam_center + ((pmtpos - AV1_cross).Mag() + (AV1_cross - sourcepos).Mag()) * n_water /c        
                #AV2 reflection time off the inside of the AV after crossing the detector
                AV_ref2_time = beam_center + (((pmtpos - AV2_cross).Mag() + (AV2_cross - sourcepos).Mag() - water_path)*n_scint + water_path*n_water) /c 
                #PSUP reflection time
                PSUP_ref_time = beam_center + (((pmtpos - PSUP_cross).Mag() + scint_path - water_path)*n_scint + 2*water_path*n_water) /c

                #count total number of hits detected and fill histograms, apply multiple hits correction
                total += 1*c_mh
                hist.t_res.Fill(time-beam_center,c_mh)
                hist.angle_time.Fill(time-beam_center,alpha_mc)
                hist.z_time.Fill(time-beam_center,z)
                hist.theta_phi.Fill(phi,theta)
                hist.h_theta.Fill(theta,c_mh)
                hist.h_phi.Fill(phi,c_mh)

                ##### apply cuts, count photons and fill histograms for each each cut, including multiple hits correction ####

                #apply direct beam cuts
                if alpha_mc_rad<=(maxBeam/180.)*math.pi and z < z_beam_max and z > z_beam_min and time < Beam_time+tbeam and (pmt_time - PathTime - beam_center) <= beam_tres:
                    beam += 1*c_mh
                    
                    hist.t_res_beam.Fill(time-beam_center,c_mh)
                    hist.angle_time_beam.Fill(time-beam_center,alpha_mc)
                    hist.z_time_beam.Fill(time-beam_center,z)    
                    hist.theta_phi_beam.Fill(phi,theta)
                    hist.h_theta_beam.Fill(theta,c_mh)
                    hist.h_phi_beam.Fill(phi,c_mh)

                #apply late pulse cuts
                elif alpha_mc_rad<=(maxBeam/180.)*math.pi and z < z_beam_max and z > z_beam_min and time < Beam_time+tbeam and (pmt_time - PathTime - beam_center) > beam_tres and (pmt_time - PathTime - beam_center) < 50:
                    double_refl += 1*c_mh
                    
                    hist.t_res_double.Fill(time-beam_center,c_mh)
                    hist.angle_time_double.Fill(time-beam_center,alpha_mc)
                    hist.z_time_double.Fill(time-beam_center,z)    
                    hist.theta_phi_double.Fill(phi,theta)
                    hist.h_theta_double.Fill(theta,c_mh)
                    hist.h_phi_double.Fill(phi,c_mh)

                else:
                    #apply cuts on outer (1st) AV reflections
                    if time < AV_ref1_time+tAV1 and alpha_mc_rad > (alpha_min/180.)*math.pi and alpha_mc_rad < (alpha_max/180.)*math.pi and (pmt_time -  PathTime - beam_center) < t and z < z_avout_max and z > z_avout_min: 
                        avout += 1*c_mh

                        hist.t_res_avout.Fill(time-beam_center,c_mh)
                        hist.angle_time_avout.Fill(time-beam_center,alpha_mc)
                        hist.z_time_avout.Fill(time-beam_center,z)
                        hist.theta_phi_avout.Fill(phi,theta)
                        hist.h_theta_avout.Fill(theta,c_mh)
                        hist.h_phi_avout.Fill(phi,c_mh)

                    #apply cuts on scattered events
                    elif time < AV_ref2_time-tAV:
                        scatt += 1*c_mh
                        
                        hist.t_res_scatt.Fill(time-beam_center,c_mh)
                        hist.angle_time_scatt.Fill(time-beam_center,alpha_mc)
                        hist.z_time_scatt.Fill(time-beam_center,z)
                        hist.theta_phi_scatt.Fill(phi,theta)
                        hist.h_theta_scatt.Fill(theta,c_mh)
                        hist.h_phi_scatt.Fill(phi,c_mh)

                    #apply cuts on inner (2nd) AV reflections
                    elif time > AV_ref2_time-tAV and ((time < PSUP_ref_time-tpsup and alpha_mc_rad > (alpha_avin/180.)*math.pi and alpha_mc_rad < ((alpha_avin+15)/180.)*math.pi) or (time < PSUP_ref_time-tpsup+10 and alpha_mc_rad > ((alpha_avin+15)/180.)*math.pi and alpha_mc_rad < ((alpha_avin+20)/180.)*math.pi) or (time < PSUP_ref_time-tpsup+20 and alpha_mc_rad > ((alpha_avin+20)/180.)*math.pi and alpha_mc_rad < ((alpha_avin+30)/180.)*math.pi) or (time < PSUP_ref_time-tpsup+25 and alpha_mc_rad > ((alpha_avin+30)/180.)*math.pi and alpha_mc_rad < ((alpha_avin+40)/180.)*math.pi) or (time < PSUP_ref_time-tpsup+35 and alpha_mc_rad > ((alpha_avin+40)/180.)*math.pi and alpha_mc_rad < ((alpha_avin+50)/180.)*math.pi) or (time < PSUP_ref_time-tpsup+40 and alpha_mc_rad > ((alpha_avin+50)/180.)*math.pi and alpha_mc_rad < ((alpha_avin+60)/180.)*math.pi) or (time < PSUP_ref_time-tpsup+45 and alpha_mc_rad > ((alpha_avin+60)/180.)*math.pi)):
                        avin += 1*c_mh
                        
                        hist.t_res_avin.Fill(time-beam_center,c_mh)
                        hist.angle_time_avin.Fill(time-beam_center,alpha_mc)
                        hist.z_time_avin.Fill(time-beam_center,z)
                        hist.theta_phi_avin.Fill(phi,theta)
                        hist.h_theta_avin.Fill(theta,c_mh)
                        hist.h_phi_avin.Fill(phi,c_mh)

                    #apply cuts on PSUP reflections
                    elif time >  AV_ref2_time-tAV and time < PSUP_ref_time+tmulti: 
                        psup += 1*c_mh

                        hist.t_res_psup.Fill(time-beam_center,c_mh)
                        hist.angle_time_psup.Fill(time-beam_center,alpha_mc)
                        hist.z_time_psup.Fill(time-beam_center,z)
                        hist.theta_phi_psup.Fill(phi,theta)
                        hist.h_theta_psup.Fill(theta,c_mh)
                        hist.h_phi_psup.Fill(phi,c_mh)

                    #apply cuts on multiple effects
                    elif time > PSUP_ref_time+tmulti:
                        multi += 1*c_mh

                        hist.t_res_multi.Fill(time-beam_center,c_mh)
                        hist.angle_time_multi.Fill(time-beam_center,alpha_mc)
                        hist.z_time_multi.Fill(time-beam_center,z)
                        hist.theta_phi_multi.Fill(phi,theta)
                        hist.h_theta_multi.Fill(theta,c_mh)
                        hist.h_phi_multi.Fill(phi,c_mh)

    #save histograms to root file
    outputroot.Write()
    outputroot.Close()

    #calculate poissonian uncertainty on all cut region counts
    total_error = math.sqrt(float(total))
    scatt_error = math.sqrt(float(scatt))
    beam_error = math.sqrt(float(beam))
    double_refl_error = math.sqrt(float(double_refl))
    avin_error = math.sqrt(float(avin))
    avout_error = math.sqrt(float(avout))
    psup_error = math.sqrt(float(psup))
    multi_error = math.sqrt(float(multi))

    #calculate the ratio of scattered/in-beam events to the total number of events and the uncertainty on the ratio
    ratio_scatt = float(scatt)/float(total)
    ratio_scatt_error = math.sqrt(math.pow(scatt_error/float(total),2) + math.pow(float(scatt)*total_error/math.pow(float(total),2),2))
    ratio_beam = float(beam)/float(total)
    ratio_beam_error = math.sqrt(math.pow(beam_error/float(total),2) + math.pow(float(beam)*total_error/math.pow(float(total),2),2))

    #save all values to a txt file
    outputfile.write(wl)
    outputfile.write('\t')
    outputfile.write(ratio.replace('p','.'))
    outputfile.write('\t scattering ratio: ')
    outputfile.write(str(ratio_scatt))
    outputfile.write('\t +/- ')
    outputfile.write(str(ratio_scatt_error))
    outputfile.write('\t beam ratio: ')
    outputfile.write(str(ratio_beam))
    outputfile.write('\t +/- ')
    outputfile.write(str(ratio_beam_error))
    outputfile.write('\n \n \n \n')

    outputfile.write("total: " + str(total) + " +/-" + str(total_error) + "\n")
    outputfile.write("beam: " + str(beam) + " +/-" + str(beam_error) + "\n")
    outputfile.write("double_refl: " + str(double_refl) + " +/-" + str(double_refl_error) + "\n")
    outputfile.write("avin: " + str(avin) + " +/-" + str(avin_error) + "\n")
    outputfile.write("avout: " + str(avout) + " +/-" + str(avout_error) + "\n")
    outputfile.write("scatt: " + str(scatt) + " +/-" + str(scatt_error) + "\n")
    outputfile.write("psup: " + str(psup) + " +/-" + str(psup_error) + "\n")
    outputfile.write("multi: " + str(multi) + " +/-" + str(multi_error) + "\n")

    outputfile.close()

if __name__ == '__main__':
    """Script runs through all files in the directory given by the -d flag. It only selects the 
       files for the given fibre-wavelength combination (-f and -w flags). It runs either over 
       the MC branch (simultions only!) or the EV branch of the RAT root files. The root file
       names have to have the structure fibre_wavelength_ratio.root
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

                    #run through all files for the chosen fibre-wavelength combination
                    if str(options.fibre) == str(components[0]) and str(options.wavelength) == str(components[1]):

                        if len(components) > 2:
                            ratio = components[2] 

                        #run over MC branch
                        if options.type == "sim":
                            generate_mc_output(infile,options.fibre,options.wavelength,ratio)

                        #run over EV branch
                        elif options.type == "data":
                            generate_data_output(infile,options.fibre,options.wavelength,ratio)
                        else:
                            print "Invalid data type!"
                            raise Exception

                    else:
                        continue
