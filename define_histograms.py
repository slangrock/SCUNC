#!/usr/bin/env python
import ROOT

class DefineHistograms:
    """Class to define all histograms to be filled for SMELLIE analysis

    Attributes: 
        t_res (TH1D) : time residual histogram for all events
        t_res_beam (TH1D) : time residual histogram for direct beam light
        t_res_double (TH1D) : time residual histogram for late pulses
        t_res_psup (TH1D) : time residual histogram for PSUP reflections
        t_res_avin (TH1D) : time residual histogram for far side AV reflections
        t_res_avout (TH1D) : time residual histogram for near side AV reflections
        t_res_scatt (TH1D) : time residual histogram for scattered photons
        t_res_multi (TH1D) : time residual histogram for multiple optical interactions
        angle_time (TH2D) : PMT angle with respect to the fibre direction vs time 
                            residual for all events
        angle_time_beam (TH2D) : PMT angle with respect to the fibre direction 
                                 vs time residual for direct beam light
        angle_time_double (TH2D) : PMT angle with respect to the fibre direction 
                                 vs time residual for late pulses
        angle_time_scatt (TH2D) : PMT angle with respect to the fibre direction 
                                  vs time residual for scattered photons
        angle_time_avin (TH2D) : PMT angle with respect to the fibre direction
                                 vs time residual for far side AV reflections
        angle_time_avout (TH2D) : PMT angle with respect to the fibre direction
                                  vs time residual for near side AV reflections
        angle_time_psup (TH2D) : PMT angle with respect to the fibre direction
                                 vs time residual for PSUP reflections
        angle_time_multi (TH2D) : PMT angle with respect to the fibre direction
                                  vs time residual for multiple optical interactions
        z_time (TH2D) : PMT z coordinate vs time residual for all events
        z_time_beam (TH2D) : PMT z coordinate vs time residual for direct beam light
        z_time_double (TH2D) : PMT z coordinate vs time residual for late pulses
        z_time_scatt (TH2D) : PMT z coordinate vs time residual for scattered photons
        z_time_avin (TH2D) : PMT z coordinate vs time residual for far side AV reflections
        z_time_avout (TH2D) : PMT z coordinate vs time residual for near side AV reflecitons
        z_time_psup (TH2D) : PMT z coordinate vs time residual for PSUP reflections
        z_time_multi (TH2D) : PMT z coordinate vs time residual for multiple optical interactions
        theta_phi (TH2D) : PMT theta coordinate vs PMT phi coordinate for all events
        theta_phi_beam (TH2D) : PMT theta coordinate vs PMT phi coordinate for direct beam light
        theta_phi_double (TH2D) : PMT theta coordinate vs PMT phi coordinate for late pulses
        theta_phi_scatt (TH2D) : PMT theta coordinate vs PMT phi coordinate for scattered photons
        theta_phi_avin (TH2D) : PMT theta coordinate vs PMT phi coordinate for far side AV reflections
        theta_phi_avout (TH2D) : PMT theta coordinate vs PMT phi coordinate for near side AV reflections
        theta_phi_psup (TH2D) : PMT theta coordinate vs PMT phi coordinate for PSUP reflections
        theta_phi_multi (TH2D) : PMT theta coordinate vs PMT phi coordinate for multiple interactions
        h_theta (TH1D) : PMT theta coordinate for all events
        h_theta_beam (TH1D) : PMT theta coordinate for direct beam light
        h_theta_double (TH1D) : PMT theta coordinate for late pulses
        h_theta_scatt (TH1D) : PMT theta coordinate for scattered photons
        h_theta_avin (TH1D) : PMT theta coordinate for far side AV reflections
        h_theta_avout (TH1D) : PMT theta coordinate for near side AV reflections
        h_theta_psup (TH1D) : PMT theta coordinate for PSUP reflections
        h_theta_multi (TH1D) : PMT theta coordinate for multiple optical interactions
        h_phi (TH1D) : PMT phi coordinate for all events
        h_phi_beam (TH1D) : PMT phi coordinate for direct beam light
        h_phi_double (TH1D) : PMT phi coordinate for late pulses
        h_phi_scatt (TH1D) : PMT phi coordinate for scattered photons
        h_phi_avin (TH1D) : PMT phi coordinate for far side AV reflections
        h_phi_avout (TH1D) : PMT phi coordinate for near side AV reflections
        h_phi_psup (TH1D) : PMT phi coordinate for PSUP reflections
        h_phi_multi (TH1D) : PMT phi coordinate for multiple optical interactions
        h_nhits (TH1D) : Number of hits for all events

    """


    def __init__(self):
        self.t_res = ROOT.TH1D("time_residual", "", 600, -20.0, 500.0)
        self.t_res_beam = ROOT.TH1D("time_residual_beam", "", 500, -20.0, 300.0)
        self.t_res_double = ROOT.TH1D("time_residual_double", "", 500, -20.0, 300.0)
        self.t_res_psup = ROOT.TH1D("time_residual_psup", "", 500, -20.0, 300.0)
        self.t_res_avin = ROOT.TH1D("time_residual_avin", "", 500, -20.0, 300.0)
        self.t_res_avout = ROOT.TH1D("time_residual_avout", "", 500, -20.0, 300.0)
        self.t_res_scatt = ROOT.TH1D("time_residual_scatt", "", 500, -20.0, 300.0)
        self.t_res_multi = ROOT.TH1D("time_residual_multi", "", 500, -20.0, 300.0)
        self.angle_time_beam = ROOT.TH2D("angle_time_beam", "", 150, -20.0, 300.0, 100, 0.0, 200.0)
        self.angle_time_double = ROOT.TH2D("angle_time_double", "", 150, -20.0, 300.0, 100, 0.0, 200.0)
        self.angle_time = ROOT.TH2D("angle_time", "", 150, -20.0, 300.0, 100, 0.0, 200.0)
        self.angle_time_scatt = ROOT.TH2D("angle_time_scatt", "", 150, -20.0, 300.0, 100, 0.0, 200.0)
        self.angle_time_avin = ROOT.TH2D("angle_time_avin", "", 150, -20.0, 300.0, 100, 0.0, 200.0)
        self.angle_time_avout = ROOT.TH2D("angle_time_avout", "", 150, -20.0, 300.0, 100, 0.0, 200.0)
        self.angle_time_psup = ROOT.TH2D("angle_time_psup", "", 150, -20.0, 300.0, 100, 0.0, 200.0)
        self.angle_time_multi = ROOT.TH2D("angle_time_multi", "", 150, -20.0, 300.0, 100, 0.0, 200.0)
        self.z_time = ROOT.TH2D("z_time", "", 150, -20.0, 300.0, 1000, -10000.0, 10000.0)
        self.z_time_beam = ROOT.TH2D("z_time_beam", "", 150, -20.0, 300.0, 1000, -10000.0, 10000.0)
        self.z_time_double = ROOT.TH2D("z_time_double", "", 150, -20.0, 300.0, 1000, -10000.0, 10000.0)
        self.z_time_scatt = ROOT.TH2D("z_time_scatt", "", 150, -20.0, 300.0, 1000, -10000.0, 10000.0)
        self.z_time_avin = ROOT.TH2D("z_time_avin", "", 150, -20.0, 300.0, 1000, -10000.0, 10000.0)
        self.z_time_avout = ROOT.TH2D("z_time_avout", "", 150, -20.0, 300.0, 1000, -10000.0, 10000.0)
        self.z_time_psup = ROOT.TH2D("z_time_psup", "", 150, -20.0, 300.0, 1000, -10000.0, 10000.0)
        self.z_time_multi = ROOT.TH2D("z_time_multi", "", 150, -20.0, 300.0, 1000, -10000.0, 10000.0)
        self.theta_phi = ROOT.TH2D("theta_phi", "", 100, -3.5, 3.5, 100, 0, 3.5)
        self.theta_phi_beam = ROOT.TH2D("theta_phi_beam", "", 100, -3.5, 3.5, 100, 0, 3.5)
        self.theta_phi_double = ROOT.TH2D("theta_phi_double", "", 100, -3.5, 3.5, 100, 0, 3.5)
        self.theta_phi_scatt = ROOT.TH2D("theta_phi_scatt", "", 100, -3.5, 3.5, 100, 0, 3.5)
        self.theta_phi_avin = ROOT.TH2D("theta_phi_avin", "", 100, -3.5, 3.5, 100, 0, 3.5)
        self.theta_phi_avout = ROOT.TH2D("theta_phi_avout", "", 100, -3.5, 3.5, 100, 0, 3.5)
        self.theta_phi_psup = ROOT.TH2D("theta_phi_psup", "", 100, -3.5, 3.5, 100, 0, 3.5)
        self.theta_phi_multi = ROOT.TH2D("theta_phi_multi", "", 100, -3.5, 3.5, 100, 0, 3.5)
        self.h_theta =  ROOT.TH1D("theta", "", 100, 0.0, 3.5)
        self.h_theta_beam =  ROOT.TH1D("theta_beam", "", 100, 0.0, 3.5)
        self.h_theta_double =  ROOT.TH1D("theta_double", "", 100, 0.0, 3.5)
        self.h_theta_scatt =  ROOT.TH1D("theta_scatt", "", 100, 0.0, 3.5)
        self.h_theta_avin =  ROOT.TH1D("theta_avin", "", 100, 0.0, 3.5)
        self.h_theta_avout =  ROOT.TH1D("theta_avout", "", 100, 0.0, 3.5)
        self.h_theta_psup =  ROOT.TH1D("theta_psup", "", 100, 0.0, 3.5)
        self.h_theta_multi =  ROOT.TH1D("theta_multi", "", 100, 0.0, 3.5)
        self.h_phi =  ROOT.TH1D("phi", "", 100, -3.5, 3.5)
        self.h_phi_beam =  ROOT.TH1D("phi_beam", "", 100, -3.5, 3.5)
        self.h_phi_double =  ROOT.TH1D("phi_double", "", 100, -3.5, 3.5)
        self.h_phi_scatt =  ROOT.TH1D("phi_scatt", "", 100, -3.5, 3.5)
        self.h_phi_avin =  ROOT.TH1D("phi_avin", "", 100, -3.5, 3.5)
        self.h_phi_avout =  ROOT.TH1D("phi_avout", "", 100, -3.5, 3.5)
        self.h_phi_psup =  ROOT.TH1D("phi_psup", "", 100, -3.5, 3.5)
        self.h_phi_multi =  ROOT.TH1D("phi_multi", "", 100, -3.5, 3.5)
        self.h_nhits = ROOT.TH1D("number_of_hits","",200,0.0,1000.0)



