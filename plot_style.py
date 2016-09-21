#!/usr/bin/env python
import ROOT       

class PlotStyle:
        """Class to define the style for all plots produced by the code framework.

         Attributes:
          hipStyle (TStyle) : object to define special plot styles
        """
    
        def __init__(self):
            self.hipStyle = ROOT.TStyle("clearRetro","HIP plots style for publications")

        def set_style(self):
            """Function to define plot style for the SMELLIE analysis.

             Returns:
              plot style
            """
            # use plain black on white colors
            self.hipStyle.SetFrameBorderMode(0)
            self.hipStyle.SetCanvasBorderMode(0)
            self.hipStyle.SetPadBorderMode(0)
            self.hipStyle.SetPadBorderSize(0)
            self.hipStyle.SetPadColor(0)
            self.hipStyle.SetCanvasColor(0)
            self.hipStyle.SetTitleColor(0)
            self.hipStyle.SetStatColor(0)
            self.hipStyle.SetFillColor(0)

            # use bold lines 
            self.hipStyle.SetHistLineWidth(2)
            self.hipStyle.SetLineWidth(2)

            # no title, stats box or fit as default
            self.hipStyle.SetOptTitle(0)
            self.hipStyle.SetOptStat(0)
            self.hipStyle.SetOptFit(0)

            # postscript dashes
            self.hipStyle.SetLineStyleString(2,"[12 12]") # postscript dashes

            # text style and size
            self.hipStyle.SetTextFont(61)
            self.hipStyle.SetTextSize(0.24)
            self.hipStyle.SetLabelFont(61,"x")
            self.hipStyle.SetLabelFont(61,"y")
            self.hipStyle.SetLabelFont(61,"z")
            #self.hipStyle.SetTitleFont(61,"x")
            #self.hipStyle.SetTitleFont(61,"y")
            self.hipStyle.SetTitleFont(61,"z")
            self.hipStyle.SetLabelSize(0.04,"x")
            self.hipStyle.SetTitleSize(0.05,"x")
            self.hipStyle.SetTitleColor(1,"x")
            self.hipStyle.SetLabelSize(0.04,"y")
            self.hipStyle.SetTitleSize(0.05,"y")
            self.hipStyle.SetTitleColor(1,"y")
            self.hipStyle.SetLabelSize(0.04,"z")
            self.hipStyle.SetTitleSize(0.05,"z")
            self.hipStyle.SetTitleColor(1,"z")
  
            # AXIS OFFSETS
            self.hipStyle.SetTitleOffset(0.8,"x")
            self.hipStyle.SetTitleOffset(0.8,"y")
            self.hipStyle.SetTitleOffset(0.8,"z")

            # Legends
            self.hipStyle.SetLegendBorderSize(1)
            # graphs - set default martker to cross, rather than .
            self.hipStyle.SetMarkerStyle(2)  # cross +

            return self.hipStyle
