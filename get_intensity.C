#include <TH1D.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TFile.h>
#include <TF1.h>
#include <TMarker.h>

#include <iostream>
#include <fstream>
#include <string>
#include <math.h> 
using namespace std;

//functions to evaluate the distance between two functions for each point of the function - used to determine the intersection between two functions
TF1 *f1, *f2, *f1_up, *f1_down, *f2_data_up, *f2_data_down;
double finter(double *x, double*par) {
   return TMath::Abs(f1->EvalPar(x,par) - f2->EvalPar(x,par));
}

double finter_up(double *x, double*par) {
   return TMath::Abs(f1_up->EvalPar(x,par) - f2->EvalPar(x,par));
}

double finter_down(double *x, double*par) {
   return TMath::Abs(f1_down->EvalPar(x,par) - f2->EvalPar(x,par));
}

double finter_data_up(double *x, double*par) {
   return TMath::Abs(f1->EvalPar(x,par) - f2_data_up->EvalPar(x,par));
}

double finter_data_down(double *x, double*par) {
   return TMath::Abs(f1->EvalPar(x,par) - f2_data_down->EvalPar(x,par));
}


void main(string fibre,string wl) {
   /* 
      Function to measure the intensity of SMELLIE data using the txt files returned by intensity_measurement.py

      Arguments:
       fibre (string) : fibre label
       wl (dtring) : wavelength the run was produced at
   */

   string filename;
   string filename_sim;
   string filename_data;
   string line;

   //define location and name of the data and simulation text files returned by intensity_measurement.py
   filename += "/data/langrock/rat-5.0-SMELLIE_analysis/";
   filename += fibre;
   filename += "/";
   filename += wl;
   filename_sim += filename;
   filename_sim += "_sim_intensity_fake.txt";

   cout << filename_sim << " ";

   filename_data += filename;
   filename_data += "_data_intensity_fake.txt";

   cout << filename_data << endl;

   //open simulation text file
   ifstream simfile (filename_sim.c_str());
   if (simfile.is_open())
   {
     while ( getline (simfile,line) )
       {
         //get the fit parameters determined by intensity_measurement.py
         stringstream   linestream(line);
         string fibrename;
         double wavelength;
         double p0;
         double p1;
         double p2;
         double p0_up;
         double p1_up;
         double p2_up;
         double p0_down;
         double p1_down;
         double p2_down;

         linestream >> fibrename >> wavelength;
         linestream >> p0 >> p1 >> p2;
         linestream >> p0_up >> p1_up >> p2_up;
         linestream >> p0_down >> p1_down >> p2_down;

         //define functions and set parameters to the fit parameters from file
         f1 = new TF1("f1","[0]*log([1]+[2]*x)",0,10000);
         //f1 = new TF1("f1","[0]*x+[1]",0,1000);
         f1->SetParameter(0,p0); 
         f1->SetParameter(1,p1); 
         f1->SetParameter(2,p2); 

         f1_up = new TF1("f1_up","[0]*log([1]+[2]*x)",0,10000);
         //f1_up = new TF1("f1_up","[0]*x+[1]",0,1000);
         f1_up->SetParameter(0,p0_up); 
         f1_up->SetParameter(1,p1_up); 
         f1_up->SetParameter(2,p2_up); 

         f1_down = new TF1("f1_down","[0]*log([1]+[2]*x)",0,10000);
         //f1_down = new TF1("f1_down","[0]**x+[1]",0,1000);
         f1_down->SetParameter(0,p0_down); 
         f1_down->SetParameter(1,p1_down); 
         f1_down->SetParameter(2,p2_down); 

         //define output file the measured intensities are saved to
         ofstream output;
         string output_name;
         output_name += "/users/langrock/plotting_macros/SMELLIE/SMELLIE_analysis_framework/datatables/intensity_";
         output_name += fibre;
         output_name += "_";
         output_name += wl;
         output_name += "_fake.txt";

         //open data text file
         output.open(output_name.c_str());

         string line_data;
         ifstream datafile (filename_data.c_str());
         if (datafile.is_open())
           {
             while ( getline (datafile,line_data) )
               {
                 //get nhits and standard deviation as determined by intensity_measurement.py
                 stringstream   linestream(line_data);
                 string fibrename_data;
                 double wavelength_data;
                 string intensity;
                 double n_hits;
                 double rms;

                 linestream >> fibrename_data >> wavelength_data;
                 linestream >> intensity >> n_hits >> rms;

                 cout << fibrename_data << " " << wavelength_data << " " << intensity << " " << n_hits << " " << rms << endl;

                 double n_hits_up;
                 double n_hits_down;

                 n_hits_up = n_hits + rms;
                 n_hits_down = n_hits - rms;

                 //define linear functions with the nhits, nhits + standard deviation and nhits - standard deviation as the sole parameter 
                 f2 = new TF1("f2","[0]",0,10000);
                 f2->SetParameter(0,n_hits);

                 f2_data_up = new TF1("f2_data_up","[0]",0,10000);
                 f2_data_up->SetParameter(0,n_hits_up);

                 f2_data_down = new TF1("f2_data_down","[0]",0,10000);
                 f2_data_down->SetParameter(0,n_hits_down);

                 //draw fit parameter function and nhits function
                 TCanvas *c1 = new TCanvas();
                 c1->Draw();

                 f1->SetTitle("");
                 f1->GetXaxis()->SetTitle("Number of photons per beam pulse");
                 f1->GetYaxis()->SetTitle("Mean number of hits");
                 f1->SetLineColor(50);
                 f1->SetLineWidth(2);
                 f2->SetLineWidth(2);
                 f1->Draw();
                 f2->Draw("same");

                 TLegend *l1 = new TLegend(0.53,0.45,0.8955,0.7,"","brNDC");
                 l1->AddEntry(f1, "Fit on MC simulations","l");
                 l1->AddEntry(f2, "Fake data","l");

                 l1->Draw("same");

                 //get intersection between the two functions
                 TF1 *fint = new TF1("fint",finter,0,10000,0);
                 double xint = fint->GetMinimumX();
                 TMarker *m = new TMarker(xint,f1->Eval(xint),24);
                 m->SetMarkerColor(kRed);
                 m->SetMarkerSize(3);
                 m->Draw();

                 ostringstream int_string;
                 int_string << intensity;
                 string output_plot;
                 output_plot = "/data/langrock/rat-5.0-SMELLIE_analysis/";
                 output_plot += fibre;
                 output_plot += "/plots/";
                 output_plot += wl;
                 output_plot += "new_";
                 output_plot += int_string.str();
                 output_plot += "_intensity_fit_water_log.pdf";
                
                 //save canvas  
                 c1->Print(output_plot.c_str());

                 c1->Close();

                 //draw fit parameter up function and nhits function
                 TCanvas *c2 = new TCanvas();
                 c2->Draw();
                   
                 f1_up->SetTitle("");
                 f1_up->GetXaxis()->SetTitle("Number of photons per beam pulse");
                 f1_up->GetYaxis()->SetTitle("Mean number of hits");
                 f1_up->SetLineColor(50);
                 f1_up->SetLineWidth(2);
                 f2->SetLineWidth(2);
                 f1_up->Draw();
                 f2->Draw("same");

                 l1->Draw("same");

                 //get intersection between the two functions
                 TF1 *fint_up = new TF1("fint_up",finter_up,0,10000,0);
                 double xint_up = fint_up->GetMinimumX();
                 TMarker *m_up = new TMarker(xint_up,f1_up->Eval(xint_up),24);
                 m_up->SetMarkerColor(kRed);
                 m_up->SetMarkerSize(3);
                 m_up->Draw();

                 c2->Close();

                 //draw fit parameter down function and nhits function
                 TCanvas *c3 = new TCanvas();
                 c3->Draw();
                   
                 f1_down->SetTitle("");
                 f1_down->GetXaxis()->SetTitle("Number of photons per beam pulse");
                 f1_down->GetYaxis()->SetTitle("Mean number of hits");
                 f1_down->SetLineColor(50);
                 f1_down->SetLineWidth(2);
                 f2->SetLineWidth(2);
                 f1_down->Draw();
                 f2->Draw("same");

                 l1->Draw("same");

                 //get intersection between the two functions
                 TF1 *fint_down = new TF1("fint_down",finter_down,0,10000,0);
                 double xint_down = fint_down->GetMinimumX();
                 TMarker *m_down = new TMarker(xint_down,f1_down->Eval(xint_down),24);
                 m_down->SetMarkerColor(kRed);
                 m_down->SetMarkerSize(3);
                 m_down->Draw();

                 c3->Close();

                 //draw fit parameter function and nhits up function
                 TCanvas *c4 = new TCanvas();
                 c4->Draw();

                 f1->Draw();
                 f2_data_up->Draw("same");

                 l1->Draw("same");

                 //get intersection between the two functions
                 TF1 *fint_data_up = new TF1("fint_data_up",finter_data_up,0,10000,0);
                 double xint_data_up = fint_data_up->GetMinimumX();
                 TMarker *m_data_up = new TMarker(xint_data_up,f1->Eval(xint_data_up),24);
                 m_data_up->SetMarkerColor(kRed);
                 m_data_up->SetMarkerSize(3);
                 m_data_up->Draw();
 
                 c4->Close();

                 //draw fit parameter function and nhits down function
                 TCanvas *c5 = new TCanvas();
                 c5->Draw();

                 f1->Draw();
                 f2_data_down->Draw("same");

                 l1->Draw("same");

                 //get intersection between the two functions
                 TF1 *fint_data_down = new TF1("fint_data_down",finter_data_down,0,10000,0);
                 double xint_data_down = fint_data_down->GetMinimumX();
                 TMarker *m_data_down = new TMarker(xint_data_down,f1->Eval(xint_data_down),24);
                 m_data_down->SetMarkerColor(kRed);
                 m_data_down->SetMarkerSize(3);
                 m_data_down->Draw();

                 c5->Close();

                 //calculate uncertainty on the measurement
                 double sigma_up;
                 double sigma_down;
                 double sigma_data_up;
                 double sigma_data_down;

                 sigma_up = (xint - xint_up)*(xint - xint_up);
                 sigma_down = (xint - xint_down)*(xint - xint_down);
                 sigma_data_up = (xint - xint_data_up)*(xint - xint_data_up);
                 sigma_data_down = (xint - xint_data_down)*(xint - xint_data_down);

                 double sigma_sum;
                 sigma_sum = (sigma_up + sigma_down + sigma_data_up + sigma_data_down)/4;
                 
                 double sigma;
                 sigma = sqrt (sigma_sum);

                 //save measured intensity and uncertainty for that data run
                 output << intensity << " \t" << xint << " +/- " << sigma << "\n";

               }
             datafile.close();
             output.close();
           }
         else cout << "Unable to open data file! \n"; 

       }
     simfile.close();
   }
   else cout << "Unable to open sim file! \n"; 
}

