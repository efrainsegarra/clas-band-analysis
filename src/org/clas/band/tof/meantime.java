package org.clas.band.tof;
import java.io.Reader;

import org.jlab.groot.data.H1F;
import org.jlab.groot.data.H2F;
import org.jlab.groot.ui.TCanvas;

import org.jlab.jnp.hipo.data.HipoEvent;
import org.jlab.jnp.hipo.data.HipoGroup;
import org.jlab.jnp.hipo.data.HipoNode;
import org.jlab.jnp.hipo.io.HipoReader;

import org.jlab.clas.physics.*;
import org.jlab.io.base.DataBank;
import org.jlab.io.base.DataEvent;
import org.jlab.io.evio.*;
import org.jlab.io.hipo.HipoDataSource;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Scanner;

public class meantime {

	// =========================================================================================================
	public static void main(String[] args){		


		// ----------------------------------------------------------------------------------
		// Useful variables
		double mp      = 0.93827; //GeV
		double mtar    = mp;
		double rad2deg = 180./3.14159;
		double Ebeam = 10.6; //GeV
		double c = 29.9792;

		// ----------------------------------------------------------------------------------
		// Event selection cuts
		double cut_ep      =     2; //GeV
		double cut_chi2pid =     5;
		double cut_min_vz  =   -15; //cm
		double cut_max_vz  =    10; //cm
		double cut_W       =     0; //GeV
		double cut_uvw     =    15; //cm
		double cut_Epcal   = 0.060; //GeV (60 MeV)
		double cut_tof_e   =    10; //ns


		// ----------------------------------------------------------------------------------
		// Declaring histograms
		H1F h1_meantime_fadc = new H1F("h1_meantime_fadc","",50,125,225); 	PrettyH1F(h1_meantime_fadc ,"(L+R)/2 for FADC"  ,"Counts",1);
		H1F h1_meantime_tdc  = new H1F("h1_meantime_tdc","",250,400,900); 	PrettyH1F(h1_meantime_tdc ,"(L+R)/2 for TDC"  ,"Counts",1);
		H1F h1_tdiff_fadc = new H1F("h1_tdiff_fadc","",100,-40,40); 		PrettyH1F(h1_tdiff_fadc ,"L-R for FADC"  ,"Counts",1);
		H1F h1_tdiff_tdc  = new H1F("h1_tdiff_tdc","",100,-40,40);			PrettyH1F(h1_tdiff_tdc ,"L-R for TDC"  ,"Counts",1);
		H1F h1_adcL = new H1F("h1_adcL","",1000,0,30000);					PrettyH1F(h1_adcL ,"L ADC"  ,"Counts",1);
		H1F h1_adcR = new H1F("h1_adcR","",1000,0,30000);					PrettyH1F(h1_adcR ,"R ADC"  ,"Counts",1);
		H1F h1_ToF_fadc = new H1F("h1_ToF_fadc","",225,-50,400); 			PrettyH1F(h1_ToF_fadc ,"(L+R)/2 - RF for FADC"  ,"Counts",1);
		H1F h1_ToF_tdc = new H1F("h1_ToF_tdc","",750,0,1500); 				PrettyH1F(h1_ToF_tdc ,"(L+R)/2 - RF for TDC"  ,"Counts",1);
		
		H1F h1_bar_nPho = new H1F("h1_bar_nPho","",	652,110,652);	PrettyH1F(h1_bar_nPho,"Bar ID","Number Hits Between 5-15 ns",1);


		// ----------------------------------------------------------------------------------
		// Opening HIPO file
		HipoReader reader = new HipoReader();
		for(String arg: args) {
			String dataFile = arg;

			//String dataFile = "/Users/efrainsegarra/Documents/band/prod_data/out_clas_006164.evio.00000.hipo"; // target lH2
			reader.open(dataFile);
			GenericKinematicFitter fitter = new GenericKinematicFitter(Ebeam);
			int cnt = 0;
			// ----------------------------------------------------------------------------------
			// Loop over events and print them on the screen
			while(reader.hasNext()){
				HipoEvent event = (HipoEvent) reader.readNextEvent();
				cnt++;
				if( cnt % 10000 == 0 ) System.out.println("Working on event: "+cnt);
				if(     (event.hasGroup("REC::Particle"    ))&&
						(event.hasGroup("REC::Calorimeter" ))&&
						(event.hasGroup("REC::Scintillator"))&&
						(event.hasGroup("REC::Event"       ))&&
						(event.hasGroup("BAND::hits"	   ))
				  ){	

					HipoGroup bank_particle     = event.getGroup("REC::Particle"    );
					HipoGroup bank_calorimeter  = event.getGroup("REC::Calorimeter" );
					HipoGroup bank_scintillator = event.getGroup("REC::Scintillator");
					HipoGroup bank_event        = event.getGroup("REC::Event"       );


					// first particle variables

					// Particle bank
					int pid0       = bank_particle    .getNode("pid"    ).getInt  (0);		// electron candidate id assigned by clas
					int chr0       = bank_particle    .getNode("charge" ).getInt  (0);		// electron candidate charge
					double chi2pid = bank_particle    .getNode("chi2pid").getFloat(0);		// electron candidate goodness of PID from EB
					double epx     = bank_particle    .getNode("px"     ).getFloat(0);		// electron candidate momentum x-component [GeV]
					double epy     = bank_particle    .getNode("py"     ).getFloat(0);		// electron candidate momentum y-component [GeV]
					double epz     = bank_particle    .getNode("pz"     ).getFloat(0);		// electron candidate momentum z-component [GeV]
					double evx     = bank_particle    .getNode("vx"     ).getFloat(0);		// electron candidate vertex x coordinate [cm]
					double evy     = bank_particle    .getNode("vy"     ).getFloat(0);		// electron candidate vertex y coordinate [cm]
					double evz     = bank_particle    .getNode("vz"     ).getFloat(0);		// electron candidate vertex z coordinate [cm]

					// Event bank
					double t_vtx   = bank_event       .getNode("STTime" ).getFloat(0);		// event vertex time

					// Scintillator bank
					double t_e     = bank_scintillator.getNode("time"   ).getFloat(0);		// electron candidate time at FTOF

					short cal_id = 0;
					double Epcal = 0;
					double Eecin = 0;
					double Eecou = 0;

					// Calorimeter bank
					int nCal = bank_calorimeter.getNode("pindex").getDataSize();

					double lU      = bank_calorimeter.getNode("lu"    ).getFloat(0);	// electron candidate distance on U-side [cm?]
					double lV      = bank_calorimeter.getNode("lv"    ).getFloat(0);	// electron candidate distance on V-side [cm?]
					double lW      = bank_calorimeter.getNode("lw"    ).getFloat(0);	// electron candidate distance on W-side [cm?]

					for(int cal=0 ; cal < nCal ; cal++) {
						cal_id = bank_calorimeter.getNode("pindex").getShort(cal);
						int layer = bank_calorimeter.getNode("layer").getByte(cal);
						if     (cal_id==0&&layer==1) Epcal = bank_calorimeter.getNode("energy").getFloat(cal);	// PCal
						else if(cal_id==0&&layer==4) Eecin = bank_calorimeter.getNode("energy").getFloat(cal);	// ECinner
						else if(cal_id==0&&layer==7) Eecou = bank_calorimeter.getNode("energy").getFloat(cal);	// ECouter
					}

					double Ee      = Epcal + Eecin + Eecou;	// electron candidate energy from calorimeter [GeV]

					// calculated variables
					double ep     = Math.sqrt(epx*epx + epy*epy + epz*epz);		// electron candidate momentum magnitude [GeV]
					Vector3 v3_ep = new Vector3(epx,epy,epz);					// electron candidate momentum vector [GeV]
					double th_e   = v3_ep.theta();//Math.acos(epz/ep);			// electron candidate theta [rad]
					double phi_e  = v3_ep.phi();								// electron candidate phi [rad]
					double Q2     = 4*ep*Ebeam*Math.pow(Math.sin(th_e/2.),2);	// Q-squared [GeV^2]
					double nu     = Ebeam - ep;									// Transfer energy [GeV]
					double W2     = mtar*mtar-Q2+2*nu*mtar;						// Invariant mass ^2 [GeV]
					double xB     = Q2/2./mp/nu;								// Bjorken-x

					double tof_e  = t_e - t_vtx;								// electron candidate time-of-flight [ns]

					// Transfer variables
					double qx = - epx;
					double qy = - epy;
					double qz = Ebeam - epz;

					// -------------------------------------------------------------------------
					// Only keep events for which the first particle is an electron
					if(     (pid0!=11            )||
							(chr0!=-1            )
							//(chi2pid>=cut_chi2pid)||
							//(ep<=cut_ep          )||
							//(ep>=Ebeam           )||
							//(evz>cut_max_vz      )||
							//(evz<cut_min_vz      )||
							//(lU<cut_uvw          )||
							//(lV<cut_uvw          )||
							//(lW<cut_uvw          )||
							//(Epcal<cut_Epcal     )||
							//(Math.sqrt(W2)<=cut_W)||
							//(tof_e<cut_tof_e     )
					  ) continue;






					HipoGroup band_hits    = event.getGroup("BAND::hits");


					int nHits = band_hits.getNode("id").getDataSize();
					if( nHits > 1 ) continue;
					for(int hit = 0; hit < nHits; hit++) {
						int sector = 	(int)band_hits.getNode("sector").getByte(hit);
						int layer = 	(int)band_hits.getNode("layer").getByte(hit);
						int component = (int)band_hits.getNode("component").getShort(hit);

						float meantimeTdc	= band_hits.getNode("meantimeTdc").getFloat(hit);
						float meantimeFadc 	= band_hits.getNode("meantimeFadc").getFloat(hit);
						float difftimeTdc 	= band_hits.getNode("difftimeTdc").getFloat(hit);
						float difftimeFadc 	= band_hits.getNode("difftimeFadc").getFloat(hit);


						float adcLcorr = band_hits.getNode("adcLcorr").getFloat(hit);
						float adcRcorr = band_hits.getNode("adcRcorr").getFloat(hit);

						if( adcLcorr < 4000 || adcRcorr < 4000 ) continue;
						if( sector < 3 || sector > 4 ) continue;

						float tFadcLcorr = band_hits.getNode("tFadcLcorr").getFloat(hit);
						float tFadcRcorr = band_hits.getNode("tFadcRcorr").getFloat(hit);


						float tTdcLcorr = band_hits.getNode("tTdcLcorr").getFloat(hit);
						float tTdcRcorr = band_hits.getNode("tTdcRcorr").getFloat(hit);

						float x = band_hits.getNode("x").getFloat(hit);
						float y = band_hits.getNode("y").getFloat(hit);
						float z = band_hits.getNode("z").getFloat(hit);

						float ux = band_hits.getNode("ux").getFloat(hit);
						float uy = band_hits.getNode("uy").getFloat(hit);
						float uz = band_hits.getNode("uz").getFloat(hit);

						//if( sector == 3 || sector == 4 ) continue;

						double dL = Math.sqrt( Math.pow(x,2) + Math.pow(y,2) + Math.pow(z,2) );
						// Fill histograms
						h1_meantime_fadc.fill(meantimeFadc);
						h1_meantime_tdc.fill(meantimeTdc);
						h1_tdiff_fadc.fill(difftimeFadc);
						h1_tdiff_tdc.fill(difftimeTdc);
						h1_adcL.fill(adcLcorr);
						h1_adcR.fill(adcRcorr);
						h1_ToF_fadc.fill(meantimeFadc-t_vtx);
						h1_ToF_tdc.fill(meantimeTdc-t_vtx);

						if( Math.abs( meantimeFadc-t_vtx-40 - 10 ) < 5 ){
							h1_bar_nPho.fill( layer*100+sector*10+component);
						}

					}// end loop over hits in event

				}//end bank if
			}// end file

		}//end loop over all files

		TCanvas c0 = new TCanvas("c0", 800, 600);
		c0.divide(2, 4);
		c0.cd(0);	c0.draw(h1_adcL);
		c0.cd(1);	c0.draw(h1_adcR);
		c0.cd(2);	c0.draw(h1_tdiff_fadc);
		c0.cd(3);	c0.draw(h1_tdiff_tdc);
		c0.cd(4);	c0.draw(h1_meantime_fadc);
		c0.cd(5);	c0.draw(h1_meantime_tdc);
		c0.cd(6);	c0.draw(h1_ToF_fadc);
		c0.cd(7);	c0.draw(h1_bar_nPho);

	}
	// =========================================================================================================
	public static void PrettyH1F(H1F h1,String titx,String tity,int color) {
		h1.setTitleX(titx);
		h1.setTitleY(tity);
		h1.setLineColor(color);
		h1.setLineWidth(3);
	}
	public static void PrettyH2F(H2F h2,String titx,String tity) {
		h2.setTitleX(titx);
		h2.setTitleY(tity);
	}
}
