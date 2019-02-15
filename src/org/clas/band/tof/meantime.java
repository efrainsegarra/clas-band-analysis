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
		
		
		
		double Ebeam = 10.6;//GeV
		double mp = 0.93827; //GeV
		double mtar = mp;
		double rad2deg = 180./3.14159;

		double cut_ep = 2; //GeV
		double cut_chi2pid = 5;
		double cut_min_vz = -15;
		double cut_max_vz = 10;
		double cut_W = 0;

		
		// ----------------------------------------------------------------------------------
		// Declaring histograms
		H1F h1_meantime_fadc = new H1F("h1_meantime_fadc","",50,0,500); 	PrettyH1F(h1_meantime_fadc ,"(L+R)/2 for FADC"  ,"Counts",1);
		H1F h1_meantime_tdc  = new H1F("h1_meantime_tdc","",50,400,900); 	PrettyH1F(h1_meantime_tdc ,"(L+R)/2 for TDC"  ,"Counts",1);
		H1F h1_tdiff_fadc = new H1F("h1_tdiff_fadc","",100,-40,40); 		PrettyH1F(h1_tdiff_fadc ,"L-R for FADC"  ,"Counts",1);
		H1F h1_tdiff_tdc  = new H1F("h1_tdiff_tdc","",100,-40,40);			PrettyH1F(h1_tdiff_tdc ,"L-R for TDC"  ,"Counts",1);
		H1F h1_adcL = new H1F("h1_adcL","",1000,0,30000);					PrettyH1F(h1_adcL ,"L ADC"  ,"Counts",1);
		H1F h1_adcR = new H1F("h1_adcR","",1000,0,30000);					PrettyH1F(h1_adcR ,"R ADC"  ,"Counts",1);
		H1F h1_ToF_fadc = new H1F("h1_ToF_fadc","",60,-100,500); 			PrettyH1F(h1_ToF_fadc ,"(L+R)/2 - RF for FADC"  ,"Counts",1);
		H1F h1_ToF_tdc = new H1F("h1_ToF_tdc","",60,300,900); 				PrettyH1F(h1_ToF_tdc ,"(L+R)/2 - RF for TDC"  ,"Counts",1);

		
		// ----------------------------------------------------------------------------------
		// Opening HIPO file
		HipoReader  reader = new HipoReader();
		for(String arg: args) {
			String dataFile = arg;

			//String dataFile = "/Users/efrainsegarra/Documents/band/prod_data/out_clas_006164.evio.00000.hipo"; // target lH2
			reader.open(dataFile);
			GenericKinematicFitter fitter = new GenericKinematicFitter(Ebeam);
		
			// ----------------------------------------------------------------------------------
			// Loop over events and print them on the screen
			while(reader.hasNext()){
				HipoEvent event = (HipoEvent) reader.readNextEvent();
	
				if(		(event.hasGroup("BAND::hits"))&&
						(event.hasGroup("RUN::rf"))
						){
					
					HipoGroup rfBank = event.getGroup("RUN::rf");
					
					double trf = 0.0;
					for (int rfIdx = 0; rfIdx < rfBank.getNode("id").getDataSize(); rfIdx++) {
						if (rfBank.getNode("id").getShort(rfIdx) == 1) {
							trf = rfBank.getNode("time").getFloat(rfIdx);
						}
					}
					
					
					HipoGroup band_hits    = event.getGroup("BAND::hits");
					
				
					int nHits = band_hits.getNode("id").getDataSize();
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
						
						if( adcLcorr < 2000 || adcRcorr < 2000 ) continue;
						/*
						float tFadcLcorr = band_hits.getNode().getFloat("tFadcLcorr", hit);
						float tFadcRcorr = band_hits.getNode().getFloat("tFadcRcorr", hit);
						
						
						float tTdcLcorr = band_hits.getFloat("tTdcLcorr", hit);
						float tTdcRcorr = band_hits.getFloat("tTdcRcorr", hit);
						
						float x = band_hits.getFloat("x", hit);
						float y = band_hits.getFloat("y", hit);
						float z = band_hits.getFloat("z", hit);
						
						float ux = band_hits.getFloat("ux", hit);
						float uy = band_hits.getFloat("uy", hit);
						float uz = band_hits.getFloat("uz", hit);
						*/				
						
						if( sector != 2 || layer != 5 ) continue;
						
						// Fill histograms
						h1_meantime_fadc.fill(meantimeFadc);
						h1_meantime_tdc.fill(meantimeTdc);
						h1_tdiff_fadc.fill(difftimeFadc);
						h1_tdiff_tdc.fill(difftimeTdc);
						h1_adcL.fill(adcLcorr);
						h1_adcR.fill(adcRcorr);
						h1_ToF_fadc.fill(meantimeFadc-trf);
						h1_ToF_tdc.fill(meantimeTdc-trf);
						
						
						
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
		c0.cd(7);	c0.draw(h1_ToF_tdc);
 		
	}
	// =========================================================================================================
		public static void PrettyH1F(H1F h1,String titx,String tity,int color) {
			h1.setTitleX(titx);
			h1.setTitleY(tity);
			h1.setLineColor(color);
			h1.setLineWidth(3);
		}
}
