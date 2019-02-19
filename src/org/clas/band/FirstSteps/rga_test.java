package org.clas.band.FirstSteps;
import java.io.Reader;

import org.jlab.groot.data.H1F;
import org.jlab.groot.data.H2F;
import org.jlab.groot.ui.TCanvas;

import org.jlab.jnp.hipo.data.HipoEvent;
import org.jlab.jnp.hipo.data.HipoGroup;
import org.jlab.jnp.hipo.data.HipoNode;
import org.jlab.jnp.hipo.io.HipoReader;
import org.jlab.jnp.hipo.io.HipoWriter;

import org.jlab.clas.physics.*;

import org.jlab.io.evio.*;

import java.util.ArrayList;
import java.util.Scanner;

import javax.swing.Box.Filler;

public class rga_test {

	// =========================================================================================================
	public static void main(String[] args){		

		// ----------------------------------------------------------------------------------
		// Useful variables
		double mp      = 0.93827; //GeV
		double mtar    = mp;
		double rad2deg = 180./3.14159;

		// ----------------------------------------------------------------------------------
		// Event selection cuts
		double cut_ep      =     2; //GeV
		double cut_chi2pid =     5;
		double cut_min_vz  =   -15; //cm
		double cut_max_vz  =    10; //cm
		double cut_W       =     0; //GeV
		double cut_uvw     =    15; //cm
		double cut_Epcal   = 0.060; //GeV (60 MeV)

		// ----------------------------------------------------------------------------------
		// Declaring histograms
		// 1D histograms
		H1F h1_e_vz  = new H1F("h1_e_vz"  ,"h1_e_vz"  ,100, -50, 50);	PrettyH1F(h1_e_vz  ,"v_{z} [cm]"           ,"Counts",4);

		H1F h1_e_px  = new H1F("h1_e_px"  ,"h1_e_px"  ,100,  -2,  2);	PrettyH1F(h1_e_px  ,"electron p_{x} [GeV]" ,"Counts",4);
		H1F h1_e_py  = new H1F("h1_e_py"  ,"h1_e_py"  ,100,  -2,  2);	PrettyH1F(h1_e_py  ,"electron p_{y} [GeV]" ,"Counts",4);
		H1F h1_e_pz  = new H1F("h1_e_pz"  ,"h1_e_pz"  ,100,   0,  9);	PrettyH1F(h1_e_pz  ,"electron p_{z} [GeV]" ,"Counts",4);
		H1F h1_e_p   = new H1F("h1_e_p"   ,"h1_e_p"   ,100,   0,  9);	PrettyH1F(h1_e_p   ,"electron |p| [GeV]"   ,"Counts",4);

		H1F h1_e_th  = new H1F("h1_e_th"  ,"h1_e_th"  ,100,   0, 30);	PrettyH1F(h1_e_th  ,"#theta_e [deg]"       ,"Counts",4);
		H1F h1_e_phi = new H1F("h1_e_phi" ,"h1_e_phi" ,100,-190,190);	PrettyH1F(h1_e_phi ,"#phi_e [deg]"         ,"Counts",4);
		H1F h1_e_lu  = new H1F("h1_e_lu"  ,"h1_e_lu"  ,100,   0,500);	PrettyH1F(h1_e_lu  ,"distance on U-side"   ,"Counts",4);
		H1F h1_e_lv  = new H1F("h1_e_lv"  ,"h1_e_lv"  ,100,   0,500);	PrettyH1F(h1_e_lv  ,"distance on V-side"   ,"Counts",4);
		H1F h1_e_lw  = new H1F("h1_e_lw"  ,"h1_e_lw"  ,100,   0,500);	PrettyH1F(h1_e_lw  ,"distance on W-side"   ,"Counts",4);

		H1F h1_p_vz  = new H1F("h1_p_vz"  ,"h1_p_vz"  ,100, -50, 50);	PrettyH1F(h1_p_vz  ,"v_{z} [cm]"           ,"Counts",2);
		H1F h1_p_num = new H1F("h1_p_num" ,"h1_p_num" , 20,   0, 10);	PrettyH1F(h1_p_num ,"proton number"        ,"Counts",4);

		H1F h1_p_px  = new H1F("h1_p_px"  ,"h1_p_px"  ,100,  -2,  2);	PrettyH1F(h1_p_px  ,"proton p_{x} [GeV]"   ,"Counts",4);
		H1F h1_p_py  = new H1F("h1_p_py"  ,"h1_p_py"  ,100,  -2,  2);	PrettyH1F(h1_p_py  ,"proton p_{y} [GeV]"   ,"Counts",4);
		H1F h1_p_pz  = new H1F("h1_p_pz"  ,"h1_p_pz"  ,100,   0,  9);	PrettyH1F(h1_p_pz  ,"proton p_{z} [GeV]"   ,"Counts",4);
		H1F h1_p_p   = new H1F("h1_p_p"   ,"h1_p_p"   ,100,   0,  9);	PrettyH1F(h1_p_p   ,"proton |p| [GeV]"     ,"Counts",4);

		H1F h1_p_th  = new H1F("h1_p_th"  ,"h1_p_th"  ,100,   0, 80);	PrettyH1F(h1_p_th  ,"#theta_p [deg]"       ,"Counts",4);
		H1F h1_p_phi = new H1F("h1_p_phi" ,"h1_p_phi" ,100,-190,190);	PrettyH1F(h1_p_phi ,"#phi_p [deg]"         ,"Counts",4);

		H1F h1_pmiss = new H1F("h1_pmiss" ,"P_{miss}" ,100,   0,  5);	PrettyH1F(h1_pmiss ,"Pm [GeV]"             ,"Counts",4);
		H1F h1_pmx   = new H1F("h1_pmx"   ,"h1_pmx"   ,100,  -2,  2);	PrettyH1F(h1_pmx   ,"Pmx [GeV]"            ,"Counts",4);
		H1F h1_pmy   = new H1F("h1_pmy"   ,"h1_pmy"   ,100,  -2,  2);	PrettyH1F(h1_pmy   ,"Pmy [GeV]"            ,"Counts",4);
		H1F h1_pmz   = new H1F("h1_pmz"   ,"h1_pmz"   ,100,  -3,  3);	PrettyH1F(h1_pmz   ,"Pmz [GeV]"            ,"Counts",4);
		H1F h1_Mmiss = new H1F("h1_Mmiss" ,"h1_Mmiss" ,100,  -1,  2);	PrettyH1F(h1_Mmiss ,"m_{miss} [GeV]"       ,"Counts",4);
		
		H1F h1_W     = new H1F("h1_W"     ,"h1_W"     ,100,   0,  4);	PrettyH1F(h1_W     ,"W [GeV]"              ,"Counts",4);
		H1F h1_xB    = new H1F("h1_xB"    ,"h1_xB"    ,100,   0,  4);	PrettyH1F(h1_xB    ,"x_B"                  ,"Counts",4);

		H1F h1_dlt_vz= new H1F("h1_dlt_vz","h1_dlt_vz",100, -20, 20);	PrettyH1F(h1_dlt_vz,"#Delta v_{z} [cm]"    ,"Counts",4);

		// 2D histograms
		H2F h2_e_Ep_p_0 = new H2F("h2_e_Ep_p_0" ,"h2_e_Ep_p_0" ,100,1   ,  7,100,  0,0.4);
		H2F h2_e_Ep_p_1 = new H2F("h2_e_Ep_p_1" ,"h2_e_Ep_p_1" ,100,1   ,  7,100,  0,0.4);	
		H2F h2_e_th_phi = new H2F("h2_e_th_phi" ,"h2_e_th_phi" ,100,-190,190,100,  0, 30);
		H2F h2_p_th_phi = new H2F("h2_p_th_phi" ,"h2_p_th_phi" ,100,-190,190,100,  0, 80);
		H2F h2_beta_p_0 = new H2F("h2_beta_p_0" ,"h2_beta_p_0" ,100,0   ,  4,100,0.1,1.1);
		H2F h2_beta_p_1 = new H2F("h2_beta_p_1" ,"h2_beta_p_1" ,100,0   ,  4,100,0.1,1.1);
		H2F h2_e_vz_phi = new H2F("h2_e_vz_phi" ,"h2_e_vz_phi" ,100,-190,190,100,-50, 50);
		H2F h2_p_vz_phi = new H2F("h2_p_vz_phi" ,"h2_p_vz_phi" ,100,-190,190,100,-50, 50);
		
		PrettyH2F(h2_e_Ep_p_0,"p_{e}"       ,"E_{e}/p_{e}"   );
		PrettyH2F(h2_e_Ep_p_1,"p_{e}"       ,"E_{e}/p_{e}"   );
		PrettyH2F(h2_e_th_phi,"#phi_e [deg]","#theta_e [deg]");
		PrettyH2F(h2_p_th_phi,"#phi_p [deg]","#theta_p [deg]");
		PrettyH2F(h2_beta_p_0,"p [GeV]"     ,"#beta"         );
		PrettyH2F(h2_beta_p_1,"p [GeV]"     ,"#beta"         );
		PrettyH2F(h2_e_vz_phi,"#phi_e [deg]","e v_{z} [cm]"  );
		PrettyH2F(h2_p_vz_phi,"#phi_p [deg]","p v_{z} [cm]"  );
				
		// ----------------------------------------------------------------------------------
		// Opening input HIPO file
		HipoReader reader = new HipoReader();
		
		// RGA example
		String dataFile = "/Users/Reynier/WORK/CLAS12/data/cooked/rga/out_clas_003842.evio.295.hipo"; // target lH2
		double Ebeam = 6.4; //GeV
		
		// RGB example
		//String dataFile = "/Users/Reynier/WORK/CLAS12/data/cooked/rgb/out_clas_006164.evio.00266.hipo"; // target deuterium
		//double Ebeam = 10.6; //GeV
		
		
		reader.open(dataFile);

		// ----------------------------------------------------------------------------------
		// Creating output HIPO file
		HipoWriter writer = reader.createWriter();
		writer.open("/Users/Reynier/WORK/CLAS12/data/cooked/rga/filter_output.hipo");

		int event_counter = 0;
		// ----------------------------------------------------------------------------------
		// Loop over events and print them on the screen
		while(reader.hasNext()==true){
			HipoEvent event = reader.readNextEvent();
			//event.show();

			if(event_counter%10000==0) System.out.println("event: "+event_counter);
			event_counter++;

			if(     (event.hasGroup("REC::Particle"))&&
					(event.hasGroup("REC::Calorimeter"))
					){	

				HipoGroup bank_particle    = event.getGroup("REC::Particle"   );
				HipoGroup bank_calorimeter = event.getGroup("REC::Calorimeter");

				//HipoGroup rEC_event = event.getGroup("REC::Event");
				//HipoGroup bank_FTOF = event.getGroup("FTOF::hbhits"); 

				//bank_particle.show();
				//bank_calorimeter.show();

				// -------------------------------------------------------------------------
				// Electron PID from Dan Carman
				// - (DONE)	pid=11 from EB
				// - (DONE)	p > 2 GeV
				// - (DONE)	p < Ebeam
				// - 		TOF > 10 ns //(need bank 330) event.json
				// - (DONE)	vz: -15 to 10 cm
				// - (DONE)	W > 0 GeV
				// - (DONE-ish)		Sampling fraction +/-3sigma about E/p vs. p mean
				// - (DONE)	~15 cm fiducial cuts on U, V, W to contain full shower (need bank 332 lu, lv, lw)
				// - (DONE)	abs(chisq PID) < 5 (goodness of PID from EB)
				// - (DONE)	PCAL > 60 MeV (to remove min-i) (bank 332 layers 1(PCAL), 4(EC inner), 7(EC outter))
				// -------------------------------------------------------------------------

				// first particle variables
				int pid0       = bank_particle.getNode("pid"    ).getInt  (0);		// electron candidate id assigned by clas
				int chr0       = bank_particle.getNode("charge" ).getInt  (0);		// electron candidate charge
				double chi2pid = bank_particle.getNode("chi2pid").getFloat(0);		// electron candidate goodness of PID from EB
				double epx     = bank_particle.getNode("px"     ).getFloat(0);		// electron candidate momentum x-component [GeV]
				double epy     = bank_particle.getNode("py"     ).getFloat(0);		// electron candidate momentum y-component [GeV]
				double epz     = bank_particle.getNode("pz"     ).getFloat(0);		// electron candidate momentum z-component [GeV]
				double evx     = bank_particle.getNode("vx"     ).getFloat(0);		// electron candidate vertex x coordinate [cm]
				double evy     = bank_particle.getNode("vy"     ).getFloat(0);		// electron candidate vertex y coordinate [cm]
				double evz     = bank_particle.getNode("vz"     ).getFloat(0);		// electron candidate vertex z coordinate [cm]

				short cal_id = 0;
				double Epcal = 0;
				double Eecin = 0;
				double Eecou = 0;

				int nCal = bank_calorimeter.getNode("pindex").getDataSize();
				for(int cal=0 ; cal < nCal ; cal++) {
					cal_id = bank_calorimeter.getNode("pindex").getShort(cal);
					int layer = bank_calorimeter.getNode("layer").getByte(cal);
					if     (cal_id==0&&layer==1) Epcal = bank_calorimeter.getNode("energy").getFloat(cal);	// PCal
					else if(cal_id==0&&layer==4) Eecin = bank_calorimeter.getNode("energy").getFloat(cal);	// ECinner
					else if(cal_id==0&&layer==7) Eecou = bank_calorimeter.getNode("energy").getFloat(cal);	// ECouter
				}

				double Ee      = Epcal + Eecin + Eecou;	// electron candidate energy from calorimeter [GeV]
				double lU      = bank_calorimeter.getNode("lu"    ).getFloat(0);	// electron candidate distance on U-side [cm?]
				double lV      = bank_calorimeter.getNode("lv"    ).getFloat(0);	// electron candidate distance on V-side [cm?]
				double lW      = bank_calorimeter.getNode("lw"    ).getFloat(0);	// electron candidate distance on W-side [cm?]

				// calculated variables
				double ep     = Math.sqrt(epx*epx + epy*epy + epz*epz);		// electron candidate momentum magnitude [GeV]
				Vector3 v3_ep = new Vector3(epx,epy,epz);					// electron candidate momentum vector [GeV]
				double th_e   = v3_ep.theta();//Math.acos(epz/ep);			// electron candidate theta [rad]
				double phi_e  = v3_ep.phi();								// electron candidate phi [rad]
				double Q2     = 4*ep*Ebeam*Math.pow(Math.sin(th_e/2.),2);	// Q-squared [GeV^2]
				double nu     = Ebeam - ep;									// Transfer energy [GeV]
				double W2     = mtar*mtar-Q2+2*nu*mtar;						// Invariant mass ^2 [GeV]
				double xB     = Q2/2./mp/nu;								// Bjorken-x

				// Transfer variables
				double qx = - epx;
				double qy = - epy;
				double qz = Ebeam - epz;

				// -------------------------------------------------------------------------
				// Fill some histograms before cutting on good electrons
				h2_e_Ep_p_0.fill(ep,Ee/ep);
				h2_e_vz_phi.fill(rad2deg*phi_e, evz);
				
				// -------------------------------------------------------------------------
				// Only keep events for which the first particle is an electron
				if(     (pid0!=11            )||
						(chr0!=-1            )||
						(chi2pid>=cut_chi2pid)||
						(ep<=cut_ep          )||
						(ep>=Ebeam           )||
						(evz>cut_max_vz      )||
						(evz<cut_min_vz      )||
						(lU<cut_uvw          )||
						(lV<cut_uvw          )||
						(lW<cut_uvw          )||
						(Epcal<cut_Epcal     )||
						(Math.sqrt(W2)<=cut_W)
						) continue;

				//if(Ee/ep>0.30||Ee/ep<0.20) continue;

				// QE cut
				//if((Math.sqrt(W2)>1.097+0.1)||(Math.sqrt(W2)<1.097-0.1)) continue;
				//if((xB>1+0.2)||(xB<1-0.2)) continue;
				if((xB>1+0.5)||(xB<1-0.5)) continue;

				// Filling electron histograms
				h1_e_lu .fill(lU           );
				h1_e_lv .fill(lV           );
				h1_e_lw .fill(lW           );
				h1_e_px .fill(epx          );
				h1_e_py .fill(epy          );
				h1_e_pz .fill(epz          );
				h1_e_p  .fill(ep           );
				h1_e_vz .fill(evz          );
				h1_W    .fill(Math.sqrt(W2));
				h1_xB   .fill(xB           );
				h1_e_th .fill(rad2deg*th_e );
				h1_e_phi.fill(rad2deg*phi_e);

				h2_e_Ep_p_1.fill(ep,Ee/ep);
				h2_e_th_phi.fill(rad2deg*phi_e,rad2deg*th_e);

				// Saving these events in output file
				writer.writeEvent(event);
				
				// -------------------------------------------------------------------------
				// Loop over the remaining particles and require a proton
				int nProtons = 0;
				int tmp_fast_p_idx  = 0;	// index of the fastest proton
				double tmp_fast_p_p = 0;	// momentum of the fastest proton

				int nParticles = bank_particle.getNode("pid").getDataSize();
				for(int par = 1; par < nParticles; par++) {
					int pidi      = bank_particle.getNode("pid"   ).getInt  (par);		// proton candidate id assigned by clas
					int chri      = bank_particle.getNode("charge").getInt  (par);		// proton candidate charge
					double beta_p = bank_particle.getNode("beta"  ).getFloat(par);		// proton candidate beta (v/c)
					double ppx    = bank_particle.getNode("px"    ).getFloat(par);		// proton candidate momentum x-component [GeV]
					double ppy    = bank_particle.getNode("py"    ).getFloat(par);		// proton candidate momentum y-component [GeV]
					double ppz    = bank_particle.getNode("pz"    ).getFloat(par);		// proton candidate momentum z-component [GeV]
					double pvx    = bank_particle.getNode("vx"    ).getFloat(par);		// proton candidate vertex x coordinate [cm]
					double pvy    = bank_particle.getNode("vy"    ).getFloat(par);		// proton candidate vertex y coordinate [cm]
					double pvz    = bank_particle.getNode("vz"    ).getFloat(par);		// proton candidate vertex z coordinate [cm]

					double pp     = Math.sqrt(ppx*ppx + ppy*ppy + ppz*ppz);		

					// Fill some histograms before checking the particles are protons
					h2_beta_p_0.fill(pp, beta_p);

					// ----------------------------------------------------------------------
					// If there are protons, find the fastest one and consider that (e,e'p)X
					// (Add proton event selection cuts below)
					if(     (pidi==2212)&&
							(chri==1   )
							) {
						nProtons++;
						if(pp>tmp_fast_p_p) {
							tmp_fast_p_p = pp;
							tmp_fast_p_idx=par;
						}
					}
					// ----------------------------------------------------------------------
				}
				if(nProtons>0) {

					int pidi      = bank_particle.getNode("pid"   ).getInt  (tmp_fast_p_idx);		// proton candidate id assigned by clas
					int chri      = bank_particle.getNode("charge").getInt  (tmp_fast_p_idx);		// proton candidate charge
					double beta_p = bank_particle.getNode("beta"  ).getFloat(tmp_fast_p_idx);		// proton candidate beta (v/c)
					double ppx    = bank_particle.getNode("px"    ).getFloat(tmp_fast_p_idx);		// proton candidate momentum x-component [GeV]
					double ppy    = bank_particle.getNode("py"    ).getFloat(tmp_fast_p_idx);		// proton candidate momentum y-component [GeV]
					double ppz    = bank_particle.getNode("pz"    ).getFloat(tmp_fast_p_idx);		// proton candidate momentum z-component [GeV]
					double pvx    = bank_particle.getNode("vx"    ).getFloat(tmp_fast_p_idx);		// proton candidate vertex x coordinate [cm]
					double pvy    = bank_particle.getNode("vy"    ).getFloat(tmp_fast_p_idx);		// proton candidate vertex y coordinate [cm]
					double pvz    = bank_particle.getNode("vz"    ).getFloat(tmp_fast_p_idx);		// proton candidate vertex z coordinate [cm]

					double pp     = Math.sqrt(ppx*ppx + ppy*ppy + ppz*ppz);
					double Ep     = Math.sqrt(pp*pp+mp*mp);
					Vector3 v3_pp = new Vector3(ppx,ppy,ppz);							// proton candidate momentum vector [GeV]		
					double th_p   = v3_pp.theta();										// proton candidate theta [rad]
					double phi_p  = v3_pp.phi();										// proton candidate phi [rad]

					double delt_vz= pvz - evz;
					
					// Missing momentum components
					double pmx = ppx - qx;
					double pmy = ppy - qy;
					double pmz = ppz - qz;
					double Pm = Math.sqrt(pmx*pmx + pmy*pmy + pmz*pmz);

					// Missing mass
					double E_mmiss = Ebeam + mtar - ep - Ep;
					double Mmiss = Math.sqrt(E_mmiss*E_mmiss - Pm*Pm);
					
					h1_p_vz    .fill(pvz          );
					h1_dlt_vz  .fill(delt_vz      );
					h1_p_px    .fill(ppx          );
					h1_p_py    .fill(ppy          );
					h1_p_pz    .fill(ppz          );
					h1_p_p     .fill(pp           );
					h1_p_th    .fill(rad2deg*th_p );
					h1_p_phi   .fill(rad2deg*phi_p);
					h1_pmx     .fill(pmx          );
					h1_pmy     .fill(pmy          );
					h1_pmz     .fill(pmz          );
					h1_pmiss   .fill(Pm           );
					h1_Mmiss   .fill(Mmiss        );
					
					h2_p_th_phi.fill(rad2deg*phi_p, rad2deg*th_p);
					h2_p_vz_phi.fill(rad2deg*phi_p, pvz         );
					h2_beta_p_1.fill(pp           , beta_p      );
					
				}

				h1_p_num.fill((double)(nProtons));

			}

		} // End of while loop over events

		writer.close();

		// -------------------------------------------------------------------------------------------
		// Drawing histograms
		TCanvas c0 = new TCanvas("c0", 800, 600);
		c0.divide(2,1);
		c0.cd(0);
		c0.draw(h1_e_vz);
		c0.draw(h1_p_vz,"same");
		c0.cd(1);
		c0.draw(h1_dlt_vz);

		TCanvas c1 = new TCanvas("c1", 800, 600);
		c1.divide(2, 2);
		c1.cd(0);	c1.draw(h1_pmx);
		c1.cd(1);	c1.draw(h1_pmy);
		c1.cd(2);	c1.draw(h1_pmz);
		c1.cd(3);	c1.draw(h1_pmiss);

		TCanvas c2 = new TCanvas("c2", 800, 600);
		c2.divide(2, 2);
		c2.cd(0);	c2.draw(h1_e_th);
		c2.cd(1);	c2.draw(h2_e_th_phi);
		c2.cd(3);	c2.draw(h1_e_phi);

		TCanvas c3 = new TCanvas("c3", 800, 600);
		c3.draw(h1_W);

		TCanvas c4 = new TCanvas("c4", 800, 600);
		c4.draw(h1_xB);

		TCanvas c5 = new TCanvas("c5", 800, 600);
		c5.divide(2, 1);
		c5.cd(0);	c5.draw(h2_e_Ep_p_0);
		c5.cd(1);	c5.draw(h2_e_Ep_p_1);

		TCanvas c6 = new TCanvas("c6", 800, 600);
		c6.divide(2, 2);
		c6.cd(0);	c6.draw(h1_e_px);
		c6.cd(1);	c6.draw(h1_e_py);
		c6.cd(2);	c6.draw(h1_e_pz);
		c6.cd(3);	c6.draw(h1_e_p );

		TCanvas c7 = new TCanvas("c7", 800, 600);
		c7.divide(2, 2);
		c7.cd(0);	c7.draw(h1_p_px);
		c7.cd(1);	c7.draw(h1_p_py);
		c7.cd(2);	c7.draw(h1_p_pz);
		c7.cd(3);	c7.draw(h1_p_p );

		TCanvas c8 = new TCanvas("c8", 800, 600);
		c8.divide(2, 1);
		c8.cd(0);	c8.draw(h2_beta_p_0);
		c8.cd(1);	c8.draw(h2_beta_p_1);

		TCanvas c9 = new TCanvas("c9", 800, 600);
		c9.divide(2, 2);
		c9.cd(0);	c9.draw(h1_p_th);
		c9.cd(1);	c9.draw(h2_p_th_phi);
		c9.cd(3);	c9.draw(h1_p_phi);

		TCanvas c10 = new TCanvas("c10", 800, 600);
		c10.divide(2, 2);
		c10.cd(0);	c10.draw(h1_e_lu);
		c10.cd(1);	c10.draw(h1_e_lv);
		c10.cd(2);	c10.draw(h1_e_lw);

		TCanvas c11 = new TCanvas("c11", 800, 600);
		c11.draw(h1_p_num);
		
		TCanvas c12 = new TCanvas("c12", 800, 600);
		c12.draw(h1_Mmiss);
		
		TCanvas c13 = new TCanvas("c13", 800, 600);
		c13.divide(2, 1);
		c13.cd(0);	c13.draw(h2_e_vz_phi);
		c13.cd(1);	c13.draw(h2_p_vz_phi);

	}
	// =========================================================================================================
	public static void PrettyH1F(H1F h1,String titx,String tity,int color) {
		h1.setTitleX(titx);
		h1.setTitleY(tity);
		h1.setLineColor(color);
		h1.setLineWidth(3);
	}
	// =========================================================================================================
	public static void PrettyH2F(H2F h2,String titx,String tity) {
		h2.setTitleX(titx);
		h2.setTitleY(tity);
	}
}
