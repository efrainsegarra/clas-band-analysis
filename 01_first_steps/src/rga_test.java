import java.io.Reader;

import org.jlab.groot.data.H1F;
import org.jlab.groot.data.H2F;
import org.jlab.groot.ui.TCanvas;

import org.jlab.jnp.hipo.data.HipoEvent;
import org.jlab.jnp.hipo.data.HipoGroup;
import org.jlab.jnp.hipo.data.HipoNode;
import org.jlab.jnp.hipo.io.HipoReader;

import org.jlab.clas.physics.*;

import org.jlab.io.evio.*;

import org.jlab.rec.band.constants.Parameters;
import org.jlab.rec.band.hit.BandHit;
import org.omg.CORBA.PRIVATE_MEMBER;

import com.sun.javafx.runtime.SystemProperties;

import jdk.internal.dynalink.beans.StaticClass;
import sun.security.krb5.internal.crypto.crc32;

import java.util.ArrayList;
import java.util.Scanner;

public class code {
	
	// =========================================================================================================
	public static void main(String[] args) {
		
		double Ebeam = 6.4;//GeV
		double mp = 0.93827; //GeV
		double mtar = mp;
		
		int ctr=0;
		double rad2deg = 180./3.14159;
		
		// ----------------------------------------------------------------------------------
		// Declaring histograms
		H1F h1_pmiss = new H1F("h1_pmiss","h1_pmiss",100, 0,5);
		H1F h1_pmx   = new H1F("h1_pmx"  ,"h1_pmx"  ,100,-2,2);
		H1F h1_pmy   = new H1F("h1_pmy"  ,"h1_pmy"  ,100,-2,2);
		H1F h1_pmz   = new H1F("h1_pmz"  ,"h1_pmz"  ,100,-3,3);
		
		H1F h1_theta = new H1F("h1_theta","h1_theta",100,0,90);
		H1F h1_W2    = new H1F("h1_W2"   ,"h1_W2"   ,100,0,4);
		H1F h1_xB    = new H1F("h1_xB"   ,"h1_xB"   ,100,0,4);
		// ----------------------------------------------------------------------------------
		// Opening HIPO file
		HipoReader reader = new HipoReader();
		String dataFile = "../../../data/cooked/rga/out_clas_003842.evio.295.hipo"; // target lH2
		reader.open(dataFile);
		 GenericKinematicFitter fitter = new GenericKinematicFitter(Ebeam);
	
		// ----------------------------------------------------------------------------------
		// Loop over events and print them on the screen
		while(reader.hasNext()==true){
			HipoEvent event = reader.readNextEvent();
			
			/*
			EventFilter filter = new EventFilter("11:2212"); // exclusive e-,p,p
			
			PhysicsEvent recEvent = fitter.getPhysicsEvent(event);
			
			if(filter.isValid(event)==true){
				
			}
			*/
			
			//PhysicsEvent recEvent = fitter.getPhysicsEvent(event);
				
			//PhysicsEvent recEvent = fitter.getPhysicsEvent(event);
			//event.show();
			
			/*
			if(event.hasGroup("RECHB::Track")){	
				HipoGroup bank_temp = event.getGroup("RECHB::Track");
				bank_temp.show();
			}
			*/
			
			if(event.hasGroup("RECHB::Particle")){	
				HipoGroup bank_particle = event.getGroup("RECHB::Particle");

				//bank_particle.show();
				
				// Only keep events for which the first particle is an electron
				int pid0 = bank_particle.getNode("pid").getInt(0);
				int chr0 = bank_particle.getNode("charge").getInt(0);
				if(     (pid0!=11)||
						(chr0!=-1)
						) continue;
			
				// Loop over the remaining particles and require a proton
				int nParticles = bank_particle.getNode("pid").getDataSize();
				for(int par = 1; par < nParticles; par++) {
					int pidi = bank_particle.getNode("pid").getInt(par);
					int chri = bank_particle.getNode("charge").getInt(par);
					if(     (pidi==2212)&&
							(chri==1)
							) {
						double epx = bank_particle.getNode("px").getFloat(0);
						double epy = bank_particle.getNode("py").getFloat(0);
						double epz = bank_particle.getNode("pz").getFloat(0);
						double ep = Math.sqrt(epx*epx + epy*epy + epz*epz);
						
						double th_e = Math.acos(epz/ep);
						double Q2 = 4*ep*Ebeam*Math.pow(Math.sin(th_e/2.),2);
						double nu = Ebeam - ep;
						
						double W2 = mtar*mtar-Q2+2*nu*mtar;
						h1_W2.fill(Math.sqrt(W2));
						
						double xB = Q2/2./mp/nu;
						h1_xB.fill(xB);
						
						double qx = - epx;
						double qy = - epy;
						double qz = Ebeam - epz;
						
						double ppx = bank_particle.getNode("px").getFloat(par);
						double ppy = bank_particle.getNode("py").getFloat(par);
						double ppz = bank_particle.getNode("pz").getFloat(par);
						
						double pmx = ppx - qx;
						double pmy = ppy - qy;
						double pmz = ppz - qz;
						double Pm = Math.sqrt(pmx*pmx + pmy*pmy + pmz*pmz);
						
						h1_theta.fill(rad2deg*th_e);
						
						h1_pmx.fill(pmx);
						h1_pmy.fill(pmy);
						h1_pmz.fill(pmz);
						h1_pmiss.fill(Pm);
					}
				}
				//	bank_particle.show();
				//}
			}
			
		}
		// Drawing histograms
		TCanvas c1 = new TCanvas("c1", 800, 600);
		c1.divide(2, 2);
		c1.cd(0);
		PrettyHistogram(h1_pmx,"Pmx [GeV]");
		c1.draw(h1_pmx);
		c1.cd(1);
		PrettyHistogram(h1_pmy,"Pmy [GeV]");
		c1.draw(h1_pmy);
		c1.cd(2);
		PrettyHistogram(h1_pmz,"Pmz [GeV]");
		c1.draw(h1_pmz);
		c1.cd(3);
		PrettyHistogram(h1_pmiss,"Pm [GeV]");
		c1.draw(h1_pmiss);
	
		TCanvas c2 = new TCanvas("c2", 800, 600);
		PrettyHistogram(h1_theta,"#theta_e [deg]");
		c2.draw(h1_theta);
		
		TCanvas c3 = new TCanvas("c3", 800, 600);
		PrettyHistogram(h1_W2,"W [GeV]");
		c3.draw(h1_W2);
		
		TCanvas c4 = new TCanvas("c4", 800, 600);
		PrettyHistogram(h1_xB,"x_B");
		c4.draw(h1_xB);
		
	}
	// =========================================================================================================
	public static void PrettyHistogram(H1F h1,String titx) {
		h1.setTitleX(titx);
		h1.setLineColor(2);
		h1.setLineWidth(2);
	}
}
