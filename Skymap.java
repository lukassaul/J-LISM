
/**
* Simulating response at hypothetical imaging platform
*
*  Create interstellar medium and perform integration
*
*/
public class Skymap {

	public int mcN = 8000;  // number of iterations per 3D integral
	public double smartMin = 1.0;	// lowest flux to continue integration of sector
	public String outFileName = "skymap_out4.txt";
	public String filePrefix = "sky_add_2fine_";
	String dir = "SPRING";
	//String dir = "FALL";

	public static final double EV = 1.60218 * Math.pow(10, -19);
	public static double MP = 1.672621*Math.pow(10.0,-27.0);
	public static double AU = 1.49598* Math.pow(10,11); //meters
	public static double EARTHSPEED = AU*2.0*3.14/365.0/24.0/3600.0;



	public double he_ion_rate = 1.0*Math.pow(10,-7);
	public double o_ion_rate = 8.0*Math.pow(10,-7);


	// define angles and energies of the instrument here

	// H mode (standard) energy bins in eV


	// angular acceptance - assume rectangular windows..
	// chaged from IBEX to higher resolution "perfect" imager
	// these are half widths, e.g. +- for acceptance
	public double angleWidth = 0.5*Math.PI/180.0;

	public file outF;
	public MultipleIntegration mi;

	/**
	*
	*
	*
	*/
	public Skymap() {

		o("earthspeed: " + EARTHSPEED);
		mi = new MultipleIntegration();
		mi.smartMin=smartMin;

		// set up output file
		outF = new file(outFileName);


		// BUILD INTERSTERLLAR MEDIUM - helium and oxygen here, primary and secondary oxygen..
		//
		HelioVector bulkHE = new HelioVector(HelioVector.SPHERICAL,26000.0, (74+180)*Math.PI/180.0, 95.5*Math.PI/180);
		HelioVector bulkO1 = new HelioVector(HelioVector.SPHERICAL,26000.0, (74+180)*Math.PI/180.0, 95.5*Math.PI/180);
		HelioVector bulkO2 = new HelioVector(HelioVector.SPHERICAL,21000.0, (74+180)*Math.PI/180.0, 95.5*Math.PI/180);
		// implement test 1 - a vector coming in along x axis
		//HelioVector bulkHE = new HelioVector(HelioVector.SPHERICAL,26000.0, Math.PI, Math.PI/2.0);
		//HelioVector bulkO1 = new HelioVector(HelioVector.SPHERICAL,25000.0, Math.PI, Math.PI/2.0);
		//HelioVector bulkO2 = new HelioVector(HelioVector.SPHERICAL,21000.0, Math.PI, Math.PI/2.0);
		// standard hot helium

		GaussianVLISM gv1 = new GaussianVLISM(bulkHE,0.015*100.0*100.0*100.0,6300.0);
		GaussianOxVLISM gv2 = new GaussianOxVLISM(bulkO1, bulkO2, 0.00005*100*100*100, 1000.0, 90000.0, 0.0);
		// last is fraction in 2ndry

		final NeutralDistribution ndHe = new NeutralDistribution(gv1, 0.0, he_ion_rate);
		final NeutralDistribution ndO = new NeutralDistribution(gv2, 0.0, o_ion_rate);

		ndHe.debug=false;
		ndO.debug=false;
		// DONE CREATING MEDIUM

		final double activeArea = 0.0001; // one square cm, assuming a nice detector here, detects everything that hits
		int index = 1;
		for (double pos = 0.0*Math.PI/180.0; pos<360.0*Math.PI/180.0; pos+= 4.0*Math.PI/180.0) {
			// position of spacecraft
			final HelioVector posVec = new HelioVector(HelioVector.SPHERICAL, 1*AU, pos, Math.PI/2.0);
			final HelioVector scVelVec = new HelioVector(HelioVector.SPHERICAL, EARTHSPEED, pos+Math.PI/2.0, Math.PI/2.0);

			// helium neutral distribution
			FunctionIII he_func = new FunctionIII() {
				public double function3(double v, double p, double t) {
					double tbr = activeArea*ndHe.dist(posVec, new HelioVector(HelioVector.SPHERICAL,v,p,t).sum(scVelVec))*v*v*Math.sin(t)*v;  // extra v for flux
					// that is the integrand of our v-space integration
					return tbr;
				}
			};
		/*	FunctionIII o_func = new FunctionIII() {   // oxygen neutral distribution
				public double function3(double v, double p, double t) {
					double tbr = activeArea*ndO.dist(posVec, new HelioVector(HelioVector.SPHERICAL,v,p,t))*v*v*Math.sin(t)*v;  // extra v for flux
					// that is the integrand of our v-space integration
					return tbr;
				}
			};

		*/
			outF=new file(filePrefix+""+(1000+index)+".txt");
			outF.initWrite(false);
			// now make the image
			for (double theta=Math.PI/3.0; theta<2.0*Math.PI/3.0; theta+=1.0*Math.PI/180.0) {
				for (double phi=0.0; phi<2.0*Math.PI; phi+=1.0*Math.PI/180) {
					// set limits for integration..
					double[] lims = {0.0, 100000.0, phi-angleWidth, phi+angleWidth, theta-angleWidth, theta+angleWidth};

					// time * spin_duty * energy_duty * area_factor
					//double o_ans =  mi.mcIntegrateSS(o_func, lims, mcN);
					double he_ans = mi.mcIntegrateSS(he_func, lims, mcN);
					//System.out.println(phi+"\t"+theta+"\t"+mi.lastNP+"\t"+he_ans);
					//outF.write(pos +"\t"+theta+"\t"+phi+"\t"+he_ans+"\t"+o_ans+"\n");
					outF.write(pos*180.0/Math.PI +"\t"+theta*180.0/Math.PI +"\t"+phi*180.0/Math.PI+"\t"+he_ans+"\n");
				}
				System.out.println("Theta= " + theta);
			}
			outF.closeWrite();
			index++;
		}
	}

	public static final void main(String[] args) {
		Skymap il = new Skymap();
	}

	public static void o(String s) {
		System.out.println(s);
	}
}

