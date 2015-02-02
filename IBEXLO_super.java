
/**
* Simulating response at IBEX_LO with this class
*
*  Create interstellar medium and perform integration
*
*  This one assumes spacecraft can look all around the world and makes 2D maps for each spacecraft location.
*/
public class IBEXLO_super {

	public int mcN = 20000;  // number of iterations per 3D integral
	public String outFileName = "ibex_lo_out17.txt";
	String dir = "SPRING";
	//String dir = "FALL";

	public static final double EV = 1.60218 * Math.pow(10, -19);
	public static double MP = 1.672621*Math.pow(10.0,-27.0);
	public static double AU = 1.49598* Math.pow(10,11); //meters
	public static double EARTHSPEED = AU*2.0*3.14/365.0/24.0/3600.0;

	// define angles and energies of the instrument here

	// H mode (standard) energy bins in eV
	public double[] minEnergies = {8.8, 17.6, 35.8, 74.2, 153.5, 330.6, 663.3, 1298.6};
	public double[] maxEnergies = {23.2, 46.4, 94.3, 195.8, 404.6, 871.5, 1748.7, 3423.5};
	public double[] minSpeedsHe, maxSpeedsHe, minSpeedsO, maxSpeedsO;

	// Interstellar mode energies
	// { spring min, spring max, fall min, fall max }
	public double[] heEnergies = { 73.7, 194.3, 4.4, 11.6 };
	public double[] oEnergies = {293.7, 774.3, 18.15, 47.9 };

	public double[] heSpeeds;
	public double[] oSpeeds;

	// O mode efficiencies
	public double heSpringEff = 1.0*Math.pow(10,-7);
	public double heFallEff = 1.6*Math.pow(10,-6);
	public double oSpringEff = 0.004;
	public double oFallEff = 0.001;

	// angular acceptance - assume rectangular windows..
	// these are half widths, e.g. +- for acceptance
	public double spinWidth = 4.0*Math.PI/180.0;
	public double hrWidth = 3.5*Math.PI/180.0;
	public double lrWidth = 7.0*Math.PI/180.0;

	public double eightDays = 8.0*24.0*60.0*60.0;
	public double oneHour = 3600.0;

	public double upWind = 74.0 * 3.14 / 180.0;
	public double downWind = 254.0 * 3.14 / 180.0;

	public double startPerp = 180*3.14/180;  // SPRING
	public double endPerp = 240*3.14/180;    // SPRING

	public double he_ion_rate = 1.0*Math.pow(10,-7);
	public double o_ion_rate = 8.0*Math.pow(10,-7);


	public double activeArea = 100.0/100.0/100.0;  // in square meters now!!
	public file outF;
	public MultipleIntegration mi;


	/**
	*
	*
	*
	*/
	public IBEXLO_super() {

		if (dir.equals("SPRING")) {
			startPerp = 175*3.14/180;
			endPerp = 235*3.14/180;
		}

		else if (dir.equals("FALL")) {
			startPerp = 290.0*3.14/180.0;  // FALL
			endPerp = 360.0*3.14/180.0;    // FALL
		}

		o("earthspeed: " + EARTHSPEED);
		mi = new MultipleIntegration();

		// set up output file
		outF = new file(outFileName);
		outF.initWrite(false);

		// calculate speed ranges with v = sqrt(2E/m)
		heSpeeds = new double[heEnergies.length];
		oSpeeds = new double[oEnergies.length];
		o("calculating speed range");

		for (int i=0; i<4; i++) {
			heSpeeds[i]=Math.sqrt(2.0*heEnergies[i]/4.0/MP*EV);
			oSpeeds[i]=Math.sqrt(2.0*oEnergies[i]/16.0/MP*EV);
			System.out.println(heEnergies[i] +" = " + heSpeeds[i] + " , " + oEnergies[i] + " = " + oSpeeds[i]);
		}

		minSpeedsHe = new double[8];
		maxSpeedsHe = new double[8];
		minSpeedsO = new double[8];
		maxSpeedsO = new double[8];

		// still setting up speeds to integrate over..
		for (int i=0; i<8; i++) {
			minSpeedsHe[i]=Math.sqrt(2.0*minEnergies[i]/4.0/MP*EV);
			maxSpeedsHe[i]=Math.sqrt(2.0*maxEnergies[i]/4.0/MP*EV);
			minSpeedsO[i]=Math.sqrt(2.0*minEnergies[i]/16.0/MP*EV);
			maxSpeedsO[i]=Math.sqrt(2.0*maxEnergies[i]/16.0/MP*EV);
			System.out.println(minSpeedsHe[i] +" = " + maxSpeedsHe[i] + " , " + minSpeedsO[i] + " = " + maxSpeedsO[i]);
		}

		System.out.println(minSpeedsHe[0] + "\t" + maxSpeedsHe[0] +"\t" + heSpeeds[2] + "\t" + heSpeeds[3]);
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

		GaussianVLISM gv1 = new GaussianVLISM(bulkHE,0.015*100*100*100,6300.0);
		GaussianOxVLISM gv2 = new GaussianOxVLISM(bulkO1, bulkO2, 0.00005*100*100*100, 1000.0, 90000.0, 0.0); // last is fraction in 2ndry

		final NeutralDistribution ndHe = new NeutralDistribution(gv1, 0.0, he_ion_rate);
		final NeutralDistribution ndO = new NeutralDistribution(gv2, 0.0, o_ion_rate);

		ndHe.debug=false;
		ndO.debug=false;
		// DONE CREATING MEDIUM

		//  PASSAGE SIMULATION - IN ECLIPTIC - SPRING or FALL  -  use String dir = "FALL" or "SPRING" to set all pointing params properly
		double initPos = startPerp;
		double finalPos = endPerp;
		double step = 8.0*360.0/365.0*3.14/180.0;  // these are orbit adjusts - every eight days
		//double step2 = step/8.0/2.0;               // here we take smaller integration timest than an orbit.. current: 12 hours
		// sweep this now around 360

		// determine entrance velocity vector and integration range
		double midPhi = initPos-Math.PI/2.0;
		double midThe = Math.PI/2.0;

		// move spacecraft
		for (double pos = initPos; pos<finalPos; pos+= step) {

			// inner loop moves pos2 - pos is what gives us pointing, not position!
			for (double pos2 = pos; pos2<pos+step; pos2+= step2) {

				// here's the now position
				final HelioVector posVec   = new HelioVector(HelioVector.SPHERICAL, 1*AU, pos2, Math.PI/2.0);
				final HelioVector scVelVec = new HelioVector(HelioVector.SPHERICAL, EARTHSPEED, pos2+Math.PI/2.0, Math.PI/2.0);

				FunctionIII he_func = new FunctionIII() {  // helium neutral distribution
					public double function3(double v, double p, double t) {
						double tbr = activeArea*ndHe.dist(posVec, new HelioVector(HelioVector.SPHERICAL,v,p,t).difference(scVelVec) )*v*v*Math.sin(t)*v;  // extra v for flux
						// that is the integrand of our v-space integration
						return tbr;
					}
				};

				FunctionIII o_func = new FunctionIII() {   // oxygen neutral distribution
					public double function3(double v, double p, double t) {
						double tbr = activeArea*ndO.dist(posVec, new HelioVector(HelioVector.SPHERICAL,v,p,t).difference(scVelVec) )*v*v*Math.sin(t)*v;  // extra v for flux
						// that is the integrand of our v-space integration
						return tbr;
					}
				};

				// set look direction and limits -
				if (dir.equals("FALL"))
					midPhi = pos-Math.PI/2+(4.0*3.1416/180.0); // change look direction at downwind point

				else if (dir.equals("SPRING"))
					midPhi = pos+Math.PI/2+(4.0*3.1416/180.0); // change look direction at downwind point

				// now step through the spin - out of ecliptic
				double theta_step = 6.0*Math.PI/180.0;

				for (double pos3 = Math.PI/2.0 - 3.0*theta_step; pos3<= (Math.PI/2.0 + 4.0*theta_step); pos3 += theta_step) {
					midThe = pos3;

					// set limits for integration, fall is first here..
					//if (dir.equals("FALL")) {
					double[] he_hr_lims = { heSpeeds[2], heSpeeds[3], midPhi-lrWidth, midPhi+lrWidth, midThe-lrWidth-spinWidth, midThe+lrWidth+spinWidth };
					double[] he_lr_lims = { minSpeedsHe[0], maxSpeedsHe[0], midPhi-lrWidth, midPhi+lrWidth, midThe-lrWidth-spinWidth, midThe+lrWidth+spinWidth };
					double[] o_hr_lims = { oSpeeds[2], oSpeeds[3], midPhi-hrWidth, midPhi+hrWidth, midThe-hrWidth-spinWidth, midThe+hrWidth+spinWidth };
					double[] o_lr_lims = { minSpeedsO[1], maxSpeedsO[1], midPhi-lrWidth, midPhi+lrWidth, midThe-lrWidth-spinWidth, midThe+lrWidth+spinWidth };
					//}

					if (dir.equals("SPRING")) {
						he_hr_lims[0]=heSpeeds[0];  he_hr_lims[1]=heSpeeds[1]; //, midPhi-hrWidth, midPhi+hrWidth, midThe-hrWidth-spinWidth, midThe+hrWidth+spinWidth };
						he_lr_lims[0]=minSpeedsHe[3]; he_lr_lims[1]=maxSpeedsHe[3]; //, midPhi-lrWidth, midPhi+lrWidth, midThe-lrWidth-spinWidth, midThe+lrWidth+spinWidth };
						o_hr_lims[0]=oSpeeds[0]; o_hr_lims[1]=oSpeeds[1]; //, midPhi-hrWidth, midPhi+hrWidth, midThe-hrWidth-spinWidth, midThe+hrWidth+spinWidth };
						o_lr_lims[0]=minSpeedsO[5]; o_lr_lims[1]=maxSpeedsO[5]; //, midPhi-lrWidth, midPhi+lrWidth, midThe-lrWidth-spinWidth, midThe+lrWidth+spinWidth };
					}

					// time * spin_duty * energy_duty * area_factor
					double o_ans_hr =  12.0*oneHour *0.017 *7.0/8.0 *0.25 *   mi.mcIntegrate(o_func, o_hr_lims, mcN);
					double o_ans_lr =  12.0*oneHour *0.017 *1.0/8.0 *1.0  *   mi.mcIntegrate(o_func, o_lr_lims, mcN);
					//double he_ans_hr = 12.0*oneHour *0.017 *1.0/8.0 *1.0 *   mi.mcIntegrate(he_func, he_hr_lims, mcN);
					double he_ans_lr = 12.0*oneHour *0.017 *1.0/8.0 *1.0  *   mi.mcIntegrate(he_func, he_lr_lims, mcN);

					if (dir.equals("SPRING")) {
						o_ans_hr*=oSpringEff;
						o_ans_lr*=oSpringEff;
						//he_ans_hr*=heSpringEff;
						he_ans_lr*=heSpringEff;
					}
					else if (dir.equals("FALL")) {
						o_ans_hr*=oFallEff;
						o_ans_lr*=oFallEff;
						//he_ans_hr*=heFallEff;
						he_ans_lr*=heFallEff;
					}

				//	double doy = (pos2*180.0/3.14)*365.0/360.0;
					outF.write(pos2*180.0/3.14 + "\t" + pos3 + "\t" + o_ans_hr  + "\t" + o_ans_lr +"\t" + he_ans_lr + "\n");
					o("pos: " + pos2 +"\t"+ "the: " + pos3*180.0/3.14 + "\t" + o_ans_hr  + "\t" + o_ans_lr + "\t" + he_ans_lr);
				}
			}
		}
		outF.closeWrite();



			//double[] limits = {heSpeeds[0], heSpeeds[1], midPhi-hrWidth, midPhi+hrWidth, midThe-hrWidth-spinWidth, midThe+hrWidth+spinWidth};
			//  THIS TEST is moving the look angle with the Earth at the supposed proper position.
			//double step = Math.abs((finalPos-initPos)/11.0);
	/*		initPos = 0.0;
			finalPos = 2.0*Math.PI;
			step = Math.abs((finalPos-initPos)/60.0);
			for (double pos = initPos; pos<finalPos; pos+= step) {
				// here's the now position
				final HelioVector posVec = new HelioVector(HelioVector.SPHERICAL, 100*AU, 0.0*Math.PI, Math.PI/2.0);
				midPhi = pos;
				midThe = Math.PI/2.0;
				double[] limits = { 30000.0, 300000.0, midPhi-hrWidth, midPhi+hrWidth, midThe-hrWidth-spinWidth, midThe+hrWidth+spinWidth };

				o("doing pos: " + pos);
				FunctionIII f3 = new FunctionIII() {
					public double function3(double v, double p, double t) {
						double tbr = ndH.dist(posVec, new HelioVector(HelioVector.SPHERICAL,v,p,t))*v*v*Math.sin(t);  // extra v for flux
						// that is the integrand of our v-space integration
						return tbr;
					}
				};
				outF.write(pos + "\t" + mi.mcIntegrate(f3, limits, 10000) + "\n");

			}
			outF.closeWrite();
	*/
			// 8 days, break to 12 hr

	}

	public static final void main(String[] args) {
		IBEXLO_super il = new IBEXLO_super();
	}

	public static void o(String s) {
		System.out.println(s);
	}
}

