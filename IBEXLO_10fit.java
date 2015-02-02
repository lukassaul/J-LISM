import java.util.StringTokenizer;


/**
* Simulating response at IBEX_LO with this class
*
*  Create interstellar medium and perform integration
*
*/
public class IBEXLO_09fit {

	public int mcN = 100000;  // number of iterations per 3D integral
	public String outFileName = "lo_fit_09.txt";

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
//	public double heSpringEff = 1.0*Math.pow(10,-7); // this was the original..  now going to adjust to match data including sputtering etc.

	public double heSpringEff =  1.0*Math.pow(10,-7)/7.8; // nominal efficiency to "roughly" match data..


	public double heFallEff = 1.6*Math.pow(10,-6);
	public double oSpringEff = 0.004;
	public double oFallEff = 0.001;

	// angular acceptance - assume rectangular windows..
	// these are half widths, e.g. +- for acceptance
	public double spinWidth = 4.0*Math.PI/180.0;
	public double hrWidth = 3.5*Math.PI/180.0;
	public double lrWidth = 10.0*Math.PI/180.0;
	public double eightDays = 8.0*24.0*60.0*60.0;
	public double oneHour = 3600.0;
	public double upWind = 74.0 * 3.14 / 180.0;
	public double downWind = 254.0 * 3.14 / 180.0;
	public double startPerp = 180*3.14/180;  // SPRING
	public double endPerp = 240*3.14/180;    // SPRING
	//public double startPerp = 260.0*3.14/180.0;  // FALL
	//public double endPerp = 340.0*3.14/180.0;    // FALL
	public double he_ion_rate = 6.00*Math.pow(10,-8);
	public GaussianVLISM gv1;

	// temporarily make them very small
	//public double spinWidth = 0.5*Math.PI/180.0;
	//public double hrWidth = 0.5*Math.PI/180.0;
	//public double lrWidth = 0.5*Math.PI/180.0;
	public double activeArea = 100.0/100.0/100.0;  // 100 cm^2 in square meters now!!
	public file outF;
	public MultipleIntegration mi;
	public NeutralDistribution ndHe;
	public double look_offset=Math.PI;
	public double currentLongitude, currentSpeed, currentDens, currentTemp;
	public HelioVector bulkHE;


	/**
	*
	*
	*
	*/
	public IBEXLO_09fit() {

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
		bulkHE = new HelioVector(HelioVector.SPHERICAL, 26000.0, (74.68)*Math.PI/180.0, 95.5*Math.PI/180.0);
		//HelioVector bulkO1 = new HelioVector(HelioVector.SPHERICAL,26000.0, (74+180)*Math.PI/180.0, 95.5*Math.PI/180);
		GaussianVLISM gv1 = new GaussianVLISM(bulkHE,0.015*100.0*100.0*100.0,6000.0);


		//KappaVLISM gv3 = new KappaVLISM(bulkO1, 0.00005*100*100*100, 6300, 1.6);
		ndHe = new NeutralDistribution(gv1, 0.0, he_ion_rate);  // mu is 0.0 for he
		ndHe.debug=false;
		currentLongitude = 74.68;
		currentSpeed = 26000.0;
		currentDens = 0.015;
		currentTemp = 6000.0;
		// DONE CREATING MEDIUM

		o("earthspeed: " + EARTHSPEED);
		mi = new MultipleIntegration();
		o("done creating IBEXLO_09 object");

		// DONE SETUP
		//-
		//-
		// now lets run the test, see how this is doing
/*

		// moving curve fitting routine here for exhaustive and to make scan plots
		file ouf = new file("scanResults1.txt");

		file f = new file("cleaned_ibex_data.txt");
		int numLines = 413;
		f.initRead();
		double[] x = new double[numLines];
		double[] y = new double[numLines];
		String line = "";

		for (int i=0; i<numLines; i++) {
			line = f.readLine();
			StringTokenizer st = new StringTokenizer(line);
			x[i]=Double.parseDouble(st.nextToken());
			//x[i]-=1.0; // IMPORTANT _ SHIFT BY 1.0 TO MOVE TO HELIOCENTRIC FRAME!!
			y[i]=Double.parseDouble(st.nextToken());
			//System.out.println("row: " + i + " x: " + x[i] + " y: " + y[i]);
		}
		f.closeRead();

		// first scan:  V vs Long.
//		for (int i=0; i<
*/


		// LOAD THE REAL DATA
/*		file df = new file("second_pass.txt");
		int numLines2 = 918;
		df.initRead();
		double[] xD = new double[numLines2];
		double[] yD = new double[numLines2];
		String line = "";
		for (int i=0; i<numLines2; i++) {
			line = df.readLine();
			StringTokenizer st = new StringTokenizer(line);
			xD[i]=Double.parseDouble(st.nextToken());
			//x[i]-=1.0; // IMPORTANT _ SHIFT BY 1.0 TO MOVE TO HELIOCENTRIC FRAME!!
			yD[i]=Double.parseDouble(st.nextToken());
			System.out.println("row: " + i + " x: " + xD[i] + " y: " + yD[i]);
			//line = f.readLine();
			//line = f.readLine();
		}
		df.closeRead();
		// END LOADING REAL DATA
*/

		// MAKE SPIN PHASE SIMULATION
		file spFile = new file("spin_phase_sim.txt");
		spFile.initWrite();
		for (int tt = 0; tt<360; tt++) {
			// tt is the theta angle
			double ans = getRate(40.0, tt*Math.PI/360.0);
			spFile. write(tt+"\t"+ans+"\n");
		}
		spFile.closeWrite();



/*
		// MAKE SIMULATION FROM REAL DATA.. STEP WITH TIME AND PARAMS FOR MAKING 2D PARAM PLOT
		// loop through parameters here
		// standard params:
		double lamda = 74.7;
		//double temp = 6000.0;
		double v = 26000.0;
		int res = 30;
		double lamdaWidth = 20;
		//double tWidth = 10000;
		double vWidth = 10000;
		double lamdaMin = 65.0;
		//double tMin = 1000.0;
		double vMin=21000.0;

		double lamdaDelta = lamdaWidth/res;
		//double tDelta = tWidth/res;
		double vDelta = vWidth/res;

		lamda = lamdaMin;
		//temp=tMin;
		v = vMin;
		//file outF = new file("testVL_pass2_Out.txt");
		//outF.initWrite(false);
		for (int i=0; i<res; i++) {
			v=vMin;
			for (int j=0; j<res; j++) {

				file outf = new file("fitP2_vl"+lamda+"_"+v+".txt");
				outf.initWrite(false);
				setParams(lamda,v,0.015,6000.0);
				System.out.println(" now calc for: fitTTT_"+lamda+"_"+v+".txt");
				double [] doyD = new double[110];
				double [] rateD = new double[110];
				in doyIndex = 0;
				for (double doy=10.0 ; doy<65.0 ; doy+=0.5) {

					outf.write(doy+"\t");
					doyD[doyIndex]=doy;
					doyIndex++;
					//for (double off = 0.0 ; off<=3.0*Math.PI/2.0; off+=Math.PI/2.0) {
						//look_offset=off;
						look_offset=3.0*Math.PI/2.0;
						double rr = getRate(doy);
						System.out.println("trying: " + doy + "  got "+rr);

						outf.write(rr +"\t");
					//}
					//double rr = getRate(doy);
					//System.out.println("trying: " + doy + "  got "+rr);
					outf.write("\n");
				}
				outf.closeWrite();

				v+=vDelta;

			}
			//temp+=tDelta;
			lamda+=lamdaDelta;

		}

		//outF.closeWrite();
*/





/*

		for (int temptemp = 1000; temptemp<20000; temptemp+=2000) {
			file outf = new file("fitB_"+temptemp+"_7468_26000.txt");
			outf.initWrite(false);
			setParams(74.68,26000,0.015,(double)temptemp);
			for (double doy=10.0 ; doy<65.0 ; doy+=0.1) {
				outf.write(doy+"\t");

				//for (double off = 0.0 ; off<=3.0*Math.PI/2.0; off+=Math.PI/2.0) {
					//look_offset=off;
					look_offset=3.0*Math.PI/2.0;
					double rr = getRate(doy);
					System.out.println("trying: " + doy + "  got "+rr);
					outf.write(rr +"\t");
				//}
				//double rr = getRate(doy);
				//System.out.println("trying: " + doy + "  got "+rr);
				outf.write("\n");
			}
			outf.closeWrite();
		}
*/

	}
	// END CONSTRUCTOR


	public void setParams(double longitude, double speed, double dens, double temp) {

		currentLongitude = longitude;
		currentSpeed = speed;
		currentDens = dens;
		currentTemp = temp;
		// BUILD INTERSTERLLAR MEDIUM - helium and oxygen here, primary and secondary oxygen..
		//
		bulkHE = new HelioVector(HelioVector.SPHERICAL, speed, longitude*Math.PI/180.0, 95.5*Math.PI/180.0);
		//HelioVector bulkO1 = new HelioVector(HelioVector.SPHERICAL,26000.0, (74+180)*Math.PI/180.0, 95.5*Math.PI/180);
		GaussianVLISM gv1 = new GaussianVLISM(bulkHE,dens*100.0*100.0*100.0,temp);
		//KappaVLISM gv3 = new KappaVLISM(bulkO1, 0.00005*100*100*100, 6300, 1.6);
		ndHe = new NeutralDistribution(gv1, 0.0, he_ion_rate);  // mu is 0.0 for he
		ndHe.debug=false;
		// DONE CREATING MEDIUM
	}


	public void setParams(double dens) {
		currentDens = dens;
		// BUILD INTERSTERLLAR MEDIUM - helium and oxygen here, primary and secondary oxygen..
		//
		//HelioVector bulkHE = new HelioVector(HelioVector.SPHERICAL, currentSpeed, currentLongitude*Math.PI/180.0, 95.5*Math.PI/180.0);
		//HelioVector bulkO1 = new HelioVector(HelioVector.SPHERICAL,26000.0, (74+180)*Math.PI/180.0, 95.5*Math.PI/180);
		GaussianVLISM gv1 = new GaussianVLISM(bulkHE,dens*100.0*100.0*100.0,currentTemp);
		//KappaVLISM gv3 = new KappaVLISM(bulkO1, 0.00005*100*100*100, 6300, 1.6);
		ndHe = new NeutralDistribution(gv1, 0.0, he_ion_rate);  // mu is 0.0 for he
		ndHe.debug=false;
		// DONE CREATING MEDIUM
	}



	public void setParams(double dens, double temp) {

		currentDens = dens;
		currentTemp = temp;
		// BUILD INTERSTERLLAR MEDIUM - helium and oxygen here, primary and secondary oxygen..
		//
		//HelioVector bulkHE = new HelioVector(HelioVector.SPHERICAL, currentSpeed, currentLongitude*Math.PI/180.0, 95.5*Math.PI/180.0);
		//HelioVector bulkO1 = new HelioVector(HelioVector.SPHERICAL,26000.0, (74+180)*Math.PI/180.0, 95.5*Math.PI/180);
		GaussianVLISM gv1 = new GaussianVLISM(bulkHE,dens*100.0*100.0*100.0,temp);
		//KappaVLISM gv3 = new KappaVLISM(bulkO1, 0.00005*100*100*100, 6300, 1.6);
		ndHe = new NeutralDistribution(gv1, 0.0, he_ion_rate);  // mu is 0.0 for he
		ndHe.debug=false;
		// DONE CREATING MEDIUM
	}

	public void setParamsN(double dens, double ion) {

		currentDens = dens;
		currentTemp = temp;
		// BUILD INTERSTERLLAR MEDIUM - helium and oxygen here, primary and secondary oxygen..
		//
		//HelioVector bulkHE = new HelioVector(HelioVector.SPHERICAL, currentSpeed, currentLongitude*Math.PI/180.0, 95.5*Math.PI/180.0);
		//HelioVector bulkO1 = new HelioVector(HelioVector.SPHERICAL,26000.0, (74+180)*Math.PI/180.0, 95.5*Math.PI/180);
		GaussianVLISM gv1 = new GaussianVLISM(bulkHE,dens*100.0*100.0*100.0,temp);
		//KappaVLISM gv3 = new KappaVLISM(bulkO1, 0.00005*100*100*100, 6300, 1.6);
		ndHe = new NeutralDistribution(gv1, 0.0, he_ion_rate);  // mu is 0.0 for he
		ndHe.debug=false;
		// DONE CREATING MEDIUM
	}


	public double getRate(double dens, double doy) {
		if ((currentDens!=dens)) {
			System.out.println("new params: " + dens);
			setParams(dens);
		}
		return getRate(doy);
	}

	public double getRate(double dens, double temp, double doy) {
		if ((currentDens!=dens)|(currentTemp!=temp)) {
			System.out.println("new params: " + dens + " " + temp);
			setParams(dens,temp);
		}
		return getRate(doy);
	}

	public double getRate(double longitude, double speed, double dens, double temp, double doy) {
		if ((currentLongitude!=longitude)|(currentSpeed!=speed)|(currentDens!=dens)|(currentTemp!=temp)) {
			System.out.println("new params: " + longitude + " " + speed + " " + dens + " " + temp);
			setParams(longitude,speed,dens,temp);
		}
		return getRate(doy);
	}


	/* use this to set the look direction in latitude */
	double midThe = Math.PI/2.0;

	/**
	* Use this one to set the look direction in spin phase and then get the rate
	*
	*/
	public double getRate(double doy, double theta) {
		midThe = theta;
		getRate(doy);
	}


	/**
	* We want to call this from an optimization routine..
	*
	*/
	public double getRate(double doy) {


		// spacecraft position ..  assume at earth at doy
		//double pos = ((doy-80.0)/365.0);
		double pos = (doy-265.0)*(360.0/365.0)*Math.PI/180.0;// - Math.PI;

		// set look direction
		double midPhi = 0.0;
		midThe = Math.PI/2.0;
		double appogee = 0;

		// SET LOOK DIRECTION
		// THIS SET IS FOR FIRST PASS
		/*if (doy>10 && doy<15.5) {
			//orbit 13
			appogee = 13.1;
		}
		else if (doy>15.5 && doy<24.1) {
			//orbit 14
			appogee = 20.5;
		}
		else if (doy>24.1 && doy<32.0) {
			//orbit 15
			appogee = 28.1;
		}
		else if (doy>32.0 && doy<39.7) {
			//orbit 16
			appogee = 35.8;
		}
		else if (doy>39.7 && doy<47.5) {
			//orbit 17
			appogee = 43.2;
		}
		else if (doy>47.5 && doy<55.2) {
			//orbit 18
			appogee = 51.3;
		}
		else if (doy>55.2 && doy<62.6) {
			//orbit 19
			appogee = 58.8;
		}
		else if (doy>62.6 && doy<70.2) {
			//orbit 20
			appogee = 66.4;
		}*/


		// SECOND PASS
		if (doy>10.5 && doy<18.0) {
			//orbit 13
			appogee = 14.0;
		}
		else if (doy>18.0 && doy<25.6) {
			//orbit 14
			appogee = 21.6;
		}
		else if (doy>25.6 && doy<33.3) {
			//orbit 15
			appogee = 29.25;
		}
		else if (doy>33.3 && doy<40.85) {
			//orbit 16
			appogee = 36.8;
		}
		else if (doy>40.85 && doy<48.85) {
			//orbit 17
			appogee = 44.6;
		}
		else if (doy>48.85 && doy<56.4) {
			//orbit 18
			appogee = 52.4;
		}
		else if (doy>56.4 && doy<64.0) {
			//orbit 19
			appogee = 60.0;
		}
		else if (doy>64.0 && doy<71.6) {
			//orbit 20
			appogee = 67.6;
		}



		//orbit 21 - apogee= 73.9
		//orbit 22 - apogee = 81.8
		else {
			appogee=doy;
		}

		//doy = doy-183.0;
		//appogee = appogee-183.0;

		double he_ans_lr;
		midPhi = (appogee-3.0-265.0)*360.0/365.0*Math.PI/180.0 + look_offset;
		final double lookPhi = midPhi;

		// centered at look diretion..  convert doy+10 (dec 21-jan1) to angle
		//midPhi = (appogee+3.0-80.0)*360.0/365.0*Math.PI/180.0 + Math.PI/2;
		//midPhi = midPhi;


		// here's the now position
		final HelioVector posVec = new HelioVector(HelioVector.SPHERICAL, 1*AU, pos, Math.PI/2.0);
		final HelioVector scVelVec = new HelioVector(HelioVector.SPHERICAL, EARTHSPEED, pos+Math.PI/2.0, Math.PI/2.0);
		// placeholder for zero
		//final HelioVector scVelVec = new HelioVector(HelioVector.SPHERICAL, 0.0, pos+Math.PI/2.0, Math.PI/2.0);
		FunctionIII he_func = new FunctionIII() {  // helium neutral distribution
			public double function3(double v, double p, double t) {
				double tbr = activeArea*ndHe.dist(posVec, new HelioVector(HelioVector.SPHERICAL,v,p,t).sum(scVelVec) )*v*v*Math.sin(t)*v;  // extra v for flux
				tbr *= angleResponse(lookPhi,p);
				// that is the integrand of our v-space integration
				return tbr;
			}
		};


		// set limits for integration - spring only here
		//double[] he_lr_lims = { minSpeedsHe[0], maxSpeedsHe[0], midPhi-lrWidth, midPhi+lrWidth, midThe-lrWidth-spinWidth, midThe+lrWidth+spinWidth };
		// old limits for including all spin..
		double[] he_lr_lims = { minSpeedsHe[1], maxSpeedsHe[7], midPhi-lrWidth, midPhi+lrWidth, 3.0*Math.PI/8.0, 5.0*Math.PI/8.0 };

		// new limits for including just a single point in spin phase
		double[] he_lr_lims = { minSpeedsHe[1], maxSpeedsHe[7], midPhi-lrWidth, midPhi+lrWidth, midThe-lrWidth, midThe+lrWidth };


		he_lr_lims[0]=minSpeedsHe[3]; he_lr_lims[1]=maxSpeedsHe[5]; //, midPhi-lrWidth, midPhi+lrWidth, midThe-lrWidth-spinWidth, midThe+lrWidth+spinWidth };

		//PERFORM INTEGRATION
		he_ans_lr = 12.0*oneHour *1.0/8.0 *1.0  *   mi.mcIntegrate(he_func, he_lr_lims, mcN);
		he_ans_lr*=heSpringEff;
		return he_ans_lr;
	}

	public	double sigma=6.1;
	public double angleResponse(double centerA, double checkA) {
		double cA = centerA*180.0/Math.PI;
		double chA = checkA*180.0/Math.PI;
		double tbr = Math.exp(-(cA-chA)*(cA-chA)/sigma/sigma);
		return tbr;
	}


	public static final void main(String[] args) {
		IBEXLO_09fit il = new IBEXLO_09fit();
	}


	public static void o(String s) {
		System.out.println(s);
	}

	/*
	* This is to load a given model from a file and give a proper normalization to match the ibex data
	* for using a one parameter (normalization) fit.
	*/
	public class FitData {
		public StringTokenizer st;
		public double[] days;
		public double[] rates;
		public Vector daysV;
		public Vector ratesV;
		public FitData(String filename) {
			daysV = new Vector();
			ratesV = new Vector();
			file f = new file(filename);
			f.initRead();
			String line = "";
			while ((line=f.readLine())!=null) {
				st = new StringTokenizer(line);
				daysV.add(st.nextToken());
				ratesV.add(st.nextToken());
			}
			// time to fix the arrays
			days = new double[daysV.size()];
			rates = new double[daysV.size()];
			for (int i=0; i<days.length; i++) {
				days[i]=Double.parseDouble((String)daysV.elementAt(i));
				rates[i]=Double.parseDouble((String)ratesV.elementAt(i));
			}
		}

		// we are going to interpolate here
		public double getRate(double day) {
			for (int i=0; i<days.length; i++) {
				if (day<days[i]) { // this is where we want to be
					return (rates[i]+rates[i+1])/2;
				}
			}
			return 0;
		}
	}
}


/*
We need to keep track of Earth's Vernal Equinox
and use J2000 Ecliptic Coordinates!!!

March
2009 	20 	11:44
2010 	20 	17:32
2011 	20 	23:21
2012 	20 	05:14
2013 	20 	11:02
2014 	20 	16:57
2015 	20 	22:45
2016 	20 	04:30
2017 	20 	10:28

This gives the location of the current epoch zero ecliptic longitude..
However


/*

/*
Earth peri and aphelion

		perihelion			aphelion
2007 	January 3 	20:00 	July 7 	00:00
2008 	January 3 	00:00 	July 4 	08:00
2009 	January 4 	15:00 	July 4 	02:00
2010 	January 3 	00:00 	July 6 	12:00
2011 	January 3 	19:00 	July 4 	15:00
2012 	January 5 	01:00 	July 5 	04:00
2013 	January 2 	05:00 	July 5 	15:00
2014 	January 4 	12:00 	July 4 	00:00
2015 	January 4 	07:00 	July 6 	20:00
2016 	January 2 	23:00 	July 4 	16:00
2017 	January 4 	14:00 	July 3 	20:00
2018 	January 3 	06:00 	July 6 	17:00
2019 	January 3 	05:00 	July 4 	22:00
2020 	January 5 	08:00 	July 4 	12:00

*/