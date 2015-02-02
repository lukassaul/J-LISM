import java.util.StringTokenizer;
import java.util.Vector;
import java.util.StringTokenizer;
import java.util.Date;
import java.util.*;


/**
*
* Simulating response at IBEX_LO with this class
*
* Create interstellar medium and perform integration
*
*/
public class IBEXWind {

	public int mcN = 100000;  // number of iterations per 3D integral
	public String outFileName = "lo_fit_09.txt";

	public static final double EV = 1.60218 * Math.pow(10, -19);
	public static double MP = 1.672621*Math.pow(10.0,-27.0);
	public static double AU = 1.49598* Math.pow(10,11); //meters
	public static double EARTHSPEED = AU*2.0*3.14/365.0/24.0/3600.0;
	public static double Ms = 1.98892 * Math.pow(10,30);  // kg
	//public static double Ms = 0.0;  // kg
	public static double G = 6.673 * Math.pow(10,-11);  // m^3/s^2/kg
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

	// angular acceptance - assume rectangular windows.. this is for integration limits only
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

	public double low_speed, high_speed;
	public int yearToAnalyze = 2009;

	public EarthIBEX ei = new EarthIBEX(yearToAnalyze);
	public TimeZone tz = TimeZone.getTimeZone("UTC");

	public file logFile = new file("IBEXWind_logfile2.txt");


	/**
	*   General constructor to do IBEX simulation
	*/
	public IBEXWind() {

		logFile.initWrite(true);

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
			//System.out.println(minSpeedsHe[i] +" = " + maxSpeedsHe[i] + " , " + minSpeedsO[i] + " = " + maxSpeedsO[i]);
		}

		System.out.println(minSpeedsHe[0] + "\t" + maxSpeedsHe[0] +"\t" + heSpeeds[2] + "\t" + heSpeeds[3]);

		// lets calculate the speed of the exact cold gas at 1AU
		// potential energy
		double pot = 4.0*MP*Ms*G/AU;
		o("pot: "+ pot);
		double energy1= 4.0*MP/2.0*28000*28000;
		o("energy1: " + energy1);
		double kEn = energy1+pot;
		double calcHeSpeed = Math.sqrt(2.0*kEn/4.0/MP);
		calcHeSpeed += EARTHSPEED;
		o("calcHeSpeed: "+ calcHeSpeed);





		// BUILD INTERSTERLLAR MEDIUM - helium and oxygen here, primary and secondary oxygen..
		//
		// test distribution, coming in here:
		bulkHE = new HelioVector(HelioVector.SPHERICAL, 28000.0, (74.68)*Math.PI/180.0, 85.0*Math.PI/180.0);

		// here's a test distribution putting in helium from x axis (vel in -x)
		//bulkHE = new HelioVector(HelioVector.CARTESIAN, 28000.0, 0.0, 0.0);

		//HelioVector.invert(bulkHE);
		// invert it, we want the real velocity vector out there at infinity :D

        // LISM PARAMS
		GaussianVLISM gv1 = new GaussianVLISM(bulkHE,0.015*100.0*100.0*100.0, 6000.0);
		KappaVLISM kv1 = new KappaVLISM(bulkHE,0.015*100.0*100.0*100.0, 6000.0, 2.4);
		//KappaVLISM gv3 = new KappaVLISM(bulkO1, 0.00005*100*100*100, 6300, 1.6);
		ndHe = new NeutralDistribution(kv1, 0.0, he_ion_rate);  // mu is 0.0 for he

		//ndHe = new NeutralDistribution(gv1, 0.0, 0.0);  // mu is 0.0 for he
		ndHe.debug=false;
		currentLongitude = 74.0;
		currentSpeed = 28000.0;
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
		//for (double vx = 0;


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
		// test here to match the ibex data and determine normalization
/*		file spFile = new file("doy_sim_test_31.txt");
		spFile.initWrite(false);


		//for (int ss = 20000; ss<100000; ss+=5000) {

		//	low_speed = ss;
		//	high_speed = ss+5000;
			//spFile.write("\n\n ss= "+ss);

		for (double tt = 68.0; tt<113.0; tt+=2) {
			//tt is the theta angle

			double ans = getRate(30.00, tt*Math.PI/180.0);
			spFile. write(tt+"\t"+Math.abs(ans)+"\n");
			o("trying tt = " + tt + " : "+ Math.abs(ans));
		}
		//
		spFile.closeWrite();

*/
		// OR MAKE TIME SERIES SIMULATION
/*		file outFF = new file("sp_sim_01.txt");
		outFF.initWrite(false);
		for (double doy=10.0 ; doy<70.0 ; doy+=1) {

			outFF.write(doy+"\t");
			double rr = getRate(doy);
			System.out.println("trying: " + doy + "  got "+rr);

			outFF.write(9*rr +"\t");  // TESTED WITH 2010 data...  factor 9 is about right for starting

			outFF.write("\n");
		}
		outFF.closeWrite();
*/




		// MAKE SIMULATION FROM REAL DATA.. STEP WITH TIME AND PARAMS FOR MAKING 2D PARAM PLOT
		// April 2014 for kappa sim now
		// loop through parameters here
		// standard params:

		int res = 30;
		double lamdaWidth = 20;
		double vWidth = 8000;
		double lamdaMin = 65.0;
		//double tMin = 1000.0;
		double vMin=22000.0;

		double lamdaDelta = lamdaWidth/res;
		double vDelta = vWidth/res;

		double lamda = lamdaMin;
		//temp=tMin;
		double v = vMin;

		// this time use a curvefitter in real time to fit the model to the data and report min error
		//CurveFitter cf = new CurveFitter();
		//cf.setFitData("2010_fit_time.txt");

		file outF = new file("kappa_test_2_4.txt");
		outF.initWrite(false);
		
		//for (int i=0; i<res; i++) {
		//	v=vMin;
		//	for (int j=0; j<res; j++) {
//
//				setParams(lamda,v,0.015,6000.0);
//				System.out.println(" now calc for params landa: "+lamda+" and v: "+v+".txt");

		double [] doyD = new double[51];
		double [] rateD = new double[51];

		int doyIndex = 0;
		for (double doy=12.0 ; doy<=62.0 ; doy+=1.0) {

			doyD[doyIndex]=doy;
			double rr = getRate(doy);
			outF.write(doy+"\t"+rr+"\n");
			System.out.println("trying: " + doy + "  got "+rr);
			rateD[doyIndex] = rr;
			doyIndex++;
		}

				// ok we have a set of model data, now let's bring the curvefitter out and go buck wild
			//	cf.setData(doyD,rateD);
			//	cf.doFit(CurveFitter.ONE_PARAM_FIT);
			//	System.out.println("did fit.. results: " + cf.minimum_pub + " at param " + cf.bestParams[0]);
			//	outF.write(lamda+"\t"+v+ "\t"+cf.minimum_pub+ "\t"+cf.bestParams[0]+ "\n");

			//	v+=vDelta;
			//}
			//temp+=tDelta;
			//lamda+=lamdaDelta;
		//}

		outF.closeWrite();


		// MAKE SIMULATION FROM REAL DATA.. STEP WITH SPIN PHASE AND PARAMS FOR MAKING 2D PARAM PLOT
		// feb 2011
		// loop through parameters here
		// standard params:
/*
		int res = 30;
		double tempWidth = 7000;
		double vWidth = 10000;
		double tempMin = 1000.0;
		//double tMin = 1000.0;
		double vMin=20000.0;

		double tempDelta = tempWidth/res;
		double vDelta = vWidth/res;

		double temp = tempMin;
		//temp=tMin;
		double v = vMin;

		// this time use a curvefitter in real time to fit the model to the data and report min error
		CurveFitter cf = new CurveFitter();
		cf.setFitData("2010_fit_sp.txt");

		file outF = new file("tVSv_75_85p7_1m.txt");
		outF.initWrite(false);
		for (int i=0; i<res; i++) {
			v=vMin;
			for (int j=0; j<res; j++) {

				bulkHE = new HelioVector(HelioVector.SPHERICAL, v, 75.0*Math.PI/180.0, 85.7*Math.PI/180.0);
				gv1 = new GaussianVLISM(bulkHE,100.0*100.0*100.0,temp);
				ndHe = new NeutralDistribution(gv1, 0.0, he_ion_rate);  // mu is 0.0 for he
				ndHe.debug=false;
				System.out.println(" now calc for params t: "+temp+" and v: "+v+".txt");

				double [] thetaD = new double[25];
				double [] rateD = new double[25];

				int doyIndex = 0;
				for (double theta=65.0 ; theta<=108.0 ; theta+=1.75) {

					thetaD[doyIndex]=theta;
					double rr = getRate(39.15, theta*Math.PI/180.0);
					System.out.println("trying: " + theta + "  got "+rr);
					rateD[doyIndex] = rr;
					doyIndex++;
				}

				System.out.println("doy index should be 25: "+ doyIndex);

				// ok we have a set of model data, now let's bring the curvefitter out and go buck wild
				cf.setData(thetaD,rateD);
				cf.doFit(CurveFitter.ONE_PARAM_FIT);
				System.out.println("did fit.. results: " + cf.minimum_pub + " at param " + cf.bestParams[0]);
				outF.write(temp+"\t"+v+ "\t"+cf.minimum_pub+ "\t"+cf.bestParams[0]+ "\n");

				v+=vDelta;
			}
			//temp+=tDelta;
			temp+=tempDelta;
		}

		outF.closeWrite();
*/


		// MAKE SIMULATION FROM REAL DATA.. STEP WITH SPIN PHASE AND PARAMS FOR MAKING 2D PARAM PLOT
		//  THIS TIME WE USE BOTH SPIN PHASE AND TIME SERIES TO GET CRAZY ON YO ASS
		// feb 2011
		// loop through parameters here
		// standard params:
/*
		int res = 30;
		double lamdaWidth = 20.0;
		double vWidth = 12000;
		double lamdaMin = 65.0;
		//double tMin = 1000.0;
		double vMin=20000.0;

		double lamdaDelta = lamdaWidth/res;
		double vDelta = vWidth/res;

		double lamda = lamdaMin;
		//temp=tMin;
		double v = vMin;

		// this time use a curvefitter in real time to fit the model to the data and report min error
		CurveFitter cf = new CurveFitter();
		cf.setFitData("2010_fit_sp.txt");

		CurveFitter cf2 = new CurveFitter();
		cf2.setFitData("2010_fit_time.txt");


		file outF = new file("lVSv_6000_86p5_200k.txt");
		outF.initWrite(false);
		for (int i=0; i<res; i++) {
			v=vMin;
			for (int j=0; j<res; j++) {

				bulkHE = new HelioVector(HelioVector.SPHERICAL, v, lamda*Math.PI/180.0, 86.5*Math.PI/180.0);
				gv1 = new GaussianVLISM(bulkHE,100.0*100.0*100.0,6000.0);
				ndHe = new NeutralDistribution(gv1, 0.0, he_ion_rate);  // mu is 0.0 for he
				ndHe.debug=false;
				System.out.println(" now calc for params lam: "+lamda+" and v: "+v+".txt");

				// first create the spin phase profile
				double [] thetaD = new double[25];
				double [] rateD = new double[25];

				int spIndex = 0;
				for (double theta=65.0 ; theta<=108.0 ; theta+=1.75) {

					thetaD[spIndex]=theta;
					double rr = getRate(39.15, theta*Math.PI/180.0);
					System.out.println("trying: " + theta + "  got "+rr);
					rateD[spIndex] = rr;
					spIndex++;
				}

				cf.setData(thetaD,rateD);
				cf.doFit(CurveFitter.ONE_PARAM_FIT);
				System.out.println("did fit.. results: " + cf.minimum_pub + " at param " + cf.bestParams[0]);
				//outF.write(temp+"\t"+v+ "\t"+cf.minimum_pub+ "\t"+cf.bestParams[0]+ "\n");

				// now create the time profile


				double [] doyC = new double[51];
				double [] rateC = new double[51];

				int doyIndex = 0;
				for (double doy=12.0 ; doy<=62.0 ; doy+=1.0) {

					doyC[doyIndex]=doy;
					double rr = getRate(doy);
					System.out.println("trying: " + doy + "  got "+rr);
					rateC[doyIndex] = rr;
					doyIndex++;
				}

				System.out.println("doy index should be 25: "+ doyIndex);

				// ok we have a set of model data, now let's bring the curvefitter out and go buck wild
				cf2.setData(doyC,rateC);
				cf2.doFit(CurveFitter.ONE_PARAM_FIT);
				System.out.println("did fit.. results: " + cf2.minimum_pub + " at param " + cf2.bestParams[0]);
				outF.write(lamda+"\t"+v+ "\t"+(cf.minimum_pub+cf2.minimum_pub)+"\n");

				v+=vDelta;
			}
			//temp+=tDelta;
			lamda+=lamdaDelta;
		}

		outF.closeWrite();
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
/*
	public void setParamsN(double dens, double ion) {

		currentDens = dens;
		//currentTemp = temp;
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
*/

/*	public double getRate(double dens, double doy) {
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
	*  We want to enable testing the rate for either spin phase or DOY
	*  thus we are going to add 10000 to the spin phase data, then we know
	*  which test we are doing
	*/
	public double getRateMulti(double n1, double n2, double t, double v, double phi, double theta, double doy_or_sp) {  // n,t,v, phi, theta
		logFile.initWrite(true);
		double tbr = 0.0;
        // params are density_series, density_sp, temp, v, lamda, theta

		bulkHE = new HelioVector(HelioVector.SPHERICAL, v, phi*Math.PI/180.0, theta*Math.PI/180.0);

		GaussianVLISM gv1 = new GaussianVLISM(bulkHE,100.0*100.0*100.0,t);
		//KappaVLISM gv3 = new KappaVLISM(bulkO1, 0.00005*100*100*100, 6300, 1.6);
		ndHe = new NeutralDistribution(gv1,0.0, he_ion_rate);  // mu is 0.0 for he
		ndHe.debug=false;

		if (doy_or_sp >= 10000) {
			// we are looking for spin phase at doy 41
			doSP = true;
			tbr = n2*getRate(39.15, (doy_or_sp-10000.0)*Math.PI/180.0);
		}

		else {
			doSP = false;
			midThe = Math.PI/2.0;
			tbr = n1*getRate(doy_or_sp); // a real DOY
		}

		logFile.write(n1+"\t"+n2+"\t"+t+"\t"+v+"\t"+phi+"\t"+theta+"\t"+doy_or_sp+"\t"+tbr+"\n");
		logFile.closeWrite();
		return tbr;

	}


	public boolean doSP = false;
	/**
	* Use this one to set the look direction in spin phase and then get the rate
	*
	*/
	public double getRate(double doy, double theta) {
		midThe = theta;
		doSP = true;
		return getRate(doy);
	}


	/**
	* We want to call this from an optimization routine..
	*
	*/
	public double getRate(double doy) {
		int day = (int)doy;
		double fract = doy-day;
		//System.out.println(fract);
		int hour = (int)(fract*24.0);


		Calendar c = Calendar.getInstance();
		c.setTimeZone(tz);
		c.set(Calendar.YEAR, yearToAnalyze);
		c.set(Calendar.DAY_OF_YEAR, day);
		c.set(Calendar.HOUR_OF_DAY, hour);
		//c.set(Calendar.MINUTE, 0);
		//c.set(Calendar.SECOND, 0);
		Date ourDate = c.getTime();
		//System.out.println("trying getRate w/ " + doy + " the test Date: " +  ourDate.toString());

		if (!doSP) midThe = Math.PI/2.0;  // we look centered if we aren't doing a scan

		// spacecraft position ..  assume at earth at doy
		final HelioVector pointVect = ei.getIbexPointing(ourDate);
		// test position for outside of heliosphere
		//final HelioVector pointVect = new HelioVector(HelioVector.CARTESIAN, 0.0,1.0,0.0);

		//if (doSP)
		final HelioVector lookVect = new HelioVector(HelioVector.SPHERICAL,1.0, pointVect.getPhi()+Math.PI/2.0, midThe);
		//else final HelioVector lookVect = new HelioVector(HelioVector.SPHERICAL,1.0, pointVect.getPhi()+Math.PI/2.0, pointVect.getTheta()); // TESTED SERIES OK !/11
		//final HelioVector lookVect = new HelioVector(HelioVector.SPHERICAL,1.0, pointVect.getPhi()-Math.PI/2.0, midThe);

		// here's the now position (REAL)
		final HelioVector posVec = ei.getEarthPosition(ourDate);
		final HelioVector scVelVec = ei.getEarthVelocity(ourDate);

		// temp! set it to zero for debug
//		final HelioVector scVelVec = new HelioVector(HelioVector.CARTESIAN, 0.0, 0.0, 0.0);

		//System.out.println("trying getRate w/ " + doy + " the test Date: " +  ourDate.toString() + " posvec: " + posVec.toAuString() +
		//	" scvelvec: " + scVelVec.toKmString() + " midthe " + midThe);

		// use a test position that puts the Earth out in the VLISM .. lets go with 200AU
		//final HelioVector posVec = new HelioVector(HelioVector.CARTESIAN, 200*AU, 0.0, 0.0); // due upwind
		//final HelioVector scVelVec = new HelioVector(HelioVector.CARTESIAN, 0.0, 0.0, 0.0); // not moving



		double midPhi = lookVect.getPhi();
		//final HelioVector posVec = new HelioVector(HelioVector.SPHERICAL, 1*AU, pos, Math.PI/2.0);
		//final HelioVector scVelVec = new HelioVector(HelioVector.SPHERICAL, EARTHSPEED, pos+Math.PI/2.0, Math.PI/2.0);

		// placeholder for zero
		//final HelioVector scVelVec = new HelioVector(HelioVector.SPHERICAL, 0.0, pos+Math.PI/2.0, Math.PI/2.0);

		// HERES WHAT WE NEED TO INTEGRATE
		final FunctionIII he_func = new FunctionIII() {  // helium neutral distribution
			public double function3(double v, double p, double t) {
				HelioVector testVect = new HelioVector(HelioVector.SPHERICAL,v,p,t);
				HelioVector.difference(testVect,scVelVec); // TESTED SERIES OK 1/11
				double tbr = activeArea*ndHe.dist(posVec, testVect)*v*v*Math.sin(t)*v;  // extra v for flux
				tbr *= angleResponse(testVect,lookVect);
				//tbr *= angleResponse(lookPhi,p);
				// that is the integrand of our v-space integration
				return tbr;
			}
		};


		// set limits for integration - spring only here
		//double[] he_lr_lims = { minSpeedsHe[0], maxSpeedsHe[0], midPhi-lrWidth, midPhi+lrWidth, midThe-lrWidth-spinWidth, midThe+lrWidth+spinWidth };
		// old limits for including all spin..
		//double[] he_lr_lims = { minSpeedsHe[1], maxSpeedsHe[7], midPhi-lrWidth, midPhi+lrWidth, 3.0*Math.PI/8.0, 5.0*Math.PI/8.0 };

		// new limits for including just a single point in spin phase

		double[] he_lr_lims = { 0.0, maxSpeedsHe[5], midPhi-lrWidth, midPhi+lrWidth, 0.0, Math.PI }; //TESTED SERIES OK 1/11

		if (doSP) {
			he_lr_lims[4]=midThe-2*lrWidth;
			he_lr_lims[5]=midThe+2*lrWidth;
			//lookVect = new HelioVector(HelioVector.SPHERICAL,lookVect.getR(), lookVect.getPhi(), midThe);
		}



		// this needs to change for the case of spin phase simulation.. w


		//double[] he_lr_lims = { minSpeedsHe[2], maxSpeedsHe[5], midPhi-lrWidth, midPhi+lrWidth, 0, Math.PI };


		//he_lr_lims[0]=low_speed; he_lr_lims[1]=high_speed; //, midPhi-lrWidth, midPhi+lrWidth, midThe-lrWidth-spinWidth, midThe+lrWidth+spinWidth };

		//PERFORM INTEGRATION
		double he_ans_lr = 12.0*oneHour *1.0/8.0 *1.0  *   mi.mcIntegrate(he_func, he_lr_lims, mcN);
		he_ans_lr*=heSpringEff;
		doSP = false;  //
		return he_ans_lr;
	}

	public	double sigma=6.1;

	/**
	*  From calibrated response
	*
	*/
	public double angleResponse(double centerA, double checkA) {

		double cA = centerA*180.0/Math.PI;
		double chA = checkA*180.0/Math.PI;


		// This is the Gaussian Response
		double tbr = Math.exp(-(cA-chA)*(cA-chA)/sigma/sigma);
		return tbr;
	}

	public double angleResponse(HelioVector center, HelioVector look) {
		double angle = center.getAngle(look);
		angle = Math.abs(angle * 180.0/Math.PI);
		// This is the Triangle Response
		double tbr = 0.0;
		if (angle<7.0) {
			double m = -1.0/7.0;
			tbr = m*angle+1.0;
		}

		// This is the Gaussian Response
		//double tbr = Math.exp(-angle*angle/sigma/sigma);
		return tbr;
	}

	public static final void main(String[] args) {
		IBEXWind il = new IBEXWind();
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
