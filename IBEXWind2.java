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
public class IBEXWind2 {

	public int mcN = 10000;  // number of iterations per 3D integral
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


	/**
	*   830 - zurri carlie wednesday -
	*   patrick, melissa, tupac  Tuesday
	*  lara..  sophie..  vero..
	*/
	public IBEXWind2() {

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
		//bulkHE = new HelioVector(HelioVector.SPHERICAL, 28000.0, (180.0 + 74.68)*Math.PI/180.0, 95.5*Math.PI/180.0);
		//HelioVector bulkO1 = new HelioVector(HelioVector.SPHERICAL,26000.0, (74+180)*Math.PI/180.0, 95.5*Math.PI/180);
		//GaussianVLISM gv1 = new GaussianVLISM(bulkHE,0.015*100.0*100.0*100.0,10000.0);


		// test distribution, coming in here:

		//ndHe = new NeutralDistribution(gv1, 0.0, 0.0);  // mu is 0.0 for he
		currentLongitude = 74.0;
		currentSpeed = 28000.0;
		currentDens = 0.015;
		currentTemp = 100.0;
		// DONE CREATING MEDIUM

		o("earthspeed: " + EARTHSPEED);
		mi = new MultipleIntegration();
		o("done creating EARTHSPEED object");

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


		// LOAD THE REAL DATA - AZIMUTH or TIME SERIES
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

		// LOAD THE REAL DATA - Spin Phase Info - Orbit 61 for now!!!
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
/*		for (double ttt = 3000; ttt<

		file spFile = new file("doy_sim.txt");
		spFile.initWrite(false);


		//for (int ss = 20000; ss<100000; ss+=5000) {

		//	low_speed = ss;
		//	high_speed = ss+5000;
			//spFile.write("\n\n ss= "+ss);
		bulkHE = new HelioVector(HelioVector.SPHERICAL, 26000.0, (74.0)*Math.PI/180.0, 96.5*Math.PI/180.0);
		HelioVector.invert(bulkHE);
		// invert it, we want the real velocity vector out there at infinity :D

		GaussianVLISM gv1 = new GaussianVLISM(bulkHE,0.015*100.0*100.0*100.0, 3900.0);
		ndHe = new NeutralDistribution(gv1, 0.0, he_ion_rate);  // mu is 0.0 for he

		for (double tt = 50.0; tt<120; tt+=1) {
			//tt is the theta angle

			double ans = getRate(41.0, tt*Math.PI/180.0);
			spFile. write(tt+"\t"+Math.abs(ans)+"\n");
			o("trying tt = " + tt + " : "+ Math.abs(ans));
		}
		//}
		spFile.closeWrite();
*/

		// MAKE SPIN PHASE PARAM PLOT

		double v = 23000.0;
		double vwidth = 8000.0;
		double t = 3000.0;
		double twidth = 7000.0;
		int res = 30;  // resolution, 30 by 30 is enough apparently



		double tdelta = twidth/res;
		double vdelta = vwidth/res;

		double vmin = v;
		double tmin = t;

		// this one is for checking that we are comparing the right data
		file sampleFile = new file("sample_data.txt");
		sampleFile.initWrite(false);

		file spFile = new file("sp_param_65_74_6_5.txt");
		spFile.initWrite(false);
		// initialze file with header of x and y axis numbers
		for (int i=0; i<res; i++) {
			spFile.write(v+"\n");
			v+=vdelta;
		}
		v=vmin;
		for (int i=0; i<res; i++) {
			spFile.write(t+"\n");
			t+=tdelta;
		}
		t=tmin;

		for (int i=0; i<res; i++) {  // v loop
			t=tmin;
			for (int j=0; j<res; j++) {  // t loop
				o("trying temp: "+ t + " vel: " + v);

				// set up the interstellar medium flow model
				bulkHE = new HelioVector(HelioVector.SPHERICAL, v, (74.0)*Math.PI/180.0, 96.5*Math.PI/180.0);
				HelioVector.invert(bulkHE);
				// invert it, we want the real velocity vector out there at infinity :D

				GaussianVLISM gv1 = new GaussianVLISM(bulkHE,0.015*100.0*100.0*100.0, t);
				ndHe = new NeutralDistribution(gv1, 0.0, he_ion_rate);  // mu is 0.0 for he
				ndHe.debug=false;

				// store the results here.. fit them later and determine how well they fit the data
				double [] ans = new double[28];
				double ansMax = 0.0;
				double theta = 59;
				for (int k = 0; k<28; k++) {
					//tt is the theta angle

					// start of orbit 65 is doy 41 !!  that's where the good data is.. .
					ans[k] = getRate(41.0, theta*Math.PI/180.0);
					if (ans[k]>ansMax) ansMax=ans[k];
					//spFile. write(tt+"\t"+Math.abs(ans)+"\n");
					//o("trying theta = " + theta + " : "+ Math.abs(ans[k]) + " i: " + i + " j: " + j);
					theta+= 2.0;
				}

				// ok we have the answers, now lets compare to the loaded data
				double nfactor = (double)spDataMax/ansMax;

				double error = 0.0;
				sampleFile.write(t+ "\t" + v + "\n");
				for (int k = 0; k<28; k++) {
					ans[k]*=nfactor;
					error += Math.abs(ans[k]-spData[k]);
					sampleFile.write(ans[k] + "\t" + spData[k] + "\n");
				}

				// we got the error..  output to file
				//spFile.write(t + "\t" + v + "\t" + nfactor);
				spFile.write(error + "\n");
				o("nfactor data/model : "+ nfactor + " error: " + error);


				t += tdelta;


			}
			v += vdelta;
		}
		sampleFile.closeWrite();
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


	/*public double getRate(double dens, double temp, double doy) {
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
*/

	/* use this to set the look direction in latitude */
	double midThe = Math.PI/2.0;

	/**
	* Use this one to set the look direction in spin phase and then get the rate
	*
	*/
	public double getRate(double doy, double theta) {
		midThe = theta;
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


		// spacecraft position ..  assume at earth at doy
		final HelioVector pointVect = ei.getIbexPointing(ourDate);
		// test position for outside of heliosphere
		//final HelioVector pointVect = new HelioVector(HelioVector.CARTESIAN, 0.0,1.0,0.0);

		final HelioVector lookVect = new HelioVector(HelioVector.SPHERICAL,1.0, pointVect.getPhi()-Math.PI/2.0, pointVect.getTheta());
		//final HelioVector lookVect = new HelioVector(HelioVector.SPHERICAL,1.0, pointVect.getPhi()-Math.PI/2.0, midThe);

		// here's the now position (REAL)
		final HelioVector posVec = ei.getEarthPosition(ourDate);
		final HelioVector scVelVec = ei.getEarthVelocity(ourDate);

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
		FunctionIII he_func = new FunctionIII() {  // helium neutral distribution
			public double function3(double v, double p, double t) {
				HelioVector testVect = new HelioVector(HelioVector.SPHERICAL,v,p,t);
				HelioVector.sum(testVect,scVelVec);
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
		double[] he_lr_lims = { 0.0, maxSpeedsHe[4], midPhi-lrWidth, midPhi+lrWidth, midThe-lrWidth, midThe+lrWidth };


		//double[] he_lr_lims = { minSpeedsHe[2], maxSpeedsHe[5], midPhi-lrWidth, midPhi+lrWidth, 0, Math.PI };


		//he_lr_lims[0]=low_speed; he_lr_lims[1]=high_speed; //, midPhi-lrWidth, midPhi+lrWidth, midThe-lrWidth-spinWidth, midThe+lrWidth+spinWidth };

		//PERFORM INTEGRATION
		double he_ans_lr = 12.0*oneHour *1.0/8.0 *1.0  *   mi.mcIntegrate(he_func, he_lr_lims, mcN);
		he_ans_lr*=heSpringEff;
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

	// starts at 59, goes up by 2 !!
	int [] spData = {4,12,8,24,51,75,133,196,388,464,799,848,1221,1255,1048,1125,704,706,420,307,143,120,51,22,17,9,7,7,4,9,3,6};
	int spDataMax = 1255;

// spin phase data from 1 day start of orbit 65
/*
	59.0000136	4
	60.9999972	12
	63.0000024	8
	65.000004	24
	67.00002	51
	69	75
	71.0000016	133
	72.999996	196
	75.000012	388
	77.0000172	464
	79.000008	799
	81.000006	848
	83.000004	1221
	85.00002	1255
	87	1048
	89.0000052	1125
	91.0000104	704
	93.000012	706
	94.999992	420
	97.000008	307
	99.0000096	143
	101.0000148	120
	102.9999984	51
	105	22
	107.000016	17
	109.000014	9
	110.9999976	7
	113.0000028	7
	115.000008	4
	117.000024	9
	119.000004	3
	121.000002	6
*/






	public static final void main(String[] args) {
		IBEXWind2 il = new IBEXWind2();
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
