/**
* Use for calculations of SW scatter to Neutrals from Lunar Surface
*
*/
public class LunarScatter {

	/**
	*
	*/
	public LunarScatter() {

		MultipleIntegration mi = new MultipleIntegration();;
		// HERES WHAT WE NEED TO INTEGRATE

		final FunctionIII lunarFunc = new FunctionIII() {  // helium neutral distribution
			public double function3(double phi, double theta, double sza) {
				return f(phi,theta,sza)*sin(theta)*sin(sza);
			}
		};

		// limits of integration - full sphere
		double[] lunar_lims = { 0.0, 360.0, 0.0, 90.0, 0.0, 90.0 }; //TESTED SERIES OK 1/11

		//PERFORM INTEGRATION
		double answer = mi.mcIntegrate(lunarFunc, lunar_lims, 1000000)*Math.PI/2.0;
		System.out.println("Full sphere integral: " + answer);
		 answer = mi.mcIntegrate(lunarFunc, lunar_lims, 1000000)*Math.PI/2.0;
		System.out.println("Full sphere integral: " + answer);
		 answer = mi.mcIntegrate(lunarFunc, lunar_lims, 1000000)*Math.PI/2.0;
		System.out.println("Full sphere integral: " + answer);

		double[] lunar_lims2 = { 180-45, 180+45, 0.0, 90.0, 0.0, 90.0 }; //TESTED SERIES OK 1/11

		//PERFORM INTEGRATION
		answer = mi.mcIntegrate(lunarFunc, lunar_lims2, 1000000,10);
		System.out.println("Hemisphere Answer: " + answer);

		// test the function
		file outf = new file("scatter.txt");
		outf.initWrite(false);

		// look at three fixed szas, take phi=0, graph theta 0 to 90
		for (double theta=0.0; theta<90.0; theta+=1.0) {
			outf.write(theta+"\t"+f(0.0,theta,0.0)+"\t"+f(0.0,theta,45.0)+"\t"+f(0.0,theta,85.0)+"\n");
		}

		outf.write("\n\n\n\n\n\n\n");
		//  take phi=0,90,180 theta = 0 sza 0 to 90
		for (double sza=0.0; sza<90.0; sza+=1.0) {
			outf.write(sza+"\t"+f(0.0,0.0,sza)+"\t"+f(90.0,0.0,sza)+"\t"+f(180.0,0.0,sza)+"\n");
		}

		outf.write("\n\n\n\n\n\n\n");

		//  take theta=0,45,85 sza = 45 phi 0 to 360
		for (double phi=0.0; phi<360.0; phi+=2.0) {
			outf.write(phi+"\t"+f(phi,0.0,45.0)+"\t"+f(phi,45.0, 45.0)+"\t"+f(phi,85.0,45.0)+"\n");
		}

		outf.write("\n\n\n\n\n\n\n");
		//  make a polar plot, fix SZA=45, go up at phi=0, down at phi=180
		for (double theta=0.0; theta<90.0; theta+=1.0) {
			outf.write(f(0.0,theta,45.0)*cos(theta)+"\t"+f(0.0,theta,45.0)*sin(theta)+"\n");
		}
		for (double theta=90.0; theta>0.0; theta-=1.0) {
			outf.write(f(180.0,theta,45.0)*cos(theta)+"\t"+f(180.0,theta,45.0)*sin(theta)+"\n");
		}


		// this is IBEX integral
		final FunctionII lunarFunc2 = new FunctionII() {  // helium neutral distribution
			public double function2(double chi, double xi) {
				return f(180.0-chi+xi,90.0-chi,chi)*sin(chi);
			}
		};

		// limits of integration - full sphere
		double[] lunar_lims22 = { 0.0, 90.0, 0.0, 90.0 }; //TESTED SERIES OK 1/11

		//PERFORM INTEGRATION
		answer = mi.mcIntegrate(lunarFunc2, lunar_lims22, 1000000);
		System.out.println("ibex integral: " + answer);


		// this is normalization integral
		final FunctionIII lunarFuncNorm = new FunctionIII() {  // helium neutral distribution
			public double function3(double phi, double theta,double xi) {
				return f(phi,theta,xi)*sin(xi)*sin(theta);
			}
		};

		// limits of integration - full sphere
		double[] lunar_limsNorm = { 0.0, 360.0, 0.0, 90.0, 0.0, 90.0 }; //TESTED SERIES OK 1/11

		//PERFORM INTEGRATION
		answer = mi.mcIntegrate(lunarFuncNorm, lunar_limsNorm, 1000000,10);
		System.out.println("norm integral: " + answer);
		answer = mi.mcIntegrate(lunarFuncNorm, lunar_limsNorm, 1000000);
		System.out.println("norm integral: " + answer);
		answer = mi.mcIntegrate(lunarFuncNorm, lunar_limsNorm, 1000000,10);
		System.out.println("norm integral: " + answer);

		outf.closeWrite();
	}

	// z functions, show dependence on solar zenith angle GIVE IN DEGREES
	public double z0(double sza) { return 1.24*cos(1.44*sza - 48.70); }
	public double z1(double sza) { return Math.PI/180.0*(0.3*sza+1.72); }
	public double z2(double sza) { return 0.24*cos(74.48-1.52*sza); }
	public double z3(double sza) { return 90.0-1.03*sza; }


	// f functions, uniform component, azimuth dependent component, sunward-antisunward assym. component - polar angle piece
	public double f0(double phi, double theta, double sza) { return 0.009*z0(sza); }
	public double f1(double phi, double theta, double sza) { return z1(sza)*cos(2.0*phi)+(1.0-z1(sza)); }
	public double f2(double phi, double theta, double sza) { return z2(sza)*cos(phi) + (1.0-z2(sza)); }
	public double f3(double phi, double theta, double sza) { return sin(theta+z3(sza))*(1.0-z3(sza)/90.0) + z3(sza)/90.0; }


	/**
	* Scattering function
	* arguments in degrees
	*/
	public double f(double phi, double theta, double sza) {
		return f0(phi,theta,sza)*f1(phi,theta,sza)*f2(phi,theta,sza)*f3(phi,theta,sza);
	}

	// this sucks sorry convert to radians for Math.sin args
	public double sin(double a) { return Math.sin(a*Math.PI/180.0); }
	public double cos(double a) { return Math.cos(a*Math.PI/180.0); }

	public static final void main(String[] ags) {
		LunarScatter ls = new LunarScatter();
	}
}