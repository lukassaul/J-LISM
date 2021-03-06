import java.lang.Math;
import java.util.Date;
import java.util.Random;
import flanagan.integration.*;
//import drasys.or.nonlinear.*; // our friend mr. simpson resides here

/**
*  This class should take care of doing multi-dimensional numeric integrations.
*  This will be slow!!
*
*   Actually not half bad...
*
*
*  Lukas Saul
*  UNH Physics
*  May, 2001
*
*   Updated Aug. 2003 - only works with 2D integrals for now...
*
*   Updated Aug. 2004 - 3D integrals OK!  Did that really take 1 year??
*
*   About 2.1 seconds to integrate e^(-x^2-y^2-z^2) !!
*
*  Testing in jan 2008 again..
*
*    16 points seems to work well enough for 4 digit precision of 3D gaussian integral in < 20ms
*/
public class MultipleIntegration {

	public static double PLUS_INFINITY = Double.POSITIVE_INFINITY;
	public static double MINUS_INFINITY = Double.NEGATIVE_INFINITY;
	private double result;
	public int npoints=128;

	Random r = new Random();

	/**
	* for testing - increments at every 1D integration performed
	*/
	public long counter, counter2;

	/**
	* Constructor just sets up the 1D integration routine for running
	*  after which all integrations may be called.
	*/
	public MultipleIntegration() {
		counter = 0;
		counter2 = 0;
		result = 0.0;
	}

	/**
	* sets counter to zero
	*/
	public void reset() {
		counter = 0;
	}
	public void reset2() {
		counter2 = 0;
	}

	/**
	* Set accuracy of integration
	*/
	public void setEpsilon(double d) {
		//s.setErrMax(d);
	}


	/**
	* Deal with infinit limits here and call Flanagan gaussian quadrature for numeric integration
	*
	*/
	public double gaussQuad(final FunctionI fI, final double min_l, final double max_l, final int np) {
		try{ // if inifinite limits, replace arg x by x/1-x and integrate -1 to 1

			if (min_l!=MINUS_INFINITY && max_l!=PLUS_INFINITY)
				return Integration.gaussQuad(fI, min_l, max_l, np);

			else if (min_l==MINUS_INFINITY && max_l==PLUS_INFINITY) {
				// lets try recursion instead...

				return
					gaussQuad(fI, min_l, -10000.0, np/4) +
					gaussQuad(fI, -10000.0, 0.0, np/4) +
					gaussQuad(fI, 0.0, 10000.0, np/4) +
					gaussQuad(fI, 10000.0, max_l, np/4);


				/*FunctionI newFI = new FunctionI () {
					public double function(double x) {
						return (fI.function((1.0-x)/x)+fI.function((-1.0+x)/x))/x/x;
					}
				};
				return Integration.gaussQuad(newFI,0.0,1.0,np);
				*/
			}
			else if (min_l==MINUS_INFINITY) {
				FunctionI newFI = new FunctionI () {
					public double function(double x) {
						return fI.function(max_l-(1.0-x)/x)/x/x;
					}
				};
				return Integration.gaussQuad(newFI,0.0,1.0,np);
			}
			else if (max_l==PLUS_INFINITY) {
				FunctionI newFI = new FunctionI () {
					public double function(double x) {
						return fI.function(min_l+(1.0-x)/x)/x/x;
					}
				};
				return Integration.gaussQuad(newFI,0.0,1.0,np);
			}
		}
		catch (Exception e) {
			e.printStackTrace();
			return 0.0;
		}
		return 0.0;
	}

	/**
	* Here's the goods, for 3D integrations.
	*
	* Limits are in order as folows:  zlow, zhigh, ylow, yhigh, xlow, xhigh
	*
	*/
	public double integrate(final FunctionIII f, final double[] limits, final int np) {
		reset();
	//	System.out.println("Called 3D integrator");
	//	System.out.println("Integrating from: \n"+limits[0]+" "+limits[1]+
	//				"\n"+limits[2]+" "+limits[3]+"\n"+limits[4]+" "+limits[5]);
		double[] nextLims = new double[4];
		for (int i=0; i<4; i++) {
			nextLims[i] = limits[i+2];
		}

		final double[] nextLimits = nextLims;

		FunctionI funcII = new FunctionI () {
			public double function(double x) {
				return integrate(f.getFunctionII(2,x), nextLimits, np);
			}
		};
		result = gaussQuad(funcII,limits[0],limits[1],np);
		return result;

	}

	/**
	* Here's the goods, for 2D integrations
	*/
	public double integrate(final FunctionII f, final double[] limits, final int np) {
		FunctionI f1 = new FunctionI() {
			public double function(double x) {
				return integrate(f.getFunction(1,x),limits[2],limits[3], np);
			}
		};
		result = gaussQuad(f1, limits[0], limits[1], np);
		return result;
		// each call to f1 by the intgrator does an integration
	}


	/**
	* Here's the simple goods, for 1D integrations
	*   courtesy of our friends at drasys.or.nonlinear
	*/
	public double integrate(final FunctionI f, double lowLimit, double highLimit, int np) {
		//System.out.println("Called 1D integrator");
		counter2++;
		counter++;
		//if (counter%10000==1) System.out.println("Counter: " + counter);
		//s.setInterval(lowLimit,highLimit);
		//return s.getResult(f);
		return gaussQuad(f,lowLimit,highLimit,np);
		//return
	}

	/**
	* Monte-Carlo 1D Integration (not a good idea but for test)
	*/
	public double mcIntegrate(final FunctionI f, double lowLimit, double highLimit, int np) {
		if (lowLimit == MINUS_INFINITY && highLimit == PLUS_INFINITY) {
			FunctionI newFI = new FunctionI () {
				public double function(double x) {
					return (f.function((1.0-x)/x)+f.function((-1.0+x)/x))/x/x;
				}
			};
			return mcIntegrate(newFI, 0.0, 1.0, np);
		}
		else {
			double tbr = 0.0;
			for (int i=0; i<np; i++) {
				// select random point
				double point = r.nextDouble()*(highLimit-lowLimit) + lowLimit;
				tbr += f.function(point);
			}
			tbr*=(highLimit-lowLimit);
			tbr/=np;
			return tbr;
		}
	}
	public double mcIntegrate(final FunctionI f, double center, double sigma,
						double lowLimit, double highLimit, int np) {

		double tbr = 0.0;
		for (int i=0; i<np; i++) {
			// select random point
			double point=lowLimit-1.0;
			while (point<lowLimit || point>highLimit) {
				point = r.nextGaussian()*sigma+center;
			}
			tbr += f.function(point)/gauss(point,center,sigma);
		}
		//tbr *= (highLimit-lowLimit);
		tbr/=np;
		return tbr;

	}
	public double ntemp = Math.sqrt(2.0*Math.PI);
	public double gauss(double x, double c, double s) {
		return 1/s/ntemp*Math.exp(-(x-c)*(x-c)/2/s/s);
	}

	/**
	* Monte carlo 2D integral
	*/
	public double mcIntegrate(final FunctionII f, double[] limits, int np) {
		r.setSeed((new Date()).getTime());
		double w1 = limits[1]-limits[0];
		double w2 = limits[3]-limits[2];
		double tbr = 0.0;
		double xpoint,ypoint,zpoint;
		System.out.println("widths: " + w1 + " " + w2);
		for (int i=0; i<np; i++) {
			// select random point
			xpoint = r.nextDouble()*w1 + limits[0];
			ypoint = r.nextDouble()*w2 + limits[2];
			tbr += f.function2(xpoint,ypoint);
		}
		tbr*=w1*w2;
		tbr=tbr/(double)np;
		return tbr;
	}

	/**
	* Monte carlo 3D integral
	*/
	public double mcIntegrate(final FunctionIII f, double[] limits, int np) {

		if (limits[0] == MINUS_INFINITY && limits[1] == PLUS_INFINITY &&
			limits[2] == MINUS_INFINITY && limits[3] == PLUS_INFINITY &&
			limits[4] == MINUS_INFINITY && limits[5] == PLUS_INFINITY) {

			System.out.println("doing infinite MC");
			FunctionIII newFI = new FunctionIII () {
				public double function(double x,double y, double z) {
					double xp = (1.0-x)/x;
					double yp = (1.0-y)/y;
					double zp = (1.0-z)/z;
					return (f.function3(xp,yp,zp)+f.function3(-xp,yp,zp)+
							f.function3(xp,-yp,zp)+f.function3(-xp,-yp,zp)+
							f.function3(xp,yp,-zp)+f.function3(-xp,yp,-zp)+
							f.function3(xp,-yp,-zp)+f.function3(-xp,-yp,-zp))/x/x/y/y/z/z;
				}
			};
			double [] newLims = {0.0,1.0, 0.0,1.0, 0.0,1.0};
			return mcIntegrate(newFI, newLims, np);
		}
		else {
			r.setSeed((new Date()).getTime());
			double w1 = limits[1]-limits[0];
			double w2 = limits[3]-limits[2];
			double w3 = limits[5]-limits[4];
			double tbr = 0.0;
			double xpoint,ypoint,zpoint;
			System.out.println("widths: " + w1 + " " + w2 + " " + w3);
			for (int i=0; i<np; i++) {
				// select random point
				xpoint = r.nextDouble()*w1 + limits[0];
				ypoint = r.nextDouble()*w2 + limits[2];
				zpoint = r.nextDouble()*w3 + limits[4];
				tbr += f.function3(xpoint,ypoint,zpoint);
			}
			tbr*=w1*w2*w3;
			tbr=tbr/(double)np;
			return tbr;
		}
	}





	/**
	* Just for testing only here!!!
	*
	* seems to work - may 3, 2001
	*
	*  lots more testing for limits of 3d , polar vs. cartesian -  Oct. 2004
	*/
	public static final void main(String[] args) {

		MultipleIntegration mi = new MultipleIntegration();

		// some functions to integrate!!
		FunctionI testf1 = new FunctionI() {
			public double function(double x) {
				return Math.exp(-(x)*(x));
			}
		};
		FunctionIII testf3 = new FunctionIII () {
			public double function3(double x, double y, double z) {
				return (Math.exp(-(x*x+y*y+z*z)));
			}
		};

		//System.out.println("f1 gauss: "+mi.integrate(testf1,MINUS_INFINITY,PLUS_INFINITY,10000));
		//System.out.println("f1 gauss mc: "+mi.mcIntegrate(testf1,MINUS_INFINITY,PLUS_INFINITY,1000000));
		System.out.println("sqrt pi: "+ Math.sqrt(Math.PI));
		System.out.println("mc 10000: " + mi.mcIntegrate(testf1,-10.0,10.0,1000));
		System.out.println("mc weight: " + mi.mcIntegrate(testf1,0.0,1.0,-10.0,10.0,1000));

		//System.out.println(""+mi.integrate(testf1,MINUS_INFINITY,0.0,512));
		//System.out.println(""+mi.integrate(testf1,-100.0,0.0,512));
		double[] lims = {MINUS_INFINITY, PLUS_INFINITY, MINUS_INFINITY, PLUS_INFINITY, MINUS_INFINITY, PLUS_INFINITY};
		double[] lims2 = {-10.0,10.0,-10.0,10.0,-10.0,10.0};
	//	double test3 = mi.integrate(testf3,lims2,64);
	//	double test32 = mi.mcIntegrate(testf3,lims,1000000);
	//	double test323 = mi.mcIntegrate(testf3,lims2,1000000);
	//	System.out.println("test3 reg to 10: " + test3);
	//	System.out.println("test3 mc: " + test32);
	//	System.out.println("test3 mc to 10: " + test323);
	//	System.out.println("pi to 3/2: "+ Math.pow(Math.PI,3.0/2.0));

		/*
		// FOR TESTING 3D INTEGRALS
		file f = new file("int_test.txt");
		f.initWrite(false);

		for (int i=1; i<256; i*=2) {
			//System.out.println("lim: " + i);
			Date d1 = new Date();
			double test1 = mi.integrate(testf1, -100,100,i);
			Date d2 = new Date();
			//System.out.println("1D integrations: " + mi.counter);
			double[] lims = {MINUS_INFINITY, PLUS_INFINITY, MINUS_INFINITY, PLUS_INFINITY, MINUS_INFINITY, PLUS_INFINITY};
			double test3 = mi.integrate(testf3,lims,i);
			Date d3 = new Date();
			//System.out.println("3D integrations: " + mi.counter);

			//System.out.println("Answer single gauss: " + test1);
			//System.out.println("Answer trip gauss: " + test3);
			//System.out.println("took: " + (d2.getTime()-d1.getTime()));
			//System.out.println("took 3: " + (d3.getTime()-d2.getTime()));
			f.write(i+"\t"+test3+"\t"+(d3.getTime()-d2.getTime()) + "\n");
		}
		f.closeWrite();
		*/



		/*
		double[] lims2 = new double[4];
		lims2[0]=-100.0;
		lims2[1]=100.0;
		lims2[2]=-100.0;
		lims2[3]=100.0;

		Date d3 = new Date();
		double test2 = mi.integrate(testf2,lims2);
		Date d4 = new Date();

		System.out.println("Answer frm 2d testf2: " + test2);
		System.out.println("took: " + (d4.getTime()-d3.getTime()));
		System.out.println("1D integrations: " + mi.counter);


		d3 = new Date();
		test2 = mi.integrate(testf2a,lims2);
		d4 = new Date();
		System.out.println("Answer frm 2d testf2a: " + test2);
		System.out.println("took: " + (d4.getTime()-d3.getTime()));
		System.out.println("1D integrations: " + mi.counter);



		System.out.println("trying polar 2d now");
		lims2 = new double[4];
		lims2[0]=0;
		lims2[1]=2*Math.PI;
		lims2[2]=0;
		lims2[3]=10;

		d3 = new Date();
		double ttest = mi.integrate(testf2b,lims2);
		d4 = new Date();

		System.out.println("2d polar Answer: " + ttest);
		System.out.println("took: " + (d4.getTime()-d3.getTime()));
		System.out.println("1D integrations: " + mi.counter);





		System.out.println("trying 3d now... ");
		// basic limit test here,
		double[] lims = new double[6];
		lims[0]=0.0;
		lims[1]=3.00;
		lims[2]=0.0;
		lims[3]=1.0;
		lims[4]=0.0;
		lims[5]=2.0;

		Date d1 = new Date();
		double test = mi.integrate(testlims,lims);
		Date d2 = new Date();

		System.out.println("Answer: " + test);
		System.out.println("took: " + (d2.getTime()-d1.getTime()));
		System.out.println("1D integrations: " + mi.counter);
		System.out.println("answer: " + 8*81/2/3/4);



		lims = new double[6];
		lims[0]=-10.0;
		lims[1]=10.00;
		lims[2]=-10.0;
		lims[3]=10.0;
		lims[4]=-10.0;
		lims[5]=10.0;

		d1 = new Date();
		test = mi.integrate(testf,lims);
		d2 = new Date();

		System.out.println("Answer: " + test);
		System.out.println("took: " + (d2.getTime()-d1.getTime()));
		System.out.println("1D integrations: " + mi.counter);
		System.out.println("test^2/3: " + Math.pow(test,2.0/3.0));
		System.out.println("trying 3d now... ");



		// 3d Function integration working in spherical coords??
	    lims = new double[6];
		lims[0]=0;
		lims[1]=Math.PI;
		lims[2]=0;
		lims[3]=2*Math.PI;
		lims[4]=0;
		lims[5]=10.0;

		d1 = new Date();
		test = mi.integrate(testfs,lims);
		d2 = new Date();

		System.out.println("Answer: " + test);
		System.out.println("took: " + (d2.getTime()-d1.getTime()));
		System.out.println("1D integrations: " + mi.counter);
		System.out.println("test^2/3: " + Math.pow(test,2.0/3.0));


		*/

	}
}