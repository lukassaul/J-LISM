
import java.util.*;

/**
*  Lukas Saul, UNH Space Science Dept.
*  Curvefitting utilities here.
*
*    Adapted to use simplex method from flanagan java math libraries - March 2006
*
*
*  >
*   --  added step function and pulled out of package IJ -- lukas saul, 11/04
*
*   -- adding EXHAUSTIVE SEARCH ability..  for nonlinear fits
*
*
*   -- adding HEMISPHERIC -  analytic fit to hemispheric model w/ radial field
*      see Hemispheric.java for info
*
*   -- adding quadratic that passes through origin for PUI / Np correlation
*
*   -- adding SEMI - analytic fit assuming constant adiabatic acceleration
*      that uses legendre polynomials for pitch angle diffusion.  see SemiModel.java
*      March, 2006
*
*/
public class CurveFitter implements MinimisationFunction {


    public static final int STRAIGHT_LINE=0,POLY2=1,POLY3=2,POLY4=3,
    EXPONENTIAL=4,POWER=5,LOG=6,RODBARD=7,GAMMA_VARIATE=8,
    LOG2=9, STEP=10, HEMISPHERIC=11, HTEST=12, QUAD_ORIGIN=13,
    GAUSSIAN=14, SEMI=15, SEMI_1=16, IBEX=17, /** this one is custom, change at will */ IBEXFIT=18, IBEX1=19,
    SPINPHASE = 20,COMBO = 21, DOUBLE_NORM = 22, ONE_PARAM_FIT = 23, A2DPOLY = 24, A3DPOLY = 25;

    public static final String[] fitList = {"Straight Line","2nd Degree Polynomial",
    	"3rd Degree Polynomial", "4th Degree Polynomial","Exponential","Power",
    	"log","Rodbard", "Gamma Variate", "y = a+b*ln(x-c)", "Step Function",
    	"Hemispheric", "H-test", "Quadratic Through Origin","Gaussian","Semi","Semi-one-param",
    	"IBEX_sim","IBEX_fitfit","1 param IBEX fit"," T,V,lat fit ", "T, v,lamda,beta fit",
    	"two norm", "one norm","2D polynomial fit","3D polynomial fit"};

    public static final String[] fList = {"y = a+bx","y = a+bx+cx^2",
    	"y = a+bx+cx^2+dx^3", "y = a+bx+cx^2+dx^3+ex^4","y = a*exp(bx)","y = ax^b",
    	"y = a*ln(bx)", "y = d+(a-d)/(1+(x/c)^b)", "y = a*(x-b)^c*exp(-(x-b)/d)",
    	"y = a+b*ln(x-c)", "y = a*step(x-b)", "y=hem.f(a,x)",
    	"y = norm*step()..","y=ax+bx^2", "y=a*EXP(-(x-b)^2/c)","y=semi.f(a,b,x)","y=ibex.f(temp,norm,pos,doy(a,b,c,x))",
    	"y = a*fit(x)", "y=ibex.f(norm)", "y=ibex.f(t,v,lat,sp)", "y=ibex.f(t,v,lamda,beta",
    	"y=a*stuff orb*stuff", "y=a*stuff", "y=sum(Aijk)(x^i+y}","y=sum(Aijk)(x^i+y^j+z^k)"};

    private static final double root2 = Math.sqrt(2.0);//1.414214; // square root of 2

    private int fit;                // Number of curve type to fit
    private double[] xData, yData;  // x,y data to fit

	public double[] bestParams; // take this after doing exhaustive for the answer
	public double minimum_pub;
	private double[] gridDataX, gridDataY, gridDataZ; // these are the value on the ticks of axis
	private double[][] zData; // this is the data to fit for our multiD fit
	private double[][][] zzData;



	private int ord = 5; // order for polynomial fit

	//public Hemispheric hh;
	public SemiModel sm;
	public AdiabaticModel am;
	public IBEXLO_09fit ib;
	public FitData ff;
	public IBEXWind iw, ibw;

	public double[] start, step;
	public double ftol;

	public double gamma = 1.5;
	public Minimisation min;
	public String fitFile;
	public FitData fitData;

	public CurveFitter() {
		min = new Minimisation();
		// beware to set data before using this one!
	}

	/**
	* build it with floats and cast
	*/
    public CurveFitter (float[] xDat, float[] yDat) {
		xData = new double[xDat.length];
		yData = new double[yDat.length];
		for (int i=0; i<xData.length; i++) {
			xData[i]=(double)xDat[i];
			yData[i]=(double)yDat[i];
		}
		min = new Minimisation();
		//numPoints = xData.length;
	}

	public CurveFitter (double[] xData, double[] yData, String fname) {
		fitFile = fname;
        this.xData = xData;
        this.yData = yData;
		min = new Minimisation();
	}


	public CurveFitter (double[] xData, double[] yData, FitData fname) {
		fitData = fname;
		this.xData = xData;
		this.yData = yData;
		min = new Minimisation();
	}


    /** Construct a new CurveFitter. */
    public CurveFitter (double[] xData, double[] yData) {
        this.xData = xData;
        this.yData = yData;
		min = new Minimisation();
        //numPoints = xData.length;
    }

    /**
    * Construct a new CurveFitter.
    *
    * this time with some params for the ibex fitter
    */
    public CurveFitter (double[] xData, double[] yData, double lam, double vv) {
        this.xData = xData;
        this.yData = yData;
		min = new Minimisation();
		ib = new IBEXLO_09fit();
		// public void setParams(double longitude, double speed, double dens, double temp) {
		ib.setParams(lam, vv, ib.currentDens, ib.currentTemp);

        //numPoints = xData.length;
    }

    /**
	* Construct a new CurveFitter.
	*
	* this time for MULTI DIMENSIONAL FIT
	*/
	public CurveFitter (double[] gridDataX, double[] gridDataY, double[][] zData) {
		this.gridDataX = gridDataX;
		this.gridDataY = gridDataY;
		this.zData = zData;

		// place holders to avoid null pointer exceptions.. o
		xData = gridDataX;
		yData = zData[0];

		min = new Minimisation();
    }

    /**
	* Construct a new CurveFitter.
	*
	* this time for MULTI DIMENSIONAL FIT X#  3  3   3D
	*/
	public CurveFitter (double[] gridDataX, double[] gridDataY, double[] gridDataZ, double[][][] zzData) {
		this.gridDataX = gridDataX;
		this.gridDataY = gridDataY;
		this.gridDataZ = gridDataZ;
		this.zzData = zzData;

		// place holders to avoid null pointer exceptions.. o
		xData = gridDataX;
		zData = zzData[0];
		yData = zData[0];

		min = new Minimisation();
    }


    /**
    *  Set the actual data that we are fitting
    *
    */
    public void setData(double[] xx, double[] yy) {
		xData = xx;
		yData = yy;
	}

    /**
    *  Set a curve to use for the one parameter fit method
    *
    */
    public void setFitData(String fname) {
		fitData = new FitData(fname);
	}


	/**
	*      Perform curve fitting with the simplex method
	*/
    public void doFit(int fitType) {

		fit = fitType;
        //if (fit == HEMISPHERIC) hh=new Hemispheric();
        if (fit == SEMI | fit==SEMI_1) am = new AdiabaticModel();
        if (fit == IBEX) ib = new IBEXLO_09fit();
        //if (fit == IBEX1) ib = new IBEXLO_09fit();
        //if (fit == IBEXFIT) ff = new FitData(fitFile);
		if (fit == IBEXFIT) ff = fitData;
		if (fit == SPINPHASE) iw = new IBEXWind();
		if (fit == COMBO) ibw = new IBEXWind();


		// SIMPLEX INITIALIZTION
		initialize();
		ftol = 1e-9;
		System.out.println("started an optimize");
		// Nelder and Mead minimisation procedure
		min.nelderMead(this, start, step, ftol, 100000);
		System.out.println("completed a Simplex ");
		//min.nelderMead(this, start, step, 1000);
		// get the minimum value
		minimum_pub = min.getMinimum();

		// get values of y and z at minimum
		double[] param = min.getParamValues();
		bestParams = param;


		// Output the results to screen
		for (int i=0; i<param.length; i++) {
			System.out.println("param_" + i + " : " + param[i]);
		}
		if (fit==A3DPOLY) {
			int index =0;
			for (int i=0; i<ord; i++) {
				for (int j=0; j<ord; j++) {
					System.out.println(i+"\t"+j+"\t"+param[index]+"\n");
					index++;
				}
			}
		}
		System.out.println("Minimum = " + min.getMinimum());
        System.out.println("Error: " + getError(param));
        //if (fit == SEMI) System.out.println("sm tests: "+ sm.counter + " "+ sm.counter/25);

		file fil = new file("3dfittestB.txt");
		fil.initWrite(false);
		for (int i=0; i<gridDataX.length; i++) {
			fil.write(gridDataX[i]+"\n");
		}
		fil.write("\n");
		for (int i=0; i<zData[0].length; i++) {
			for (int j=0; j<zData[0].length; j++) {
				fil.write(zData[i][j]+"\n");
			}
		}
		fil.closeWrite();


		fil = new file("3dfittestB2.txt");
		fil.initWrite(false);
		for (int i=0; i<gridDataX.length; i++) {
			fil.write(gridDataX[i]+"\n");
		}
		fil.write("\n");
		for (int i=0; i<zData[0].length; i++) {
			for (int j=0; j<zData[0].length; j++) {
				fil.write(f2(1,param,gridDataX[i],gridDataY[j])+"\n");
			}
		}
		fil.closeWrite();

		/*
        String tbr = "";
        tbr+=("Minimum = " + min.getMinimum() + "\n");
        tbr+=("Error: " + getError(param) + "\n");
		for (int i=0; i<param.length; i++) {
			System.out.println("param_" + i + " : " + param[i]);
			tbr+=("param_" + i + " : " + param[i] + "\n");
		}
        //if (fit == SEMI) System.out.println("sm tests: "+ sm.counter + " "+ sm.counter/25);
        tbr+="xdata \t ydata \t fitdata \n";
		//outF.closeWrite();
		for (int i=0; i<xData.length; i++) {
			tbr+=(xData[i]+"\t"+yData[i]+"\t"+f(fit,param,xData[i])+"\n");
		}
        file f = new file("curvefit_results.txt");
        f.saveShit(tbr);
        */

    }



	/**
	*  Here we cannot check our error and move toward the next local minima in error..
	*
	*   we must compare every possible fit.
	*
	*  For now, only 2d fits here..
	*/
    public void doExhaustiveFit(int fitType, double[] param_limits, int[] steps) {

		fit = fitType;
       // if (fit == HEMISPHERIC) hh=new Hemispheric();
        if (fit == SEMI) sm = new SemiModel();
        if (fit == SEMI_1) sm = new SemiModel();

		int numParams = getNumParams();
		if (numParams!=2)
			throw new IllegalArgumentException("Invalid fit type");

		double[] params = new double[numParams];
		double bestFit = Math.pow(10,100.0);
		bestParams = new double[numParams];
		double[] startParams = new double[numParams];
		double[] deltas = new double[numParams];

		for (int i=0; i<params.length; i++) {

			// set miminimum of space to search
			params[i]=param_limits[2*i];

			// find deltas of space to search
			deltas[i]=(param_limits[2*i+1]-param_limits[2*i])/steps[i];

			// set the best to the first for now
			bestParams[i]=params[i];
			startParams[i]=params[i];
		}

		System.out.println("deltas: " + deltas[0] + " " + deltas[1]);
		System.out.println("steps: " + steps[0] + " " + steps[1]);
		System.out.println("startParams: " + startParams[0] + " " + startParams[1]);
		//start testing them
		for (int i=0; i<steps[0]; i++) {
			for (int j=0; j<steps[1]; j++) {
				params[0]=startParams[0]+i*deltas[0];
				params[1]=startParams[1]+j*deltas[1];
				double test = getError(params);
			//	System.out.println("prms: " + params[0] + " " + params[1]+ " er: " + test);

				if (test<bestFit) {
					//System.out.println("found one: " + test + " i:" + params[0] + " j:" + params[1]);
					bestFit = test;
					bestParams[0] = params[0];
					bestParams[1] = params[1];
				}
			}
		}
	}


    /**
    *
    *Initialise the simplex
    *
    * Here we put starting point for search in parameter space
    */
    void initialize() {
	   start = new double[getNumParams()];
	   step = new double[getNumParams()];
	   java.util.Arrays.fill(step,1.0);
	   java.util.Arrays.fill(start,1.0);


	   double firstx = xData[0];
	   double firsty = yData[0];
	   double lasty = yData[xData.length-1];
	   double lastx = xData[xData.length-1];
	   	   double xmean = (firstx+lastx)/2.0;
	   double ymean = (firsty+lasty)/2.0;
	   double slope;
	   if ((lastx - firstx) != 0.0)
		   slope = (lasty - firsty)/(lastx - firstx);
	   else
		   slope = 1.0;
       double yintercept = firsty - slope * firstx;


       switch (fit) {
            case STRAIGHT_LINE:
                start[0] = yintercept;
                start[1] = slope;
                break;
            case POLY2:
                start[0] = yintercept;
                start[1] = slope;
                start[2] = 0.0;
                break;
            case POLY3:
                start[0] = yintercept;
                start[1] = slope;
                start[2] = 0.0;
                start[3] = 0.0;
                break;
            case POLY4:
                start[0] = yintercept;
                start[1] = slope;
                start[2] = 0.0;
                start[3] = 0.0;
                start[4] = 0.0;
                break;
            case EXPONENTIAL:
                start[0] = 0.1;
                start[1] = 0.01;
                break;
            case POWER:
                start[0] = 0.0;
                start[1] = 1.0;
                break;
            case LOG:
                start[0] = 0.5;
                start[1] = 0.05;
                break;
            case RODBARD:
                start[0] = firsty;
                start[1] = 1.0;
                start[2] = xmean;
                start[3] = lasty;
                break;
            case GAMMA_VARIATE:
                //  First guesses based on following observations:
                //  t0 [b] = time of first rise in gamma curve - so use the user specified first limit
                //  tm = t0 + a*B [c*d] where tm is the time of the peak of the curve
                //  therefore an estimate for a and B is sqrt(tm-t0)
                //  K [a] can now be calculated from these estimates
                start[0] = firstx;
                double ab = xData[getMax(yData)] - firstx;
                start[2] = Math.sqrt(ab);
                start[3] = Math.sqrt(ab);
                start[1] = yData[getMax(yData)] / (Math.pow(ab, start[2]) * Math.exp(-ab/start[3]));
                break;
            case LOG2:
                start[0] = 0.5;
                start[1] = 0.05;
                start[2] = 0.0;
                break;
            case STEP:
            	start[0] = yData[getMax(yData)];
            	start[1] = 2.0;
            	break;
            case HEMISPHERIC:
            	start[0] = 1.0;
            	start[1] = yData[getMax(yData)];
            	break;
            case QUAD_ORIGIN:
            	start[0] = yData[getMax(yData)]/2;
            	start[1] = 1.0;
            	break;
            case GAUSSIAN:
            	start[0] = yData[getMax(yData)];
            	start[1] = xData[getMax(yData)];
            	start[2] = Math.abs((lastx-firstx)/2);
            	System.out.println(""+start[0]+" "+start[1]+" "+start[2]);
            	break;
            case SEMI:
            	start[0] = 1;
            	start[1] = 1.5;
            	//start[2] = 1.5;
            	step[0] = 100;
            	step[1] = 0.1;
            	//step[2] = 1;
            //	min.addConstraint(0,-1,0);
            	min.addConstraint(1,-1,0);
            //	min.addConstraint(1,1,0.001);
            	//min.addConstraint(2,-1,0);
            	//start[2] = 1.5;
            	//step[2] = 1.0;
            	break;
            case SEMI_1:
            	start[0] = 1;
            	step[0] = 100;
            	break;
            case IBEX:
            	// new params are density and temp

            	start[0] = 0.015;
            	start[1] = 6300;
            	step[0] = 0.008;
            	step[1] = 6000;

            	// params are longitude (deg), speed(m/s), density(cm3), temp(k)
            	// old start params for 4 parameter fit..  now we use two parameter fit
            	/*start[0] = 75.0;
            	start[1] = 26000.0;
            	start[2] = 0.015;
            	start[3] = 6300;
            	step[0] = 5;
            	step[1] = 5000;
            	step[2] = 0.003;
            	step[3] = 1000;*/
            	break;
            case IBEXFIT:

            	start[0] = 6000.0;
            	step[0] = 20.0;
            	break;

            case IBEX1:

            	start[0] = 0.0001;
            	step[0] = 0.001;
            	break;

            case SPINPHASE:
            	// params are density, temp, v, theta
            	start[0] = 0.015;
            	start[1] = 5000.0;
            	start[2] = 25000.0;
            	start[3] = 96.5;

            	step[0] = 0.008;
            	step[1] = 3000.0;
            	step[2] = 5000.0;
            	step[3] = 4.0;
            	break;

            case COMBO:
            	// params are density_series, density_sp, temp, v, lamda, theta
            	start[0] = 0.10;
            	start[1] = 0.37;
            	start[2] = 5800.0;
            	start[3] = 30000.0;
            	//start[4] = 73.4;
            	start[4] = 85.0;

            	step[0] = 0.03;
            	step[1] = 0.05;
            	step[2] = 500.0;
            	step[3] = 1000.0;
            	//step[4] = 2.0;
            	step[4] = 1.0;

				// temp constraint
				min.addConstraint(2,-1,3000);
				min.addConstraint(2,+1,12000);

				// v constraint
				min.addConstraint(3,-1,20000);
				min.addConstraint(3,+1,30000);

				// norm constraints
				min.addConstraint(0,-1,0);
				min.addConstraint(1,-1,0);

				// longitude constraint
				//min.addConstraint(4,-1,60);
				//min.addConstraint(4,+1,85);


            	break;

            case ONE_PARAM_FIT:
            	start[0] = 1.0;
            	break;

            case A2DPOLY:
            	int sign = 1;
            	for (int i=0; i<getNumParams(); i++) {
					start[i]=0.0;
					step[i]=-1.0*sign;
					sign*=-1;
				}
				start[0]=10.0;

            case A3DPOLY:
            	int sign3 = 1;
            	for (int i=0; i<getNumParams(); i++) {
					start[i]=0.0;
					step[i]=-1.0*sign3;
					sign3*=-1;
				}
				start[0]=10.0;

        }
    }


    /** Get number of parameters for current fit function */
    public int getNumParams() {
        switch (fit) {
            case STRAIGHT_LINE: return 2;
            case POLY2: return 3;
            case POLY3: return 4;
            case POLY4: return 5;
            case EXPONENTIAL: return 2;
            case POWER: return 2;
            case LOG: return 2;
            case RODBARD: return 4;
            case GAMMA_VARIATE: return 4;
            case LOG2: return 3;
            case STEP: return 2;
            case HEMISPHERIC: return 2;
            case SEMI: return 2;
            case QUAD_ORIGIN: return 2;
            case GAUSSIAN: return 3;
            case SEMI_1: return 1;
            case IBEX: return 2;
            case IBEXFIT: return 1;
            case IBEX1: return 1;
            case SPINPHASE: return 4;
            case COMBO: return 5;
            case ONE_PARAM_FIT: return 1;
            case A2DPOLY: return ord*ord;
            case A3DPOLY: return ord*ord*ord;
        }
        return 0;
    }

    /**
    *Returns "fit" function value for parametres "p" at "x,y"
    *
    *  Define 2D function to fit to here!!
    */
    public double f2(int fit, double[] p, double xx, double yy) {
		//if (fit==A3DPOLY) {
			// first lets get the params into matrix form
			int index = 0;
			double[][] paraMatrix = new double[ord][ord];
			for (int i=0; i<ord; i++) {
				for (int j=0; j<ord; j++) {
					paraMatrix[i][j] = p[index];
					index++;
				}
			}
			//System.out.println("final index making array: " + index);

			// now calculate the value of our polynomial at the point x[][]
			double tbr = 0.0;
			double xTemp = 1.0;
			for (int i=0; i<ord; i++) {
				double yTemp = 1.0;
				for (int j=0; j<ord; j++) {
					tbr += paraMatrix[i][j]*xTemp*yTemp;
					yTemp*=yy;
				}
				xTemp*=xx;
			}

			// if this were a single power series:
			/*tbr = 0.0
			tempx = x;
			for (int i=0; i<5; i++) {
				tbr += A[i]*x;
				x*=x;
			}
			*/
			return tbr;
		//}
		//else return Double.MAX_VALUE;
	}

	/**
	*Returns "fit" function value for parametres "p" at "x,y"
	*
	*  Define 3D function to fit to here!!
	*/
	public double f3(int fit, double[] p, double xx, double yy, double zz) {
		//if (fit==A3DPOLY) {
			// first lets get the params into matrix form
			int index = 0;
			double[][][] paraMatrix = new double[ord][ord][ord];
			for (int i=0; i<ord; i++) {
				for (int j=0; j<ord; j++) {
					for (int k=0; k<ord; k++) {
						paraMatrix[i][j][k] = p[index];
						index++;
					}
				}
			}
			//System.out.println("final index making array: " + index);

			// now calculate the value of our polynomial at the point x[][]
			double tbr = 0.0;
			double xTemp = 1.0;
			for (int i=0; i<ord; i++) {
				double yTemp = 1.0;
				for (int j=0; j<ord; j++) {
					double zTemp = 1.0;
					for (int k=0; k<ord; k++) {
						tbr += paraMatrix[i][j][k]*xTemp*yTemp*zTemp;
						zTemp*=zz;
					}
					yTemp*=yy;
				}
				xTemp*=xx;
			}

			// if this were a single power series:
			/*tbr = 0.0
			tempx = x;
			for (int i=0; i<5; i++) {
				tbr += A[i]*x;
				x*=x;
			}
			*/
			return tbr;
		//}
		//else return Double.MAX_VALUE;
	}


    /**
    *Returns "fit" function value for parametres "p" at "x"
    *
    *  Define function to fit to here!!
    */
    public double f(int fit, double[] p, double x) {
        switch (fit) {
            case STRAIGHT_LINE:
                return p[0] + p[1]*x;
            case POLY2:
                return p[0] + p[1]*x + p[2]* x*x;
            case POLY3:
                return p[0] + p[1]*x + p[2]*x*x + p[3]*x*x*x;
            case POLY4:
                return p[0] + p[1]*x + p[2]*x*x + p[3]*x*x*x + p[4]*x*x*x*x;
            case EXPONENTIAL:
                return p[0]*Math.exp(p[1]*x);
            case POWER:
                if (x == 0.0)
                    return 0.0;
                else
                    return p[0]*Math.exp(p[1]*Math.log(x)); //y=ax^b
            case LOG:
                if (x == 0.0)
                    x = 0.5;
                return p[0]*Math.log(p[1]*x);
            case RODBARD:
                double ex;
                if (x == 0.0)
                    ex = 0.0;
                else
                    ex = Math.exp(Math.log(x/p[2])*p[1]);
                double y = p[0]-p[3];
                y = y/(1.0+ex);
                return y+p[3];
            case GAMMA_VARIATE:

                if (p[0] >= x) return 0.0;
                if (p[1] <= 0) return -100000.0;
                if (p[2] <= 0) return -100000.0;
                if (p[3] <= 0) return -100000.0;

                double pw = Math.pow((x - p[0]), p[2]);
                double e = Math.exp((-(x - p[0]))/p[3]);
                return p[1]*pw*e;
            case LOG2:
                double tmp = x-p[2];
                if (tmp<0.001) tmp = 0.001;
                return p[0]+p[1]*Math.log(tmp);
            case STEP:
            	if (x>p[1]) return 0.0;
            	else return p[0];
            case HEMISPHERIC:
            	//return hh.eflux(p[0],p[1],x);
            	return 0.0;
            case SEMI:
            	return am.f(p[0],p[1],x);
            case SEMI_1:
            	return am.f(p[0],gamma,x);
            	//return sm.f(p[0],p[1],p[2],x);
            case QUAD_ORIGIN:
            	return p[0]*x + p[1]*x*x;
            case GAUSSIAN:
            	return p[0]/p[2]/Math.sqrt(Math.PI*2)*Math.exp(-(x-p[1])*(x-p[1])/p[2]/p[2]/2);
            case IBEX:
            	//return ib.getRate(p[0],p[1],x);
            	return 0.0;
            case IBEX1:
            	return 0.0;
            	//return ib.getRate(p[0],x);
            case IBEXFIT:
            	return p[0]*ff.getRate(x);
            case SPINPHASE:
            	//density, temperature, v, theta, x (sp)
            	return ib.getRate(p[0], p[1], p[2], p[3],x);
            case COMBO:
            	//return ibw.getRateMulti(p[0],p[1],p[2],p[3],p[4],p[5],x);
            	return ibw.getRateMulti(p[0],p[1],p[2],p[3],74.68,p[4],x);
            case ONE_PARAM_FIT:
            	return p[0]*fitData.getRate(x);
            default:
                return 0.0;
        }
    }


	file f_log = new file("curve_fit_log_9.txt");
	public int num = 0;

	/**
	* Returns sum of squares of residuals
	*   er huh?  returns the sume of the square errors at each point  for least square fittig
	*/
	public double getError(double[] params_) {
		num++;
		//f_log.initWrite(true);
		double tbr = 0.0;
		if (fit!=A3DPOLY) {
			for (int i = 0; i < xData.length; i++) {
				tbr += sqr(yData[i] - f(fit, params_, xData[i]));
			}
		}
		else if (fit==A2DPOLY) {  // looking for a square here
			for (int i=0; i<gridDataX.length; i++) {
				for (int j=0; j<gridDataY.length; j++) {
					tbr += sqr(zData[i][j] - f2(fit, params_, gridDataX[i], gridDataY[j]));
				}
			}
		}
		else if (fit==A3DPOLY) { // looking for a cube here
			for (int i=0; i<gridDataX.length; i++) {
				for (int j=0; j<gridDataY.length; j++) {
					for (int k=0; k<gridDataZ.length; k++) {
						tbr += sqr(zzData[i][j][k] - f3(fit, params_, gridDataX[i], gridDataY[j], gridDataZ[k]));
					}
				}
			}
		}

		//for (int i=0; i<params_.length; i++) {
		//	f_log.write(params_[i]+"\t" );
		//}
		//System.out.println("made it here : " + num);
		//f_log.write(tbr + "\n");
		//f_log.closeWrite();
		if (num%1000==0) System.out.println("num: " + num + " error: " + tbr +"\n");
		return tbr;
    }

    /** Here's the one to minimize!!  */
    public double function(double[] params) {
		return getError(params);
	}

    public static double sqr(double d) { return d * d; }


    /**
	* Gets index of highest value in an array.
	*
	* @param              Double array.
	* @return             Index of highest value.
	*/
	public static int getMax(double[] array) {
		double max = array[0];
		int index = 0;
		for(int i = 1; i < array.length; i++) {
			if(max < array[i]) {
				max = array[i];
				index = i;
			}
		}
		return index;
    }


    /**
    *  Use this to fit a curve or
    * for testing...
    *
    *  Reads a file for data and fits to a curve at command line
    *
    */
    public static final void main(String[] args) {

		int numLines = 62;
		int linesToSkip = 0;


		if (args.length==1) {
			//double[] x = {0,1,2,3,4,5,6,7,8,9};
			//double[] y = {1,3,5,7,9,11,13,15,17,19};
			//CurveFitter cf = new CurveFitter(x,y);
			//cf.doFit(CurveFitter.STRAIGHT_LINE);

			file f = new file(args[0]);
			f.initRead();
			double[] x = new double[numLines];
			double[] y = new double[numLines];
			String line = "";

			for (int i=0; i<linesToSkip; i++) {
				line = f.readLine();
			}

			for (int i=0; i<numLines; i++) {
				line = f.readLine();
				StringTokenizer st = new StringTokenizer(line);
				x[i]=Double.parseDouble(st.nextToken());
				//x[i]-=1.0; // IMPORTANT _ SHIFT BY 1.0 TO MOVE TO HELIOCENTRIC FRAME!!
				y[i]=Double.parseDouble(st.nextToken());
				//System.out.println("row: " + i + " x: " + x[i] + " y: " + y[i]);
			}
			f.closeRead();

			CurveFitter cf = new CurveFitter(x,y);
			System.out.println("starting doFit routine");
			cf.doFit(CurveFitter.COMBO);
		}

		//if (args.length==2) {
		else if (args.length==2) {
			file f = new file(args[0]);
			f.initRead();
			double[] x = new double[numLines];
			double[] y = new double[numLines];
			String line = "";
			for (int i=0; i<linesToSkip; i++) {
				line = f.readLine();
			}

			for (int i=0; i<numLines; i++) {
				line = f.readLine();
				StringTokenizer st = new StringTokenizer(line);
				x[i]=Double.parseDouble(st.nextToken());
				x[i]-=1.0; // IMPORTANT _ SHIFT BY 1.0 TO MOVE TO HELIOCENTRIC FRAME!!
				y[i]=Double.parseDouble(st.nextToken());
				System.out.println("row: " + i + " x: " + x[i] + " y: " + y[i]);
			}
			f.closeRead();

			for (double gamma =0.8; gamma<2.2; gamma+=0.1) {
				CurveFitter cf = new CurveFitter(x,y);
				cf.gamma = gamma;
				System.out.println("gamma: " + gamma);
				cf.doFit(CurveFitter.SEMI_1);
				//ystem.out.println(
			}
		}
		// testing here!!
		else if (args.length==0) {
			int size = 10;
			double[] x = new double[size];
			double[] y = new double[size];
			double[] z = new double[size];
			double[][][] zz = new double[size][size][size];
			for (int i=0; i<size; i++) {
				x[i]=i*20.0/size-10.0;
				y[i]=i*20.0/size-10.0;
				z[i]=i*20.0/size-10.0;
			}
			for (int i=0; i<size; i++) {
				for (int j=0; j<size; j++) {
					for (int k=0; k<size; k++) {
						zz[i][j][k]=10.0*Math.exp( -(x[i]*x[i] + y[j]*y[j] + z[k]*z[k]) );
					}
				}
			}
			CurveFitter cf = new CurveFitter(x,y,z,zz);
			System.out.println("Starting 3D fit test... ");
			cf.doFit(CurveFitter.A3DPOLY);
		}


		//double[] lims = {1E6,1E10, 0.00000001, 0.0001};
		//int[] steps = {4,4};
		//cf.doExhaustiveFit(CurveFitter.SEMI,lims,steps);
		//System.out.println("param[0]: " + cf.bestParams[0]);
		//System.out.println("param[1]: " + cf.bestParams[1]);


		//cf.setRestarts(100);
		//Date d1 = new Date();
		//cf.doFit(type);
		//Date d2 = new Date();

//		cf.print();

//		System.out.println(cf.getResultString()+"\n\n");
		//System.out.println("param[0]: " + cf.getParams()[0]);
		//System.out.println("param[1]: " + cf.getParams()[1]);
		//System.out.println("took: " + (d2.getTime()-d1.getTime()));


		//int max = getMax(y);
		//System.out.println("y_max: " + y[max]);
		//double[] lims = {0.0 , 10.0, y[max]/2.0, 4*y[max]};
		//System.out.println("lims: " + lims[0] + " " + lims[1] + " " + lims[2] + " " + lims[3]);
		//int[] steps = {256,128};
		//Date d3 = new Date();
		//cf.doExhaustiveFit(CurveFitter.QUAD_ORIGIN,lims,steps);
		//cf.doFit(type);
		//Date d4 = new Date();
		//double[] ans = cf.bestParams;
		//System.out.println("a[0]: " + ans[0]);
		//System.out.println("a[1]: " + ans[1]);
		//System.out.println(cf.getResultString());
		//System.out.println("took: " + (d4.getTime()-d3.getTime()));

	}

/*
	public double getSPError(double[] spModel) {
		file f = new file("2010_fit_sp.txt");
		int numLines = f.readShitNumLines();
		f.initRead();
		double[] x = new double[numLines];
		double[] y = new double[numLines];
		String line = "";

		for (int i=0; i<linesToSkip; i++) {
			line = f.readLine();
		}

		for (int i=0; i<numLines; i++) {
			line = f.readLine();
			StringTokenizer st = new StringTokenizer(line);
			x[i]=Double.parseDouble(st.nextToken());
			y[i]=Double.parseDouble(st.nextToken());
		}
		f.closeRead();

	}

	public double getTimeError(double[] timeModel) {
		file f = new file("2010_fit_time.txt");
		int numLines = f.readShitNumLines();
		f.initRead();
		double[] x = new double[numLines];
		double[] y = new double[numLines];
		String line = "";

		for (int i=0; i<linesToSkip; i++) {
			line = f.readLine();
		}

		for (int i=0; i<numLines; i++) {
			line = f.readLine();
			StringTokenizer st = new StringTokenizer(line);
			x[i]=Double.parseDouble(st.nextToken());
			//x[i]-=1.0; // IMPORTANT _ SHIFT BY 1.0 TO MOVE TO HELIOCENTRIC FRAME!!
			y[i]=Double.parseDouble(st.nextToken());
			//System.out.println("row: " + i + " x: " + x[i] + " y: " + y[i]);
		}
		f.closeRead();

	}

	public double getNormalizedError(double[] timeModel, double[] spModel) {
		// first load the fit set to compare the model to
		file f = new file("2010_fit_set.txt");
		f.initRead();
		double[] x = new double[numLines];
		double[] y = new double[numLines];
		String line = "";

		for (int i=0; i<linesToSkip; i++) {
			line = f.readLine();
		}

		for (int i=0; i<numLines; i++) {
			line = f.readLine();
			StringTokenizer st = new StringTokenizer(line);
			x[i]=Double.parseDouble(st.nextToken());
			//x[i]-=1.0; // IMPORTANT _ SHIFT BY 1.0 TO MOVE TO HELIOCENTRIC FRAME!!
			y[i]=Double.parseDouble(st.nextToken());
			//System.out.println("row: " + i + " x: " + x[i] + " y: " + y[i]);
		}
		f.closeRead();

		// ok we have the data to fit our result to, and the data itself is passed in here..

*/





	/*
	* This is to load a given model from a file and give a proper normalization to match the ibex data
	* for using a one parameter (normalization) fit.
	*/
	public class FitData {
		public StringTokenizer st;
		public double[] days;
		public double[] rates;
		public Vector daysV;
		public Vector ratesV;  // we call this days and rates, but could be any x and y
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
			for (int i=0; i<days.length-1; i++) {
				if (day<days[i]) { // this is where we want to be
					return (rates[i]+rates[i+1])/2;
				}
			}
			return 0;
		}
	}


}

