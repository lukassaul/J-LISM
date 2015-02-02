

//import ij.*;
//import ij.gui.*;
import java.util.*;
import flanagan.math.*;

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
*/
public class CurveFitter {
    public static final int STRAIGHT_LINE=0,POLY2=1,POLY3=2,POLY4=3,
    EXPONENTIAL=4,POWER=5,LOG=6,RODBARD=7,GAMMA_VARIATE=8,
    LOG2=9, STEP=10, HEMISPHERIC=11, HTEST=12, QUAD_ORIGIN=13, GAUSSIAN=14, SEMI=15;

    public static final String[] fitList = {"Straight Line","2nd Degree Polynomial",
    	"3rd Degree Polynomial", "4th Degree Polynomial","Exponential","Power",
    	"log","Rodbard", "Gamma Variate", "y = a+b*ln(x-c)", "Step Function",
    	"Hemispheric", "H-test", "Quadratic Through Origin","Gaussian","Semi"};

    public static final String[] fList = {"y = a+bx","y = a+bx+cx^2",
    	"y = a+bx+cx^2+dx^3", "y = a+bx+cx^2+dx^3+ex^4","y = a*exp(bx)","y = ax^b",
    	"y = a*ln(bx)", "y = d+(a-d)/(1+(x/c)^b)", "y = a*(x-b)^c*exp(-(x-b)/d)",
    	"y = a+b*ln(x-c)", "y = a*step(x-b)", "y=hem.f(a,x)",
    	"y = norm*step()..","y=ax+bx^2", "y=a*EXP(-(x-b)^2/c)","y=semi.f(a,b,x)"};

    private static final double root2 = 1.414214; // square root of 2

    private int fit;                // Number of curve type to fit
    private double[] xData, yData;  // x,y data to fit

	public double[] bestParams; // take this after doing exhaustive for the answer
	public Hemispheric hh;
	public SemiModel sm;

	/**
	* build it with floats and cast
	*/
    public CurveFitter (float[] xData, float[] ydata) {
		this.xData = new double[xData.length];
		this.yData = new double[xData.length];
		for (int i=0; i<xData.length; i++) {
			this.xData[i]=(float)xData[i];
			this.yData[i]=(float)yData[i];
		}
		numPoints = xData.length;
	}


   /** Construct a new CurveFitter. */
    public CurveFitter (double[] xData, double[] yData) {
        this.xData = xData;
        this.yData = yData;
        numPoints = xData.length;
    }

	 /**
	 *      Perform curve fitting with the simplex method
	 */
    public void doFit(int fitType) {
        doFit(fitType, false);
    }

	/**
	*  Here we cannot check our error and move toward the next local minima in error..
	*
	*   we must compare every possible fit.
	*
	*  For now, only 2d fits here..
	*/
    public void doExhaustiveFit(int fitType, double[] param_limits, int[] steps) {
        if (fitType < STRAIGHT_LINE || fitType > GAUSSIAN)
            throw new IllegalArgumentException("Invalid fit type");

		fit = fitType;
        if (fit == HEMISPHERIC) hh=new Hemispheric();
        if (fit == SEMI) semi = new SemiModel();

		numParams = getNumParams();
		if (numParams*2!=param_limits.length)
			throw new IllegalArgumentException("Invalid fit type");

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



    public void doFit(int fitType, boolean showSettings) {
        if (fitType < STRAIGHT_LINE || fitType > GAUSSIAN)
            throw new IllegalArgumentException("Invalid fit type");
        fit = fitType;

        if (fit == HEMISPHERIC) hh=new Hemispheric();

        initialize();
        //if (showSettings) settingsDialog();
        restart(0);

        numIter = 0;
        boolean done = false;
        double[] center = new double[numParams];  // mean of simplex vertices
        while (!done) {
            numIter++;
            for (int i = 0; i < numParams; i++) center[i] = 0.0;
            // get mean "center" of vertices, excluding worst
            for (int i = 0; i < numVertices; i++)
                if (i != worst)
                    for (int j = 0; j < numParams; j++)
                        center[j] += simp[i][j];
            // Reflect worst vertex through centre
            for (int i = 0; i < numParams; i++) {
                center[i] /= numParams;
                next[i] = center[i] + alpha*(simp[worst][i] - center[i]);
            }
            sumResiduals(next);
            // if it's better than the best...
            if (next[numParams] <= simp[best][numParams]) {
                newVertex();
                // try expanding it
                for (int i = 0; i < numParams; i++)
                    next[i] = center[i] + gamma * (simp[worst][i] - center[i]);
                sumResiduals(next);
                // if this is even better, keep it
                if (next[numParams] <= simp[worst][numParams])
                    newVertex();
            }
            // else if better than the 2nd worst keep it...
            else if (next[numParams] <= simp[nextWorst][numParams]) {
                newVertex();
            }
            // else try to make positive contraction of the worst
            else {
                for (int i = 0; i < numParams; i++)
                    next[i] = center[i] + beta*(simp[worst][i] - center[i]);
                sumResiduals(next);
                // if this is better than the second worst, keep it.
                if (next[numParams] <= simp[nextWorst][numParams]) {
                    newVertex();
                }
                // if all else fails, contract simplex in on best
                else {
                    for (int i = 0; i < numVertices; i++) {
                        if (i != best) {
                            for (int j = 0; j < numVertices; j++)
                                simp[i][j] = beta*(simp[i][j]+simp[best][j]);
                            sumResiduals(simp[i]);
                        }
                    }
                }
            }
            order();

            double rtol = 2 * Math.abs(simp[best][numParams] - simp[worst][numParams]) /
            (Math.abs(simp[best][numParams]) + Math.abs(simp[worst][numParams]) + 0.0000000001);

            if (numIter >= maxIter) done = true;
            else if (rtol < maxError) {
                System.out.print(getResultString());
                restarts--;
                if (restarts < 0) {
                    done = true;
                }
                else {
                    restart(best);
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
        // Calculate some things that might be useful for predicting parametres
        numParams = getNumParams();
        numVertices = numParams + 1;      // need 1 more vertice than parametres,
        simp = new double[numVertices][numVertices];
        next = new double[numVertices];

        double firstx = xData[0];
        double firsty = yData[0];
        double lastx = xData[numPoints-1];
        double lasty = yData[numPoints-1];
        double xmean = (firstx+lastx)/2.0;
        double ymean = (firsty+lasty)/2.0;
        double slope;
        if ((lastx - firstx) != 0.0)
            slope = (lasty - firsty)/(lastx - firstx);
        else
            slope = 1.0;
        double yintercept = firsty - slope * firstx;
        maxIter = IterFactor * numParams * numParams;  // Where does this estimate come from?
        restarts = 1;
        maxError = 1e-9;
        switch (fit) {
            case STRAIGHT_LINE:
                simp[0][0] = yintercept;
                simp[0][1] = slope;
                break;
            case POLY2:
                simp[0][0] = yintercept;
                simp[0][1] = slope;
                simp[0][2] = 0.0;
                break;
            case POLY3:
                simp[0][0] = yintercept;
                simp[0][1] = slope;
                simp[0][2] = 0.0;
                simp[0][3] = 0.0;
                break;
            case POLY4:
                simp[0][0] = yintercept;
                simp[0][1] = slope;
                simp[0][2] = 0.0;
                simp[0][3] = 0.0;
                simp[0][4] = 0.0;
                break;
            case EXPONENTIAL:
                simp[0][0] = 0.1;
                simp[0][1] = 0.01;
                break;
            case POWER:
                simp[0][0] = 0.0;
                simp[0][1] = 1.0;
                break;
            case LOG:
                simp[0][0] = 0.5;
                simp[0][1] = 0.05;
                break;
            case RODBARD:
                simp[0][0] = firsty;
                simp[0][1] = 1.0;
                simp[0][2] = xmean;
                simp[0][3] = lasty;
                break;
            case GAMMA_VARIATE:
                //  First guesses based on following observations:
                //  t0 [b] = time of first rise in gamma curve - so use the user specified first limit
                //  tm = t0 + a*B [c*d] where tm is the time of the peak of the curve
                //  therefore an estimate for a and B is sqrt(tm-t0)
                //  K [a] can now be calculated from these estimates
                simp[0][0] = firstx;
                double ab = xData[getMax(yData)] - firstx;
                simp[0][2] = Math.sqrt(ab);
                simp[0][3] = Math.sqrt(ab);
                simp[0][1] = yData[getMax(yData)] / (Math.pow(ab, simp[0][2]) * Math.exp(-ab/simp[0][3]));
                break;
            case LOG2:
                simp[0][0] = 0.5;
                simp[0][1] = 0.05;
                simp[0][2] = 0.0;
                break;
            case STEP:
            	simp[0][0] = yData[getMax(yData)];
            	simp[0][1] = 2.0;
            	break;
            case HEMISPHERIC:
            	simp[0][0] = 1.0;
            	simp[0][1] = yData[getMax(yData)];
            	break;
            case QUAD_ORIGIN:
            	simp[0][0] = yData[getMax(yData)]/2;
            	simp[0][1] = 1.0;
            	break;
            case GAUSSIAN:
            	simp[0][0] = yData[getMax(yData)];
            	simp[0][1] = xData[getMax(yData)];
            	simp[0][2] = Math.abs((lastx-firstx)/2);
            	System.out.println(""+simp[0][0]+" "+simp[0][1]+" "+simp[0][2]);
            	break;
            case SEMI:
            	simp[0][0] = 1.0;
            	simp[0][2] = y.Data[getMax(yData)];
            	break;
        }
    }

    /** Restart the simplex at the nth vertex */
    void restart(int n) {
        // Copy nth vertice of simplex to first vertice
        for (int i = 0; i < numParams; i++) {
            simp[0][i] = simp[n][i];
        }
        sumResiduals(simp[0]);          // Get sum of residuals^2 for first vertex
        double[] step = new double[numParams];
        for (int i = 0; i < numParams; i++) {
            step[i] = simp[0][i] / 2.0;     // Step half the parametre value
            if (step[i] == 0.0)             // We can't have them all the same or we're going nowhere
                step[i] = 0.01;
        }
        // Some kind of factor for generating new vertices
        double[] p = new double[numParams];
        double[] q = new double[numParams];
        for (int i = 0; i < numParams; i++) {
            p[i] = step[i] * (Math.sqrt(numVertices) + numParams - 1.0)/(numParams * root2);
            q[i] = step[i] * (Math.sqrt(numVertices) - 1.0)/(numParams * root2);
        }
        // Create the other simplex vertices by modifing previous one.
        for (int i = 1; i < numVertices; i++) {
            for (int j = 0; j < numParams; j++) {
                simp[i][j] = simp[i-1][j] + q[j];
            }
            simp[i][i-1] = simp[i][i-1] + p[i-1];
            sumResiduals(simp[i]);
        }
        // Initialise current lowest/highest parametre estimates to simplex 1
        best = 0;
        worst = 0;
        nextWorst = 0;
        order();
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
            case SEMI return 2;
            case QUAD_ORIGIN: return 2;
            case GAUSSIAN: return 3;
        }
        return 0;
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
            	return hh.eflux(p[0],p[1],x);
            case SEMI:
            	return sm.f(p[0],p[1],x);
            case QUAD_ORIGIN:
            	return p[0]*x + p[1]*x*x;
            case GAUSSIAN:
            	return p[0]/p[2]/Math.sqrt(Math.PI*2)*Math.exp(-(x-p[1])*(x-p[1])/p[2]/p[2]/2);
            default:
                return 0.0;
        }
    }

    /** Get the set of parameter values from the best corner of the simplex */
    public double[] getParams() {
        order();
        return simp[best];
    }

    /** Returns residuals array ie. differences between data and curve. */
    public double[] getResiduals() {
        double[] params = getParams();
        double[] residuals = new double[numPoints];
        for (int i = 0; i < numPoints; i++)
            residuals[i] = yData[i] - f(fit, params, xData[i]);
        return residuals;
    }



	/** Returns sum of squares of residuals */
	public double getError(double[] params_) {
		double tbr = 0.0;
		for (int i = 0; i < numPoints; i++) {
			tbr += sqr(yData[i] - f(fit, params_, xData[i]));
		}
		//System.out.println("error: " + tbr);
		return tbr;
    }

    /* Last "parametre" at each vertex of simplex is sum of residuals
     * for the curve described by that vertex
     */
    public double getSumResidualsSqr() {
        double sumResidualsSqr = (getParams())[getNumParams()];
        return sumResidualsSqr;
    }



    /**  SD = sqrt(sum of residuals squared / number of params+1)
     */
    public double getSD() {
        double sd = Math.sqrt(getSumResidualsSqr() / numVertices);
        return sd;
    }

    /**  Get a measure of "goodness of fit" where 1.0 is best.
     *
     */
    public double getFitGoodness() {
        double sumY = 0.0;
        for (int i = 0; i < numPoints; i++) sumY += yData[i];
        double mean = sumY / numVertices;
        double sumMeanDiffSqr = 0.0;
        int degreesOfFreedom = numPoints - getNumParams();
        double fitGoodness = 0.0;
        for (int i = 0; i < numPoints; i++) {
            sumMeanDiffSqr += sqr(yData[i] - mean);
        }
        if (sumMeanDiffSqr > 0.0 && degreesOfFreedom != 0)
            fitGoodness = 1.0 - (getSumResidualsSqr() / degreesOfFreedom) * ((numParams) / sumMeanDiffSqr);

        return fitGoodness;
    }

    /** Get a string description of the curve fitting results
     * for easy output.
     */
    public String getResultString() {
        StringBuffer results = new StringBuffer("\nNumber of iterations: " + getIterations() +
        "\nMaximum number of iterations: " + getMaxIterations() +
        "\nSum of residuals squared: " + getSumResidualsSqr() +
        "\nStandard deviation: " + getSD() +
        "\nGoodness of fit: " + getFitGoodness() +
        "\nParameters:");
        char pChar = 'a';
        double[] pVal = getParams();
        for (int i = 0; i < numParams; i++) {
            results.append("\n" + pChar + " = " + pVal[i]);
            pChar++;
        }
        return results.toString();
    }

    public static double sqr(double d) { return d * d; }

    /** Adds sum of square of residuals to end of array of parameters */
    void sumResiduals (double[] x) {
        x[numParams] = 0.0;
        for (int i = 0; i < numPoints; i++) {
            x[numParams] = x[numParams] + sqr(f(fit,x,xData[i])-yData[i]);
            //        if (IJ.debugMode) ij.IJ.log(i+" "+x[n-1]+" "+f(fit,x,xData[i])+" "+yData[i]);
        }
    }

    /** Keep the "next" vertex */
    void newVertex() {
        for (int i = 0; i < numVertices; i++)
            simp[worst][i] = next[i];
    }

    /** Find the worst, nextWorst and best current set of parameter estimates */
    void order() {
        for (int i = 0; i < numVertices; i++) {
            if (simp[i][numParams] < simp[best][numParams]) best = i;
            if (simp[i][numParams] > simp[worst][numParams]) worst = i;
        }
        nextWorst = best;
        for (int i = 0; i < numVertices; i++) {
            if (i != worst) {
                if (simp[i][numParams] > simp[nextWorst][numParams]) nextWorst = i;
            }
        }
        //        IJ.write("B: " + simp[best][numParams] + " 2ndW: " + simp[nextWorst][numParams] + " W: " + simp[worst][numParams]);
    }

    /** Get number of iterations performed */
    public int getIterations() {
        return numIter;
    }

    /** Get maximum number of iterations allowed */
    public int getMaxIterations() {
        return maxIter;
    }

    /** Set maximum number of iterations allowed */
    public void setMaxIterations(int x) {
        maxIter = x;
    }

    /** Get number of simplex restarts to do */
    public int getRestarts() {
        return restarts;
    }

    /** Set number of simplex restarts to do */
    public void setRestarts(int x) {
        restarts = x;
    }

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
		int numLines = 17;
		int linesToSkip = 0;
		int type = CurveFitter.QUAD_ORIGIN;

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

		CurveFitter cf = new CurveFitter(x,y);
		/*cf.setRestarts(100);
		Date d1 = new Date();
		cf.doFit(type);
		Date d2 = new Date();

		System.out.println(cf.getResultString()+"\n\n");
		System.out.println("param[0]: " + cf.getParams()[0]);
		System.out.println("param[1]: " + cf.getParams()[1]);
		System.out.println("took: " + (d2.getTime()-d1.getTime()));
		*/

		//int max = getMax(y);
		//System.out.println("y_max: " + y[max]);
		//double[] lims = {0.0 , 10.0, y[max]/2.0, 4*y[max]};
		//System.out.println("lims: " + lims[0] + " " + lims[1] + " " + lims[2] + " " + lims[3]);
		//int[] steps = {256,128};
		//Date d3 = new Date();
		//cf.doExhaustiveFit(CurveFitter.QUAD_ORIGIN,lims,steps);
		cf.doFit(type);
		//Date d4 = new Date();
		//double[] ans = cf.bestParams;
		//System.out.println("a[0]: " + ans[0]);
		//System.out.println("a[1]: " + ans[1]);
		System.out.println(cf.getResultString());
		//System.out.println("took: " + (d4.getTime()-d3.getTime()));

	}

}

