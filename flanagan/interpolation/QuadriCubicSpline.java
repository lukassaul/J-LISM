/**********************************************************
*
*   QuadriCubicSpline.java
*
*   Class for performing an interpolation on the tabulated
*   function y = f(x1,x2,x3,x4) using a natural quadricubic spline
*   Assumes second derivatives at end points = 0 (natural spine)
*
*   WRITTEN BY: Dr Michael Thomas Flanagan
*
*   DATE:	May 2003
*   UPDATE  July 2007, 4 December 2007
*
*   DOCUMENTATION:
*   See Michael Thomas Flanagan's Java library on-line web page:
*   http://www.ee.ucl.ac.uk/~mflanaga/java/QuadriCubicSpline.html
*   http://www.ee.ucl.ac.uk/~mflanaga/java/
*
*   Copyright (c) May 2003, December 2007   Michael Thomas Flanagan
*
*   PERMISSION TO COPY:
*   Permission to use, copy and modify this software and its documentation for
*   NON-COMMERCIAL purposes is granted, without fee, provided that an acknowledgement
*   to the author, Michael Thomas Flanagan at http://www.ee.ucl.ac.uk/~mflanaga, appears in all copies.
*
*   Dr Michael Thomas Flanagan makes no representations about the suitability
*   or fitness of the software for any or for a particular purpose.
*   Michael Thomas Flanagan shall not be liable for any damages suffered
*   as a result of using, modifying or distributing this software or its derivatives.
*
***************************************************************************************/

package flanagan.interpolation;

import flanagan.math.Fmath;

public class QuadriCubicSpline{

    	private int nPoints = 0;   	    // no. of x1 tabulated points
    	private int mPoints = 0;   	    // no. of x2 tabulated points
    	private int lPoints = 0;   	    // no. of x3 tabulated points
    	private int kPoints = 0;   	    // no. of x4 tabulated points

    	private double[][][][] y = null;    // y=f(x1,x2,x3,x4) tabulated function
    	private double[] x1 = null;   	    // x1 in tabulated function f(x1,x2,x3,x4)
    	private double[] x2 = null;   	    // x2 in tabulated function f(x1,x2,x3,x4)
    	private double[] x3 = null;   	    // x3 in tabulated function f(x1,x2,x3,x4)
    	private double[] x4 = null;   	    // x4 in tabulated function f(x1,x2,x3,x4)
    	private double[] xMin = new double[4];   // minimum values of x1, x2, x3 and x4
    	private double[] xMax = new double[4];   // maximum values of x1, x2, x3 and x4

    	private TriCubicSpline[] tcsn = null;   // nPoints array of TriCubicSpline instances
    	private CubicSpline csm = null;         // CubicSpline instance
    	private double[][][][] d2ydx2 = null;   // inner matrix of second derivatives
        private boolean derivCalculated = false;    // = true when the called triicubic spline derivatives have been calculated
        private String subMatrixIndices = " ";      // String of indices of the submatrices that have called QuadriCubicSpline if called from higher dimension interpolation
        private boolean averageIdenticalAbscissae = false;  // if true: the the ordinate values for identical abscissae are averaged
                                                            // If false: the abscissae values are separated by 0.001 of the total abscissae range;


    	// Constructor
    	public QuadriCubicSpline(double[] x1, double[] x2, double[] x3, double[] x4, double[][][][] y){
        	this.nPoints=x1.length;
        	this.mPoints=x2.length;
        	this.lPoints=x3.length;
        	this.kPoints=x4.length;
        	if(this.nPoints!=y.length)throw new IllegalArgumentException("Arrays x1 and y-row are of different length " + this.nPoints + " " + y.length);
        	if(this.mPoints!=y[0].length)throw new IllegalArgumentException("Arrays x2 and y-column are of different length "+ this.mPoints + " " + y[0].length);
        	if(this.lPoints!=y[0][0].length)throw new IllegalArgumentException("Arrays x3 and y-column are of different length "+ this.mPoints + " " + y[0][0].length);
        	if(this.kPoints!=y[0][0][0].length)throw new IllegalArgumentException("Arrays x4 and y-column are of different length "+ this.kPoints + " " + y[0][0][0].length);
          	if(this.nPoints<3 || this.mPoints<3 || this.lPoints<3 || this.kPoints<3)throw new IllegalArgumentException("The tabulated 4D array must have a minimum size of 3 X 3 X 3 X 3");

        	this.csm = new CubicSpline(this.nPoints);
        	this.tcsn = TriCubicSpline.oneDarray(this.nPoints, this.mPoints, this.lPoints, this.kPoints);
        	this.x1 = new double[this.nPoints];
        	this.x2 = new double[this.mPoints];
        	this.x3 = new double[this.lPoints];
        	this.x4 = new double[this.kPoints];

        	this.y = new double[this.nPoints][this.mPoints][this.lPoints][this.kPoints];
        	this.d2ydx2 = new double[this.nPoints][this.mPoints][this.lPoints][this.kPoints];
        	for(int i=0; i<this.nPoints; i++){
            		this.x1[i]=x1[i];
        	}
        	this.xMin[0] = Fmath.minimum(this.x1);
        	this.xMax[0] = Fmath.maximum(this.x1);

        	for(int j=0; j<this.mPoints; j++){
            		this.x2[j]=x2[j];
        	}
        	this.xMin[1] = Fmath.minimum(this.x2);
        	this.xMax[1] = Fmath.maximum(this.x2);

        	for(int j=0; j<this.lPoints; j++){
            		this.x3[j]=x3[j];
        	}
        	this.xMin[2] = Fmath.minimum(this.x3);
        	this.xMax[2] = Fmath.maximum(this.x3);

        	for(int j=0; j<this.kPoints; j++){
            		this.x4[j]=x4[j];
        	}
        	this.xMin[3] = Fmath.minimum(this.x4);
        	this.xMax[3] = Fmath.maximum(this.x4);

        	for(int i =0; i<this.nPoints; i++){
            		for(int j=0; j<this.mPoints; j++){
            		    for(int k=0; k<this.lPoints; k++){
            		        for(int l=0; l<this.kPoints; l++){
                		        this.y[i][j][k][l]=y[i][j][k][l];
                		    }
                		}
            		}
        	}
    	}

    	//  METHODS

    	// Reset the default handing of identical abscissae with different ordinates
        // from the default option of separating the two relevant abscissae by 0.001 of the range
        // to avraging the relevant ordinates
    	public void averageIdenticalAbscissae(){
    	    this.averageIdenticalAbscissae = true;
    	    for(int i=0; i<this.tcsn.length; i++)this.tcsn[i].averageIdenticalAbscissae();
    	    this.csm.averageIdenticalAbscissae();
    	}

   	    // Get minimum limits
    	public double[] getXmin(){
    	    return this.xMin;
    	}

    	// Get maximum limits
    	public double[] getXmax(){
    	    return this.xMax;
    	}

    	// Get limits to x
    	public double[] getLimits(){
    	    double[] limits = {xMin[0], xMax[0], xMin[1], xMax[1], xMin[2], xMax[2], xMin[3], xMax[3]};
    	    return limits;
    	}

    	// Display limits to x
    	public void displayLimits(){
    	    System.out.println(" ");
    	    for(int i=0; i<2; i++){
    	        System.out.println("The limits to the x array " + i + " are " + xMin[i] + " and " + xMax[i]);
    	    }
    	    System.out.println(" ");
    	}

    	//	Returns an interpolated value of y for values of x1, x2, x3 and x4
    	//  	from a tabulated function y=f(x1,x2,x3,x4)
    	public double interpolate(double xx1, double xx2, double xx3, double xx4){

            double[][][] yTempml = new double[this.mPoints][this.lPoints][this.kPoints];
            for(int i=0; i<this.nPoints; i++){
                String workingIndices = new String(subMatrixIndices);
	        	workingIndices += "QuadriCubicSpline rows  " + i;
	        	for(int j=0; j<this.mPoints; j++){
	        	    workingIndices += ",  " + j;
	        	    for(int k=0; k<this.lPoints; k++){
	        	        workingIndices += ",  " + k + ": ";
            		    for(int l=0; l<this.kPoints; l++){
	        	            yTempml[j][k][l]=y[i][j][k][l];
	        	        }
	        	    }
	        	}
	        	this.tcsn[i].setSubMatrix(workingIndices);
	        	this.tcsn[i].resetData(x2,x3,x4,yTempml);
	    	}
	    	double[] yTempm = new double[nPoints];

	    	for (int i=0;i<nPoints;i++){
	    	    if(this.derivCalculated)this.tcsn[i].setDeriv(d2ydx2[i]);
		    	yTempm[i]=this.tcsn[i].interpolate(xx2, xx3, xx4);
		    	if(!this.derivCalculated)d2ydx2[i] = this.tcsn[i].getDeriv();
	    	}
	    	derivCalculated = true;

	    	String workingIndices = new String(subMatrixIndices);
            workingIndices += "QuadriCubicSpline interpolated column:  ";
	    	this.csm.setSubMatrix(workingIndices);
	    	this.csm.resetData(x1,yTempm);
	    	return this.csm.interpolate(xx1);
    	}

    	// Return second derivatives
    	public double[][][][] getDeriv(){
    	    return this.d2ydx2;
    	}
}

