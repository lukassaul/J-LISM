/**********************************************************
*
*   PolyCubicSpline.java
*
*   Class for performing an interpolation on the tabulated
*   function y = f(x1,x2, x3 .... xn) using a natural cubic splines
*   Assumes second derivatives at end points = 0 (natural spines)
*
*   WRITTEN BY: Dr Michael Thomas Flanagan
*
*   DATE:	    15 February 2006,
*   UPDATES:    9 June 2007, 27 July 2007, 4 December 2007
*
*   DOCUMENTATION:
*   See Michael Thomas Flanagan's Java library on-line web page:
*   http://www.ee.ucl.ac.uk/~mflanaga/java/PolyCubicSpline.html
*   http://www.ee.ucl.ac.uk/~mflanaga/java/
*
*   Copyright (c) February 2006, July 2007  Michael Thomas Flanagan
*
*   PERMISSION TO COPY:
*   Permission to use, copy and modify this software and its documentation for
*   NON-COMMERCIAL purposes is granted, without fee, provided that an acknowledgement
*   to the author, Michael Thomas Flanagan at www.ee.ucl.ac.uk/~mflanaga, appears in all copies.
*
*   Dr Michael Thomas Flanagan makes no representations about the suitability
*   or fitness of the software for any or for a particular purpose.
*   Michael Thomas Flanagan shall not be liable for any damages suffered
*   as a result of using, modifying or distributing this software or its derivatives.
*
***************************************************************************************/

package flanagan.interpolation;

import flanagan.math.Fmath;
import java.lang.reflect.Array;

public class PolyCubicSpline{

    	private int nDimensions = 0;   	    // number of the dimensions of the tabulated points array, y=f(x1,x2,x3 . . xn), i.e. n
    	private Object fOfX = null;         // tabulated values of y = f(x1,x2,x3 . . fn)
    	                                    // as a multidimensional array of double [x1 length][x2 length] ... [xn length]
    	private Object internalArray = null;// Object to store second derivatives
    	                                    // as a multidimensional array of Objects the innermost three layers being double arrays
    	private Object xArrays = null;      // The variable arrays x1, x2, x3 . . . xn
    	                                    // packed as an Object as a multidimensional array of double [][]
    	                                    // where xArrays[0] = array of x1 values, xArrays[1] = array of x2 values etc
    	private double[][] xArray = null;   // The variable arrays x1, x2, x3 . . . xn
        private double yValue = 0.0D;       // returned interpolated value
    	private double[] xMin = null;       // minimum values of the x arrays
    	private double[] xMax = null;       // maximum values of the x arrays
    	private boolean calculationDone = false;    // = true when derivatives calculated
        private String subMatrixIndices = "PolyCubicSpline submatrices ";      // String of indices of the submatrices that have called lower dimension interpolations
        private boolean averageIdenticalAbscissae = false;  // if true: the the ordinate values for identical abscissae are averaged
                                                            // If false: the abscissae values are separated by 0.001 of the total abscissae range;


    	// Constructor
    	public PolyCubicSpline(Object xArrays, Object fOfX){

    	    this.fOfX = Fmath.copyObject(fOfX);
    	    this.xArrays = Fmath.copyObject(xArrays);

    	    // Calculate fOfX array dimension number
    	    Object internalArray = Fmath.copyObject(fOfX);
    	    this.nDimensions = 1;
            while(!((internalArray  =  Array.get(internalArray, 0)) instanceof Double))this.nDimensions++;

            // Repack xArrays as 2 dimensional array if entered a single dimensioned array for a simple cubic spline
            if(this.xArrays instanceof double[] && this.nDimensions == 1){
                double[][] xArraysTemp = new double[1][];
                xArraysTemp[0] = (double[])this.xArrays;
                this.xArrays = (Object)xArraysTemp;
            }
            else{
               if(!(this.xArrays instanceof double[][]))throw new IllegalArgumentException("xArrays should be a two dimensional array of doubles");
            }

            // x -arrays and their limits
            this.xArray = (double[][])this.xArrays;
            this.limits();
        }

        // Limits to x arrays
        private void limits(){
            this.xMin = new double[this.nDimensions];
            this.xMax = new double[this.nDimensions];
            for(int i=0; i<this.nDimensions; i++){
                this.xMin[i] = Fmath.minimum(xArray[i]);
                this.xMax[i] = Fmath.maximum(xArray[i]);
            }
        }

   	    // Get minimum limits
    	public double[] getXmin(){
    	    return this.xMin;
    	}

    	// Get maximum limits
    	public double[] getXmax(){
    	    return this.xMax;
    	}

    	// Get number of dimensions
    	public int getNumberOfDimensions(){
    	    return this.nDimensions;
    	}

    	// Get limits to x
    	public double[] getLimits(){
    	    double[] limits = {xMin[0], xMax[0], xMin[1], xMax[1], xMin[2], xMax[2], xMin[3], xMax[3]};
    	    return limits;
    	}

    	// Display limits to x
    	public void displayLimits(){
    	    System.out.println(" ");
    	    for(int i=0; i<this.nDimensions; i++){
    	        System.out.println("The limits to the x array " + i + " are " + xMin[i] + " and " + xMax[i]);
    	    }
    	    System.out.println(" ");
    	}

    	// Reset the default handing of identical abscissae with different ordinates
        // from the default option of separating the two relevant abscissae by 0.001 of the range
        // to avraging the relevant ordinates
    	public void averageIdenticalAbscissae(){
    	    this.averageIdenticalAbscissae = true;
    	}

    	// Set sub-matrix indices - for use in the recursion
    	public void setSubMatrix(String subMatrixVector){
    	    this.subMatrixIndices = subMatrixIndices;
    	}

    	//  Interpolation method
    	public double interpolate(double[] unknownCoord){

    	    int nUnknown = unknownCoord.length;
    	    if(nUnknown!=this.nDimensions)throw new IllegalArgumentException("Number of unknown value coordinates, " + nUnknown + ", does not equal the number of tabulated data dimensions, " + this.nDimensions);

            int kk = 0;
            switch(this.nDimensions){
                case 0: throw new IllegalArgumentException("data array must have at least one dimension");
                case 1: // If fOfX is one dimensional perform simple cubic spline
                        CubicSpline cs = new CubicSpline(this.xArray[0], (double[])this.fOfX);
                        if(this.calculationDone)cs.setDeriv((double[])this.internalArray);
                        if(this.averageIdenticalAbscissae)cs.averageIdenticalAbscissae();
                        this.yValue = cs.interpolate(unknownCoord[0]);
                        if(!this.calculationDone){
                            double[] deriv = cs.getDeriv();
                            this.internalArray = (Object)deriv;
                            this.calculationDone = true;
                        }
                        break;
                case 2: // If fOfX is two dimensional perform bicubic spline
                        BiCubicSpline bcs = new BiCubicSpline(this.xArray[0], this.xArray[1], (double[][])this.fOfX);
                        if(this.calculationDone)bcs.setDeriv((double[][])this.internalArray);
                        if(this.averageIdenticalAbscissae)bcs.averageIdenticalAbscissae();
                        this.yValue = bcs.interpolate(unknownCoord[0], unknownCoord[1]);
                        if(!this.calculationDone){
                            double[][] deriv = bcs.getDeriv();
                            this.internalArray = (Object)deriv;
                            this.calculationDone = true;
                        }
                        break;
                case 3: // If fOfX is three dimensional perform tricubic spline
                        TriCubicSpline tcs = new TriCubicSpline(xArray[0], xArray[1], xArray[2], (double[][][])this.fOfX);
                        if(this.calculationDone)tcs.setDeriv((double[][][])this.internalArray);
                        if(this.averageIdenticalAbscissae)tcs.averageIdenticalAbscissae();
                        this.yValue = tcs.interpolate(unknownCoord[0], unknownCoord[1], unknownCoord[2]);
                        if(!this.calculationDone){
                            double[][][] deriv = tcs.getDeriv();
                            this.internalArray = (Object)deriv;
                            this.calculationDone = true;
                        }
                        break;
                default:// If fOfX is greater than three dimensional, recursively call PolyCubicSpline
                        //  with, as arguments, the n1 fOfX sub-arrays, each of (number of dimensions - 1) dimensions,
                        //  where n1 is the number of x1 variables.
                        Object obj = fOfX;
                        int dimOne = Array.getLength(obj);
                        double[] csArray = new double [dimOne];
                        double[][] newXarrays= new double[this.nDimensions-1][];
                        double[] newCoord = new double[this.nDimensions-1];
                        for(int i=0; i<this.nDimensions-1; i++){
                            newXarrays[i] = xArray[i+1];
                            newCoord[i] = unknownCoord[i+1];
                        }
                        Object[] objDeriv = new Object[dimOne];
                        if(calculationDone)objDeriv = (Object[])this.internalArray;
                        for(int i=0; i<dimOne; i++){
                            Object objT = (Object)Array.get(obj, i);
                            PolyCubicSpline pcs = new PolyCubicSpline(newXarrays, objT);
                            if(this.averageIdenticalAbscissae)pcs.averageIdenticalAbscissae();
                            String workingIndices = new String(subMatrixIndices);
	        	            workingIndices += "" + i + ", ";
	        	            pcs.setSubMatrix(workingIndices);
                            if(this.calculationDone)pcs.setDeriv(objDeriv[i]);
                            csArray[i] = pcs.interpolate(newCoord);
                            if(!this.calculationDone)objDeriv[i] = pcs.getDeriv();
                        }
                        this.internalArray = (Object)objDeriv;
                        this.calculationDone = true;


                        // Perform simple cubic spline on the array of above returned interpolates
                        CubicSpline ncs = new CubicSpline(xArray[0], csArray);
                        String workingIndices = "PolyCubic Spline interpolated column: ";
         	            ncs.setSubMatrix(workingIndices);
            	    	this.yValue = ncs.interpolate(unknownCoord[0]);
            }

            return this.yValue;
    	}

    	// Set derivatives (internal array)
    	public void setDeriv(Object internalArray){
    	    this.internalArray = internalArray;
    	}

    	// Get derivatives (internal array)
    	public Object getDeriv(){
    	    return this.internalArray;
    	}


}

