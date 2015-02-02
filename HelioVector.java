import java.lang.Math;

/**  Utility class for transforming coordinates, etc.
*
*  (this is a general vector class, formulated for use in heliosphere)
*  (GSE not really implemented yet)
*
* Lukas Saul - Warsaw - July 2000
*  Last updated May, 2001
*   renamed to HelioVector,  Oct 2002
*
*/
public class HelioVector {

	/*
	* These static ints are for discerning coordinate systems
	*/
	public static final int CARTESIAN = 1;  // heliocentric
	public static final int SPHERICAL = 2; // heliocentric
	public static final int PLANAR = 3;     // ??
	public static final int GSE = 4;        //??
	public static final int NAHV = 5;     //  not a heliovector
	public static double AU = 1.49598* Math.pow(10,11); //meters

	private double x,y,z;
	private double r,theta,phi;

	/**
	* these booleans are so we know when we need to calculate the conversion
	*    just to make sure we don't do extra work here
	*/
	private boolean haveX, haveY, haveZ, haveR, havePhi, haveTheta;

	/**
	* Default constructor - creates zero vector
	*/
	public HelioVector() {
		haveX=true; haveY=true; haveZ=true; haveR=true; haveTheta=true; havePhi=true;
		x=0; y=0; z=0; r=0; theta=0; phi=0;
	}

	/**
	* constructor sets up SPHERICAL and CARTESIAN coordinates
	* CARTESIAN heliocentric coordinates. x axis = jan 1 00:00:00.
	*/
	public HelioVector(int type, double a, double b, double c) {
		if (type == CARTESIAN) {
			haveX=true; haveY=true; haveZ=true; haveR=false; haveTheta=false; havePhi=false;
			x = a; y = b; z = c;
		}
		if (type == SPHERICAL) {
			haveX=false; haveY=false; haveZ=false; haveR=true; haveTheta=true; havePhi=true;
			r = a;	phi = b; theta = c;  //syntax=  give azimuth first!
		}
	}

	/**
	* Use this to save overhead by eliminating "new" statements,
	*   otherwise same as constructors
	*/
	public void setCoords(int type, double a, double b, double c) {
		if (type == CARTESIAN) {
			haveX=true; haveY=true; haveZ=true; haveR=false; haveTheta=false; havePhi=false;
			x = a; y = b; z = c;
		}
		if (type == SPHERICAL) {
			haveX=false; haveY=false; haveZ=false; haveR=true; haveTheta=true; havePhi=true;
			r = a;	phi = b; theta = c;  //syntax=  give azimuth first!
		}
	}

	/**
	* Use this to save overhead by eliminating "new" statements
	*/
	public void setCoords(HelioVector hp) {
		haveX=true; haveY=true; haveZ=true; haveR=false; haveTheta=false; havePhi=false;
		x = hp.getX(); y = hp.getY(); z = hp.getZ();
	}

	/**
	* Access methods - use these to get coordinates of this point
	*/
	public double getX() {
		if (haveX) return x;
		else {
			x = getR()*Math.cos(getPhi())*Math.sin(getTheta());
			haveX = true;
			return x;
		}
	}

	/**
	* Access methods - use these to get coordinates of this point
	*/
	public double getY(){
		if (haveY) return y;
		else {
			y = getR()*Math.sin(getPhi())*Math.sin(getTheta());
			haveY = true;
			return y;
		}
	}

	/**
	* Access methods - use these to get coordinates of this point
	*/
	public double getZ() {
		if (haveZ) return z;
		else {
			z = getR()*Math.cos(getTheta());
			haveZ = true;
			return z;
		}
	}


	/**
	* Access methods - use these to get coordinates of this point
	*/
	public double getR() {
		if (haveR) return r;
		else {
			r = Math.sqrt(getX()*x + getY()*y + getZ()*z);
			haveR = true;
			return r;
		}
	}

	/**
	* Access methods - use these to get coordinates of this point
	*/
	public double getPhi() {
		if (havePhi) return phi;
		else {
			// this stuff is because Math.atrig(arg) returns -pi/2<theta<pi/2 :
			// we want 0<phi<2pi
			if (getX()>0) phi = Math.atan(getY()/x);
			else if (x<0) phi = Math.atan(y/x) + Math.PI;
			else if (x==0 & y>0) phi = Math.PI/2;
			else if (x==0 & y<0) phi = 3*Math.PI/2;
			else if (x==0 & y==0) phi = 0;
			havePhi = true;
			if (phi<0) phi+=Math.PI*2;
			else if (phi>Math.PI*2) phi-=Math.PI*2;

			return phi;
		}
	}

	/**
	* Access methods - use these to get coordinates of this point
	*/
	public double getTheta() {
		if (haveTheta) return theta;
		else if(getR()==0) { theta=0; haveTheta=true; return 0; }
		else {
			// we want theta>=0 & theta<= PI
			theta = Math.acos(getZ()/getR());
			//else if (z<0) theta = Math.PI + Math.acos(z/r);
			//else if (z==0) theta = Math.PI/2;
			haveTheta = true;
			//if (theta<0) theta+=|(theta>Math.PI)) System.out.println("theta? " + theta);
			return theta;
		}
	}

	// SOME VECTOR OPERATIONS
	// each operation comes in two forms, one creating a new object
	// and one changing an argument

	/**
	* Add two vectors easily
	*/
	public HelioVector sum(HelioVector hp) {
		return new HelioVector( CARTESIAN,getX()+hp.getX(),getY()+hp.getY(),getZ()+hp.getZ() );
	}
	public static void sum(HelioVector x, HelioVector y) {
		x.setCoords(CARTESIAN, x.getX()+y.getX(), x.getY()+y.getY(), x.getZ()+y.getZ());
	}

	/**
	* Subtract two vectors easily
	*/
	public HelioVector difference(HelioVector hp) {
		return new HelioVector( CARTESIAN,getX()-hp.getX(),getY()-hp.getY(),getZ()-hp.getZ() );
	}
	public static void difference(HelioVector x, HelioVector y) {
		x.setCoords( CARTESIAN, x.getX()-y.getY(), x.getY()-y.getY(), x.getZ()-y.getZ() );
	}

	/**
	* Multiply a vector by a scalar
	*/
	public HelioVector product(double d) {
		return new HelioVector( SPHERICAL, d*getR(), getPhi(), getTheta() );
	}
	public static void product(HelioVector x, double d) {
		x.setCoords( SPHERICAL, d*x.getR(), x.getPhi(), x.getTheta() );
	}
	public HelioVector product(int d) {
		return new HelioVector( SPHERICAL, d*getR(), getPhi(), getTheta() );
	}

	/**
	* returns vector of same size pointing in opposite direction -
	*    ADDITIVE inverse NOT MULTIPLICATIVE
	*
	*/
	public HelioVector invert() {
		return new HelioVector(CARTESIAN, -getX(), -getY(), -getZ());
	}
	public static void invert(HelioVector x) {
		x.setCoords( CARTESIAN, -x.getX(), -x.getY(), -x.getZ() );
	}

	/**
	* returns unit vector in same direction
	*
	*/
	public HelioVector normalize() {
		return new HelioVector(SPHERICAL, 1, getPhi(), getTheta());
	}
	public static void normalize(HelioVector x) {
		x.setCoords( SPHERICAL, 1, x.getPhi(), x.getTheta() );
	}

	/**
	*  returns standard dot product
	*/
	public double dotProduct(HelioVector hp) {
		return getX()*hp.getX() + getY()*hp.getY() + getZ()*hp.getZ();
	}

	/**
	*  returns standard cross product
	*/
	public HelioVector crossProduct(HelioVector hp) {
		return new HelioVector( CARTESIAN,
					(getY()*hp.getZ() - getZ()*hp.getY()),
					(getZ()*hp.getX() - getX()*hp.getZ()),
					(getX()*hp.getY() - getY()*hp.getX())   );
	}
	public static void crossProduct(HelioVector x, HelioVector y) {
		x.setCoords( CARTESIAN,
					(x.getY()*y.getZ() - x.getZ()*y.getY()),
					(x.getZ()*y.getX() - x.getX()*y.getZ()),
					(x.getX()*y.getY() - x.getY()*y.getX())   );
	}

	/**
	*   this returns the angle between our vector and inflow as expressed
	*   with elevation and azimuth
	*
	*     (assuming inflow is at phi, theta)
	*    (Range = from 0 to PI)
	*/
	public double angleToInflow(double phi0,double theta0) {
		HelioVector lism = new HelioVector(HelioVector.SPHERICAL, 1, phi0, theta0);
		return getAngle(lism);
	}

	/**
	* This calculates the angle between two vectors.
	* Should range from 0 to PI
	*/
	public double getAngle(HelioVector hp) {
		double cosAngle = this.dotProduct(hp)/(getR()*hp.getR());
		return Math.acos(cosAngle);
	}


	/**
	*  HOPEFULLY DEPRECATED BY NOW
	* these routines for converting to a cylindrical coordinate system
	* note: axis must be non-zero vector!
	*/
	public double cylCoordRad(HelioVector axis) {
		o("Deprecated: cylCoordRad in HelioVector");
		double tbr = this.dotProduct(axis) / axis.getR(); // ref plane doesnt matter here
		return tbr;
	}

	public double cylCoordTan(HelioVector axis) {
		o("Deprecated: cylCoordTan in HelioVector");
		return this.crossProduct(axis).getR() / axis.getR(); // lenght of cross product = Vtangent?
	}

	public double cylCoordPhi(HelioVector axis, HelioVector refPlaneVector) {
		o("Deprecated: cylCoordPhi in HelioVector");
		// a bit trickier here...  using convention from Judge and Wu (sin Psi)
		HelioVector newThis = moveZAxis(axis);
		HelioVector newRefPlaneVector = moveZAxis(axis);
		return (newThis.getPhi() - newRefPlaneVector.getPhi() - Math.PI/2);
	}

	/**
	* this routine returns a new helioPoint expressed in terms of new z axis
	* lets try to take the cross product with current z axis and rotate around that.
	*/
	public HelioVector moveZAxis(HelioVector newAxis) {
		HelioVector cp = new HelioVector(CARTESIAN, getY(), -getX(), 0);
		cp = cp.normalize();
		return rotateAroundArbitraryAxis(cp, newAxis.getTheta()).invert();
	}

	/**
	*  This routine uses a matrix transformation to
	*   return a vector (HelioVector) which has been rotated about
	*   the axis (axis) by some degree (t).
	*/
	public HelioVector rotateAroundArbitraryAxis(HelioVector axis, double t) {
		x=getX(); y=getY(); z=getZ(); // make sure we have coords we need
		double s = Math.sin(t);
		double c=Math.cos(t);// useful quantities
		double oc = 1-c ;  // this is unneccessary but I think I have an extra 16bits around somewhere
		double u = axis.getX();
		double v = axis.getY();
		double w = axis.getZ();

		// this is a matrix transformation
		return new HelioVector(CARTESIAN,
				x*(c + u*u*oc) + y*(-w*s + u*v*oc) + z*(v*s + u*w*oc),  //x
				x*(w*s + u*v*oc) + y*(c + v*v*oc) + z*(-u*s + v*w*oc),  //y
				x*(-v*s + u*w*oc) + y*(u*s + v*w*oc) + z*(c + w*w*oc)   //z
					);
	}

	public String o() {
		return new String("X=" + getX()/AU + " Y=" + getY()/AU + " Z=" + getZ()/AU);
	}

	public String toAuString() {
		String tbr = new String("\n");
		tbr += "X = " + getX()/AU + "\n";
		tbr += "Y = " + getY()/AU + "\n";
		tbr += "Z = " + getZ()/AU + "\n";
		tbr += "r = " + getR()/AU + "\n";
		tbr += "phi = " + getPhi() + "\n";
		tbr += "theta = " + getTheta() + "\n";
		return tbr;
	}

	public String toKmString() {
		String tbr = new String("\n");
		tbr += "X = " + getX()/1000 + "\n";
		tbr += "Y = " + getY()/1000 + "\n";
		tbr += "Z = " + getZ()/1000 + "\n";
		tbr += "r = " + getR()/1000 + "\n";
		tbr += "phi = " + getPhi() + "\n";
		tbr += "theta = " + getTheta() + "\n";
		return tbr;
	}

	public String toString() {
		String tbr = new String("\n");
		tbr += "X = " + getX() + "\n";
		tbr += "Y = " + getY() + "\n";
		tbr += "Z = " + getZ() + "\n";
		tbr += "r = " + getR() + "\n";
		tbr += "phi = " + getPhi() + "\n";
		tbr += "theta = " + getTheta() + "\n";
		return tbr;
	}

	/**
	* for testing!!
	*/
	public static final void main(String[] args) {
		// let's test rotation routine
		//
		//  THese are all checked and good to go May 2004
		/*HelioVector hp1 = new HelioVector(CARTESIAN, 1,0,0);
		HelioVector hp2 = new HelioVector(CARTESIAN, 0,1,0);
		HelioVector hp3 = new HelioVector(CARTESIAN, 0,0,1);
		HelioVector hp4 = new HelioVector(CARTESIAN, -1,0,0);
		HelioVector hp5 = new HelioVector(CARTESIAN, 0,-1,0);
		HelioVector hp6 = new HelioVector(CARTESIAN, 0,0,-1);
		HelioVector hp7 = new HelioVector(CARTESIAN, 1,1,0);
		HelioVector hp8 = new HelioVector(CARTESIAN, 0,1,-1);
		HelioVector hp9 = new HelioVector(CARTESIAN, 3,1,0);
		HelioVector hp10 = new HelioVector(CARTESIAN, 0,4,3);
		*/
		HelioVector hp11 = new HelioVector(CARTESIAN, 1.0*AU, 0.01*AU,0.01*AU);
		HelioVector hp12 = new HelioVector(SPHERICAL, 1.0*AU, 0.0, Math.PI/2.0);

		/*o("hp1: " + hp1);
		o("hp2: " + hp2);
		o("hp3: " + hp3);
		o("hp4: " + hp4);
		o("hp5: " + hp5);
		o("hp6: " + hp6);
		o("hp4: " + hp7);
		o("hp5: " + hp8);
		o("hp6: " + hp9);
		o("hp6: " + hp10);
		*/
		o("hp6: " + hp11);
		o("hp12: "+ hp12);



		//HelioVector hp3 = hp1.rotateAroundArbitraryAxis(hp2, Math.PI/2);
		//System.out.println(hp3.getX()+" "+hp3.getY()+" "+hp3.getZ());

		// ok, let's test this cylindrical coordinate bullshit
		//HelioVector hp4 = new HelioVector(CARTESIAN, 0,0,1);
		//HelioVector hp5 = new HelioVector(CARTESIAN, -1,0,0);
		//System.out.println(hp5.cylCoordRad(hp4)+"");
		//System.out.println(hp5.cylCoordTan(hp4)+"");
		//System.out.println(hp5.cylCoordPhi(hp4,new HelioVector(CARTESIAN,1,0,0)));

	}

	public static void o(String s) {
		System.out.println(s);
	}
}