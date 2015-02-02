
	/**
	*    this routine calculates the distribution function (6D)
	*   however based on 5D inputs as there is symmetry about inflow axis
	*
	*  NOTE:  parameters of distribution must be given in coordinates
    *	    relative to inflow direction for this routine
    *	    Damby & Camm 1957, Mostly from Judge & Wu 1979
    *
    *  Most of the math of the calcuation is in this method!
    *
    * r = dist. from sun
    * theta = angle from inflow
    *
    *
	*/
	public double computeReferenceCoords(
		double r, double theta, double vr, double vt, double vPsi) {

		// we need to make sure that vt is within an appropriate range
		//  this amounts to assuming no non-hyperbolic trajectories
		if (Math.abs(vt) < vtLimit(vr,r)) return 0;

		else {

			Q = Math.sqrt((1-mu)*G*Ms/r);
			//System.out.println(""+(vr*vr+vt*vt-2*Q*Q));

			V0 = Math.sqrt(vr*vr + vt*vt - 2*Q*Q);
			//System.out.println("computed easy numbers "+Q+" "+V0);

			F1 = 2*Vb*V0/(Vtemp*Vtemp) * Math.cos(theta)*(V0 - vr)*(V0 - vr)/(V0*(V0-vr)+Q*Q) -
					((Vb*Vb)+(V0*V0)+2*V0*Vb*Math.cos(theta))/(Vtemp*Vtemp);

			F2 = 2*V0*Vb*vt/(Vtemp*Vtemp) * Math.sin(theta)*(V0-vr)/ (vr*(V0-vr) + (Q*Q));
			//System.out.println("computed F1 & F2 "+F1+" "+F2);

			return  sp.compute(r,theta,vr,vt)*
				N * test1 * Math.exp(F1 + F2*Math.sin(vPsi));

		}
	}

	/**
	* This routine does the analytic integration over psi (V psi)
	*  using the appropriate modified Bessel funciton.
	*
	*  This routine takes coordinates in the reference frame
	*  Useful for computing density quicker.
	*/
	public double computeReference2D(double r, double theta, double vr, double vt)  {
		if (Math.abs(vt) < vtLimit(vr,r)) return 0;

		else {
			Q = Math.sqrt((1-mu)*G*Ms/r);

			V0 = Math.sqrt(vr*vr + vt*vt - 2*Q*Q);

			F1 = 2*Vb*V0/(Vtemp*Vtemp) * Math.cos(theta)*(V0 - vr)*(V0 - vr) / (V0*(V0-vr)+Q*Q) -
							((Vb*Vb)+(V0*V0)+2*V0*Vb*Math.cos(theta))/(Vtemp*Vtemp);

			F2 = 2*V0*Vb/(Vtemp*Vtemp) * Math.sin(theta)*(V0-vr)*vt/ (V0*(V0-vr) + (Q*Q));

			return sp.compute(r,theta,vr,vt) * N * test1 *
				Math.exp(F1)*Bessel.i0(F2);
		}
	}


	/**
	* This tells what the limit on tangential velocity is, based on r and radial velocity
	*  (returns an absolute value that vt must be greater than)
	*/
	public double vtLimit(double vr, double r) {
		double QQ = Math.sqrt((1-mu)*G*Ms/r);
		double t211 = 0;
		if ((t211 = 2*QQ*QQ - vr*vr)<0) return 0;
		else return Math.sqrt(t211) + 0.01;   // this is in m/s so we're OK with this factor
	}


	/**
	* This takes parameters in spherical (heliospheric) coordinates
	*
	* Does appropriate transformation and calculates computeReferenceCoords(args)
	*/

	/*public double compute(double r, double theta, double phi,
						   double vmag, double vtheta, double vphi) {

		hpr = new HelioVector(HelioVector.SPHERICAL, r, phi, theta);
		hpv = new HelioVector(HelioVector.SPHERICAL, vmag, vphi, vtheta);
		//hpin = new HelioVector(HelioVector.Spherical, 1, phi0, theta0); // length doesn't matter here

		return computeReferenceCoords(r, hpr.angleToInflow(phi0,theta0),
		             hpv.cylCoordRad(hpr), hpv.cylCoordTan(hpr), hpv.cylCoordPhi(hpr,hpin));
	}*/

	/**
	* This takes parameters in spherical (heliospheric) coordinates
	*
	* Does appropriate transformation and calculates computeRefernce2D()
	*/

	/*public double compute2D(double r, double theta, double phi,
						   double vmag, double vtheta) {

		hpr = new HelioVector(HelioVector.SPHERICAL, r, phi, theta);
		hpv = new HelioVector(HelioVector.SPHERICAL, vmag, 0, vtheta);
		//hpin = new HelioVector(HelioVector.Spherical, 1, phi0, theta0); // length doesn't matter here

		return computeReference2D(r, hpr.angleToInflow(phi0,theta0),
		             hpv.cylCoordRad(hpr), hpv.cylCoordTan(hpr));
	}*/


	/**
	* This takes coords in a more general frame, as in compute() above
	*  Returns a 3D function in velocity space, with args. vr, vtheta, vphi
	*
	*/
	public FunctionIII getVelocityDistribution(final double r, final double theta, final double phi) {
		return new FunctionIII () {
			public double function3(double x, double y, double z) {
				return dist(r,theta,phi,x,y,z);
			}
		};
	}

	/**
	* This takes coords in a more general frame, as in compute() above
	*  Returns a 2D function in velocity space, with args. vr, vtheta
	*   (we have analytically integrated the Vpsi coordinate over all angles)
	*/
	public FunctionII get2DVelocityDistribution(final double r, final double theta, final double phi) {
		return new FunctionII () {
			public double function2(double x, double y) {
				return computeReference2D(r,theta,x,y);
			}
		};
	}



	/**
	* This takes coords relative to inflow, as in computeRefernceCoords above
	*  Returns a 3D function in velocity space, with args. vr, vtheta, vphi
	*
	*/
	public FunctionIII getReferenceVelocityDistribution(final double r, final double theta) {
		return new FunctionIII () {
			public double function3(double x, double y, double z) {
				return computeReferenceCoords(r,theta,x,y,z);
			}
		};
	}


	/**
	* This takes coords in the reference frame.
	*  Returns a 2D function in velocity space, with args. vr, vtheta.
	*  The function returned is the analytic integration of all psi (Vpsi)
	*/
	public FunctionII get2DReferenceVelocityDistribution(final double r, final double theta) {
		return new FunctionII () {
			public double function2(double x, double y) {
				return computeReference2D(r,theta,x,y);
			}
		};
	}

	/**
	*  This returns the density at a point in the heliosphere by doing
	*   a numerical 3D integration in all velocity space
	*   See MultipleIntegration.java for details
	*   (default is at .001 accuracy)
	*/
	public double density(double r, double theta, double phi) {
		MultipleIntegration mi = new MultipleIntegration();

		// old fashioned way
		double [] lims = new double[6];
		lims[0]=-3*Math.pow(10,5);
		lims[1]=3*Math.pow(10,5); // m/s
		lims[2]=-3*Math.pow(10,5);
		lims[3]=3*Math.pow(10,5);
		lims[4]=-3*Math.pow(10,5);
		lims[5]=3*Math.pow(10,5);
		// do the integral over + and then - velocity space to avoid the zero
		/*double [] lims = new double[6];
		lims[0]=-3*Math.pow(10,5);
		lims[1]=0; // m/s
		lims[2]=-3*Math.pow(10,5);
		lims[3]=0;
		lims[4]=-3*Math.pow(10,5);
		lims[5]=0;

		double [] lims2 = new double[6];
		lims2[0]=0;
		lims2[1]=3*Math.pow(10,5);
		lims2[2]=0;
		lims2[3]=3*Math.pow(10,5);
		lims2[4]=0;
		lims2[5]=3*Math.pow(10,5);*/

		return mi.integrate(getVelocityDistribution(r,theta, phi) , lims);
			   // + mi.integrate(getVelocityDistribution(r,theta, phi) ,lims2);

		// the integration routine will give us stats on the operation.
	}


	/**
	*  This returns the density at a point in the heliosphere by doing
	*   a numerical 2D integration, plus the
	*   analytic integratoin of the PSI coordinate with modified bessel function
	*   See MultipleIntegration.java for details
	*   (default is at .001 accuracy)
	*/
	public double density2D (double r, double theta, double phi) {
		MultipleIntegration mi = new MultipleIntegration();
		double [] lims = new double[4];
		lims[0]=-3*Math.pow(10,7);
		lims[1]=3*Math.pow(10,7);
		lims[2]=-3*Math.pow(10,7);
		lims[3]=3*Math.pow(10,7);

		return mi.integrate(get2DReferenceVelocityDistribution(r,theta),lims);

		// the integration routine will give us stats on the operation.
	}

