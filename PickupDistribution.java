import java.lang.Math;
import java.util.*;


/*  Attempt to model pickup distribution including Vincident and B in 3D

 	Warsaw - July 2000 - Lukas Saul


    Parameters from PickMain (Moebius et al) :

    "Solar wind: Velocity [km/s], Azimuth[deg],Elevation[deg], Density [cm^(-3)]y"
	     540.000    177.000  0.    5.00000
	 "IMF:        Magnitude [nT],  Azimuth [deg], Elevation [deg]"
	    10.00000    135.000  0.
	 "Neutr.gas: Density [cm^(-3)],Velocity [km/s],Loss rate [1/s],Ioniz. rate [1/s]"
	     7.00000E-03    22.5000    6.90000E-08    6.50000E-08
	 "Pickup ions:	Mass [amu],	Charge [e],	Mean Free Path [AU]"
	     4.00000    1.00000    1.00000E-01
	 "Case studies: Interpl. cond., Neutral gas , Velocity Evolution Model."
	     1.00000    3.00000    1.00000
	 "Distr.func.: Conv.crit.(dF/F), Diff. Param. Chi0,  Adiab. Exp. GAMMA"
	     1.00000E-02   0.500000    2.00000
	 "SC attitude (north ecl.pole: elev.=90 deg) and position: "
	 "             Elevation [deg], Azimuth [deg],	Sun distance [AU]"
	     90.0000  0.    1.00000
	 "Vel. distr.  F(u,V)- calc. grid: N (<72)-angle, M(<28)-velocity"
	     18.0000    32.0000
	 "Accumulation grid in instr. acceptance: E-steps, Ang.steps"
    10.00000   10.00000


 */

// Distribution in a given "parcel" - region of constant vsw, and B field.
public class PickupDistribution {
	public NeutralDistribution NHe;
	public SurvivalProbability LHe;
	public double br,btheta,bphi,vsw,vb,phi0,theta0,T,N0,mu,betaPlus,betaMinus;
	public PickupDistribution(    // all the parameters:

	// sw parameters:  assume radial
			double br, double btheta, double bphi, double vsw
	// vlism parameters:
			double bulkVelocity, double phi0 double theta0 double temperature, double density,
	// ionization parameters:
			double radiationToGravityRatio, double lossRate1AU, double pickupRate1AU)
	{

		this.br = br;
		this.btheta = btheta;
		this.bphi = bphi;
		this.vsw = vsw;
		vb = bulkVelocity;
		this.phi0 = phi0;
		this.theta0 = theta0;
		T = temperature;
		N0 = density;
		mu = radiationToGravityRation;
		betaMinus = lossRate1AU;
		betaPlus = pickupRate1AU;

		NHe = new NeutralDistribution(vb, phi0, theta0, T, N0, mu);
		LHe = new SurvivalProbability(betaMinus);
	}
	// done constructor

	public compute(r,theta,phi, vr, vtheta, vphi) {



}