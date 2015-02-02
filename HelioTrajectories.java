import java.util.Random;

/**
* do some time stepping and track trajectories
*
*/
public class HelioTrajectories {

	public int numToRun = 10000;
	public int timeStep = 10;  // seconds

	public static double KB = 1.38065 * Math.pow(10,-23);
	public static double Ms = 1.98892 * Math.pow(10,30);  // kg
	//public static double Ms = 0.0;  // kg
	public static double G = 6.673 * Math.pow(10,-11);  // m^3/s^2/kg
	public static double AU = 1.49598* Math.pow(10,11); //meters
	public static double NaN = Double.NaN;
	public static double MAX = Double.MAX_VALUE;


   // photon momentum = hv/c
    public static double LA = 2.47 * Math.pow(10,15);   // 1/s
    public static double PC = 6.626 * Math.pow(10,-34);  // J*sec
    public static double CC = 2.998 * Math.pow(10,8);
    public static double MOM = PC*LA/CC;

    public HelioVector cPos;  // current position
    public HelioVector cVel; // current velocity

    public double auFlux = 3.5*Math.pow(10,11); //  photons per cm^2
    public double crossSection = 1.0*Math.pow(10,-14);
    public double auFluxOnAtom = auFlux/crossSection;

    public HelioVector radVec = new HelioVector();

    public HelioTrajectories () {
		r=new Random();

		o("one photon momentum: " + MOM);

		// start them all at the same location..
		public HelioVector initPos = new HelioVector(HelioVector.CARTESIAN,-100*AU,0,0);
		public HelioVector initSpeed = new HelioVector(HelioVector.CARTESIAN,28000.0,0,0);
		cPos = initPos;
		cVel = initVel;
		for (int i=0; i<numToRun; i++) {

			int numSteps = 100000;
			for (int step = 0; step<numSteps; step++) {
				// move the particle according to current velocity
				cPos.sum(cVel.product(timeStep));

				// add to the velocity according to gravity and photon pressure
				radVec = cPos.normalize().invert(); // points to sun
				cVel.sum(radVec.product(timeStep*Ms/G/cPos.getR()/cPos.getR()));  // ad

				// give them a speed according to 6300 degrees kelvin and 28 m/s ba
				// bulk speed is in +x direction
				//  actually lets go with cold model for now.. keep it simple
			}

		while (


	}



	public static final void main(String[] args) {
		HelioTrajectories ht = new HelioTrajectories();
	}

	public void o(String s) {
		System.out.println(s);
	}
}