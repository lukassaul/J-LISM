/**
* Test particle motion
* for antimatter study
* Jan. 2012
* TGV Oct. 2012 continuing study
*/

public class GravSim {

	public static final double EV = 1.60218 * Math.pow(10, -19);
	public static double MP = 1.672621*Math.pow(10.0,-27.0);
	public static double AU = 1.49598* Math.pow(10,11); //meters
	public static double EARTHSPEED = AU*2.0*3.14/365.0/24.0/3600.0;
	public static double Ms = 1.98892 * Math.pow(10,30);  // kg
	public static double G = 6.673 * Math.pow(10,-11);  // m^3/s^2/kg

	public double timeStep;
	public double currentMass;
	public double density;
	public double currentRadius;


	//public HelioVector initPos;
	//public HelioVector initVel;
	public HelioVector currentPos;
	public HelioVector currentVel;

	//public HelioVector deltaV;
	public HelioVector aCurrentPos;
	public HelioVector aCurrentVel;

	public GravSim() {

		file outf = new file("gravsim_out.txt");
		outf.initWrite(false);

		currentPos = new HelioVector(HelioVector.CARTESIAN, -100.0*AU, 0.01*AU, 0.01*AU);
		currentVel = new HelioVector(HelioVector.CARTESIAN, 0.0, EARTHSPEED, 0.0);

		currentRadius = 1000.0; // lets start with a km size body
		double vsw = 400.0;
		currentMass = 4.0*Math.PI/3.0*currentRadius*currentRadius*currentRadius;

		timeStep = 30000; // seconds
		int numSteps = 1000000;
		int numOutput = 500;
		int counterLimit = numSteps/numOutput;
		int counter = 0;

		for (int i=0; i<numSteps; i++) {

			if (counter==counterLimit) {
				counter=0;
				System.out.println(i);
				outf.write(timeStep*i/365.0/24.0/3600.0+"\t"+currentPos.getX()/AU+"\t"+currentPos.getY()/AU+"\t"+currentPos.getZ()/AU+"\t"+
											currentVel.getX()/1000.0+"\t"+currentVel.getY()/1000.0+"\t"+currentVel.getZ()/1000.0+"\n");
			}
			else {
				counter++;
			}

			// adjust pos by current vel*timestep
			HelioVector.sum(currentPos,currentVel.product(timeStep));

			// here we use the force equation
			// dv/dt= GM/r/r
			HelioVector.sum(currentVel,currentPos.product(timeStep*-1.0*G*Ms/currentPos.getR()/currentPos.getR()/currentPos.getR()));

			// lets add a force term here proportional to sw incident flux
			// the force vector will be solar wind velocity - object velocity
			double factor = 0.0;
			HelioVector.sum(currentVel, HelioVector.difference(currentVel, new Heliovector(Heliovector.RADIAL, vsw, 0.0,0.0)

		outf.closeWrite();


	}

	public static final void main(String[] args) {
		GravSim gs = new GravSim();
	}
}