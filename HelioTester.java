import java.util.Date;
import java.lang.Math;
public class HelioTester {

	public static void main(String[] args) {

		// THis tests the moveZAxis routine
		Date d1 = new Date();
		//HelioVector h1 = new HelioVector(HelioVector.Spherical, 1,Math.PI/2, Math.PI/4);
		//HelioVector h2 = new HelioVector(HelioVector.Spherical, 1,Math.PI/2, Math.PI/4);



		HelioVector h1 = new HelioVector(HelioVector.Cartesian, 1,0,1);
		HelioVector h2 = new HelioVector(HelioVector.Cartesian, 1,0,0);
		//HelioVector h3 = h2.moveZAxis(h1);
		System.out.println("cylCOordPhi = " + h1.cylCoordPhi(h1,h2));
		//System.out.println("new point: "+h3.getX()+" "+h3.getY()+" "+h3.getZ());


		HelioVector h4 = new HelioVector(HelioVector.Cartesian, 0,1,0);
		HelioVector h5 = new HelioVector(HelioVector.Cartesian, 0,1,0);
		HelioVector h6 = h5.moveZAxis(h1);
		System.out.println("new point: "+h6.getX()+" "+h6.getY()+" "+h6.getZ());


		HelioVector h7 = new HelioVector(HelioVector.Cartesian, 2,0,0);
		HelioVector h8 = new HelioVector(HelioVector.Cartesian, 0,1,0);
		HelioVector h9 = h7.rotateAroundArbitraryAxis(h8,3.14/2);
		System.out.println("new point: "+h9.getX()+" "+h9.getY()+" "+h9.getZ());


		Date d2 = new Date();
		System.out.println(d2.getTime() - d1.getTime()+"");
	}
}