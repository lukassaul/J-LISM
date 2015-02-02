import java.util.Date;
import java.util.Random;

public class MCThread extends Thread {
	public boolean finished = false;
	private double[] limits;
	private int np;
	private final FunctionIII f;
	public double tbr;
	private Random r;
	private double w1,w2,w3;


	public MCThread(FunctionIII funk, double[] lims, int num) {
		limits = lims;
		w1 = limits[1]-limits[0];
		w2 = limits[3]-limits[2];
		w3 = limits[5]-limits[4];
		np = num;
		f = funk; // do you feel it
		tbr = 0;
		r = new Random();
	}

	/**
	* Just call the function a certain number of times and add the results together
	* this is the dartboard
	*/
	public void run() {

		tbr = 0.0;
		double xpoint,ypoint,zpoint;
	//	System.out.println("widths: " + w1 + " " + w2 + " " + w3);
		for (int i=0; i<np; i++) {
			// select random point
			xpoint = r.nextDouble()*w1 + limits[0];
			ypoint = r.nextDouble()*w2 + limits[2];
			zpoint = r.nextDouble()*w3 + limits[4];
			tbr += f.function3(xpoint,ypoint,zpoint);
		}
		finished = true;
	}


	/**
	* For testing - use GVLISM
	*/
	public static final void main(String[] args) {
		MultipleIntegration mi = new MultipleIntegration();
		final GaussianVLISM gvli = new GaussianVLISM(new HelioVector(HelioVector.CARTESIAN,
						-27000.0,0.0,0.0), 100.0, 6000.0);
		System.out.println("Created GVLISM" );

		FunctionIII distZw = new FunctionIII () {
			public double function3(double r, double p, double t) {
				HelioVector hv = new HelioVector(HelioVector.SPHERICAL,r,p,t);
				return gvli.heliumDist(hv)*r*r*Math.sin(t);
			}
		};
		double[] limitsZw = {0.0,100000.0,0.0,2*Math.PI,0,Math.PI};

		Date d1 = new Date();
		double ans = mi.mcIntegrate(distZw, limitsZw,5000000);
		Date d2 = new Date();
		o("5000000 answer: " + ans + " Time: " + (d2.getTime() - d1.getTime()));

		d1 = new Date();
		ans = mi.mcIntegrate(distZw, limitsZw,5000000,1);
		d2 = new Date();
		o("5000000,1 answer: " + ans + " Time: " + (d2.getTime() - d1.getTime()));

		d1 = new Date();
		ans = mi.mcIntegrate(distZw, limitsZw,5000000,2);
		d2 = new Date();
		o("5000000,2 answer: " + ans + " Time: " + (d2.getTime() - d1.getTime()));

		d1 = new Date();
		ans = mi.mcIntegrate(distZw, limitsZw,5000000,3);
		d2 = new Date();
		o("5000000,3 answer: " + ans + " Time: " + (d2.getTime() - d1.getTime()));


		System.out.println("answer in spherical coords = " + ans );
	}


	public static void o(String s) {
		System.out.println(s);
	}

}