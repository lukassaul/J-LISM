import java.util.Vector;


/**
* DIsplay the results of a NeutralDensity calculation
*   -  lukas saul Nov. 2004

*
*/
public class DensityDisplayAsym {
	private file f;
	private JColorGraph jcg;
	private double[] x,y,z;
	public static double AU = 1.49598* Math.pow(10,11); //meters

	public DensityDisplayAsym() {
		this("testoutput3.dat");
	}

	public DensityDisplayAsym(String filename) {

		f = new file(filename);
		f.initRead();

		Vector xv = new Vector();
		Vector yv = new Vector();
		String line = "";

		// stop when we hit a blank line
		while ( (line=f.readLine()).length()>2 ) {
			xv.add(new Double(Double.parseDouble(line)));
		}

		// stop when we hit a blank line
		while ( (line=f.readLine()).length()>2 ) {
			yv.add(new Double(Double.parseDouble(line)));
		}



		x = new double[xv.size()];
		y = new double[yv.size()];
		z = new double[x.length * y.length];

		// cast vector to array
		for (int i=0; i<x.length; i++) {
			x[i] = ((Double)xv.elementAt(i)).doubleValue();
			//if (x[i]>AU/10)
			x[i] = x[i]/AU;

			System.out.println("i: " + i + " x[i]: " + x[i]);
		}
		for (int i=0; i<x.length; i++) {
			y[i] = ((Double)yv.elementAt(i)).doubleValue();
			//if (x[i]>AU/10)
			y[i] = y[i]/AU;

			System.out.println("i: " + i + " y[i]: " + y[i]);
		}


		for (int i=0; i<z.length; i++) {
			line=f.readLine();
			z[i]=Double.parseDouble(line);
		}

		jcg = new JColorGraph(x,y,z);
		jcg.xLabel1 = "x";
		jcg.xLabel2 = "AU";
		jcg.yLabel1 = "y";
		jcg.yLabel2 = "AU";
		jcg.run();
		jcg.showIt();
	}


	/**
	* run this bad boy
	*
	*/
	public static final void main (String[] args) {
		if (args.length==1)	{
			DensityDisplayAsym dd = new DensityDisplayAsym(args[0]);
		}
		else { DensityDisplayAsym dd = new DensityDisplayAsym();}
	}
}