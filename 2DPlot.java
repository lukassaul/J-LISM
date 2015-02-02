
import java.util.StringTokenizer;

public class TwoDPlot {

	public CurveFitter cf;

	public TwoDPlot() {

		file f = new file("second_pass_ascii.txt");
		f.initRead();

		double[] x = new double[30];
		double[] y = new double[30];
		double[] z = new double[900];

		for (int i=0; i<x.length; i++) 	x[i]=Double.parseDouble(f.readLine());
		for (int i=0; i<y.length; i++) 	y[i]=Double.parseDouble(f.readLine());
		for (int i=0; i<z.length; i++) 	z[i]=Double.parseDouble(f.readLine());

		JColorGraph jcg = new JColorGraph(x,y,z);

		String unitString = "log (sum square model error)";
		jcg.setLabels("IBEX-LO","2010",unitString);

		jcg.run();
		jcg.showIt();





	}

	public static boolean contains(double[] set, double test) {
		for (int i=0; i<set.length; i++) {
			if (set[i]==test) return true;
		}
		return false;
	}

	public static void o(String s) {
		System.out.println(s);
	}

	public static final void main(String[] args) {
		TwoDPlot theMain = new TwoDPlot();
	}
}