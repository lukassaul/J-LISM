import java.util.StringTokenizer;
import java.util.Vector;
/*
* This is to load a given model from a file and give a proper normalization to match the ibex data
* for using a one parameter (normalization) fit.
*/
public class FitData {
	public StringTokenizer st;
	public double[] days;
	public double[] rates;
	public Vector daysV;
	public Vector ratesV;

	public FitData(double[] qdays, double[] qrates) {
		days = qdays;
		rates = qrates;
	}

	public FitData(String filename) {
		daysV = new Vector();
		ratesV = new Vector();
		file f = new file(filename);
		f.initRead();
		String line = "";
		while ((line=f.readLine())!=null) {
			st = new StringTokenizer(line);
			daysV.add(st.nextToken());
			ratesV.add(st.nextToken());
		}
		// time to fix the arrays
		days = new double[daysV.size()];
		rates = new double[daysV.size()];
		for (int i=0; i<days.length; i++) {
			days[i]=Double.parseDouble((String)daysV.elementAt(i));
			rates[i]=Double.parseDouble((String)ratesV.elementAt(i));
		}
	}

	// we are going to interpolate here
	public double getRate(double day) {
		for (int i=0; i<days.length; i++) {
			if (day<days[i]) { // this is where we want to be
				return (rates[i]+rates[i+1])/2;
			}
		}
		return 0;
	}
}