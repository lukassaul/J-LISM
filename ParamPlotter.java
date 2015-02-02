import java.util.StringTokenizer;

public class ParamPlotter {

	public ParamPlotter() {

		file f = new file("lVSv_6000_86p5_200k.txt");

		double[] x = new double[30];
		double[] y = new double[30];
		double[] z = new double[900];

		// we are going to set the x and y graph params by hand here..
		// easier than reading from a file because we created the x and y with code previously
		// get the params from IBEXWind.java :

		double res = 30;

		/*
		double lamdaWidth = 20;
		double vWidth = 8000;
		double lamdaMin = 65.0;
		//double tMin = 1000.0;
		double vMin=22000.0;
		double lamdaDelta = lamdaWidth/res;
		double vDelta = vWidth/res;
		*/

/*
		double tempWidth = 7000;
		double vWidth = 10000;
		double tempMin = 1000.0;
		//double tMin = 1000.0;
		double vMin=20000.0;
		double tempDelta = tempWidth/res;
		double vDelta = vWidth/res;
*/

		double lamdaWidth = 20.0;
		double vWidth = 12000;
		double lamdaMin = 65.0;
		//double tMin = 1000.0;
		double vMin=20000.0;

		double lamdaDelta = lamdaWidth/res;
		double vDelta = vWidth/res;

		double lamda = lamdaMin;
		//temp=tMin;
		double v = vMin;



		for (int i=0; i<30; i++) {
			x[i]=lamdaMin + i*lamdaDelta;
			y[i]=vMin  + i*vDelta;
		}


		String line = "";
		f.initRead();

		int zIndex = 0;
		while ((line=f.readLine())!=null) {
			StringTokenizer st = new StringTokenizer(line);
			String garb = st.nextToken();
			garb = st.nextToken();
			try {
				z[zIndex]=Math.log(Double.parseDouble(st.nextToken()));
			}
			catch (Exception e) {
				e.printStackTrace();
			}
			zIndex++;
		}

		System.out.println("going to throw up a plot now ");
		JColorGraph jcg;
		//jcg = new JColorGraph(x,y,z,false);
		jcg = new JColorGraph(x,y,z);

		String unitString = "log sum of squares difference";
		jcg.setLabels("IBEX-LO","2010",unitString);
		jcg.setXAxisLabels("Interstellar He Longitude", "deg");
		jcg.setYAxisLabels("Interstellar He Temperature", "deg K");


//		sg.setXMetaData(new SGTMetaData("Interstellar He Wind Longitude", "deg"));
//		sg.setYMetaData(new SGTMetaData("Interstellar He Wind Temperature", "deg K"));

		jcg.run();
		jcg.showIt();
	}


	public static final void main(String[] args) {
		ParamPlotter pp = new ParamPlotter();
	}

}