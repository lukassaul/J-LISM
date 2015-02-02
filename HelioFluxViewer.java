import java.util.StringTokenizer;

/**
* Adapted for heliosphere display
* this one for showing flux as a skymap
*/
public class HelioFluxViewer {

	public static double AU = 1.49598* Math.pow(10,11); //meters
	public HelioFluxViewer (String filename) {

		file f = new file(filename);
		//int num2 = 101;
		double[] xArray = new double[60];  // x is theta
		double[] yArray = new double[361];   // y is phi

		// set up x and y axes for reading file
		for (int i=0; i<xArray.length; i+=1) {
			xArray[i]=(double)i+60.0;
		}
		for (int i=0; i<yArray.length; i+=1) {
			yArray[i]=(double)i;
		}
		double[] zArray = new double[xArray.length*yArray.length];

		for (int i=0; i<zArray.length; i++) zArray[i]=0.0;

		f.initRead();
		String line = "";  //f.readLine(); // garbage - header
		int index = 0;
		String garbage = "";
		while ((line=f.readLine())!=null) {

			double x=0;
			double y=0;
			double num=0;
			StringTokenizer st = new StringTokenizer(line);
			try {
				garbage = st.nextToken();
				x = Double.parseDouble(st.nextToken());
				y = Double.parseDouble(st.nextToken());
				//if (numS.equals("Infinity")) numS="0.001";
				num = Double.parseDouble(st.nextToken());
				//if (num==0.0) num = 0.001;
				//num = num+1.0;
				//num = Math.log(num);

			}
			catch (Exception e) {e.printStackTrace();}
			zArray[index]=num;
			//System.out.println(index + " " + zArray[index]);
			index++;
		}

		System.out.println("index: " + index + " zArray.length: " + zArray.length);

		// get min and max of zarray
		double minz = 1000.0;
		double maxz = 0.0;
		for (int j=0; j<zArray.length; j++) {
			if (zArray[j]>maxz) maxz=zArray[j];
			if (zArray[j]<minz) minz=zArray[j];
		}
		System.out.println("minz: "+minz+" maxz: "+maxz);


		/*	int index = 0;
		for (int i=0; i<x.length; i++) {
			x[i]=(tpTimes[i]-GOODS)/24.0/60.0/60.0;
			for (int j=0; j<y.length; j++) {
				y[j]=tpHist[i].label[j];
				if (tpHist[i].data[j]==0) z[index]=0.0;
				else {
					z[index]=tpHist[i].data[j];
					//o("nonzero z: " + z[index]);
				}
				index++;
			}
		}*/


		JColorGraph jcg = new JColorGraph(xArray,yArray,zArray,true);
		jcg.fileName = filename;

		String unitString = "Counts";
		//if (histType == "Flux") unitString = "Diff. Flux (1/cm^2/s/sr/keV)";
		//if (histType == "Distribution Function") unitString = "Dist. Function (s^3/km^6)";
		//if (histType == "Log Distribution Function") unitString = "log Dist. Function (s^3/km^6)";
		jcg.setLabels("Heliosphere Simulation","UniBern 2008",unitString);

		jcg.run();
		jcg.showIt();
	}

	public static final void main(String[] args) {
		HelioFluxViewer  m = new HelioFluxViewer (args[0]);
	}
}
