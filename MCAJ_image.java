import java.util.StringTokenizer;

/**
* Adapted for heliosphere display
*
*/
public class MCAJ_image {

	public static double AU = 1.49598* Math.pow(10,11); //meters
	public MCAJ_image(String filename) {

		file f = new file(filename);
		int num2 = 101;
		double[] xArray = new double[num2];
		double[] yArray = new double[num2];
		int in = 0;
		for (double pos = -5.0; pos<5.0; pos+=10.0/((double)num2)) {
			xArray[in]=pos;
			yArray[in]=pos;
			in++;
			System.out.println("pos: " + pos + " next in: " + in);
		}
		double[] zArray = new double[num2*num2];

		for (int i=0; i<zArray.length; i++) zArray[i]=0;

		f.initRead();
		String line = f.readLine(); // garbage - header
		int index = 0;

		while ((line=f.readLine())!=null) {

			double x=0;
			double y=0;
			double num=0;
			StringTokenizer st = new StringTokenizer(line);
			try {

				x = Double.parseDouble(st.nextToken());
				y = Double.parseDouble(st.nextToken());
				String numS = st.nextToken();
				//if (numS.equals("Infinity")) numS="0.001";
				num = Double.parseDouble(numS);
				if (num==0.0) num = 0.001;
				//num = num+1.0;
				//num = Math.log(num);

			}
			catch (Exception e) {e.printStackTrace();}
			zArray[index]=num;
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

		String unitString = "Density";
		//if (histType == "Flux") unitString = "Diff. Flux (1/cm^2/s/sr/keV)";
		//if (histType == "Distribution Function") unitString = "Dist. Function (s^3/km^6)";
		//if (histType == "Log Distribution Function") unitString = "log Dist. Function (s^3/km^6)";
		jcg.setLabels("Heliosphere Simulation","UniBern 2008",unitString);

		jcg.run();
		jcg.showIt();
	}

	public static final void main(String[] args) {
		MCAJ_image m = new MCAJ_image(args[0]);
	}
}
