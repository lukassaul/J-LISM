import java.util.StringTokenizer;
import java.util.Vector;

/**
* Adapted for heliosphere display-
*
* Lets generalize this to read from an unknown 3d (color +2d) field
* where we don't know how many points there are.
*
*/
public class IBEX_image {

	public static double AU = 1.49598* Math.pow(10,11); //meters


	public IBEX_image(String filename) {

		file f = new file(filename);
		//int num2 = 101;
		Vector xVec = new Vector();
		Vector yVec = new Vector();
		Vector zVec = new Vector();
		String garbage = "garbage";
		f.initRead();
		String line = ""; // garbage - header
		int index = 0;
		while ((line=f.readLine())!=null) {

			int x=0;
			int y=0;
			int z=0;
			double num=0;
			StringTokenizer st = new StringTokenizer(line);
			try {

				x = (int)(1000.0*(double)Float.parseFloat(st.nextToken()));
				y = (int)(1000.0*(double)Float.parseFloat(st.nextToken()));
				garbage = st.nextToken();
				z = (int)(1000.0*(double)Float.parseFloat(st.nextToken()));
				Integer x_d = new Integer(x);
				Integer y_d = new Integer(y);
				Integer z_d = new Integer(z);
				if (!xVec.contains(x_d)) xVec.add(x_d);
				if (!yVec.contains(y_d)) yVec.add(y_d);
				zVec.add(z_d);

			}
			catch (Exception e) {e.printStackTrace();}
			index ++;
		}

		// lets take a look at what we got..  not sure I like "contains"..
		for (int i=0; i<yVec.size(); i++) 	o("y_"+i+" : " + (Integer)yVec.elementAt(i));
		for (int i=0; i<xVec.size(); i++) 	o("x_"+i+" : " + (Integer)xVec.elementAt(i));


		System.out.println("x leng: " + xVec.size());
		System.out.println("y leng: " + yVec.size());
		System.out.println("z leng: " + zVec.size());

		double[] xArray = new double[xVec.size()];
		double[] yArray = new double[yVec.size()];
		double[] zArray = new double[zVec.size()];


		for (int i = 0; i<xVec.size(); i++) {
			xArray[i]=(double)(Integer)xVec.elementAt(i);
			xArray[i]/=1000.0;
		}

		for (int i = 0; i<yVec.size(); i++) {
			yArray[i]=(double)(Integer)yVec.elementAt(i);
			yArray[i]/=1000.0;
		}

		for (int i = 0; i<zVec.size(); i++) {
			zArray[i]=(double)(Integer)zVec.elementAt(i);
			zArray[i]/=1000.0;
			//zArray[i] = 10.0+Math.log(zArray[i]);
		}

		// get min and max of zarray
		double minz = 1000.0;
		double maxz = 0.0;
		for (int j=0; j<zArray.length; j++) {
			if (zArray[j]>maxz) maxz=zArray[j];
			if (zArray[j]<minz) minz=zArray[j];
		}
		System.out.println("minz: "+minz+" maxz: "+maxz);

		double[][] zArray2 = new double[xVec.size()][yVec.size()];
		for (int i=0; i<xArray.length; i++) {
			for (int j=0; j<yArray.length; j++) {
				zArray2[i][j]= minz;
			}
		}

		// read file again and assign z value in double array to fix any problems in grid
		f.initRead();
		line = ""; // garbage - header
		index = 0;
		while ((line=f.readLine())!=null) {
			int x=0;
			int y=0;
			double z=0.0;
			double num=0;
			StringTokenizer st = new StringTokenizer(line);
			try {
				x = (int)(1000.0*(double)Float.parseFloat(st.nextToken()));
				y = (int)(1000.0*(double)Float.parseFloat(st.nextToken()));

				garbage = st.nextToken();
				z = Double.parseDouble(st.nextToken());
				//z = 10.0+Math.log(Double.parseDouble(st.nextToken()));
				Integer x_d = new Integer(x);
				Integer y_d = new Integer(y);
				int x_ind = xVec.indexOf(x_d);
				int y_ind = yVec.indexOf(y_d);
				if (z<1.0) zArray2[x_ind][y_ind]=0.0;
				else zArray2[x_ind][y_ind]=z;
			}
			catch (Exception e) {e.printStackTrace();}
			index ++;
		}


		zArray = new double[xArray.length*yArray.length];
		int iii = 0;
		for (int i=0; i<xArray.length; i++) {
			for (int j=0; j<yArray.length; j++) {
				zArray[iii]=zArray2[i][j];
				iii++;
			}
		}

		for (int i = 0; i<yVec.size(); i++) {
			yArray[i]=yArray[i]*180.0/Math.PI;
		}

		JColorGraph jcg = new JColorGraph(xArray,yArray,zArray,true);

		String unitString = "Counts";
		//if (histType == "Flux") unitString = "Diff. Flux (1/cm^2/s/sr/keV)";
		//if (histType == "Distribution Function") unitString = "Dist. Function (s^3/km^6)";
		//if (histType == "Log Distribution Function") unitString = "log Dist. Function (s^3/km^6)";
		jcg.setLabels("IBEX-LO Simulation","UniBern 2008",unitString);
		jcg.run();
		jcg.showIt();
	}

	public static final void main(String[] args) {
		IBEX_image m = new IBEX_image(args[0]);
	}

	public static final void o(String arg) {
		System.out.println(arg);
	}
}