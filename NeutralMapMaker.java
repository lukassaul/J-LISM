

/**
* use this in conjunction with other 0LISM classes to make map of neutral densities or temperatures etc.
*
*/
public class NeutralMapMaker {

	public NeutralMapMaker() {

		HelioVector bulk = new HelioVector(HelioVector.SPHERICAL,26000.0,(74+180)*Math.PI/180.0,5.5*Math.PI/180);

		// inflow velocity, density, temperature
		GaussianVLISM gv = new GaussianVLISM(bulk,1,6300.0);

		// neutral distribution, mu, ionization rate
		NeutralDistribution ndH = new NeutralDistribution(gv, 0.0, 3.0*Math.pow(10,-7));

		MultipleIntegration mi = new MultipleIntegration();

		// set up the coordinates
		int res = 100;
		double[] x = new double[res];
		double[] y = new double[res];
		double[] z = new double[res*res];

		//double xMin =
		//double yMin =

		//for (int i=0; i<res; i++)

		int index = 0;
		for (int i=0; i<res; i++) {
			for (int j=0; j<res; j++) {
				z[index]=
				index++;

		JColorGraph jcg;
			if (theMain.monochromeCheckBox.isSelected()) jcg = new JColorGraph(x,y,z,false);
			else jcg = new JColorGraph(x,y,z);

			String unitString = "Diff. EFlux (1/cm^2/s/sr)";
			if (histType == "Flux") unitString = "Diff. Flux (1/cm^2/s/sr/keV)";
			if (histType == "Distribution Function") unitString = "Dist. Function (s^3/km^6)";
			if (histType == "Log Distribution Function") unitString = "log Dist. Function (s^3/km^6)";
			jcg.setLabels("SOHO CTOF He+","1996",unitString);

			jcg.run();
			jcg.showIt();\





	public static final void main(String[] args) {
		NeutralMapMaker nmm = new NeutralMapMaker();
	}
}