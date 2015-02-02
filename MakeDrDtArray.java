public class MakeDrDtArray {

	public static void main(String[] args) {
		VIonBulk vib;
		file myFile;
		double AU = 1.5* Math.pow(10,11);

		myFile = new file("testArray2.txt");
		myFile.initWrite(false);
		for (double theta = -3.14; theta<= 3.14; theta=theta+3.14/30) {
			vib = new VIonBulk(AU, 28000, theta);
			myFile.write(theta+" "+vib.Vr()+" " + vib.Vperp()+System.getProperty("line.separator"));
		}
		myFile.closeWrite();
	}
}

