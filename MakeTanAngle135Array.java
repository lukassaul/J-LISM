public class MakeTanAngle135Array {

	public static void main(String[] args) {
		VIonBulk vib;
		file myFile;
		double AU = 1.5* Math.pow(10,11);

		myFile = new file("135array.txt");
		myFile.initWrite(false);
		for (int v0=20000; v0<40000; v0=v0+1000) {
			vib = new VIonBulk(AU, 28000, 135*3.14159/180);
			double ourAngle = Math.atan(vib.Vr()/vib.Vperp());
			myFile.write(v0+" "+ourAngle+System.getProperty("line.separator"));
		}
		myFile.closeWrite();
	}
}

