import java.lang.Math;
import java.util.Date;

/**
*  Lukas Saul - March 27, 2001
*  Taken from C routine Four1, from Numerical Recipes
*/
public class JFFT {
	public double[] dOutReal;
	public double[] dOutImaginary;
	private double[] dInReal;

	//private long nn;
	//private int iSign;
	//private double wr, wi, wpr, wpi, wtemp, theta;

	// for testing here...
	public static void main(String[] args) {
		double[] testData = new double[10000];
		for (int m=0; m<testData.length; m++) {
			testData[m]=m*m+2*m+3; // random function
		}
		JFFT j = new JFFT(testData);

		//for (int q=0; q<testData.length; q++) {
		//	System.out.println(j.dOutReal[q] + " " + j.dOutImaginary[q]);
		//}
	}


	public JFFT(double[] _dIn) {
		Date d1 = new Date();
		dInReal = _dIn;
		int size = dInReal.length;
		System.out.println(size);

		// what if dIn.length is not a power of 2?
		int test=1;
		int i = 0;
		while (test != 0) {
			i++;
			test = size >> i; // divide by 2^i (binary shift)
			//System.out.println("test = " + test);
		}

		i=i-1;
		System.out.println("final i = " + i);
		int dif = size - (int)Math.pow(2,i);
		System.out.println("dif = " + dif);

		if (dif > 0) { // In this case we need to pad the array with zeroes

			double[] dIn2 = new double[(int)Math.pow(2.0,i+1)]; //this could be faster
			//copy by hand
			for (int k=0; k<dInReal.length; k++) {
				dIn2[k]=dInReal[k];
			}
			for (int j=dInReal.length; j<dIn2.length; j++) {
				dIn2[j]=0;
			}
			// the new data:
			dInReal = dIn2;
			System.out.println("new size " + dInReal.length);
		}
		Date d2 = new Date();
		System.out.println("setup took: " + (d2.getTime() - d1.getTime()) );

		// Set up our discrete complex array
		double[] ourData = new double[dInReal.length*2];
		for (int j=0; j<dInReal.length; j++) {
			ourData[2*j] = dInReal[j];
			ourData[2*j+1] = 0; // set imaginary part to zero
		}

		// Ok, let's do the FFT on the ourData array:
		four1(ourData, dInReal.length, 1);

		dOutReal = new double[dInReal.length];
		dOutImaginary = new double[dInReal.length];

		for (int j=0; j<dInReal.length; j++) {
			dOutReal[j] = ourData[2*j];
			dOutImaginary[j] = ourData[2*j+1];
		}
		Date d3 = new Date();
		System.out.println("fft took: " + (d3.getTime() - d1.getTime()) );
	}

	public void four1(double[] data, int nn, int isign) {
		System.out.println("computing fft with args: " + data.length + " " + nn);
		long n, mmax, m, j, istep, i;
		double wtemp, wr, wpr, wpi, wi, theta;
		float tempr, tempi;

		n = nn/2;
		System.out.println("compare " + n + " " + nn + " " + nn/2);
		j=1;

		for (i=1; i<n; i+=2) {
			if (j>i) {
				swap(data[(int)j], data[(int)i]);
				swap(data[(int)j+1], data[(int)i+1]);
			}
			m=n/2;
			while (m >= 2 && j>m)  {
				j -= m;
				m=m/2;
			}
			j += m;
		}

		// Here begins Danielson-Lanczos section of routine
		mmax=2;
		while (n > mmax) {
			istep = mmax << 1;
			theta = isign*(6.28318530717959/mmax);
			wtemp = Math.sin(0.5*theta);
			wpr = -2.0*wtemp*wtemp;
			wpi = Math.sin(theta);
			wr=1.0;
			wi=0.0;
			for (m=1; m<mmax; m+=2) {
				for (i=m; i<=n; i+=istep) {
					j=i+mmax;
					tempr = (float)(wr*data[(int)j] -
						wi*data[(int)j+1]);
					tempi = (float)(wr*data[(int)j+1] - wi*data[(int)j]);
					data[(int)j] = data[(int)i] - tempr;
					data[(int)j+1] = data[(int)i+1] - tempi;
					data[(int)i] += tempr;
					data[(int)i+1] += tempi;
				}
				wr = (wtemp=wr)*wpr - wi*wpi + wr;
				wi = wi*wpr + wtemp*wpi + wi;
			}
			mmax = istep;
		}
	}

	private void swap(double a, double b) {
		double temp = a;
		a=b;
		b=temp;
	}
}