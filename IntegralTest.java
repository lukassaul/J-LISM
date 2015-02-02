
//import cern.jet.math.Bessel;
import java.lang.Math;
import java.util.Date;
import drasys.or.nonlinear.*;

public class IntegralTest {

	public static void main(String[] args) {
		Simpsons s = new Simpsons();
		MyF myF = new MyF();
		Integration i = new Integration();
		i.setMaxError(.001);
		//s.setEpsilon(.001);


		try {
			Date d1 = new Date();
			//System.out.println(""+s.getEpsilon());
			System.out.println("Simp: "+s.integrate(myF,0.0,Math.PI));
			Date d2 = new Date();
			System.out.println("" + (d2.getTime() - d1.getTime()));
		}
		catch(Exception e) {}

		i.setFunction(myF);
		i.setStartDivision(4);
		Date d3 = new Date();
		System.out.println("MyGuy: "+i.integrate(0.0,Math.PI)+"");
		Date d4 = new Date();
		System.out.println("" + (d4.getTime() - d3.getTime()));

		/*try {
			s.setEpsilon(.000001);
			System.out.println(""+s.getEpsilon());
			System.out.println(""+s.integrate(myF,1,10));
		}
		catch(Exception e) {}try {
			s.setEpsilon(.5);
			System.out.println(""+s.getEpsilon());
			System.out.println(""+s.integrate(myF,1,10));
		}
		catch(Exception e) {}try {
			s.setEpsilon(.001);
			System.out.println(""+s.getEpsilon());
			System.out.println(""+s.integrate(myF,1,10));
		}
		catch(Exception e) {}try {
			s.setEpsilon(.9);
			System.out.println(""+s.getEpsilon());
			System.out.println(""+s.integrate(myF,1,10));
		}
		catch(Exception e) {}*/

	}
}

class MyF extends Function  {
	public double function(double x) {
		return Math.sin(x);
	}
}


