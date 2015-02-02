
//
import flanagan.integration.*;
//import drasys.or.nonlinear.*;

/*
* Utility class for passing around functions
*
*  implements the drasys "FunctionI" interface to enable Simpsons method
*   outsource.
*/
public abstract class Function implements IntegralFunction {
	public double function(double var) {
		return 0;
	}
}
