
//import drasys.or.nonlinear.*;
import gaiatools.numeric.function.Function;

/*
* Utility class for passing around functions
*
*  implements the drasys "FunctionI" interface to enable Simpsons method
*   outsource.
*/
public abstract class FunctionI extends Function {
	public double value(double var) {
		return 0;
	}
}
