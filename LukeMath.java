/** This is a wrapper for some static functions I wrote
*
*/

public class LukeMath {

	// recursively computes legendre polynomial of order n, argument x
	public static double compute(int n,double x) {
			if (n == 0) return 1;
			else if (n == 1) return x;
			return ( (2*n-1)*x*compute(n-1,x) - (n-1)*compute(n-2,x) ) / n;
		}
	}


	// Guess we don't need this stuff...

	/*public static long factorial(int m) {
			if (m==1) return 1;
			else return m*factorial(m-1);
		}
	}*/// L

	/*C***    LEGENDRE POLYNOM P(X) OF ORDER N
	C***	EXPLICIT EXPRESSION

		   POLLEG=0.
		   NEND=N/2
			DO M=0,NEND
		        	N2M2=N*2-M*2
				 NM2=N-M*2
				  NM=N-M
			    TERM=X**NM2*(-1.)**M
			    TERM=TERM/NFAK(M)*NFAK(N2M2)
			    TERM=TERM/NFAK(NM)/NFAK(NM2)
			  POLLEG=POLLEG+TERM
			END DO
		  POLLEG=POLLEG/2**N

		 RETURN

		END
	C*/

// I love fortran!