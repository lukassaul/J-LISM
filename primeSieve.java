import SimpleInput;
import SimpleOutput;


public class primeSieve {
    // This program implements the Sieve of Eratosthenes

    public static void main (String args []) {

    // declare basic objects
    SimpleInput in = new SimpleInput();
    SimpleOutput out = new SimpleOutput();
    int M;   // size of sieve
    int i, number;   // indices

    // Determine size of sieve
    out.println ("Use of the Sieve of Eratosthenes to compute primes");
    out.println ("Enter how large an integer to process:  ");
    M = in.readInt();

    // set up sieve
    boolean [] crossedOut = new boolean [M+1];
        // use size M+1 for elemens 0 through M
    for (i = 0; i<= M; i++)  // i++ is an abbreviation for i = i + 1
        crossedOut[i] = false;

    // follow cross-out process
    for (number = 2; number < M; number++)
        {if (!crossedOut[number])
            { // cross out multiples of number
                i = 2*number;
                while (i <= M)
                    { crossedOut[i] = true;
                      i = i + number;
                    }
            }
        }

    // print list of primes
    out.println ("The following prime numbers have been found");
    int numberOnLine = 0;

    for (i = 2; i <= M; i++)
        {if (!crossedOut[i])
            {
                out.print(i + "\t");  // add tab after each number
                numberOnLine++;
                if (numberOnLine == 8)
                    {
                        out.println(); // start a new line every 8 numbers
                        numberOnLine = 0;
                    }
            }
        }
    out.println();
    }
}
