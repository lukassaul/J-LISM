/*
*   Class   ComplexPoly
*
*   Defines a complex polynomial
*   y = a[0] + a[1].x + a[2].x^2 + a[3].3 + . . . + a[n].x^n
*   where x and all a[i] may be real or complex
*   and deg is the degree of the polynomial, i.e. n,
*   and includes the methods associated with polynomials,
*   e.g. complex root searches
*
*   WRITTEN BY: Dr Michael Thomas Flanagan
*
*   See class Complex for standard complex arithmetic
*
*   DATE:    February 2002
*   UPDATED: 22 June 2003, 19 January 2005, 12 May 2005, 11 October 2005, 30 April 2007, 22 November 2007
*
*   DOCUMENTATION:
*   See Michael Thomas Flanagan's Java library on-line web pages:
*   http://www.ee.ucl.ac.uk/~mflanaga/java/ComplexPoly.html
*   http://www.ee.ucl.ac.uk/~mflanaga/java/
*
*
*   Copyright (c) October 2005, April 2007
*
*   PERMISSION TO COPY:
*   Permission to use, copy and modify this software and its documentation for
*   NON-COMMERCIAL purposes is granted, without fee, provided that an acknowledgement
*   to the author, Michael Thomas Flanagan at www.ee.ucl.ac.uk/~mflanaga, appears in all copies.
*
*   Dr Michael Thomas Flanagan makes no representations about the suitability
*   or fitness of the software for any or for a particular purpose.
*   Michael Thomas Flanagan shall not be liable for any damages suffered
*   as a result of using, modifying or distributing this software or its derivatives.
*
***************************************************************************************/

package flanagan.complex;

import flanagan.io.FileOutput;

public class ComplexPoly{

        private int deg = 0;            // Degree of the polynomial
        private int degwz = 0;          // Degree of the polynomial with zero roots removed
        private Complex[] coeff;        // Array of polynomial coefficients
        private Complex[] coeffwz;      // Array of polynomial coefficients with zero roots removed

        // CONSTRUCTORS
        public ComplexPoly(int n){
                this.deg = n;
                this.coeff = Complex.oneDarray(n+1);
        }

        // Coefficients are complex
        public ComplexPoly(Complex[] aa){
                this.deg =aa.length-1;
                this.coeff = Complex.oneDarray(this.deg+1);
                for(int i=0; i<=this.deg; i++){
                        this.coeff[i]=Complex.copy(aa[i]);
                }
        }

        // Coefficients are real
        public ComplexPoly(double[] aa){
                this.deg =aa.length-1;
                coeff = Complex.oneDarray(this.deg+1);
                for(int i=0; i<=deg; i++){
                        this.coeff[i].reset(aa[i], 0.0);
                }
        }

        // Single constant -  complex
        // y = aa
        // needed in class Loop
        public ComplexPoly(Complex aa){
                this.deg = 0;
                coeff = Complex.oneDarray(1);
                this.coeff[0]=Complex.copy(aa);
        }

        // Single constant -  double
        // y = aa
        // needed in class Loop
        public ComplexPoly(double aa){
                this.deg = 0;
                coeff = Complex.oneDarray(1);
                this.coeff[0].reset(aa, 0.0);
        }

        // Straight line - coefficients are complex
        // y = aa + bb.x
        public ComplexPoly(Complex aa, Complex bb){
                this.deg = 1;
                coeff = Complex.oneDarray(2);
                this.coeff[0]=Complex.copy(aa);
                this.coeff[1]=Complex.copy(bb);
        }

        // Straight line - coefficients are real
        // y = aa + bb.x
        public ComplexPoly(double aa, double bb){
                this.deg = 1;
                coeff = Complex.oneDarray(2);
                this.coeff[0].reset(aa, 0.0);
                this.coeff[1].reset(bb, 0.0);
        }

        // Quadratic - coefficients are complex
        // y = aa + bb.x + cc.x^2
        public ComplexPoly(Complex aa, Complex bb, Complex cc){
                this.deg = 2;
                coeff = Complex.oneDarray(3);
                this.coeff[0]=Complex.copy(aa);
                this.coeff[1]=Complex.copy(bb);
                this.coeff[2]=Complex.copy(cc);
        }

        // Quadratic - coefficients are real
        // y = aa + bb.x + cc.x^2
        public ComplexPoly(double aa, double bb, double cc){
                this.deg = 2;
                coeff = Complex.oneDarray(3);
                this.coeff[0].reset(aa, 0.0);
                this.coeff[1].reset(bb, 0.0);
                this.coeff[2].reset(cc, 0.0);
        }

        // Cubic - coefficients are complex
        // y = aa + bb.x + cc.x^2 + dd.x^3
        public ComplexPoly(Complex aa, Complex bb, Complex cc, Complex dd){
                this.deg = 3;
                coeff = Complex.oneDarray(4);
                this.coeff[0]=Complex.copy(aa);
                this.coeff[1]=Complex.copy(bb);
                this.coeff[2]=Complex.copy(cc);
                this.coeff[3]=Complex.copy(dd);
        }

        // Cubic - coefficients are real
        // y = aa + bb.x + cc.x^2 + dd.x^3
        public ComplexPoly(double aa, double bb, double cc, double dd){
                this.deg = 3;
                coeff = Complex.oneDarray(4);
                this.coeff[0].reset(aa, 0.0);
                this.coeff[1].reset(bb, 0.0);
                this.coeff[2].reset(cc, 0.0);
                this.coeff[3].reset(dd, 0.0);
        }

        // METHODS

        // Returns a ComplexPoly given the polynomial's roots
        public static ComplexPoly rootsToPoly(Complex[] roots){
            int pdeg = roots.length;

            Complex[] rootCoeff = Complex.oneDarray(2);
            rootCoeff[0] = roots[0].times(Complex.minusOne());
            rootCoeff[1] = Complex.plusOne();
            ComplexPoly rPoly = new ComplexPoly(rootCoeff);
            for(int i=1; i<pdeg; i++){
                    rootCoeff[0] = roots[i].times(Complex.minusOne());
                    ComplexPoly cRoot = new ComplexPoly(rootCoeff);
                    rPoly = rPoly.times(cRoot);
            }
            return rPoly;
        }

        // Reset the polynomial
        public void resetPoly(Complex [] aa){
                if((this.deg+1)!=aa.length)throw new IllegalArgumentException("array lengths do not match");
                for(int i=0; i<this.deg; i++){
                        this.coeff[i] = Complex.copy(aa[i]);
                }
        }

        // Reset a coefficient
        public void resetCoeff(int i, Complex aa){
                this.coeff[i] = Complex.copy(aa);
        }

        // Return a copy of this ComplexPoly [instance method]
        public ComplexPoly copy(){
            if(this==null){
                return null;
            }
            else{
                ComplexPoly aa = new ComplexPoly(this.deg);
                for(int i=0; i<=this.deg; i++){
                        aa.coeff[i] = Complex.copy(this.coeff[i]);
                }
                aa.deg=this.deg;
                return aa;
            }
        }

        // Return a copy of this ComplexPoly [static]
        public static ComplexPoly copy(ComplexPoly bb){
            if(bb==null){
                return null;
            }
            else{
                ComplexPoly aa = new ComplexPoly(bb.deg);
                for(int i=0; i<=bb.deg; i++){
                        aa.coeff[i] = Complex.copy(bb.coeff[i]);
                }
                aa.deg=bb.deg;
                return aa;
            }
        }

        // Clone a ComplexPoly
        public Object clone(){
            if(this==null){
                return null;
            }
            else{
                ComplexPoly aa = new ComplexPoly(this.deg);
                for(int i=0; i<=this.deg; i++){
                        aa.coeff[i] = Complex.copy(this.coeff[i]);
                }
                aa.deg=this.deg;
                return (Object) aa;
            }
        }

        // Return a copy of the polynomial
        public Complex[] polyNomCopy(){
                Complex[] aa = Complex.oneDarray(this.deg+1);
                for(int i=0; i<=this.deg; i++){
                        aa[i] = Complex.copy(this.coeff[i]);
                }
                return aa;
        }

        // Return a reference to the polynomial
        public Complex[] polyNomReference(){
                return this.coeff;
        }
        // Return a reference to the polynomial
        public Complex[] polyNomPointer(){
                return this.coeff;
        }

        // Return a copy of a coefficient
        public Complex coeffCopy(int i){
                return Complex.copy(this.coeff[i]);
        }

        // Return a reference to a coefficient
        public Complex coeffReference(int i){
                return this.coeff[i];
        }

        // Return a reference to a coefficient
        public Complex coeffPointer(int i){
                return this.coeff[i];
        }

        // Return the degree
        public int getDeg(){
                return this.deg;
        }

        // Sets the representation of the square root of minus one to j in Strings
        public void setj(){
                 Complex.setj();
        }

        // Sets the representation of the square root of minus one to i in Strings
        public void seti(){
                Complex.seti();
        }

        // Convert to a String of the form (a+jb)[0] + (a+jb)[1].x + (a+jb)[2].x^2  etc.
        public String toString(){
                String ss = "";
                ss =  ss + this.coeffCopy(0).toString();
                if(this.deg>0)ss = ss + " + (" + this.coeffCopy(1).toString() + ").x";
                for(int i=2; i<=this.deg; i++){
                    ss = ss + " + (" + this.coeffCopy(i).toString() + ").x^" + i;
                }
                return ss;
        }

        // Print the polynomial to screen
        public void print(){
                System.out.print(this.toString());
        }

        // Print the polynomial to screen with line return
        public void println(){
                System.out.println(this.toString());
        }

        // Print the polynomial to a text file with title
        public void printToText(String title){
                title = title + ".txt";
                FileOutput fout = new FileOutput(title, 'n');

                fout.println("Output File for a ComplexPoly");
                fout.dateAndTimeln();
                fout.println();
                fout.print("Polynomial degree is ");
                fout.println(this.deg);
                fout.println();
                fout.println("The coefficients are ");

                for(int i=0;i<=this.deg;i++){
                        fout.println(this.coeff[i]);
                }
                fout.println();
                fout.println("End of file.");
                fout.close();
        }

        // Print the polynomial to a text file without a given title
        public void printToText(){
                String title = "ComplexPolyOut";
                printToText(title);
        }

        // LOGICAL TESTS
        // Check if two polynomials are identical
        public boolean equals(ComplexPoly cp){
            return isEqual(cp);
        }

        public boolean isEqual(ComplexPoly cp){
            boolean ret = false;
            int nDegThis = this.getDeg();
            int nDegCp = cp.getDeg();
            if(nDegThis==nDegCp){
                boolean test = true;
                int i=0;
                while(test){
                    if(!this.coeff[i].isEqual(cp.coeffReference(i))){
                        test = false;
                    }
                    else{
                        i++;
                        if(i>nDegCp){
                            test = false;
                            ret = true;
                        }
                    }
                }
            }
            return ret;
        }

        // Check if two polynomials are identical (static)
        public static boolean isEqual(ComplexPoly cp1, ComplexPoly cp2){
            boolean ret = false;
            int nDegCp1 = cp1.getDeg();
            int nDegCp2 = cp2.getDeg();
            if(nDegCp1==nDegCp2){
                boolean test = true;
                int i=0;
                while(test){
                    if(!cp1.coeffReference(i).isEqual(cp2.coeffReference(i))){
                        test = false;
                    }
                    else{
                        i++;
                        if(i>nDegCp1){
                            test = false;
                            ret = true;
                        }
                    }
                }
            }
            return ret;
        }

        // ADDITION OF TWO POLYNOMIALS
        // Addition,  instance method
        public ComplexPoly plus(ComplexPoly b){
                int n = Math.max(this.deg, b.deg);
                ComplexPoly c = new ComplexPoly(n);
                if(n==this.deg){
                        for(int i=0; i<=n; i++)c.coeff[i]=Complex.copy(this.coeff[i]);
                        for(int i=0; i<=b.deg; i++)c.coeff[i]=this.coeff[i].plus(b.coeff[i]);
                }
                else{
                        for(int i=0; i<=n; i++)c.coeff[i]=Complex.copy(b.coeff[i]);
                        for(int i=0; i<=this.deg; i++)c.coeff[i]=this.coeff[i].plus(b.coeff[i]);
                }
                return c;
        }

        // Addition,  static method
        public static ComplexPoly plus(ComplexPoly a, ComplexPoly b){
                int n = Math.max(a.deg, b.deg);
                ComplexPoly c = new ComplexPoly(n);
                if(n==a.deg){
                        for(int i=0; i<=n; i++)c.coeff[i]=Complex.copy(a.coeff[i]);
                        for(int i=0; i<=b.deg; i++)c.coeff[i]=a.coeff[i].plus(b.coeff[i]);
                }
                else{
                        for(int i=0; i<=n; i++)c.coeff[i]=Complex.copy(b.coeff[i]);
                        for(int i=0; i<=a.deg; i++)c.coeff[i]=a.coeff[i].plus(b.coeff[i]);
                }
                return c;
        }

        // SUBTRACTION OF TWO POLYNOMIALS
        // Subtraction,  instance method
        public ComplexPoly minus(ComplexPoly b){
                int n = Math.max(this.deg, b.deg);
                ComplexPoly c = new ComplexPoly(n);
                if(n==this.deg){
                        for(int i=0; i<=n; i++)c.coeff[i]=Complex.copy(this.coeff[i]);
                        for(int i=0; i<=b.deg; i++)c.coeff[i]=this.coeff[i].minus(b.coeff[i]);
                }
                else{
                        for(int i=0; i<=n; i++)c.coeff[i]=(b.coeff[i]).times(Complex.minusOne());
                        for(int i=0; i<=this.deg; i++)c.coeff[i]=b.coeff[i].plus(this.coeff[i]);
                }
                return c;
        }

        // Subtraction,  static method
        public static ComplexPoly minus(ComplexPoly a, ComplexPoly b){
                int n = Math.max(a.deg, b.deg);
                ComplexPoly c = new ComplexPoly(n);
                if(n==a.deg){
                        for(int i=0; i<=n; i++)c.coeff[i]=Complex.copy(a.coeff[i]);
                        for(int i=0; i<=b.deg; i++)c.coeff[i]=a.coeff[i].minus(b.coeff[i]);
                }
                else{
                        for(int i=0; i<=n; i++)c.coeff[i]=(b.coeff[i]).times(Complex.minusOne());
                        for(int i=0; i<=a.deg; i++)c.coeff[i]=b.coeff[i].plus(a.coeff[i]);
                }
                return c;
        }

        // MULTIPLICATION OF TWO POLYNOMIALS
        // Multiplication,  instance method
        public ComplexPoly times(ComplexPoly b){
                int n = this.deg + b.deg;
                ComplexPoly c = new ComplexPoly(n);
                for(int i=0; i<=this.deg; i++){
                        for(int j=0; j<=b.deg; j++){
                                c.coeff[i+j].plusEquals(this.coeff[i].times(b.coeff[j]));
                        }
                }
                return c;
        }

        // Multiplication,  static method
        public static ComplexPoly times(ComplexPoly a, ComplexPoly b){
                int n = a.deg + b.deg;
                ComplexPoly c = new ComplexPoly(n);
                for(int i=0; i<=a.deg; i++){
                        for(int j=0; j<=b.deg; j++){
                                c.coeff[i+j].plusEquals(a.coeff[i].times(b.coeff[j]));
                        }
                }
                return c;
        }

        // DERIVATIVES
        // Return the coefficients, as a new ComplexPoly,  of the nth derivative
        public ComplexPoly nthDerivative(int n){
                ComplexPoly dnydxn;

                if(n>this.deg){
                        dnydxn = new ComplexPoly(0.0);
                }
                else{
                        dnydxn = new ComplexPoly(this.deg-n);
                        Complex[] nc = Complex.oneDarray(this.deg - n + 1);

                        int k = this.deg - n;
                        for(int i=this.deg; i>n-1; i--){
                                nc[k]=Complex.copy(this.coeff[i]);
                                for(int j=0; j<n; j++){
                                        nc[k]=Complex.times(nc[k], i-j);
                                }
                                k--;
                        }
                        dnydxn = new ComplexPoly(nc);
                }
                return dnydxn;
        }

        // EVALUATION OF A POLYNOMIAL AND ITS DERIVATIVES
        // Evaluate the polynomial
        public Complex evaluate(Complex x){
                Complex y = new Complex();
                if(this.deg==0){
                        y=Complex.copy(this.coeff[0]);
                }
                else{
                        y=Complex.copy(this.coeff[deg]);
                        for(int i=deg-1; i>=0; i--){
                                y=Complex.plus(Complex.times(y, x),this.coeff[i]);
                        }
                }
                return y;
        }

        public Complex evaluate(double xx){
                Complex x =new Complex(xx,0.0);
                Complex y = new Complex();
                if(deg==0){
                        y=Complex.copy(this.coeff[0]);
                }
                else{
                        y=Complex.copy(this.coeff[deg]);
                        for(int i=deg-1; i>=0; i--){
                                y=Complex.plus(Complex.times(y, x),this.coeff[i]);
                        }
                }
                return y;
        }

        // Evaluate the nth derivative of the polynomial
        public Complex nthDerivEvaluate(int n, Complex x){
                Complex dnydxn = new Complex();
                Complex[] nc = Complex.oneDarray(this.deg+1);

                if(n==0)
                {
                        dnydxn=this.evaluate(x);
                        System.out.println("n = 0 in ComplexPoly.nthDerivative");
                        System.out.println("polynomial itself evaluated and returned");
                }
                else{
                        ComplexPoly nthderiv = this.nthDerivative(n);
                        dnydxn = nthderiv.evaluate(x);
                }
                return dnydxn;
        }

        public Complex nthDerivEvaluate(int n, double xx){
                Complex x = new Complex(xx, 0.0);
                return nthDerivEvaluate(n, x);
        }

        // ROOTS OF POLYNOMIALS
        // For general details of root searching and a discussion of the rounding errors
        // see Numerical Recipes, The Art of Scientific Computing
        // by W H Press, S A Teukolsky, W T Vetterling & B P Flannery
        // Cambridge University Press,   http://www.nr.com/

        // Calculate the roots (real or complex) of a polynomial (real or complex)
        // polish = true ([for deg>3 see laguerreAll(...)]
        // initial root estimates are all zero [for deg>3 see laguerreAll(...)]
        public Complex[] roots(){
                boolean polish=true;
                Complex estx = new Complex(0.0, 0.0);
                return roots(polish, estx);
        }

        // Calculate the roots (real or complex) of a polynomial (real or complex)
        // initial root estimates are all zero [for deg>3 see laguerreAll(...)]
        // for polish  see laguerreAll(...)[for deg>3]
        public Complex[] roots(boolean polish){
                Complex estx = new Complex(0.0, 0.0);
                return roots(polish, estx);
        }

        // Calculate the roots (real or complex) of a polynomial (real or complex)
        // for estx  see laguerreAll(...)[for deg>3]
        // polish = true  see laguerreAll(...)[for deg>3]
        public Complex[] roots(Complex estx){
                boolean polish=true;
                return roots(polish, estx);
        }

        // Calculate the roots (real or complex) of a polynomial (real or complex)
        public Complex[] roots(boolean polish, Complex estx){
                if(this.deg==0){
                    System.out.println("degree of the polynomial is zero in the method ComplexPoly.roots");
                    System.out.println("null returned");
                    return null;
                }

                // check for zero roots
                boolean testzero=true;
                int ii=0, nzeros=0;
                while(testzero){
                    if(this.coeff[ii].isZero()){
                        nzeros++;
                        ii++;
                    }
                    else{
                        testzero=false;
                    }
                }
                if(nzeros>0){
                    this.degwz = this.deg - nzeros;
                    this.coeffwz = Complex.oneDarray(this.degwz+1);
                    for(int i=0; i<=this.degwz; i++)this.coeffwz[i] = this.coeff[i+nzeros].copy();
                }
                else{
                    this.degwz = this.deg;
                    this.coeffwz = Complex.oneDarray(this.degwz+1);
                    for(int i=0; i<=this.degwz; i++)this.coeffwz[i] = this.coeff[i].copy();
                }

                // calculate non-zero roots
                Complex[] roots = Complex.oneDarray(this.deg);
                Complex[] root = Complex.oneDarray(this.degwz);

                switch(this.degwz){
                        case 1: root[0]=Complex.negate(this.coeffwz[0].over(this.coeffwz[1]));
                                break;
                        case 2: root=quadratic(this.coeffwz[0],this.coeffwz[1],this.coeffwz[2]);
                                break;
                        case 3: root=cubic(this.coeffwz[0],this.coeffwz[1],this.coeffwz[2], this.coeffwz[3]);
                                break;
                        default: root=laguerreAll(polish, estx);
                }

                for(int i=0; i<this.degwz; i++){
                        roots[i]=root[i].copy();
                }
                if(nzeros>0){
                    for(int i=this.degwz; i<this.deg; i++){
                        roots[i]=Complex.zero();
                    }
                }

                return roots;
        }

        // ROOTS OF A QUADRATIC EQUATION
        // ax^2 + bx + c = 0
        // roots returned in root[]
        // 4ac << b*b accomodated by these methods
        public static Complex[] quadratic(Complex c, Complex b, Complex a){
                double  qsign=1.0;
                Complex qsqrt = new Complex();
                Complex qtest = new Complex();
                Complex bconj = new Complex();
                Complex[] root = Complex.oneDarray(2);

                bconj = b.conjugate();
                qsqrt = Complex.sqrt((Complex.square(b)).minus((a.times(c)).times(4)));

                qtest = bconj.times(qsqrt);

                if( qtest.getReal() < 0.0 ) qsign = -1.0;

                qsqrt = ((qsqrt.over(qsign)).plus(b)).over(-2.0);
                root[0] = Complex.over(qsqrt, a);
                root[1] = Complex.over(c, qsqrt);

                return root;
        }

        public static Complex[] quadratic(double c, double b, double a){
                Complex aa = new Complex(a, 0.0);
                Complex bb = new Complex(b, 0.0);
                Complex cc = new Complex(c, 0.0);

                return quadratic(cc, bb, aa);
        }

        // ROOTS OF A CUBIC EQUATION
        // ddx^3 + ccx^2 + bbx + aa = 0
        // roots returned in root[]
        public static Complex[] cubic(Complex aa, Complex bb, Complex cc, Complex dd){

            Complex r = aa.over(dd);
            Complex q = bb.over(dd);
            Complex p = cc.over(dd);

            Complex a = ((q.times(3)).minus(p.times(p))).over(3);
            Complex b = ((p.times(p).times(p).times(2)).minus(p.times(q).times(9)).plus(r.times(27))).over(27);

            Complex bOver2 = b.over(2);
            Complex aOver3 = a.over(3);
            Complex minusbOver2 = bOver2.negate();
            Complex bOver2squared = bOver2.times(bOver2);
            Complex aOver3cubed = aOver3.times(aOver3).times(aOver3);
            Complex sqrtTerm = Complex.sqrt(bOver2squared.plus(aOver3cubed));
            Complex bigA = minusbOver2.plus(sqrtTerm);
            bigA = Complex.pow(bigA, 1.0/3.0);
            Complex bigB = minusbOver2.minus(sqrtTerm);
            bigB = Complex.pow(bigB, 1.0/3.0);
            Complex plusAplusB = bigA.plus(bigB);
            Complex minusAplusB = plusAplusB.negate();
            Complex plusAminusB = bigA.minus(bigB);
            Complex minusAminusB = plusAminusB.negate();
            Complex sqrtMinus3 = Complex.sqrt(new Complex(-3.0, 0.0));
            Complex[] roots = Complex.oneDarray(3);
            roots[0] = plusAplusB;
            roots[1] = (minusAplusB.plus(plusAminusB.times(sqrtMinus3))).over(2.0);
            roots[2] = (minusAplusB.plus(minusAminusB.times(sqrtMinus3))).over(2.0);

            return roots;
        }


        public static Complex[] cubic(double d, double c, double b, double a){
                Complex aa = new Complex(a, 0.0);
                Complex bb = new Complex(b, 0.0);
                Complex cc = new Complex(c, 0.0);
                Complex dd = new Complex(d, 0.0);

                return cubic(dd, cc, bb, aa);
        }

        // LAGUERRE'S METHOD FOR COMPLEX ROOTS OF A COMPLEX POLYNOMIAL

        // Laguerre method for one of the roots
        // Following the procedure in Numerical Recipes for C [Reference above]
        // estx     estimate of the root
        // coeff[]  coefficients of the polynomial
        // m        degree of the polynomial
        public static Complex laguerre(Complex estx, Complex[] pcoeff, int m){
                double  eps = 1e-7;     // estimated fractional round-off error
                int     mr = 8;         // number of fractional values in Adam's method of breaking a limit cycle
                int     mt = 1000;      // number of steps in breaking a limit cycle
                int     maxit = mr*mt;  // maximum number of iterations allowed
                int     niter = 0;      // number of iterations taken

                // fractions used to break a limit cycle
                double  frac[]={0.5, 0.25, 0.75, 0.13, 0.38, 0.62, 0.88, 1.0};

                Complex root = new Complex();    // root
                Complex b   = new Complex();
                Complex d   = new Complex();
                Complex f   = new Complex();
                Complex g   = new Complex();
                Complex g2  = new Complex();
                Complex h   = new Complex();
                Complex sq  = new Complex();
                Complex gp  = new Complex();
                Complex gm  = new Complex();
                Complex dx  = new Complex();
                Complex x1  = new Complex();
                Complex temp1  = new Complex();
                Complex temp2  = new Complex();

                double  abp = 0.0D, abm = 0.0D;
                double  err = 0.0D, abx = 0.0D;

                for(int i=1; i<=maxit; i++){
                        niter=i;
                        b=Complex.copy(pcoeff[m]);
                        err=Complex.abs(b);
                        d=f=Complex.zero();
                        abx=Complex.abs(estx);
                        for(int j=m-1; j>=0;j--)
                        {
                                // Efficient computation of the polynomial and its first two derivatives
                                f=Complex.plus(Complex.times(estx, f),  d);
                                d=Complex.plus(Complex.times(estx, d),  b);
                                b=Complex.plus(Complex.times(estx, b),  pcoeff[j]);
                                err=Complex.abs(b)+abx*err;
                        }
                        err*=eps;

                        // Estimate of round-off error in evaluating polynomial
                        if(Complex.abs(b)<=err)
                        {
                                root=Complex.copy(estx);
                                niter=i;
                                return root;
                        }
                        // Laguerre formula
                        g=Complex.over(d, b);
                        g2=Complex.square(g);
                        h=Complex.minus(g2, Complex.times(2.0, Complex.over(f, b)));
                        sq=Complex.sqrt(Complex.times((double)(m-1), Complex.minus(Complex.times((double)m, h), g2)));
                        gp=Complex.plus(g, sq);
                        gm=Complex.minus(g, sq);
                        abp=Complex.abs(gp);
                        abm=Complex.abs(gm);
                        if( abp < abm ) gp = gm;
                        temp1.setReal((double)m);
                        temp2.setReal(Math.cos((double)i));
                        temp2.setImag(Math.sin((double)i));
                        dx=((Math.max(abp, abm) > 0.0 ? Complex.over(temp1, gp):Complex.times(Math.exp(1.0+abx),temp2)));
                        x1=Complex.minus(estx, dx);
                        if(Complex.isEqual(estx, x1))
                        {
                                root=Complex.copy(estx);
                                niter=i;
                                return root;     // converged
                        }
                        if ((i % mt)!= 0){
                                estx=Complex.copy(x1);
                        }
                        else{
                                // Every so often we take a fractional step to break any limit cycle
                                // (rare occurence)
                                estx=Complex.minus(estx, Complex.times(frac[i/mt], dx));
                        }
                        niter=i;
                }
                // exceeded maximum allowed iterations
                root=Complex.copy(estx);
                System.out.println("Maximum number of iterations exceeded in laguerre");
                System.out.println("root returned at this point");
                return root;
        }

        // Finds all roots of a complex polynomial by successive calls to laguerre
        // Following the procedure in Numerical Recipes for C [Reference above]
        // Initial estimates are all zero, polish=true
        public Complex[] laguerreAll(){
                Complex estx = new Complex(0.0, 0.0);
                boolean polish = true;
                return laguerreAll(polish, estx);
        }

        //  Initial estimates estx, polish=true
        public Complex[] laguerreAll(Complex estx){
                boolean polish = true;
                return laguerreAll(polish, estx);
        }

        //  Initial estimates are all zero.
        public Complex[] laguerreAll(boolean polish){
                Complex estx = new Complex(0.0, 0.0);
                return laguerreAll(polish, estx);
        }

        // Finds all roots of a complex polynomial by successive calls to laguerre
        //  Initial estimates are estx
        public Complex[] laguerreAll(boolean polish, Complex estx){
                // polish boolean variable
                // if true roots polished also by Laguerre
                // if false roots returned to be polished by another method elsewhere.
                // estx estimate of root - Preferred default value is zero to favour convergence
                //   to smallest remaining root

                int     m = this.degwz;
                double  eps = 2.0e-6;  // tolerance in determining round off in imaginary part

                Complex x = new Complex();
                Complex b = new Complex();
                Complex c = new Complex();
                Complex[] ad = new Complex[m+1];
                Complex[] roots = new Complex[m+1];

                // Copy polynomial for successive deflation
                for(int j=0; j<=m; j++) ad[j]=Complex.copy(this.coeffwz[j]);

                // Loop over each root found
                for(int j=m; j>=1; j--){
                        x=Complex.copy(estx);   // Preferred default value is zero to favour convergence to smallest remaining root
                                                // and find the root
                        x=laguerre(x, ad, j);
                        if(Math.abs(x.getImag())<=2.0*eps*Math.abs(x.getReal())) x.setImag(0.0);
                        roots[j]=Complex.copy(x);
                        b=Complex.copy(ad[j]);
                        for(int jj=j-1; jj>=0; jj--){
                                c=Complex.copy(ad[jj]);
                                ad[jj]=Complex.copy(b);
                                b=(x.times(b)).plus(c);
                        }
                }

                if(polish){
                        // polish roots using the undeflated coefficients
                        for(int j=1; j<=m; j++){
                                roots[j]=laguerre(roots[j], this.coeffwz, m);
                        }
                }

                // Sort roots by their real parts by straight insertion
                for(int j=2; j<=m; j++){
                        x=Complex.copy(roots[j]);
                        int i=0;
                        for(i=j-1; i>=1; i--){
                                if(roots[i].getReal() <= x.getReal()) break;
                                roots[i+1]=Complex.copy(roots[i]);
                        }
                        roots[i+1]=Complex.copy(x);
                }
                // shift roots to zero initial index
                for(int i=0; i<m; i++)roots[i]=Complex.copy(roots[i+1]);
                return roots;
        }
}

