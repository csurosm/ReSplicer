/*
 * The MIT License
 *
 * Copyright 2016 Mikl&oacute;s Cs&#369;r&ouml;s.
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */

package splice.model;

/**
 * 
 * @author Mikl&oacute;s Cs&#369;r&ouml;s
 */
public class FunctionMinimization 
{
    

    private FunctionMinimization(){}
    
    private static double sign(double a, double b)
    {
        if (b>=0.0)
            return Math.abs(a);
        else
            return -Math.abs(a);
    }
    
    
    private static final double R=0.5*(Math.sqrt(5.0)-1.0); // The golden ratios.
    private static final double C=(1.0-R);
    
    private static final int BRENT_ITMAX=200;
    private static final double ZEPS=1.0e-10;
    // Here ITMAX is the maximum allowed number of iterations; 
    // ZEPS is a small number that protects against trying to achieve fractional accuracy for a minimum that
    // happens to be exactly zero.
    
    /**
     * Brent's method of function minimization
     *
     * 
     * Based on Numerical Recipes 10.2.
     *
     * Given a function f, and given a bracketing triplet of abscissas ax, bx, cx (such that bx is
     * between ax and cx, and f(bx) is less than both f(ax) and f(cx)), this routine isolates
     * the minimum to a fractional precision of about tol using Brent's method. 
     * 
     * @param ax bracket interval endpoint
     * @param bx bracket midpoint 
     * @param cx bracket interval endpoint
     * @param func the function to be minimized
     * @param tol relative length for the smaller bracket length for stopping the bracket shrinking (set to sqrt(precision))
     * @return [x,y] where y is minimum function value, and x is its abscissa
     */

    public static double[] brent(double ax, double bx, double cx, OneParameterFunction func, double tol)
    {
        double x,w,v,fx,fw,fv,u,fu;
        double e=0.0;
        double a=(ax < cx ? ax : cx); // a and b must be in ascending order, but input abscissas need not be. 
        double b=(ax > cx ? ax : cx);
        x=w=v=bx; 
        fw=fv=fx=func.eval(x);

        double d=0.0;
        for (int iter=1;iter<=BRENT_ITMAX;iter++) { // Main program loop.
            double xm=0.5*(a+b);
            double tol1;//=tol*Math.abs(x)+ZEPS;
            double tol2=2.0*(tol1=tol*Math.abs(x)+ZEPS);
            double threshold=tol2-0.5*(b-a);
            //System.out.println("#**B.brent "+iter+" x "+x+", xm "+xm+"; tol1 "+tol1+", tol2 "+tol2+"; a "+a+", b "+b+"; thresh "+threshold+", xdiff "+Math.abs(x-xm));
            if (Math.abs(x-xm) <= threshold) { // Test for done here.
                double[] retval=new double[2];
                retval[0]=x;
                retval[1]=fx;
                return retval;
            }
            if (Math.abs(e) > tol1) { // Construct a trial parabolic fit.
                double r=(x-w)*(fx-fv);
                double q=(x-v)*(fx-fw);
                double p=(x-v)*q-(x-w)*r;
                q=2.0*(q-r);
                if (q > 0.0) p = -p;
                q=Math.abs(q);
                double etemp=e;
                e=d;
                if (Math.abs(p) >= Math.abs(0.5*q*etemp) || p <= q*(a-x) || p >= q*(b-x))
                    d=C*(e=(x >= xm ? a-x : b-x));
                    // The above conditions determine the acceptability of the parabolic fit. Here we
                    // take the golden section step into the larger of the two segments.
                else {
                    d=p/q;//  Take the parabolic step.
                    u=x+d;
                    if (u-a < tol2 || b-u < tol2)
                    d=sign(tol1,xm-x);
                }
            } else {
                e=(x >= xm ? a-x : b-x);
                d=C*e;
            }
            u=(Math.abs(d) >= tol1 ? x+d : x+sign(tol1,d));
            fu=func.eval(u);
            // This is the one function evaluation per iteration.
            if (fu <= fx) { // Now decide what to do with our function evaluation. 
                //System.out.println("#**B.brent "+iter+" fu<=fx u "+u+", fu "+fu+"; fx "+fx+"; v "+v+", w "+w);
                if (u >= x) a=x; else b=x;
                v=w; w=x; x=u;
                fv=fw; fw=fx; fx=fu;
            } else {
                //System.out.println("#**B.brent "+iter+" fu>>fx u "+u+", fu "+fu+"; fx "+fx+"; v "+v+", w "+w);
                if (u < x) a=u; else b=u;
                if (fu <= fw || w == x) {
                    v=w; w=u;
                    fv=fw; fw=fu;
                } else if (fu <= fv || v == x || v == w) {
                    v=u;
                    fv=fu;
                }
            } // Done with housekeeping. Back for another iteration. 
        }
        throw new OptimizationException("Too many iterations in brent");
    }
    /**
     * Exception thrown when too many iteration in one of the routines.
     */
    public static class OptimizationException extends RuntimeException 
    {
        private OptimizationException(String message){
            super(message);
        }
    }
    
    public interface OneParameterFunction 
    {
        /**
         * Calculates the value of a function <var>f</var> at a point.
         * @param param <var>x</var> coordinate
         * @return function value <var>f</var>(<var>x</var>)
         */
        public double eval(double param);

    }
    
}
