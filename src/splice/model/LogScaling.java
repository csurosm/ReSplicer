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
 * Arithmetics for log-likelihood scores.
 * Let <var>p</var><sub>1</sub> be the observation's likelihood by hypothesis 1 and 
 * <var>p</var><sub>0</sub> by the null hypothesis. The score is then
 * ln(<var>p</var><sub>1</sub>/<var>p</var><sub>0</sub>) / <var>s</var>, where <var>s</var> 
 * is the scaling value (specified at instantiation). Integer-valued scores are 
 * obtained by floor-rounding.
 * 
 * @author Mikl&oacute;s Cs&#369;r&ouml;s
 */
public final class LogScaling 
{
    /**
     * Scoring with default scale (base e: nats). 
     */
    public LogScaling()
    {
        this(1.0);
    }
    /** 
     * Scoring with arbitrary scale.
     * 
     * @param scale non-0 positive scaling factor (with respect to nats)
     */
    public LogScaling(double scale)
    {
        this.scale = scale;
        //this.scaled_eps = (-(int)(-scaleDouble(LOG_EPS))); // <0
        // log_z (1+z^q) < 1 
        // (1+z^q) < z 
        // z^q < z-1
        // q < log_z (z-1)
        // q < ln (z-1)/ ln(z)
        // q < ln (exp(s)-1) / s
        // q < (s+ln (1-exp(-s))) /s 
        // q < 1 + ln (1-exp(-s))/s
        double log1p_max = 1+Math.log1p(-Math.exp(scale))/scale;
        int log1p_array_size = ((int)-log1p_max)+1;
        scaled_log1plusx = new int[log1p_array_size];
        precompute();
    }
    private final double scale;
    //private int scaled_eps;
    /**
     * Precomputed values for log (1+exp(-<var>i</var>))
     */
    private final int[] scaled_log1plusx;
    
    /**
     * Phred scale: 10 log<sub>10</sub> <var>x</var>
     */
    public static final double SCALE_PHRED = 0.1*Math.log(10.0);
    /**
     * Binary scale: log<sub>2</sub> <var>x</var>
     */
    public static final double SCALE_BITS = Math.log(2.0);
    /**
     * Decimal scale: log<sub>10</sub> <var>x</var>
     */
    public static final double SCALE_BANS = Math.log(10.0);
    /**
     * Stretched Phred scale: 1000 log<sub>10</sub> <var>x</var>
     */
    public static final double SCALE_MILLIBANS = 0.001*Math.log(10.0);
    
    /**
     * Scaled value for given nats
     * 
     * @param nats log-likelihood ratio in nats (using <code>Math.log</code>)
     * @return scaled value (nats/{@link #getScaleNats() }) 
     */
    public double scaleDouble(double nats)
    {
        return nats / scale;
    }
    
    /**
     * Scaling factor in nats. 
     * @return scaling factor
     */
    public final double getScaleNats()
    {
        return scale;
    }
    
    private double exp(int q)
    {
        return Math.exp(q*scale);
    }
    
    /**
     * Fills up {@link #scaled_log1plusx} array.
     */
    private void precompute()
    {
        for (int i=0; i<scaled_log1plusx.length; i++)
        {
            double p = exp(-i);
            double m = Math.log1p(p);
            int g = (int)scaleDouble(m);
            scaled_log1plusx[i]=g;
        }
    }
    
    
    /**
     * Array of scores for given empirical counts.
     * 
     * @param count1 counts in hypothesis 1
     * @param count0 null hypothesis frequencies
     * @param laplace_pseudocount Pseudocount added to count0[[] entries before normalization
     * @return unscaled scores (in nats)
     */
    public double[] getUnscaledScore(int[] count1, int[] count0, int laplace_pseudocount)
    {
        for (int i=0; i<count0.length; i++)
            count0[i] += laplace_pseudocount;
        double[] s = getUnscaledScore(count1, count0);
        for (int i=0; i<count0.length; i++)
            count0[i] -= laplace_pseudocount;
        return s;
    }
    
    /**
     * Array of scores for given empirical counts.
     * 
     * @param count1 counts in hypothesis 1
     * @param count0 null hypothesis frequencies
     * @return unscaled scores (in nats)
     */
    public static double[] getUnscaledScore(int[] count1, int[] count0)
    {
        assert (count1.length==count0.length);
        int sum1=0;
        for (int c: count1){ sum1+= c;}
        int sum0=0;
        for (int c: count0){ sum0+=c;}
        double[] llr = new double[count1.length];
        for (int idx=0; idx<count1.length; idx++)
        {
            double n1 = count1[idx];
            double p1 = n1/sum1;
            double n0 = count0[idx];
            double p0 = n0/sum0;
            double r = getUnscaledLogRatio(p1, p0);
            llr[idx] = r;
        }
        return llr;
    }
    /**
     * Kullback-Leibler divergence between two empirical distributions.
     * Probabilities are calculated from counts as 
     * <var>p'</var><sub><var>i</var></sub>=<var>p</var><sub><var>i</var></sub> / sum <var>p</var><sub><var>i</var></sub>.
     * 
     * @param p array of empirical counts for  first distribution 
     * @param q array of empirical counts for second distribution
     * @return KL divergence in nats b/w the two distributions 
     */
    public static double getKLDivergence(int[] p, int[] q)
    {
        assert (p.length==q.length);
        double sump=0.0;
        for (int c: p){ sump+= c;}
        double sumq=0.0;
        for (int c: q){ sumq+=c;}
        double kld = 0.0;
        for (int idx=0; idx<p.length; idx++)
        {
            double p_i = p[idx]/sump;
            double q_i = q[idx]/sumq;
            if (p_i!=0)
            {
                double r = getUnscaledLogRatio(p_i, q_i);
                kld += p_i * r;
                
                //System.out.println("#*LS.gKLD "+idx+" p "+p_i+" q "+q_i+" r "+r+" kld "+kld);
            }
        }
        return kld;
    }
    
    /**
     * Array of scores for given empirical counts.
     * 
     * @param count1 counts in hypothesis 1
     * @param count0 null hypothesis frequencies
     * @return scaled scores (by {@link #getScaleNats() scale})
     */
    public int[] getScaledScore(int[] count1, int[] count0)
    {
        assert (count1.length==count0.length);
        int sum1=0;
        for (int c: count1){ sum1+= c;}
        int sum0=0;
        for (int c: count0){ sum0+=c;}
        int[] llr = new int[count1.length];
        for (int idx=0; idx<count1.length; idx++)
        {
            double n1 = count1[idx];
            double n0 = count0[idx];
            if (n1==0.0)
            {
                llr[idx]=(n0==n1?0:Integer.MIN_VALUE);
            }
            else
            {
                double p1 = n1/sum1;
                double p0 = n0/sum0;
                double r = Math.log(p1/p0); // n1*sum0 / n0 /sum1
                double s = scaleDouble(r);
                llr[idx] = (int)s;
            }
        }
        return llr;
    }
    
    /**
     * Log-likelihood ratio in nats.
     * 
     * @param p1 likelihood for hypothesis 1 (probability, between 0 and 1)
     * @param p0 likelihood for null hypothesis (probability, between 0 and 1)
     * @return log(p1/p0)
     */
    public static double getUnscaledLogRatio(double p1, double p0)
    {
        double r = (p1<=0.0 && p0==p1)?0.0:Math.log(p1/p0);
        return r;
    }
    
    /**
     * Log-likelihood ratio on the chosen scale.
     * 
     * @param p1 likelihood for hypothesis 1 (probability, between 0 and 1)
     * @param p0 likelihood for null hypothesis (probability, between 0 and 1)
     * @return log(p1/p0)
     */
    public int getScaledLogRatio(double p1, double p0)
    {
        if (p1==0.0)
        {
            return (p0==p1?0:Integer.MIN_VALUE);
        } else
        {
                double r = Math.log(p1/p0);
                double s = scaleDouble(r);
                return (int)s;
        }
    }
    
//    /**
//     * Calculates log<sub><var>z</var></sub>(1+<var>z</var><sup><var>-x</var></sup>) for a given 
//     * log-scaled value <var>x</var>&ge;0. <var>z</var> is the scaling base 
//     * (equals e<sup>{@link #getScaleNats() }</sup>). 
//     * 
//     * 
//     * @param x log-transformed value on our scale
//     * @return log-transformed value 
//     */
//    private final int log1_plus_exp_minus_x(int x)
//    {
//        assert (x>=0);
//        final int i=x;
//        return (i<scaled_log1plusx.length)
//                ?scaled_log1plusx[i]
//                :0;
//    }
    
    /**
     * Addition for two log-scaled values.
     * 
     * @param a log <var> x</var>
     * @param b log <var> y</var>
     * @return log (<var>x</var>+<var>y</var>)
     */
    public int add(int a, int b)
    {
        int larger;
        int diff;
        if (a<b)
        {
            larger = b;
            diff = b-a;
        } else
        {
            larger = a;
            diff = a-b;
        }
        int y = (diff < scaled_log1plusx.length)
                ?scaled_log1plusx[diff]
                :0;
        return larger+y;
        // a+b = a*(1+b/a)
    }
    
    // 2^{-53} is our cutoff
    private static final double LOG_EPS = 53*Math.log(0.5);
    
}
