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

package splice.empirical;

import java.util.Arrays;
import java.util.Collection;
import java.util.List;
import java.util.Set;
import splice.annotation.AnnotatedGenomes;
import splice.annotation.GenePred;
import splice.model.FunctionMinimization;
import splice.model.Functions;
import splice.model.LogScaling;

/**
 *
 * Parametric model for intorn length: Poisson + geometric.
 * 
 * @author Mikl&oacute;s Cs&#369;r&ouml;s
 */
public class IntronLengthDistribution 
{
    private IntronLengthDistribution()
    {
    }
    private double poisson_lambda;
    private double poisson_weight;
    private double geom_fail_prob;
    private int geom_shift;

    public IntronLengthDistribution(Empirical E)
    {
        init(E);
        fit(E);
    }
    
    /**
     * Natural logarithm for the probability of a given intron length, assuming 
     * short length (Poisson) distribution.
     * 
     * @param ilen
     * @return 
     */
    public double getShortIntronLengthLL(int ilen)
    {
        double prob_p = Functions.Poisson_ln(poisson_lambda, ilen);
        return prob_p;
    }
    
    public double getLongIntronLengthLL(int ilen)
    {
        if (ilen < geom_shift)
            return Double.NEGATIVE_INFINITY;
        else 
        {
            double prob_g = Math.log1p(-geom_fail_prob)
                    +((ilen>geom_shift)
                        ?(ilen-geom_shift)*Math.log(geom_fail_prob)
                        :0.0);
            return prob_g;
        }
    }
    
    public double getShortIntronProbability()
    {
        return poisson_weight;
    }
    
    public double getShortIntronMean()
    {
        return poisson_lambda;
    }
    
    public double getLongIntronMean()
    {
        return (geom_shift-1.0)+(1.0/(1.0-geom_fail_prob));
    }
    
    public abstract class IncrementalLengthScoring
    {
        private final LogScaling scale;
    
        IncrementalLengthScoring(LogScaling scale)
        {
            this.scale = scale;
        }
        /**
         * Score for a length-1 intron
         * 
         * @return score
         */
        abstract double getIntronOpenScoreNats();
        /**
         * Score for intron length extension 
         * 
         * @param len new length
         * @return score difference for <var>len</var>-1 &rarr; <var>len</var>
         */
        abstract double getIntronAdvanceScoreNats(int len);
        public double getScaledIntronOpenScore()
        {
            return scale.scaleDouble(getIntronOpenScoreNats());
        }
        public double getScaledIntronAdvanceScore(int len)
        {
            return scale.scaleDouble(getIntronAdvanceScoreNats(len));
        }
        
        public abstract double getScaledLengthScore(int ilen);
    }
    
    public IncrementalLengthScoring getShortIntronScoring(LogScaling scale)
    {
        return new ShortIntronScoring(scale);
    }
    
    public IncrementalLengthScoring getLongIntronScoring(LogScaling scale)
    {
        return new LongIntronScoring(scale);
    }
    
    private class ShortIntronScoring extends IncrementalLengthScoring
    {
        ShortIntronScoring(LogScaling scale)
        {
            super(scale);
            this.log_lambda = Math.log(poisson_lambda);
            this.log_wt = Math.log(poisson_weight);
        }

        private final double log_lambda;
        private final double log_wt;
        
        @Override
        double getIntronOpenScoreNats() 
        {
            double log_score 
                    = log_wt // short intron prob
                    - poisson_lambda;          // e^{-l} * l^k / k! // k=0
            return log_score + getIntronAdvanceScoreNats(1); 
        }

        @Override
        double getIntronAdvanceScoreNats(int len) 
        {
            double log_diff = log_lambda - Math.log(len);
            return log_diff;
        }
        
        @Override
        public double getScaledLengthScore(int ilen)        
        {
            return getShortIntronLengthLL(ilen);
        }
    }
    
    private class LongIntronScoring extends IncrementalLengthScoring
    {
        LongIntronScoring(LogScaling scale)
        {
            super(scale);
            this.log_wt = Math.log1p(-poisson_weight);
            this.log_succ = Math.log1p(-geom_fail_prob);
            this.log_fail = Math.log(geom_fail_prob);
        }
        private final double log_wt;
        private final double log_succ;
        private final double log_fail;

        @Override
        double getIntronOpenScoreNats() 
        {
            double log_score 
                    = log_wt // long intron prob
                    + log_succ;
            return log_score;
        }

        @Override
        double getIntronAdvanceScoreNats(int len) 
        {
            return log_fail;
        }
        
        @Override 
        public double getScaledLengthScore(int ilen)
        {
            return getLongIntronLengthLL(ilen);
        }
        
    }
    
    private static final int MIN_FIT_LENGTH = 1;
    
    private void go(String org, AnnotatedGenomes annotations)
    {
        IntronLengthDistribution.Empirical E = new IntronLengthDistribution.Empirical();
        
        Collection<GenePred> genes = annotations.getAllAnnotations(org);
        E.countIntrons(annotations.getAllAnnotations(org));
        init(E);
        fit(E);
        
        System.out.println("# "+org+"\tfitted "+toString());
        
        int[] cnt = E.getAllCounts();
        int most_frequent_length=0;
        int highest_frequency=0;
        double tot_intron_length=0.0;
        int num_introns=0;
        for (int ilen=1; ilen<cnt.length; ilen++)
        {
            int n = cnt[ilen];
            double sh_ll = getShortIntronLengthLL(ilen);
            double ln_ll = getLongIntronLengthLL(ilen);
            System.out.println(org+"\t"+ilen+"\t"+n+"\t"+sh_ll+"\t"+ln_ll);
            if (n>highest_frequency)
            {
                highest_frequency=n;
                most_frequent_length=ilen;
            }
            tot_intron_length+=n*ilen;
            num_introns+=n;
        }
        
        double tot_exon_length=0.0;
        int num_exons=0;
        for (GenePred G: genes)
        {
            for (int i=0; i<G.getExonCount(); i++)
            {
                int estart = G.getExonStart(i);
                int eend = G.getExonEnd(i);
                num_exons++;
                tot_exon_length += eend-estart;
            }
        }
        tot_exon_length /= 3;
        double avg_exon_length= (num_exons==0?0.0:tot_exon_length/ num_exons);
        double avg_length = (num_introns==0?0.0:tot_intron_length/num_introns);
        double avg_cnt = num_introns/((double)genes.size());
        int ph0 = (int)(E.getPhaseCount(0)*100.0/num_introns+0.5);
        int ph1 = (int)(E.getPhaseCount(1)*100.0/num_introns+0.5);
        int ph2 = (int)(E.getPhaseCount(2)*100.0/num_introns+0.5);
        
                
        System.out.println("#TALLY\t"+org+"\tngenes "+genes.size()+"\tnintrons "+num_introns
                +" (ph1 "+E.getPhaseCount(1)+"="+ph1+"%"
                +", ph2 "+E.getPhaseCount(2)+"="+ph2+"%"
                +", ph0 "+E.getPhaseCount(0)+"="+ph0+"%"
                +")"
                +"\tintr/gene "+avg_cnt
                +"\ttot_intron_len "+tot_intron_length+"\tavg_ilen "+avg_length+"\ttypical "+most_frequent_length+"\ttot_exon_len "+tot_exon_length+"\tnexons "+num_exons+"\tavg_elen "+avg_exon_length);
        for (int prev=-1; prev<3; prev++)
        {
            int t0 = E.getTransitionCount(prev, 0);
            int t1 = E.getTransitionCount(prev, 1);
            int t2 = E.getTransitionCount(prev, 2);
            double tot = t0+t1+t2;
            int pct0 = (int) (t0*100.0/tot+0.5);
            int pct1 = (int) (t1*100.0/tot+0.5);
            int pct2 = (int) (t2*100.0/tot+0.5);
            System.out.print("#TRANSITION\t"+org+"\t"+prev);
            System.out.print("\tt0 "+t0+" ("+pct0+"%)");
            System.out.print("\tt1 "+t1+" ("+pct1+"%)");
            System.out.print("\tt2 "+t2+" ("+pct2+"%)");
            System.out.println();
        }
    }
    
    public static Empirical countIntrons(Collection<GenePred> annotations)
    {
        IntronLengthDistribution.Empirical E = new IntronLengthDistribution.Empirical();
        E.countIntrons(annotations);
        return E;
    }
    
    /**
     * Collects statistics about intron lengths in gene structure annotations. 
     * 
     */
    public static class Empirical
    {
        private final Counter<Integer> intron_count;
        private final Counter<Integer> phase_count;
        private final Counter<Integer> phase_transition_count;
        
        private Empirical()
        {
            intron_count = new Counter<>();
            phase_count = new Counter<>();
            phase_transition_count = new Counter<>();
        }
        
        /**
         * Increments the length counters for a single gene's introns.
         * 
         * @param gp 
         */
        private void countIntrons(GenePred gp)
        {
            int num_exons = gp.getExonCount();
            int cds_offset = 0;
            int prev_phase = -1;
            if (gp.isReverseStrand())
            { 
                for (int eidx=num_exons-1; eidx>0; eidx--)
                {
                    // <<<EEE..........EEE<<<
                    //       ^         ^
                    //       |         | 
                    //       ie        is
                    int is = gp.getExonStart(eidx); // exclusive
                    int ie = gp.getExonEnd(eidx-1); // inclusive
                    if (gp.getCDSStart()<=is && gp.getCDSEnd()>is)
                    {
                        int cds_delta = Math.min(gp.getCDSEnd(), gp.getExonEnd(eidx))-is;
                        assert (cds_delta>=0);
                        cds_offset += cds_delta;
                        int cur_phase = cds_offset % 3;
                        phase_count.increment(cur_phase);
                        phase_transition_count.increment(prev_phase*4+cur_phase);
                        prev_phase = cur_phase;
                    }
                    
                    int ilen = is-ie;
                    intron_count.increment(ilen);
                }
            } else
            {
                for (int eidx=1; eidx<num_exons; eidx++)
                {
                    int is = gp.getExonEnd(eidx-1); // inclusive
                    int ie = gp.getExonStart(eidx); // exclusive
                    
                    //
                    // >>>EEE.........EEE>>>
                    //       ^        ^
                    //       |        |
                    //       is       ie
                    if (gp.getCDSStart()<is && gp.getCDSEnd()>=is)
                    {
                        int cds_delta = is-Math.max(gp.getCDSStart(), gp.getExonStart(eidx-1));
                        assert (cds_delta>=0);
                        cds_offset += cds_delta;
                        int cur_phase = cds_offset % 3;
                        phase_count.increment(cur_phase);
                        phase_transition_count.increment(prev_phase*4+cur_phase);
                        prev_phase = cur_phase;
                    }
                    
                    int ilen = ie-is;
                    intron_count.increment(ilen);
                }
            }
        }
        
        /**
        /**
         * Increments the intron length counters for a set of genes.
         * 
         * @param genome 
         */
        public void countIntrons(Iterable<GenePred> genome)
        {
            for (GenePred gp: genome)
                countIntrons(gp);
        }
        
        
        /**
         * Array of observed intron legnths, shortest first.
         * 
         * @return an array that contains all observed intron lengths once, in increasing order
         */
        int[] getObservedIntronLengths()
        {
            Set<Integer> ilen_set = intron_count.keySet();
            
            int nc = ilen_set.size();
            int[] ilen_array = new int[nc];
            int arr_idx=0;
            for (int i:ilen_set)
            {
                ilen_array[arr_idx++] = i;
            }
            Arrays.sort(ilen_array);
            return ilen_array;
        }
        
        /**
         * Number of introns with given length. 
         * 
         * @param intron_length
         * @return 0 if no such introns 
         */
        int getIntronCount(int intron_length)
        {
            return intron_count.getCount(intron_length);
        }
        
        int getPhaseCount(int phase)
        {
            return phase_count.getCount(phase);
        }
        
        int getTransitionCount(int prev_phase, int next_phase)
        {
            return phase_transition_count.getCount(prev_phase*4+next_phase);
        }
                
        /**
         * Array of intron length frequencies. 
         * 
         * @return array where <var>j</var>-th element contains the number of introns with length <var>j</var>.
         */
        int[] getAllCounts()
        {
            int[] ilen_arr = getObservedIntronLengths();
            int ilen_max = ilen_arr[ilen_arr.length-1];
            int[] counts = new int[ilen_max+1];
            for (int ilen: ilen_arr)
            {
                counts[ilen]=intron_count.getCount(ilen);
            }
            return counts;
        }
    }
    
    /**
     * Rough guess for parameters {@link #poisson_lambda}, {@link #poisson_weight}, {@link #geom_fail_prob} and {@link #geom_shift}. 
     * Call before {@link #fit(splice.empirical.IntronLengthDistribution.Empirical) }.
     * 
     * @param E 
     */
    private void init(Empirical E)
    {
        int[] ilen_counts = E.getAllCounts();
        {
            int tot_ilen = 0;
            int tot_count = 0;
            int ilen_mode = 0;
            int count_before_mode = 0;
            for (int i=0; i<ilen_counts.length; i++)
            {
                int cnt = ilen_counts[i];
                if (cnt>ilen_counts[ilen_mode])
                {
                    ilen_mode = i;
                    count_before_mode = tot_count;
                }
                tot_count += cnt;
                tot_ilen += i*cnt;
            }
            double ilen_mean = ((double)tot_ilen)/((double)tot_count);
            this.poisson_lambda = ilen_mode;
            double poi_cumul = Functions.Poisson_cumulative(poisson_lambda, ilen_mode);
            double w = ((double)count_before_mode)/tot_count;
            poisson_weight = w/poi_cumul ;

            this.geom_shift = 1; //(int)poisson_lambda;
            this.geom_fail_prob = 1.0-(1.0-poisson_weight)/(ilen_mean-poisson_weight*poisson_lambda-(1.0-poisson_weight)*(geom_shift-1));

            // w/(1-f) = L
            // 1-f = w/L
            // f = 1-w/L
            
            // w*lm + (1-w)(s-1+1/(1-f)) = L
            // f = 1-(1-w)/(L-w*lm-(1-w)*(s-1))


//            System.out.println("#*ILD.init\ttot "+tot_count+"\tcumul "+poi_cumul+"\tw "+w+"\tbef "+count_before_mode+"\t// "+this);
        }
    }

    /**
     * Fits distribution parameters to observed intron length frequencies.
     * Parameter fitting is done by likelihood maximization, excluding  
     * very short and very long introns (length below {@link #MIN_FIT_LENGTH}
     * or above average + 3*sdev). 
     * 
     * @param E 
     */
    private void fit(Empirical E)
    {
        int[] ilen_counts = E.getAllCounts();

        double tot_ilen = 0.0;
        int tot_count = 0;
        double tot_ilensq = 0.0;
        for (int i=0; i<ilen_counts.length; i++)
        {
            int cnt = ilen_counts[i];
            tot_count += cnt;
            tot_ilen += i*cnt;
            tot_ilensq += i*i*cnt;
        }
        double ilen_mean = tot_ilen/tot_count;
        double ilensq_mean = tot_ilensq/tot_count;
        double ilen_var = ilensq_mean - ilen_mean*ilen_mean;
        double ilen_sd = Math.sqrt(ilen_var);
        int ilen_truncate_at = 1+(int)(ilen_mean + 3.0*ilen_sd);

        int[] trunc = new int[ilen_truncate_at];
        System.arraycopy(ilen_counts, 0, trunc, 0, ilen_truncate_at);
        ilen_counts = trunc;

        Arrays.fill(trunc, 0, MIN_FIT_LENGTH, 0);


        OptimizeParameter opt_poi = new OptimizePoisson(ilen_counts);
        OptimizeParameter opt_geo = new OptimizeGeometric(ilen_counts);
        OptimizeParameter opt_weight = new OptimizeWeight(ilen_counts);

        int REP=5;
        double tol = 1e-7;

        for (int rep=0; rep<REP; rep++)
        {
            {
                double[] zpoi = FunctionMinimization.brent(1.0, poisson_lambda, ilen_counts.length/2, opt_poi, tol);
                opt_poi.set(zpoi[0]);
                //System.out.println("#*ILD.PG.fit/"+Integer.toString(1+rep)+"\tpoi "+zpoi[0]+"\tll "+zpoi[1]+"\t// "+this.toString()+"\tll "+logLikelihood(ilen_counts));
            }
            {
                double[] zgeo = FunctionMinimization.brent(0.0, geom_fail_prob, 1.0, opt_geo, tol);
                opt_geo.set(zgeo[0]);
                //System.out.println("#*ILD.PG.fit/"+Integer.toString(1+rep)+"\tgeo "+zgeo[0]+"\tll "+zgeo[1]+"\t// "+this.toString()+"\tll "+logLikelihood(ilen_counts));
            }
            {
                double[] zwt = FunctionMinimization.brent(0.0, poisson_weight, 1.0, opt_weight, tol);
                opt_weight.set(zwt[0]);
                //System.out.println("#*ILD.PG.fit/"+Integer.toString(1+rep)+"\twt "+zwt[0]+"\tll "+zwt[1]+"\t// "+this.toString()+"\tll "+logLikelihood(ilen_counts));
            }
        }

//        System.out.println("#*ILD.fit trunc "+MIN_FIT_LENGTH+".."+ilen_truncate_at+"\tmean "+ilen_mean+"\tsdev "+ilen_sd+"\tfit "+this.toString());
    }

    @Override
    public String toString()
    {
        StringBuilder sb = new StringBuilder("ILD[");
        sb.append("poi ");
        sb.append(poisson_lambda);
        sb.append(" *");
        sb.append(poisson_weight);
        sb.append("; geo ");
        sb.append(geom_fail_prob);
        sb.append(" + ");
        sb.append(geom_shift);
        sb.append(" (avglen ");
        sb.append(geom_shift+1.0/(1.0-geom_fail_prob));
        sb.append(")");
        sb.append("]");
        return sb.toString();
    }


    public void writeData(java.io.PrintStream out)
    {
        out.println("#INTRON_LENGTH");
        out.println(this.poisson_weight+"\t"+this.poisson_lambda+"\t"+this.geom_shift+"\t"+this.geom_fail_prob);
    }
    
    public static IntronLengthDistribution readData(java.io.BufferedReader in) throws java.io.IOException
    {
        String[] fields = Aggregator.skipCommentLines(in);
        IntronLengthDistribution ILD = new IntronLengthDistribution();
        ILD.poisson_weight = Double.parseDouble(fields[0]);
        ILD.poisson_lambda = Double.parseDouble(fields[1]);
        ILD.geom_shift = Integer.parseInt(fields[2]);
        ILD.geom_fail_prob = Double.parseDouble(fields[3]);
        return ILD;
    }
    
    
//        double[] getGeometricDistribution(int n)
//        {
//            double[] d = new double[n+1];
//            if (geom_shift<=n)
//            {
//                int ilen = geom_shift;
//                double prev = d[ilen] = 1.0-geom_fail_prob;
//                while (ilen<n)
//                {
//                    ilen++;
//                    prev = d[ilen] = prev*geom_fail_prob;
//                }
//            }
//            return d;
//        }
//        
//        double[] getPoissonDistribution(int n)
//        {
//            //System.out.println("#**Poisson  r="+r);
//            double[] d = new double[n+1];
//            d[0] = Math.exp(-poisson_lambda);
//            double prev = d[0];
//            //System.out.println("#**Poisson[0]\t"+d[0]);
//            for (int j=1; j<=n; j++)
//            {
//                double f = poisson_lambda/j;
//                d[j] = prev*f;
//                //System.out.println("#**Poisson["+j+"]\t"+d[j]);
//                prev = d[j];
//            }
//            return d;
//        }        

    double logLikelihood(int[] observed_counts)
    {
        double ll=0.0;
        double prob_g = Double.NEGATIVE_INFINITY;
        double wt_p = Math.log(poisson_weight);
        double wt_g = Math.log(1.0-poisson_weight);
        double fail_log = Math.log(geom_fail_prob);

        for (int ilen=0; ilen<observed_counts.length; ilen++)
        {
            if (ilen==geom_shift)
            {
                prob_g = Math.log1p(-geom_fail_prob);
            } else if (ilen>geom_shift)
            {
                prob_g = prob_g + fail_log;

            }
            int cnt = observed_counts[ilen];
            if (cnt != 0)
            {
                double prob_p = Functions.Poisson_ln(poisson_lambda, ilen);
                double wp = wt_p+prob_p;
                double wg = wt_g+prob_g;
                double q;
                if (Double.isInfinite(prob_g))
                {
                    q = wp;
                } else if (wp<wg)
                {
                    q = wg + Math.log1p(Math.exp(wp-wg));
                } else
                {
                    q = wp + Math.log1p(Math.exp(wg-wp));
                }
                //System.out.println("#*ILD.lL "+ilen+"\t"+cnt+"\t"+q+"\tpoi "+prob_p+"\tgeo "+prob_g);
                ll += cnt * q;
            }
        }
        return ll;
    }

    private abstract class OptimizeParameter implements FunctionMinimization.OneParameterFunction
    {
        OptimizeParameter(int[] counts)
        {
            this.counts = counts;
        }

        private final int[] counts;

        abstract void set(double param);
        abstract double get();

        @Override
        public double eval(double x)
        {
            set(x);
            double ll = logLikelihood(counts);
            //System.out.println("#*ILD.OP "+this.getClass().getSimpleName()+"\t"+x+"\tll "+ll);
            return -ll;
        }
    }


    private class OptimizePoisson extends OptimizeParameter
    {
        OptimizePoisson(int[] counts)
        {
            super(counts);
        }

        @Override
        void set(double lm)
        {
            poisson_lambda = lm;
            //geom_shift = (int)lm;
        }

        @Override 
        double get()
        {
            return poisson_lambda;
        }
    }

    private class OptimizeWeight extends OptimizeParameter
    {
        OptimizeWeight(int[] counts)
        {
            super(counts);
        }

        @Override
        void set(double w)
        {
            poisson_weight = w;
        }

        @Override 
        double get()
        {
            return poisson_weight;
        }
    }

    private class OptimizeGeometric extends OptimizeParameter
    {
        OptimizeGeometric(int[] counts)
        {
            super(counts);
        }

        @Override
        void set(double f)
        {
            geom_fail_prob = f;
        }

        @Override 
        double get()
        {
            return geom_fail_prob;
        }
    }
    
    public static void main(String[] args) throws Exception 
    {
        if (args.length==0 || args[0].equals("-h"))
        {
            System.err.println("Tests: call as $0 org1=g1.fa,org2=g2.fa... <annotation file>");
            System.err.println("Reads gene annotations from file.");
            System.exit(9);
        }
        
        int arg_idx=0;
        String genome_fa_list = args[arg_idx++];
        String annotation_file = args[arg_idx++];
        AnnotatedGenomes annotations = new AnnotatedGenomes();
        List<String> wanted_organisms = annotations.readMultipleGenomes(genome_fa_list);
        annotations.readAnnotations(annotation_file);  
        
        IntronLengthDistribution ILD = new IntronLengthDistribution();
        for (String org: wanted_organisms)
        {
            ILD.go(org, annotations);
        }
    }
}