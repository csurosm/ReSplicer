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

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import splice.GeneralizedFileReader;
import splice.annotation.AnnotatedGenomes;
import splice.annotation.GenePred;
import splice.model.LogScaling;
import splice.sequence.DNASequence;
import splice.sequence.ProteinSequence;

/**
 * For each species: 
 * <ul>
 * <li>codon distribution in coding sequences (null distribution?)</li>
 * <li>donor and acceptor site composition at introns, and null distribution</li>
 * <li>codon substitutions in codon alignments</li>
 * </ul>
 * 
 * 
 * @author Mikl&oacute;s Cs&#369;r&ouml;s
 */
public class Aggregator
{
    private static final long VERSION_ID = 20151110;
    
    private static final int DEFAULT_DONOR_MOTIF_LENGTH = 4;
    private static final int DEFAULT_DONOR_BOUNDARY_LENGTH = 2;
    private static final int DEFAULT_ACCEPTOR_MOTIF_LENGTH = 3;
    private static final int DEFAULT_ACCEPTOR_BOUNDARY_LENGTH = 15;
    
    /**
     * Instantiation and computation of statistics on splice sites and codon composition.
     * 
     * @param wanted_organisms list of organisms used in the analysis
     * @param annotations 
     */
    public Aggregator(List<String> wanted_organisms, AnnotatedGenomes annotations)
    {
        this(wanted_organisms);
        initDataStructures(annotations);
    }
    
    /**
     * Instantiation without computing statistics.
     * 
     * @param wanted_organisms list of organisms used in the analysis
     */
    public Aggregator(List<String> wanted_organisms)
    {
        this.wanted_organisms = wanted_organisms;
//        this.annotations = annotations;
        this.splice_site_stats = new HashMap<>();
        this.codon_stats = new HashMap<>();
        this.intron_length = new HashMap<>();
        this.codon_trans = new CodonTransitions(wanted_organisms);
    }
    
    private void initDataStructures(AnnotatedGenomes annotations)
    {
        initSpliceSiteStats(annotations);
        initIntronLengthStats(annotations);
        initCodonStats(annotations);
    }
    
//    private final AnnotatedGenomes annotations;
    private final List<String> wanted_organisms;
    private final Map<String,SpliceSiteStatistics> splice_site_stats;
    private final Map<String,CodonStatistics> codon_stats;
    private final CodonTransitions codon_trans;
    private final LogScaling scoring_scale = new LogScaling(LogScaling.SCALE_MILLIBANS); // millibans = Phred*100
    private final Map<String, IntronLengthDistribution> intron_length;
            
    public LogScaling getScoringScale(){ return scoring_scale;}
    
    private void initSpliceSiteStats(AnnotatedGenomes annotations)
    {
        collectSpliceSiteStats(annotations,
                DEFAULT_DONOR_MOTIF_LENGTH,DEFAULT_DONOR_BOUNDARY_LENGTH,
                DEFAULT_ACCEPTOR_MOTIF_LENGTH,DEFAULT_ACCEPTOR_BOUNDARY_LENGTH);
    }
    
    public List<String> getOrganisms(){ return wanted_organisms;}
    
    // Best: 4 2 3 15
    private void collectSpliceSiteStats(AnnotatedGenomes annotations,
            int donor_motif_length,
            int donor_boundary_length,
            int acceptor_motif_length,
            int acceptor_boundary_length)
    {
        splice_site_stats.clear();
        for (String org : wanted_organisms) 
        {
            //System.out.println("#A.cSSS organism "+org);
            
            SpliceSiteStatistics stat = new SpliceSiteStatistics(donor_motif_length, donor_boundary_length, acceptor_motif_length, acceptor_boundary_length);
            stat.init(annotations.getAllAnnotations(org));
            splice_site_stats.put(org,stat);
            
            //stat.reportScores(System.out);
        }
    }
    
    private void initIntronLengthStats(AnnotatedGenomes annotations)
    {
        intron_length.clear();
        for (String org: wanted_organisms)
        {
            IntronLengthDistribution.Empirical E = IntronLengthDistribution.countIntrons(annotations.getAllAnnotations(org));
            IntronLengthDistribution ILD = new IntronLengthDistribution(E);
            intron_length.put(org, ILD);
        }
    }
    
    private void initCodonStats(AnnotatedGenomes annotations)
    {
        codon_stats.clear();
        for (String org: wanted_organisms)
        {
            //System.out.println("#A.iCS organism "+org);
            CodonStatistics stat = new CodonStatistics();
            stat.init(annotations.getAllAnnotations(org));
            codon_stats.put(org, stat);
            //getCodonScoring(org).reportScores(System.out);
            
        }
    }
    
    public void countAlignment(AnnotatedGenomes annotations, ProteinSequence[] alignment) 
    {
        this.codon_trans.countAlignment(annotations, alignment);
    }
    
    
    
    private static class CodonStatistics
    {
        final CodonCounter codon_cnt;
        final CodonCounter codon_cnt_null;
        
        CodonStatistics(CodonCounter codon_cnt, CodonCounter codon_null)
        {
            this.codon_cnt = codon_cnt;
            this.codon_cnt_null = codon_null;
        }
        
        CodonStatistics()
        {
            codon_cnt = new CodonCounter();
            codon_cnt_null = new CodonCounter();
        }
        
        void init(Iterable<GenePred> gene_set)
        {
            codon_cnt.countExons(gene_set);
            codon_cnt_null.countIntrons(gene_set);
        }
    }
    
    private Map<String, CodonCounter.Scoring> memoized_codon_scoring;
    
    public CodonCounter.Scoring getCodonScoring(String org)
    {
        if (memoized_codon_scoring == null)
            memoized_codon_scoring = new HashMap<>();

        if (memoized_codon_scoring.containsKey(org))
            return memoized_codon_scoring.get(org);
        else
        {
            CodonStatistics stats = codon_stats.get(org);
            CodonCounter.Scoring sc = stats.codon_cnt.getScoring(scoring_scale, stats.codon_cnt_null);

            memoized_codon_scoring.put(org, sc);
            return sc;
        }
    }
    
    private Map<String, Map<String, CodonTransitions.Scoring>> memoized_codon_subst_scoring;
    
    public CodonTransitions.Scoring getCodonSubstScoring(String org_from, String org_to)
    {
        if (memoized_codon_subst_scoring==null)
            memoized_codon_subst_scoring = new HashMap<>();
        if (!memoized_codon_subst_scoring.containsKey(org_from))
            memoized_codon_subst_scoring.put(org_from, new HashMap<>());

        if (memoized_codon_subst_scoring.get(org_from).containsKey(org_to))
            return memoized_codon_subst_scoring.get(org_from).get(org_to);
        else
        {
            CodonStatistics stats_from = codon_stats.get(org_from);
            CodonStatistics stats_to = codon_stats.get(org_to);

            CodonTransitions.Scoring sc = codon_trans.getScoring(org_from, org_to, scoring_scale, stats_from.codon_cnt_null, stats_to.codon_cnt_null);
            memoized_codon_subst_scoring.get(org_from).put(org_to, sc);
            return sc;
        }
    }
    
    private Map<String, IntronLengthDistribution.IncrementalLengthScoring> memoized_short_intron_length_scoring;
    public IntronLengthDistribution.IncrementalLengthScoring getShortIntronLengthScoring(String org)
    {
        if (memoized_short_intron_length_scoring==null)
            memoized_short_intron_length_scoring = new HashMap<>();
        if (memoized_short_intron_length_scoring.containsKey(org))
            return memoized_short_intron_length_scoring.get(org);
        else
        {
            IntronLengthDistribution ILD = intron_length.get(org);
            IntronLengthDistribution.IncrementalLengthScoring sc = ILD.getShortIntronScoring(scoring_scale);
            memoized_short_intron_length_scoring.put(org, sc);
            return sc;
        }
    }
    
    public IntronLengthDistribution getIntronLengthDistribution(String org)
    {
        return intron_length.get(org);
    }

    private static class SpliceSiteStatistics
    {
        final SpliceSiteComposition ss5;
        final SpliceSiteComposition ss3;
        final SpliceSiteComposition ss5_null;
        final SpliceSiteComposition ss3_null;
        
        
        SpliceSiteStatistics(int ss5_motif_length, int ss5_boundary, int ss3_motif_length, int ss3_boundary)
        {
            ss5 = new SpliceSiteComposition(ss5_motif_length, ss5_boundary, true);
            ss5_null = new SpliceSiteComposition(ss5_motif_length, ss5_boundary, true);
            ss3 = new SpliceSiteComposition(ss3_motif_length, ss3_boundary,false);
            ss3_null = new SpliceSiteComposition(ss3_motif_length, ss3_boundary,false);
        }
        
        SpliceSiteStatistics(SpliceSiteComposition ss5, SpliceSiteComposition ss5_null, SpliceSiteComposition ss3, SpliceSiteComposition ss3_null)
        {
            this.ss5 = ss5;
            this.ss5_null = ss5_null;
            this.ss3 = ss3;
            this.ss3_null = ss3_null;
        }
        
        void init(Iterable<GenePred> gene_set)
        {
            for (GenePred gene: gene_set)
            {
                ss5.countIntronSites(gene);
                ss5_null.countNonSites(gene, false);
                ss3.countIntronSites(gene);
                ss3_null.countNonSites(gene, false);
            }
        }
        
        void reportScores(java.io.PrintStream out)
        {
            LogScaling scale = new LogScaling(LogScaling.SCALE_MILLIBANS);
            int[] ss5_mot1 = ss5.getAllMotifCounts();
            int[] ss5_mot0 = ss5_null.getAllMotifCounts();
            int[] ss5_mot_score = scale.getScaledScore(ss5_mot1, ss5_mot0);
            double ss5_mot_kl = LogScaling.getKLDivergence(ss5_mot1, ss5_mot0);
            
            int[][] ss5_bnd1 = ss5.getBoundaryCounts();
            int[][] ss5_bnd0 = ss5_null.getBoundaryCounts();
            
            int[][] ss5_bnd_score = new int[ss5_bnd1.length][];
            double[] ss5_bnd_kl = new double[ss5_bnd_score.length];
            
            double kl5 = ss5_mot_kl;
            for (int j=0; j<ss5_bnd_score.length; j++)
            {
                ss5_bnd_score[j] = scale.getScaledScore(ss5_bnd1[j], ss5_bnd0[j]);
                ss5_bnd_kl[j] = LogScaling.getKLDivergence(ss5_bnd1[j], ss5_bnd0[j]);
                kl5 += ss5_bnd_kl[j];
            }
            int[] ss3_mot1 = ss3.getAllMotifCounts();
            int[] ss3_mot0 = ss3_null.getAllMotifCounts();
            int[] ss3_mot_score = scale.getScaledScore(ss3_mot1, ss3_mot0);
            double ss3_mot_kl = LogScaling.getKLDivergence(ss3_mot1, ss3_mot0);
            
            int[][] ss3_bnd1 = ss3.getBoundaryCounts();
            int[][] ss3_bnd0 = ss3_null.getBoundaryCounts();
            
            int[][] ss3_bnd_score = new int[ss3_bnd1.length][];
            double[] ss3_bnd_kl = new double[ss3_bnd_score.length];
            double kl3 = ss3_mot_kl;
            for (int j=0; j<ss3_bnd_score.length; j++)
            {
                ss3_bnd_score[j] = scale.getScaledScore(ss3_bnd1[j], ss3_bnd0[j]);
                ss3_bnd_kl[j] = LogScaling.getKLDivergence(ss3_bnd1[j], ss3_bnd0[j]);
                kl3 += ss3_bnd_kl[j];
            }
            
            for (int j=0; j<ss5_mot_score.length; j++)
            {
                if (ss5_mot1[j]>0)
                {
                    out.println("#SPLICE 5'mot "+ss5.decodeMotif(j)+"\t"+ss5_mot1[j]+"\t"+ss5_mot_score[j]);
                }
            }
            out.println("#SPLICE 5'mot KL "+ss5_mot_kl);
            for (int bidx=0; bidx<ss5_bnd_score.length; bidx++)
            {
                out.print("#SPLICE 5'++"+bidx);
                for (int nuc=0; nuc<4; nuc++)
                {
                    out.print("\t"+DNASequence.toChar(nuc)+" "+ss5_bnd1[bidx][nuc]+"/"+ss5_bnd_score[bidx][nuc]);
                }
                out.print("\tKL "+ss5_bnd_kl[bidx]);
                out.println();
            }
            
            for (int j=0; j<ss3_mot_score.length; j++)
            {
                if (ss3_mot1[j]>0)
                {
                    out.println("#SPLICE 3'mot "+DNASequence.reverseComplement(ss3.decodeMotif(j))+"\t"+ss3_mot1[j]+"\t"+ss3_mot_score[j]);
                }
            }
            out.println("#SPLICE 3'mot KL "+ss3_mot_kl);
            for (int bidx=0; bidx<ss3_bnd_score.length; bidx++)
            {
                out.print("#SPLICE 3'--"+bidx);
                for (int nuc=0; nuc<4; nuc++)
                {
                    out.print("\t"+DNASequence.toChar(nuc)+" "+ss3_bnd1[bidx][nuc]+"/"+ss3_bnd_score[bidx][nuc]);
                }
                out.print("\tKL "+ss3_bnd_kl[bidx]);
                out.println();
            }
            
            kl5 /= LogScaling.SCALE_MILLIBANS;
            kl3 /= LogScaling.SCALE_MILLIBANS;
            out.println("#SPLICE 5'KLmbans "+kl5+"\t3'KLmbans "+kl3);
            
        }
    }
    
    public static final String GENERIC_SPLICE_SITE = "this is not a real genome";
    private Map<String, SpliceSiteComposition.Scoring> memoized_donor_site_scoring;
    /**
     * Scoring system appropriate for a genome
     * 
     * @param org use {@link #GENERIC_SPLICE_SITE} for generic scoring (canonical, non-canonical and U12 splicing motifs)
     * @return scoring scheme
     */
    public SpliceSiteComposition.Scoring getDonorSiteScoring(String org)
    {
        if (memoized_donor_site_scoring==null)
            memoized_donor_site_scoring = new HashMap<>();
        if (memoized_donor_site_scoring.containsKey(org))
            return memoized_donor_site_scoring.get(org);
        else
        {
            SpliceSiteComposition.Scoring sc ;
            if (GENERIC_SPLICE_SITE.equals(org))
            {
                sc = SpliceSiteComposition.genericDonorScoring(scoring_scale);
            } else
            {  
                SpliceSiteStatistics stats = splice_site_stats.get(org);
                sc = stats.ss5.getScoring(scoring_scale, stats.ss5_null);
            }
            memoized_donor_site_scoring.put(org, sc);
            return sc;
        }
    }
    
    private Map<String, SpliceSiteComposition.Scoring> memoized_acceptor_site_scoring;
    
    public SpliceSiteComposition.Scoring getAcceptorSiteScoring(String org)
    {
        if (memoized_acceptor_site_scoring==null)
            memoized_acceptor_site_scoring = new HashMap<>();
        if (memoized_acceptor_site_scoring.containsKey(org))
            return memoized_acceptor_site_scoring.get(org);
        else
        {
            SpliceSiteComposition.Scoring sc;
            if (GENERIC_SPLICE_SITE.equals(org))
            {
                sc = SpliceSiteComposition.genericAcceptorScoring(scoring_scale);
            } else
            {  
                SpliceSiteStatistics stats =splice_site_stats.get(org);
                sc = stats.ss3.getScoring(scoring_scale, stats.ss3_null);
            }
            memoized_acceptor_site_scoring.put(org, sc);
            return sc;
        }
    }
    
    static void writeOrganisms(java.io.PrintStream out, List<String> wanted_organisms)
    {
        out.println("#ORGANISMS");
        for (int oidx=0; oidx<wanted_organisms.size(); oidx++)
        {
            String org = wanted_organisms.get(oidx);
            if (oidx>0)
                out.print("\t");
            out.print(org);
        }
        out.println();
    }
    
    static List<String> readOrganisms(java.io.BufferedReader in) throws java.io.IOException
    {
        String[] fields = skipCommentLines(in);
        List<String> listed_org = new ArrayList<>();
        for (String o: fields)
            listed_org.add(o);
        return listed_org;
        
    }
    
    public void writeData(java.io.PrintStream out)
    {
        out.println("#"+getClass().getName()+"\t"+VERSION_ID);
        writeOrganisms(out, wanted_organisms);
        out.println("#SUBST");
        codon_trans.writeData(out);
        
        for (String org: wanted_organisms)
        {
            SpliceSiteStatistics S = splice_site_stats.get(org);
            out.println("#SS5 "+org);
            S.ss5.writeData(out);
            out.println("#SS5/null "+org);
            S.ss5_null.writeData(out);
            out.println("#SS3 "+org);
            S.ss3.writeData(out);
            out.println("#SS3/null "+org);
            S.ss3_null.writeData(out);

            IntronLengthDistribution ILD = intron_length.get(org);
            out.println("#ILEN "+org);
            ILD.writeData(out);
            
            CodonStatistics C = codon_stats.get(org);
            out.println("#CODON "+org);
            C.codon_cnt.writeData(out);
            out.println("#CODON/null "+org);
            C.codon_cnt_null.writeData(out);
        }
        out.println("#END");
    }
    
    public static String[] skipCommentLines(java.io.BufferedReader in) throws java.io.IOException
    {
        String line;
        do
        {
            line = in.readLine();
        } while (line != null && line.startsWith("#"));
        return (line==null?null:line.split("\\t"));
    }

    public static Aggregator readData(java.io.BufferedReader in) throws java.io.IOException
    {
        List<String> listed_org = readOrganisms(in);
        Aggregator agg = new Aggregator(listed_org);
        agg.codon_trans.readData(in);
        
        for (int oidx=0; oidx<listed_org.size(); oidx++)
        {
            String org = listed_org.get(oidx);
            SpliceSiteComposition ss5 = SpliceSiteComposition.readData(in);
            SpliceSiteComposition ss5_null = SpliceSiteComposition.readData(in);
            SpliceSiteComposition ss3 = SpliceSiteComposition.readData(in);
            SpliceSiteComposition ss3_null = SpliceSiteComposition.readData(in);
            SpliceSiteStatistics S = new SpliceSiteStatistics(ss5, ss5_null, ss3, ss3_null);

            IntronLengthDistribution ILD = IntronLengthDistribution.readData(in);
            
            CodonCounter codon_cnt = CodonCounter.readData(in);
            CodonCounter codon_cnt_null = CodonCounter.readData(in);
            CodonStatistics C = new CodonStatistics(codon_cnt, codon_cnt_null);
            agg.splice_site_stats.put(org, S);
            agg.codon_stats.put(org, C);
            agg.intron_length.put(org, ILD);
        }
        return agg;
    }
    
    public static void main(String[] args) throws Exception
    {
        if (args.length==0 || args[0].equals("-h")){
            System.err.println("Tests: call as $0 org1=g1.fa,org2=g2.fa annot.txt ali1.fa ali2.fa ...");
            System.err.println("Outputs splice site compositions.");
            System.exit(9);
        }
        
        if (args.length == 1)
        {
            java.io.BufferedReader R = new GeneralizedFileReader(args[0]);
            Aggregator agg = Aggregator.readData(R);
            agg.writeData(System.out);
            R.close();
        } else
        {
            int arg_idx=0;
            String genome_fa_list = args[arg_idx++];
            String annotation_file = args[arg_idx++];
    //        int donor_motif_length = Integer.parseInt(args[arg_idx++]);
    //        int donor_boundary_length = Integer.parseInt(args[arg_idx++]);
    //        int acceptor_motif_length = Integer.parseInt(args[arg_idx++]);
    //        int acceptor_boundary_length = Integer.parseInt(args[arg_idx++]);

            AnnotatedGenomes annotations = new AnnotatedGenomes();
            List<String> wanted_organisms = annotations.readMultipleGenomes(genome_fa_list);
            annotations.readAnnotations(annotation_file);  

            Aggregator agg = new Aggregator(wanted_organisms, annotations);
            
            
            while (arg_idx<args.length)
            {
                String alignment_file = args[arg_idx++];
                agg.codon_trans.countAlignment(annotations, ProteinSequence.readMultiFasta(alignment_file));
            }
            agg.writeData(System.out);
        } 
                
    }    
}
