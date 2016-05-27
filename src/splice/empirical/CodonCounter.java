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

import java.io.PrintStream;
import splice.annotation.AnnotatedGenomes;
import splice.annotation.GenePred;
import splice.model.LogScaling;
import splice.sequence.Codon;
import splice.sequence.DNASequence;
import splice.sequence.TranslationTable;

/**
 *
 * @author Mikl&oacute;s Cs&#369;r&ouml;s
 */
public class CodonCounter 
{
    private final int[] first_codons;
    private final int[] middle_codons;
    private final int[] last_codons;
    
    public CodonCounter()
    {
        first_codons = new int[Codon.NUM_CODONS];
        middle_codons = new int[Codon.NUM_CODONS];
        last_codons = new int[Codon.NUM_CODONS];
    }
    
    /**
     * Array of first, middle and last codon counts.
     * 
     * @return 3x64 size array of codon counts
     */
    public int[][] getAllCodonCounts()
    {
        int[][] retval = new int[3][];
        retval[0]=getFirstCodonCounts();
        retval[1] = getMiddleCodonCounts();
        retval[2] = getLastCodonCounts();
        
        return retval;
    }

    /**
     * Sum of codon counts in all positions.
     * 
     * @return 64-element array of codon counts
     */
    public int[] getCombinedCodonCounts()
    {
        int[][] all_cnt = getAllCodonCounts();
        for (int cidx=0; cidx<Codon.NUM_CODONS; cidx++)
        {
            all_cnt[0][cidx] += all_cnt[1][cidx]+all_cnt[2][cidx];
        }
        return all_cnt[0];
    }
    
    public int[] getFirstCodonCounts()
    {
        return (int[])first_codons.clone();
    }

    public int[] getMiddleCodonCounts()
    {
        return (int[])middle_codons.clone();
    }

    public int[] getLastCodonCounts()
    {
        return (int[])last_codons.clone();
    }
    
    /**
     * Counts the combined occurrence of nucleotides in first, second, and third codon positions.
     * @param codon_counts 64-element array of codon counts
     * @return 4-element array of nucleotide counts, indexed  by IDX_ constants from {@link splice.sequence.NucleotideEncoding}
     */
    public static int[] getNucleotideCounts(int[] codon_counts)
    {
        int[] nuc_counts = new int[4];
        
        for (int codon_idx=0;codon_idx<codon_counts.length; codon_idx++)
        {
            Codon codon = Codon.getCodon(codon_idx);
            int n= codon_counts[codon_idx];
            nuc_counts[DNASequence.toIndex(codon.nucleotideAt(1))]+=n;
            nuc_counts[DNASequence.toIndex(codon.nucleotideAt(2))]+=n;
            nuc_counts[DNASequence.toIndex(codon.nucleotideAt(3))]+=n;
        }
        return nuc_counts;
    }
    
    /**
     * Counts the codon distribution in exons for a set of genes.
     * 
     * @param genes set of gene structures
     */
    public void countExons(Iterable<GenePred> genes)
    {
        for (GenePred gp: genes)
        {
            String gs = gp.getSequenceCDS();
            boolean truncated_start;
            boolean truncated_end;
            if (gp.isReverseStrand())
            {
                truncated_start = (gp.getExonEnd(gp.getExonCount()-1)==gp.getChrom().getLength());
                truncated_end = (gp.getExonStart(0)==0);
            } else
            {
                truncated_start = (gp.getExonStart(0)==0);
                truncated_end = (gp.getExonEnd(gp.getExonCount()-1)==gp.getChrom().getLength());
            }
//            if (truncated_start || truncated_end)
//                System.out.println("#*CDS.E.cE truncated gene structure annotation "+gp.shortDesc());
            countSequence(gs, truncated_start, truncated_end);
        }
    }

    /**
     * Calculates a null distribution for a set of genes, using the 
     * triplets observed in annotated introns.
     * 
     * @param genes set of gene structures
     */
    public void countIntrons(Iterable<GenePred> genes)
    {
        for (GenePred gp:genes)
        {
            int num_exons = gp.getExonCount();
            for (int eidx=1; eidx<num_exons; eidx++)
            {
                String iseq = gp.getSequenceIntron(eidx-1);
                countSequence(iseq, false, false);
                countSequence(DNASequence.reverseComplement(iseq), false, false);
            }
        }
    }    

    private void countSequence(String gs, boolean truncated_start, boolean truncated_end)
    {
        char[] gc = gs.toCharArray();
        int glen = gc.length;
        if (glen<3)
            return;
        int[] current_codon_set = (truncated_start?middle_codons:first_codons);
        int cidx = 0;
        do
        {
            if (3*cidx+2<glen)
            {
                char c1=gc[3*cidx];
                char c2=gc[3*cidx+1];
                char c3=gc[3*cidx+2];
                if (!DNASequence.isAmbiguous(c1)
                    && !DNASequence.isAmbiguous(c2)
                        && !DNASequence.isAmbiguous(c3))
                {
                    Codon C = new Codon(c1,c2,c3);
                    current_codon_set[C.getIndex()]++;
                    current_codon_set = middle_codons;
                }
            } 
            if (3*cidx+3>=glen && !truncated_end)
            {
                char c1=gc[glen-3];
                char c2=gc[glen-2];
                char c3=gc[glen-1];
                if (!DNASequence.isAmbiguous(c1)
                    && !DNASequence.isAmbiguous(c2)
                        && !DNASequence.isAmbiguous(c3))
                {
                    Codon C = new Codon(c1,c2,c3);
                    last_codons[C.getIndex()]++;
                }
            }
            ++cidx;
        } while (3*cidx<glen);
    }        
    
    public Scoring getScoring(LogScaling scoring_scale, CodonCounter codon_null)
    {
        return new Scoring(scoring_scale, codon_null);
    }
    
    public class Scoring
    {
        private final int[] first_codon_score;
        private final int[] middle_codon_score;
        private final int[] last_codon_score;
        private Scoring(LogScaling scoring_scale, CodonCounter codon_null)
        {
            int[] cnt_null = codon_null.getCombinedCodonCounts();
            first_codon_score = scoring_scale.getScaledScore(first_codons, cnt_null);
            middle_codon_score = scoring_scale.getScaledScore(middle_codons, cnt_null);
            last_codon_score = scoring_scale.getScaledScore(last_codons, cnt_null);
        }
        
        public int getFirstCodonScore(Codon cod)
        {
            return first_codon_score[cod.getIndex()];
        }
        
        public int getMiddleCodonScore(Codon cod)
        {
            return middle_codon_score[cod.getIndex()];
        }
        
        public int getLastCodonScore(Codon cod)
        {
            return last_codon_score[cod.getIndex()];
        }
        
        public void reportScores(java.io.PrintStream out)
        {
            for (Codon cod: Codon.ALL_CODONS)
            {
                out.println("#CODON\t"+cod.toString(TranslationTable.STANDARD_TABLE)
                        +"\t"+getFirstCodonScore(cod)
                        +"\t"+getMiddleCodonScore(cod)
                        +"\t"+getLastCodonScore(cod));
            }
        }
    }
    
    
    private static void reportCounts(Iterable<GenePred> genes, PrintStream out)
    {
        CodonCounter counter_e = new CodonCounter();
        counter_e.countExons(genes);
        CodonCounter counter_i = new CodonCounter();
        counter_i.countIntrons(genes);
        
        LogScaling scale = new LogScaling(LogScaling.SCALE_MILLIBANS);
        
        double[] tot_e = new double[3];
        double tot_i = 0.0;
        
        int[][] num_e = counter_e.getAllCodonCounts();
        int[] num_i = counter_i.getMiddleCodonCounts();
        
        for (int cidx=0; cidx<Codon.NUM_CODONS; cidx++)
        {
            for (int pos_idx=0; pos_idx<3; pos_idx++)
            {
                tot_e[pos_idx] += num_e[pos_idx][cidx];
            }
            tot_i += num_i[cidx];
        }
        
        for (int cidx=0; cidx<Codon.NUM_CODONS; cidx++)
        {
            Codon codon = Codon.getCodon(cidx);
            double[] p_e = new double[3];
            double p_i = num_i[cidx]/tot_i;
            
            for (int pos_idx=0; pos_idx<3; pos_idx++)
            {
                p_e[pos_idx] = num_e[pos_idx][cidx]/tot_e[pos_idx];
            }
            
            out.printf("%s\t%d\t%d\t%d\t%d\t%.2f%%\t%.2f%%\t%.2f%%\t%.2f%%\t%d\t%d\t%d\n", 
                    codon.toString(TranslationTable.STANDARD_TABLE),
                    num_e[0][cidx], num_e[1][cidx], num_e[2][cidx],
                    num_i[cidx],
                    100.0*p_e[0], 100.0*p_e[1], 100.0*p_e[2], 100.0*p_i,
                    scale.getScaledLogRatio(p_e[0], p_i),
                    scale.getScaledLogRatio(p_e[1], p_i),
                    scale.getScaledLogRatio(p_e[2], p_i)
                    );
        }
        
        out.printf("# Expected score for codons in coding sequence wrt null :\t%f\t%f\t%f\n",
                scale.scaleDouble(LogScaling.getKLDivergence(num_e[0], num_i)),
                scale.scaleDouble(LogScaling.getKLDivergence(num_e[1],num_i)),
                scale.scaleDouble(LogScaling.getKLDivergence(num_e[2], num_i))
                );
    }
    
    public void writeData(java.io.PrintStream out)
    {
        out.println("#CODON");
        for (int cidx=0; cidx<Codon.NUM_CODONS; cidx++)
        {
            out.println(first_codons[cidx]+"\t"+middle_codons[cidx]+"\t"+last_codons[cidx]
                +"\t# "+Codon.getCodon(cidx).toString(TranslationTable.STANDARD_TABLE));
        }
    }

//    private static String[] skipCommentLines(java.io.BufferedReader in) throws java.io.IOException
//    {
//        String line;
//        do
//        {
//            line = in.readLine();
//        } while (line != null && line.startsWith("#"));
//        return (line==null?null:line.split("\\t"));
//    }

    public static CodonCounter readData(java.io.BufferedReader in) throws java.io.IOException
    {
        CodonCounter cc = new CodonCounter();
        for (int cidx=0; cidx<Codon.NUM_CODONS; cidx++)
        {
            String[] fields = Aggregator.skipCommentLines(in);
            cc.first_codons[cidx] = Integer.parseInt(fields[0]);
            cc.middle_codons[cidx] = Integer.parseInt(fields[1]);
            cc.last_codons[cidx] = Integer.parseInt(fields[2]);
        }
        return cc;
    }
    
    /**
     * Reads a genome sequence and a set of gene annotations; outputs 
     * the empirical codon counts.
     * 
     * @param args command-line arguments
     * @throws Exception if something goes wrong
     */
    public static void main(String[] args) throws Exception
    {
        if (args.length==0 || args[0].equals("-h")){
            System.err.println("Tests: call as $0 genome.fa annot.txt organism");
            System.err.println("Reads genome, and gene annotations for a given organism, outputs empirical codon counts");
            System.exit(9);
        }
        int arg_idx=0;
        String genome_fa = args[arg_idx++];
        String annotation_file = args[arg_idx++];
        String target_org = args[arg_idx++];
        
        AnnotatedGenomes A = new AnnotatedGenomes();
        DNASequence[] refseq = DNASequence.readMultiFasta(genome_fa);
        A.addGenome(target_org, refseq);
        A.readAnnotations(annotation_file);
        
        reportCounts(A.getAllAnnotations(target_org), System.out);
    }
}

