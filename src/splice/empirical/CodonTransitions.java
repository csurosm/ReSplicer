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
import java.util.Arrays;
import java.util.List;
import splice.GeneralizedFileReader;
import splice.annotation.AnnotatedGenomes;
import splice.annotation.GenePred;
import splice.model.CodonSubstitutionModel;
import splice.model.LogScaling;
import splice.sequence.AminoAcid;
import splice.sequence.Codon;
import splice.sequence.ProteinSequence;
import splice.sequence.TranslationTable;

/**
 *
 * @author Mikl&oacute;s Cs&#369;r&ouml;s
 */
public class CodonTransitions 
{
    private CodonTransitions()
    {
        this.wanted_organisms = new ArrayList<>();
    }
    
    public CodonTransitions(List<String> wanted_organisms)
    {
        this.wanted_organisms = wanted_organisms;
    }
    
    private final List<String> wanted_organisms;
    private final Counter<CodonColumn> column_counter = new Counter<>();
    
    private Codon[][] getCodonSequences(ProteinSequence[] aligned_sequences, AnnotatedGenomes annotations)
    {
        Codon[][] aligned_codons = new Codon[aligned_sequences.length][];
        for (int seq_idx=0; seq_idx<aligned_sequences.length; seq_idx++)
        {
            ProteinSequence seq_aa = aligned_sequences[seq_idx];
            GenePred gp = annotations.getAnnotationForGene(seq_aa);
            
            int seqlen_aa = seq_aa.getLength();
            Codon[] seq_cod = aligned_codons[seq_idx] = new Codon[seqlen_aa];
            String cds = gp.getSequenceCDS();
            
            {
                int pos_aa=0;
                int pos_nt=0;
                
                while (pos_aa<seqlen_aa)
                {
                    AminoAcid aa = seq_aa.getResidueAt(pos_aa);
//                    if (pos_nt+3>cds.length())
//                    {
//                        System.out.println("#**CT.dA different lengths "+gp.shortDesc(true)+"\tcdslen "+cds.length()+"\taalen "+seqlen_aa+"\tscaf "+gp.getChrom().getLength());
//                        System.out.println(seq_aa.toFasta());
//                        System.out.println(cds);
//                        System.out.flush();
//                    }
                    
                    if (aa.isRegularResidue() && pos_nt+3<=cds.length())
                    {
                        char nuc1 = cds.charAt(pos_nt++);
                        char nuc2 = cds.charAt(pos_nt++);
                        char nuc3 = cds.charAt(pos_nt++);
                        Codon codon = Codon.getCodon(nuc1, nuc2, nuc3);
                        seq_cod[pos_aa] = codon;
                    } else
                    {
                        seq_cod[pos_aa] = null;
                    }
                    pos_aa++;
                }
            }
        }
        return aligned_codons;
    }
    
    /**
     * Sorts the protein sequences by the given organism order.
     * 
     * @param wanted_organisms desired order
     * @param aligned_aa aligned sequences, one per organism
     * @param annotations 
     */
    private void sortSequences(ProteinSequence[] aligned_aa, List<String> wanted_organisms, AnnotatedGenomes annotations)
    {
        ProteinSequence[] retval = new ProteinSequence[aligned_aa.length];
        for (int seq_idx=0; seq_idx<aligned_aa.length; seq_idx++)
        {
            ProteinSequence seq_aa = aligned_aa[seq_idx];
            GenePred gene = annotations.getAnnotationForGene(seq_aa);
            String organism = gene.getChrom().getOrganism(); // safer than looking for organism tag in the protein seq defline
            int org_idx = wanted_organisms.indexOf(organism);
            retval[org_idx] = seq_aa;
        }
        for (int seq_idx=0; seq_idx<aligned_aa.length; seq_idx++)
        {
            aligned_aa[seq_idx] = retval[seq_idx];
        }
    }
    
    private static String encode(CodonColumn col)
    {
        char[] retval = new char[col.codons.length];
        for (int cod_idx=0; cod_idx<col.codons.length; cod_idx++)
        {
            Codon cod = col.codons[cod_idx];
            char codon_code = (cod==null?CODON_CHAR_OFFSET:(char)(CODON_CHAR_OFFSET+1+cod.getIndex()));
            retval[cod_idx] = codon_code;
        }
        return new String(retval);
    }
    
    private static CodonColumn decode(String column_code)
    {
        Codon[] codons = new Codon[column_code.length()];
        for (int cod_idx=0; cod_idx<codons.length; cod_idx++)
        {
            char c = column_code.charAt(cod_idx);
            if (c==CODON_CHAR_OFFSET)
            {
                codons[cod_idx] = null;
            } else
            {
                int codon_idx = c-1-CODON_CHAR_OFFSET;
                codons[cod_idx] = Codon.getCodon(codon_idx);
            }
        }
        return new CodonColumn(codons);
    }
    
    private static class CodonColumn implements Comparable<CodonColumn>
    {
        CodonColumn(Codon[] codons)
        {
            this.codons = codons;
        }
        
        private Codon[] codons;
        
        @Override
        public int hashCode()
        {
            return Arrays.hashCode(codons);
        }
        
        @Override
        public boolean equals(Object o)
        {
            if (o instanceof CodonColumn)
            {
                CodonColumn cc = (CodonColumn) o;
                return Arrays.equals(codons, cc.codons); // handles null appropriately
            } else
                return super.equals(o); // false presumably
        }
        
        @Override
        public String toString()
        {
            StringBuilder sb = new StringBuilder('{');
            for (int j=0; j<codons.length; j++)
            {
                if (j>0) sb.append(',');
                if (codons[j]==null)
                    sb.append("---/-");
                else
                    sb.append(codons[j].toString(TranslationTable.STANDARD_TABLE));
            }
            sb.append('}');
            return sb.toString();
        }
        
        int getCodonIdx(int org_idx)
        { 
            Codon cod = codons[org_idx];
            if (cod==null)
                return Codon.NUM_CODONS;
            else
                return cod.getIndex();
        }

        @Override
        public int compareTo(CodonColumn o) 
        {
            int cmp =0;
            for (int col_idx=0; col_idx<codons.length && cmp==0; col_idx++)
            {
                if (codons[col_idx]==null)
                {
                    cmp = (o.codons[col_idx]==null?0:-1);
                } else
                {
                    cmp = (o.codons[col_idx]==null?1:codons[col_idx].compareTo(o.codons[col_idx]));
                }
            }
            return cmp;
        }
    }
    
    /**
     * Selects one column of aligned codons.
     * @param aligned_codons array of aligned codon columns (rows = sequences)
     * @param column_idx column index
     * @return column constructed from <var>aligned_codons</var>[*][<var>column_idx</var>]
     */
    private CodonColumn newColumn(Codon[][] aligned_codons, int column_idx)
    {
        Codon[] col = new Codon[aligned_codons.length];
        for (int seq_idx=0; seq_idx<aligned_codons.length; seq_idx++)
        {
            col[seq_idx]  = aligned_codons[seq_idx][column_idx];
        }
        
        return new CodonColumn(col);
    }
    
    /**
     * Array of aligned codon pairs in two sequences.
     * 
     * @param org_idx_row organism index in {@link CodonColumn} for rows
     * @param org_idx_col organism index in {@link CodonColumn} for columns
     * @return a 65*65 array, last row/column corresponds to gap
     */
    private int[][] getPairwiseCounts(int org_idx_row, int org_idx_col)
    {
        int[][] codon_trans_counts = new int[Codon.NUM_CODONS+1][Codon.NUM_CODONS+1];
        for (CodonColumn col: column_counter.keySet())
        {
            int codon_from_idx = col.getCodonIdx(org_idx_row);
            int codon_to_idx = col.getCodonIdx(org_idx_col);
            int cnt = column_counter.getCount(col);
            codon_trans_counts[codon_from_idx][codon_to_idx] += cnt;// stop codons are counted here as regular ones ...
        }
        return codon_trans_counts;
    }
    
    private static final double CODON_SUBST_PSEUDOCOUNT = 1.0;
    private static final TranslationTable GENETIC_CODE = TranslationTable.STANDARD_TABLE;
    
    private class Substitution implements CodonSubstitutionModel
    {
        Substitution(String org_from, String org_to)
        {
            idx_from=wanted_organisms.indexOf(org_from);
            idx_to = wanted_organisms.indexOf(org_to);
            this.codon_to_frequencies = new double[Codon.NUM_CODONS];
            this.codon_subst_probs = new double[Codon.NUM_CODONS][Codon.NUM_CODONS];
            computeProbabilities();
        }
        private final int idx_from;
        private final int idx_to;
        private final double[] codon_to_frequencies;
        private final double[][] codon_subst_probs;
        private void computeProbabilities() 
        {
            int[][] codon_subst_counts = getPairwiseCounts(idx_from, idx_to);
            
            int[] codon_to_counts = new int[Codon.NUM_CODONS];
            double codon_to_tot = 0.0;
            for (int j=0; j<Codon.NUM_CODONS; j++)
                if (!GENETIC_CODE.isStopCodon(Codon.getCodon(j)))
                {
                    for (int i=0; i<Codon.NUM_CODONS; i++)
                        if (!GENETIC_CODE.isStopCodon(Codon.getCodon(i)))
                        {
                            codon_to_counts[j] += codon_subst_counts[i][j]; 
                        }
                    codon_to_tot += codon_to_counts[j];
                }
            for (int j=0; j<Codon.NUM_CODONS; j++)
            {
                codon_to_frequencies[j] = codon_to_counts[j]/codon_to_tot;
            }
            
            for (int i=0; i<Codon.NUM_CODONS; i++)
                if (!GENETIC_CODE.isStopCodon(Codon.getCodon(i)))
                {
                    double row_tot = 0.0;
                    for (int j=0; j<Codon.NUM_CODONS; j++)
                        if (!GENETIC_CODE.isStopCodon(Codon.getCodon(j)))
                        {
                            double cnt = codon_subst_counts[i][j] + CODON_SUBST_PSEUDOCOUNT*codon_to_frequencies[j];
                            codon_subst_probs[i][j] = cnt;
                            row_tot += codon_subst_counts[i][j];
                        }
                    row_tot += CODON_SUBST_PSEUDOCOUNT;
                    for (int j=0; j<Codon.NUM_CODONS; j++)
                    {
                        codon_subst_probs[i][j] /= row_tot;
                    }
                }
            
            for (int i=0; i<Codon.NUM_CODONS; i++)
                if (GENETIC_CODE.isStopCodon(Codon.getCodon(i)))
                    codon_subst_probs[i][i] = 1.0;

            System.out.println("#*CT.S.cP transitions ============= "+wanted_organisms.get(idx_from)+" -> "+wanted_organisms.get(idx_to));
            for (int codon_from_idx=0; codon_from_idx<Codon.NUM_CODONS; codon_from_idx++)
            {
                double row_tot = 0.0;
                for (int codon_to_idx=0; codon_to_idx<Codon.NUM_CODONS; codon_to_idx++)
                {
                    int cnt = codon_subst_counts[codon_from_idx][codon_to_idx];
                    double p = codon_subst_probs[codon_from_idx][codon_to_idx];
                    System.out.print("#*CT.S.cP "+Codon.getCodon(codon_from_idx).toString(TranslationTable.STANDARD_TABLE)
                                +"\t-> "
                                +Codon.getCodon(codon_to_idx).toString(TranslationTable.STANDARD_TABLE));
                    System.out.printf("\t%.4f\t%d\n", p, cnt);
                    row_tot += p;
                }
                System.out.println("#*CT.S.cP "+Codon.getCodon(codon_from_idx)+"\t-> *** "+row_tot);
            }

            System.out.println("#*CT.S.cP frequencies ============= "+wanted_organisms.get(idx_from)+" -> "+wanted_organisms.get(idx_to));
            double col_tot = 0.0;
            for (int codon_to_idx=0; codon_to_idx<Codon.NUM_CODONS; codon_to_idx++)
            {
                double p = codon_to_frequencies[codon_to_idx];
                System.out.println("#*CT.S.cP "+Codon.getCodon(codon_to_idx).toString(TranslationTable.STANDARD_TABLE)+"\t"+p);
                col_tot += p;
            }
            System.out.println("#*CT.S.cP total "+col_tot);
            
            
        }
        
        @Override
        public double getStationaryProbability(Codon cod)
        {
            return codon_to_frequencies[cod.getIndex()];
        }
        @Override
        public double getTransitionProbability(Codon codon_from, Codon codon_to)
        {
            return codon_subst_probs[codon_from.getIndex()][codon_to.getIndex()];
        }
        
    }
    
    public Scoring getScoring(String org_from, String org_to, LogScaling scoring_scale, 
                CodonCounter from_null,
                CodonCounter to_null)
    {
        return new Scoring(org_from, org_to, scoring_scale, from_null, to_null);
    }
    
    // homology prob: p(x,y) = p(x) * p(y|x)
    // null prob: q(x)*q(y)
    // p(x) * p(y|x) / q(x) / q(y) = p(x)/q(x) * p(y|x) / q(y)
    
    public class Scoring 
    {
        private Scoring(String org_from, String org_to, LogScaling scoring_scale, 
                CodonCounter from_null,
                CodonCounter to_null)
        {
            int org_from_idx = wanted_organisms.indexOf(org_from);
            int org_to_idx = wanted_organisms.indexOf(org_to);
            
            int[][] codon_trans_cnt = getPairwiseCounts(org_from_idx, org_to_idx);
            int[] codon_from_cnt = new int[Codon.NUM_CODONS];
            for (int i=0; i<Codon.NUM_CODONS; i++)
                if (!GENETIC_CODE.isStopCodon(Codon.getCodon(i)))
                    for (int j=0; j<Codon.NUM_CODONS; j++)
                        if (!GENETIC_CODE.isStopCodon(Codon.getCodon(j)))
                        {
                            codon_from_cnt[i] += codon_trans_cnt[i][j];
                        }
            
            codon_score_from = scoring_scale.getScaledScore(codon_from_cnt, from_null.getCombinedCodonCounts());
            
            codon_score_to = new int[Codon.NUM_CODONS][];
            int[] codon_to_null = to_null.getCombinedCodonCounts();
            double tot_counts = 0.0;
            double exp_score = 0.0;
            for (int i=0; i<Codon.NUM_CODONS; i++)
            {
                if (!GENETIC_CODE.isStopCodon(Codon.getCodon(i)))
                {
                    int[] subst_cnt = new int[Codon.NUM_CODONS];
                    for (int j=0; j<Codon.NUM_CODONS; j++)
                        if (!GENETIC_CODE.isStopCodon(Codon.getCodon(j)))
                        {
                            subst_cnt[j] = codon_trans_cnt[i][j];
                            tot_counts += subst_cnt[j];
                        }
                    codon_score_to[i] = scoring_scale.getScaledScore(subst_cnt, codon_to_null);
                    for (int j=0; j<Codon.NUM_CODONS; j++)
                        if (!GENETIC_CODE.isStopCodon(Codon.getCodon(j)))
                        {
                            if (codon_score_to[i][j]!=Integer.MIN_VALUE)
                            {
                                codon_score_to[i][j]+=codon_score_from[i];
                                exp_score += codon_trans_cnt[i][j]*codon_score_to[i][j];
                            }
                        }
                }
            }
            this.expected_codon_score = exp_score / tot_counts;
            
        }
        
        private final int[] codon_score_from;
        private final int[][] codon_score_to;
        
        private final double expected_codon_score;
        
        public int getScore(Codon codon_from, Codon codon_to)
        {
            if (GENETIC_CODE.isStopCodon(codon_from))
            {
                return GENETIC_CODE.isStopCodon(codon_to)?0:Integer.MIN_VALUE;
            } else
            {
                return (GENETIC_CODE.isStopCodon(codon_to))
                        ?Integer.MIN_VALUE
                        :codon_score_to[codon_from.getIndex()][codon_to.getIndex()];
            }
        }
        
        /**
         * Expected score of aligned codon pairs
         * 
         * @return expected alignment score per codon position
         */
        public double getExpectedScore()
        {
            return this.expected_codon_score;
        }
        
        public void reportScores(java.io.PrintStream out)
        {
            out.print("#codon\tscore0");
            for (int i=0; i<Codon.NUM_CODONS; i++)
            {
                out.print("\t ");
                out.print(Codon.getCodon(i).toString(GENETIC_CODE));
            }
            out.println();
            for (int i=0; i<Codon.NUM_CODONS; i++)
            {
                out.print(Codon.getCodon(i).toString(GENETIC_CODE));
                out.print(" \t");
                out.printf("%6d", codon_score_from[i]);
                for (int j=0; j<Codon.NUM_CODONS; j++)
                {
                    out.printf("\t%6d", getScore(Codon.getCodon(i), Codon.getCodon(j)));
                }
                out.println();
            }
            out.println("#codon\texpected: "+this.expected_codon_score);
        }
        
    }
    

    private static final char CODON_CHAR_OFFSET = '!';
    
    public void writeData(java.io.PrintStream out)
    {
//        out.println("#SUBST");
//        out.println(wanted_organisms.size()+"\t# number of organisms");
//        for (int oidx=0; oidx<wanted_organisms.size(); oidx++)
//        {
//            if (oidx>0)
//                out.print('\t');
//            out.print(wanted_organisms.get(oidx));
//            
//        }
//        out.println("\t# organism order");
        CodonColumn[] all_columns = column_counter.keySet().toArray(new CodonColumn[0]);
        out.println(all_columns.length+"\t# number of alignment column types");
        Arrays.sort(all_columns);
        for (CodonColumn col:all_columns)
        {
            out.print("|\t");
            out.print(encode(col));
            out.print("\t");
            out.println(column_counter.getCount(col));
        }
    }
    
//    /**
//     * Writes the content to an output stream (more economical than standard Serializable).
//     * 
//     * @param out
//     * @throws java.io.IOException 
//     */
//    private void writeDataOOS(java.io.ObjectOutputStream out)
//         throws java.io.IOException
//    {
//        out.writeInt(wanted_organisms.size());
//        for (String org: wanted_organisms)
//            out.writeObject(org);
//        Set<CodonColumn> all_columns = column_counter.keySet();
//        out.writeInt(all_columns.size());
//        for (CodonColumn col: all_columns)
//        {
//            out.writeObject(encode(col));
//            out.writeInt(column_counter.getCount(col));
//        }
//    }

//    private static String[] skipCommentLines(java.io.BufferedReader in) throws java.io.IOException
//    {
//        String line;
//        do
//        {
//            line = in.readLine();
//        } while (line != null && line.startsWith("#"));
//        return (line==null?null:line.split("\\t"));
//    }
    

    public void readData(java.io.BufferedReader in) throws java.io.IOException
    {
//        System.out.print("#*CT.rD ");
//        Aggregator.writeOrganisms(System.out, wanted_organisms);
        String[] fields = Aggregator.skipCommentLines(in);
        int num_col = Integer.parseInt(fields[0]);
        column_counter.clear();
        for (int col_idx=0; col_idx<num_col; col_idx++)
        {
            fields = Aggregator.skipCommentLines(in);
            CodonColumn col = decode(fields[1]);
            int cnt = Integer.parseInt(fields[2]);
            column_counter.set(col, cnt);
        }
    }
    
//    /**
//     * Reinitializes using a previous state saved with {@link #writeDataOOS(java.io.ObjectOutputStream) }.
//     * @param in
//     * @throws java.io.IOException
//     * @throws ClassNotFoundException 
//     */
//    private void readDataOIS(java.io.ObjectInputStream in)
//        throws java.io.IOException, ClassNotFoundException
//    {
//        wanted_organisms.clear();
//        int num_org = in.readInt();
//        for (int oidx=0; oidx<num_org; oidx++)
//        {
//            String org = (String) in.readObject();
//            wanted_organisms.add(org);
//        }
//        
//        column_counter.clear();
//        int num_col = in.readInt();
//        for (int col_idx=0; col_idx<num_col; col_idx++)
//        {
//            String column_code = (String) in.readObject();
//            CodonColumn col = decode(column_code);
//            int cnt = in.readInt();
//            column_counter.set(col, cnt);
//        }
//    }
    
    public void countAlignment(AnnotatedGenomes annotations, ProteinSequence[] alignment) 
    {
        ProteinSequence[] aligned_aa;
        {
            List<ProteinSequence> L = new ArrayList<>();
            for (ProteinSequence PS: alignment)
            {
                if (annotations.getAnnotationForGene(PS)!=null)
                    L.add(PS);
            }
            aligned_aa = L.toArray(new ProteinSequence[0]);
            sortSequences(aligned_aa, wanted_organisms, annotations);
        }
        Codon[][] aligned_codons = getCodonSequences(aligned_aa, annotations);

//            System.out.println("#*CT.dA ali "+alignment_file+"\tnseq "+aligned_codons.length+"\talilen "+aligned_codons[0].length);
//            for (int sidx=0; sidx<aligned_codons.length; sidx++)
//            {
//                System.out.println("#*CT.dA ali "+alignment_file+"\t\tsi "+sidx+"\tseq "+aligned_aa[sidx].getDefLine());
//            }

        for (int col_idx=0; col_idx<aligned_codons[0].length; col_idx++)
        {
            CodonColumn col = newColumn(aligned_codons, col_idx);
            column_counter.increment(col);
        }        
    }
    
    
    /**
     * Collects statistics about column occurrences.
     * 
     * @param args command-line arguments
     * @throws Exception 
     */
    private void uzsgyi(String[] args) throws Exception
    {
        if (args.length==0 || args[0].equals("-h")){
            System.err.println("Tests: call as "+getClass().getName()+" org1=g1.fa,org2=g2.fa annot.txt outfile ali1.fa ali2.fa ...");
            System.err.println("Outputs aligned codon frequencies.");
            System.exit(9);
        }
        
        int arg_idx=0;
        String genome_fa_list = args[arg_idx++];
        String annotation_file = args[arg_idx++];
        String outfile = args[arg_idx++];
//        java.io.ObjectOutputStream oos = new java.io.ObjectOutputStream(new java.io.FileOutputStream(outfile));

        AnnotatedGenomes annotations = new AnnotatedGenomes();
        wanted_organisms.addAll(annotations.readMultipleGenomes(genome_fa_list));
        annotations.readAnnotations(annotation_file);
        
        while (arg_idx<args.length)
        {
            String alignment_file = args[arg_idx];
            countAlignment(annotations, ProteinSequence.readMultiFasta(alignment_file));
            ++arg_idx;
        }

        java.io.PrintStream out = ("-".equals(outfile)?System.out:new java.io.PrintStream(outfile));
        Aggregator.writeOrganisms(out, wanted_organisms);
        writeData(out);
        out.close();
//        writeDataOOS(oos);
//        oos.close();
        
//        if (arg_idx<args.length)
//        {
//            Substitution subst = new Substitution("hyal", "albu");
//        }
        
    }
    
    private void printCounts()
    {
        CodonColumn[] all_columns = column_counter.toArray(new CodonColumn[0]);
                
//                = column_counter.keySet().toArray(new CodonColumn[0]);
//        Arrays.sort(all_columns, new Comparator<CodonColumn>() {
//            @Override
//            public int compare(CodonColumn o1, CodonColumn o2) 
//            {
//                return Integer.compare(column_counter.getCount(o2), column_counter.getCount(o1));
//            }
//        });
        
        int num_multiple_occ=0;
        int tot_len = 0;
        for (CodonColumn col: all_columns)
        {
            int cnt = column_counter.getCount(col);
            if (cnt>1)
            {
                num_multiple_occ++;
                System.out.println("COUNTS "+cnt+"\t"+col);
            }
            tot_len += cnt;
        }
        System.out.println("COUNTS ---- total length: "+tot_len+"\tcolumn types: "+all_columns.length+"\t(multi,single) "+num_multiple_occ+", "+(all_columns.length-num_multiple_occ));
    }
    
    public static void main(String[] args) throws Exception
    {
        CodonTransitions CT;
        if (args.length>1)
        {
            CT = new CodonTransitions();
            CT.uzsgyi(args);
        } else
        {
            java.io.BufferedReader in = new GeneralizedFileReader(args[0]);
            CT = new CodonTransitions(Aggregator.readOrganisms(in));
            CT.readData(in);
            in.close();
//            java.io.ObjectInputStream ois = new java.io.ObjectInputStream(new java.io.FileInputStream(args[0]));
//            CT.readDataOIS(ois);
//            ois.close();
        }        
        //CT.printCounts();
    }
    
}
