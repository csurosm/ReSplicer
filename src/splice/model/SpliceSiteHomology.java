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

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import splice.GeneralizedFileReader;
import splice.annotation.AnnotatedGenomes;
import splice.annotation.GenePred;
import splice.empirical.Aggregator;
import splice.empirical.CodonCounter;
import splice.empirical.CodonTransitions;
import splice.empirical.Counter;
import splice.empirical.SpliceSiteComposition;
import splice.model.SpliceSiteHistory;
import splice.sequence.AminoAcid;
import splice.sequence.Codon;
import splice.sequence.DNASequence;
import splice.sequence.DNAStrand;
import splice.sequence.MultiHomology;
import splice.sequence.ProteinSequence;
import splice.sequence.TranslationTable;

/**
 *
 * @author Mikl&oacute;s Cs&#369;r&ouml;s
 */
public class SpliceSiteHomology 
{
    private static int NUM_MATCHING_NUCLEOTIDES = 4;
    
    private static boolean DEBUG_TRACK_BRACKETING_ALIGNMENT = false;
    private static boolean DEBUG_TRACK_BRACKETING_BOUNDARIES = false;
    
    private static final String OPTION_NUM_MATCHING = "-anchorlen";
    private static final String OPTION_TRACK_ALIGNMENT = "-showali";
    private static final String OPTION_TRACK_BOUNDARY = "-showbrackets";
    
    private SpliceSiteHomology()
    {
        
    }
    
    public SpliceSiteHomology(Aggregator agg)
    {
        this();
        setScoring(agg);
    }
    
    
    private void setScoring(Aggregator agg){ this.agg=agg;}
    private Aggregator agg;
    
    private ProteinSequence[] aligned_sequences;
    private GenePred[] homolog_genes;
    private DNAStrand[] transcript_strands;
    
    /**
     * <var>aligned_genomic_positions</var>: for each gene, gives the vector of the aligned codons' first genomic position. 
     * <var>aligned_genomic_positions</var>[<var>i</var>][<var>j</var>] is the first nucleotide's position for the query_codon in 
     * sequence <var>i</var> in alignment column <var>j</var>. If the alignment column contains a gap 
     * in that column, then the value is {@link #ALIGNMENT_GAP}.
     */
    private int[][] aligned_genomic_positions;
    /**
     * <var>aligned_codons</var>: for each gene, gives the vector of aligned codons. 
     * <var>aligned_codons</var>[<var>i</var>][<var>j</var>] is the query_codon in 
 sequence <var>i</var> in alignment column <var>j</var>. If the alignment has a gap there, or if 
     * a covered nucleotide is ambiguous, the value is <code>null</code>.
     */
    private Codon[][] aligned_codons;
    /**
     * <var>intron_columns</var>: for each gene, gives the intron columns.
     * <var>intron_columns</var>[<var>i</var>][<var>j</var>] is the column where intron <var>j</var> falls in 
     * when projected onto the alignment's sequence <var>i</var>. <var>j</var>=0 is the first intron 
     * in the direction of the transcription, <var>j</var>=1 is the second etc.
     */
    private int[][] intron_columns;
    
    /**
     * Intron phase mapped into an alignment column: for an alignment column, gives 
     * organism &rarr; phase mapping <var>ph</var>=0,1,2; intron falls after <var>ph</var> nucleotides of 
 the column's query_codon. 
     */
    private Map<Integer, Map<String,Integer>> projected_intron_phase;
    
    public static final int ALIGNMENT_GAP = -1;
    
    
    public void reportGenes(java.io.PrintStream out)
    {
        for (int gene_idx=0; gene_idx<aligned_sequences.length; gene_idx++)
        {
            ProteinSequence P = aligned_sequences[gene_idx];
            GenePred G = homolog_genes[gene_idx];
            out.println("#GENE\t"+gene_idx+"\torganism "+P.getOrganism()+"\t"+G.shortDesc(true)+"\tdefline "+P.getDefLine());
        }
    }
    /**
     * Initializes 
     * {@link #aligned_sequences}, 
     * {@link #homolog_genes} and 
     * {@link #transcript_strands}.
     * 
     * @param annotations
     * @param aligned_sequences
     */
    private void initDataStructures(AnnotatedGenomes annotations, ProteinSequence[] all_aligned_sequences)
    {
        {
            List<ProteinSequence> kept_sequences=new ArrayList<>();
            for (ProteinSequence pseq: all_aligned_sequences)
            {
                if (annotations.getAnnotationForGene(pseq)!=null) kept_sequences.add(pseq);
            }
            this.aligned_sequences = kept_sequences.toArray(new ProteinSequence[0]);
        }
        
        
        this.homolog_genes = new GenePred[aligned_sequences.length];
        this.transcript_strands = new DNAStrand[aligned_sequences.length];
        
        for (int gene_idx=0; gene_idx<aligned_sequences.length; gene_idx++)
        {
            ProteinSequence P = aligned_sequences[gene_idx];
            GenePred G = annotations.getAnnotationForGene(P);
            homolog_genes[gene_idx] = G;
            transcript_strands[gene_idx] = new DNAStrand(G.getChrom(), G.getStrand());
            
        }
        
        aligned_genomic_positions = new int[aligned_sequences.length][];
        aligned_codons = new Codon[aligned_sequences.length][];
        intron_columns = new int[aligned_sequences.length][];
        projected_intron_phase=new HashMap<>();        
    }
    
    /**
     * Alignment length. (Works after {@link #initDataStructures(splice.annotation.AnnotatedGenomes, splice.sequence.ProteinSequence[]) } was called.)
     * 
     * @return number of alignment columns
     */
    public int getNumAlignmentColumns()
    {
        if (aligned_sequences.length==0) return 0;
        else return aligned_sequences[0].getLength();
    }
    
    /**
     * Genomic coordinate of the aligned codon's first nucleotide.
     * 
     * @param seq_idx sequnce index
     * @param column alignment column
     * @return {@link #ALIGNMENT_GAP} if a gap, otherwise the position of the first aligned nucleotide
     */
    public int getAlignedGenomicPosition(int seq_idx, int column)
    {
        return this.aligned_genomic_positions[seq_idx][column];
    }
    
    /**
     * Combines gene structure data and protein alignment. 
     * Sets {@link #aligned_genomic_positions}, {@link #aligned_codons}, {@link #intron_columns},
     * {@link #projected_intron_phase}.
     * 
     * @param annotations genome annotations
     * @param all_aligned_sequences aligned sequences in arbitrary order
     * @return false if two introns fall within a single codon
     */
    public boolean computeAlignmentColumns(AnnotatedGenomes annotations, ProteinSequence[] all_aligned_sequences)
    {
        initDataStructures(annotations, all_aligned_sequences);
        
        int aligned_seqlen = getNumAlignmentColumns(); 

        assert (aligned_genomic_positions.length == aligned_sequences.length);
        assert (aligned_codons.length == aligned_sequences.length);
        assert (intron_columns.length == aligned_sequences.length);
        
        for (int seq_idx=0; seq_idx<aligned_sequences.length; seq_idx++)
        {
            aligned_genomic_positions[seq_idx] = new int[aligned_seqlen];
            aligned_codons[seq_idx] = new Codon[aligned_seqlen];
            intron_columns[seq_idx] = new int[homolog_genes[seq_idx].getExonCount()-1];
            GenePred G = homolog_genes[seq_idx];
            DNASequence chrom = G.getChrom();
            int chrom_length = chrom.getLength();
            int pos = G.isReverseStrand()
                    ? G.getCDSEnd()-1
                    : G.getCDSStart();
            int eidx=0;
            while (eidx<G.getExonCount() && G.getExonEnd(eidx)<=pos) eidx++;

            byte[] codon_nuc = new byte[3];
            
            for (int column_idx=0; column_idx<aligned_seqlen; column_idx++)
            {
                AminoAcid aa = aligned_sequences[seq_idx].getResidueAt(column_idx);
                
                if (AminoAcid.GAP.equals(aa))
                {
                    aligned_genomic_positions[seq_idx][column_idx] = ALIGNMENT_GAP;
                    aligned_codons[seq_idx][column_idx] = null;
                } else
                {   
                    aligned_genomic_positions[seq_idx][column_idx] = pos;
                    if (G.isReverseStrand())
                    {
                        for (int nuc_idx=0; nuc_idx<3; nuc_idx++)
                        {
                            if (pos<0)
                                break;
                            codon_nuc[nuc_idx] = DNASequence.complement(chrom.getNucleotideAt(pos));
                            if (pos==G.getExonStart(eidx) && eidx!=0)
                            {
                                // e0 e1 e2 ; exoncount=3
                                intron_columns[seq_idx][G.getExonCount()-1-eidx] = column_idx;
                                Map<String, Integer> intron_phase_pos = projected_intron_phase.get(column_idx);
                                if (intron_phase_pos==null)
                                    intron_phase_pos = new HashMap<>();
                                if (intron_phase_pos.containsKey(G.getName()))
                                {
                                    System.err.println("### gene "+G.shortDesc()
                                            +": two introns within same codon (phase "+intron_phase_pos.get(G.getName())
                                            +", "+nuc_idx+") in column "+column_idx
                                            +", exon "+eidx);
//                                    throw new java.lang.RuntimeException("### gene "+G.shortDesc()
//                                            +": two introns within same codon (phase "+intron_phase_pos.get(G.getName())
//                                            +", "+nuc_idx+") in column "+column_idx
//                                            +", exon "+eidx);
                                    return false;
                                }
                                intron_phase_pos.put(G.getName(), 1+nuc_idx);
                                projected_intron_phase.put(column_idx, intron_phase_pos);
                                --eidx;
                                pos = G.getExonEnd(eidx);
                            }
                            --pos;
                        }
                    } else
                    {
                        for (int nuc_idx=0; nuc_idx<3; nuc_idx++)
                        {
                            if (pos>=chrom_length)
                                break;
                            codon_nuc[nuc_idx] = chrom.getNucleotideAt(pos);
                            pos++;
                            if (pos==G.getExonEnd(eidx) && eidx!=G.getExonCount()-1)
                            {
                                intron_columns[seq_idx][eidx] = column_idx;
                                Map<String, Integer> intron_phase_pos = projected_intron_phase.get(column_idx);
                                if (intron_phase_pos==null)
                                    intron_phase_pos = new HashMap<>();
                                if (intron_phase_pos.containsKey(G.getName()))
                                {
                                    System.err.println("### gene "+G.shortDesc()
                                            +": two introns within same codon (phase "+intron_phase_pos.get(G.getName())
                                            +", "+nuc_idx+") in column "+column_idx
                                            +", exon "+eidx);
                                    return false;
//                                    throw new java.lang.RuntimeException("### gene "+G.shortDesc()
//                                            +": two introns within same codon (phase "+intron_phase_pos.get(G.getName())
//                                            +", "+nuc_idx+") in column "+column_idx
//                                            +", exon "+eidx);
                                }
                                intron_phase_pos.put(G.getName(), 1+nuc_idx);
                                projected_intron_phase.put(column_idx, intron_phase_pos);
                                eidx++;
                                pos = G.getExonStart(eidx);
                            }
                        }
                    }
                    if (DNASequence.isAmbiguous(codon_nuc[0])
                            || DNASequence.isAmbiguous(codon_nuc[1])
                            || DNASequence.isAmbiguous(codon_nuc[2]))
                        aligned_codons[seq_idx][column_idx] = null;
                    else
                        aligned_codons[seq_idx][column_idx] = Codon.getCodon(codon_nuc[0], codon_nuc[1], codon_nuc[2]);
                }
            }
        }
        return true;
    }
    
    public void printAlignment(java.io.PrintStream out)
    {
        int aligned_seqlen = this.getNumAlignmentColumns();
        out.print("#ALIGNMENT\tcolumn"); // header
        for (int seq_idx=0; seq_idx<aligned_sequences.length; seq_idx++)
        {
            String org = homolog_genes[seq_idx].getChrom().getOrganism();
            out.print("\t");
            out.print(org);
        }
        out.println();
        for (int col_idx=0; col_idx<aligned_seqlen; col_idx++)
        {
            out.print("#ALIGNMENT\t"+col_idx);
            for (int seq_idx=0; seq_idx<aligned_sequences.length; seq_idx++)
            {
                out.print("\t");
                out.print(describeAlignedCodonPosition(seq_idx, col_idx));
            }
            out.println();
            for (int seq_idx=0; seq_idx<homolog_genes.length; seq_idx++)
                if (getIntronPhase(homolog_genes[seq_idx], col_idx)!=NO_INTRON)
                {
                    out.println("#INTRON\t"+getOrganism(seq_idx)+"\t"+getIntronSequence(seq_idx, col_idx));
                }
        }
        
        for (int seq_idx=0; seq_idx<aligned_sequences.length; seq_idx++)
        {
            out.println("#PROJECTION\t"+homolog_genes[seq_idx].getChrom().getOrganism()+":\tintron columns\t"+Arrays.toString(intron_columns[seq_idx]));
        }
    }
    
    
    public static final int NO_INTRON = -1;
    
    /**
     * Intron phase.
     * 
     * @param gene one of the aligned genes
     * @param column_idx alignment column
     * @return {@link #NO_INTRON} if no intron there, or intron phase (1,2,3)
     */
    public int getIntronPhase(GenePred gene, int column_idx)
    {
        if (projected_intron_phase.containsKey(column_idx))
        {
            return projected_intron_phase.get(column_idx).getOrDefault(gene.getName(), NO_INTRON);
        } else
            return NO_INTRON;
    }
    public Integer[] getIntronColumns()
    {
        Integer[] columns_with_introns = projected_intron_phase.keySet().toArray(new Integer[0]);
        Arrays.sort(columns_with_introns);        
        return columns_with_introns;
    }
    
    public String getOrganism(int seq_idx)
    {
        return homolog_genes[seq_idx].getChrom().getOrganism();
    }
    
    public int getNumSequences()
    {
        return homolog_genes.length;
    }
    
    public GenePred getGene(int seq_idx)
    {
        return homolog_genes[seq_idx];
    }
    
    public DNAStrand getTranscriptStrand(int seq_idx)
    {
        return transcript_strands[seq_idx];
    }
    
    private final class BracketedHomology
    {
        BracketedHomology(int query_idx, int ref_idx, int intron_column)
        {
            this();
            this.left_bracket_column = intron_column;
            this.right_bracket_column = intron_column;
            calculateBracketing(query_idx, ref_idx, intron_column, true);
            calculateBracketing(query_idx, ref_idx, intron_column, false);
        }
        
        private BracketedHomology()
        {
            this.segments = new int[homolog_genes.length][];
            this.left_bracket_score = Integer.MIN_VALUE;
            this.right_bracket_score = Integer.MIN_VALUE;
        }
        
        private int left_bracket_score;
        private int left_bracket_column;
        private int right_bracket_score;
        private int right_bracket_column;
        private final int[][] segments;

        boolean hasValidBrackets()
        {
            return left_bracket_score > 0 && right_bracket_score > 0;
        }
        
        @Override
        public boolean equals(Object o)
        {
            if (o instanceof BracketedHomology)
            {
                BracketedHomology other = (BracketedHomology) o;
                return (left_bracket_score > 0 == other.left_bracket_score>0)
                        && (left_bracket_column == other.left_bracket_column)
                        && (right_bracket_score>0 == other.right_bracket_score>0)
                        && (right_bracket_column == other.right_bracket_column);
            } else
                return super.equals(o); // false    
        }

        @Override
        public int hashCode()
        {
            int hash = (this.left_bracket_score>0?1:0)
                    + 2*(this.right_bracket_score>0?2:0);
            hash = 37 * hash + this.left_bracket_column;
            hash = 37 * hash + this.right_bracket_column;
            return hash;
        }

        private void calculateBracketing(int query_idx, int ref_idx, int intron_column, boolean scan_upstream)
        {
            GenePred query_gene = homolog_genes[query_idx];
            GenePred ref_gene = homolog_genes[ref_idx];
            CodonTransitions.Scoring subst_scoring = agg.getCodonSubstScoring(getOrganism(query_idx), getOrganism(ref_idx));

            int max_score_threshold = (int)(NUM_MATCHING_NUCLEOTIDES*subst_scoring.getExpectedScore());
            int query_phase = getIntronPhase(query_gene, intron_column);

            int scanned_column_idx = intron_column; // used in scanning
            if (scan_upstream)
            {
                if (query_phase == 3) 
                {
                    scanned_column_idx++;
                }
            } else // downstream
            {
                scanned_column_idx++;
            }

            int best_start_col = -1;
            int best_end_col = -1;
            int best_score = 0;

            if (scan_upstream)
            {
                int score = 0;
                int end_col = scanned_column_idx;
                int start_col = scanned_column_idx;

                while(start_col>0)
                {
                    start_col--;
                    if (aligned_genomic_positions[query_idx][start_col]==ALIGNMENT_GAP)
                    {
                        if (aligned_genomic_positions[ref_idx][start_col]==ALIGNMENT_GAP)
                        {
                            // both have gaps here
                            // nothing to do
                        } else
                        {
                            // query has a gap but not ref
                            score = 0;
                        }
                    } else // query has no gap here
                    {
                        if (aligned_genomic_positions[ref_idx][start_col]==ALIGNMENT_GAP)
                        {
                            // ref has a gap but query has none
                            score = 0;
                        } else
                        {
                            // check if they have introns here
                            query_phase = getIntronPhase(query_gene, start_col);
                            int ref_phase = getIntronPhase(ref_gene, start_col);

                            boolean query_nointron = query_phase == -1 || (query_phase == 3 && end_col==start_col+1);
                            boolean ref_nointron = ref_phase == -1 || (ref_phase == 3 && end_col==start_col+1);
                            // get codons, get subst score
                            if (query_nointron && ref_nointron)
                            {
                                Codon query_codon = aligned_codons[query_idx][start_col];
                                Codon ref_codon = aligned_codons[ref_idx][start_col];
                                if (query_codon != null && ref_codon != null) // exclude codons with ambiguous nucleotides
                                {
                                    score += subst_scoring.getScore(query_codon, ref_codon);
                                }
                            } else
                            {
                                score = 0;
                            }
                        }
                    }
                    if (score <= 0)
                    {
                        end_col = start_col;
                        score = 0;
                    } else if (score>best_score)
                    {
                        best_score = score;
                        best_start_col =start_col;
                        best_end_col = end_col;
                        if (best_score>max_score_threshold)
                            break;
                    }
                } // while start_codon>0

                if (DEBUG_TRACK_BRACKETING_ALIGNMENT)
                {
                    if (DEBUG_TRACK_BRACKETING_BOUNDARIES)
                    {
                        System.out.println("#*SSH.IC.cC col/up "+scanned_column_idx+"\t"+getOrganism(query_idx)+"/"+getOrganism(ref_idx)+"\tsc "+best_score
                            +((best_score>max_score_threshold)?"*":"")
                            +"\t"+best_start_col+".."+best_end_col
                            +"\t// "+max_score_threshold);
                    }
                    if (best_score>0)
                    {
                        for (int c=Math.max(best_end_col-12,best_start_col); c<best_end_col; c++)
                        {
                            int sc =  (aligned_codons[query_idx][c]!=null && aligned_codons[ref_idx][c]!=null)
                                    ?subst_scoring.getScore(aligned_codons[query_idx][c], aligned_codons[ref_idx][c]):0;
                            System.out.println("#*SSH.IC.cC "+c+"|\t"+describeAlignedCodonPosition(query_idx, c)+"\t"+describeAlignedCodonPosition(ref_idx, c)+"\t"+sc);
                        }
                        for (int c=best_end_col; c<=scanned_column_idx; c++)
                        {
                            int sc =  (aligned_codons[query_idx][c]!=null && aligned_codons[ref_idx][c]!=null)
                                    ?subst_scoring.getScore(aligned_codons[query_idx][c], aligned_codons[ref_idx][c]):0;
                            System.out.println("#*SSH.IC.cC "+c+".\t"+describeAlignedCodonPosition(query_idx, c).toLowerCase()+"\t"+describeAlignedCodonPosition(ref_idx, c).toLowerCase()+"\t"+sc);
                        }
                    }
                }
            } else
            { // scan downstream
                int score = 0;
                int start_col = scanned_column_idx;
                int end_col = scanned_column_idx;
                int ali_len = getNumAlignmentColumns();


                while (end_col<ali_len)
                {
                    query_phase = -1;
                    int ref_phase = -1;
                    if (aligned_genomic_positions[query_idx][end_col]==ALIGNMENT_GAP)
                    {
                        if (aligned_genomic_positions[ref_idx][end_col]==ALIGNMENT_GAP)
                        {
                            // both have gaps here
                            // nothing to do
                        } else
                        {
                            // query has a gap but not ref
                            score = 0;
                        }
                    } else // query has no gap here
                    {
                        if (aligned_genomic_positions[ref_idx][end_col]==ALIGNMENT_GAP)
                        {
                            // ref has a gap but query has none
                            score = 0;
                        } else
                        {
                            // check if they have introns here
                            query_phase = getIntronPhase(query_gene, end_col);
                            ref_phase = getIntronPhase(ref_gene, end_col);

                            boolean query_intron = query_phase == 1 || query_phase==2;
                            boolean ref_intron = ref_phase == 1 || ref_phase==2;

                            // get codons, get subst score
                            if (!query_intron && !ref_intron)
                            {
                                Codon query_codon = aligned_codons[query_idx][end_col];
                                Codon ref_codon = aligned_codons[ref_idx][end_col];
                                if (query_codon != null && ref_codon != null) // exclude codons with ambiguous nucleotides
                                {
                                    score += subst_scoring.getScore(query_codon, ref_codon);
                                }
                            } else
                            {
                                score = 0;
                            }
                        }
                    }
                    end_col++;

                    if (score>best_score)
                    {
                        best_score = score;
                        best_start_col =start_col;
                        best_end_col = end_col;
                        if (best_score>max_score_threshold)
                            break;
                    }
                    if (score <= 0 || ref_phase==3 || query_phase==3)
                    {
                        start_col = end_col;
                        score = 0;
                        //if (best_score>max_score_threshold)
                        //    break;
                    } 
                } // while end<alignment length

                if (DEBUG_TRACK_BRACKETING_ALIGNMENT)
                {
                    if (DEBUG_TRACK_BRACKETING_BOUNDARIES)
                    {
                        System.out.println("#*SSH.IC.cC col/dn "+scanned_column_idx+"\t"+getOrganism(query_idx)+"/"+getOrganism(ref_idx)+"\tsc "+best_score
                            +((best_score>max_score_threshold)?"*":"")
                            +"\t"+best_start_col+".."+best_end_col
                            +"\t// "+max_score_threshold);
                    }
                    if (best_score>0)
                    {
                        for (int c=scanned_column_idx; c<best_start_col; c++)
                        {
                            int sc =  (aligned_codons[query_idx][c]!=null && aligned_codons[ref_idx][c]!=null)
                                    ?subst_scoring.getScore(aligned_codons[query_idx][c], aligned_codons[ref_idx][c]):0;
                            System.out.println("#*SSH.IC.cC "+c+".\t"+describeAlignedCodonPosition(query_idx, c).toLowerCase()+"\t"+describeAlignedCodonPosition(ref_idx, c).toLowerCase()+"\t"+sc);
                        }
                        for (int c=best_start_col; c<Math.min(best_start_col+12,best_end_col); c++)
                        {
                            int sc =  (aligned_codons[query_idx][c]!=null && aligned_codons[ref_idx][c]!=null)
                                    ?subst_scoring.getScore(aligned_codons[query_idx][c], aligned_codons[ref_idx][c]):0;
                            System.out.println("#*SSH.IC.cC "+c+"|\t"+describeAlignedCodonPosition(query_idx, c)+"\t"+describeAlignedCodonPosition(ref_idx, c)+"\t"+sc);
                        }
                    }
                }
            } // if downstream


            if (scan_upstream)
            {
                this.left_bracket_score =  (best_score<max_score_threshold?0:best_score); 
                this.left_bracket_column = best_end_col;
//                if (best_end_col<IntronContext.this.first_column_idx) // best_end_col+1 <= first_column_idx
//                    first_column_idx = best_end_col+1;
            } else
            {
                this.right_bracket_score = (best_score<max_score_threshold?0:best_score);
                this.right_bracket_column = best_start_col;
//                if (best_start_col>OldIntronContext.this.last_column_idx)
//                    last_column_idx = best_start_col;
            }                
        }

        private int[] getCodingSegments(int seq_idx)
        {
            assert(hasValidBrackets());
            
            GenePred gene = homolog_genes[seq_idx];
            //StringBuilder sb = new StringBuilder();

            List<Integer> coding_segments = new ArrayList<>();

            int num_coding_nuc=0;
            //int lcol = left_bracket_column;
            //int rcol = right_bracket_column;

            int col= left_bracket_column-1;
            if (getIntronPhase(gene, col)==3)
            {
                //sb.append(num_coding_nuc).append(BRACKET_INTRON_SIGN);
                coding_segments.add(num_coding_nuc); // which is 0
            }
            col++;

            while (col<right_bracket_column)
            {
                int pos=aligned_genomic_positions[seq_idx][col];
                int intron_phase = getIntronPhase(gene, col);
                if (intron_phase == -1)
                {
                    num_coding_nuc += (pos==ALIGNMENT_GAP?0:3);
                } else
                {
                    num_coding_nuc += intron_phase;
                    //sb.append(num_coding_nuc).append(BRACKET_INTRON_SIGN);
                    coding_segments.add(num_coding_nuc);
                    num_coding_nuc = 3-intron_phase;
                }
                col++;
            }
            //sb.append(num_coding_nuc);
            coding_segments.add(num_coding_nuc);
            if (DEBUG_TRACK_BRACKETING_ALIGNMENT)
            {
                String org = getOrganism(seq_idx);
                for (int c=left_bracket_column; c<right_bracket_column; c++)
                {
                    System.out.print("\n#*SSH.IC.BH.gCS "+org+"\t"+c+"\t"+describeAlignedCodonPosition(seq_idx, c));
                }
                System.out.println();
            }
            //sb.append("(").append(first_column_idx).append("..").append(last_column_idx).append(")");

            int[] coding_segments_array = new int[coding_segments.size()];
            for (int ci=0; ci<coding_segments.size(); ci++)
                coding_segments_array[ci] = coding_segments.get(ci);

//            System.out.println("#*SSH.ICS.gCS seq "+seq_idx+"\t"+lcol+".."+rcol+"\t"+Arrays.toString(coding_segments_array));

            return coding_segments_array;
        }  

        BracketedHomology fuseBrackets(BracketedHomology other)
        {
            BracketedHomology fused = new BracketedHomology();
            if (this.left_bracket_score == Integer.MIN_VALUE) // unlikely
            {
                fused.left_bracket_column = other.left_bracket_column;
                fused.left_bracket_score = other.left_bracket_score;
            } else
            {
                fused.left_bracket_column = this.left_bracket_column;
                fused.left_bracket_score = this.left_bracket_score;
            }
            if (other.right_bracket_score == Integer.MIN_VALUE) // unlikely
            {
                fused.right_bracket_column = this.right_bracket_column;
                fused.right_bracket_score = this.right_bracket_score;
            } else
            {
                fused.right_bracket_column = other.right_bracket_column;
                fused.right_bracket_score = other.right_bracket_score;
            }
            return fused;
        }
        
        
        int getRegionStart(int seq_idx)
        {
            int pos_up = left_bracket_column-1;
            int seq_pos = aligned_genomic_positions[seq_idx][pos_up];
            if (seq_pos==ALIGNMENT_GAP)
                return seq_pos;
            else
                return seq_pos + (homolog_genes[seq_idx].isReverseStrand()?-3:+3);
        }

        int getRegionEnd(int seq_idx)
        {
            int pos_dn = right_bracket_column;

            int region_end = aligned_genomic_positions[seq_idx][pos_dn];
            return region_end;
        }
        
        int getSegmentCount(int seq_idx)
        { 
            if (segments[seq_idx]==null)
            {
                segments[seq_idx] = getCodingSegments(seq_idx);
            }
            return segments[seq_idx].length;
        }

        int getSegmentLength(int seq_idx, int segment_idx)
        {
            if (segments[seq_idx]==null)
            {
                segments[seq_idx] = getCodingSegments(seq_idx);
            }
            return segments[seq_idx][segment_idx];
        }
        
        
        void addSites(SpliceSiteHistory site_history, int qidx, int ridx)
        {
            DNAStrand query_transcript = transcript_strands[qidx];
            DNAStrand ref_transcript = transcript_strands[ridx];
            
            int query_region_start = getRegionStart(qidx);
            int query_region_end = getRegionEnd(qidx);
            int query_region_length = query_transcript.sub(query_region_end, query_region_start);
            int ref_region_start = getRegionStart(ridx);
            int ref_region_end = getRegionEnd(ridx);
            int ref_region_length = ref_transcript.sub(ref_region_end, ref_region_start);
                            
            int query_donor_offset = getSegmentLength(qidx, 0);
            int ref_donor_offset = getSegmentLength(ridx, 0);
            int donor_shift = ref_donor_offset - query_donor_offset;
            int query_acceptor_offset = getSegmentLength(qidx, getSegmentCount(qidx)-1);
            int ref_acceptor_offset = getSegmentLength(ridx, getSegmentCount(ridx)-1);
            int acceptor_shift = ref_acceptor_offset - query_acceptor_offset;
            
            if (getSegmentCount(ridx)==1)
            {
                site_history.addNoIntron(ref_transcript, true);
            } else 
            {
                int intron_phase = ref_donor_offset % 3;
                if (intron_phase==0) intron_phase = 3;
                assert(intron_phase>0 && intron_phase<=3);
                site_history.addDonorSite(ref_transcript, ref_transcript.add(ref_region_start, ref_donor_offset), intron_phase, true);
                site_history.addAcceptorSite(ref_transcript, ref_transcript.add(ref_region_end, -ref_acceptor_offset-1), intron_phase, true);
            }
            if (getSegmentCount(qidx)==1)
            {
                site_history.addNoIntron(query_transcript, true);
            } else
            {
                int intron_phase = query_donor_offset % 3;
                if (intron_phase==0) intron_phase = 3;
                assert(intron_phase>0 && intron_phase<=3);
                site_history.addDonorSite(query_transcript, query_transcript.add(query_region_start, query_donor_offset), intron_phase, true);
                site_history.addAcceptorSite(query_transcript, query_transcript.add(query_region_end, -query_acceptor_offset-1), intron_phase, true);
            }
            
            if (getSegmentCount(qidx)==2 && getSegmentCount(ridx)==1)
            {
                // intron loss
                int query_coding_length = query_donor_offset + query_acceptor_offset;
                int query_insertion_length = query_region_length-query_coding_length;

                
                int ref_insertion_length = ref_region_length-query_coding_length;
                if (query_region_length%3==0) // query_coding_length %3 == 0 for sure
                {
                    CodonCounter.Scoring qry_codon_scoring = agg.getCodonScoring(getOrganism(qidx));
                    int codon_sc = 0;
                    int donor_codon_pos = query_transcript.add(query_region_start,3*(query_donor_offset/3)); // inclusive
                    int acceptor_codon_pos = query_transcript.add(query_region_end,-3*(query_acceptor_offset/3)); // exclusive
                    int pos=donor_codon_pos;
//                    //  
//                    //     EEE
//                    //     ^ 
//                    //     |
//                    //     qend
//                    // offset
//                    //      0   0
//                    //      1   0
//                    //      2   0
//                    //      3   3
                    while (pos!=acceptor_codon_pos)
                    {
                        Codon query_codon = transcript_strands[qidx].getCodonAt(pos);
                        if (query_codon != null)
                        {
                            int s = qry_codon_scoring.getMiddleCodonScore(query_codon);
                            if (s==Integer.MIN_VALUE)
                            {
                                codon_sc=s; break;
                            }
                            codon_sc += s;
                        }
                        pos += query_transcript.getDownstreamStep()*3;
                    }
                    
                    if (codon_sc != Integer.MIN_VALUE && ref_insertion_length > query_insertion_length/3)
                    {
                        // exonify this intron?
                        site_history.addNoIntron(query_transcript, false);
                    }
                }
            }
            
            if (getSegmentCount(ridx)== 1 || (donor_shift != 0 && query_donor_offset <ref_region_length))
            {
                SpliceSiteComposition.Scoring ss5_scoring = agg.getDonorSiteScoring(Aggregator.GENERIC_SPLICE_SITE);
                int intron_phase = query_donor_offset % 3;
                if (intron_phase==0) intron_phase = 3;
                assert(intron_phase>0 && intron_phase<=3);
                
                if (ref_region_length>=ss5_scoring.getSiteLength()+query_donor_offset)
                {
                    int projected_pos = ref_transcript.add(ref_region_start, query_donor_offset);
                    int sc5 = ss5_scoring.getScore(ref_transcript.getDNA(), projected_pos, !ref_transcript.isReverseStrand());
                    if (sc5 != Integer.MIN_VALUE)
                    {
                        site_history.addDonorSite(ref_transcript, projected_pos, intron_phase, false);
                    } else
                    {
//                        ss5_scoring = agg.getDonorSiteScoring(getOrganism(ridx));
//                        if (ref_region_length>=ss5_scoring.getSiteLength()+query_donor_offset)
//                        {
//                            sc5 = ss5_scoring.getScore(ref_transcript.getDNA(), projected_pos, !ref_transcript.isReverseStrand());
//                            if (sc5 != Integer.MIN_VALUE)
//                            {
//                                site_history.addDonorSite(ref_transcript, projected_pos, intron_phase, false);
//                            }
//                        }
                    }
                }
            }         
            
            if (getSegmentCount(ridx)== 1 || (acceptor_shift != 0 && query_acceptor_offset<ref_region_length))
            {
                SpliceSiteComposition.Scoring ss3_scoring = agg.getAcceptorSiteScoring(Aggregator.GENERIC_SPLICE_SITE);
                if (ref_region_length >= ss3_scoring.getSiteLength()+query_acceptor_offset)
                {
                    int intron_phase = 3-(query_acceptor_offset % 3);
                    assert(intron_phase>0 && intron_phase<=3);
                    // offset phase
                    //  0     3
                    //  1     2
                    //  2     1
                    //  3     3
                    int projected_pos= ref_transcript.add(ref_region_end, -(1+query_acceptor_offset));
                    int sc3 = ss3_scoring.getScore(ref_transcript.getDNA(), projected_pos, !ref_transcript.isReverseStrand());
//                    System.out.println("#*SSH.BH.aS "+getOrganism(qidx)+"->"+getOrganism(ridx)+"\t"+projected_pos+":"+intron_phase+"\t"
//                            +DNASequence.toChar(ref_transcript.getNucleotideAt(ref_transcript.add(projected_pos, -1)))
//                            +DNASequence.toChar(ref_transcript.getNucleotideAt(projected_pos))
//                            +"/"+sc3+"\t// "+query_region_end+"->"+ref_region_end+"-"+query_acceptor_offset);
                    if (sc3 != Integer.MIN_VALUE)
                    {
                        site_history.addAcceptorSite(ref_transcript, projected_pos, intron_phase, false);
                    } else
                    { // try it again with the genome-specific scoring scheme
//                        ss3_scoring=agg.getAcceptorSiteScoring(getOrganism(ridx));
//                        if (ref_region_length >= ss3_scoring.getSiteLength()+query_acceptor_offset)
//                        {
//                            sc3 =ss3_scoring.getScore(ref_transcript.getDNA(), projected_pos, !ref_transcript.isReverseStrand());
//                            if (sc3 != Integer.MIN_VALUE)
//                            {
//                                site_history.addAcceptorSite(ref_transcript, projected_pos, intron_phase, false);
//                            }
//                        }
                    }
                }
            }
        }
        
        private void reportPairing(int qidx, int ridx)
        {
            GenePred query_gene = homolog_genes[qidx];
            GenePred ref_gene = homolog_genes[ridx];
            
            int query_region_start = getRegionStart(qidx);
            int query_region_end = getRegionEnd(qidx);
            int query_region_length = 
                    query_gene.isReverseStrand()?(query_region_start-query_region_end):(query_region_end-query_region_start);
            int ref_region_start = getRegionStart(ridx);
            int ref_region_end = getRegionEnd(ridx);
            int ref_region_length =
                    ref_gene.isReverseStrand()?(ref_region_start-ref_region_end):(ref_region_end-ref_region_start);
                            
            int query_donor_offset = getSegmentLength(qidx, 0);
            int ref_donor_offset = getSegmentLength(ridx, 0);
            int donor_shift = ref_donor_offset - query_donor_offset;
            int query_acceptor_offset = getSegmentLength(qidx, getSegmentCount(qidx)-1);
            int ref_acceptor_offset = getSegmentLength(ridx, getSegmentCount(ridx)-1);
            int acceptor_shift = ref_acceptor_offset - query_acceptor_offset;
                            
            StringBuilder classification = new StringBuilder();

            if (getSegmentCount(qidx)>2)
            {
                classification.append("\tQmul").append(getSegmentCount(qidx)-1);
            }
            
            if (getSegmentCount(ridx)>2)
            {
                classification.append("\tRmul").append(getSegmentCount(ridx)-1);
            }
            
            if (getSegmentCount(qidx)==2 && getSegmentCount(ridx)==1)
            {
                // intron loss
                int query_coding_length = query_donor_offset + query_acceptor_offset;
                int query_insertion_length = query_region_length-query_coding_length;

                assert (getSegmentLength(2,0)==ref_region_length);
                
                int ref_insertion_length = ref_region_length-query_coding_length;
                if (ref_insertion_length==0)
                {
                    classification.append("\tLexa(").append(query_insertion_length).append(")");
                } else if (ref_insertion_length>0)
                {
                    classification.append("\tLins(").append(ref_insertion_length).append(")");
                } else if (ref_insertion_length<0)
                {
                    classification.append("\tLdel(").append(-ref_insertion_length).append(")");
                }
                if (query_region_length%3==0) // query_coding_length %3 == 0 for sure
                {
                    CodonCounter.Scoring qry_codon_scoring = agg.getCodonScoring(getOrganism(qidx));
                    CodonTransitions.Scoring subst_scoring = agg.getCodonSubstScoring(getOrganism(ridx), getOrganism(qidx));
                    int codon_sc = 0;
                    int subst_sc = 0;
                    int donor_codon_pos = query_region_start+(query_donor_offset/3)*(query_gene.isReverseStrand()?-3:+3);
                    int acceptor_codon_pos = query_region_end-(query_acceptor_offset/3)*(query_gene.isReverseStrand()?-3:+3);
                    int pos=donor_codon_pos;
                    int ref_pos = ref_region_start + (query_donor_offset/3)*(ref_gene.isReverseStrand()?-3:3);
                    //  
                    //     EEE
                    //     ^ 
                    //     |
                    //     qend
                    // offset
                    //      0   0
                    //      1   0
                    //      2   0
                    //      3   3
                    while (pos!=acceptor_codon_pos)
                    {
                        Codon query_codon = transcript_strands[qidx].getCodonAt(pos);
                        if (query_codon != null)
                        {
                            int s = qry_codon_scoring.getMiddleCodonScore(query_codon);
                            if (s==Integer.MIN_VALUE)
                            {
                                codon_sc=s; break;
                            }
                            codon_sc += s;
                            Codon ref_codon = transcript_strands[ridx].getCodonAt(ref_pos);
                            if (ref_codon != null)
                            {
                                int r = subst_scoring.getScore(ref_codon, query_codon);
                                if (r==Integer.MIN_VALUE)
                                    subst_sc = r;
                                else 
                                    subst_sc = subst_sc+(subst_sc==Integer.MIN_VALUE?0:r);
                            }
                        }
                        pos += query_gene.isReverseStrand()?-3:+3;
                        ref_pos += ref_gene.isReverseStrand()?-3:+3;
                    }
                    classification.append(", Q3n(").append(query_region_length-query_coding_length);
                    classification.append(", cod ").append(codon_sc);
                    classification.append(", subst ").append(subst_sc).append(")");
                    if (codon_sc>0 && ref_insertion_length > query_insertion_length/3)
                    {
                        int query_donor_pos = query_region_start 
                        + (query_gene.isReverseStrand()?-1:+1)*query_donor_offset;
                        int query_acceptor_pos = query_region_end
                            -(query_gene.isReverseStrand()?-1:1)*(1+query_acceptor_offset);
                        // EDIT: exonify query intron
                    }
                } else
                {
                    classification.append(", qIlen(").append(query_region_length-query_coding_length).append(")");
                }
            }
            if (getSegmentCount(ridx)>1)
            {
                if (donor_shift==0)
                {
                    classification.append("\tM5(").append(query_donor_offset).append(")");
                    int projected_pos = ref_region_start 
                        + (ref_gene.isReverseStrand()?-1:+1)*ref_donor_offset;
                    // EDIT: donor match
                } else
                {
                    if (donor_shift<0)
                    {
                        classification.append("\tI5(").append(query_donor_offset).append(donor_shift).append(")");
                    } else 
                    {
                        classification.append("\tE5(").append(query_donor_offset).append("+").append(donor_shift).append(")");
                    }
                }
            }
            
            if (getSegmentCount(ridx)== 1 || (donor_shift != 0 && query_donor_offset <ref_region_length))
            {
                SpliceSiteComposition.Scoring ss5_scoring = agg.getDonorSiteScoring(Aggregator.GENERIC_SPLICE_SITE);
                if (ref_region_length>=ss5_scoring.getSiteLength()+query_donor_offset)
                {
                    int projected_pos = ref_region_start 
                            + (ref_gene.isReverseStrand()?-1:+1)*query_donor_offset;

                    int sc5 = ss5_scoring.getScore(ref_gene.getChrom(), projected_pos, !ref_gene.isReverseStrand());
                    classification.append(", sc5 ").append(sc5).append(" @ ").append(projected_pos);
                    if (sc5 != Integer.MIN_VALUE)
                    {
                        // EDIT: donor shift @projected_pos for ridx
                    }
                }
            }
            if (getSegmentCount(ridx)>1)
            {
                if (acceptor_shift==0)
                {
                    classification.append("\tM3(").append(query_acceptor_offset).append(")");
                    int projected_pos= ref_region_end
                            -(ref_gene.isReverseStrand()?-1:1)*(1+ref_acceptor_offset);
                    // EDIT: acceptor match
                } else
                {
                    if (acceptor_shift<0)
                    {
                        classification.append("\tI3(").append(query_acceptor_offset).append(acceptor_shift).append(")");
                    } else
                    {
                        classification.append("\tE3(").append(query_acceptor_offset).append("+").append(acceptor_shift).append(")");
                    }
                }
            }
            if (getSegmentCount(ridx)== 1 || (acceptor_shift != 0 && query_acceptor_offset<ref_region_length))
            {
                SpliceSiteComposition.Scoring ss3_scoring = agg.getAcceptorSiteScoring(Aggregator.GENERIC_SPLICE_SITE);
                if (ref_region_length >= ss3_scoring.getSiteLength()+query_acceptor_offset)
                {
                    int projected_pos= ref_region_end
                            -(ref_gene.isReverseStrand()?-1:1)*(1+query_acceptor_offset);
                    int sc3 = ss3_scoring.getScore(ref_gene.getChrom(), projected_pos, !ref_gene.isReverseStrand());
                    classification.append(", sc3 ").append(sc3).append(" @ ").append(projected_pos);
                    if (sc3 != Integer.MIN_VALUE)
                    {
                        // EDIT: acceptor shift @projected_pos for ridx
                    }
                }
            }
            System.out.println("#*SSH.BH.rP "+toString()
                    +"\t"+getOrganism(qidx)+"/"+getOrganism(ridx)+":\t"
                    +Arrays.toString(getCodingSegments(qidx))
                    +"/"
                    +Arrays.toString(getCodingSegments(ridx))
                    +classification
                    +"\t\t// "+query_region_start+".."+query_region_end+":"+query_region_length
                    +"/"+ref_region_start+".."+ref_region_end+":"+ref_region_length);
        }                
        

        @Override
        public String toString()
        {
            StringBuilder sb = new StringBuilder(); //getClass().getSimpleName());
            sb.append("[").append(left_bracket_column).append("..").append(right_bracket_column);
            sb.append(", sc ").append(left_bracket_score).append("/").append(right_bracket_score);
            sb.append("]");
            return sb.toString();
        }

    }
    
    public final class IntronContext
    {
        private final BracketedHomology[][] pairwise_homologies;
        private int first_column;
        private int last_column;
        
        private final Map<BracketedHomology,BracketedHomology> brackets;
        
        private IntronContext(int intron_column)
        {
            pairwise_homologies = new BracketedHomology[homolog_genes.length][homolog_genes.length];
            this.brackets = new HashMap<>();
            this.first_column = intron_column;
            this.last_column = intron_column;
            
            Map<String, Integer> projected_introns = projected_intron_phase.get(intron_column);
            for (int qidx=0; qidx<homolog_genes.length; qidx++) 
            {
                GenePred query_gene = homolog_genes[qidx];
                if (projected_introns.containsKey(query_gene.getName()))
                {
                    calculateBrackets(qidx, intron_column);
                }
            }
        }

        /**
         * Leftmost alignment column between the brackets (=minimum of {@link BracketedHomology#right_bracket_column} with positive scores) 
         * 
         * @return starting column for context (inclusive)
         */
        int getFirstColumn(){ return first_column;}

        /**
         * Rightmost alignment column between brackets (=maximum of {@link BracketedHomology#left_bracket_column} with positive scores)
         * 
         * @return end column for context (exclusive)
         */
        int getLastColumn(){ return last_column;}
        
        private boolean hasValidBrackets(int qidx, int ridx)
        {
            BracketedHomology H = pairwise_homologies[qidx][ridx];
            return H!=null && H.hasValidBrackets();
        }
        
        private void calculateBrackets(int query_idx, int column_idx)
        {
            for (int ref_idx=0; ref_idx<homolog_genes.length; ref_idx++)
            {
                if (query_idx != ref_idx)
                {
                    BracketedHomology H = new BracketedHomology(query_idx, ref_idx, column_idx);
                    if (brackets.containsKey(H))
                    {
                        pairwise_homologies[query_idx][ref_idx] = brackets.get(H);
                        //System.out.println("#*SSH.IC.cB "+getOrganism(query_idx)+"/"+getOrganism(ref_idx)+"\treusing "+pairwise_homologies[query_idx][ref_idx]);
                    } else
                    {
                        pairwise_homologies[query_idx][ref_idx] = H;
                        brackets.put(H,H);
                        //System.out.println("#*SSH.IC.cB "+getOrganism(query_idx)+"/"+getOrganism(ref_idx)+"\tnew "+pairwise_homologies[query_idx][ref_idx]);
                    }
                    if (H.hasValidBrackets())
                    {
                        if (H.left_bracket_column<first_column)
                            first_column = H.left_bracket_column;
                        if (H.right_bracket_column>last_column)
                            last_column = H.right_bracket_column;
                    }
                }
            }
        }
        
        
        void fuseContexts(IntronContext other)
        {
            int fused_first = this.first_column;
            int fused_last = this.last_column;
            
            brackets.clear();
            for (int i=0; i<homolog_genes.length; i++)
            {
                for (int j=0; j<homolog_genes.length; j++)
                {
                    if (i!=j)
                    {
                        BracketedHomology this_homo = pairwise_homologies[i][j];
                        BracketedHomology that_homo = other.pairwise_homologies[i][j];
                        BracketedHomology fused_homo;
                        if (this_homo==null)
                        {
                            fused_homo=that_homo;
                        } else
                        {
                            if (that_homo==null)
                            {
                                fused_homo=this_homo;
                            } else
                            {
                                fused_homo = this_homo.fuseBrackets(that_homo);
                            }
                        }
                        if (fused_homo != null)
                        {
//                            System.out.println("#*SSH.IC.fC "+getFirstColumn()+".."+getLastColumn()+" + "+other.getFirstColumn()+".."+other.getLastColumn()
//                                    +"\t"+homolog_genes[i].getChrom().getOrganism()+"/"+homolog_genes[j].getChrom().getOrganism()
//                                    +"\tthis "+this_homo+"\tthat "+that_homo+"\tfused "+fused_homo);
                            
                            if (brackets.containsKey(fused_homo))
                            {
                                fused_homo = brackets.get(fused_homo);
                            } else
                            {
                                brackets.put(fused_homo, fused_homo);
                            }
                            pairwise_homologies[i][j]=fused_homo;
                            if (fused_homo.hasValidBrackets())
                            {
                                if (fused_homo.left_bracket_column<fused_first)
                                    fused_first = fused_homo.left_bracket_column;
                                if (fused_homo.right_bracket_column>fused_last)
                                    fused_last = fused_homo.right_bracket_column;
                            }
                        }
                    }
                }
            } // for i
            this.first_column = fused_first;
            this.last_column = fused_last;
        }
        
        void reportPairings()
        {
            System.out.println("\n\n#*SSH.IC.rP "+this.toString());
            for (int i=0; i<homolog_genes.length; i++)
                for (int j=0; j<homolog_genes.length; j++)
                {
                    if (i!=j && hasValidBrackets(i,j))
                    {
                        pairwise_homologies[i][j].reportPairing(i,j);
                    }
                }
        }
        
        int getIntronCount(int seq_idx)
        {
            int num_introns = 0;
            GenePred gene = homolog_genes[seq_idx];
            int col= first_column-1;
            if (getIntronPhase(gene, col)==3)
            {
                num_introns++;
            }
            col++;
            while (col<last_column)
            {
                int ph = getIntronPhase(gene, col);
                if (ph!=-1)
                    num_introns++;
                col++;
            }
            return num_introns;
        }
        

        public SpliceSiteHistory collectSites()
        {
            MultiHomology[] donor_homologies = new MultiHomology[transcript_strands.length];
            MultiHomology[] acceptor_homologies = new MultiHomology[transcript_strands.length];
            
            Map<MultiHomology, Boolean> all_donor_homo=new HashMap<>();
            Map<MultiHomology, MultiHomology> donor_acceptor = new HashMap<>();
            
            int max_intron_count = 0;
            for (int seq_idx=0; seq_idx<transcript_strands.length; seq_idx++)
            {
                int num_introns = getIntronCount(seq_idx);
                max_intron_count = Math.max(num_introns, max_intron_count);
                if (num_introns>0)
                {
                    if (donor_homologies[seq_idx]==null)
                    {
                        MultiHomology donor = new MultiHomology(transcript_strands[seq_idx]);
                        MultiHomology acceptor = new MultiHomology(transcript_strands[seq_idx]);
                        donor_homologies[seq_idx] = donor;
                        acceptor_homologies[seq_idx] = acceptor;
                        boolean consistent_homo = collectMappings(seq_idx, donor_homologies, acceptor_homologies);
                        all_donor_homo.put(donor, consistent_homo);
                        donor_acceptor.put(donor, acceptor);
                    }
                }
            }
            int num_consistent = 0;
            for (MultiHomology MH: all_donor_homo.keySet())
            {
//                //int i=0; while (donor_homologies[i]!=MH) i++;
//                System.out.println("#*SSH.cG "
//                        +"\t"+getFirstColumn()+".."+getLastColumn()
//                        +"\tintrons "+max_intron_count
//                        +"\tnorgs "+MH.size()
//                        +"\t"+(all_donor_homo.get(MH)?"good":"bad")
//                        +"\tdonor "+MH);//+"\tacceptor "+acceptor_homologies[i]
                if (all_donor_homo.get(MH))
                    num_consistent++;
            }
            
            SpliceSiteHistory site_history = null;
            if (num_consistent==1 && max_intron_count==1)
            {
                for (MultiHomology donor: all_donor_homo.keySet())
                    if (all_donor_homo.get(donor) && donor.size()>homolog_genes.length/2)
                    {
                        MultiHomology acceptor = donor_acceptor.get(donor);
                        site_history = new SpliceSiteHistory(donor, acceptor);
                        for (int i=0; i<homolog_genes.length; i++)
                        {
                            if (donor.contains(transcript_strands[i]))
                            {
                                for (int j=0; j<homolog_genes.length; j++)
                                {
                                    if (i!=j && donor.contains(transcript_strands[j]) && hasValidBrackets(i,j))
                                    {
                                        pairwise_homologies[i][j].addSites(site_history, i, j);
                                    }
                                }
                            }
                        }
//                        for (DNAStrand seq: donor.getSequences())
//                        {
//                            System.out.println("#*SSH.cG "+seq.getDNA().getOrganism()+seq.getStrand()+"\tsites "+site_history.getSiteDescription(seq));
//                        }
                    }
            }
            return site_history;
        }
        
        private boolean collectMappings(int query_idx, MultiHomology[] donor_homologies, MultiHomology[] acceptor_homologies)
        {
//            System.out.println("#*SSH.cM traverse "+getOrganism(query_idx)+"\t"+donor_homologies[query_idx]);
            boolean mappings_ok = true;
            DNAStrand query_strand = transcript_strands[query_idx];
            for (int ref_idx=0; ref_idx<transcript_strands.length; ref_idx++)
            {
                if (hasValidBrackets(query_idx, ref_idx))
                {
                    DNAStrand ref_strand = transcript_strands[ref_idx];
                    BracketedHomology B = pairwise_homologies[query_idx][ref_idx];
                    int query_start_pos = B.getRegionStart(query_idx);
                    int ref_start_pos = B.getRegionStart(ref_idx);
                    boolean ref_donor_ok = donor_homologies[query_idx].addMatch(query_strand, query_start_pos, ref_strand, ref_start_pos);
                    if (!ref_donor_ok)
                    {
//                        System.out.println("#*SSH.cM bad donor map "+B+"\t"+getOrganism(query_idx)+"/"+getOrganism(ref_idx)+"\tqstart "+query_start_pos+"\trstart "+ref_start_pos+"\twant "+donor_homologies[query_idx].projectCoordinate(query_strand, query_start_pos, ref_strand));
                    }
                    int query_end_pos = B.getRegionEnd(query_idx);
                    int ref_end_pos = B.getRegionEnd(ref_idx);
                    boolean ref_acceptor_ok = acceptor_homologies[query_idx].addMatch(query_strand, query_end_pos, ref_strand, ref_end_pos);
                    if (!ref_acceptor_ok)
                    {
//                        System.out.println("#*SSH.cM bad accceptor map "+B+"\t"+getOrganism(query_idx)+"/"+getOrganism(ref_idx)+"\tqend "+query_end_pos+"\trend "+ref_end_pos+"\twant "+acceptor_homologies[query_idx].projectCoordinate(query_strand, query_end_pos, ref_strand));
                    }
                    mappings_ok = mappings_ok && ref_donor_ok && ref_acceptor_ok;
                    if (donor_homologies[ref_idx]==null)
                    {
                        // depth-first traversal
                        donor_homologies[ref_idx] = donor_homologies[query_idx];
                        acceptor_homologies[ref_idx] = acceptor_homologies[query_idx];
                        boolean connections_ok = collectMappings(ref_idx, donor_homologies, acceptor_homologies);
                        mappings_ok = mappings_ok && connections_ok;
                    }
                }
            }
            return mappings_ok;
        }
        
        private List<List<Integer>> getGrouping()
        {
            int[] group_ident=new int[homolog_genes.length];
            for (int i=0; i<homolog_genes.length; i++)
                group_ident[i]=i;
            for (int i=0; i<homolog_genes.length; i++)
            {
                int id = group_ident[i];
                while (group_ident[id]!=id)
                {
                    id = group_ident[id]=group_ident[group_ident[id]]; // path halving
                }

                for (int j=0; j<homolog_genes.length; j++)
                    if (i!=j && hasValidBrackets(i,j))
                    {
                        group_ident[j] = id;
                    }
            }
            Counter<Integer> group_size = new Counter<>();
            for (int i=0;i<homolog_genes.length; i++)
            { // compression
                int id = group_ident[i];
                while (id!=group_ident[id]) id=group_ident[id];
                group_ident[i]=id;
                group_size.increment(id);
            }
            Integer[] group_identifiers = group_size.toArray(new Integer[0]); // decreasing order of size
            List<List<Integer>> groupings = new ArrayList<>();
            for (int id: group_identifiers) groupings.add(new ArrayList<>());
            
            for (int i=0; i<homolog_genes.length; i++)
            {
                int gidx=0; while (group_identifiers[gidx]!=group_ident[i]) gidx++;
                groupings.get(gidx).add(i);
            }
            return groupings;
        }        
        

        @Override
        public String toString()
        {
            StringBuilder sb = new StringBuilder("C"); //getClass().getSimpleName());
            sb.append("{").append(getFirstColumn()).append("..").append(getLastColumn());
            sb.append(" ").append(brackets.keySet());
//            sb.append(" {");
//            for (BracketedHomology H: brackets.keySet())
//            {
//                sb.append(" ").append(H);
//            }
//            sb.append("}");
            sb.append("}");
            return sb.toString();
        }
    }
    
    public IntronContext[] calculateIntronContexts()
    {
        List<IntronContext> contexts = new ArrayList<>();
        Integer[] columns_with_introns = getIntronColumns();
        for (int col_idx: columns_with_introns)
        {
            IntronContext current_context = new IntronContext(col_idx);
            
            int previous_context_idx = contexts.size();
            while (previous_context_idx>0)
            {
                previous_context_idx--;
                IntronContext previous_context = contexts.get(previous_context_idx);
                if (previous_context.getLastColumn()>=current_context.getFirstColumn())
                {
                    //System.out.println("#*SSH.cIC fuse "+previous_context+"\t"+current_context);
                    
                    previous_context.fuseContexts(current_context);
                    contexts.remove(previous_context_idx);
                    current_context = previous_context;
                } else
                    break;
            }
            //System.out.println("#*SSH.cIC "+contexts.size()+"\tcol "+col_idx+"\t"+current_context);
            contexts.add(current_context);
        }
        
        return contexts.toArray(new IntronContext[0]);
        
    }

    
    /**
     * String representation of an alignment cell, with intron info.
     * 
     * @param seq_idx sequence index
     * @param col_idx column index
     * @return string giving aligned query_codon's genomic position, query_codon with amino acid translation, and phased intron (if any)
     */
    private String describeAlignedCodonPosition(int seq_idx, int col_idx)
    {
        Codon C = aligned_codons[seq_idx][col_idx];
        GenePred gene = homolog_genes[seq_idx];

        Map<String, Integer> intron_phase_pos = projected_intron_phase.get(col_idx);

        String codon_desc = (C==null?"---/-":C.toString(TranslationTable.STANDARD_TABLE));
        if (intron_phase_pos!=null && intron_phase_pos.containsKey(gene.getName()))
        { // intron inserted according  to phase
            int phase = intron_phase_pos.get(homolog_genes[seq_idx].getName());
            codon_desc = codon_desc.substring(0,phase)+"^"+Integer.toString(phase)+
                    codon_desc.substring(phase);
        }
        int pos = aligned_genomic_positions[seq_idx][col_idx];
        
        // padding for pretty-print if gap, otherwise aligned query_codon's position
        String pos_info = 
                pos == ALIGNMENT_GAP
                ?"            ".substring(0,Integer.toString(gene.getCDSEnd()).length()+1)
                :(Integer.toString(pos)+":");
        return pos_info+codon_desc;
    }
    
    private String getIntronSequence(int seq_idx, int col_idx)
    {
        GenePred gene = homolog_genes[seq_idx];
        int pos = aligned_genomic_positions[seq_idx][col_idx];
        int phase = getIntronPhase(gene, col_idx);
        
        int wanted_site_context_length = 60;

        String intronseq = null;
        
        if (phase != NO_INTRON)
        {
            if (gene.isReverseStrand())
            {
                pos -= phase;
                // <<<<<  iiEEEEE
                //         ^
                //         |
                //        pos
                int eidx=0; while (eidx<gene.getExonCount() && gene.getExonStart(eidx)-1<pos) eidx++;
                if (eidx<gene.getExonCount() && eidx>0)
                {
                    intronseq = gene.getSequenceIntron(eidx-1).toLowerCase();
                } 
            } else
            {
                pos += phase;
                // EEEEEiii
                //      ^
                //      |
                //     pos
                int eidx=0; while (eidx<gene.getExonCount() && gene.getExonEnd(eidx)<pos) eidx++;
                if (eidx<gene.getExonCount()-1)
                {
                    intronseq = gene.getSequenceIntron(eidx).toLowerCase();
                }
            }
        }
        if (intronseq==null)
        {
            return "(no intron at "+pos+"in column "+col_idx+"; gene "+gene.shortDesc(true)+")";
        } 
        if (intronseq.length()>2*wanted_site_context_length+10)
        {
            StringBuilder sb = new StringBuilder();
            sb.append(intronseq.substring(0,wanted_site_context_length));
            sb.append("..").append(intronseq.length()-2*wanted_site_context_length).append("..");
            sb.append(intronseq.substring(intronseq.length()-wanted_site_context_length));
            intronseq = sb.toString();
        } 
        return intronseq;
    }

    
    
    private void neki(String[] args) throws Exception
    {
        if (args.length==0 || args[0].equals("-h")){
            System.err.println("Tests: call as "+getClass().getName()+" org1=g1.fa,org2=g2.fa annot.txt agg.data aa-alignment.fa");
            System.err.println("Checks splice site homologies.");
            System.exit(9);
        }
        
        int arg_idx=0;
        
        while (arg_idx<args.length && args[arg_idx].startsWith("-"))
        {
            String option = args[arg_idx++];
            String value = args[arg_idx++];
            
            if (OPTION_NUM_MATCHING.equals(option))
            {
                NUM_MATCHING_NUCLEOTIDES = Integer.parseInt(value);
            } else if (OPTION_TRACK_ALIGNMENT.equals(option))
            {
                DEBUG_TRACK_BRACKETING_ALIGNMENT = "true".equalsIgnoreCase(value) || "yes".equalsIgnoreCase(value);
            } else if (OPTION_TRACK_BOUNDARY.equals(option))
            {
                DEBUG_TRACK_BRACKETING_BOUNDARIES = "true".equalsIgnoreCase(value) || "yes".equalsIgnoreCase(value);
            } else
            {
                throw new IllegalArgumentException("Undefined option "+option+"; known options are "+OPTION_NUM_MATCHING+", "+OPTION_TRACK_ALIGNMENT+", "+OPTION_TRACK_BOUNDARY);
            }
        }
        
        String genome_fa_list = args[arg_idx++];
        String annotation_file = args[arg_idx++];
        String scoring_file = args[arg_idx++];
        String alignment_file = args[arg_idx++];
        
        AnnotatedGenomes annotations = new AnnotatedGenomes();
        List<String> wanted_organisms = annotations.readMultipleGenomes(genome_fa_list);
        annotations.readAnnotations(annotation_file);  

        setScoring(Aggregator.readData(new GeneralizedFileReader(scoring_file)));
        if (!computeAlignmentColumns(annotations, ProteinSequence.readMultiFasta(alignment_file)))
            System.exit(1999);
        
        printAlignment(System.out);

//        Integer[] columns_with_introns = projected_intron_phase.keySet().toArray(new Integer[0]);
//        Arrays.sort(columns_with_introns);
//        for (int col_idx: columns_with_introns)
//            reportIntronColumn(col_idx);
        
        IntronContext[] contexts = calculateIntronContexts();
        for (IntronContext C: contexts)
        {
            C.reportPairings();
            C.collectSites();
        }
    }
    
    public static void main(String[] de_nagyon) throws Exception
    {
        
        SpliceSiteHomology adj = new SpliceSiteHomology();
        adj.neki(de_nagyon);
        
    }
    
    
}
