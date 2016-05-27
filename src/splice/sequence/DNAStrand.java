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
package splice.sequence;

import java.util.Arrays;
import java.util.Objects;

/**
 *  Class for DNA strand and cached codon translation.
 * 
 * @author Mikl&oacute;s Cs&#369;r&ouml;s
 */
public final class DNAStrand 
{
    /**
     * Character constant for reverse strand (<code>'-'</code>).
     */
    public static final char REVERSE = '-';
    /**
     * Character constant for forward strand (<code>'-'</code>).
     */
    public static final char FORWARD = '+';
    /**
        * Character constant for unknown strandedness (<code>'.'</code>).
     */
    public static final char UNKNOWN = '.';
    
    /**
     * Instantiation with DNA sequence and strand. 
     * 
     * @param seq DNA sequence for forward or reverse strand
     * @param strand which strand is used here
     */
    public DNAStrand(DNASequence seq, char strand)
    {
        this.seq = seq;
        this.strand = strand;
        this.cache = new CodonCache();//seq, strand);
    }
    
//    public DNAStrand(GenePred gene)
//    {
//        this(gene.getChrom(), gene.getStrand());
//    }
    
    private final DNASequence seq;
    private final char strand;
    private final CodonCache cache;
    
    public boolean isReverseStrand()
    {
        return strand == REVERSE;
    }
    
    public Codon getCodonAt(int target_pos)
    {
        int codon_idx = cache.getCodonIdx(target_pos);
        return (codon_idx==AMBIGUOUS_CODON?null:Codon.getCodon(codon_idx));
    }
    
    private static final int AMBIGUOUS_CODON = Codon.NUM_CODONS;
    public int getCodonIdxAt(int pos)
    {
        return cache.getCodonIdx(pos);
    }
    
    /**
     * DNA sequence used at instantiation (not inverted for reverse strand). 
     * 
     * @return DNA sequence 
     */
    public DNASequence getDNA()
    {
        return seq;
    }

    /**
     * Nucleotide (complemented if necessary) at a given position.
     * 
     * @param pos non-negative 0-based coordinate from start of DNA sequence
     * @return complemented nucleotide if reverse strand
     */
    public byte getNucleotideAt(int pos)
    {
        byte b = seq.getNucleotideAt(pos);
        return (isReverseStrand()
                ?DNASequence.complement(b)
                :b);
    }

    public char getStrand(){ return strand;}

    /**
     * Downstream offset determined by strandedness (+1 or -1). 
     * 
     * @return -1 for reverse strand, +1 otherwise
     */
    public int getDownstreamStep()
    {
        return isReverseStrand()?-1:1;
    }

    
    
    /**
     * Equality test
     * 
     * @param other the other DNA strand
     * @return true if same underlying DNA sequence and strand
     */
    @Override
    public boolean equals(Object other)
    {
        if (other instanceof DNAStrand)
        {
            DNAStrand o = (DNAStrand)other;
            return this.getDNA().equals(o.getDNA())
                    && this.getStrand() == o.getStrand();
        } else
            return super.equals(other);
    }

    @Override
    public int hashCode() 
    {
        int hash = 5;
        hash = 17 * hash + Objects.hashCode(this.seq);
        hash = 17 * hash + this.strand;
        return hash;
    }
    
    @Override
    public String toString()
    {
        return seq.getIdent()+"/"+strand;
    }
           
    public SlidingWindow getWindow(int window_size)
    {
        return new SlidingWindow(window_size);
    }
    
    
    public int add(int pos, int offset)
    {
        return isReverseStrand()?pos-offset:pos+offset;
    }
    
    public int sub(int end_pos, int start_pos)
    {
        return isReverseStrand()?(start_pos-end_pos):(end_pos-start_pos);
    }
    
//    public Arithmetics getCoordinateArithmetics()
//    {
//        if (isReverseStrand())
//        {
//            return new ReverseStrandArithmetics();
//        } else
//        {
//            return new ForwardStrandArithmetics();
//        }
//    }
//    
//    public static Arithmetics getCoordinateArithmetics(char strand)
//    {
//        if (strand==REVERSE)
//            return new ReverseStrandArithmetics();
//        else
//            return new ForwardStrandArithmetics();
//    }
//    
//    public interface Arithmetics
//    {
//        public abstract int add(int position,int offset);
//        public abstract int sub(int end_position, int start_position);
//    }
//    
//    private static final class ForwardStrandArithmetics implements Arithmetics
//    {
//        @Override
//        public int add(int position, int offset)
//        {
//            return position+offset;
//        }
//        
//        @Override
//        public int sub(int end_position, int start_position)
//        {
//            return end_position-start_position;
//        }
//    }
//    
//    private static final class ReverseStrandArithmetics implements Arithmetics
//    {
//        @Override
//        public int add(int position, int offset)
//        {
//            return position-offset;
//        }
//        
//        @Override
//        public int sub(int end_position, int start_position)
//        {
//            return start_position-end_position;
//        }
//    }
    
    /**
     * Sliding window respecting strand orientation.
     */
    public final class SlidingWindow extends DNAWindow
    {
        /**
         * 
         * @param width window size
         */
        private SlidingWindow(int width)
        {
            super(width);
        }
        
        private int next_position;
        
        /**
         * Sets the window position and initializes window content.
         * 
         * @param pos position of the first nucleotide in the window (right end on reverse strand, left end on forward strand)
         */
        public void setPosition(int pos)
        {
            next_position = pos;
            super.reset();
            for (int j=0; j<getWindowLength(); j++)
            {
                byte nuc = getNucleotideAt(next_position);
                nextNucleotide(nuc);
                next_position += getDownstreamStep();
            }
        }
        
        /** 
         * Slides the window by one position downstream.
         */
        public void stepForward()
        {
            byte nuc = getNucleotideAt(next_position);
            nextNucleotide(nuc);
            next_position += getDownstreamStep();
        }
        
        /**
         * Coordinate of the window's 5' end. 
         * 
         * @return position forthe first nucleotide covered by this window
         */
        public int getPosition()
        {
            return next_position - getDownstreamStep()*getWindowLength();
        }        
    }
    
    /**
     *
     * Cached access to translated DNA sequence on forward or reverse strand.
     * 
     * @author Mikl&oacute;s Cs&#369;r&ouml;s
     */
    private class CodonCache 
    {
        private static final int BLOCK_SIZE= 256; 

        private final byte[][][] codon_indices;

//        private final DNASequence seq;
//        private final char strand;

    //    private CodonCache(int sequence_length) 
    //    {
    //        int num_blocks = ((sequence_length+2)/3 + BLOCK_SIZE - 1) / BLOCK_SIZE;
    //        codon_indices = new byte[3][num_blocks][];
    //        seq = null;
    //        strand = '.';
    //    }

//        /**
//         * Instantiation with a sequence and strand specification.
//         * 
//         * @param seq 
//         * @param strand 
//         */
//        CodonCache(DNASequence seq, char strand)
        CodonCache()
        {
            int sequence_length = seq.getLength();
            int num_blocks = ((sequence_length+2)/3 + BLOCK_SIZE - 1) / BLOCK_SIZE;
            codon_indices = new byte[3][num_blocks][];
//            this.seq = seq;
//            this.strand = strand;
        }

        private int get(int pos) 
        {
            int phase = pos % 3;
            int phase_pos = pos / 3;
            byte[][] phase_indices = codon_indices[phase];
            int block_idx = phase_pos / BLOCK_SIZE;

            if (phase_indices[block_idx] == null) 
            {
                return -1;
            } else 
            {
                int offset = phase_pos % BLOCK_SIZE;
                return (int) (phase_indices[block_idx][offset]);
            }
        }

        private void put(int pos, byte codon_idx) 
        {
            int phase = pos % 3;
            int phase_pos = pos / 3;
            byte[][] phase_indices = codon_indices[phase];
            int block_idx = phase_pos / BLOCK_SIZE;
            if (phase_indices[block_idx] == null) {
                phase_indices[block_idx] = new byte[BLOCK_SIZE];
                Arrays.fill(phase_indices[block_idx], (byte) -1);
            }
            int offset = phase_pos % BLOCK_SIZE;
            phase_indices[block_idx][offset] = codon_idx;
        }
        

        /**
         * Codon index for a given position along the DNA sequence {@link #seq}.
         * 
         * @param pos position for first nucleotide in the codon
         * @return 
         */
        public int getCodonIdx(int pos)
        {
            int codon_idx = get(pos);
            if (codon_idx<0)
            {
                if (strand == REVERSE)
                {
                    byte b1 = DNASequence.complement(seq.getNucleotideAt(pos));
                    byte b2 = DNASequence.complement(seq.getNucleotideAt(pos-1));
                    byte b3 = DNASequence.complement(seq.getNucleotideAt(pos-2));
                    
                    Codon cod = Codon.getCodon(b1, b2, b3);
                    codon_idx = (cod==null?AMBIGUOUS_CODON:cod.getIndex());
//                    
//                    if (DNASequence.isAmbiguous(b1)
//                       || DNASequence.isAmbiguous(b2)
//                       || DNASequence.isAmbiguous(b3))
//                    {
//                        codon_idx = AMBIGUOUS_CODON;
//                    } else
//                    {
//                        Codon codon = new Codon(
//                            DNASequence.toChar(b1),
//                            DNASequence.toChar(b2),
//                            DNASequence.toChar(b3));
//                        codon_idx = codon.getIndex();
//                        //System.out.println("#*CS.RS.gS "+pos+"\t"+codon.toString(SimpleTranslationTable.STANDARD_TABLE)+"\t"+getCodonScore(codon)+"\t// "+getCodonScore(codon_idx));
//                    } 
                } else
                {
                    byte b1 = seq.getNucleotideAt(pos);
                    byte b2 = seq.getNucleotideAt(pos+1);
                    byte b3 = seq.getNucleotideAt(pos+2);
                    Codon cod = Codon.getCodon(b1, b2, b3);
                    codon_idx = (cod==null?AMBIGUOUS_CODON:cod.getIndex());
//                    if (DNASequence.isAmbiguous(b1)
//                       || DNASequence.isAmbiguous(b2)
//                       || DNASequence.isAmbiguous(b3))
//                    {
//                        codon_idx = AMBIGUOUS_CODON;
//                    } else
//                    {
//                        Codon codon = new Codon(
//                            DNASequence.toChar(b1),
//                            DNASequence.toChar(b2),
//                            DNASequence.toChar(b3));
//                        codon_idx = codon.getIndex();
//                        //System.out.println("#*CS.FS.gS "+pos+"\t"+codon.toString(SimpleTranslationTable.STANDARD_TABLE)+"\t"+getCodonScore(codon)+"\t// "+getCodonScore(codon_idx));
//                    } 
                }
                put(pos, (byte)codon_idx);
                
            }
            return codon_idx;
        }
    }
    
    
    
}