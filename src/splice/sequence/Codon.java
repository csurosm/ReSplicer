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

/**
* Class for representing codons. 
 *
 * @author Mikl&oacute;s Cs&#369;r&ouml;s
 */
public class Codon implements NucleotideEncoding, Comparable<Codon>
{

    public static final int NUM_CODONS = 64;
    
    
    /** Creates a new instance of Codon 
     * @param nuc1 nucleotide in first position
     * @param nuc2 nucleotide in second position
     * @param nuc3 nucleotide in third position
     */
    public Codon(char nuc1, char nuc2, char nuc3) 
    {
//        int i1=DNASequence.toIndex(nuc1);
//        int i2=DNASequence.toIndex(nuc2);
//        int i3=DNASequence.toIndex(nuc3);
        value=codonIdx(nuc1, nuc2, nuc3);
        
        nucleotides=new char[3];
        nucleotides[0]=nuc1;
        nucleotides[1]=nuc2;
        nucleotides[2]=nuc3;
    }
    
    private static final int codonIdx(char nuc1, char nuc2, char nuc3)
    {
        int i1=DNASequence.toIndex(nuc1);
        int i2=DNASequence.toIndex(nuc2);
        int i3=DNASequence.toIndex(nuc3);
        
        int idx=i1*16+i2*4+i3;
        return idx;
    }
    
    /**
     * Instantiation with codon index
     * 
     * @param idx 0..63 (only lower 6 bits are used for other values)
     */
    public Codon(int idx)
    {
        this(DNASequence.toChar((idx & 48)>>4), DNASequence.toChar((idx & 12)>>2), DNASequence.toChar(idx & 3));
    }
    
    public static final Codon[] ALL_CODONS;
    static 
    {
        ALL_CODONS = new Codon[NUM_CODONS];
        for (int cidx=0; cidx<NUM_CODONS; cidx++)
            ALL_CODONS[cidx] = new Codon(cidx);
    }
    
    /**
     * @param cidx of a codon in lexicographic ordering
     * @return a codon for which getIndex() gives cidx
     * @throws ArrayIndexOutOfBoundsException if codon index is out of range (0..{@link #NUM_CODONS}-1)
     */
    public static Codon getCodon(int cidx)
    {
        return ALL_CODONS[cidx];
    }
            
    public static Codon getCodon(char nuc1, char nuc2, char nuc3)
    {
        if (DNASequence.isAmbiguous(nuc1)||DNASequence.isAmbiguous(nuc2)||DNASequence.isAmbiguous(nuc3))
            return null;
        else
            return getCodon(codonIdx(nuc1, nuc2, nuc3));
    }
    
    public static Codon getCodon(byte nuc1, byte nuc2, byte nuc3)
    {
        if (DNASequence.isAmbiguous(nuc1)||DNASequence.isAmbiguous(nuc2)||DNASequence.isAmbiguous(nuc3))
            return null;
        else
            return getCodon(codonIdx(DNASequence.toChar(nuc1), DNASequence.toChar(nuc2), DNASequence.toChar(nuc3)));
    }
    
    protected final int value;
    private final char[] nucleotides;
    /**
     * @return a value 0..63 by the lexicographical ordering of codons
     */
    public int getIndex(){return value;}
    
    /**
     * @param i 1,2,3 for position within the codon
     * @return A,C,G,T
     */
    public char nucleotideAt(int i)
    {
        return nucleotides[i-1];
    }
    
    /**
     * @return a three-letter String for this codon
     */
    @Override
    public String toString()
    {
        return new String(nucleotides);
    }
    

    /**
     * Description of this codon and its translation.
     * 
     * @param T codon-to-amino-acid translation table
     * @return a String of "nnn/a" where nnn is the codon String and a is its amino acid translation
     */
    public String toString(TranslationTable T)
    {
        StringBuilder sb=new StringBuilder();
        sb.append(nucleotides[0]);
        sb.append(nucleotides[1]);
        sb.append(nucleotides[2]);
        sb.append('/');
        sb.append(T.toAA(this));
        return sb.toString();
    }
    
//    /**
//     * @param idx of a codon in lexicographic ordering
//     * @return a codon for which getIndex() gives idx
//     */
//    public static Codon getCodon(int idx)
//    {
//        return new Codon(idx);
//    }
    

    @Override
    public boolean equals(Object o)
    {
        if (o instanceof Codon)
            return equals((Codon)o);
        else
            return false;
    }
    
    public boolean equals(Codon c)
    {
        if (c==null)
            return false;
        else
            return (nucleotides[0] == c.nucleotides[0]
                && nucleotides[1] == c.nucleotides[1]
                && nucleotides[2] == c.nucleotides[2]);
    }

    @Override
    public int hashCode()
    {
        return value;
    }
    
    /**
     * Number of different nucleotides between two codons. 
     * 
     * @param c the other codon
     * @return Hamming distance
     */
    public int distance(Codon c)
    {
        int d = 0;
        if (nucleotides[0] != c.nucleotides[0]) d++;
        if (nucleotides[1] != c.nucleotides[1]) d++;
        if (nucleotides[2] != c.nucleotides[2]) d++;
        return d;
    }

    @Override
    public int compareTo(Codon o) 
    {
        return Integer.compare(value, o.value);
    }
}