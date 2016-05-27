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
 *
 * @author Mikl&oacute;s Cs&#369;r&ouml;s
 */
public class TranslationTable implements NucleotideEncoding 
{
    

    /**
     * This translation table is instantiated by a string of 
     * translated amino acids, where codons are enumerated in
     * lexicographic order. The lexicographic order itself 
     * depends on the chosen order of nucleotides: nuc1&lt;nuc2&lt;nuc3&lt;nuc4.
     * Stop codon is denoted by <code>*</code>.
     * 
     * @param aa is a 64-lettter String of amino acid 1-character codes (stop=*)
     *        in the lexicographic order of codons
     * @param nuc1 first nucleotide in the ordering
     * @param nuc2 second nucleotide in the ordering
     * @param nuc3 third nucleotide in the ordering 
     */
    public TranslationTable(String aa, char nuc1, char nuc2, char nuc3) 
    {
        init(aa, nuc1, nuc2, nuc3);    
    }
    
    /**
     * This translation table is instantiated by a string of 
     * translated amino acids, where codons are enumerated in
     * lexicographic order, assuming T C A G order as at NCBI. 
     * @param aa is a 64-lettter String of amino acid 1-character codes (stop=*),
     *        in the lexicographic order of codons.
     */
    public TranslationTable(String aa)
    {
        this(aa, CHAR_T, CHAR_C, CHAR_A);
    }
    
    private void init(String aa, char nuc1, char nuc2, char nuc3)
    {
        nuc_order=new int[4];
        nuc_order[0]=nuc_order[1]=nuc_order[2]=nuc_order[3]=3;
        nuc_order[DNASequence.toIndex(nuc1)]=0;
        nuc_order[DNASequence.toIndex(nuc2)]=1;
        nuc_order[DNASequence.toIndex(nuc3)]=2;
        
        table=aa.toCharArray();       
        for (int i=0; i<63; i++)
        {
            Codon C=Codon.getCodon(i);
            //System.out.println("#STT.i codon "+C.toString(this));
        }
    }
    
    private int[] nuc_order;
    private char[] table;
    
    /**
     * Whether a codon represents a STOP signal in translation.
     * 
     * @param codon codon
     * 
     * @return true if codon encodes for stop
     */
    public boolean isStopCodon(Codon codon) 
    {
        return table[getIndex(codon)]==CHAR_STOP;
    }
    
    /**
     * @param codon the codon that needs to be translated
     * @return a one-letter code of an amino acid if a sense codon, otherwise a String with CHAR_STOP  
     */
    public char toAA(Codon codon) 
    {
        return table[getIndex(codon)];
    }
    
    private int getIndex(Codon codon)
    {
        int i1=nuc_order[DNASequence.toIndex(codon.nucleotideAt(1))];
        int i2=nuc_order[DNASequence.toIndex(codon.nucleotideAt(2))];
        int i3=nuc_order[DNASequence.toIndex(codon.nucleotideAt(3))];
        int idx=i1*16+i2*4+i3;
        
        return idx;
    }
    
    public static final char CHAR_STOP='*';
    
    public static final TranslationTable STANDARD_TABLE
        =new TranslationTable(
            "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG");
}


