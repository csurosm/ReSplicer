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
 * One-character, short and long names for amino acids.
 * 
 * @author Mikl&oacute;s Cs&#369;r&ouml;s
 */
public final class AminoAcid implements Comparable<AminoAcid>
{

    /**
     * Nobody should start introducing new amino acids.
     */
    private AminoAcid(int index, char oneCharacterName, String shortName, String longName){
        this.index = index;
        this.code = oneCharacterName;
        this.shortName = shortName;
        this.longName = longName;
    }

    private final int index;

    private final char code;
    private final String shortName;
    private final String longName;


    public static final AminoAcid A    = new AminoAcid( 0, 'A', "Ala", "Alanine");
    public static final AminoAcid C    = new AminoAcid( 1, 'C', "Cys", "Cysteine");
    public static final AminoAcid D    = new AminoAcid( 2, 'D', "Asp", "Aspartic acid");
    public static final AminoAcid E    = new AminoAcid( 3, 'E', "Glu", "Glutamic acid");
    public static final AminoAcid F    = new AminoAcid( 4, 'F', "Phe", "Phenylalanine");
    public static final AminoAcid G    = new AminoAcid( 5, 'G', "Gly", "Glycine");
    public static final AminoAcid H    = new AminoAcid( 6, 'H', "His", "Histidine");
    public static final AminoAcid I    = new AminoAcid( 7, 'I', "Ile", "Isoleucine");
    public static final AminoAcid K    = new AminoAcid( 8, 'K', "Lys", "Lysine");
    public static final AminoAcid L    = new AminoAcid( 9, 'L', "Leu", "Leucine");
    public static final AminoAcid M    = new AminoAcid(10, 'M', "Met", "Methionine");
    public static final AminoAcid N    = new AminoAcid(11, 'N', "Asn", "Asparagine");
    public static final AminoAcid P    = new AminoAcid(12, 'P', "Pro", "Proline");
    public static final AminoAcid Q    = new AminoAcid(13, 'Q', "Gln", "Glutamine");
    public static final AminoAcid R    = new AminoAcid(14, 'R', "Arg", "Arginine");
    public static final AminoAcid S    = new AminoAcid(15, 'S', "Ser", "Serine");
    public static final AminoAcid T    = new AminoAcid(16, 'T', "Thr", "Threonine");
    public static final AminoAcid V    = new AminoAcid(17, 'V', "Val", "Valine");
    public static final AminoAcid W    = new AminoAcid(18, 'W', "Trp", "Tryptophan");
    public static final AminoAcid Y    = new AminoAcid(19, 'Y', "Tyr", "Tyrosine");
    // ambiguous characters
    public static final AminoAcid B    = new AminoAcid(20, 'B', "Asx", "Aspartic acid or Asparagine"); // N or D
    public static final AminoAcid Z    = new AminoAcid(21, 'Z', "Glx", "Glutamine or Glutamic acid"); // Q or E
    public static final AminoAcid X    = new AminoAcid(22, 'X', "Xaa", "Any amino acid");
    public static final AminoAcid STOP = new AminoAcid(23, '*', "*", "Stop");

    // indel characters
    public static final AminoAcid GAP  = new AminoAcid(24, '-', "-", "GAP");
    public static final AminoAcid NULL = new AminoAcid(25, '.', "???", "error");

    public static final AminoAcid[] REGULAR_RESIDUES = {A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y,B,Z,X,STOP};
    public static final AminoAcid[] ALL_RESIDUES = {A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y,B,Z,X,STOP,GAP,NULL};

    public static final int getIndex(AminoAcid aa){return aa.getIndex();}
    public final int getIndex(){return index;}

    public static final AminoAcid getResidue(int index){return ALL_RESIDUES[index];}

    /**
     * With the exception of null and gap, all others are OK
     * @return false if this is gap or null
     */
    public final boolean isRegularResidue()
    {
        return isRegularResidue(this);
    }

    /**
     * With the exception of null and gap, all others are OK
     * @param res an amino acid
     * @return false if gap or null
     */
    public static final boolean isRegularResidue(AminoAcid res)
    {
        return (res != GAP && res != NULL && res != null);

    }

    private static final AminoAcid decode[];

    static
    {
        decode = new AminoAcid[256];
        java.util.Arrays.fill(decode, NULL);
        for (AminoAcid aa:ALL_RESIDUES)
        {
            decode[aa.code]=aa;
        }
    }
    /**
     * Finds out which amino acid (or gap or null) corresponds to the code
     * @param code IUBMB common alpha-amino acid code, gap (<code>-</code>), stop (<code>*</code>), or NULL (<code>.</code>)
     * @return the corresponding amino acid
     */
    public static AminoAcid getResidue(char code)
    {
        return decode[code];
    }

    /**
     * String description of this object
     *
     * @return string description
     */
    @Override
    public String toString()
    {
        StringBuilder sb = new StringBuilder("AA");
        sb.append("['");
        sb.append(code);
        sb.append("' idx ");
        sb.append(index);
        sb.append(']');
        return sb.toString();
    }

    /**
     * One-letter code for the amino acid
     *
     * @return one-letter code
     */
    public char getCode(){return code;}
    /**
     * Three-letter code for the amino acid
     * @return three-letter code
     */
    public String getShortName(){return shortName;}
    /**
     * Full name of the amino acid
     * @return full name
     */
    public String getLongName(){return longName;}

    /**
     * Alphabetical ordering by on one-letter code.
     * @param obj an AminoAcid object
     * @return an integer that can be used to order amino acids
     */
    @Override
    public int compareTo(AminoAcid obj)
    {
        return code - obj.code;
    }    
    
}