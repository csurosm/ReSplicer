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
 * Defines AGCT&lt;-&gt;0..3 encoding and useful constants for character-, byte- and integer-encoding of nucleotides.
 *
 * There are three types of constants:
 * <ul>
 * <li>int constants starting with IDX_ used for encoding of indexes into arrays: range of 0..4, with {@link #IDX_GAP} included.
 * <li>byte constants starting with BYTE_: four-bit encoding for possibly ambiguous nucleotides for storing large sequences in byte[] arrays;
 *         lower two bits are for unambiguous nucleotides: same encoding as IDX_
 * <li>char constants starting with CHAR_: character transcription of nucleotides; upper case only
 * </ul>
 *
 * @author Mikl&oacute;s Cs&#369;r&ouml;s
 */
public interface NucleotideEncoding
{

    public static final int IDX_A=0; // A-G and C-T have the same first bit: transition (usually more frequent that transversion) changes the second bit only
    public static final int IDX_C=2;
    public static final int IDX_G=1;
    public static final int IDX_T=3;
    public static final int IDX_GAP=4;

    public static final char CHAR_A= 'A';
    public static final char CHAR_C= 'C';
    public static final char CHAR_G= 'G';
    public static final char CHAR_T= 'T';
    public static final char CHAR_N= 'N';
    public static final char CHAR_X= 'X';

    public static final char CHAR_R= 'R'; // Purine: A or G
    public static final char CHAR_Y= 'Y'; // Pyrimidine: C or T
    public static final char CHAR_W= 'W'; // Weak: A or T
    public static final char CHAR_S= 'S'; // Strong: C or G
    public static final char CHAR_K= 'K'; // G or T
    public static final char CHAR_M= 'M'; // A or C

    public static final char CHAR_H= 'H'; // not G
    public static final char CHAR_B= 'B'; // not A
    public static final char CHAR_V= 'V'; // not T
    public static final char CHAR_D= 'D'; // not D

    public static final char CHAR_U= 'U'; // synonym for T

    public static final char CHAR_GAP = '-';
    public static final char CHAR_INS = ' ';

    // can verify ambiguity by (b & 12 == 4)
    // byte 0 is left empty : it is always a bug if you get that byte.
    public static final byte BYTE_A=(byte) IDX_A+0x04;
    public static final byte BYTE_C=(byte) IDX_C+0x04;
    public static final byte BYTE_G=(byte) IDX_G+0x04;
    public static final byte BYTE_T=(byte) IDX_T+0x04;
    public static final byte BYTE_R=(byte) 0x08; // Purine: A or G
    public static final byte BYTE_Y=(byte) 0x09; // Pyrimidine: C or T
    public static final byte BYTE_W=(byte) 0x0a; // Weak: A or T
    public static final byte BYTE_S=(byte) 0x0b; // Strong: C or G
    public static final byte BYTE_K=(byte) 0x0c; // G or T
    public static final byte BYTE_M=(byte) 0x0d; // A or C
    public static final byte BYTE_N=(byte) 0x0f;
    public static final byte BYTE_B=(byte) 0x01; // not A
    public static final byte BYTE_D=(byte) 0x02; // not C
    public static final byte BYTE_H=(byte) 0x03; // not G
    public static final byte BYTE_V=(byte) 0x0e; // not T
    public static final byte BYTE_ERROR=(byte) 0;
    
}
