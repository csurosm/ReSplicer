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
 * All-purpose methods for a sliding DNA window. Provides
 * encoding/decoding on forward and reverse strands.
 * 
 * @author Mikl&oacute;s Cs&#369;r&ouml;s
 */
public class DNAWindow 
{
    /**
     * Window length
     */
    private int window_length;

    /**
     * Content of the window on forward strand
     */
    private long forward_content;

    /**
     * Encoding of the window content on the reverse strand: this is
     * the complement inverse of the forward strand content.
     */
    private long reverse_content;

    /**
     * Bit mask used to exclude the most significant bits in the encoding.
     */
    protected long bit_mask;

    /**
     * Number of rightmost positions without ambiguous nucleotides
     */
    private int valid_positions;

    public DNAWindow(int window_length)
    {
        if (window_length<1 || window_length>31)
            throw new IllegalArgumentException("Window length must be between 1 and 31");
        this.window_length = window_length;
        bit_mask=0L;
        for (int i=0; i<window_length-1; i++)
            bit_mask = (bit_mask << 2) | 3L;
        reset();
    }


    /**
     * Initialization of the window at the start of a new sequence.
     * @param content is an array of which the rightmost (larger index) bytes determine the window content
     */
    protected void setContent(byte[] content)
    {
        reset();
        for (int j=0; j<content.length; j++)
        {
            byte nuc = content[j];
            nextNucleotide(nuc);
        }
    }

    /**
     * Encoding of the current window content.
     * The encoding is by using two bits per nucleotide
     * (<code>IDX_</code><var>c</var> constants in
     * {@link NucleotideEncoding}), and little endian
     * combination into the <code>long</code> value.
     * The first nucleotide is encoded on the lowest two bits,
     * the next nucleotide on the two next bits, etc.
     *
     * @return <code>long</code> encoding of the window content on the forward strand
     */
    public long getForwardContent()
    {
        return forward_content;
    }

    /**
     * Encoding of the current window content.
     * The encoding is by using two bits per nucleotide
     * (<code>IDX_</code><var>c</var> constants in
     * {@link NucleotideEncoding}), and little endian
     * combination into the <code>long</code> value.
     * The first nucleotide is encoded on the lowest two bits,
     * the next nucleotide on the two next bits, etc.
     *
     * @return <code>long</code> encoding of the window content on the reverse strand
     */
    public long getReverseContent()
    {
        return reverse_content;
    }

    protected final void reset()
    {
        forward_content=reverse_content=0L;
        valid_positions = 0;
    }

    /**
     * Window length
     * @return window length
     */
    public int getWindowLength()
    {
        return window_length;
    }


    /**
     * Translates an encoded window content into a DNA sequence.
     *
     * @param content window content
     * @return String description of the content
     */
    public final String toDNASequence(long content)
    {
        return new String(toDNASequence(content, window_length));
    }

    private static char[] toDNASequence(long content, int length)
    {
        char[] retval = new char[length];
        for (int i=length-1; i>=0; i--)
        {
            int idx = (int)(content & 3L);
            char nuc_code = DNASequence.toChar(idx);
            retval[length-1-i]=nuc_code;
            content = content >>> 2;
        }
        return retval;
    }


    /**
     * Updates window content when the window slides forward by one position.
     *
     * @param nuc the new nucleotide covered by the window
     */
    protected void nextNucleotide(byte nuc)
    {
        if (DNASequence.isAmbiguous(nuc))
        {
            forward_content = reverse_content=0L;
            valid_positions = 0;
        }
        else
        {
            int f = DNASequence.toIndex(nuc);
            int r = DNASequence.complement(f);
            forward_content = (forward_content >>> 2) | ((long)f << 2*(window_length-1));
            reverse_content = ((reverse_content & bit_mask)<< 2) | r;
            if (valid_positions<window_length) valid_positions++;
        }
    }


    /**
     * Whether the window content is <em>valid</em>. 
     * The window is valid if it is filled and 
     * all nucleotides in the window are non-ambiguous.
     *
     * @return true if the window is valid
     */
    public boolean isValid()
    {
        return valid_positions == window_length;
    }    
}
