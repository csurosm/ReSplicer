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
 * Common methods for Fasta-encoded sequences. 
 * 
 * @author Mikl&oacute;s Cs&#369;r&ouml;s
 */
public class FastaSequence 
{
    private final String defline;
    public FastaSequence(String defline)
    {
        this.defline = defline;
        setOrganism (getTagValue(DEFLINE_TAG_ORGANISM));
    }
    
    private String organism;
    public String getOrganism(){return organism;}
    public void setOrganism(String s){this.organism = s;}

    public String getDefLine()
    {
        return defline;
    }
    
    protected void append(String s)
    {}
    
    /**
     * The sequence identifier comes from the defline up to the first space character. 
     * @return identifier for this sequence
     */
    public String getIdent()
    {
        int first_white_space = defline.indexOf(' ');
        if (first_white_space<0)
            return defline;
        else 
            return defline.substring(0,first_white_space);
    }
    
    /**
     * Extracts the value for tags in the syntax of /tag=value
     * 
     * @param tag is what we are interested in, without the leading `/'
     *
     * @return null if there is no such tag
     */
    public final String getTagValue(String tag)
    {
        int tag_loc = defline.indexOf("/"+tag+"=");
        if (tag_loc==-1)
            return null;
        tag_loc += 2+tag.length();
        StringBuilder sb=new StringBuilder();
        int name_length = defline.length();
        char quote = '0';
        for (int i=0; tag_loc+i<name_length; i++)
        {
            char c = defline.charAt(tag_loc+i);
            if (i==0 && (c=='\'' || c=='"'))
                quote = c;
            else 
            {
                if (quote == '0')
                {
                    if (Character.isWhitespace(c))
                        break;
                    else 
                        sb.append(c);
                } else
                {
                    if (c == quote)
                        break;
                    else
                        sb.append(c);
                }
            }
        }
        return sb.toString();
    }
    
    /**
     * Generic pretty print for fasta files.
     * 
     * @param chars_per_line characters per line 
     * @param defline definition line, without starting &gt; (will be prepended)
     * @param sequence sequence to be printed, any characters are allowed
     * @return formatted string  
     */
    public static String toFasta(int chars_per_line, String defline, String sequence)
    {
        StringBuilder sb = new StringBuilder(">");
        sb.append(defline);
        sb.append("\n");
        int pos = 0;
        int len = sequence.length();
        while (pos<len)
        {
            
            int chunk_length = Math.min(chars_per_line, len-pos);
            String ss = sequence.substring(pos, pos+chunk_length);
            sb.append(ss);
            pos += chunk_length;
            if (pos<len)
                sb.append("\n");
        }
        return sb.toString();
    }    
    
    /**
     * Tag used to define the organism in the defline (<code>/organism=xxx</code>).
     */
    public static final String DEFLINE_TAG_ORGANISM = "organism";
        
    
}
