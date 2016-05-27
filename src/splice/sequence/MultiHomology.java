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

import java.util.HashMap;
import java.util.Map;
import java.util.Set;

/**
 *
 * Collection of multiple diagonals for alignment. 
 * Checks consistency of pairwise diagonals: 
 * basically, if x maps to y and y maps to z, then x's mapping to z is implied. 
 * 
 * @author Mikl&oacute;s Cs&#369;r&ouml;s
 */
public class MultiHomology 
{
    /**
     * All coordinates are calculated relative to a single (arbitrary) sequence stored in {@link #origin}
     */
    private DNAStrand origin=null;
    /** 
     * Diagonal offsets with respect to {@link #origin} coordinates
     */
    private final Map<DNAStrand, Integer> strand_coordinate_offset;
    
    public MultiHomology()
    {
        strand_coordinate_offset = new HashMap<>();
    }
    
    public MultiHomology(DNAStrand origin)
    {
        this();
        setOrigin(origin);
    }
    
    public int size()
    {
        return strand_coordinate_offset.size();
    }
    
    public Set<DNAStrand> getSequences()
    {
        return strand_coordinate_offset.keySet();
    }
    
    public boolean contains(DNAStrand seq)
    {
        return strand_coordinate_offset.containsKey(seq);
    }
    
    public int getAbsolutePosition(DNAStrand seq, int pos)
    {
        return projectCoordinate(seq, pos, origin);
    }
    
    public int projectAbsolutePosition(DNAStrand seq_to, int abs_pos)
    {
        return projectCoordinate(origin, abs_pos, seq_to);
    }
    
    /**
     * Called when the first sequence is known.
     * 
     * @param seq selected origin sequence (no other sequences were added yet)
     */
    private void setOrigin(DNAStrand seq)
    {
        assert (origin==null);
        
        this.origin = seq;
        strand_coordinate_offset.put(origin, 0);
    }
    
    /**
     * Adds a pair of matching positions. Verifies if known mappings do not contradict this mapping.  
     * 
     * @param seq1 first sequence with strand
     * @param pos1 first sequence position
     * @param seq2 second sequence with strand (order is immaterial between the sequences)
     * @param pos2 second sequence position
     * @return whether this match can be added consistently to the set of matching positions already known
     */
    public boolean addMatch(DNAStrand seq1, int pos1, DNAStrand seq2, int pos2)
    {
        if (origin==null)
        {
            setOrigin(seq1);
        }
        if (seq1.equals(origin))
        {
            return addProjection(pos1, seq2, pos2);
        } else if (seq2.equals(origin))
        {
            return addProjection(pos2, seq1, pos1);
        } else
        {
            if (strand_coordinate_offset.containsKey(seq1))
            {
                int pos0 = projectCoordinate(seq1, pos1, origin);
                return addProjection(pos0, seq2, pos2);
            } else if (strand_coordinate_offset.containsKey(seq2))
            {
                int pos0 = projectCoordinate(seq2, pos2, origin);
                return addProjection(pos0, seq1, pos1);
            } else
                throw new IllegalArgumentException("addMatch called with two novel sequences ("+seq1+", "+seq2+"): known sequences are "+strand_coordinate_offset.keySet());
        }
    }
    
    /**
     * Match between origin and another sequence
     * 
     * @param pos_origin position on {@link #origin}
     * @param seq_to other sequence
     * @param pos position on other sequence that aligns with <var>pos_origin</var>
     * @return whether this mapping can be added consistently
     */
    private boolean addProjection(int pos_origin, DNAStrand seq_to, int pos_to)
    {
        if (strand_coordinate_offset.containsKey(seq_to)) // this we know already
        {
            int projected_to = projectCoordinate(origin, pos_origin, seq_to);
            return (projected_to==pos_to);
        } else
        {
            int diagonal_offset = seq_to.getStrand()==origin.getStrand()
                    ?pos_to-pos_origin
                    :pos_to+pos_origin;
            strand_coordinate_offset.put(seq_to, diagonal_offset);
            return true;
        }
    }
    
//    /**
//     * Offset for mapping positions from one sequence to another. 
//     * 
//     * @param seq_from sequence from which the positions are mapped
//     * @param seq_to sequence to which the positions are mapped
//     * @return an offset value used in {@link #projectCoordinate(splice.sequence.DNAStrand, int, splice.sequence.DNAStrand) }
//     */
//    private int getDiagonalOffset(DNAStrand seq_from, DNAStrand seq_to)
//    {
//        // have (x,y)  and (x',z') : want y->z
//        //
//        //  X   Y   Z   X=Y X=Z Y=Z df
//        //  F   F   F   t   t   t   +-   z'' = (z'-x')+x''   = (z'-x')+((x-y)+y'')   = (z'-x'+x-y) + y''
//        //  F   F   R   t   f   f   --   z'' = (z'+x')-x''   = (z'+x')-((x-y)+y'')   = (z'+x'-x+y) - y''
//        //  F   R   F   f   t   f   ++   z'' = (z'-x')+x''   = (z'-x')+((x+y)-y'')   = (z'-x'+x+y) - y''
//        //  F   R   R   f   f   t   -+   z'' = (z'+x')-x''   = (z'+x')-((x+y)-y'')   = (z'+x'-x-y) + y''
//        //  R   F   F   f   f   t   -+   z'' = (z'+x')-x''   = (z'+x')-((x+y)-y'')   = (z'+x'-x-y) + y''
//        //  R   F   R   f   t   f   ++   z'' = (z'-x')+x''   = (z'-x')+((x+y)-y'')   = (z'-x'+x+y) - y''
//        //  R   R   F   t   f   f   --   z'' = (z'+x')-x''   = (z'+x')-((x-y)+y'')   = (z'+x'-x+y) - y''
//        //  R   R   R   t   t   t   +-   z'' = (z'-x')+x''   = (z'-x')+((x-y)+y'')   = (z'-x'+x-y) + y''
//        //
//        // d_to = d(x,z) = {z'-x' / z'+x'}
//        // d_from = d(x,y) = {y-x / y+x } 
//
//        
//        //  feq  teq
//        //  t   t   dt-df
//        //  t   f   dt+df
//        //  f   t   dt+df
//        //  f   f   dt-df
//        
//        int d_to = strand_coordinate_offset.get(seq_to);  // == 0 if seq_to==origin
//        int d_from = strand_coordinate_offset.get(seq_from);    // == 0 if seq_from==origin
//        
//        
//        return (seq_to.getStrand()==seq_from.getStrand())
//                ?d_to-d_from
//                :d_to+d_from;
//    
//    }
    
    /**
     * Mapping coordinates from one sequence to another. 
     * 
     * @param seq_from sequence from which position is mapped
     * @param pos_from genomic coordinate 
     * @param seq_to sequence to which the position is mapped
     * @return position on the second sequence that aligns with <var>pos_from</var> according to known diagonal mappings
     */
    public int projectCoordinate(DNAStrand seq_from, int pos_from, DNAStrand seq_to)
    {
        int d_to = strand_coordinate_offset.get(seq_to);  // == 0 if seq_to==origin
        int d_from = strand_coordinate_offset.get(seq_from);    // == 0 if seq_from==origin
        
        
        return (seq_to.getStrand()==seq_from.getStrand())
                ?d_to-d_from+pos_from
                :d_to+d_from-pos_from;
        
//        // if seq_from==origin, then d_from==0, and we return d_to
//        // if seq_to==origin, then d_to==0, and we return either -d_from (same strandedness) or d_from (opposite strandedness)
//
//        // projecting on (x,y) diagonal 
//        //  X   Y   projection
//        //  F   F   y' = y+(x'-x) = (y-x)+x'    x' = (x-y)+y' 
//        //  F   R   y' = y-(x'-x) = (y+x)-x'    x' = (x+y)-y'
//        //  R   F   y' = y+(x-x') = (y+x)-x'    x' = (x+y)-y'
//        //  R   R   y' = y-(x-x') = (y-x)+x'    x' = (x-y)+y'
//        return (seq_from.getStrand()==seq_to.getStrand()?diagonal_offset+pos_from:diagonal_offset-pos_from);
    }
    
    @Override
    public String toString()
    {
        StringBuilder sb = new StringBuilder("MH");
        sb.append("[");
        sb.append("orig ").append(origin.getDNA().getOrganism()).append(origin.getStrand());
        sb.append(", size ").append(size());
        if (size()>1)
        {
            sb.append(" {");
            for (DNAStrand member: getSequences())
                if (!member.equals(origin))
                {
                    sb.append(member.getDNA().getOrganism()).append(member.getStrand()).append("/").append(strand_coordinate_offset.get(member)).append(" ");
                }
            sb.append("}");
        }
        sb.append("]");
        return sb.toString();
    }
}
