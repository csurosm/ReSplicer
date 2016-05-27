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

import java.util.Arrays;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Map;
import java.util.Set;

/**
 *
 * Hashtable-based counter.
 * 
 * @author Mikl&oacute;s Cs&#369;r&ouml;s
 * @param <O> type of items being counted
 */
public final class Counter<O> 
{
    private final Map<O, Integer> counts;
    private int tot_count;
    
    public Counter()
    {
        this.counts = new HashMap<>();
        this.tot_count = 0;
    }
    
    /**
     * Registers one occurrence of an item.
     * 
     * @param item occurrence
     */
    public void increment(O item)
    {
        increment(item, 1);
    }
    
    public void increment(O item, int diff)
    {
        int n = getCount(item);
        counts.put(item, n+diff);
        tot_count += diff;
    }
    
    /**
     * Number of times the item was seen (i.e., {@link #increment(java.lang.Object) } was called).
     * 
     * @param item occurrence
     * @return 0 if never seen, otherwise number of occurrences
     */
    public int getCount(O item)
    {
        return (counts.containsKey(item)?counts.get(item):0);
    }
    
    /**
     * Resets all counters to zero.
     * 
     */
    public void clear()
    {
        counts.clear();
        tot_count = 0;
    }

    /**
     * Sets the count for a given item.
     * 
     * @param item occurrence
     * @param cnt count
     */
    public void set(O item, int cnt)
    {
        assert (cnt>=0);
        if (counts.containsKey(item))
            tot_count-=counts.get(item);
        counts.put(item, cnt);
        tot_count += cnt;
    }
    
    /**
     * Array of items sorted by decreasing order of frequency (most frequent first).
     * 
     * @param proxy sets the return type
     * @return sorted array of items
     */
    public O[] toArray(O[] proxy)
    {
        O[] sorted_keys = keySet().toArray(proxy);
        Arrays.sort(sorted_keys, new Comparator<O>() {
            @Override
            public int compare(O o1, O o2) 
            {
                return Integer.compare(getCount(o2), getCount(o1));
            }
        });
        return sorted_keys;
    }
    
    /**
     * Set of items seen (different arguments for {@link #increment(java.lang.Object) }).
     * 
     * @return set of items
     */
    public Set<O> keySet()
    {
        return counts.keySet();
    }

    /**
     * Number of different items counted.
     * 
     * @return number of different items
     */
    public int keyCount()
    {
        return counts.size();
    }
    
    /**
     * Total of all counts.
     * Takes O(1) time.
     * 
     * @return sum of all {@link #getCount(java.lang.Object) } values
     */
    public int totalCounts()
    {
        return tot_count;
    }
}
